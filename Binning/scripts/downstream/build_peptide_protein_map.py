#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd

# ---------------------------
# CLI
# ---------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=("Build per-sample protein coverage from an annotated peptide crosstab. "
                     "Identity is PeptideFlanked; coverage marks only the core peptide span, "
                     "requiring flank matches in the protein sequence.")
    )
    p.add_argument("-i", "--tsv_file", required=True,
                   help="Input annotated crosstab TSV (must include Peptide, PeptideFlanked, Protein, and sample columns).")
    p.add_argument("-db", "--database-dir", required=True,
                   help="Directory containing exactly one FASTA file (*.fasta or *.faa).")
    p.add_argument("-p", "--min-num-pep", type=int, default=1,
                   help="Minimum (flanked) peptides per protein to keep (default: 1).")
    p.add_argument("--mode", choices=["unique_only", "requires_unique", "all_matches"],
                   default="all_matches",
                   help="Peptide inclusion mode. 'unique_only' keeps Unique==True only; "
                        "'requires_unique' keeps proteins with any unique peptide; "
                        "'all_matches' keeps all.")
    p.add_argument("-f", "--filter-mode", choices=["per_sample", "global"], default="per_sample",
                   help="Apply min peptide threshold per sample or globally across all samples.")
    p.add_argument("--outdir", default="map_files",
                   help="Output directory for per-sample mapfiles, merged coverage TSV, and heatmap (default: map_files).")
    # Column names (fixed for this pipeline)
    p.add_argument("--peptide-col", default="Peptide", help="Core peptide column name (default: Peptide).")
    p.add_argument("--peptide-flanked-col", default="PeptideFlanked",
                   help="Flanked peptide column name (default: PeptideFlanked).")
    p.add_argument("--protein-col", default="Protein", help="Protein column name (default: Protein).")
    return p.parse_args()

# ---------------------------
# FASTA utilities
# ---------------------------

def find_single_fasta(db_dir: Path) -> Path:
    if not db_dir.is_dir():
        raise FileNotFoundError(f"Database directory not found: {db_dir}")
    candidates = list(db_dir.glob("*.fasta")) + list(db_dir.glob("*.faa"))
    if len(candidates) == 0:
        raise FileNotFoundError(f"No FASTA found in {db_dir} (expect one .fasta or .faa).")
    if len(candidates) > 1:
        names = ", ".join(sorted(p.name for p in candidates))
        raise RuntimeError(f"Multiple FASTAs found in {db_dir}: {names} — keep exactly one.")
    return candidates[0]

def load_fasta_sequences(fasta_path: Path) -> dict:
    """
    Return {protein_id -> sequence}, where protein_id is the header's first token (up to first whitespace).
    """
    seqs = {}
    acc = None
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                header = line[1:].strip()
                acc = header.split()[0]
                seqs.setdefault(acc, [])
            else:
                if acc is not None:
                    seqs[acc].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}

# ---------------------------
# Flank-aware mapping helpers
# ---------------------------

def parse_peptide_flanked(pflanked: str):
    """
    Parse 'N.PEPTIDE.C' → (N, PEPTIDE, C). '-' means terminus (no residue).
    Returns (n_flank_or_None, core, c_flank_or_None).
    """
    parts = str(pflanked).split(".")
    if len(parts) != 3:
        return None, str(pflanked), None
    n, core, c = parts
    n = None if n in ("", "-") else n.upper()
    c = None if c in ("", "-") else c.upper()
    return n, core, c

def find_flanked_positions(protein_seq: str, core: str, n_flank: str | None, c_flank: str | None):
    """
    Return 0-based start indices where 'core' occurs AND adjacent residues satisfy flank constraints.
    Flanks are a gate only; they are not counted as covered.
    """
    seq = (protein_seq or "").upper()
    core_u = (core or "").upper()
    L = len(core_u)
    if L == 0:
        return []
    hits, start = [], 0
    while True:
        idx = seq.find(core_u, start)
        if idx == -1:
            break
        n_ok = True if n_flank is None else (idx > 0 and seq[idx - 1] == n_flank)
        c_ok = True if c_flank is None else (idx + L < len(seq) and seq[idx + L] == c_flank)
        if n_ok and c_ok:
            hits.append(idx)
        start = idx + 1
    return hits

def mark_core_coverage(protein_seq: str, core: str, starts: list[int], cov: np.ndarray):
    """Set cov[start : start+len(core)] = True for each valid start."""
    L = len(core or "")
    if L == 0:
        return
    for s in starts:
        cov[s:s+L] = True

def build_coverage_with_flanks(protein_seq: str, rows_for_protein: pd.DataFrame) -> np.ndarray:
    """
    Given a protein sequence and the subset of rows for that protein containing
    ['PeptideFlanked','Peptide'], return a boolean coverage vector that marks
    the CORE peptide spans only, but requires matching flanks when provided.
    """
    cov = np.zeros(len(protein_seq or ""), dtype=bool)
    if rows_for_protein.empty:
        return cov
    for pflanked, core in rows_for_protein[["PeptideFlanked", "Peptide"]].drop_duplicates().itertuples(index=False):
        n_flank, core_seq, c_flank = parse_peptide_flanked(pflanked)
        if not core_seq:
            core_seq = core
        starts = find_flanked_positions(protein_seq, core_seq, n_flank, c_flank)
        mark_core_coverage(protein_seq, core_seq, starts, cov)
    return cov

# ---------------------------
# I/O checks & sample detection
# ---------------------------

RESERVED_NON_SAMPLE = {
    "Peptide", "PeptideFlanked", "Protein",
    "Cluster", "Unique", "Unique to Cluster", "Shared",
    "Gene", "Function", "nested"
}

def ensure_required_columns(df: pd.DataFrame, peptide_col: str, peptide_flanked_col: str, protein_col: str):
    missing = [c for c in (peptide_col, peptide_flanked_col, protein_col) if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. "
                         f"Expected at least {peptide_col}, {peptide_flanked_col}, {protein_col}.")
    for c in (peptide_col, peptide_flanked_col, protein_col):
        df[c] = df[c].astype(str).str.strip()

def detect_sample_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in RESERVED_NON_SAMPLE]

# ---------------------------
# Per-sample prep & coverage
# ---------------------------

def prepare_subtable_for_sample(df: pd.DataFrame, sample_col: str, mode: str) -> pd.DataFrame:
    """
    Build per-sample subtable with ['PeptideFlanked','Peptide','Protein', sample_col],
    apply intensity > 0, optional uniqueness gating, and dedup by (Protein, PeptideFlanked).
    """
    cols = ["PeptideFlanked", "Peptide", "Protein", sample_col]
    sub = df[cols].copy()
    sub[sample_col] = pd.to_numeric(sub[sample_col], errors="coerce").fillna(0)
    sub = sub[sub[sample_col] > 0]

    if mode == "unique_only":
        if "Unique" not in df.columns:
            raise ValueError("Requested mode 'unique_only' but 'Unique' column is missing.")
        sub = sub.loc[df.loc[sub.index, "Unique"].values]

    elif mode == "requires_unique":
        if "Unique" not in df.columns:
            raise ValueError("Requested mode 'requires_unique' but 'Unique' column is missing.")
        keep_by_prot = df.groupby("Protein")["Unique"].any()
        sub = sub[sub["Protein"].map(keep_by_prot).fillna(False)]

    if not sub.empty:
        sub = (sub.sort_values(by=sample_col, ascending=False)
                  .drop_duplicates(subset=["Protein", "PeptideFlanked"]))
    return sub

def compute_coverage_for_sample(sub: pd.DataFrame,
                                protein_seqs: dict,
                                min_num_pep: int,
                                filter_mode: str,
                                sample_col: str) -> pd.DataFrame:
    """
    Return per-protein DataFrame with columns:
      Protein, num_pep (distinct PeptideFlanked count), and Coverage for this sample.
    """
    rows = []
    for prot_id, rows_for_prot in sub.groupby("Protein"):
        seq = protein_seqs.get(str(prot_id))
        if not seq:
            # Protein not in FASTA → skip (or warn)
            continue
        cov_vec = build_coverage_with_flanks(seq, rows_for_prot[["PeptideFlanked", "Peptide"]])
        coverage_pct = 100.0 * float(cov_vec.sum()) / max(1, len(seq))
        num_pep = int(rows_for_prot["PeptideFlanked"].nunique())
        rows.append((prot_id, num_pep, coverage_pct))

    out = pd.DataFrame(rows, columns=["Protein", "num_pep", "Coverage"])

    if filter_mode == "per_sample" and min_num_pep > 1 and not out.empty:
        out = out[out["num_pep"] >= min_num_pep]

    return out.rename(columns={"Coverage": sample_col})

# ---------------------------
# Main
# ---------------------------

def main():
    args = parse_args()

    in_path = Path(args.tsv_file)
    db_dir = Path(args.database_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not in_path.is_file():
        print(f"ERROR: Input table not found: {in_path}", file=sys.stderr)
        sys.exit(1)

    try:
        fasta_path = find_single_fasta(db_dir)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    protein_seqs = load_fasta_sequences(fasta_path)

    # Load input and validate columns
    df = pd.read_csv(in_path, sep="\t", low_memory=False, dtype=str)
    ensure_required_columns(df, args.peptide_col, args.peptide_flanked_col, args.protein_col)

    sample_cols = detect_sample_columns(df)
    if not sample_cols:
        print("ERROR: No sample columns detected (only ID/annotation columns present).", file=sys.stderr)
        sys.exit(1)

    # Per-sample coverage tables and (optional) global counts
    per_sample_tables = []   # for merging: each has only ["Protein", <sample_col>]
    global_counts = {}       # Protein -> set of PeptideFlanked across all samples

    for sample_col in sample_cols:
        sub = prepare_subtable_for_sample(df, sample_col, args.mode)

        if not sub.empty:
            # update global counts
            for prot, pfl in sub[["Protein", "PeptideFlanked"]].drop_duplicates().itertuples(index=False):
                global_counts.setdefault(prot, set()).add(pfl)

        cov_tbl = compute_coverage_for_sample(
            sub, protein_seqs, args.min_num_pep, args.filter_mode, sample_col
        )

        # 1) Write the per-sample mapfile with num_pep + Coverage
        out_map = outdir / f"{sample_col}_mapfile.tsv"
        map_tbl = cov_tbl.rename(columns={sample_col: "Coverage"})  # Protein, num_pep, Coverage
        map_tbl.to_csv(out_map, sep="\t", index=False)

        # 2) Keep ONLY Protein + coverage for merging across samples
        merge_tbl = cov_tbl[["Protein", sample_col]].copy()
        per_sample_tables.append(merge_tbl)

    # Global peptide threshold gating (applied after collecting all samples)
    merged = None
    for tbl in per_sample_tables:
        merged = tbl.copy() if merged is None else merged.merge(tbl, on="Protein", how="outer")

    if merged is None or merged.empty:
        print("No proteins passed filters; nothing to write.", file=sys.stderr)
        sys.exit(0)

    if args.filter_mode == "global" and args.min_num_pep > 1:
        keep = {prot for prot, s in global_counts.items() if len(s) >= args.min_num_pep}
        merged = merged[merged["Protein"].isin(keep)]

    # Tidy & write merged matrix
    sample_cols_sorted = sorted([c for c in merged.columns if c != "Protein"])
    merged = merged[["Protein"] + sample_cols_sorted]
    for c in sample_cols_sorted:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0)

    out_matrix = outdir / "grouped_coverage.tsv"
    merged.to_csv(out_matrix, sep="\t", index=False)

    # Always generate heatmap (top 100 by max coverage if large)
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        data = merged.set_index("Protein")[sample_cols_sorted]
        if len(data) > 100:
            top_idx = data.max(axis=1).sort_values(ascending=False).head(100).index
            data = data.loc[top_idx]
        plt.figure(figsize=(10, max(6, len(data) * 0.15)))
        sns.heatmap(data, cmap="viridis")
        plt.tight_layout()
        pdf_path = outdir / "coverage_heatmap.pdf"
        plt.savefig(pdf_path, bbox_inches="tight")
        plt.close()
        print(f"Wrote per-sample mapfiles: {outdir}")
        print(f"Wrote merged coverage matrix: {out_matrix}  (proteins: {len(merged)})")
        print(f"Wrote heatmap: {pdf_path}")
    except Exception as e:
        print(f"Heatmap generation skipped due to error: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()

