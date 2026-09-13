#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import pandas as pd

def is_decoy(protein: str) -> bool:
    return isinstance(protein, str) and protein.startswith("XXX_")

def load_safe_tsv(path: Path, sep: str = "\t") -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=sep)
    except Exception as e:
        raise RuntimeError(f"Failed reading {path}: {e}")

def detect_scan_col(df: pd.DataFrame) -> str:
    """
    Return the scan column name used in df: 'ScanNum' or 'Scan'.
    Raise if neither present.
    """
    for c in ("ScanNum", "Scan"):
        if c in df.columns:
            return c
    raise ValueError("Missing required scan column (ScanNum or Scan).")

def pick_template_rows(fdr: pd.DataFrame, group_cols=("Peptide",)) -> pd.DataFrame:
    """
    Pick one template row per grouping key (lowest QValue if available).
    We use this to copy spectrum/metadata columns onto appended rows.
    """
    tmp = fdr.copy()
    if "QValue" in tmp.columns:
        tmp["__Q__"] = pd.to_numeric(tmp["QValue"], errors="coerce")
        tmp = tmp.sort_values(list(group_cols) + ["__Q__"], kind="mergesort")
    else:
        tmp["__Q__"] = pd.NA
        tmp = tmp.sort_values(list(group_cols), kind="mergesort")
    best = tmp.groupby(list(group_cols), as_index=False).head(1).drop(columns="__Q__", errors="ignore")
    return best

def append_syn_to_fdr(fdr_path: Path, syn_path: Path, out_suffix: str = "_withsyn.tsv") -> tuple[int, int]:
    fdr = load_safe_tsv(fdr_path, sep="\t")
    syn = load_safe_tsv(syn_path, sep="\t")

    # Required columns
    for col in ("Peptide", "Protein"):
        if col not in fdr.columns:
            raise ValueError(f"{fdr_path.name} missing required column: {col}")
        if col not in syn.columns:
            raise ValueError(f"{syn_path.name} missing required column: {col}")

    # Detect scan columns (can differ across tables)
    fdr_scan = detect_scan_col(fdr)
    syn_scan = detect_scan_col(syn)

    # Light normalization to avoid accidental misses due to stray spaces
    for df in (fdr, syn):
        df["Peptide"] = df["Peptide"].astype(str).str.strip()
        df["Protein"] = df["Protein"].astype(str).str.strip()

    # Build unique triples
    fdr_triples = fdr[[fdr_scan, "Peptide", "Protein"]].drop_duplicates()
    syn_triples = syn[[syn_scan, "Peptide", "Protein"]].drop_duplicates()

    # Only consider SYN rows whose (scan, peptide) exist in FDR
    syn_triples = syn_triples.merge(
        fdr[[fdr_scan, "Peptide"]].drop_duplicates().assign(_present=True),
        left_on=[syn_scan, "Peptide"],
        right_on=[fdr_scan, "Peptide"],
        how="inner"
    ).drop(columns=["_present", fdr_scan])

    # New triples = in SYN but not already present in FDR
    new_triples = syn_triples.merge(
        fdr_triples.assign(_in_fdr=True),
        left_on=[syn_scan, "Peptide", "Protein"],
        right_on=[fdr_scan, "Peptide", "Protein"],
        how="left"
    )
    # Keep only non-existing triples; keep syn_scan for now
    new_triples = new_triples[new_triples["_in_fdr"].isna()].drop(columns=["_in_fdr", fdr_scan], errors="ignore")

    if new_triples.empty:
        out_path = fdr_path.with_name(fdr_path.stem.replace("_fdrstats", "") + out_suffix)
        fdr.to_csv(out_path, sep="\t", index=False)
        return (0, fdr_triples.shape[0])

    # Template per (scan, peptide) from FDR (so we copy spectrum-level fields from the correct instance)
    template = pick_template_rows(fdr, group_cols=(fdr_scan, "Peptide"))

    # Merge new triples with their (scan, peptide) template rows
    new_rows = new_triples.merge(
        template,
        left_on=[syn_scan, "Peptide"],
        right_on=[fdr_scan, "Peptide"],
        how="left",
        suffixes=("", "_tmpl")
    )

    # Ensure the scan column name matches FDR's column naming in output
    # Copy syn scan values into the FDR scan column name, then drop the SYN scan column
    if syn_scan != fdr_scan:
        new_rows[fdr_scan] = new_rows[syn_scan]
    # Replace the template protein with the SYN protein (already in 'Protein').
    # Drop all *_tmpl helper columns and the SYN scan column to avoid duplicates.
    drop_these = [c for c in new_rows.columns if c.endswith("_tmpl")] + ([syn_scan] if syn_scan in new_rows.columns else [])
    new_rows = new_rows.drop(columns=drop_these, errors="ignore")

    # Compute IsDecoy for appended rows (existing FDR rows are untouched)
    if "IsDecoy" in new_rows.columns:
        new_rows["IsDecoy"] = new_rows["Protein"].astype(str).map(is_decoy)
    else:
        # If template lacked IsDecoy, compute it freshly
        new_rows["IsDecoy"] = new_rows["Protein"].astype(str).map(is_decoy)

    # Align appended-row columns to FDR column order (preserve only FDR's columns)
    new_rows = new_rows[[c for c in fdr.columns if c in new_rows.columns]]

    # Append; do not drop duplicates globally (we've already enforced triple-uniqueness)
    combined = pd.concat([fdr, new_rows], ignore_index=True)

    # Write output alongside the original FDR file
    out_name = fdr_path.stem.replace("_fdrstats", "") + out_suffix
    out_path = fdr_path.with_name(out_name)
    combined.to_csv(out_path, sep="\t", index=False)

    # Report: how many appended and now how many unique (scan, peptide, protein) triples
    uniq_cols = [fdr_scan, "Peptide", "Protein"]
    return (len(new_rows), combined[uniq_cols].drop_duplicates().shape[0])

def find_pairs(workdir: Path, fdr_dir: str = "fdr_esti", syn_dir: str = "results/PHRPOut"):
    fdr_root = (workdir / fdr_dir)
    syn_root = (workdir / syn_dir)
    if not fdr_root.is_dir():
        raise FileNotFoundError(f"Missing FDR directory: {fdr_root}")
    if not syn_root.is_dir():
        raise FileNotFoundError(f"Missing SYN directory: {syn_root}")

    pairs = []
    for fdr_file in fdr_root.glob("*_fdrstats.tsv"):
        base = fdr_file.name[:-len("_fdrstats.tsv")]
        syn_file = syn_root / f"{base}_syn.txt"
        if syn_file.exists():
            pairs.append((fdr_file, syn_file))
        else:
            print(f"[WARN] No matching SYN for {fdr_file.name} (looked for {syn_file.name})", file=sys.stderr)
    return pairs

def main():
    ap = argparse.ArgumentParser(
        description="Append peptide→protein combos from PHRP SYN into FDRStats, keyed by (Scan/ScanNum, Peptide, Protein)."
    )
    ap.add_argument("-w", "--workdir", type=Path, default=Path("."), help="Working directory containing fdr_esti and results/PHRPOut")
    ap.add_argument("--fdr-dir", default="fdr_esti", help="Relative path to FDRStats directory")
    ap.add_argument("--syn-dir", default="results/PHRPOut", help="Relative path to SYN directory")
    ap.add_argument("--suffix", default="_withsyn.tsv", help="Suffix for output files")
    args = ap.parse_args()

    pairs = find_pairs(args.workdir, args.fdr_dir, args.syn_dir)
    if not pairs:
        print("No file pairs found. Nothing to do.")
        return

    total_added = 0
    for fdr_path, syn_path in pairs:
        added, unique_triples = append_syn_to_fdr(fdr_path, syn_path, out_suffix=args.suffix)
        print(f"[OK] {fdr_path.name}  +{added} rows  -> unique (scan,peptide,protein) now: {unique_triples}")
        total_added += added

    print(f"\nDone. Total new rows appended across all files: {total_added}")

if __name__ == "__main__":
    main()
