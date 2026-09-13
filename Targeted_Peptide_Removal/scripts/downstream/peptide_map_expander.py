#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd


# ----------------------------
# Helpers
# ----------------------------
def is_decoy(protein: str) -> bool:
    # Keep your existing convention
    return isinstance(protein, str) and protein.startswith("XXX_")


def detect_scan_col(df: pd.DataFrame) -> str:
    for c in ("ScanNum", "Scan"):
        if c in df.columns:
            return c
    raise ValueError("Missing required scan column (ScanNum or Scan).")


def canonicalize_peptide(peptide: str) -> str:
    """
    Convert typical PHRP/MS-GF+ peptide strings into the digest-map key.

    Handles:
      - Flanked form: X.PEPTIDE.Y   -> PEPTIDE
      - Unflanked: PEPTIDE         -> PEPTIDE
      - Removes common mod markers by stripping everything except A–Z
    """
    if peptide is None:
        return ""
    s = str(peptide).strip()

    # If flanked with two dots, pull the middle
    if s.count(".") >= 2:
        first = s.find(".")
        last = s.rfind(".")
        if 0 <= first < last:
            s = s[first + 1:last]

    # Strip anything that isn't a capital letter (mods, charges, brackets, etc.)
    s = re.sub(r"[^A-Z]", "", s)
    return s


def protein_id_from_header(header: str) -> str:
    """
    Extract a stable record identifier from a FASTA header.

    We do NOT assume any specific prefix like 'T_'.
    Convention: use the first whitespace-delimited token.
    """
    return str(header).strip().split()[0] if header is not None else ""


def read_fasta_sequences(path: Path) -> Dict[str, str]:
    """
    Minimal FASTA reader returning {record_key -> sequence}.

    Keys stored:
      - first token of header (protein_id_from_header)
      - full header (minus '>') as an alias key

    This makes the lookup resilient to future header style changes, and to whether
    the 'Protein' column contains a bare accession/id or the full header string.
    """
    if not path.exists():
        raise FileNotFoundError(f"FASTA not found: {path}")

    seqs: Dict[str, str] = {}
    header_full: Optional[str] = None
    header_id: Optional[str] = None
    chunks: List[str] = []

    def flush() -> None:
        nonlocal header_full, header_id, chunks
        if header_full is None:
            return
        seq = "".join(chunks).strip().upper()
        chunks = []
        if not seq:
            return

        hid = header_id or ""
        hf = header_full or ""

        if hid:
            seqs.setdefault(hid, seq)
        if hf:
            seqs.setdefault(hf, seq)

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header_full = line[1:].strip()
                header_id = protein_id_from_header(header_full)
            else:
                chunks.append(line)

    flush()
    return seqs


def infer_parent_fasta_path(
    workdir: Path,
    inputfile_name: str = "inputfile.tsv",
    database_dirname: str = "database",
) -> Optional[Path]:
    """
    Infer the peptide FASTA path from the pipeline's inputfile.tsv.

    Expected inputfile.tsv columns include 'database' (as produced by your automator).
    The FASTA is assumed to live at: workdir / database_dirname / <database_filename>

    Returns:
      Path if inference succeeds and file exists; otherwise None.
    """
    input_path = workdir / inputfile_name
    if not input_path.exists():
        return None

    try:
        df = pd.read_csv(input_path, sep="\t", dtype=str)
    except Exception:
        return None

    if "database" not in df.columns:
        return None

    db_vals = [str(x).strip() for x in df["database"].dropna().tolist() if str(x).strip()]
    if not db_vals:
        return None

    uniq = sorted(set(db_vals))
    if len(uniq) != 1:
        # Ambiguous: multiple database FASTAs referenced in one run.
        return None

    db_name = uniq[0]
    candidate = workdir / database_dirname / db_name
    return candidate if candidate.exists() else None


def load_peptide_map(map_path: Path) -> Dict[str, List[str]]:
    """
    Read peptide_protein_digest_map.tsv and return:
      { canonical_peptide -> sorted list of protein IDs }

    Expected columns: Peptide, ProteinCount, Proteins
    Proteins field is ';'-separated.
    """
    df = pd.read_csv(map_path, sep="\t", dtype=str, keep_default_na=False)

    expected = ["Peptide", "ProteinCount", "Proteins"]
    for c in expected:
        if c not in df.columns:
            raise ValueError(f"{map_path.name} missing required column: {c}")

    df["Peptide"] = df["Peptide"].fillna("").astype(str)
    df["Proteins"] = df["Proteins"].fillna("").astype(str)

    pep2prots: Dict[str, List[str]] = {}

    for _, row in df.iterrows():
        pep = canonicalize_peptide(row["Peptide"])
        prots_field = str(row["Proteins"]).strip()

        if not pep:
            continue

        if not prots_field:
            pep2prots.setdefault(pep, [])
            continue

        prots = [p.strip() for p in prots_field.split(";") if p.strip()]
        prots = sorted(set(prots))
        pep2prots[pep] = prots

    return pep2prots


def expand_fdr_by_map(
    fdr_path: Path,
    pep2prots: Dict[str, List[str]],
    out_path: Path,
    peptide_col: str = "Peptide",
    protein_col: str = "Protein",
    recompute_isdecoy: bool = True,
    dedup_triples: bool = False,
    parent_fasta_seqs: Optional[Dict[str, str]] = None,
    enable_parent_fallback: bool = True,
) -> Dict[str, int | float]:
    """
    For each row in fdrstats:
      - canonicalize peptide
      - look up proteins in map
      - if found and non-empty: emit N rows (one per protein)
      - else:
          - optionally try a "parent peptide" fallback (peptide-FASTA mode):
              * treat the row's Protein field as a FASTA record key
              * fetch the parent peptide sequence from the searched peptide FASTA
              * map the parent peptide via pep2prots and use those proteins
          - if still not found: emit original row unchanged

    Returns:
      Dict of metrics (see code)
    """
    fdr = pd.read_csv(fdr_path, sep="\t", dtype=str)

    for col in (peptide_col, protein_col):
        if col not in fdr.columns:
            raise ValueError(f"{fdr_path.name} missing required column: {col}")

    scan_col = detect_scan_col(fdr)

    # Normalize key cols (avoid whitespace mismatches)
    fdr[peptide_col] = fdr[peptide_col].astype(str).str.strip()
    fdr[protein_col] = fdr[protein_col].astype(str).str.strip()

    out_rows = []
    direct_mapped_input_rows = 0
    fallback_mapped_input_rows = 0
    # Output rows emitted for mapped input rows (sum of proteins emitted per mapped row)
    output_rows_from_mapped = 0

    for _, r in fdr.iterrows():
        pep_key = canonicalize_peptide(r[peptide_col])
        mapped = pep2prots.get(pep_key)

        # Primary mapping path
        if mapped:
            direct_mapped_input_rows += 1
            output_rows_from_mapped += len(mapped)
            for prot in mapped:
                r2 = r.copy()
                r2[protein_col] = prot
                out_rows.append(r2)
            continue

        # Parent-peptide fallback (peptide FASTA mode)
        if ( 
            enable_parent_fallback 
            and parent_fasta_seqs
            and not is_decoy(r[protein_col])
        ):
            prot_key = str(r[protein_col]).strip()
            parent_seq = parent_fasta_seqs.get(prot_key)

            # If Protein contains extra description, try first token too
            if parent_seq is None and prot_key:
                parent_seq = parent_fasta_seqs.get(protein_id_from_header(prot_key))

            if parent_seq:
                parent_pep_key = canonicalize_peptide(parent_seq)
                parent_mapped = pep2prots.get(parent_pep_key)
                if parent_mapped:
                    fallback_mapped_input_rows += 1
                    output_rows_from_mapped += len(parent_mapped)
                    for prot in parent_mapped:
                        r2 = r.copy()
                        r2[protein_col] = prot
                        out_rows.append(r2)
                    continue

        # Default: keep row as-is
        out_rows.append(r)

    out_df = pd.DataFrame(out_rows, columns=fdr.columns)

    # Recompute IsDecoy if present or if you want it guaranteed correct after swapping proteins
    if recompute_isdecoy:
        if "IsDecoy" in out_df.columns:
            out_df["IsDecoy"] = out_df[protein_col].astype(str).map(is_decoy)

    if dedup_triples:
        out_df = out_df.drop_duplicates(subset=[scan_col, peptide_col, protein_col], keep="first")

    out_df.to_csv(out_path, sep="\t", index=False)

    rows_in = int(len(fdr))
    rows_out = int(len(out_df))
    mapped_input_rows = int(direct_mapped_input_rows + fallback_mapped_input_rows)
    unmapped_input_rows = int(rows_in - mapped_input_rows)
    # Extra rows attributable to expansion (vs keeping 1 row per mapped input row)
    # If a peptide maps to 1 protein, it contributes 0 extra rows.
    extra_rows_from_expansion = int(output_rows_from_mapped - mapped_input_rows)
    avg_proteins_per_mapped_row = (output_rows_from_mapped / mapped_input_rows) if mapped_input_rows else 0.0

    return {
        "RowsIn": rows_in,
        "RowsOut": rows_out,
        "MappedInputRows": mapped_input_rows,
        "UnmappedInputRows": unmapped_input_rows,
        "DirectMappedInputRows": int(direct_mapped_input_rows),
        "FallbackMappedInputRows": int(fallback_mapped_input_rows),
        "OutputRowsFromMapped": int(output_rows_from_mapped),
        "ExtraRowsFromExpansion": extra_rows_from_expansion,
        "AvgProteinsPerMappedRow": float(avg_proteins_per_mapped_row),
    }


def find_fdr_files(workdir: Path, fdr_dir: str = "fdr_esti") -> List[Path]:
    root = workdir / fdr_dir
    if not root.is_dir():
        raise FileNotFoundError(f"Missing FDR directory: {root}")
    return sorted(root.glob("*_fdrstats.tsv"))


# ----------------------------
# CLI
# ----------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Expand *_fdrstats.tsv so each peptide has one row per mapped protein from "
            "peptide_protein_digest_map.tsv (no SYN merge). Output naming matches peptide_data_merger. "
            "Includes an optional parent-peptide fallback for peptide-FASTA searches."
        )
    )
    ap.add_argument("-w", "--workdir", type=Path, default=Path("."), help="Working directory (contains fdr_esti and peptide_protein_digest_map.tsv)")
    ap.add_argument("--fdr-dir", default="fdr_esti", help="Relative path to FDRStats directory")
    ap.add_argument("--map-name", default="peptide_protein_digest_map.tsv", help="Digest map TSV filename in workdir")
    ap.add_argument("--suffix", default="_withsyn.tsv", help="Output suffix (keep default to preserve pipeline expectations)")
    ap.add_argument("--peptide-col", default="Peptide", help="Peptide column name in fdrstats (default: Peptide)")
    ap.add_argument("--protein-col", default="Protein", help="Protein column name in fdrstats (default: Protein)")
    ap.add_argument("--no-recompute-isdecoy", action="store_true", help="Do not recompute IsDecoy after protein replacement")
    ap.add_argument("--dedup-triples", action="store_true", help="Drop duplicate (scan, peptide, protein) triples in output")
    ap.add_argument("--audit", default="peptide_map_expansion_audit.tsv", help="Audit TSV filename (written in workdir)")

    # Parent-peptide fallback configuration
    ap.add_argument(
        "--parent-peptide-fasta",
        type=Path,
        default=None,
        help=(
            "Optional: explicitly provide the peptide FASTA used in the MS-GF+ search. "
            "If omitted, the script will attempt to infer it from inputfile.tsv and "
            "assume it is located at workdir/database/<database_fasta_name>."
        ),
    )
    ap.add_argument(
        "--inputfile-name",
        default="inputfile.tsv",
        help="Name of the pipeline inputfile TSV (default: inputfile.tsv)",
    )
    ap.add_argument(
        "--database-dirname",
        default="database",
        help="Database directory under workdir (default: database)",
    )
    ap.add_argument(
        "--no-parent-fallback",
        action="store_true",
        help="Disable parent-peptide fallback even if a peptide FASTA is found/provided.",
    )

    args = ap.parse_args()

    map_path = args.workdir / args.map_name
    if not map_path.is_file():
        print(f"ERROR: Missing map file: {map_path}", file=sys.stderr)
        return 1

    pep2prots = load_peptide_map(map_path)

    # Determine parent peptide FASTA (explicit path wins; otherwise infer from inputfile.tsv)
    parent_fasta_path: Optional[Path] = None
    if args.parent_peptide_fasta is not None:
        parent_fasta_path = args.parent_peptide_fasta
    else:
        parent_fasta_path = infer_parent_fasta_path(
            workdir=args.workdir,
            inputfile_name=args.inputfile_name,
            database_dirname=args.database_dirname,
        )

    parent_seqs: Optional[Dict[str, str]] = None
    if parent_fasta_path is not None and (not args.no_parent_fallback):
        try:
            parent_seqs = read_fasta_sequences(parent_fasta_path)
            print(f"[INFO] Parent-peptide fallback enabled using FASTA: {parent_fasta_path}")
        except Exception as e:
            print(f"[WARN] Could not read inferred/provided peptide FASTA for parent fallback: {e}", file=sys.stderr)
            parent_seqs = None
    else:
        if args.no_parent_fallback:
            print("[INFO] Parent-peptide fallback disabled (--no-parent-fallback).")
        else:
            print("[INFO] No parent-peptide FASTA provided/found; parent fallback will be skipped.")

    fdr_files = find_fdr_files(args.workdir, args.fdr_dir)
    if not fdr_files:
        print("No *_fdrstats.tsv files found. Nothing to do.")
        return 0

    audit_rows = []
    total_in = total_out = 0
    total_mapped = total_unmapped = 0
    total_direct = total_fallback = 0
    total_output_from_mapped = 0
    total_extra_rows = 0

    for fdr_path in fdr_files:
        out_name = fdr_path.stem.replace("_fdrstats", "") + args.suffix
        out_path = fdr_path.with_name(out_name)

        metrics = expand_fdr_by_map(
            fdr_path=fdr_path,
            pep2prots=pep2prots,
            out_path=out_path,
            peptide_col=args.peptide_col,
            protein_col=args.protein_col,
            recompute_isdecoy=(not args.no_recompute_isdecoy),
            dedup_triples=args.dedup_triples,
            parent_fasta_seqs=parent_seqs,
            enable_parent_fallback=(not args.no_parent_fallback),
        )

        audit_rows.append({
            "File": fdr_path.name,
            "RowsIn": metrics["RowsIn"],
            "RowsOut": metrics["RowsOut"],
            "MappedInputRows": metrics["MappedInputRows"],
            "UnmappedInputRows": metrics["UnmappedInputRows"],
            "DirectMappedInputRows": metrics["DirectMappedInputRows"],
            "FallbackMappedInputRows": metrics["FallbackMappedInputRows"],
            "OutputRowsFromMapped": metrics["OutputRowsFromMapped"],
            "ExtraRowsFromExpansion": metrics["ExtraRowsFromExpansion"],
            "AvgProteinsPerMappedRow": metrics["AvgProteinsPerMappedRow"],
            "OutFile": out_path.name,
        })

        total_in += metrics["RowsIn"]
        total_out += metrics["RowsOut"]
        total_mapped += metrics["MappedInputRows"]
        total_unmapped += metrics["UnmappedInputRows"]
        total_direct += metrics["DirectMappedInputRows"]
        total_fallback += metrics["FallbackMappedInputRows"]
        total_output_from_mapped += metrics["OutputRowsFromMapped"]
        total_extra_rows += metrics["ExtraRowsFromExpansion"]

        if parent_seqs:
            fb = metrics["FallbackMappedInputRows"]
            fb_str = f"  fallback-mapped-input-rows: {fb}"
        else:
            fb_str = ""
        print(
            f"[OK] {fdr_path.name} -> {out_path.name}  rows: {metrics['RowsIn']} -> {metrics['RowsOut']}"
            f"  mapped-input-rows: {metrics['MappedInputRows']} (direct: {metrics['DirectMappedInputRows']}, fallback: {metrics['FallbackMappedInputRows']})"
            f"  extra-rows-from-expansion: {metrics['ExtraRowsFromExpansion']}  avg-proteins/mapped-row: {metrics['AvgProteinsPerMappedRow']:.4f}"
            f"{fb_str}"
        )

    audit_path = args.workdir / args.audit
    pd.DataFrame(audit_rows).to_csv(audit_path, sep="\t", index=False)

    print(f"\nDone.")
    print(f"Total rows in : {total_in}")
    print(f"Total rows out: {total_out}")
    print(f"Total mapped input rows: {total_mapped} (direct: {total_direct}, fallback: {total_fallback})")
    print(f"Total unmapped input rows: {total_unmapped}")
    print(f"Total extra rows from expansion: {total_extra_rows}")
    if total_mapped:
        print(f"Avg proteins per mapped row: {total_output_from_mapped/total_mapped:.4f}")
    print(f"Audit written: {audit_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
