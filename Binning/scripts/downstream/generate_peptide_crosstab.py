#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected (True/False).")

def parse_args():
    p = argparse.ArgumentParser(
        description="Build a peptide–protein crosstab across samples. "
                    "Consolidate intensities per (PeptideFlanked, Protein); "
                    "retain both PeptideFlanked and core Peptide columns."
    )
    p.add_argument(
        "-i", "--input_dir", required=True,
        help="Directory containing per-sample *_filtered.tsv files."
    )
    p.add_argument(
        "-o", "--output_tsv", required=True,
        help="Path to write the peptide–protein crosstab TSV."
    )
    p.add_argument(
        "--only-highest-intensity", type=str2bool, default=False,
        help="True → keep only the highest ParentIonIntensity per (PeptideFlanked, Protein). "
             "False → sum intensities per (PeptideFlanked, Protein)."
    )
    p.add_argument(
        "--peptide-col", default="Peptide",
        help="Column name for core peptide sequences (default: Peptide)."
    )
    p.add_argument(
        "--peptide-flanked-col", default="PeptideFlanked",
        help="Column name for flanked peptide sequences (default: PeptideFlanked)."
    )
    p.add_argument(
        "--protein-col", default=None,
        help="Column name for protein IDs. If not provided, auto-detects 'Protein' or 'ProteinID'."
    )
    p.add_argument(
        "--intensity-col", default="ParentIonIntensity",
        help="Column name for intensities (default: ParentIonIntensity)."
    )
    return p.parse_args()

def detect_protein_col(df, user_col):
    if user_col is not None:
        if user_col in df.columns:
            return user_col
        else:
            raise ValueError(f"--protein-col='{user_col}' not found in columns: {list(df.columns)}")
    for cand in ("Protein", "ProteinID"):
        if cand in df.columns:
            return cand
    raise ValueError("Could not find a protein column. Tried 'Protein' and 'ProteinID'. "
                     "Pass --protein-col to specify explicitly.")

def load_and_reduce_one(tsv_path, peptide_col, peptide_flanked_col, protein_col_hint, intensity_col, use_max):
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    protein_col = detect_protein_col(df, protein_col_hint)

    # Required columns
    required = (peptide_col, peptide_flanked_col, protein_col, intensity_col)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns {missing} in {tsv_path}. "
            f"Have: {list(df.columns)}"
        )

    # Sanitize keys & coerce numeric
    for c in (peptide_col, peptide_flanked_col, protein_col):
        df[c] = df[c].astype(str).str.strip()
    df[intensity_col] = pd.to_numeric(df[intensity_col], errors="coerce")

    # Drop rows with nulls in key/intensity cols
    df = df.dropna(subset=[peptide_col, peptide_flanked_col, protein_col, intensity_col])

    # Group by (PeptideFlanked, Protein)
    gb_keys = [peptide_flanked_col, protein_col]
    gb = df.groupby(gb_keys, as_index=False, sort=False)
    if use_max:
        red = gb[intensity_col].max()
    else:
        red = gb[intensity_col].sum()

    # Attach the core Peptide (first observed per group) so we can retain it downstream
    core_map = (
        df[[peptide_flanked_col, protein_col, peptide_col]]
        .drop_duplicates(subset=[peptide_flanked_col, protein_col, peptide_col])
    )
    # If multiple cores ever appear for the same flanked/protein (unlikely), pick the most frequent:
    core_first = (
        df.groupby(gb_keys + [peptide_col])[intensity_col]
          .size()
          .reset_index(name="n")
          .sort_values(["PeptideFlanked", protein_col, "n"], ascending=[True, True, False])
          .drop_duplicates(subset=gb_keys)
          [[peptide_flanked_col, protein_col, peptide_col]]
    )
    # Prefer frequency-based, fallback to simple drop_duplicates if needed
    try:
        red = red.merge(core_first, on=gb_keys, how="left")
    except Exception:
        red = red.merge(core_map, on=gb_keys, how="left")

    # Standardize output columns for this sample-reduced table
    red = red.rename(columns={
        peptide_flanked_col: "PeptideFlanked",
        protein_col: "Protein",
        peptide_col: "Peptide",
        intensity_col: "Intensity"
    })

    return red

def merge_wide(wide_df, red, sample_name):
    # red: columns PeptideFlanked, Protein, Peptide, Intensity (single sample)
    red = red.copy()
    red = red.rename(columns={"Intensity": sample_name})

    if wide_df is None:
        return red

    # Outer merge on composite key (PeptideFlanked, Protein, Peptide)
    out = wide_df.merge(
        red,
        on=["PeptideFlanked", "Protein", "Peptide"],
        how="outer"
    )
    return out

def main():
    args = parse_args()
    in_dir = Path(args.input_dir)
    if not in_dir.is_dir():
        print(f"ERROR: Input directory not found: {in_dir}", file=sys.stderr)
        sys.exit(1)

    # Only pull *_filtered.tsv files
    tsvs = sorted([p for p in in_dir.iterdir() if p.is_file() and p.name.endswith("_filtered.tsv")])
    if not tsvs:
        print(f"ERROR: No files matching '*_filtered.tsv' in {in_dir}", file=sys.stderr)
        sys.exit(1)

    wide = None
    for pth in tsvs:
        sample = pth.stem.replace("_filtered","")  # e.g., 'S_P_UV_filtered' -> this whole stem becomes the column header
        red = load_and_reduce_one(
            pth,
            peptide_col=args.peptide_col,
            peptide_flanked_col=args.peptide_flanked_col,
            protein_col_hint=args.protein_col,
            intensity_col=args.intensity_col,
            use_max=args.only_highest_intensity
        )
        wide = merge_wide(wide, red, sample)

    # Fill missing with 0 (no observation in that sample)
    numeric_cols = [c for c in wide.columns if c not in ("PeptideFlanked", "Protein", "Peptide")]
    wide[numeric_cols] = wide[numeric_cols].fillna(0)

    # Sort for readability (Protein, then PeptideFlanked, then Peptide)
    wide = wide.sort_values(by=["Protein", "PeptideFlanked", "Peptide"], kind="mergesort").reset_index(drop=True)
    
    # Reorder columns: Peptide, PeptideFlanked, Protein, then samples
    sample_cols = [c for c in wide.columns if c not in ("Peptide", "PeptideFlanked", "Protein")]
    wide = wide[["Peptide", "PeptideFlanked", "Protein"] + sample_cols]
    
    # Write output
    out_path = Path(args.output_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote crosstab: {out_path}  (rows: {len(wide)}, samples: {len(numeric_cols)})")

if __name__ == "__main__":
    main()
