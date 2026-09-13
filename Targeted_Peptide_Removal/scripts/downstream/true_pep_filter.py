#!/usr/bin/env python3
import pandas as pd
import sys
import re
import os
import time
import argparse
import warnings
warnings.filterwarnings('ignore')

import plot_functions
import fdr_estimator

DEFAULT_THRESHOLDS = {
    'FDR': 0.01,
    'MSMSScore': 10.0,
    'QValue': 0.01
}

# ---------------------- Utilities ----------------------
def _none_if_string_none(val):
    if val is None:
        return None
    if isinstance(val, str) and val.strip().lower() in ("none", ""):
        return None
    return val

def _parse_optional_float(val):
    val = _none_if_string_none(val)
    if val is None:
        return None
    return float(val)

def ensure_isdecoy(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure IsDecoy exists; compute from Protein prefix 'XXX_' if missing."""
    if 'IsDecoy' not in df.columns:
        df['IsDecoy'] = df['Protein'].astype(str).str.startswith('XXX_')
    else:
        if df['IsDecoy'].dtype != bool:
            df['IsDecoy'] = df['IsDecoy'].astype(str).str.lower().isin(['true', '1', 't', 'yes'])
    return df

_pep_core_pat = re.compile(r'^[A-Z\-]\.(.+)\.[A-Z\-]$')
def peptide_core(peptide: str) -> str:
    """
    Extract core sequence from PHRP-style 'K.ABCDE.R' and strip mod markers.
    Keeps only A–Z (drops '*', digits, parentheses, +/-). If no flanks, strip non-letters.
    """
    s = str(peptide)
    m = _pep_core_pat.match(s)
    core = m.group(1) if m else s
    core = re.sub(r'[^A-Z]', '', core.upper())
    return core

def add_peptide_core(df: pd.DataFrame) -> pd.DataFrame:
    if 'PeptideCore' not in df.columns:
        df['PeptideCore'] = df['Peptide'].map(peptide_core)
    return df

def stringent_filter(df):
    """Preset: FDR <= 0.01, MSMSScore >= 10, QValue <= 0.05 (applies only to present columns)."""
    conds = []
    if 'FDR' in df.columns:         conds.append(df['FDR'] <= 0.01)
    if 'MSMSScore' in df.columns:   conds.append(df['MSMSScore'] >= 10)
    if 'QValue' in df.columns:      conds.append(df['QValue'] <= 0.05)
    if not conds:
        return df
    mask = conds[0]
    for c in conds[1:]:
        mask = mask & c
    return df[mask]

def apply_metric_filter(df, metric, threshold):
    """
    Apply a single metric filter:
      - FDR/QValue: keep rows <= threshold
      - MSMSScore:  keep rows >= threshold
      - stringent_filter: ignore threshold and apply preset
      Missing metric column => pass-through.
    """
    if metric == 'stringent_filter':
        return stringent_filter(df)

    if threshold is None:
        threshold = DEFAULT_THRESHOLDS[metric]

    if metric == 'FDR' and 'FDR' in df.columns:
        return df[df['FDR'] <= threshold]
    elif metric == 'QValue' and 'QValue' in df.columns:
        return df[df['QValue'] <= threshold]
    elif metric == 'MSMSScore' and 'MSMSScore' in df.columns:
        return df[df['MSMSScore'] >= threshold]
    else:
        return df

def peptide_level_view(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse to one row per peptide core for peptide-level stats:
      - Best (min) QValue per core (if present)
      - Peptide IsDecoy = True only if all assignments are decoy (conservative)
    """
    g = df.groupby('PeptideCore', as_index=False)
    out = g.agg({
        'QValue': (lambda x: pd.to_numeric(x, errors='coerce').min()) if 'QValue' in df.columns else 'first',
        'IsDecoy': 'all'
    })
    return out.rename(columns={'PeptideCore':'Peptide'})

# ---------------------- Core processing ----------------------
def process_directory(input_directory, output_directory,
                      metric1, threshold1,
                      metric2=None, threshold2=None):
    """
    1) Read *_withsyn.tsv files
    2) Apply primary metric filter (defaults allowed)
    3) Optional second filter if BOTH provided (or 'stringent_filter')
    4) Add PeptideCore; set PeptideFlanked = original Peptide; set Peptide = PeptideCore
    5) Only drop truly identical rows (no (Peptide,Protein) collapsing)
    6) Update cumulative stats on peptide cores & unique proteins
    7) Final trim: drop decoys + intensity floor; write as *_filtered.tsv (hide PeptideCore)
    """
    total_files_processed = 0
    rows_processed = 0
    peptides_processed = 0
    proteins_processed = 0
    target_peptides_written = 0

    all_assignments = []
    all_peptides_for_stats = []

    start_time = time.time()
    os.makedirs(output_directory, exist_ok=True)

    # Normalize optional second-filter args
    metric2 = _none_if_string_none(metric2)
    threshold2 = _parse_optional_float(threshold2)
    apply_second = (metric2 is not None) and (threshold2 is not None) if metric2 != 'stringent_filter' else True

    input_files = [f for f in os.listdir(input_directory) if f.endswith('_withsyn.tsv')]
    if not input_files:
        input_files = [f for f in os.listdir(input_directory) if f.endswith('.tsv')]

    for filename in sorted(input_files):
        input_file = os.path.join(input_directory, filename)

        # Output name: *_filtered.tsv
        stem = (filename[:-len('_withsyn.tsv')] if filename.endswith('_withsyn.tsv')
                else filename[:-len('.tsv')] if filename.endswith('.tsv')
                else filename)
        output_file = os.path.join(output_directory, f"{stem}_filtered.tsv")

        df = pd.read_csv(input_file, sep='\t', dtype=str)

        for col in ('Peptide','Protein'):
            if col not in df.columns:
                raise ValueError(f"{filename} missing required column: {col}")

        # coerce numerics if present
        for c in ('FDR','QValue','MSMSScore','ParentIonIntensity'):
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors='coerce')

        df = ensure_isdecoy(df)
        df = df[~df['Protein'].astype(str).str.contains('Contaminant', case=False, na=False)]

        # ---- primary filter ----
        df_filtered = apply_metric_filter(df, metric1, threshold1)

        # ---- optional second filter ----
        if apply_second:
            if metric2 not in (None, 'FDR', 'MSMSScore', 'QValue', 'stringent_filter'):
                raise ValueError("Invalid --filter2: choose FDR, MSMSScore, QValue, or stringent_filter (or 'None')")
            if metric2 == 'stringent_filter':
                df_filtered = stringent_filter(df_filtered)
            elif metric2 is not None and threshold2 is not None:
                df_filtered = apply_metric_filter(df_filtered, metric2, threshold2)

        # ---- peptide relabeling (ONLY change you asked for) ----
        df_filtered = add_peptide_core(df_filtered)                # adds PeptideCore
        df_filtered["PeptideFlanked"] = df_filtered["Peptide"]     # keep original flanked
        df_filtered["Peptide"] = df_filtered["PeptideCore"]        # use core as the peptide label

        # ---- exact-duplicate removal only (no (Peptide,Protein) collapsing) ----
        df_filtered = df_filtered.drop_duplicates()

        # ---- accumulate for stats (use PeptideCore) ----
        all_assignments.append(df_filtered[['Peptide','PeptideCore','Protein','IsDecoy','QValue']].copy())
        pep_stat_view = peptide_level_view(df_filtered)
        all_peptides_for_stats.append(pep_stat_view)

        total_files_processed += 1
        rows_processed += len(df_filtered)

        # ---- write per-file output: targets only + intensity floor; hide PeptideCore ----
        out_df = df_filtered[df_filtered['IsDecoy'] == False].copy()
        if 'ParentIonIntensity' in out_df.columns:
            out_df = out_df[out_df['ParentIonIntensity'] >= 5e3]

        # hide internal stats column, keep PeptideFlanked for reference
        if 'PeptideCore' in out_df.columns:
            out_df = out_df.drop(columns=['PeptideCore'])

        # count unique peptides written using the (core) Peptide column
        target_peptides_written += out_df['Peptide'].nunique()


        # Sort: MSMSScore (desc), then Protein (A→Z) as tiebreak
        if 'MSMSScore' in out_df.columns and 'Protein' in out_df.columns:
            out_df = out_df.sort_values(
                by=['MSMSScore', 'Protein'],
                ascending=[False, True],
                na_position='last',
                kind='mergesort'  # stable
            )
        elif 'Protein' in out_df.columns:
            out_df = out_df.sort_values('Protein', ascending=True, kind='mergesort')

        out_df.to_csv(output_file, sep='\t', index=False)

    # ----- cumulative stats -----
    if all_assignments:
        combined_df = pd.concat(all_assignments, ignore_index=True)
        combined_pep = pd.concat(all_peptides_for_stats, ignore_index=True).drop_duplicates(subset=['Peptide'])

        peptides_processed = combined_pep['Peptide'].nunique()
        proteins_processed = combined_df['Protein'].nunique()

        decoy_peptides = int(combined_pep['IsDecoy'].sum())
        cumulative_pep_fdr_pct = (decoy_peptides / max(1, peptides_processed)) * 100.0

        try:
            peptide_fdr = fdr_estimator.calculate_peptide_fdr(combined_pep)
        except Exception:
            peptide_fdr = round(cumulative_pep_fdr_pct, 4)

        prot_df = (combined_df
                   .assign(IsDecoy=combined_df['Protein'].astype(str).str.startswith('XXX_'))
                   .drop_duplicates(subset=['Protein'])[['Protein','IsDecoy']])
        try:
            protein_fdr = fdr_estimator.calculate_protein_fdr(prot_df)
        except Exception:
            protein_fdr = round(100.0 * prot_df['IsDecoy'].sum() / max(1, len(prot_df)), 4)

        # plots (assignment-level)
        plot_dir = os.path.join(output_directory, "combined_plots")
        try:
            plot_functions.plot_density_msms_vs_is_decoy(combined_df, plot_dir)
            plot_functions.plot_density_ppm_vs_is_decoy(combined_df, plot_dir)
        except Exception:
            pass
    else:
        combined_df = pd.DataFrame()
        peptide_fdr = protein_fdr = 0.0
        peptides_processed = proteins_processed = 0

    end_time = time.time()
    execution_time = end_time - start_time

    # ------- console summary -------
    print("#files processed:", total_files_processed)
    print("#rows (peptide–protein assignments):", rows_processed)
    print("#unique peptides:", peptides_processed, "@", f"{peptide_fdr}%")
    print("#unique proteins:", proteins_processed, "@", f"{protein_fdr}%")
    print("-" * 60)
    print()
    print("#Target peptides written (non-decoy, intensity≥5e3):", target_peptides_written)
    print("#### peptide-level FDR (conservative):", f"{(decoy_peptides / max(1, peptides_processed) * 100.0) if peptides_processed else 0.0:.2f}%")
    print()
    print("-" * 60)
    print("Execution time:", "{:.2f}".format(execution_time), "seconds")

# ---------------------- CLI ----------------------
def build_parser():
    parser = argparse.ArgumentParser(
        description="Filter MSGF+ withsyn tables; relabel Peptide=core and carry PeptideFlanked; robust peptide/protein stats."
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Input directory containing *_withsyn.tsv files (fallback: any .tsv).")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Output directory to save *_filtered.tsv files.")
    parser.add_argument("-f", "--filter_function",
                        choices=['FDR', 'MSMSScore', 'QValue', 'stringent_filter'],
                        default='FDR',
                        help="Primary filter to apply. Default: FDR.")
    parser.add_argument("-t", "--threshold", type=float,
                        help="Primary threshold. Defaults: FDR=0.01, QValue=0.01, MSMSScore=10. "
                             "Ignored for 'stringent_filter'.")
    parser.add_argument("--filter2", type=str, default="None",
                        help="Second filter: FDR | MSMSScore | QValue | stringent_filter | None. "
                             "If 'None', no second pass is applied.")
    parser.add_argument("--threshold2", type=str, default="None",
                        help="Threshold for second filter (float) or 'None'. "
                             "If 'None' (or empty), the second pass is disabled. "
                             "Ignored only if filter2 is 'None'.")
    return parser

def display_help_and_exit(parser):
    parser.print_help()
    print("\nExamples:")
    print("  python filter_withsyn.py -i fdr_esti/ -o filtered_out/ -f FDR -t 0.01 --filter2 None --threshold2 None")
    print("  python filter_withsyn.py -i fdr_esti/ -o filtered_out/ -f MSMSScore -t 10 --filter2 QValue --threshold2 0.01")
    sys.exit(0)

def main():
    if len(sys.argv) == 1:
        parser = build_parser()
        display_help_and_exit(parser)

    parser = build_parser()
    args = parser.parse_args()

    process_directory(
        input_directory=args.input_dir,
        output_directory=args.output_dir,
        metric1=args.filter_function,
        threshold1=args.threshold,
        metric2=args.filter2,
        threshold2=args.threshold2
    )

if __name__ == "__main__":
    main()
