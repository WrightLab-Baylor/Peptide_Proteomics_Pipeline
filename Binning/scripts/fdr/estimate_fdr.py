#!/usr/bin/env python3

import os
import time
import argparse
import numpy as np
import pandas as pd
#import plot_functions


def calculate_psm_fdr(msgfplus_output_file, score_column='SpecEValue', decoy_prefix='XXX_'):
    df = pd.read_csv(msgfplus_output_file, sep='\t')

    df['IsDecoy'] = df['Protein'].astype(str).str.startswith(decoy_prefix)

    # Stable sort preserves input order for exact score ties, making the
    # cumulative target/decoy calculation reproducible without introducing
    # an arbitrary secondary scientific ranking.
    df = df.sort_values(
        by=score_column,
        ascending=True,
        kind='mergesort',
    ).reset_index(drop=True)

    df['cum_decoys'] = df['IsDecoy'].cumsum()
    df['cum_targets'] = (~df['IsDecoy']).cumsum()

    # Classic target-decoy style cumulative FDR
    df['FDR'] = df['cum_decoys'] / df['cum_targets'].replace(0, np.nan)
    df['FDR'] = df['FDR'].clip(upper=1.0)
    df['FDR'] = df['FDR'][::-1].cummin()[::-1]

    df['MSMSScore'] = -np.log10(df[score_column].replace(0, 1))
    df['absPPM'] = abs(df['DelM_PPM'])

    return df[[
        'ScanNum', 'Peptide', 'Protein', 'absPPM',
        'MSGFScore', 'SpecEValue', 'EValue', 'QValue',
        'MSMSScore', 'PepQValue', 'ElutionTime',
        'ParentIonIntensity', 'StatMomentsArea',
        'IsDecoy', 'FDR'
    ]]


def calculate_peptide_fdr(msgfplus_df, score_column_name='SpecEValue'):
    if msgfplus_df.empty:
        return "0.00"

    df = msgfplus_df.sort_values(by=score_column_name, kind='mergesort').copy()
    df['IsDecoy'] = df['Protein'].astype(str).str.startswith('XXX_')

    pep_df = df.groupby('Peptide').first().reset_index()

    target_hits = pep_df[~pep_df['IsDecoy']].shape[0]
    decoy_hits = pep_df[pep_df['IsDecoy']].shape[0]

    peptide_fdr = decoy_hits / (target_hits + decoy_hits) if (target_hits + decoy_hits) > 0 else 0.0
    return "{:.2f}".format(min(peptide_fdr, 1.0) * 100)


def calculate_protein_fdr(msgfplus_df, score_column_name='SpecEValue'):
    if msgfplus_df.empty:
        return "0.00"

    df = msgfplus_df.sort_values(by=score_column_name, kind='mergesort').copy()
    df['IsDecoy'] = df['Protein'].astype(str).str.startswith('XXX_')

    prot_df = df.groupby('Protein').first().reset_index()

    target_hits = prot_df[~prot_df['IsDecoy']].shape[0]
    decoy_hits = prot_df[prot_df['IsDecoy']].shape[0]

    protein_fdr = decoy_hits / (target_hits + decoy_hits) if (target_hits + decoy_hits) > 0 else 0.0
    return "{:.2f}".format(min(protein_fdr, 1.0) * 100)


def output_name_for(filename):
    if filename.endswith("_PlusSICStats.tsv"):
        return filename.replace("_PlusSICStats.tsv", "_fdrstats.tsv")

    if filename.endswith("_consensus_culled.tsv"):
        return filename.replace("_consensus_culled.tsv", "_fdrstats.tsv")

    if filename.endswith("_classical_guided.tsv"):
        return filename.replace("_classical_guided.tsv", "_fdrstats.tsv")
        
    return None


def process_directory(directory_path, output_directory):
    total_files_processed = 0
    psms_processed = 0
    peptides_processed = 0
    proteins_processed = 0
    cumulative_decoys = 0
    cumulative_targets = 0

    start_time = time.time()

    os.makedirs(output_directory, exist_ok=True)

    combined_df = pd.DataFrame()

    for filename in os.listdir(directory_path):
        output_name = output_name_for(filename)

        if output_name is None:
            continue

        total_files_processed += 1

        input_file_path = os.path.join(directory_path, filename)
        output_file_path = os.path.join(output_directory, output_name)

        df = calculate_psm_fdr(input_file_path)
        df.to_csv(output_file_path, sep='\t', index=False)

        combined_df = pd.concat([combined_df, df], ignore_index=True)

        psms_processed += df['ScanNum'].nunique()
        peptides_processed += df['Peptide'].nunique()
        proteins_processed += df['Protein'].nunique()
        cumulative_decoys += int(df['IsDecoy'].sum())
        cumulative_targets += int((~df['IsDecoy']).sum())

    if cumulative_targets + cumulative_decoys > 0:
        cumulative_fdr_percentage = cumulative_decoys / (cumulative_targets + cumulative_decoys) * 100
    else:
        cumulative_fdr_percentage = 0.0

    cumulative_fdr_percentage = "{:.2f}".format(cumulative_fdr_percentage)

    peptide_fdr = calculate_peptide_fdr(combined_df)
    protein_fdr = calculate_protein_fdr(combined_df)

#    if not combined_df.empty:
#        plot_functions.plot_density_msms_vs_is_decoy(
#            combined_df,
#            os.path.join(output_directory, "combined_plots")
#        )
#        plot_functions.plot_density_ppm_vs_is_decoy(
#            combined_df,
#            os.path.join(output_directory, "combined_plots")
#        )

    execution_time = time.time() - start_time

    print("Execution time:", "{:.2f}".format(execution_time), "seconds")
    print()
    print("#Spectrum files:", total_files_processed)
    print("#PSMs:", psms_processed, "@", cumulative_fdr_percentage, "%")
    print("#peptides:", peptides_processed, "@", peptide_fdr, "%")
    print("#proteins:", proteins_processed, "@", protein_fdr, "%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates FDR on processed MSGFPlus/SIC output files.")
    parser.add_argument("-i", "--input_dir", help="Input directory containing processed files.", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory to save FDR files.", required=True)
    args = parser.parse_args()

    process_directory(args.input_dir, args.output_dir)