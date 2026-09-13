"""
Script Name: fdr_estimator.py

Description:
This script calculates the False Discovery Rate (FDR) for processed MSGFPlus output files. It processes each file in the input directory, calculates FDR for peptide-spectrum matches (PSMs), and generates corresponding output files with FDR statistics. Additionally, it calculates FDR for proteins based on the target-decoy competition method.

Functions:
1. calculate_fdr:
    Calculates FDR for PSMs in a MSGFPlus output file.
    Parameters:
        - msgfplus_output_file: Path to the MSGFPlus output file.
        - score_column_name (optional): Name of the column containing the scoring values (default: 'SpecEValue').
        - protein_column_name (optional): Name of the column containing protein IDs (default: 'Protein').
    Returns:
        DataFrame: DataFrame containing FDR statistics for PSMs.

2. calculate_protein_fdr:
    Calculates FDR for proteins based on the target-decoy competition method.
    Parameters:
        - msgfplus_df: DataFrame containing MSGFPlus output.
        - score_column_name (optional): Name of the column containing the scoring values (default: 'SpecEValue').
    Returns:
        float: FDR for proteins as a percentage.

3. process_directory:
    Processes all MSGFPlus output files in the input directory, calculates FDR, and saves the results.
    Parameters:
        - directory_path: Path to the input directory containing processed MSGFPlus output files.
        - output_directory: Path to the output directory to save processed files.
    Returns:
        None

Usage:
python fdr_estimator.py -i <input_directory> -o <output_directory>

Example:
python fdr_estimator.py -i processed_files/ -o fdr_results/
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import plot_functions

def calculate_psm_fdr(msgfplus_output_file, score_column='SpecEValue', decoy_prefix='XXX_'):
    # Load MSGFPlus output
    df = pd.read_csv(msgfplus_output_file, sep='\t')

    # Label decoys
    df['IsDecoy'] = df['Protein'].str.startswith(decoy_prefix)

    # Sort by score (lower SpecEValue = better)
    df = df.sort_values(by=score_column, ascending=True).reset_index(drop=True)

    # Cumulative decoys and targets
    df['cum_decoys'] = df['IsDecoy'].cumsum()
    df['cum_targets'] = (~df['IsDecoy']).cumsum()

    # FDR and q-value
    df['FDR'] = df['cum_decoys'] / df['cum_targets'].replace(0, np.nan)
    df['FDR'] = df['FDR'].clip(upper=1.0)
    df['FDR'] = df['FDR'][::-1].cummin()[::-1]

    # Optional: MSMS score
    df['MSMSScore'] = -np.log10(df[score_column].replace(0, 1))

    # Calculate MSMS score
    df['MSMSScore'] = -np.log10(df['SpecEValue'])

    # Calculate absMassError
    df['absPPM'] = (abs(df['DelM_PPM']))

    return df[['ScanNum', 'Peptide', 'Protein','absPPM','MSGFScore', 'SpecEValue', 'EValue', 'QValue', 'MSMSScore','PepQValue', 'ElutionTime', 'ParentIonIntensity', 'StatMomentsArea','IsDecoy','FDR']]

def calculate_peptide_fdr(msgfplus_df, score_column_name='SpecEValue'):
    # Sort the DataFrame by score
    df = msgfplus_df.sort_values(by=score_column_name)

    # Identify decoy proteins
    df['IsDecoy'] = df['Protein'].str.startswith('XXX_') 

    # Group by protein and keep only the best-scoring hit for each protein
    pep_df = df.groupby('Peptide').first().reset_index()

    # Count the number of target and decoy hits
    target_hits = pep_df[~pep_df['IsDecoy']].shape[0]
    decoy_hits = pep_df[pep_df['IsDecoy']].shape[0]

    # Calculate FDR using target-decoy competition method
    peptide_fdr = decoy_hits / (target_hits + decoy_hits) if target_hits > 0 else 0.0

    return "{:.2f}".format(min(peptide_fdr, 1.0) * 100)

def calculate_protein_fdr(msgfplus_df, score_column_name='SpecEValue'):
    # Sort the DataFrame by score
    df = msgfplus_df.sort_values(by=score_column_name)

    # Identify decoy proteins
    df['IsDecoy'] = df['Protein'].str.startswith('XXX_')

    # Group by protein and keep only the best-scoring hit for each protein
    prot_df = df.groupby('Protein').first().reset_index()

    # Count the number of target and decoy hits
    target_hits = prot_df[~prot_df['IsDecoy']].shape[0]
    decoy_hits = prot_df[prot_df['IsDecoy']].shape[0]

    # Calculate FDR using target-decoy competition method
    protein_fdr = decoy_hits / (target_hits + decoy_hits) if target_hits > 0 else 0.0

    # Cap FDR at 100% and return as percentage
    return "{:.2f}".format(min(protein_fdr, 1.0) * 100)

def process_directory(directory_path, output_directory):
    total_files_processed = 0
    total_psms = 0
    psms_processed = 0
    peptides_processed = 0
    proteins_processed = 0
    cumulative_decoys = 0
    cumulative_targets = 0

    start_time = time.time()
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    combined_df = pd.DataFrame()

    # Iterate over files in the directory
    for filename in os.listdir(directory_path):
        if filename.endswith("_PlusSICStats.tsv"):
            total_files_processed += 1
            input_file_path = os.path.join(directory_path, filename)
            output_file_path = os.path.join(output_directory, filename.replace("_PlusSICStats.tsv", "_fdrstats.tsv"))

            # Call calculate_fdr for each file and save the result
            df = calculate_psm_fdr(input_file_path)
            df.to_csv(output_file_path, sep='\t', index=False)  # Save DataFrame to file
            combined_df = pd.concat([combined_df, df], ignore_index=True)

            # Update counts
            total_psms += len(df)
            psms_processed += df['ScanNum'].nunique()
            peptides_processed += df['Peptide'].nunique()
            proteins_processed += df['Protein'].nunique()
            cumulative_decoys += sum(df['IsDecoy'])
            cumulative_targets += len(df) - sum(df['IsDecoy'])

    if cumulative_targets > 0:
        cumulative_fdr_percentage = (cumulative_decoys) / (cumulative_targets + cumulative_decoys) * 100
    else:
        cumulative_fdr_percentage = 0.0
    cumulative_fdr_percentage = "{:.2f}".format(cumulative_fdr_percentage)
    peptide_fdr = calculate_peptide_fdr(combined_df)
    protein_fdr = calculate_protein_fdr(combined_df)

    # Plotting combined data
    if not combined_df.empty:
        plot_functions.plot_density_msms_vs_is_decoy(combined_df, os.path.join(output_directory, "combined_plots"))  # Plot combined data
        plot_functions.plot_density_ppm_vs_is_decoy(combined_df, os.path.join(output_directory, "combined_plots"))

    end_time = time.time()
    execution_time = end_time - start_time  # Calculate execution time

    print("Execution time:", "{:.2f}".format(execution_time), "seconds")
    print()
    print("#Spectrum files:", total_files_processed)
    print("#PSMs:", psms_processed, "@", cumulative_fdr_percentage,"%")
    print("#peptides:", peptides_processed, "@", peptide_fdr,"%")
    print("#proteins:", proteins_processed, "@", protein_fdr,"%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calcultes FDR on processed MSGFPlus output files.")
    parser.add_argument("-i", "--input_dir", help="Input directory containing processed MSGFPlus output files.", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory to save processed files.", required=True)
    args = parser.parse_args()

    process_directory(args.input_dir, args.output_dir)