"""
Filter Selected Ion Chromatograms (SICs)

This script filters Selected Ion Chromatograms (SICs) TSV files generated from PHRP & MASIC Merged output files. It extracts specific columns from the input files, performs data processing, and saves the filtered data as TSV files in the specified output directory.

Usage:
    python prepare_SICs.py -i <input_directory> -o <output_directory>

Arguments:
    -i, --input_dir: Path to the input directory containing PHRP & MASIC Merged output files.
    -o, --output_dir: Path to the output directory to save processed files.

Input Format:
    - The input directory should contain SIC files with names ending in '_PlusSICStats.txt'.

Output Format:
    - Filtered TSV files are saved in the output directory with the same name as the input files, but with the '.tsv' extension.

Dependencies:
    - Python 3.x
    - pandas

Author:
    [Tulasi Rao Relangi, PhD]
"""

import os
import sys
import argparse
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

columns_list = ['Scan', 'PrecursorMZ', 'DelM_PPM', 'Peptide', 'Protein','MSGFScore','MSGFDB_SpecEValue', 'EValue', 'QValue', 'PepQValue', 'ElutionTime', 'ParentIonIntensity','PeakArea','StatMomentsArea']  # Defined columns
file_extension = '_PlusSICStats'

def preprocess_SICs(input_dir, output_dir, columns_list, file_extension, delimiter='\t'):
    """
    Extract selected columns from SIC files in a directory with specific names
    and write the results as TSV.
    """
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate over each file in the input directory
    for file_name in os.listdir(input_dir):
        # Check if the file matches the specified pattern
        if file_extension in file_name:
            print(f"Processed: {file_name} ✔️")
            # Construct the full path of the input file
            input_file_path = os.path.join(input_dir, file_name)
            # Construct the full path of the output file
            output_file_path = os.path.join(output_dir, file_name.replace('_fht_PlusSICStats.txt', '_PlusSICStats.tsv'))
            # Read the file into a DataFrame
            df = pd.read_csv(input_file_path, delimiter=delimiter)

            # Group by some identifier (e.g., Peptide) and select the row with the highest intensity as monoisotopic peak
            monoisotopic_peak_indices = df.groupby('Scan')['ParentIonIntensity'].idxmax()
            monoisotopic_peaks = df.loc[monoisotopic_peak_indices]

            # Select the desired columns
            selected_columns = monoisotopic_peaks[columns_list]

            #Renaming Scan and SpecEValue columns to be compatible with other scripts
            selected_columns.rename(columns={'MSGFDB_SpecEValue': 'SpecEValue', 'Scan': 'ScanNum'}, inplace=True)

            # Write the selected columns to the output file
            selected_columns.to_csv(output_file_path, sep='\t', index=False)
            print(f"Exported to: {output_file_path} 📄")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filters Selected Ion Chromatograms (SICs) tsv files.")
    parser.add_argument("-i", "--input_dir", help="Input directory containing PHRP & MASIC Merged output files.", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory to save processed files.", required=True)
    args = parser.parse_args()

    preprocess_SICs(args.input_dir, args.output_dir, columns_list, file_extension)
