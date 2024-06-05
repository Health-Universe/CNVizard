"""
File used to build individual reference files for plotting with CNVizard
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import argparse
import os
import pandas as pd
from cnvizard.reference_builder.data_processing import prepare_cnv_table, explode_cnv_table

def create_reference_files(path_to_input, ngs_type, path_to_output, omim_path, reference_type, starting_letter):
    # Define run name
    run_name = os.path.basename(os.path.normpath(path_to_input))

    # Load OMIM reference
    omim_all = os.path.join(omim_path, 'omim.txt')
    omim_df = pd.read_csv(omim_all, delimiter="\t")

    # Collect all relevant subdirectories
    list_of_dirs = [os.path.join(path_to_input, d) for d in os.listdir(path_to_input) if d.startswith(starting_letter)]

    # Collect individual dataframes
    list_of_dfs = []
    for path in list_of_dirs:
        current_sample = os.path.basename(path)
        if reference_type == "normal":
            total_file = os.path.join(path, 'results/CNV', f'{current_sample}.cnr') if ngs_type == "WES" else os.path.join(path, 'exome_extract/CNV', f'{current_sample}.cnr')
        else:
            total_file = os.path.join(path, 'results/CNV', f'{current_sample}_bintest.tsv') if ngs_type == "WES" else os.path.join(path, 'exome_extract/CNV', f'{current_sample}_bintest.tsv')

        current_total_df = pd.read_csv(total_file, delimiter="\t")
        list_of_dfs.append(current_total_df)

    # Concatenate all individual dataframes to a reference dataframe
    reference = pd.concat(list_of_dfs)
    exploded_reference = explode_cnv_table(reference)
    ordered_reference = prepare_cnv_table(exploded_reference, omim_df)

    # Cast the OMIMG column to string
    ordered_reference['OMIMG'] = ordered_reference['OMIMG'].astype(str).str.split('.').str[0]

    # Define export name
    export_name_parquet = os.path.join(path_to_output, f'CNV_reference_{"bintest_" if reference_type != "normal" else ""}{run_name}.parquet')

    # Export reference to parquet file
    ordered_reference.to_parquet(export_name_parquet, index=False)

def main():
    # Create command-line parser
    parser = argparse.ArgumentParser(description="Create Reference Parquet")

    # Required positional arguments
    parser.add_argument('-i', type=str, help='Input path for CNVkit .cnr files')
    parser.add_argument('-n', type=str, help='NGS method (WES or WGS)')
    parser.add_argument('-o', type=str, help='Output path')
    parser.add_argument('-r', type=str, help='Path to omim_reference')
    parser.add_argument('-t', type=str, help='Type of reference (bintest or normal)')
    parser.add_argument('-s', type=str, help='Starting letter of subdir containing the sample')

    # Get command-line arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    create_reference_files(args.i, args.n, args.o, args.r, args.t, args.s)

if __name__ == "__main__":
    main()
