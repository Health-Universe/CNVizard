"""
File used to merge previously created individual reference files for plotting with CNVizard
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import argparse
import os
import pandas as pd
import numpy as np
from cnvizard.reference_builder.data_processing import _get_call_counts, _get_frequency, _get_frequency_bintest

def merge_reference_files(path_to_input, path_to_output, path_to_bintest):
    # Load reference files
    reference_files = [os.path.join(path_to_input, f) for f in os.listdir(path_to_input) if f.endswith('.parquet')]
    reference_dfs = [pd.read_parquet(file) for file in reference_files]

    # Concatenate all reference DataFrames
    reference_df = pd.concat(reference_dfs, ignore_index=True)

    # Drop unnecessary columns
    reference_df.drop(columns=["chromosome", "start", "end", "weight", "squaredvalue", "OMIMG", "Disease", "OMIMP", "Inheritance", "comments"], inplace=True)

    # Group by gene and exon, aggregate columns
    reference_df = reference_df.groupby(["gene", "exon"]).agg({"depth": list, "log2": list, "call": list}).reset_index()

    # Apply call count and frequency functions
    reference_df["call_counts"] = reference_df["call"].apply(_get_call_counts)
    reference_df["het_del_frequency"] = reference_df["call_counts"].apply(lambda x: _get_frequency(x, 0, 4))
    reference_df["hom_del_frequency"] = reference_df["call_counts"].apply(lambda x: _get_frequency(x, 1, 4))
    reference_df["dup_frequency"] = reference_df["call_counts"].apply(lambda x: _get_frequency(x, 2, 4))

    # Calculate additional statistics
    reference_df["max_log2"] = reference_df["log2"].apply(max)
    reference_df["min_log2"] = reference_df["log2"].apply(min)
    reference_df["mean_depth"] = reference_df["depth"].apply(np.mean)
    reference_df["mean_log2"] = reference_df["log2"].apply(np.mean)
    reference_df["median_depth"] = reference_df["depth"].apply(lambda x: np.quantile(x, 0.5))
    reference_df["median_log2"] = reference_df["log2"].apply(lambda x: np.quantile(x, 0.5))
    reference_df["q1_depth"] = reference_df["depth"].apply(lambda x: np.quantile(x, 0.25))
    reference_df["q1_log2"] = reference_df["log2"].apply(lambda x: np.quantile(x, 0.25))
    reference_df["q3_depth"] = reference_df["depth"].apply(lambda x: np.quantile(x, 0.75))
    reference_df["q3_log2"] = reference_df["log2"].apply(lambda x: np.quantile(x, 0.75))
    reference_df["std_depth"] = reference_df["depth"].apply(np.std)
    reference_df["std_log2"] = reference_df["log2"].apply(np.std)
    reference_df["box_size"] = (reference_df["q3_log2"] - reference_df["q1_log2"]) * 1.5
    reference_df["actual_minimum_log2"] = reference_df["q1_log2"] - reference_df["box_size"]
    reference_df["actual_maximum_log2"] = reference_df["q3_log2"] + reference_df["box_size"]
    reference_df["box_size_depth"] = (reference_df["q3_depth"] - reference_df["q1_depth"]) * 1.5
    reference_df["actual_minimum_depth"] = reference_df["q1_depth"] - reference_df["box_size_depth"]
    reference_df["actual_maximum_depth"] = reference_df["q3_depth"] + reference_df["box_size_depth"]

    # Create new DataFrame with call counts from normal references
    ref_counts_df = reference_df[["gene", "exon", "call_counts"]]

    # Load bintest reference files
    bintest_files = [os.path.join(path_to_bintest, f) for f in os.listdir(path_to_bintest) if f.endswith('.parquet')]
    bintest_dfs = [pd.read_parquet(file) for file in bintest_files]

    # Concatenate all bintest reference DataFrames
    bintest_df = pd.concat(bintest_dfs, ignore_index=True)
    bintest_df.drop(columns=["chromosome", "start", "end", "weight", "squaredvalue", "OMIMG", "Disease", "OMIMP", "Inheritance", "comments"], inplace=True)

    # Group by gene and exon, aggregate columns
    bintest_df = bintest_df.groupby(["gene", "exon"]).agg({"depth": list, "log2": list, "call": list}).reset_index()
    bintest_df["call_counts"] = bintest_df["call"].apply(_get_call_counts)

    # Merge with bintest reference DataFrame
    bintest_df = pd.merge(bintest_df, ref_counts_df, on=["gene", "exon"], how="left")
    bintest_df["het_del_frequency"] = bintest_df.apply(lambda x: _get_frequency_bintest(x["call_counts_x"], 0, 4, x["call_counts_y"]), axis=1)
    bintest_df["hom_del_frequency"] = bintest_df.apply(lambda x: _get_frequency_bintest(x["call_counts_x"], 1, 4, x["call_counts_y"]), axis=1)
    bintest_df["dup_frequency"] = bintest_df.apply(lambda x: _get_frequency_bintest(x["call_counts_x"], 2, 4, x["call_counts_y"]), axis=1)
    bintest_df.drop(columns=["depth", "log2", "call", "call_counts_x", "call_counts_y"], inplace=True)

    # Write bintest reference to parquet file
    bintest_df.to_parquet(os.path.join(path_to_output, "cnv_reference_bintest.parquet"))

    # Drop unnecessary columns from reference DataFrame
    reference_df.drop(columns=["depth", "log2", "call", "call_counts", "max_log2", "min_log2", "box_size", "box_size_depth"], inplace=True)

    # Write merged and formatted reference to parquet file
    reference_df.to_parquet(os.path.join(path_to_output, "cnv_reference.parquet"))

def main():
    # Create command-line parser
    parser = argparse.ArgumentParser(description="Merge Reference Parquet")

    # Required positional arguments
    parser.add_argument('-i', type=str, help='Input path for reference .parquet files')
    parser.add_argument('-o', type=str, help='Output path')
    parser.add_argument('-b', type=str, help='Input path for bintest reference .parquet files')

    # Get command-line arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    merge_reference_files(args.i, args.o, args.b)

if __name__ == "__main__":
    main()
