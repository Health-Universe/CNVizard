"""
File which is used to convert panel lists (.tsv format) from the genomics england "panel app" to genelist.txt files
@author: Jeremias Krause / Matthias Begemann / Florian Kraft / Carlos Classen 
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd
import argparse

# Create Commandline parser
parser = argparse.ArgumentParser(description="Create Reference Parquet")
# Required positional arguments -i for to to .cnr files, -n for NGS-Type (WES/WGS, -o for output location and -r for path to omim.txt file)
parser.add_argument(
    "-i", type=str, help="Positional argument for input genepanel.txt file"
)
parser.add_argument("-o", type=str, help="Positional argument for output.txt file")


# Get Commandline arguments by using parser
args = parser.parse_args()

# Load parsed arguments into corresponding variables
path_to_input = args.i
path_to_output = args.o

# read panel app list
panel_app_list = pd.read_csv(path_to_input, sep="\t")
# extract gene list from panel app dataframe
gene_list = panel_app_list["Gene Symbol"]
gene_list = gene_list.drop_duplicates()
# write extracted gene list into newly created .txt file
gene_list.to_csv(path_to_output, index=False, header=False)
