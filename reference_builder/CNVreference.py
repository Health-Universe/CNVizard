"""
File which is used to build individual reference files used for plotting with the CNVizard
@author: Jeremias Krause / Matthias Begemann / Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""
import argparse
import os
import pandas as pd
import fastparquet
import pyarrow

def prepare_cnv_table(df:pd.DataFrame,df2:pd.DataFrame) -> pd.DataFrame :
    """
    Function to format and reorder the .cnr DataFrame
    Input Arguments:
        df (pandas DataFrame) -> .cnr DataFrame 
    Output Arguments:
        df (pandas DataFrame) -> reorderded and formatted DataFrame
    """
    # Drop Antitarget entries 
    df.drop(df[df['gene'].str.contains('Antitarget') == True].index, inplace = True)
    #Reverse function to log2 
    df.loc[:, 'squaredvalue'] = 2**df['log2']
    #Split gene and exon file into two seperate fields 
    df['gene'] = df['gene'].str.split('_')
    df.loc[:, 'exon'] = df['gene'].str[1].astype(int)
    df.loc[:, 'gene'] = df['gene'].str[0]
    # Reorder the columns to resemble the familiar cnv-table look 
    cols = df.columns.tolist()
    cols = cols[0:4] + cols[-1:] + cols[5:6] + cols[6:7] + cols[4:5] + cols[7:-1]
    df = df[cols]
    # Insert a call column 
    df.insert(loc=7, column='call', value='')
    # Apply calling scores, based on log2 values (values recommended by CNVkit)
    df.loc[df['log2'] <= -1.1, 'call'] = int(0)
    df.loc[(df['log2'] <= -0.4) & (df['log2'] > -1.1), 'call'] = int(1)
    df.loc[(df['log2'] <= 0.3) & (df['log2'] > -0.4), 'call'] = int(2)
    df.loc[df['log2'] > 0.3, 'call'] = int(3)
    # merge with omim df 
    df = pd.merge(df, df2, on='gene', how='left')
    df.fillna('-', inplace=True)
    df['comments'] = '.'
    return df
    
def explode_cnv_table(df:pd.DataFrame) -> pd.DataFrame :
    """
    Function to split redundant exons and explode the df
    Input Arguments:
        df (pandas DataFrame): .cnr DataFrame
    Output Arguments:
        df (pandas DataFrame): formatted .cnr DataFrame
    """
    df['gene'] = df['gene'].str.split(',')
    df = df.explode('gene')
    return df 

#Create Commandline parser
parser = argparse.ArgumentParser(description="Create Reference Parquet")
#Required positional arguments -i for to to .cnr files, -n for NGS-Type (WES/WGS, -o for output location and -r for path to omim.txt file)
parser.add_argument('-i', type=str,
                    help='Positional argument for input path (CNVkit .cnr files)')
parser.add_argument('-n', type=str,
                    help='Positional argument for NGS method (WES or WGS)')
parser.add_argument('-o', type=str,
                    help='Positional argument for output')
parser.add_argument('-r', type=str,
                    help='Positional argument for path to omim_reference')
parser.add_argument("-t", type=str,
                    help="Positional argument for Type of reference (bintest or normal)")
parser.add_argument("-s", type=str,
                    help="Positional argument for the starting letter of the subdir which contain the sample")

#Get Commandline arguments by using parser
args = parser.parse_args()

#Load parsed arguments into corresponding variables
path_to_input = args.i
NGS_Type = args.n
path_to_output = args.o
omim_path = args.r
reference_type = args.t
starting_letter = args.s

#Define a variable which contains the run name, it is assumed that it is the parentdir of the last subdir
#In case the input is given in the following form /some/pathway/ skip the last slash  
run_name = path_to_input.split('/')[-1]
if run_name == '':
    run_name = path_to_input.split('/')[-2]

#Create variable with omim reference name 
omim_all = omim_path + 'omim.txt'
#Read the omim file into a dataframe 
omim_df = pd.read_csv(omim_all, header=0,delimiter="\t")

#Create list which will store all subdirs starting with "M" "IK" or "A"
list_of_dirs = []
for dir_element in os.listdir(path_to_input):
    if dir_element.startswith(starting_letter):
        list_of_dirs.append(os.path.join(path_to_input, dir_element))

#Create list which will store the individual dataframes, which will be created for all the individuals in one run
list_of_dfs = []
if reference_type =="normal":
    for path in list_of_dirs:
        current_Mnummer = path.split('/')[-1]
        if NGS_Type == "WES":
            total_file = path + '/results/CNV/' + current_Mnummer + '.cnr'
        else:
            total_file = path + '/exome_extract/CNV/' + current_Mnummer + '.cnr'
        current_total_df = pd.read_csv(total_file, header=0, delimiter="\t")
        list_of_dfs.append(current_total_df)
else:
    for path in list_of_dirs:
        current_Mnummer = path.split('/')[-1]
        if NGS_Type == "WES":
            total_file = path + '/results/CNV/' + current_Mnummer + '_bintest.tsv'
        else:
            total_file = path + '/exome_extract/CNV/' + current_Mnummer + '_bintest.tsv'
        current_total_df = pd.read_csv(total_file, header=0, delimiter="\t")
        list_of_dfs.append(current_total_df)
#Concatenate all the individual dataframes to a reference dataframe
reference = pd.concat(list_of_dfs)
exploded_reference = explode_cnv_table(reference)
ordered_reference = prepare_cnv_table(exploded_reference,omim_df)
#Cast the OMIMG Column to string 
ordered_reference['OMIMG'] = ordered_reference['OMIMG'].astype(str)
ordered_reference['OMIMG'] = ordered_reference['OMIMG'].str.split('.').str[0]
# Define export name
if reference_type == "normal":
    export_name_parquet = path_to_output + 'CNV_reference_' + run_name + '.parquet'
else:
    export_name_parquet = path_to_output + 'CNV_reference_bintest_' + run_name + '.parquet'

# Export reference to parquet file
ordered_reference.to_parquet(export_name_parquet, index=False)
