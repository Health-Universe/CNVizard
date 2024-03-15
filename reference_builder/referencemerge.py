"""
File which is used to merge the previously created individual reference files used for plotting with the CNVizard
@author: Jeremias Krause / Matthias Begemann / Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""
import pandas as pd 
import pyarrow
import fastparquet
import numpy
import plotly
import plotly.graph_objects as go
import argparse
import os 

def _get_call_counts(call: list):
    """
    Function which will be used as a lambda function to calculate the amount of dels,dups,wild types from a call list
    Input Arguments : 
        call (list) : list inside a Pandas DataFrame column which stores all calls (per exon) inside the reference DataFrame
    Output Arguments : 
        call_list (list): list inside a pandas DataFrame column which stores the sums of all call-types (per exon) inside the reference DataFrame
    """
    #Define an empty in which the amount of heterozygous deletions, homozygous deletions, duplications, wildtypes and the total amount of calls will be stored
    call_list = []
    #Set all call types to zero
    dels_het = 0
    dels_hom = 0
    dups = 0 
    wt = 0
    all_calls = 0
    #Iterate over the call list to calculate the amount of calls per cnv type
    for x in call:
        if x == 1:
            dels_het = dels_het + 1 
        elif x == 0:
            dels_hom = dels_hom + 1
        elif x == 2:
            wt = wt + 1 
        else:
            dups = dups +1
    #Sum all calls to calculate the total amount of calls
    all_calls = all_calls + dels_het + dels_hom + wt + dups
    #Append calculated sums for call_types
    call_list.append(dels_het)
    call_list.append(dels_hom)
    call_list.append(dups)
    call_list.append(wt)
    call_list.append(all_calls)
    return call_list

def _get_frequency(call_counts: list, position1: int, position2: int):
    """
    Function to calculate the frequency of different call types inside the reference DataFrame
    This function will be applied as a lambda function onto a pandas DataFrame
    Input Arguments : 
        call_counts (list) : list which contains the previously calculated sums for each call-types
        position1 (int) : This position selects the call-type for which the frequency is to be calculated 
        position2 (int) : This position selects the total sum of individuals inside the reference DataFrame
    Output Arguments : 
        frequency (float) : frequency of a cnv-type per exon inside the reference DataFrame
    """
    #Get number of cnvs with specific cnv/type (0=heterozygous deletions, 1=homozygous deletions, 2=duplications,3=wildtype,4=total_sum)
    number_of_cnvs = call_counts[position1]
    #Get total sum
    number_of_individuals = call_counts[position2]
    #Calculate frequency :
    frequency = number_of_cnvs / number_of_individuals
    return frequency 

def _get_frequency_bintest(call_counts: list, position1: int, position2: int, call_counts_ref:list):
    """
    Function to calculate the frequency of different call types inside the reference DataFrame
    This function will be applied as a lambda function onto a pandas DataFrame
    Input Arguments : 
        call_counts (list) : list which contains the previously calculated sums for each call-types
        position1 (int) : This position selects the call-type for which the frequency is to be calculated 
        position2 (int) : This position selects the total sum of individuals inside the reference DataFrame
        call_counts_ref (list) : list which contains the previously calculated sums for each call-types (from the "normal" reference)
    Output Arguments : 
        frequency (float) : frequency of a cnv-type per exon inside the reference DataFrame
    """
    #Get number of cnvs with specific cnv/type (0=heterozygous deletions, 1=homozygous deletions, 2=duplications,3=wildtype,4=total_sum)
    number_of_cnvs = call_counts[position1]
    #Get total sum
    number_of_individuals = call_counts_ref[position2]
    #Calculate frequency :
    frequency = number_of_cnvs / number_of_individuals
    return frequency 

#Create Commandline parser
parser = argparse.ArgumentParser(description="Create Reference Parquet")
# Required positional argument -i for input path and -o for output path
parser.add_argument('-i', type=str,
                    help='Positional argument for input path (Reference .parquet files)')
parser.add_argument('-o', type=str,
                    help='Positional argument for output')
parser.add_argument('-b', type=str,
                    help='Positional argument for bintest input path (Reference .parquet files)')

#Get Commandline arguments by using parser
args = parser.parse_args()

#Load parsed arguments into corresponding variables
path_to_input = args.i
path_to_output = args.o
path_to_bintest = args.b

#List all reference parquet files 
genome_reference_list_of_files = os.listdir(path_to_input)
#Create an empty list which will store all individual reference parquet files
genome_reference_list_of_parquets=[]

#Iterate over genome_reference_list_of_files, store the individual reference inside a pandas DataFrame and append it to genome_reference_list_of_parquets
for file in genome_reference_list_of_files:
    table = pd.read_parquet(path_to_input+file)
    genome_reference_list_of_parquets.append(table)

#Concatenate all individual reference DataFrames 
genome_reference = pd.concat(genome_reference_list_of_parquets, ignore_index=True)
#Drop unnecessary columns to preserve memory
genome_reference = genome_reference.drop(["chromosome","start","end","weight","squaredvalue","OMIMG",
                                         "Disease","OMIMP","Inheritance","comments"],axis=1)

#Group by gene and exon, to aggregate the columns depth,log2 and call from all individuals into a single column
genome_reference = genome_reference.groupby(["gene","exon"]).agg({"depth": list,
                                                                 "log2": list,
                                                                 "call": list}).reset_index()

#Apply the previously defined function _get_call_counts
genome_reference["call_counts"] = genome_reference.apply(lambda x: _get_call_counts(x["call"],),axis=1,)
#Apply the previously defined function _get_frequency to get the frequency of heterozygous deletions
genome_reference["het_del_frequency"] = genome_reference.apply(lambda x: _get_frequency(x["call_counts"], 0,4,),axis=1,)
#Apply the previously defined function _get_frequency to get the frequency of homozygous deletions
genome_reference["hom_del_frequency"] = genome_reference.apply(lambda x: _get_frequency(x["call_counts"], 1,4,),axis=1,)
#Apply the previously defined function _get_frequency to get the frequency of duplications
genome_reference["dup_frequency"] = genome_reference.apply(lambda x: _get_frequency(x["call_counts"], 2,4,),axis=1,)
#Get maximal log2 value per row
genome_reference["max_log2"] = genome_reference.apply(lambda x: max(x["log2"],),axis=1,)
#Get minimal log2 value per row
genome_reference["min_log2"] = genome_reference.apply(lambda x: min(x["log2"],),axis=1,)
#Get mean of depth per row
genome_reference["mean_depth"] = genome_reference.apply(lambda x: sum(x["depth"])/len(x["depth"]),axis=1,)
#Get mean of log2 per row
genome_reference["mean_log2"] = genome_reference.apply(lambda x: sum(x["log2"])/len(x["log2"]),axis=1,)
#Get median of depth per row
genome_reference["median_depth"] = genome_reference.apply(lambda x: numpy.quantile(x["depth"],[0.5])[0],axis=1,)
#Get median of log2 per row
genome_reference["median_log2"] = genome_reference.apply(lambda x: numpy.quantile(x["log2"],[0.5])[0],axis=1,)
#Get q1 of depth per row
genome_reference["q1_depth"] = genome_reference.apply(lambda x: numpy.quantile(x["depth"],[0.25])[0],axis=1,)
#Get q1 of log2 per row
genome_reference["q1_log2"] = genome_reference.apply(lambda x: numpy.quantile(x["log2"],[0.25])[0],axis=1,)
#Get q3 of depth per row
genome_reference["q3_depth"] = genome_reference.apply(lambda x: numpy.quantile(x["depth"],[0.75])[0],axis=1,)
#Get q3 of log2 per row
genome_reference["q3_log2"] = genome_reference.apply(lambda x: numpy.quantile(x["log2"],[0.75])[0],axis=1,)
#Get std of depth per row
genome_reference["std_depth"] = genome_reference.apply(lambda x: numpy.std(x["depth"]),axis=1,)
#Get log2 of depth per row
genome_reference["std_log2"] = genome_reference.apply(lambda x: numpy.std(x["log2"]),axis=1,)

#Calculate the box_size of the log2 boxplot by calculating the distance between the q1 and q3 quartile
genome_reference["box_size"] = (genome_reference["q3_log2"] - genome_reference["q1_log2"])*1.5
#Limit the acutual minimum : the smallest value is q1_log2 - distance_between_q1_q3_quartile_times_1.5x
genome_reference["actual_minimum_log2"] = genome_reference["q1_log2"] - genome_reference["box_size"]
#Limit the acutual maximum : the smallest value is q1_log2 - distance_between_q1_q3_quartile_times_1.5x
genome_reference["actual_maximum_log2"] = genome_reference["q3_log2"] + genome_reference["box_size"]
#Calculate the box_size of the depth boxplot by calculating the distance between the q1 and q3 quartile
genome_reference["box_size_depth"] = (genome_reference["q3_depth"] - genome_reference["q1_depth"])*1.5
#Limit the acutual minimum : the smallest value is q1_depth - distance_between_q1_q3_quartile_times_1.5x
genome_reference["actual_minimum_depth"] = genome_reference["q1_depth"] - genome_reference["box_size_depth"]
#Limit the acutual maximum : the smallest value is q1_depth - distance_between_q1_q3_quartile_times_1.5x
genome_reference["actual_maximum_depth"] = genome_reference["q3_depth"] + genome_reference["box_size_depth"]

#create a new pandas dataFrame which contains gene, exon and call_counts from the "normal" references
ref_counts_dataframe = pd.DataFrame()
ref_counts_dataframe["gene"] = genome_reference["gene"]
ref_counts_dataframe["exon"] = genome_reference["exon"]
ref_counts_dataframe["call_counts"] = genome_reference["call_counts"]

#List all reference parquet files 
genome_reference_list_of_files_bintest = os.listdir(path_to_bintest)
#Create an empty lift which will store all individual reference parquet files 
genome_reference_list_of_parquets_bintest=[]

#Iterate over genome_reference_list_of_files, store the individual reference inside a pandas DataFrame and append it 
for file in genome_reference_list_of_files_bintest:
    table = pd.read_parquet(path_to_bintest+file)
    genome_reference_list_of_parquets_bintest.append(table)

#Concatenate all individual reference DataFrames 
genome_reference_bintest = pd.concat(genome_reference_list_of_parquets_bintest, ignore_index=True)
#Drop unnecessary columns to preserve memory
genome_reference_bintest = genome_reference_bintest.drop(["chromosome","start","end","weight","squaredvalue","OMIMG",
                                         "Disease","OMIMP","Inheritance","comments"],axis=1)

#Group by gene and exon, to aggregate the columns depth,log2 and call from all individuals into a single column
genome_reference_bintest = genome_reference_bintest.groupby(["gene","exon"]).agg({"depth": list,
                                                                 "log2": list,
                                                                 "call": list}).reset_index()

#Apply the previously defined function _get_call_counts 
genome_reference_bintest["call_counts"] = genome_reference_bintest.apply(lambda x: _get_call_counts(x["call"],),axis=1,)
#Merge with bintest reference dataframe 
genome_reference_bintest = pd.merge(genome_reference_bintest,ref_counts_dataframe,how="left",left_on=["gene", "exon"],right_on=["gene", "exon"])

#Apply the previously defined function _get_frequency to get the frequency of heterozygous deletions (bintest)
genome_reference_bintest["het_del_frequency"] = genome_reference_bintest.apply(lambda x: _get_frequency_bintest(x["call_counts_x"], 0,4,x["call_counts_y"],),axis=1,)
#Apply the previously defined function _get_frequency to get the frequency of homozygous deletions (bintest)
genome_reference_bintest["hom_del_frequency"] = genome_reference_bintest.apply(lambda x: _get_frequency_bintest(x["call_counts_x"], 1,4,x["call_counts_y"],),axis=1,)
#Apply the previously defined function _get_frequency to get the frequency of duplications (bintest)
genome_reference_bintest["dup_frequency"] = genome_reference_bintest.apply(lambda x: _get_frequency_bintest(x["call_counts_x"], 2,4,x["call_counts_y"],),axis=1,)

#drop unnecessary columns
genome_reference_bintest = genome_reference_bintest.drop(["depth", "log2", "call", "call_counts_x", "call_counts_y"], axis=1)
#Write genome reference for parquet 
genome_reference_bintest.to_parquet(path_to_output+"cnv_reference_bintest.parquet")

#Once all the required statistics for the boxplot have been calculated, drop the unnecessary columns
genome_reference = genome_reference.drop(["depth","log2","call","call_counts","max_log2","min_log2","box_size","box_size_depth"],axis=1)
#Export merged and formatted reference to parquet file
genome_reference.to_parquet(path_to_output + "cnv_reference.parquet")
