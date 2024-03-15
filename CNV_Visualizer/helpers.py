"""
File which contains little helpers of the CNVizard
@author: Jeremias Krause / Matthias Begemann / Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd 
import streamlit as st

def filter_tsv(tsv:pd.DataFrame,chromosome_list_cnv:list,cnv_type:list,acmg_class:list,entered_cnv_chrom:list,entered_cnv_type:list,entered_acmg_class:list):
    """
    Function which applies previously defined filters on the .tsv DataFrame 
    Note : Defined outside of the cnv_visualizer so the user can utilise the function without initialising an instance of the cnv_visualizer, which 
    requires the reference df, a .cnr df and the bintest df
    Input Arguments:
        tsv (pandas DataFrame) -> .tsv DataFrame
        chromosome_list_cnv (list) -> list which contains all the possible chromosomes (used to negate an empty filter)
        cnv_type (list) -> list which contains all the possible cnv_types (used to negate an empty filter)
        acmg_class (list) -> list which contains all the possible acmg_classes (used to negate an empty filter)
        entered_cnv_chrom (list) -> list which contains the chromosomes selected by the user 
        entered_cnv_type (list) -> list which contains the cnv types selected by the user
        entered_acmg_class (list) -> list which contains the acmg classes selected by the user
    Output Arguments:
        filtered_tsv (pandas DataFrame) -> filtered .tsv DataFrame
    """
    #Perform casts to ensure the dataframe columns contain the correct types
    tsv["SV_chrom"] = tsv["SV_chrom"].astype(str)
    tsv["SV_type"] = tsv["SV_type"].astype(str)
    tsv["ACMG_class"] = tsv["ACMG_class"].astype(int)
    #If entered_cnv_chrom is empty the empty filter is negated by assining chromosome_list_cnv
    if entered_cnv_chrom is None or entered_cnv_chrom==[]:
        entered_cnv_chrom=chromosome_list_cnv
    #If entered_cnv_type is empty the empty filter is negated by assining cnv_type
    if entered_cnv_type is None or entered_cnv_type==[]:
        entered_cnv_type=cnv_type
    #If entered_acmg_class is empty the empty filter is negated by assining acmg_class   
    if entered_acmg_class is None or entered_acmg_class==[]:
        entered_acmg_class=acmg_class
    #Apply filters
    filtered_tsv = tsv[tsv["SV_chrom"].isin(entered_cnv_chrom)]
    filtered_tsv = filtered_tsv[filtered_tsv["SV_type"].isin(entered_cnv_type)]
    filtered_tsv = filtered_tsv[filtered_tsv["ACMG_class"].isin(entered_acmg_class)]
    return filtered_tsv
