"""
File containing helper functions of the CNVizard.
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import pandas as pd
import streamlit as st


def filter_tsv(tsv: pd.DataFrame, chromosome_list_cnv: list, cnv_type: list, acmg_class: list, 
               entered_cnv_chrom: list, entered_cnv_type: list, entered_acmg_class: list) -> pd.DataFrame:
    """
    Applies filters to the .tsv DataFrame.

    Args:
        tsv (pd.DataFrame): .tsv DataFrame.
        chromosome_list_cnv (list): List of possible chromosomes.
        cnv_type (list): List of possible CNV types.
        acmg_class (list): List of possible ACMG classes.
        entered_cnv_chrom (list): List of selected chromosomes.
        entered_cnv_type (list): List of selected CNV types.
        entered_acmg_class (list): List of selected ACMG classes.

    Returns:
        pd.DataFrame: Filtered .tsv DataFrame.
    """
    tsv["SV_chrom"] = tsv["SV_chrom"].astype(str)
    tsv["SV_type"] = tsv["SV_type"].astype(str)
    tsv["ACMG_class"] = tsv["ACMG_class"].astype(int)

    entered_cnv_chrom = entered_cnv_chrom or chromosome_list_cnv
    entered_cnv_type = entered_cnv_type or cnv_type
    entered_acmg_class = entered_acmg_class or acmg_class

    filtered_tsv = tsv[tsv["SV_chrom"].isin(entered_cnv_chrom)]
    filtered_tsv = filtered_tsv[filtered_tsv["SV_type"].isin(entered_cnv_type)]
    filtered_tsv = filtered_tsv[filtered_tsv["ACMG_class"].isin(entered_acmg_class)]

    return filtered_tsv
