"""
Utility functions for data processing in CNVizard
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import pandas as pd

def prepare_cnv_table(df: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """
    Format and reorder the .cnr DataFrame.

    Parameters:
        df (pd.DataFrame): .cnr DataFrame 
        df2 (pd.DataFrame): OMIM DataFrame 

    Returns:
        pd.DataFrame: Reordered and formatted DataFrame
    """
    df = df[~df['gene'].str.contains('Antitarget')]
    df['squaredvalue'] = 2 ** df['log2']
    df[['gene', 'exon']] = df['gene'].str.split('_', expand=True)
    df['exon'] = df['exon'].astype(int)
    cols = df.columns.tolist()
    cols = cols[:4] + cols[-1:] + cols[5:6] + cols[6:7] + cols[4:5] + cols[7:-1]
    df = df[cols]
    df.insert(loc=7, column='call', value='')
    df.loc[df['log2'] <= -1.1, 'call'] = 0
    df.loc[(df['log2'] <= -0.4) & (df['log2'] > -1.1), 'call'] = 1
    df.loc[(df['log2'] <= 0.3) & (df['log2'] > -0.4), 'call'] = 2
    df.loc[df['log2'] > 0.3, 'call'] = 3
    df = pd.merge(df, df2, on='gene', how='left')
    df.fillna('-', inplace=True)
    df['comments'] = '.'
    return df

def explode_cnv_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Split redundant exons and explode the DataFrame.

    Parameters:
        df (pd.DataFrame): .cnr DataFrame

    Returns:
        pd.DataFrame: Formatted .cnr DataFrame
    """
    df['gene'] = df['gene'].str.split(',')
    df = df.explode('gene')
    return df

def _get_call_counts(call: list) -> list:
    """
    Calculate the amount of deletions, duplications, wild types from a call list.

    Parameters:
        call (list): List inside a pandas DataFrame column which stores all calls (per exon) inside the reference DataFrame.

    Returns:
        list: List inside a pandas DataFrame column which stores the sums of all call-types (per exon) inside the reference DataFrame.
    """
    dels_het, dels_hom, dups, wt = 0, 0, 0, 0
    for x in call:
        if x == 1:
            dels_het += 1
        elif x == 0:
            dels_hom += 1
        elif x == 2:
            wt += 1
        else:
            dups += 1
    all_calls = dels_het + dels_hom + wt + dups
    return [dels_het, dels_hom, dups, wt, all_calls]

def _get_frequency(call_counts: list, position1: int, position2: int) -> float:
    """
    Calculate the frequency of different call types inside the reference DataFrame.

    Parameters:
        call_counts (list): List which contains the previously calculated sums for each call-types.
        position1 (int): This position selects the call-type for which the frequency is to be calculated.
        position2 (int): This position selects the total sum of individuals inside the reference DataFrame.

    Returns:
        float: Frequency of a CNV-type per exon inside the reference DataFrame.
    """
    number_of_cnvs = call_counts[position1]
    number_of_individuals = call_counts[position2]
    frequency = number_of_cnvs / number_of_individuals
    return frequency

def _get_frequency_bintest(call_counts: list, position1: int, position2: int, call_counts_ref: list) -> float:
    """
    Calculate the frequency of different call types inside the reference DataFrame for bintest.

    Parameters:
        call_counts (list): List which contains the previously calculated sums for each call-types.
        position1 (int): This position selects the call-type for which the frequency is to be calculated.
        position2 (int): This position selects the total sum of individuals inside the reference DataFrame.
        call_counts_ref (list): List which contains the previously calculated sums for each call-types (from the "normal" reference).

    Returns:
        float: Frequency of a CNV-type per exon inside the reference DataFrame.
    """
    number_of_cnvs = call_counts[position1]
    number_of_individuals = call_counts_ref[position2]
    frequency = number_of_cnvs / number_of_individuals
    return frequency
