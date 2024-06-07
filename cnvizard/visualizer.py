"""
File which contains the main backbone of the CNVizard, that formats the input files 
@author: Jeremias Krause / Matthias Begemann / Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""
import streamlit as st
import pandas as pd
import os 
import xlsxwriter
import numpy as np
import pyarrow

class CNVVisualizer:
    """
    Class used to format the imported .cnr dataframe, bintest dataframe and the reference dataframe.
    """

    def __init__(self, reference_db: pd.DataFrame, cnr_db: pd.DataFrame, bintest_db: pd.DataFrame):
        """
        Constructor of the Class CNVVisualizer.

        Args:
            reference_db (pd.DataFrame): Contains the aggregated information of multiple .cnr files (used for frequency filtering and plots)
            cnr_db (pd.DataFrame): Contains the index patients CNV Information, created by importing the .cnr file, created by CNVkit
            bintest_db (pd.DataFrame): Contains the index patients bintest CNV Information, created by importing the bintest file, created by CNVkit
        """
        self.reference_db = reference_db
        self.cnr_db = cnr_db
        self.bintest_db = bintest_db

    def explode_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function which splits the gene column of a pandas Dataframe on a comma and subsequently explodes the column.

        Args:
            df (pd.DataFrame): DataFrame to be exploded.

        Returns:
            pd.DataFrame: Exploded DataFrame.
        """
        df['gene'] = df['gene'].str.split(',')
        df = df.explode('gene')
        return df 

    def prepare_cnv_table(self, df: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
        """
        Function to process and add relevant Information to the .cnr/bintest DataFrame.
        1. Drop antitarget Entries.
        2. Apply the reverse function of log2 for easier interpretation.
        3. Reorder Columns.
        4. Translate log2-value into call information.
        5. Merge with OMIM df.
        6. Fill None entries.

        Args:
            df (pd.DataFrame): .cnr/bintest DataFrame to be ordered.
            df2 (pd.DataFrame): Omim Dataframe created by importing the previously mentioned .txt file.

        Returns:
            pd.DataFrame: Extended and reordered pandas DataFrame.
        """
        df.drop(df[df['gene'].str.contains('Antitarget') == True].index, inplace=True)
        df.loc[:, 'CN'] = 2**df['log2']
        df['gene'] = df['gene'].str.split('_')
        df.loc[:, 'exon'] = df['gene'].str[1].astype(int)
        df.loc[:, 'gene'] = df['gene'].str[0]
        cols = df.columns.tolist()
        cols = cols[0:4] + cols[-1:] + cols[5:6] + cols[6:7] + cols[4:5] + cols[7:-1]
        df = df[cols]
        df.insert(loc=7, column='call', value='')
        df.loc[df['log2'] <= -1.1, 'call'] = int(0)
        df.loc[(df['log2'] <= -0.4) & (df['log2'] > -1.1), 'call'] = int(1)
        df.loc[(df['log2'] <= 0.3) & (df['log2'] > -0.4), 'call'] = int(2)
        df.loc[df['log2'] > 0.3, 'call'] = int(3)
        df = pd.merge(df, df2, on='gene', how='left')
        df['comments'] = '.'
        return df
    
    def prepare_parent_cnv(self, parent_df: pd.DataFrame) -> pd.DataFrame:
        """
        Slightly altered Version of prepare_cnv_table to process the index patients parental .cnr DataFrames, used for trio-visualization.

        Args:
            parent_df (pd.DataFrame): Parental .cnr DataFrame.

        Returns:
            pd.DataFrame: Processed parental .cnr DataFrame.
        """
        parent_df = self.explode_df(parent_df)
        parent_df.drop(parent_df[parent_df['gene'].str.contains('Antitarget') == True].index, inplace=True)
        parent_df.loc[:, 'CN'] = 2**parent_df['log2']
        parent_df['gene'] = parent_df['gene'].str.split('_')
        parent_df.loc[:, 'exon'] = parent_df['gene'].str[1].astype(int)
        parent_df.loc[:, 'gene'] = parent_df['gene'].str[0]
        cols = parent_df.columns.tolist()
        cols = cols[0:4] + cols[-1:] + cols[5:6] + cols[6:7] + cols[4:5] + cols[7:-1]
        parent_df = parent_df[cols]
        parent_df.insert(loc=7, column='call', value='')
        parent_df.loc[parent_df['log2'] <= -1.1, 'call'] = int(0)
        parent_df.loc[(parent_df['log2'] <= -0.4) & (parent_df['log2'] > -1.1), 'call'] = int(1)
        parent_df.loc[(parent_df['log2'] <= 0.3) & (parent_df['log2'] > -0.4), 'call'] = int(2)
        parent_df.loc[parent_df['log2'] > 0.3, 'call'] = int(3)
        return parent_df
    
    def format_df(self, omim_path: str, selected_candi_path: str):
        """
        Function used to import and subsequently preprocess the omim and candi DataFrame.

        Args:
            omim_path (str): Path to omim.txt file.
            selected_candi_path (str): Path to selected candigene.txt file.

        Returns:
            tuple: omim_df (pd.DataFrame), cand_df (pd.DataFrame), cnr_db (pd.DataFrame), bintest_db (pd.DataFrame)
        """
        omim_df = pd.read_csv(omim_path, header=0, delimiter="\t")
        candi_df = pd.read_csv(selected_candi_path, header=None, names=["gen"], delimiter="\t")
        self.cnr_db = self.explode_df(self.cnr_db)
        self.bintest_db = self.explode_df(self.bintest_db)
        self.cnr_db = self.prepare_cnv_table(self.cnr_db, omim_df)
        gene_size = self.cnr_db.groupby("gene")["gene"].size().reset_index(name="gene_size")
        self.cnr_db = pd.merge(self.cnr_db, gene_size, on="gene", how="left")
        self.bintest_db = self.prepare_cnv_table(self.bintest_db, omim_df)
        return omim_df, candi_df, self.cnr_db, self.bintest_db
    
    def filter_for_deletions_hom(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function which is used to filter for homozygously deleted exons (preset).

        Args:
            df (pd.DataFrame): .cnr DataFrame

        Returns:
            pd.DataFrame: .cnr DataFrame filtered for homozygously deleted exons.
        """
        df_del_hom = df[df['call'] == 0]
        return df_del_hom
    
    def filter_for_duplications(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function which is used to filter for duplicated exons (preset).

        Args:
            df (pd.DataFrame): .cnr DataFrame

        Returns:
            pd.DataFrame: .cnr DataFrame filtered for duplicated exons.
        """
        df_dup = df[df['call'] == 3]
        return df_dup
    
    def filter_for_deletions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function which is used to filter for heterozygously deleted exons (preset).

        Args:
            df (pd.DataFrame): .cnr DataFrame

        Returns:
            pd.DataFrame: .cnr DataFrame filtered for heterozygously deleted exons.
        """
        df_del = df[((df['call'] == 0) | (df['call'] == 1))]
        return df_del

    def prepare_filter_for_consecutive_cnvs(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function which is used to filter for consecutively deleted/duplicated exons.

        Args:
            df (pd.DataFrame): .cnr DataFrame

        Returns:
            pd.DataFrame: .cnr DataFrame filtered for consecutively deleted/duplicated exons.
        """
        df['difference_previous'] = df.groupby('gene')['exon'].diff()
        df['difference_previous'] = df.groupby('gene')['difference_previous'].fillna(method='backfill')
        df['difference_next'] = df.groupby('gene')['exon'].diff(periods=-1)
        df['difference_next'] = df.groupby('gene')['difference_next'].fillna(method='ffill')
        return df

    def filter_for_consecutive_cnvs(self, df: pd.DataFrame, del_or_dup: str, del_size: int, dup_size: int) -> pd.DataFrame:
        """
        Function which extends the prepare_filter_for_consecutive_cnvs function.
        It takes the output from the aforementioned function and selects those cnvs, which have a difference to previous/next of 1/-1.

        Args:
            df (pd.DataFrame): .cnr DataFrame annotated with consecutive cnvs
            del_or_dup (str): String that determines whether to filter for deletions or duplications
            del_size (int): Number of consecutive deletions
            dup_size (int): Number of consecutive duplications

        Returns:
            pd.DataFrame: .cnr DataFrame with applied filter for consecutive deletions/duplications.
        """
        if del_size is None or del_size == "":
            del_size = 2
        if dup_size is None or dup_size == "":
            dup_size = 2
        del_size = int(del_size)
        dup_size = int(dup_size)

        if del_or_dup == 'del':
            df_cons = df[(((df['call'] == 0) | (df['call'] == 1)) & ((df['difference_previous'] == -1) | (df['difference_previous'] == 1))) |
                          (((df['call'] == 0) | (df['call'] == 1)) & ((df['difference_next'] == -1) | (df['difference_next'] == 1)))]
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
            df_cons = df_cons[(df_cons['counts'] >= del_size) | (df_cons['gene_size'] == df_cons['counts'])]
        else:
            df_cons = df[(((df['call'] == 3)) & ((df['difference_previous'] == -1) | (df['difference_previous'] == 1))) |
                          (((df['call'] == 3)) & ((df['difference_next'] == -1) | (df['difference_next'] == 1)))]
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
            df_cons = df_cons[(df_cons['counts'] >= dup_size) | (df_cons['gene_size'] == df_cons['counts'])]
        return df_cons

    def filter_for_candi_cnvs(self, df: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
        """
        Function which filters for genes contained in the candigene list.

        Args:
            df (pd.DataFrame): .cnr DataFrame
            df2 (pd.DataFrame): candigene DataFrame

        Returns:
            pd.DataFrame: Filtered DataFrame for candigenes.
        """
        df['in_candi'] = df.gene.isin(df2.gen)
        df_filter_candi = df[(df['in_candi'] == True) & (df['call'] != 2)]
        return df_filter_candi

    def apply_filters(self, df: pd.DataFrame, start_selection: int, end_selection: int, depth_selection: float,
                      weight_selection: float, chrom_selection: list, call_selection: list, log2_selection: float,
                      gene_selection: list, chrom_list: list, call_list: list, gene_list: list,
                      het_del_selection: float, hom_del_selection: float, dup_selection: float) -> pd.DataFrame:
        """
        Function which applies the predefined filters for the .cnr file.

        Args:
            df (pd.DataFrame): .cnr DataFrame
            start_selection (int): Filter which defines a starting coordinates (only works if end coordinates are given) and only one chromosome is selected
            end_selection (int): Filter which defines a end coordinates (only works if start coordinates are given) and only one chromosome is selected
            depth_selection (float): Filter which defines a minimal depth 
            weight_selection (float): Filter which defines a minimal weight 
            chrom_selection (list): Filter which defines which chromosomes shall be displayed
            call_selection (list): Filter which defines which calls shall be displayed 
            log2_selection (float): Filter which defines a minimal log2 
            gene_selection (list): Filter which genes shall be displayed
            chrom_list (list): Previously defined list with all chromosomes (used to negate an empty filter)
            gene_list (list): Previously defined list with all genes (used to negate an empty filter)
            het_del_selection (float): Filter which defines a maximal heterozygous deletion frequency
            hom_del_selection (float): Filter which defines a maximal homozygous deletion frequency
            dup_selection (float): Filter which defines a maximal duplication frequency

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        skip_start_end = False
        filtered_df = df
        filtered_df["chromosome"] = filtered_df["chromosome"].astype(str)
        filtered_df["call"] = filtered_df["call"].astype(int)
        filtered_df["gene"] = filtered_df["gene"].astype(str)
        filtered_df["depth"] = filtered_df["depth"].astype(float)
        filtered_df["weight"] = filtered_df["weight"].astype(float)
        filtered_df["log2"] = filtered_df["log2"].astype(float)
        filtered_df["start"] = filtered_df["start"].astype(int)
        filtered_df["end"] = filtered_df["end"].astype(int)
        filtered_df["het_del_frequency"] = filtered_df["het_del_frequency"].astype(float)
        filtered_df["hom_del_frequency"] = filtered_df["hom_del_frequency"].astype(float)
        filtered_df["dup_frequency"] = filtered_df["dup_frequency"].astype(float)

        if start_selection is None or start_selection == "" or end_selection is None or end_selection == "" or len(chrom_selection) > 1:
            skip_start_end = True

        if chrom_selection is None or chrom_selection == []:
            chrom_selection = chrom_list
        if call_selection is None or call_selection == []:
            call_selection = call_list
        if gene_selection is None or gene_selection == []:
            gene_selection = gene_list
        if depth_selection is None or depth_selection == "":
            depth_selection = float(-10000)
        if weight_selection is None or weight_selection == "":
            weight_selection = float(-10000)
        if log2_selection is None or log2_selection == "":
            log2_selection = float(-10000)
        if het_del_selection is None or het_del_selection == "":
            het_del_selection = float(1)
        if hom_del_selection is None or hom_del_selection == "":
            hom_del_selection = float(1)
        if dup_selection is None or dup_selection == "":
            dup_selection = float(1)

        if skip_start_end:
            filtered_df = filtered_df[filtered_df["chromosome"].isin(chrom_selection)]
            filtered_df = filtered_df[filtered_df["call"].isin(call_selection)]
            filtered_df = filtered_df[filtered_df["gene"].isin(gene_selection)]
            filtered_df = filtered_df[filtered_df["depth"] >= depth_selection]
            filtered_df = filtered_df[filtered_df["weight"] >= weight_selection]
            filtered_df = filtered_df[filtered_df["log2"] >= log2_selection]
            filtered_df = filtered_df[filtered_df["het_del_frequency"] <= het_del_selection]
            filtered_df = filtered_df[filtered_df["hom_del_frequency"] <= hom_del_selection]
            filtered_df = filtered_df[filtered_df["dup_frequency"] <= dup_selection]
        else:
            filtered_df = filtered_df[filtered_df["chromosome"].isin(chrom_selection)]
            filtered_df = filtered_df[filtered_df["call"].isin(call_selection)]
            filtered_df = filtered_df[filtered_df["gene"].isin(gene_selection)]
            filtered_df = filtered_df[filtered_df["depth"] >= depth_selection]
            filtered_df = filtered_df[filtered_df["weight"] >= weight_selection]
            filtered_df = filtered_df[filtered_df["log2"] >= log2_selection]
            filtered_df = filtered_df[filtered_df["start"] >= start_selection]
            filtered_df = filtered_df[filtered_df["end"] <= end_selection]
            filtered_df = filtered_df[filtered_df["het_del_frequency"] <= het_del_selection]
            filtered_df = filtered_df[filtered_df["hom_del_frequency"] <= hom_del_selection]
            filtered_df = filtered_df[filtered_df["dup_frequency"] <= dup_selection]

        return filtered_df

    def apply_trio_filters(self, trio_df: pd.DataFrame, selection_index: list, selection_father: list, selection_mother: list, call_list: list) -> pd.DataFrame:
        """
        Function which applies the predefined filters onto the trio .cnr DataFrame.

        Args:
            trio_df (pd.DataFrame): DataFrame by merging the index .cnr DataFrame and the DataFrame from the parental .cnr files
            selection_index (list): Filter which selects the selected calls for the index patient
            selection_father (list): Filter which selects the selected calls for the father of the index patient
            selection_mother (list): Filter which selects the selected calls for the mother of the index patient
            call_list (list): List which contains all possible calls (used to negate an empty filter)

        Returns:
            pd.DataFrame: Filtered trio DataFrame.
        """
        if selection_index is None or selection_index == []:
            selection_index = call_list
        if selection_father is None or selection_father == []:
            selection_father = call_list
        if selection_mother is None or selection_mother == []:
            selection_mother = call_list

        filtered_df = trio_df[trio_df["call"].isin(selection_index)]
        filtered_df = filtered_df[filtered_df["call"].isin(selection_father)]
        filtered_df = filtered_df[filtered_df["call"].isin(selection_mother)]
        return filtered_df
