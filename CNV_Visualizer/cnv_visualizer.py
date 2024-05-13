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

class cnv_visualizer:
    """
    class which is used to format the imported .cnr dataframe, bintest dataframe and the reference dataframe
    """
    def __init__(self,reference_db:pd.DataFrame,cnr_db:pd.DataFrame,bintest_db:pd.DataFrame):
        """
        Constructor of the Class cnv_visualizer
        Input Arguments : 
            reference_db (pandas DataFrame) -> Contains the aggregated information of multiple .cnr files (used for frequency filtering and plots)
            cnr_db (pandas DataFrame) -> Contains the index patients CNV Information, created by importing the .cnr file, created by CNVkit
            bintest_db (pandas Dataframe) -> Contains the index patients bintest CNV Information, created by importing the bintest file, created by CNVkit
        """
        self.reference_db = reference_db
        self.cnr_db = cnr_db
        self.bintest_db = bintest_db

    def explode_df(self,df:pd.DataFrame) -> pd.DataFrame:
        """
        Function which splits the gene column of a pandas Dataframe on a comme and subsequently explodes the column
        Input Arguments:
            df (pandas DataFrame) -> DataFrame to be exploded
        Output Arguments: 
            df (pandas DataFrame) -> exploded DataFrame
        """
        df['gene'] = df['gene'].str.split(',')
        df = df.explode('gene')
        return df 
    
    def prepare_cnv_table(self,df:pd.DataFrame,df2:pd.DataFrame) -> pd.DataFrame :
        """
        Function to process and add relevant Information to the .cnr/bintest DatafFrame 
        1. Drop antitarget Entries 
        2. Apply the reverse function of log2 for easier interpretation
        3. Reorder Columns 
        4. Translate log2-value into call information
        5. merge with omim df
        6. fill None enries 

        Input Arguments:
            df (pandas DataFrame) -> .cnr/bintest DataFrame to be ordered 
            df2 (pandas DataFrame) -> Omim Dataframe created by importing the previously mentioned .txt file 

        Output Arguments:
            df (pandas DataFrame) -> extended and reordered pandas DataFrame
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
        # Apply calling scores, based on log2 values  
        df.loc[df['log2'] <= -1.1, 'call'] = int(0)
        df.loc[(df['log2'] <= -0.4) & (df['log2'] > -1.1), 'call'] = int(1)
        df.loc[(df['log2'] <= 0.3) & (df['log2'] > -0.4), 'call'] = int(2)
        df.loc[df['log2'] > 0.3, 'call'] = int(3)
        # merge with omim df 
        df = pd.merge(df, df2, on='gene', how='left')
        #df.fillna('-', inplace=True)
        df['comments'] = '.'
        return df
    
    def prepare_parent_cnv(self,parent_df:pd.DataFrame) -> pd.DataFrame :
        """
        Slightly altered Version of prepare_cnv_table to process the index patients parental .cnr DataFrames, used for trio-visualization
        Input Arguments : 
            parent_df (Pandas DataFrame) -> Parental .cnr DataFrame
        """
        #Explode parent DataFrame 
        parent_df = self.explode_df(parent_df)
        # Drop Antitarget entries  
        parent_df.drop(parent_df[parent_df['gene'].str.contains('Antitarget') == True].index, inplace = True)
        #Reverse function to log2 
        parent_df.loc[:, 'squaredvalue'] = 2**parent_df['log2']
        #Split gene and exon file into two seperate fields 
        parent_df['gene'] = parent_df['gene'].str.split('_')
        parent_df.loc[:, 'exon'] = parent_df['gene'].str[1].astype(int)
        parent_df.loc[:, 'gene'] = parent_df['gene'].str[0]
        # Reorder the columns to resemble the familiar cnv-table look 
        cols = parent_df.columns.tolist()
        cols = cols[0:4] + cols[-1:] + cols[5:6] + cols[6:7] + cols[4:5] + cols[7:-1]
        parent_df = parent_df[cols]
        # Insert a call column 
        parent_df.insert(loc=7, column='call', value='')
        # Apply calling scores, based on log2 values  
        parent_df.loc[parent_df['log2'] <= -1.1, 'call'] = int(0)
        parent_df.loc[(parent_df['log2'] <= -0.4) & (parent_df['log2'] > -1.1), 'call'] = int(1)
        parent_df.loc[(parent_df['log2'] <= 0.3) & (parent_df['log2'] > -0.4), 'call'] = int(2)
        parent_df.loc[parent_df['log2'] > 0.3, 'call'] = int(3)
        return parent_df
    
    def format_df(self,omim_path:str,selected_candi_path:str):
        """
        Function used to import and subsequently preprocess the omim and candi DataFrame
        Input Arguments : 
            omim_path (str): previously defined path to omim.txt file 
            selected_candi_path (str): previously defined path to selected candigene.txt file
        Output Arguments: 
            omim_df (pandas DataFrame) -> Imported Omim DataFrame
            cand_df (pandas DataFrame) -> Imported Candigene DataFrame
            cnr_db (pandas DataFrame) -> exploded .cnr DataFrame, which is subsequently reformatted using the function prepare_cnv_table
            bintest_db (pandas DataFrame) -> exploded bintest DataFrame
        """
        #Import Omim.txt file to pandas DataFrame
        omim_df = pd.read_csv(omim_path, header=0,delimiter="\t")
        #Import Candigene.txt file to pandas DataFrame
        candi_df = pd.read_csv(selected_candi_path,header=None,names=["gen"], delimiter="\t")
        #Explode gene column of .cnr/bintest DataFrame
        self.cnr_db = self.explode_df(self.cnr_db)
        self.bintest_db = self.explode_df(self.bintest_db)
        #Reformat the .cnr DataFrame using the prepare_cnv_table function
        self.cnr_db = self.prepare_cnv_table(self.cnr_db,omim_df)
        #Perform a groupby function to define how many exons each gene has, save it as a seperat DataFrame
        gene_size = self.cnr_db.groupby("gene")["gene"].size().reset_index(name="gene_size")
        #Merge .cnr Dataf with gene_size DataFrame
        self.cnr_db = pd.merge(self.cnr_db,gene_size,on="gene",how="left")
        #Create bintest table
        self.bintest_db = self.prepare_cnv_table(self.bintest_db,omim_df)
        return omim_df,candi_df,self.cnr_db,self.bintest_db
    
    def filter_for_deletions_hom(self,df:pd.DataFrame) -> pd.DataFrame :
        """
        Function which is used to filter for homozygously deleted exons (preset)
        Input Arguments:
            df (pandas DataFrame) -> .cnr DataFrame
        Ouput Arguments: 
            df (pandas DataFrame) -> .cnr DataFrame filtered for homozygously deleted exons
        """
        df_del_hom = df[df['call'] == 0]
        return df_del_hom
    
    def filter_for_duplications(self,df:pd.DataFrame) -> pd.DataFrame :
        """
        Function which is used to filter for duplicated exons (preset)
        Input Arguments:
            df (pandas DataFrame) -> .cnr DataFrame
        Ouput Arguments: 
            df (pandas DataFrame) -> .cnr DataFrame filtered for duplicated exons
        """
        df_dup = df[df['call'] == 3]
        return df_dup
    
    def filter_for_deletions(self,df:pd.DataFrame) -> pd.DataFrame :
        """
        Function which is used to filter for heterozygously deleted exons (preset)
        Input Arguments:
            df (pandas DataFrame) -> .cnr DataFrame
        Ouput Arguments: 
            df (pandas DataFrame) -> .cnr DataFrame filtered for heterozygously deleted exons
        """
        df_del = df[((df['call'] == 0) | (df['call'] ==1))]
        return df_del

    def prepare_filter_for_consecutive_cnvs(self,df:pd.DataFrame) -> pd.DataFrame :
        """
        Function which is used to filter for consecutively deleted/duplicated exons
        Input Arguments:
            df (pandas DataFrame) -> .cnr DataFrame
        Ouput Arguments: 
            df (pandas DataFrame) -> .cnr DataFrame filtered for consecutively deleted/duplicated exons
        """
        #Group CNVs by gene and calculate the difference between the exons in each group (consecutive exons will have 1/-1)
        # Because exons might be listed as 1,2,3 or 3,2,1 the calculation is performed both ways
        df['difference_previous'] = df.groupby('gene')['exon'].diff()
        df['difference_previous'] = df.groupby('gene')['difference_previous'].fillna(method='backfill')
        df['difference_next'] = df.groupby('gene')['exon'].diff(periods=-1)
        df['difference_next'] = df.groupby('gene')['difference_next'].fillna(method='ffill')
        return df
    
    # Filter for rows that have 1/-1 in one of the two diff rows 
    def filter_for_consecutive_cnvs(self,df:pd.DataFrame,del_or_dup:str,del_size:int,dup_size:int) -> pd.DataFrame :
        """
        Function which extends the prepare_filter_for_consecutive_cnvs function
        It takes the output from the aforementioned function and selects those cnvs, which have a difference to previous/next of 1/-1
        Input Arguments:
            df (pandas DataFrame) -> .cnr DataFrame annotated with consecutive cnvs 
            del_or_up (str) -> String that determines wether to filter for deletions or duplications 
            del_size (int) -> Int that determines the necessary amount of consecutive deletions
            dup_size (int) -> Int that determines the necessary amount of consecutive duplications
        Output Arguments:
            df_cons (pandas DataFrame) -> .cnr DataFrame with applied filter for consecutive deletions/duplications

        """
        #Filter for consecutive deletions
        if del_or_dup == 'del':
            df_cons = df[(((df['call'] == 0) | (df['call'] ==1)) & ((df['difference_previous'] == -1) | (df['difference_previous'] ==1))) | (((df['call'] == 0) | (df['call'] ==1)) & ((df['difference_next'] == -1) | (df['difference_next'] ==1)))]
            #Calculate the size of consecutive deletions per gene 
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
            #Select only those genes which contain the defined amount of deleted exons or which are smaller than the defined number
            df_cons = df_cons[(df_cons['counts']==del_size)|(df_cons['gene_size'] == df_cons['counts'])]
        #Filter for consecutive duplications
        else:
            df_cons = df[(((df['call'] == 3)) & ((df['difference_previous'] == -1) | (df['difference_previous'] ==1))) | (((df['call'] == 3)) & ((df['difference_next'] == -1) | (df['difference_next'] ==1)))]
            #Calculate the size of consecutive duplications per gene 
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
            #Select only those genes which contain the defined amount of duplicated exons or which are smaller than the defined number
            df_cons = df_cons[(df_cons['counts']==dup_size)|(df_cons['gene_size'] == df_cons['counts'])]
        return df_cons
    
    #Filter using candigene list
    def filter_for_candi_cnvs(self,df:pd.DataFrame,df2:pd.DataFrame) -> pd.DataFrame :
        """
        Function which filters for genes contained in the candigene list 
        Input Arguments : 
            df (pandas DataFrame) -> .cnr DataFrame
            df (pandas DataFrame) -> candigene DataFrame 
        """
        #Keep those genes which are contained in the candigene list 
        df['in_candi'] = df.gene.isin(df2.gen)
        #Remove non cnvs
        df_filter_candi = df[(df['in_candi'] == True)& (df['call'] != 2)]
        return df_filter_candi
    
    def apply_filters(self,df:pd.DataFrame, start_selection:int, end_selection:int,depth_selection:float,
                      weight_selection:float,chrom_selection:list,call_selection:list,log2_selection:float,
                      gene_selection:list,chrom_list:list,call_list:list,gene_list:list,
                      het_del_selection:float,hom_del_selection:float,dup_selection:float) -> pd.DataFrame:
        """
        Function which applies the previously defined filters for the .cnr file
        Input Arguments : 
            df (pandas DataFrame) -> .cnr DataFrame
            start_selection (int) -> Filter which defines a starting coordinates (only works if end coordinates are given) and only one chromosome is selected
            end_selection (int) ->  Filter which defines a end coordinates (only works if start coordinates are given) and only one chromosome is selected
            depth_deletion (float) -> Filter which defines a minimal depth 
            weight_selection (float) -> Filter which defines a minimal weight 
            chrom_selection (list) -> Filter which defines which chromosomes shall be displayed
            call_selection (list) -> Filter which defines which calls shall be displayed 
            log2_selection (float) -> Filter which defines a minimal log2 
            gene_selection (list) -> Filter which genes shall be displayed
            chrom_list (list) -> previously defined list with all chromosomes (used to negate an empty filter)
            gene_list (list) -> previously defined list with all genes (used to negate an empty filter)
            het_del_selection (float) -> Filter which defines a maximal heterozygous deletion frequency
            hom_del_selection (float) -> Filter which defines a maximal homozygous deletion frequency
            dup_selection (float) -> Filter which defines a maximal duplication frequency
        """
        #Declare a boolean which is used to determine wether the start/end filter is to be used
        skip_start_end=False
        #Perform Casts to ensure the DataFrame has the correct types
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

        # If only start or only end coordinates are provided or more then one chromosome is defined skip the end/start filter 
        if start_selection is None or start_selection == "" or end_selection is None or end_selection =="" or len(chrom_selection)>1: #or len(chrom_selection)<1:
            skip_start_end=True
        #If no chromosomes are selected negate empty filter by assining predefined chromosome list
        if chrom_selection is None or chrom_selection==[]:
            chrom_selection = chrom_list
        #If no calls are selected negate empty filter by assining predefined chromosome list
        if call_selection is None or call_selection==[]:
            call_selection = call_list
        #If no genes are selected negate empty filter by assining predefined chromosome list
        if gene_selection is None or gene_selection==[]:
            gene_selection = gene_list
        #If no depth is defined, assign a negative value to not lose values <0
        if depth_selection is None or depth_selection=="":
            depth_selection=float(-10000)
        #If no weight is defined, assign a negative value to not lose values <0
        if weight_selection is None or weight_selection=="":
            weight_selection=float(-10000)
        #If no log2 is defined, assign a negative value to not lose values <0
        if log2_selection is None or log2_selection=="":
            log2_selection=float(-10000)
        #If no maximal heterozygous deletion frequency is defined, assign a max frequency of 1 
        if het_del_selection is None or het_del_selection=="":
            het_del_selection = float(1)
         #If no maximal homozygous deletion frequency is defined, assign a max frequency of 1 
        if hom_del_selection is None or hom_del_selection=="":
            hom_del_selection = float(1)
        #If no maximal duplication frequency is defined, assign a max frequency of 1 
        if dup_selection is None or dup_selection=="":
            dup_selection = float(1)
        #This is executed, if the start/end filter is skipped
        if skip_start_end ==True:
            #Filter for chromosomes in chrom_selection
            filtered_df = filtered_df[filtered_df["chromosome"].isin(chrom_selection)]
            #Filter for calls in call_selection
            filtered_df = filtered_df[filtered_df["call"].isin(call_selection)]
            #Filter for genes in gene_selection
            filtered_df = filtered_df[filtered_df["gene"].isin(gene_selection)]
            #Filter for depth greater or equal to minimal depth
            filtered_df = filtered_df[filtered_df["depth"]>=depth_selection]
            #Filter for weight greater or equal to minimal depth
            filtered_df = filtered_df[filtered_df["weight"]>=weight_selection]
            #Filter for log2 greater or equal to minimal log2 
            filtered_df = filtered_df[filtered_df["log2"]>=log2_selection]
            #Filter for heterozygous deletion frequency smaller or equal to max frequency 
            filtered_df = filtered_df[filtered_df["het_del_frequency"]<=het_del_selection]
            #Filter for homozygous deletion frequency smaller or equal to max frequency 
            filtered_df = filtered_df[filtered_df["hom_del_frequency"]<=hom_del_selection]
            #Filter for duplication frequency smaller or equal to max frequency 
            filtered_df = filtered_df[filtered_df["dup_frequency"]<=dup_selection]
        #This is executed if the start/end filter is not skipped
        else:
            #Filter for chromosomes in chrom_selection
            filtered_df = filtered_df[filtered_df["chromosome"].isin(chrom_selection)]
            #Filter for calls in call_selection
            filtered_df = filtered_df[filtered_df["call"].isin(call_selection)]
            #Filter for genes in gene_selection
            filtered_df = filtered_df[filtered_df["gene"].isin(gene_selection)]
            #Filter for depth greater or equal to minimal depth
            filtered_df = filtered_df[filtered_df["depth"]>=depth_selection]
            #Filter for weight greater or equal to minimal depth
            filtered_df = filtered_df[filtered_df["weight"]>=weight_selection]
            #Filter for log2 greater or equal to minimal log2 
            filtered_df = filtered_df[filtered_df["log2"]>=log2_selection]
            #Filter for coordinates greater or equal to the start coordinates
            filtered_df = filtered_df[filtered_df["start"]>=start_selection]
            #Filter for coordinates smaller or equal to the end coordinates
            filtered_df = filtered_df[filtered_df["end"]<=end_selection]
            #Filter for heterozygous deletion frequency smaller or equal to max frequency 
            filtered_df = filtered_df[filtered_df["het_del_frequency"]<=het_del_selection]
            #Filter for homozygous deletion frequency smaller or equal to max frequency 
            filtered_df = filtered_df[filtered_df["hom_del_frequency"]<=hom_del_selection]
            #Filter for duplication frequency smaller or equal to max frequency 
            filtered_df = filtered_df[filtered_df["dup_frequency"]<=dup_selection]
        return filtered_df
    
    def apply_trio_filters(self,trio_df:pd.DataFrame,selection_index:list,selection_father:list,selection_mother:list,call_list:list):
        """
        Function which applies the predefined filters onto the trio .cnr DataFrame
        Input Arguments:
            trio_df (pandas DataFrame) -> DataFrame by merging the index .cnr DataFrame and the DataFrame from the parental .cnr files
            selection_index (list) -> Filter which selects the selected calls for the index patient
            selection_father (list) -> Filter which selects the selected calls for the father of the index patient
            selection_mother (list) -> Filter which selects the selected calls for the mother of the index patient
            call_list (list) -> list which contains all possible calls (used to negate an empty filter)
        """
        #If no calls are selected for the index, negate empty filter by assign value of call_list
        if selection_index is None or selection_index==[]:
            selection_index = call_list
        #If no calls are selected for the father of the index, negate empty filter by assign value of call_list
        if selection_father is None or selection_father==[]:
            selection_father = call_list
        #If no calls are selected for the mother of the index, negate empty filter by assign value of call_list
        if selection_mother is None or selection_mother==[]:
            selection_mother = call_list
        #Apply filter by selecting the defined calls
        filtered_df = trio_df[trio_df["call"].isin(selection_index)]
        filtered_df = filtered_df[filtered_df["call"].isin(selection_father)]
        filtered_df = filtered_df[filtered_df["call"].isin(selection_mother)]
        return filtered_df



    

    