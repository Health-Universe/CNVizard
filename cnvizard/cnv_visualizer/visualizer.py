"""
Main backbone of the CNVizard that formats the input files.
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import pandas as pd

class CNVVisualizer:
    """
    Class used to format the imported .cnr dataframe, bintest dataframe, and the reference dataframe.
    """

    def __init__(self, reference_db: pd.DataFrame, cnr_db: pd.DataFrame, bintest_db: pd.DataFrame):
        """
        Constructor for the CNVVisualizer class.

        Args:
            reference_db (pd.DataFrame): Contains aggregated information of multiple .cnr files (used for frequency filtering and plots).
            cnr_db (pd.DataFrame): Contains the index patient's CNV information from the .cnr file created by CNVkit.
            bintest_db (pd.DataFrame): Contains the index patient's bintest CNV information from the bintest file created by CNVkit.
        """
        self.reference_db = reference_db
        self.cnr_db = cnr_db
        self.bintest_db = bintest_db

    def explode_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Splits the gene column of a DataFrame on commas and subsequently explodes the column.

        Args:
            df (pd.DataFrame): DataFrame to be exploded.

        Returns:
            pd.DataFrame: Exploded DataFrame.
        """
        df['gene'] = df['gene'].str.split(',')
        return df.explode('gene')

    def prepare_cnv_table(self, df: pd.DataFrame, omim_df: pd.DataFrame) -> pd.DataFrame:
        """
        Processes and adds relevant information to the .cnr/bintest DataFrame.

        Args:
            df (pd.DataFrame): .cnr/bintest DataFrame to be ordered.
            omim_df (pd.DataFrame): Omim DataFrame created from the omim.txt file.

        Returns:
            pd.DataFrame: Extended and reordered DataFrame.
        """
        df = df[~df['gene'].str.contains('Antitarget')]
        df['CN'] = 2 ** df['log2']
        df[['gene', 'exon']] = df['gene'].str.split('_', expand=True)
        df['exon'] = df['exon'].astype(int)
        cols = df.columns.tolist()
        df = df[cols[:4] + cols[-1:] + cols[5:7] + cols[4:5] + cols[7:-1]]
        df.insert(loc=7, column='call', value='')
        df.loc[df['log2'] <= -1.1, 'call'] = 0
        df.loc[(df['log2'] <= -0.4) & (df['log2'] > -1.1), 'call'] = 1
        df.loc[(df['log2'] <= 0.3) & (df['log2'] > -0.4), 'call'] = 2
        df.loc[df['log2'] > 0.3, 'call'] = 3
        df = pd.merge(df, omim_df, on='gene', how='left')
        df['comments'] = '.'
        return df

    def prepare_parent_cnv(self, parent_df: pd.DataFrame) -> pd.DataFrame:
        """
        Processes the index patient's parental .cnr DataFrames for trio visualization.

        Args:
            parent_df (pd.DataFrame): Parental .cnr DataFrame.

        Returns:
            pd.DataFrame: Processed parental .cnr DataFrame.
        """
        parent_df = self.explode_df(parent_df)
        parent_df = parent_df[~parent_df['gene'].str.contains('Antitarget')]
        parent_df['CN'] = 2 ** parent_df['log2']
        parent_df[['gene', 'exon']] = parent_df['gene'].str.split('_', expand=True)
        parent_df['exon'] = parent_df['exon'].astype(int)
        cols = parent_df.columns.tolist()
        parent_df = parent_df[cols[:4] + cols[-1:] + cols[5:7] + cols[4:5] + cols[7:-1]]
        parent_df.insert(loc=7, column='call', value='')
        parent_df.loc[parent_df['log2'] <= -1.1, 'call'] = 0
        parent_df.loc[(parent_df['log2'] <= -0.4) & (parent_df['log2'] > -1.1), 'call'] = 1
        parent_df.loc[(parent_df['log2'] <= 0.3) & (parent_df['log2'] > -0.4), 'call'] = 2
        parent_df.loc[parent_df['log2'] > 0.3, 'call'] = 3
        return parent_df

    def format_df(self, omim_path: str, selected_candi_path: str):
        """
        Imports and preprocesses the omim and candi DataFrame.

        Args:
            omim_path (str): Path to omim.txt file.
            selected_candi_path (str): Path to selected candigene.txt file.

        Returns:
            tuple: Imported and formatted DataFrames.
        """
        omim_df = pd.read_csv(omim_path, delimiter="\t")
        candi_df = pd.read_csv(selected_candi_path, header=None, names=["gene"], delimiter="\t")
        self.cnr_db = self.explode_df(self.cnr_db)
        self.bintest_db = self.explode_df(self.bintest_db)
        self.cnr_db = self.prepare_cnv_table(self.cnr_db, omim_df)
        gene_size = self.cnr_db.groupby("gene").size().reset_index(name="gene_size")
        self.cnr_db = pd.merge(self.cnr_db, gene_size, on="gene", how="left")
        self.bintest_db = self.prepare_cnv_table(self.bintest_db, omim_df)
        return omim_df, candi_df, self.cnr_db, self.bintest_db

    def filter_for_deletions_hom(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filters for homozygously deleted exons.

        Args:
            df (pd.DataFrame): .cnr DataFrame.

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        return df[df['call'] == 0]

    def filter_for_duplications(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filters for duplicated exons.

        Args:
            df (pd.DataFrame): .cnr DataFrame.

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        return df[df['call'] == 3]

    def filter_for_deletions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filters for heterozygously deleted exons.

        Args:
            df (pd.DataFrame): .cnr DataFrame.

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        return df[df['call'].isin([0, 1])]

    def prepare_filter_for_consecutive_cnvs(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepares DataFrame for filtering consecutive CNVs.

        Args:
            df (pd.DataFrame): .cnr DataFrame.

        Returns:
            pd.DataFrame: Annotated DataFrame for consecutive CNVs.
        """
        df['difference_previous'] = df.groupby('gene')['exon'].diff().fillna(method='backfill')
        df['difference_next'] = df.groupby('gene')['exon'].diff(periods=-1).fillna(method='ffill')
        return df

    def filter_for_consecutive_cnvs(self, df: pd.DataFrame, del_or_dup: str, del_size: int, dup_size: int) -> pd.DataFrame:
        """
        Filters for consecutive CNVs based on the given parameters.

        Args:
            df (pd.DataFrame): Annotated DataFrame.
            del_or_dup (str): String to determine whether to filter for deletions or duplications.
            del_size (int): Minimum number of consecutive deletions.
            dup_size (int): Minimum number of consecutive duplications.

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        del_size = int(del_size or 2)
        dup_size = int(dup_size or 2)

        if del_or_dup == 'del':
            filtered_df = df[(df['call'].isin([0, 1])) & (df['difference_previous'].isin([-1, 1])) | 
                              ((df['call'].isin([0, 1])) & (df['difference_next'].isin([-1, 1])))]
            affected_size = filtered_df.groupby('gene').size().reset_index(name='counts')
            filtered_df = pd.merge(filtered_df, affected_size, on='gene', how='left')
            filtered_df = filtered_df[(filtered_df['counts'] >= del_size) | (filtered_df['gene_size'] == filtered_df['counts'])]
        else:
            filtered_df = df[(df['call'] == 3) & (df['difference_previous'].isin([-1, 1])) | 
                              ((df['call'] == 3) & (df['difference_next'].isin([-1, 1])))]
            affected_size = filtered_df.groupby('gene').size().reset_index(name='counts')
            filtered_df = pd.merge(filtered_df, affected_size, on='gene', how='left')
            filtered_df = filtered_df[(filtered_df['counts'] >= dup_size) | (filtered_df['gene_size'] == filtered_df['counts'])]

        return filtered_df

    def filter_for_candi_cnvs(self, df: pd.DataFrame, candi_df: pd.DataFrame) -> pd.DataFrame:
        """
        Filters for genes contained in the candigene list.

        Args:
            df (pd.DataFrame): .cnr DataFrame.
            candi_df (pd.DataFrame): Candigene DataFrame.

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        df['in_candi'] = df['gene'].isin(candi_df['gene'])
        return df[(df['in_candi']) & (df['call'] != 2)]

    def apply_filters(self, df: pd.DataFrame, start_selection: int, end_selection: int, depth_selection: float,
                      weight_selection: float, chrom_selection: list, call_selection: list, log2_selection: float,
                      gene_selection: list, chrom_list: list, call_list: list, gene_list: list,
                      het_del_selection: float, hom_del_selection: float, dup_selection: float) -> pd.DataFrame:
        """
        Applies predefined filters to the .cnr DataFrame.

        Args:
            df (pd.DataFrame): .cnr DataFrame.
            start_selection (int): Start coordinate for filtering.
            end_selection (int): End coordinate for filtering.
            depth_selection (float): Minimum depth for filtering.
            weight_selection (float): Minimum weight for filtering.
            chrom_selection (list): List of selected chromosomes.
            call_selection (list): List of selected calls.
            log2_selection (float): Minimum log2 value for filtering.
            gene_selection (list): List of selected genes.
            chrom_list (list): List of all chromosomes.
            call_list (list): List of all possible calls.
            gene_list (list): List of all genes.
            het_del_selection (float): Maximum heterozygous deletion frequency.
            hom_del_selection (float): Maximum homozygous deletion frequency.
            dup_selection (float): Maximum duplication frequency.

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        skip_start_end = not all([start_selection, end_selection]) or len(chrom_selection) != 1
        chrom_selection = chrom_selection or chrom_list
        call_selection = call_selection or call_list
        gene_selection = gene_selection or gene_list
        depth_selection = depth_selection or float(-10000)
        weight_selection = weight_selection or float(-10000)
        log2_selection = log2_selection or float(-10000)
        het_del_selection = het_del_selection or float(1)
        hom_del_selection = hom_del_selection or float(1)
        dup_selection = dup_selection or float(1)

        filtered_df = df.copy()
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

    def apply_trio_filters(self, trio_df: pd.DataFrame, selection_index: list, selection_father: list, 
                           selection_mother: list, call_list: list) -> pd.DataFrame:
        """
        Applies predefined filters to the trio .cnr DataFrame.

        Args:
            trio_df (pd.DataFrame): DataFrame from merging the index .cnr DataFrame and the parental .cnr files.
            selection_index (list): List of selected calls for the index patient.
            selection_father (list): List of selected calls for the father.
            selection_mother (list): List of selected calls for the mother.
            call_list (list): List of all possible calls.

        Returns:
            pd.DataFrame: Filtered DataFrame.
        """
        selection_index = selection_index or call_list
        selection_father = selection_father or call_list
        selection_mother = selection_mother or call_list

        filtered_df = trio_df.copy()
        filtered_df = filtered_df[filtered_df["call"].isin(selection_index)]
        filtered_df = filtered_df[filtered_df["call_f"].isin(selection_father)]
        filtered_df = filtered_df[filtered_df["call_m"].isin(selection_mother)]

        return filtered_df
