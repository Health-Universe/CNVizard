"""
File used to run the CNVizard, a streamlit app that visualizes germline-copy-number-variants
@author: Jeremias Krause / Matthias Begemann / Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""


import streamlit as st
import pandas as pd
import os 
from io import BytesIO
import xlsxwriter
import numpy as np
import pyarrow
from CNV_Visualizer import cnv_visualizer
from CNV_Visualizer import styler
from CNV_Visualizer import exporter
from CNV_Visualizer import plotter
from CNV_Visualizer import helpers
import plotly.express as px
import csv
import cnvlib 
from cnvlib.commands import do_scatter
from pathlib import Path 
import dotenv

#pd.set_option("styler.render.max_elements", 1264048)

st.set_page_config(layout="wide", page_title="CNVizard", page_icon="CNVizard.png")

#load env file for igv_string
dotenv.load_dotenv()
if (dotenv.load_dotenv()) == True:
    igv_string = os.environ.get("APPSETTING_IGV_OUTLINK")

#Declare current_working dir 
current_working_dir = Path.cwd()

#Declare the path to an annotation_folder 
#The folder contains an omim.txt file which which contains OMIMP/OMIMG annotations mapped to gene
#path to omim .txt file 
omim_annotation_folder = str(current_working_dir) + "/Ressources/OMIM/omim.txt"
#Declare the path to an additional annotation_folder 
#This folder contains multiple .txt files which contain a predefined candigene panel
#path to candilist folder 
candi_annotation_folder = str(current_working_dir) + "/Ressources/Candi_lists/"
#Path to AnnotSV Column config 
tsv_input = str(current_working_dir) + "/Ressources/AnnotSV_Format/Annotsv_format.txt"

#Filter_options :
#List of chromosomes
chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
              "chr20", "chr21", "chr22", "chrX", "chrY"]
#List of possible Copy-Number calls 
call_list = [0, 1, 2, 3]
#import the omim.txt file as a pandas dataframe
gene_list = pd.read_csv(omim_annotation_folder, header=0,delimiter="\t")["gene"].tolist()

#Title of the App and Unicode for the wizard emoji
st.title("CNVizard")
#General Description
st.markdown("This is a streamlit webapp which provides analyis tools for genetic copy number variants")
#Description of uploader objects and the expected format of input files 
st.subheader("Upload")
st.markdown("Please upload a reference file, an individual .cnr file and an individual .bintest file, provided by CNV_Kit")

#Declare a streamlit column object which holds the three uploader objects for the input files 
cols = st.columns(3)
#Uploader for index .cnr file, which is created by CNVkit 
entered_cnr = cols[0].file_uploader(".cnr",accept_multiple_files=False)
#Uploder for bintest file, which is created by CNVkit
entered_bintest = cols[1].file_uploader("bintest",accept_multiple_files=False)
#Select NGS_Type 
ngs_type = cols[2].radio("Choose NGS Type", ["WES", "WGS"])

#The reference file are multiple concatenated .cnr files which have been formatted (the individual information have 
#been aggregated per exon and corresponding statics have been calculated which will be used in the creation of the boxplots,
#finally the reference file contains calculated in House frequencies on an exon level)
if ngs_type =="WGS":
    entered_reference = str(current_working_dir) + "/Ressources/References/genome_cnv_reference_large.parquet"
    entered_reference_bintest = str(current_working_dir) + "/Ressources/References/genome_cnv_reference_bintest_large.parquet"
else: 
    entered_reference = str(current_working_dir) + "/Ressources/References/exome_cnv_reference_large.parquet"
    entered_reference_bintest = str(current_working_dir) + "/Ressources/References/exome_cnv_reference_bintest_large.parquet"


#Once a index .cnr file is uploaded, extract sample name, which will later be used in the title of the plots 
if entered_cnr is not None:
    sample_name = entered_cnr.name.split(".")[0]
    if igv_string is not None:
        igv_string = igv_string.replace("samplename", sample_name)

    
#Description for configurations 
st.subheader("Configurations")
st.markdown("Please provide a sample Name and define how many consecutive exons have to be deleted/duplicated to show up, using the consecutive preset (default value=2)")

#Declare a streamlit column object which holds the two text input objects, where the user can define the minimal 
#amount of deletions and duplications (in exons) which will be displayed in the filtered dataframe 
cols2 = st.columns(2)
#text input object for deletion size 
entered_del_size = cols2[0].text_input("deletion_size")
#text input object for duplication size 
entered_dup_size = cols2[1].text_input("duplication_size")

#App Description which will be depicted in the sidebar 
st.sidebar.title("About")
st.sidebar.markdown("This streamlit webapp enables you to visualize copy-number-variant-data using dataframes and plots")

#The sidebar additionally enables the user to select a candigene (which will be loaded from the previously candi_annotation_folder)
st.sidebar.subheader("Candigene Selection")
#Use os to get all the defined candigene lists (sort them alphabetically)
candigenes_options = sorted(os.listdir(candi_annotation_folder))
#Declare streamlit selection object 
selected_candi = st.sidebar.radio("Please select the desired candigenes",candigenes_options)
#Display streamlit selection object 
st.sidebar.write(selected_candi)

#Index Dataframe visualization section
st.subheader("Dataframe Visualization")
#Declare a expander which holds different streamlit selection/input objects which enable the user to filter the .cnr dataframe
with st.expander("Filter"):
    #Declare a streamlit column object which holds the filters "Chromosome", "start", "end"
    cols3 = st.columns(3)
    #Declare a streamlit multiselect object and save the input in the variable chrom_selection
    chrom_selection = cols3[0].multiselect("Chromosome",chrom_list)
    #transform the variable chrom_selection to a list, which will be used to filter the dataframe
    chrom_selection = list(chrom_selection)
    #Declare a streamlit text input object and save the input in the variable start_selection
    start_selection = cols3[1].text_input("start")
    #Once a input for "start" has been made (and therefore start is not None & not empty) cast start_selection to int
    if start_selection is not None and start_selection != "":
        start_selection = int(start_selection)
    #Declare a streamlit text input object and save the input in the variable end_selection
    end_selection = cols3[2].text_input("end")
    #Once a input for "end" has been made (and therefore end is not None & not empty) cast end_selection to int
    if end_selection is not None and end_selection != "":
        end_selection = int(end_selection)
    #Both casts are performed for filtering purposes (the column of the dataframe which will be filtered contains values of type int)
    #Declare a streamlit column object which will hold the filters "exon", "min_depth" & "min_weight"
    cols4 = st.columns(3)
    #Declare a streamlit text input object and save the input in the variable depth_selection
    depth_selection = cols4[0].text_input("min_depth")
    #Once an input has been made for depth_selection has been made (and therefore depth_selection is not None & and not empty) cast depth_selection to float
    if depth_selection is not None and depth_selection !="":
        depth_selection = float(depth_selection)
    #Declare a streamlit text input and save the input in the variable weight_selection 
    weight_selection = cols4[1].text_input("min_weight")
    #Once an input has been made for weight_selection has been made (and therefore weight_selection is not None & and not empty) cast weight_selection to float
    if weight_selection is not None and weight_selection !="":
        weight_selection = float(weight_selection)
    #Both casts are performed for filtering purposes (the column of the dataframe which will be filtered contains values of type float)
    #Declare a streamlit multiselect object and save the input in the variable call_selection
    call_selection = cols4[2].multiselect("call",call_list)
    #Declare a streamlit column objcect which holds the filters "min_log2", "gene" and "het_del_selection"
    cols5 = st.columns(3)
    #Declare a streamlit text inpout and save the input in the variable log2_selection
    log2_selection = cols5[0].text_input("min_log2")
    #Once an input has been made for log2_selection (and therefore log2_selection is not None & not empty) cast log2_selection to float
    if log2_selection is not None and log2_selection !="":
        log2_selection = float(log2_selection)
    #Cast is done for filtering purposes (the column of the dataframe which will be filtered contains values of type float)
    #Declare a streamlit multiselect input and save the input in the variable gene_selection
    gene_selection = cols5[1].multiselect("gene",gene_list)
    #Transform gene_selection to list, which will be used to filter a dataframe 
    gene_selection = list(gene_selection)
    gene_selection = [gene.upper() for gene in gene_selection]
    #Declare a streamlit text_input object and save the input in the variable het_del_selection
    het_del_selection = cols5[2].text_input("max_het_del_freq")
    #Once an input has been made for het_del_selection (and therefore het_del_selection is not None & not empty) cast het_del_selection to float
    if het_del_selection is not None and het_del_selection !="":
        het_del_selection = float(het_del_selection)
    #Cast is done for filtering purposes (the column of the dataframe which will be filtered contains values of type float)
    #Declare a streamlit column object which will hold the filters "max_hom_del_freq","max_dup_freq" (the created column object is still of size 3 for aesthetic reasons, an object of the size two would not align well with the other column objects)
    cols7 = st.columns(3)
    #Declare a streamlit text input and save the input to the variable hom_del_selection
    hom_del_selection = cols7[0].text_input("max_hom_del_freq")
    #Once an input has been made for hom_del_selection (and therefore hom_del_selection is not None & not empty) cast hom_del_selection to float
    if hom_del_selection is not None and hom_del_selection !="":
        hom_del_selection = float(hom_del_selection)
    #Cast is done for filtering purposes (the column of the dataframe which will be filtered contains values of type float)
    #Declare a streamlit text input and save the input to the variable dup_selection
    dup_selection = cols7[1].text_input("max_dup_freq")
    if dup_selection is not None and dup_selection !="":
        dup_selection = float(dup_selection)
    #Once an input has been made for dup_selection (and therefore dup_selection is not None & not empty) cast dup_selection to float
    #Cast is done for filtering purposes (the column of the dataframe which will be filtered contains values of type float)




#Define a list of possible preset dataframes from which the user can choose 
list_of_possible_dataframes = ["total", "bintest", "hom_del", "total_candi", "bintest_candi", "consecutive_del", "consecutive_dup"]
#Declare a streamlit selectbox object which contains the preset options
df_to_be_displayed=st.selectbox("Select dataframe which shall be displayed",list_of_possible_dataframes)

#Once a reference parquet file has been uploaded, import the reference parquet file into a pandas dataframe 
if entered_reference is not None:
    reference_df = pd.read_parquet(entered_reference)

#Once a .cnr file has been uploaded, import the .cnr file into a pandas dataframe 
if entered_cnr is not None:
    cnr_df = pd.read_csv(entered_cnr,header=0,delimiter="\t")

#Once a bintest file has been uploaded, import the bintest file into a pandas dataframe 
if entered_bintest is not None:
    bintest_df = pd.read_csv(entered_bintest,header=0,delimiter="\t")

#Once a reference dataframe (entered_df), a .cnr dataframe (cnr_df) & a bintest df (bintest_df) are present initialise an instance 
# of the class cnv_visualizer which will be used to process & formate the entered input files 
# Subsequently apply the filters entered by the user 
# Finally display the currently selected preset dataframe (the user can switch between the possible presets)
if (entered_reference is not None) and (entered_cnr is not None) and (entered_bintest is not None) and (entered_reference_bintest is not None):
    #Initialise an instance of the cnv_visualizer class (required arguments : reference dataframe & bintest dataframe)
    cnv_visualizer_initialised = cnv_visualizer.cnv_visualizer(reference_df,cnr_df,bintest_df)
    #Import the omim .txt file & the selected candi .txt file into a pandas dataframe and merge them with the bintest and .cnr dataframe 
    #Additionaly reorder some columns insinde the .cnr & bintest dataframe and change the format of some columns
    [omim_df,candi_df,cnr_db,bintest_db]=cnv_visualizer_initialised.format_df(omim_annotation_folder,candi_annotation_folder + selected_candi)

    #Create a new pandas dataframe which extracts the previously mentioned inhouse frequencies from the reference file 
    call_df = pd.DataFrame()
    call_df["gene"] = reference_df["gene"]
    call_df["exon"] = reference_df["exon"]
    call_df["het_del_frequency"] = reference_df["het_del_frequency"]
    call_df["hom_del_frequency"] = reference_df["hom_del_frequency"]
    call_df["dup_frequency"] = reference_df["dup_frequency"]

    #Create a new pandas dataframe which holds the inhouse frequencies for the bintest
    bintest_inhouse_df = pd.read_parquet(entered_reference_bintest)

    #merge the .cnr dataframe with the in house frequency dataframe 
    cnr_db = pd.merge(cnr_db,call_df,how="left", left_on=["gene","exon"], right_on=["gene","exon"])

    #merge the bintest dataframe with in house frequency bintest dataframe 
    bintest_db = pd.merge(bintest_db, bintest_inhouse_df, how="left", left_on=["gene", "exon"], right_on=["gene", "exon"])
    bintest_db["het_del_frequency"] = bintest_db["het_del_frequency"].fillna(0)
    bintest_db["hom_del_frequency"] = bintest_db["hom_del_frequency"].fillna(0)
    bintest_db["dup_frequency"] = bintest_db["dup_frequency"].fillna(0)

    #Apply the filters selected by the user, empty filters will be ignored. The filters only apply to the .cnr dataframe 
    cnr_db_filtered = cnv_visualizer_initialised.apply_filters(cnr_db, start_selection, end_selection, depth_selection,
                                                      weight_selection,chrom_selection,call_selection, log2_selection,
                                                      gene_selection, chrom_list, call_list, gene_list,het_del_selection,
                                                      hom_del_selection,dup_selection)
    
 
    #rename the column "chromosome" in the .cnr and bintest dataframe for aesthetic purposes 
    cnr_db_filtered = cnr_db_filtered.rename(columns={"chromosome": "chr"})
    cnr_db_filtered = cnr_db_filtered.round(2)
    bintest_db = bintest_db.rename(columns={"chromosome": "chr"})
    bintest_db = bintest_db.round(2)
    #if an env file was created, create a column for an igv outlink 
    if igv_string is not None:
        cnr_db_filtered["IGV_outlink"] = igv_string + cnr_db_filtered["chr"] + ":" + cnr_db_filtered["start"].astype(str)
        bintest_db["IGV_outlink"] = igv_string + bintest_db["chr"] + ":" + bintest_db["start"].astype(str)


    #This section responds to the users preset selection, the selection is observed by storing the selected value 
    #in the variable df_to_be_displayed
    #addtionally the currently viewed dataframe will be save in a variable download_filter and the name of the preset in the variable download_filter_name
    #download_filter and download_filter_name will be used to enable the viewer a download option for the currently selected preset
    #Display the unfiltered .cnr dataframe
    if df_to_be_displayed == "total":
        st.dataframe(cnr_db_filtered)
        download_filter = cnr_db_filtered
        download_filter_name = "total"
    #Display the unfiltered bintest dataframe, additionally a pandas styler (defined in CNV_Visualizer/styler.py) is applied 
    elif df_to_be_displayed == "bintest":
        st.dataframe(bintest_db.style.pipe(styler.make_pretty))
        download_filter = bintest_db
        download_filter_name = "bintest"
    #Display the .cnr dataframe filtered for homozygous deletions, additionally a pandas styler (defined in CNV_Visualizer/styler.py) is applied  
    elif df_to_be_displayed == "hom_del":
        hom_del_df_total = cnv_visualizer_initialised.filter_for_deletions_hom(cnr_db_filtered)
        st.dataframe(hom_del_df_total.style.pipe(styler.make_pretty))
        download_filter = hom_del_df_total
        download_filter_name = "hom_del"
    #Display the .cnr dataframe filtered with the previously defined candigene list, additionally a pandas styler (defined in CNV_Visualizer/styler.py) is applied  
    elif df_to_be_displayed == "total_candi":
        filtered_for_candi_df_total = cnv_visualizer_initialised.filter_for_candi_cnvs(cnr_db_filtered,candi_df)
        st.dataframe(filtered_for_candi_df_total.style.pipe(styler.make_pretty))
        download_filter = filtered_for_candi_df_total
        download_filter_name = "total_candi"
    #Display the bintest dataframe filtered with the previously defined candigene list, additionally a pandas styler (defined in CNV_Visualizer/styler.py) is applied  
    elif df_to_be_displayed == "bintest_candi":
        filtered_for_candi_df_bintest = cnv_visualizer_initialised.filter_for_candi_cnvs(bintest_db,candi_df)
        st.dataframe(filtered_for_candi_df_bintest.style.pipe(styler.make_pretty))
        download_filter = filtered_for_candi_df_bintest
        download_filter_name = "bintest_candi"
    #Display the .cnr dataframe filtered for consecutive deletions (the amount of required consecutive exons is previously defined in the configurations section), 
    #additionally a pandas styler (defined in CNV_Visualizer/styler.py) is applied  
    elif df_to_be_displayed == "consecutive_del":
        filtered_for_del = cnv_visualizer_initialised.filter_for_deletions(cnr_db_filtered)
        filtered_for_del = cnv_visualizer_initialised.prepare_filter_for_consecutive_cnvs(filtered_for_del)
        filtered_for_del = cnv_visualizer_initialised.filter_for_consecutive_cnvs(filtered_for_del,"del",entered_del_size,entered_dup_size)
        st.dataframe(filtered_for_del.style.pipe(styler.make_pretty))
        download_filter = filtered_for_del
        download_filter_name = "consecutive_del"
    #Display the .cnr dataframe filtered for consecutive duplications (the amount of required consecutive exons is previously defined in the configurations section), 
    #additionally a pandas styler (defined in CNV_Visualizer/styler.py) is applied 
    elif df_to_be_displayed == "consecutive_dup":
        filtered_for_dup = cnv_visualizer_initialised.filter_for_duplications(cnr_db_filtered)
        filtered_for_dup = cnv_visualizer_initialised.prepare_filter_for_consecutive_cnvs(filtered_for_dup)
        filtered_for_dup = cnv_visualizer_initialised.filter_for_consecutive_cnvs(filtered_for_dup,"dup",entered_del_size,entered_dup_size)
        st.dataframe(filtered_for_dup.style.pipe(styler.make_pretty))
        download_filter = filtered_for_dup
        download_filter_name = "consecutive_dup"
    
    #Prepare a streamlit column object which will hold the "prepare for download buttons"
    download_columns = st.columns(2)
    download_preparator_all = download_columns[0].button("prepare for download (all)")
    download_preparator_filtered = download_columns[1].button("prepare for download (individual filtered table)")

    #Prepare two additonal streamlit column objects which will hold the actual download buttons and the download button messages
    download_message_columns = st.columns(2)
    download_button_columns = st.columns(2)

    
    #Prepare the preset dataframes for an export as excel file 
    #To improve performance this is only triggered once the user clicks the "prepare for download button"
    #The required export functions are defined in /CNV_Visualizer/exporter.py     
    if download_preparator_all:
        hom_del_df_total = cnv_visualizer_initialised.filter_for_deletions_hom(cnr_db)
        filtered_for_candi_df_total = cnv_visualizer_initialised.filter_for_candi_cnvs(cnr_db,candi_df)
        filtered_for_candi_df_bintest = cnv_visualizer_initialised.filter_for_candi_cnvs(bintest_db,candi_df)
        filtered_for_del = cnv_visualizer_initialised.filter_for_deletions(cnr_db)
        filtered_for_del = cnv_visualizer_initialised.prepare_filter_for_consecutive_cnvs(filtered_for_del)
        filtered_for_del = cnv_visualizer_initialised.filter_for_consecutive_cnvs(filtered_for_del,"del",entered_del_size,entered_dup_size)
        filtered_for_dup = cnv_visualizer_initialised.filter_for_duplications(cnr_db)
        filtered_for_dup = cnv_visualizer_initialised.prepare_filter_for_consecutive_cnvs(filtered_for_dup)
        filtered_for_dup = cnv_visualizer_initialised.filter_for_consecutive_cnvs(filtered_for_dup,"dup",entered_del_size,entered_dup_size)

        table_exporter =exporter.cnv_exporter()
        to_be_exported = table_exporter.save_tables_as_excel(cnr_db,bintest_db,hom_del_df_total,filtered_for_candi_df_total,
                                                      filtered_for_candi_df_bintest,filtered_for_del,filtered_for_dup)
        download_button_columns[0].download_button(label="Download",
                           data=to_be_exported,
                           file_name = sample_name + "_df.xlsx")
    else:
        download_message_columns[0].write("Click button to prepare download")

    #Prepare the individual filtered dataframe for an export as excel file 
    #To improve performance this is only triggered once the user clicks the "prepare for download button"
    #The required export functions are defined in /CNV_Visualizer/exporter.py
    
    if download_preparator_filtered:
        table_exporter =exporter.cnv_exporter()
        to_be_exported_filtered = table_exporter.save_filtered_table_as_excel(download_filter,download_filter_name)
        download_button_columns[1].download_button(label="Download_filtered",
                                                   data=to_be_exported_filtered,
                                                   file_name = sample_name + "filtered_df.xlsx")
    else:
        download_message_columns[1].write("Click button to prepare download")
        
#In this section the user can plot Boxplots for user defined genes 
#The required statistical information are stored in the reference file 
st.subheader("Plot log2 and depth for selected genes")
entered_gene = st.multiselect("gene", gene_list,max_selections=1)
if len(entered_gene)==0:
    entered_gene=None
else:
    entered_gene = entered_gene[0]
    entered_gene = entered_gene.upper()

#Once a reference df, a .cnr df and a bintest df are present and the user entered a gene which shall be plotted 
# a boxplot for the log2 and the depth value are created 
#The boxplots are created as plotly graph objects 
#To compare the index log2 and weight from the .cnr file to the reference, the index values are plotted ontop of the 
#boxplots as a plotly scatterplot
if (entered_reference is not None) and (entered_cnr is not None) and (entered_bintest is not None) and (entered_gene is not None) and (entered_gene != ""):
    gene_plotter = plotter.cnv_plotter()
    gene_plotter.plot_log2_for_gene_precomputed(entered_gene,cnr_db,reference_df,sample_name)
    gene_plotter.plot_depth_for_gene_precomputed(entered_gene,cnr_db,reference_df,sample_name)

#In this section the user has the option to additionally upload .cnr file from the index patients parents     
st.subheader ("Load additional .cnr files from the index patients corresponding parents")
#Declare a streamlit expander which contains filter options stored in streamlit column objects 
with st.expander("Filter for Trio"):
    #Declare a streamlit column object which contains the filter "call index", "call father" and "call mother"
    #This can be used to filter for different call constellations, e.g. denovo calling (index call!=0 & mother call ==0 & father call == 0)
    cols_for_trio_filter = st.columns(3)
    #Declare a streamlit multiselect object and save the input in the variable call_selection_index
    call_selection_index = cols_for_trio_filter[0].multiselect("call index",call_list)
    #Declare a streamlit multiselect object and save the input in the variable call_selection_father
    call_selection_father = cols_for_trio_filter[1].multiselect("call father",call_list)
    #Declare a streamlit multiselect object and save the input in the variable call_selection_mother 
    call_selection_mother = cols_for_trio_filter[2].multiselect("call mother",call_list)
#Create a streamlit column object which holds the file_uploader objects, used to upload the index patients parental .cnr files
cols_trio = st.columns(3)
#Create a streamlit file_uploader object and save the input in the variable father_cnr
father_cnr = cols_trio[0].file_uploader("Father .cnr",accept_multiple_files=False)
#Create a streamlit file_uploader object and save the input in the variable mother_cnr
mother_cnr = cols_trio[1].file_uploader("Mother .cnr",accept_multiple_files=False)

#Once a reference df, a cnr df, a bintest df, a father .cnr file and a mother .cnr file are present import the parental 
#.cnr files into a pandas dataframe and format those dataframes 
if (entered_reference is not None) and (entered_cnr is not None) and (entered_bintest is not None) and (father_cnr is not None) and (mother_cnr is not None):
    
    #import the father_cnr_df and format it using the previously initialise cnv_visualizer instance and drop unnecessary columns
    father_cnr_df = pd.read_csv(father_cnr,header=0,delimiter="\t")
    father_cnr_df = cnv_visualizer_initialised.prepare_parent_cnv(father_cnr_df)
    father_cnr_df = father_cnr_df.drop(["chromosome","start","end"],axis=1)
    father_cnr_df = father_cnr_df.rename(columns={"gene": "gene_f", "exon": "exon_f", "depth": "depth_f",
                                                  "weight": "weight_f", "call": "call_f", "log2": "log2_f", "squaredvalue": "squaredvalue_f"})
    #round float values to 2 decimals
    father_cnr_df = father_cnr_df.round(2)

    #import the mother_cnr_df and format it using the previously initialise cnv_visualizer instance and drop unnecessary columns
    mother_cnr_df = pd.read_csv(mother_cnr,header=0,delimiter="\t")
    mother_cnr_df = cnv_visualizer_initialised.prepare_parent_cnv(mother_cnr_df)
    mother_cnr_df = mother_cnr_df.drop(["chromosome","start","end"],axis=1)
    mother_cnr_df = mother_cnr_df.rename(columns={"gene": "gene_m", "exon": "exon_m", "depth": "depth_m",
                                                  "weight": "weight_m", "call": "call_m", "log2": "log2_m", "squaredvalue": "squaredvalue_m"})
    #round float values to 2 decimals
    mother_cnr_df = mother_cnr_df.round(2)

    #merge the parental dataframes with the index dataframe
    trio_cnr_df = pd.merge(cnr_db,father_cnr_df,how="left",left_on=["gene","exon"],right_on=["gene_f","exon_f"])
    trio_cnr_df = pd.merge(trio_cnr_df,mother_cnr_df,how="left",left_on=["gene","exon"],right_on=["gene_m","exon_m"])

    #renmae the column "chromosome" to "chr" for aesthetic purposes 
    trio_cnr_df = trio_cnr_df.rename(columns={"chromosome": "chr"})
    #Round float values to 2 decimals
    trio_cnr_df = trio_cnr_df.round(2)
    #Apply the previously user defined filters, empty filters are ignored
    trio_cnr_df_filtered = cnv_visualizer_initialised.apply_trio_filters(trio_cnr_df,call_selection_index,call_selection_father,call_selection_mother,call_list)
    #Visualize filtered trio dataframe
    st.dataframe(trio_cnr_df_filtered)

#In this section the user can plot a genome-wide scatter plot (using cnvlibs "do_scatter" function)
st.subheader("Plot genome-wide or chromosome-wide scatter plot")
cols_cns_upload= st.columns(3)
entered_cns = cols_cns_upload[0].file_uploader(".cns file", accept_multiple_files=False)
entered_cnr_new = cols_cns_upload[1].file_uploader(".cnr file", accept_multiple_files=False)

#Add the option "All" and chrom_list to a new variable, from which the user can select to plot an individual chromosome or a genome wide scatter
scatter_options = chrom_list
scatter_options.append("All")

with st.expander("Select chromosome"):
    selected_scatter = st.radio("Please select an individual chromosome/all", scatter_options)



if (entered_cnr_new is not None) and (entered_cns is not None):
    if selected_scatter == "All":
        input_cnr = cnvlib.read(entered_cnr_new)
        input_cns = cnvlib.read(entered_cns)
        fig_scatter = do_scatter(input_cnr,segments=input_cns,y_min=-2,y_max=2)
        fig_scatter.set_figwidth(12)
        st.pyplot(fig_scatter)
    else: 
        input_cnr = cnvlib.read(entered_cnr_new)
        input_cns = cnvlib.read(entered_cns)
        fig_scatter = do_scatter(input_cnr,segments=input_cns,y_min=-2,y_max=2, show_range=selected_scatter)
        fig_scatter.set_figwidth(12)
        st.pyplot(fig_scatter)


#AnnotSV uses a differing chromosome nomenclature (without "chr") therefore create a new chromosome list which will be used for filtering
chromosome_list_cnv = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
#AnnotSV recognises only DEL/DUP as Call
cnv_type = ["DEL","DUP"]
#Create a list containing the five acmg classes, used for filtering
acmg_class =[1,2,3,4,5]
#Create a streamlit upload object, used to upload the .tsv file annotated by AnnotSV
st.subheader("Load annotated .tsv file created by AnnotSV")
entered_tsv_file = st.file_uploader(".tsv file",accept_multiple_files=False)

#Define a streamlit expander which contains multiple filters, which will be used to filter the tsv object
with st.expander("Filter for .tsv file"):
    #Define a stremlit column object which contains the filter "chromosome", "CNV_Type" and "ACMG_Class"
    cols6 = st.columns(3)
    #Create a streamlit multiselect object and save the output in the variable entered_cnv_chrom
    entered_cnv_chrom = cols6[0].multiselect("Chromosome",chromosome_list_cnv)
    #Create a streamlit multiselect object and save the output in the variable entered_cnv_type
    entered_cnv_type = cols6[1].multiselect("CNV_Type",cnv_type)
    #Create a streamlit multiselect object and save the ouput in the variable entered_acmg_class
    entered_acmg_class = cols6[2].multiselect("ACMG_Class",acmg_class)

#Once a .tsv file, annotated by AnnotSV has been uploaded 
#Import the .tsv file into a dataframe 
#The following step comes down to preference, to ensure compatibility it is assumed that 
#AnnotSV is used with all of the standard Annotations 
#For visibility purposes multiple columns are dropped and the remaining columns are rearranged
#Import Annotsv Column file (a column file, which contains a single column in which the preferred columns from Annotsv raw output are listed)
if entered_tsv_file is not None:
    if sample_name is None in globals():
        sample_name = entered_cnr.name.split(".")[0]
    tsv_df = pd.read_csv(entered_tsv_file,header=0,delimiter="\t")
    all_columns_tsv = tsv_df.columns.tolist()
    
    with open (tsv_input, "rt") as column_file:
        lines=[]
        for column in column_file:
            column = column.strip()
            lines.append(column)

    tsv_df = tsv_df[lines]
    #Create a Chromosome list which will be used to sort the pandas dataframe by turning it into a categorical
    sort_list_for_categorical = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","X","Y"]
    tsv_df["SV_chrom"] = pd.Categorical(tsv_df["SV_chrom"], sort_list_for_categorical)
    tsv_df = tsv_df.sort_values("SV_chrom")
    #After sorting the dataframe reset the index
    tsv_df = tsv_df.reset_index(drop=True)
    #Apply the user defined filters
    filtered_tsv = helpers.filter_tsv(tsv_df,chromosome_list_cnv,cnv_type,acmg_class,entered_cnv_chrom,entered_cnv_type,entered_acmg_class)
    filtered_tsv = filtered_tsv.fillna(".")
    #Display the filtered .tsv dataframe
    st.write(filtered_tsv)

    #Declate a streamlit button which the user can click to prepare the .tsv file for export as an excel file 
    #Because of Performance reasons this is only triggered if clicked by the user 
    #The export functions are defined in /CNV_Visualizer/exporter.py
    if st.button("prepare for download of annotated tsv data"):
        table_exporter = exporter.cnv_exporter()
        to_be_exported = table_exporter.save_tables_as_excel_tsv(filtered_tsv)
        st.download_button(label="Download",
                           data=to_be_exported,
                           file_name = sample_name + "_annotated_df.xlsx")
    else:
        st.write("Click button to prepare download")
    
