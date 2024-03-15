"""
File which contains plotting functions of the CNVizard
@author: Jeremias Krause / Matthias Begemann / Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd 
import numpy as np
import streamlit as st
import mpld3
import streamlit.components.v1 as components
import plotly.express as px
import plotly.graph_objects as go 

def df_to_list(df:pd.DataFrame) -> pd.DataFrame :
    """
    Functions which transfroms the log2 and depth columns to lists, used for plotting
    Input Arguments: 
        df (pandas DataFrame) -> .cnr Dataframe 
    Output Arguments:
        list_of_log2 (list) -> list containing the log2 values
        list_of_depth (list) -> list containing the depth values
    """
    col_of_log2 = df.groupby('exon')['log2'].apply(list)
    col_of_depth = df.groupby('exon')['depth'].apply(list)
    list_of_log2 = col_of_log2.tolist()
    list_of_depth = col_of_depth.tolist()
    return list_of_log2, list_of_depth

def index_ref_processor(df:pd.DataFrame,selected_gene:str):
    """
    Function used to extract the log2 and depth values from the index .cnr DataFrame for a specific gene 
    Input Arguments:
        df (pandas DataFrame) -> index .cnr DataFrame
        selected_gene (Str) -> Gene for which to filter
    Output Arguments:
        list_of_log2_index (list) -> list of log2 values from index for selected gene 
        exon_list (list) -> list of exons from index for selected gene
        list_of_depth_index (list) -> list of depth values from index for selected gene 
    """
    #Filter DataFrame for selected gene 
    index_gene_df = df[df['gene'] ==selected_gene]
    #Sort Dataframe by exon 
    index_gene_df = index_gene_df.sort_values(by=['exon'], ascending=True)
    #Extract log2/depth values 
    [list_of_log2_index, list_of_depth_index] = df_to_list(index_gene_df)
    #Get Exon names from index DataFrame
    exon_list =index_gene_df['exon'].unique().tolist()
    return list_of_log2_index,exon_list,list_of_depth_index

def plot_log2_for_gene_precomputed(gene:str,df_total:pd.DataFrame,reference_df:pd.DataFrame,sample_name:str):
    """
    Function used to create boxplot for log2 values, extracted from the reference dataFrame
    Subsequently the individual log2 values from the index cnr are plotted "on top" of the boxplot to show 
    the index log2 values compare to the reference log2 values
    Input Arguments:
        gene (str) -> selected gene 
        df_toal (pandas DataFrame) -> index .cnr DataFrame 
        reference_df (pandas DataFrame) -> containing the precalculated statistics, created by /reference_builder/referencemerge.py
        sample_name (str) -> index sample name, automatically extracted after uploading the index cnr file
    """
    #Create Plotly Graphical Object 
    fig = go.Figure()
    #Extract precomputed Statistics from the reference DataFrame
    selected_gene = reference_df[reference_df["gene"]==gene]
    listed_q1 = selected_gene["q1_log2"].tolist()
    listed_median = selected_gene["median_log2"].tolist()
    listed_q3 = selected_gene["q3_log2"].tolist()
    listed_min = selected_gene["actual_minimum_log2"].tolist()
    listed_max = selected_gene["actual_maximum_log2"].tolist()
    listed_mean = selected_gene["mean_log2"].tolist()
    listed_std = selected_gene["std_log2"].tolist()
    listed_exons = selected_gene["exon"].tolist()

    #For each exon of the selected gene plot a boxplot 
    fig.add_trace(go.Box(q1=listed_q1, median=listed_median,
                  q3=listed_q3, lowerfence=listed_min,
                  upperfence=listed_max,mean=listed_mean,x=listed_exons,showlegend=False,fillcolor="white", line= {"color": "black"}))
    #Add title to boxplot
    fig.update_layout(title="log2-plot-"+gene+"-"+sample_name)
    #Add Exon numbers on x-axis
    fig.update_xaxes(title_text="Exons",tickmode="linear")
    #Add Description for y-axis
    fig.update_yaxes(title_text="log-2 value",range=[-2,2])
    [list_of_log2_index,exon_list,list_of_depth_index] = index_ref_processor(df_total,gene)
    for log2_value,exon_name in zip(list_of_log2_index,listed_exons):
            for value in log2_value :
                y_list =[exon_name]
                x_list =[value]
                #Plot Index log2 value
                fig.add_trace(go.Scatter(y=x_list,x=y_list,mode="markers",marker=dict(color="red"),showlegend=False))
    #Add three horizontal lines to the log2-plot, which resemble the cut-off values for the cnv calls 
    fig.add_hline(y=0.3, line_width=1, line_dash="dash", line_color="blue", showlegend=True, name="CN>2")
    fig.add_hline(y=-0.4, line_width=1, line_dash="dash", line_color="red", showlegend=True, name="CN<2")
    fig.add_hline(y=-1.1, line_width=1, line_dash="dash", line_color="darkred", showlegend=True, name="CN<1")
    #Add legend which enables on&off-toggling for the horizontal lines 
    fig.update_layout(legend=dict(orientation="h", yanchor="middle",y=1.02,xanchor="left",x=0.05))
    #Visualize plotly graphical object in plotly
    if len(exon_list) > 30: 
        st.plotly_chart(fig,use_container_width=True)
    else:
        st.plotly_chart(fig)

def plot_depth_for_gene_precomputed(gene:str,df_total:pd.DataFrame,reference_df:pd.DataFrame,sample_name:str):
    """
    Function used to create boxplot for depth values, extracted from the reference dataFrame
    Subsequently the individual depth values from the index cnr are plotted "on top" of the boxplot to show 
    the index depth values compare to the reference log2 values
    Input Arguments:
        gene (str) -> selected gene 
        df_toal (pandas DataFrame) -> index .cnr DataFrame 
        reference_df (pandas DataFrame) -> containing the precalculated statistics, created by /reference_builder/referencemerge.py
        sample_name (str) -> index sample name, automatically extracted after uploading the index cnr file
    """
    #Create Plotly Graphical Object 
    fig = go.Figure()
    #Extract precomputed Statistics from the reference DataFrame
    selected_gene = reference_df[reference_df["gene"]==gene]
    listed_q1 = selected_gene["q1_depth"].tolist()
    listed_median = selected_gene["median_depth"].tolist()
    listed_q3 = selected_gene["q3_depth"].tolist()
    listed_min = selected_gene["actual_minimum_depth"].tolist()
    listed_max = selected_gene["actual_maximum_depth"].tolist()
    listed_mean = selected_gene["mean_depth"].tolist()
    listed_std = selected_gene["std_depth"].tolist()
    listed_exons = selected_gene["exon"].tolist()

    #For each exon of the selected gene plot a boxplot 
    fig.add_trace(go.Box(q1=listed_q1, median=listed_median,
                  q3=listed_q3, lowerfence=listed_min,
                  upperfence=listed_max,mean=listed_mean,x=listed_exons,showlegend=False,fillcolor="white", line= {"color": "black"}))
    #Add title to boxplot
    fig.update_layout(title="depth-plot-"+gene+"-"+sample_name)
    #Add Exon numbers on x-axis
    fig.update_xaxes(title_text="Exons",tickmode="linear")
    #Add Description for y-axis
    fig.update_yaxes(title_text="Depth value")
    [list_of_log2_index,exon_list,list_of_depth_index] = index_ref_processor(df_total,gene)
    for depth_value,exon_name in zip(list_of_depth_index,listed_exons):
            for value in depth_value :
                y_list =[exon_name]
                x_list =[value]
                #Plot index depth values
                fig.add_trace(go.Scatter(y=x_list,x=y_list,mode="markers",marker=dict(color="red"),showlegend=False))
    #Visualize plotly graphical object in plotly
    if len(exon_list) > 30: 
        st.plotly_chart(fig,use_container_width=True)
    else:
        st.plotly_chart(fig)



