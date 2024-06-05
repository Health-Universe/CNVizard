"""
Run CNVizard, a Streamlit app that visualizes germline copy number variants.
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import streamlit as st
import pandas as pd
import os
from pathlib import Path
import dotenv
from cnvizard.cnv_visualizer import styler, exporter, plotter, helpers, visualizer
import cnvlib
from cnvlib.commands import do_scatter

def main():
    # Load environment variables for IGV link
    dotenv.load_dotenv()
    igv_string = os.getenv("APPSETTING_IGV_OUTLINK")

    # Set Streamlit page configuration
    st.set_page_config(layout="wide", page_title="CNVizard", page_icon="CNVizard.png")

    # Define current working directory
    current_working_dir = Path.cwd()

    # Define paths to resources
    resources_dir = current_working_dir / "resources"
    omim_annotation_path = resources_dir / "omim/omim.txt"
    candi_annotation_dir = resources_dir / "candidate_lists"
    annotsv_format_path = resources_dir / "annotsv_format.txt"

    # Filter options
    chrom_list = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    call_list = [0, 1, 2, 3]
    gene_list = pd.read_csv(omim_annotation_path, sep="\t")["gene"].tolist()

    # Streamlit App Title
    st.title("CNVizard")
    st.markdown("This is a Streamlit web app providing analysis tools for genetic copy number variants.")

    # File uploaders
    st.subheader("Upload")
    st.markdown("Please upload a reference file, an individual .cnr file, and an individual .bintest file, provided by CNVkit.")

    cols = st.columns(3)
    entered_cnr = cols[0].file_uploader(".cnr", type=["txt", "cnr"])
    entered_bintest = cols[1].file_uploader("bintest", type=["txt", "bintest"])
    ngs_type = cols[2].radio("Choose NGS Type", ["WES", "WGS"])

    # Determine reference files based on NGS type
    if ngs_type == "WGS":
        reference_path = resources_dir / "references/genome_cnv_reference_large.parquet"
        reference_bintest_path = resources_dir / "references/genome_cnv_reference_bintest_large.parquet"
    else:
        reference_path = resources_dir / "references/exome_cnv_reference_large.parquet"
        reference_bintest_path = resources_dir / "references/exome_cnv_reference_bintest_large.parquet"

    if entered_cnr:
        sample_name = entered_cnr.name.split(".")[0]
        if igv_string:
            igv_string = igv_string.replace("samplename", sample_name)

    st.subheader("Configurations")
    st.markdown("Provide a sample name and define the number of consecutive exons to display (default value = 2).")

    cols2 = st.columns(2)
    entered_del_size = cols2[0].text_input("deletion_size", value="2")
    entered_dup_size = cols2[1].text_input("duplication_size", value="2")

    # Sidebar configuration
    st.sidebar.title("About")
    st.sidebar.markdown("This Streamlit web app enables you to visualize copy-number-variant data using dataframes and plots.")

    st.sidebar.subheader("Candigene Selection")
    candigene_options = sorted(os.listdir(candi_annotation_dir))
    selected_candi = st.sidebar.radio("Select the desired candigene", candigene_options)
    st.sidebar.write(selected_candi)

    st.subheader("Dataframe Visualization")
    with st.expander("Filter"):
        cols3 = st.columns(3)
        chrom_selection = cols3[0].multiselect("Chromosome", chrom_list)
        start_selection = cols3[1].text_input("start", value="")
        end_selection = cols3[2].text_input("end", value="")
        
        cols4 = st.columns(3)
        depth_selection = cols4[0].text_input("min_depth", value="")
        weight_selection = cols4[1].text_input("min_weight", value="")
        call_selection = cols4[2].multiselect("call", call_list)
        
        cols5 = st.columns(3)
        log2_selection = cols5[0].text_input("min_log2", value="")
        gene_selection = cols5[1].multiselect("gene", gene_list)
        het_del_selection = cols5[2].text_input("max_het_del_freq", value="")
        
        cols7 = st.columns(3)
        hom_del_selection = cols7[0].text_input("max_hom_del_freq", value="")
        dup_selection = cols7[1].text_input("max_dup_freq", value="")

    list_of_possible_dataframes = ["total", "bintest", "hom_del", "total_candi", "bintest_candi", "consecutive_del", "consecutive_dup"]
    df_to_be_displayed = st.selectbox("Select dataframe to display", list_of_possible_dataframes)

    # Load reference files
    if reference_path.exists():
        reference_df = pd.read_parquet(reference_path)

    if entered_cnr:
        cnr_df = pd.read_csv(entered_cnr, delimiter="\t")

    if entered_bintest:
        bintest_df = pd.read_csv(entered_bintest, delimiter="\t")

    if reference_df is not None and cnr_df is not None and bintest_df is not None:
        cnv_visualizer_instance = visualizer.CNVVisualizer(reference_df, cnr_df, bintest_df)
        omim_df, candi_df, cnr_db, bintest_db = cnv_visualizer_instance.format_df(
            omim_annotation_path, candi_annotation_dir / selected_candi
        )

        call_df = reference_df[["gene", "exon", "het_del_frequency", "hom_del_frequency", "dup_frequency"]]
        bintest_inhouse_df = pd.read_parquet(reference_bintest_path)

        cnr_db = pd.merge(cnr_db, call_df, on=["gene", "exon"], how="left")
        bintest_db = pd.merge(bintest_db, bintest_inhouse_df, on=["gene", "exon"], how="left")
        bintest_db.fillna(0, inplace=True)

        cnr_db_filtered = cnv_visualizer_instance.apply_filters(
            cnr_db, start_selection, end_selection, depth_selection, weight_selection,
            chrom_selection, call_selection, log2_selection, gene_selection,
            chrom_list, call_list, gene_list, het_del_selection, hom_del_selection, dup_selection
        )

        cnr_db_filtered.rename(columns={"chromosome": "chr"}, inplace=True)
        cnr_db_filtered = cnr_db_filtered.round(2)
        bintest_db.rename(columns={"chromosome": "chr"}, inplace=True)
        bintest_db = bintest_db.round(2)

        if igv_string:
            cnr_db_filtered["IGV_outlink"] = igv_string + cnr_db_filtered["chr"] + ":" + cnr_db_filtered["start"].astype(str)
            bintest_db["IGV_outlink"] = igv_string + bintest_db["chr"] + ":" + bintest_db["start"].astype(str)

        # Dataframe display logic
        display_mapping = {
            "total": cnr_db_filtered,
            "bintest": bintest_db,
            "hom_del": cnv_visualizer_instance.filter_for_deletions_hom(cnr_db_filtered),
            "total_candi": cnv_visualizer_instance.filter_for_candi_cnvs(cnr_db_filtered, candi_df),
            "bintest_candi": cnv_visualizer_instance.filter_for_candi_cnvs(bintest_db, candi_df),
            "consecutive_del": cnv_visualizer_instance.filter_for_consecutive_cnvs(
                cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.filter_for_deletions(cnr_db_filtered)
                ), "del", entered_del_size, entered_dup_size
            ),
            "consecutive_dup": cnv_visualizer_instance.filter_for_consecutive_cnvs(
                cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.filter_for_duplications(cnr_db_filtered)
                ), "dup", entered_del_size, entered_dup_size
            ),
        }

        download_filter = display_mapping[df_to_be_displayed]
        st.dataframe(display_mapping[df_to_be_displayed].style.pipe(styler.make_pretty) if df_to_be_displayed != "total" else display_mapping[df_to_be_displayed])

        # Download buttons
        download_columns = st.columns(2)
        download_preparator_all = download_columns[0].button("Prepare for download (all)")
        download_preparator_filtered = download_columns[1].button("Prepare for download (filtered)")

        download_message_columns = st.columns(2)
        download_button_columns = st.columns(2)

        if download_preparator_all:
            tables_to_export = [
                cnr_db, bintest_db,
                cnv_visualizer_instance.filter_for_deletions_hom(cnr_db),
                cnv_visualizer_instance.filter_for_candi_cnvs(cnr_db, candi_df),
                cnv_visualizer_instance.filter_for_candi_cnvs(bintest_db, candi_df),
                cnv_visualizer_instance.filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                        cnv_visualizer_instance.filter_for_deletions(cnr_db)
                    ), "del", entered_del_size, entered_dup_size
                ),
                cnv_visualizer_instance.filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                        cnv_visualizer_instance.filter_for_duplications(cnr_db)
                    ), "dup", entered_del_size, entered_dup_size
                ),
            ]
            table_exporter = exporter.CNVExporter()
            export_data = table_exporter.save_tables_as_excel(*tables_to_export)
            download_button_columns[0].download_button(label="Download", data=export_data, file_name=f"{sample_name}_df.xlsx")
        else:
            download_message_columns[0].write("Click button to prepare download")

        if download_preparator_filtered:
            table_exporter = exporter.CNVExporter()
            export_data_filtered = table_exporter.save_filtered_table_as_excel(download_filter, df_to_be_displayed)
            download_button_columns[1].download_button(label="Download filtered", data=export_data_filtered, file_name=f"{sample_name}_filtered_df.xlsx")
        else:
            download_message_columns[1].write("Click button to prepare download")

    st.subheader("Plot log2 and depth for selected genes")
    entered_gene = st.multiselect("gene", gene_list, max_selections=1)
    entered_gene = entered_gene[0].upper() if entered_gene else None

    if reference_df is not None and cnr_df is not None and bintest_df is not None and entered_gene:
        gene_plotter = plotter.CNVPlotter()
        gene_plotter.plot_log2_for_gene_precomputed(entered_gene, cnr_db, reference_df, sample_name)
        gene_plotter.plot_depth_for_gene_precomputed(entered_gene, cnr_db, reference_df, sample_name)

    st.subheader("Load additional .cnr files from the index patient's parents")
    with st.expander("Filter for Trio"):
        cols_for_trio_filter = st.columns(3)
        call_selection_index = cols_for_trio_filter[0].multiselect("call index", call_list)
        call_selection_father = cols_for_trio_filter[1].multiselect("call father", call_list)
        call_selection_mother = cols_for_trio_filter[2].multiselect("call mother", call_list)

    cols_trio = st.columns(3)
    father_cnr = cols_trio[0].file_uploader("Father .cnr", type=["txt", "cnr"])
    mother_cnr = cols_trio[1].file_uploader("Mother .cnr", type=["txt", "cnr"])

    if reference_df is not None and cnr_df is not None and bintest_df is not None and father_cnr and mother_cnr:
        father_cnr_df = cnv_visualizer_instance.prepare_parent_cnv(pd.read_csv(father_cnr, delimiter="\t"))
        mother_cnr_df = cnv_visualizer_instance.prepare_parent_cnv(pd.read_csv(mother_cnr, delimiter="\t"))

        father_cnr_df = father_cnr_df.drop(["chromosome", "start", "end"], axis=1).rename(
            columns={col: f"{col}_f" for col in father_cnr_df.columns if col not in ["gene", "exon"]}).round(2)
        mother_cnr_df = mother_cnr_df.drop(["chromosome", "start", "end"], axis=1).rename(
            columns={col: f"{col}_m" for col in mother_cnr_df.columns if col not in ["gene", "exon"]}).round(2)

        trio_cnr_df = cnr_db.merge(father_cnr_df, left_on=["gene", "exon"], right_on=["gene_f", "exon_f"], how="left").merge(
            mother_cnr_df, left_on=["gene", "exon"], right_on=["gene_m", "exon_m"], how="left")

        trio_cnr_df.rename(columns={"chromosome": "chr"}, inplace=True)
        trio_cnr_df = trio_cnr_df.round(2)
        trio_cnr_df_filtered = cnv_visualizer_instance.apply_trio_filters(
            trio_cnr_df, call_selection_index, call_selection_father, call_selection_mother, call_list)
        st.dataframe(trio_cnr_df_filtered)

    st.subheader("Plot genome-wide or chromosome-wide scatter plot")
    cols_cns_upload = st.columns(3)
    entered_cns = cols_cns_upload[0].file_uploader(".cns file", type=["txt", "cns"])
    entered_cnr_new = cols_cns_upload[1].file_uploader(".cnr file", type=["txt", "cnr"])

    scatter_options = chrom_list + ["All"]
    selected_scatter = st.radio("Select chromosome or all", scatter_options)

    if entered_cns and entered_cnr_new:
        input_cnr = cnvlib.read(entered_cnr_new)
        input_cns = cnvlib.read(entered_cns)
        fig_scatter = do_scatter(input_cnr, segments=input_cns, y_min=-2, y_max=2, show_range=None if selected_scatter == "All" else selected_scatter)
        fig_scatter.set_figwidth(12)
        st.pyplot(fig_scatter)

    chromosome_list_cnv = [str(i) for i in range(1, 23)] + ["X", "Y"]
    cnv_type = ["DEL", "DUP"]
    acmg_class = [1, 2, 3, 4, 5]
    st.subheader("Load annotated .tsv file created by AnnotSV")
    entered_tsv_file = st.file_uploader(".tsv file", type=["tsv"])

    with st.expander("Filter for .tsv file"):
        cols6 = st.columns(3)
        entered_cnv_chrom = cols6[0].multiselect("Chromosome", chromosome_list_cnv)
        entered_cnv_type = cols6[1].multiselect("CNV_Type", cnv_type)
        entered_acmg_class = cols6[2].multiselect("ACMG_Class", acmg_class)

    if entered_tsv_file:
        tsv_df = pd.read_csv(entered_tsv_file, delimiter="\t")
        with open(annotsv_format_path, "r") as column_file:
            columns_to_keep = [line.strip() for line in column_file]
        tsv_df = tsv_df[columns_to_keep]
        tsv_df["SV_chrom"] = pd.Categorical(tsv_df["SV_chrom"], categories=[str(i) for i in range(1, 23)] + ["X", "Y"])
        tsv_df.sort_values("SV_chrom", inplace=True)
        tsv_df.reset_index(drop=True, inplace=True)
        filtered_tsv = helpers.filter_tsv(tsv_df, chromosome_list_cnv, cnv_type, acmg_class, entered_cnv_chrom, entered_cnv_type, entered_acmg_class)
        st.write(filtered_tsv)

        if st.button("Prepare for download of annotated tsv data"):
            table_exporter = exporter.CNVExporter()
            to_be_exported = table_exporter.save_tables_as_excel_tsv(filtered_tsv)
            st.download_button(label="Download", data=to_be_exported, file_name=f"{sample_name}_annotated_df.xlsx")
        else:
            st.write("Click button to prepare download")

if __name__=='__main__':
    main()