import streamlit as st
import pandas as pd
import os
from pathlib import Path
import dotenv
import argparse
from cnvizard import make_pretty, CNVExporter, CNVPlotter, filter_tsv, CNVVisualizer, prepare_cnv_table, explode_cnv_table

def validate_env_file(env_file):
    """
    Validates the contents of an environment file.

    Args:
        env_file (str): Path to the environment file.

    Returns:
        tuple: A tuple containing two lists:
            - missing_keys (list): List of missing keys.
            - invalid_paths (list): List of invalid paths.
    """
    required_keys = [
        "OMIM_ANNOTATION_PATH",
        "CANDIDATE_LIST_DIR",
        "REFERENCE_FILES_DIR",
        "ANNOTS_SV_FORMAT_PATH"
    ]
    optional_keys = [
        "APPSETTING_IGV_OUTLINK"
    ]
    missing_keys = []
    invalid_paths = []

    with open(env_file, "r") as file:
        lines = file.readlines()
        env_dict = dict(line.strip().split('=') for line in lines if '=' in line)

    for key in required_keys:
        if key not in env_dict:
            missing_keys.append(key)
        elif not Path(env_dict[key]).exists():
            invalid_paths.append(env_dict[key])

    for key in optional_keys:
        if key in env_dict and not Path(env_dict[key]).exists():
            invalid_paths.append(env_dict[key])

    return missing_keys, invalid_paths

def load_and_select_env():
    """
    Streamlit function to load or create an environment file.

    This function allows users to upload an existing environment file or create a new one by providing necessary paths and settings.

    Returns:
        str: Path to the loaded or created environment file, or None if there was an error.
    """
    st.title("CNVizard Environment Selector")

    st.markdown("""
    ### Select or Create an Environment File
    """)

    uploaded_env_file = st.file_uploader("Upload environment file")

    if uploaded_env_file:
        try:
            env_file = uploaded_env_file.name
            with open(env_file, "wb") as file:
                file.write(uploaded_env_file.getbuffer())

            missing_keys, invalid_paths = validate_env_file(env_file)
            if missing_keys or invalid_paths:
                error_message = "Error in environment file:\n"
                if missing_keys:
                    error_message += f"Missing keys: {', '.join(missing_keys)}\n"
                if invalid_paths:
                    error_message += f"Invalid paths: {', '.join(invalid_paths)}"
                raise ValueError(error_message)

            dotenv.load_dotenv(env_file)
            st.success(f"Loaded environment file: {env_file}")
            return env_file
        except Exception as e:
            st.error(f"Error loading environment file: {e}")
            return None

    if st.button("Create new `.env` file"):
        st.session_state.create_env = True

    if st.session_state.get("create_env", False):
        st.markdown("#### Provide the paths for the new `.env` file")
        env_output_path = st.text_input("Output path for new environment file:", "default.env")
        
        omim_annotation_path = st.text_input("OMIM annotation path:", "./resources/omim.txt")
        candidate_list_dir = st.text_input("Candidate list directory:", "./resources/candidate_lists")
        reference_files_dir = st.text_input("Reference files directory:", "./resources/references")
        annotsv_format_path = st.text_input("Annotsv format path:", "./resources/annotsv_format.txt")
        igv_outlink = st.text_input("IGV outlink (optional):", "")

        if st.button("Save new `.env` file"):
            try:
                paths = {
                    "OMIM_ANNOTATION_PATH": omim_annotation_path,
                    "CANDIDATE_LIST_DIR": candidate_list_dir,
                    "REFERENCE_FILES_DIR": reference_files_dir,
                    "ANNOTS_SV_FORMAT_PATH": annotsv_format_path,
                }
                if igv_outlink:
                    paths["APPSETTING_IGV_OUTLINK"] = igv_outlink

                invalid_paths = [path for path in paths.values() if not Path(path).exists()]
                if invalid_paths:
                    raise ValueError(f"Invalid paths: {', '.join(invalid_paths)}")

                env_file = env_output_path if env_output_path else ".env"
                with open(env_file, "w") as file:
                    for key, value in paths.items():
                        file.write(f"{key}={value}\n")
                st.success(f"New `.env` file created: `{env_file}`")
                st.session_state.create_env = False
                return env_file
            except Exception as e:
                st.error(f"Error saving environment file: {e}")
                return None

    return None

def prepare_filter_for_consecutive_cnvs(self, df: pd.DataFrame) -> pd.DataFrame:
    """
    Function which is used to filter for consecutively deleted/duplicated exons.

    Args:
        df (pd.DataFrame): .cnr DataFrame

    Returns:
        pd.DataFrame: .cnr DataFrame filtered for consecutively deleted/duplicated exons.
    """
    df = df.copy()  # Create a copy to avoid SettingWithCopyWarning
    df['difference_previous'] = df.groupby('gene')['exon'].diff()
    df['difference_previous'] = df['difference_previous'].fillna(method='bfill')
    df['difference_next'] = df.groupby('gene')['exon'].diff(periods=-1)
    df['difference_next'] = df['difference_next'].fillna(method='ffill')
    return df


def main(env_file_path):
    if env_file_path is None:
        env_file_path = load_and_select_env()
        if not env_file_path:
            st.stop()

    # Load environment variables
    dotenv.load_dotenv(env_file_path)
    igv_string = os.getenv("APPSETTING_IGV_OUTLINK")

    # Load paths from environment variables
    omim_annotation_path = Path(os.getenv("OMIM_ANNOTATION_PATH", "./resources/omim/omim.txt"))
    candidate_list_dir = Path(os.getenv("CANDIDATE_LIST_DIR", "./resources/candidate_lists"))
    reference_files_dir = Path(os.getenv("REFERENCE_FILES_DIR", "./resources/references"))
    annotsv_format_path = Path(os.getenv("ANNOTS_SV_FORMAT_PATH", "./resources/annotsv_format.txt"))

    # Filter options
    chrom_list = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    call_list = [0, 1, 2, 3]

    # Streamlit App Title
    st.title("CNVizard")
    st.markdown("This is a Streamlit web app providing analysis tools for genetic copy number variants.")

    # File Uploads Section
    st.subheader("Upload")
    st.markdown("Please upload a reference file, an individual .cnr file, and an individual .bintest file, provided by CNVkit.")

    omim_annotation_file = None
    gene_list = []

    if omim_annotation_path.exists():
        omim_annotation_file = omim_annotation_path
        gene_list = pd.read_csv(omim_annotation_file, sep="\t")["gene"].tolist()

    cols = st.columns(3)
    entered_cnr = cols[0].file_uploader(".cnr", type=["txt", "cnr"])
    entered_bintest = cols[1].file_uploader("bintest", type=["txt", "tsv"])
    ngs_type = cols[2].radio("Choose NGS Type", ["WES", "WGS"])

    reference_path = reference_files_dir / ("genome_cnv_reference_large.parquet" if ngs_type == "WGS" else "exome_cnv_reference_large.parquet")
    reference_bintest_path = reference_files_dir / ("genome_cnv_reference_bintest_large.parquet" if ngs_type == "WGS" else "exome_cnv_reference_bintest_large.parquet")

    sample_name = ""
    reference_df = None
    cnr_df = None
    bintest_df = None
    if entered_cnr:
        sample_name = entered_cnr.name.split(".")[0]
        try:
            cnr_df = pd.read_csv(entered_cnr, delimiter="\t")
        except Exception as e:
            st.error(f"Error reading .cnr file: {e}")
            cnr_df = None
        if igv_string:
            igv_string = igv_string.replace("samplename", sample_name)
    if entered_bintest:
        try:
            bintest_df = pd.read_csv(entered_bintest, delimiter="\t")
        except Exception as e:
            st.error(f"Error reading bintest file: {e}")
            bintest_df = None
    if reference_path.exists():
        try:
            reference_df = pd.read_parquet(reference_path)
        except Exception as e:
            st.error(f"Error reading reference file: {e}")
            reference_df = None
    try:
        reference_bintest_df = pd.read_parquet(reference_bintest_path) if reference_bintest_path.exists() else pd.DataFrame()
    except Exception as e:
        st.error(f"Error reading reference bintest file: {e}")
        reference_bintest_df = pd.DataFrame()

    st.subheader("Configurations")
    st.markdown("Provide a sample name and define the number of consecutive exons to display (default value = 2).")

    cols2 = st.columns(2)
    entered_del_size = cols2[0].text_input("deletion_size", value="2")
    entered_dup_size = cols2[1].text_input("duplication_size", value="2")

    # Sidebar configuration
    st.sidebar.title("About")
    st.sidebar.markdown("This Streamlit web app enables you to visualize copy-number-variant data using dataframes and plots.")

    st.sidebar.subheader("Candidate Gene Selection")
    candidate_options = sorted(os.listdir(candidate_list_dir)) if candidate_list_dir.exists() else []
    selected_candidate = st.sidebar.radio("Select the desired candidate gene list", candidate_options)
    st.sidebar.write(selected_candidate)

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

    list_of_possible_dataframes = ["total", "bintest", "hom_del", "total_candidate", "bintest_candidate", "consecutive_del", "consecutive_dup"]
    df_to_be_displayed = st.selectbox("Select dataframe to display", list_of_possible_dataframes)

    if reference_df is not None and cnr_df is not None and bintest_df is not None:
        cnv_visualizer_instance = CNVVisualizer(reference_df, cnr_df, bintest_df)
        omim_df, candidate_df, cnr_db, bintest_db = cnv_visualizer_instance.format_df(
            omim_annotation_file, os.path.join(candidate_list_dir, selected_candidate) if candidate_list_dir.exists() else None
        )

        call_df = reference_df[["gene", "exon", "het_del_frequency", "hom_del_frequency", "dup_frequency"]]
        bintest_inhouse_df = reference_bintest_df if not reference_bintest_df.empty else pd.DataFrame()

        cnr_db = pd.merge(cnr_db, call_df, on=["gene", "exon"], how="left")
        if not bintest_inhouse_df.empty:
            bintest_db = pd.merge(bintest_db, bintest_inhouse_df, on=["gene", "exon"], how="left")
            bintest_db = bintest_db.infer_objects().fillna(0)

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
            "total_candidate": cnv_visualizer_instance.filter_for_candi_cnvs(cnr_db_filtered, candidate_df),
            "bintest_candidate": cnv_visualizer_instance.filter_for_candi_cnvs(bintest_db, candidate_df),
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
        st.dataframe(display_mapping[df_to_be_displayed].style.pipe(make_pretty) if df_to_be_displayed != "total" else display_mapping[df_to_be_displayed])

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
                cnv_visualizer_instance.filter_for_candi_cnvs(cnr_db, candidate_df),
                cnv_visualizer_instance.filter_for_candi_cnvs(bintest_db, candidate_df),
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
            table_exporter = CNVExporter()
            export_data = table_exporter.save_tables_as_excel(*tables_to_export)
            download_button_columns[0].download_button(label="Download", data=export_data, file_name=f"{sample_name}_df.xlsx")
        else:
            download_message_columns[0].write("Click button to prepare download")

        if download_preparator_filtered:
            table_exporter = CNVExporter()
            export_data_filtered = table_exporter.save_filtered_table_as_excel(download_filter, df_to_be_displayed)
            download_button_columns[1].download_button(label="Download filtered", data=export_data_filtered, file_name=f"{sample_name}_filtered_df.xlsx")
        else:
            download_message_columns[1].write("Click button to prepare download")

    st.subheader("Plot log2 and depth for selected genes")
    entered_gene = st.multiselect("gene", gene_list, max_selections=1)
    entered_gene = entered_gene[0].upper() if entered_gene else None

    if reference_df is not None and cnr_df is not None and bintest_df is not None and entered_gene:
        gene_plotter = CNVPlotter()
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
        try:
            father_cnr_df = cnv_visualizer_instance.prepare_parent_cnv(pd.read_csv(father_cnr, delimiter="\t"))
            mother_cnr_df = cnv_visualizer_instance.prepare_parent_cnv(pd.read_csv(mother_cnr, delimiter="\t"))
        except Exception as e:
            st.error(f"Error reading parent .cnr files: {e}")
            st.stop()

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
        try:
            tsv_df = pd.read_csv(entered_tsv_file, delimiter="\t")
        except Exception as e:
            st.error(f"Error reading .tsv file: {e}")
            st.stop()

        with open(annotsv_format_path, "r") as column_file:
            columns_to_keep = [line.strip() for line in column_file]

        tsv_df = tsv_df[columns_to_keep]

        # Fix ValueErrors and ensure ACMG_class is treated as int
        tsv_df["ACMG_class"] = tsv_df["ACMG_class"].apply(
            lambda x: int(x.split('=')[-1]) if isinstance(x, str) and 'full=' in x else (
                int(x) if isinstance(x, str) and x.isdigit() else -1
            )
        )
        # Remove rows with -1 in ACMG_class after conversion
        tsv_df = tsv_df[tsv_df["ACMG_class"] != -1]

        # Ensure ACMG_class is treated as int after conversion
        tsv_df["ACMG_class"] = tsv_df["ACMG_class"].astype(int)

        # Convert AnnotSV_ranking_score to numeric, handle errors by converting invalid values to NaN
        tsv_df["AnnotSV_ranking_score"] = tsv_df["AnnotSV_ranking_score"].apply(
            lambda x: pd.to_numeric(x, errors='coerce') if isinstance(x, str) else float('nan') if x == '.' else x
        )

        filtered_tsv = filter_tsv(tsv_df,chromosome_list_cnv,cnv_type,acmg_class,entered_cnv_chrom,entered_cnv_type,entered_acmg_class)
        filtered_tsv = filtered_tsv.fillna(".")
        st.write("Filtered AnnotSV DataFrame:")
        st.write(filtered_tsv)

        if st.button("Prepare for download of annotated tsv data"):
            table_exporter = CNVExporter()
            to_be_exported = table_exporter.save_tables_as_excel_tsv(filtered_tsv)
            st.download_button(label="Download", data=to_be_exported, file_name=f"{sample_name}_annotated_df.xlsx")
        else:
            st.write("Click button to prepare download")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run CNVizard Streamlit app.")
    parser.add_argument("env", nargs='?', default=None, help="Path to the .env file.")
    args = parser.parse_args()

    st.set_page_config(layout="wide", page_title="CNVizard", page_icon="CNVizard.png")

    if args.env:
        main(args.env)
    else:
        env_file_path = load_and_select_env()
        if env_file_path:
            main(env_file_path)
