# CNVizard

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.35.0-brightgreen.svg)](https://streamlit.io/)
[![Python](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)

CNVizard is a Streamlit-based application designed for the visualization and analysis of germline copy number variants (CNVs). This tool provides comprehensive features to help researchers and clinicians analyze genetic data with ease.

## Table of Contents
- [Setup](#setup)
  - [Mandatory Setup](#mandatory-setup)
  - [Optional Setup](#optional-setup)
- [Usage](#usage)
- [Creating References](#creating-references)
- [License](#license)

## Setup

### Mandatory Setup
1. **Clone this repository and change directory:**
   ```sh
   git clone https://github.com/IHGGM-Aachen/CNVizard
   cd CNVizard
   ```

2. **Create a virtual environment and activate it:**
   ```sh
   python3 -m venv venv
   source venv/bin/activate
   ```

4. **Install the CNVizard:**
   ```sh
   pip install -e .
   ```

5. **Start the application with optionally giving a path to an environment file:**
   ```sh
   streamlit run cnvizard/app.py [ENV_FILE]
   ```

### Optional Setup

1. **Create an `.env` file in the CNVizard folder (optional):**
   - Can also be used to set up IGV outlinks.
   - Example content:
     ```
     APPSETTING_IGV_OUTLINK=/path/to/samplename/samplename.cram
     OMIM_ANNOTATION_PATH=./resources/omim/omim.txt
     CANDIDATE_LIST_DIR=./resources/candidate_lists
     REFERENCE_FILES_DIR=./resources/references
     ANNOTSV_FORMAT_PATH=./resources/annotsv_format.txt
     ```

2. **Change the AnnotSV table format:**
   - Navigate to `resources/annotsv_format.txt`.
   - Modify the standard column names as needed.

3. **Add/Remove a panel-list:**
   - Navigate to `resources/candidate_lists/`.
   - Add or remove new panel `.txt` files.

4. **Modify OMIM List:**
   - Navigate to `resources/omim.txt`.
   - Modify the `.txt` file as needed.

5. **Change Reference List:**
   - Navigate to `resources/references/`.
   - Replace `bintest_ref.parquet`, `total_ref.parquet`, or other reference files as needed.

## Usage

1. **Run the application:**
   ```sh
   streamlit run cnvizard/app.py
   ```

2. **Select or upload your `.env` file:**
   - If no `.env` file is specified via the command line, you will be prompted to upload or create a new `.env` file within the application.

3. **Upload necessary files:**
   - Upload the reference file, individual `.cnr` file, and `.bintest` file provided by CNVkit.
   - Optionally, upload additional files for trio analysis or scatter plots.

4. **Configure settings:**
   - Provide sample name and define the number of consecutive exons to display.
   - Select candidate gene lists and apply filters as needed.

5. **Visualize and analyze data:**
   - View filtered dataframes and plots.
   - Download prepared data in Excel format for further analysis.

## Creating References

To create new references using the CNVizard utility functions, follow these steps:

### Create Reference Files using provided aggregated functions  
The `create_reference_files` function creates individual reference files for CNV visualization.
This function all necessary stepts until the merge step.
This function requires a predefined directory structure.
Each dir resembles an individual function.
For exome-data : .../results/CNV/...cnr & .../results/CNV/...bintest
For genome-data : .../exome_extract/CNV/...cnr & .../exome_extract/CNV/...bintest

```python
from cnvizard.reference_processing import create_reference_files

# Define parameters
path_to_input = 'path_to_input_directory'
ngs_type = 'WES'  # or 'WGS'
path_to_output = 'path_to_output_directory'
omim_path = 'path_to_omim_directory'
reference_type = 'normal'  # or 'bintest'
starting_letter = 'A'  # filter subdirectories starting with this letter

# Create reference files
create_reference_files(path_to_input, ngs_type, path_to_output, omim_path, reference_type, starting_letter)

# Define paths
path_to_input = 'path_to_individual_references'
path_to_output = 'path_to_output_directory'
path_to_bintest = 'path_to_bintest_references'

# Merge reference files and precalculate statistics
merge_reference_files(path_to_input, path_to_output, path_to_bintest)
```
### Explanation of indivual steps : 
In a first step the indivual .cnr or bintest files are loaded into pandas dataframes which are stored in a list 
```python
import pandas as pd
from cnvizard.reference_processing import prepare_cnv_table, explode_cnv_table, merge_reference_files

# Load your data / Repeat for each patient 
list_to_collect_dataframes = []
df = pd.read_csv('path_to_cnr_file.cnr', delimiter='\t')
list_to_collect_dataframes.append(df)

#Concat all samples to a large dataframe
reference = pd.concat(list_of_dfs)

#The `explode_cnv_table` function splits redundant exons and explodes the DataFrame.
exploded_df = explode_cnv_table(prepared_df)

# Prepare the CNV table
#Read omim dataframe
df2 = pd.read_csv('path_to_omim_file.txt', delimiter='\t')
#The `prepare_cnv_table` function formats and reorders the `.cnr` reference DataFrame.
prepared_df = prepare_cnv_table(df, df2)

# Export reference to parquet file
ordered_reference.to_parquet(export_name_parquet, index=False)

```

### Merge Reference Files
The `merge_reference_files` function merges previously created multiple reference files into a single reference file.
Additionally precalculates statistics required for the generation of internal frequencies and plotting.
If all samples are have been used to create one large reference dataframe, this function only calculates the necessary statistics

```python
from cnvizard.reference_processing import merge_reference_files

# Define paths
path_to_input = 'path_to_individual_references'
path_to_output = 'path_to_output_directory'
path_to_bintest = 'path_to_bintest_references'

# Merge reference files
merge_reference_files(path_to_input, path_to_output, path_to_bintest)
```

These functions will help you process and create references for your CNV analysis using CNVizard. Make sure to adjust the paths and parameters according to your specific setup and requirements.

## Convert genomics england panel app files to compatible gene lists
The `convert_genomics_england_panel_to_txt` function transfroms gene-panel files from the genomics england panel app 
to gene-list files compatible with the CNVizard.

```python
from cnvizard.reference_processing import convert_genomics_england_panel_to_txt

# Define paths
path_to_input = 'path_to_panel_app_file'
path_to_output = 'path_to_output_file'

# Merge reference files
convert_genomics_england_panel_to_txt(path_to_input, path_to_output)
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
```

This section provides clear instructions on how to use the provided utility functions to create and process reference files for CNVizard. Ensure to adjust the paths and parameters according to your specific setup and requirements.
