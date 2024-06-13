# CNVizard

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Streamlit](https://img.shields.io/badge/Streamlit-0.84.2-brightgreen.svg)](https://streamlit.io/)
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

1. **Create a virtual environment and activate it:**
   ```sh
   python3 -m venv venv
   source venv/bin/activate
   ```

2. **Navigate into the CNVizard folder:**
   ```sh
   cd CNVizard
   ```

3. **Install the required dependencies:**
   ```sh
   pip install -r requirements.txt
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

### Prepare CNV Table
The `prepare_cnv_table` function formats and reorders the `.cnr` DataFrame.

```python
import pandas as pd
from cnvizard.reference_processing import prepare_cnv_table

# Load your data
df = pd.read_csv('path_to_cnr_file.cnr', delimiter='\t')
df2 = pd.read_csv('path_to_omim_file.txt', delimiter='\t')

# Prepare the CNV table
prepared_df = prepare_cnv_table(df, df2)
```

### Explode CNV Table
The `explode_cnv_table` function splits redundant exons and explodes the DataFrame.

```python
from cnvizard.reference_processing import explode_cnv_table

# Explode the CNV table
exploded_df = explode_cnv_table(prepared_df)
```

### Merge Reference Files
The `merge_reference_files` function merges previously created individual reference files into a single reference file.

```python
from cnvizard.reference_processing import merge_reference_files

# Define paths
path_to_input = 'path_to_individual_references'
path_to_output = 'path_to_output_directory'
path_to_bintest = 'path_to_bintest_references'

# Merge reference files
merge_reference_files(path_to_input, path_to_output, path_to_bintest)
```

### Create Reference Files
The `create_reference_files` function creates individual reference files for CNV visualization.

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
```

These functions will help you process and create references for your CNV analysis using CNVizard. Make sure to adjust the paths and parameters according to your specific setup and requirements.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
```

This section provides clear instructions on how to use the provided utility functions to create and process reference files for CNVizard. Ensure to adjust the paths and parameters according to your specific setup and requirements.
