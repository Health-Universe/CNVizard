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
- [Application Structure](#application-structure)
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

5. **Start the application with optinally giving a path to an environment file:**
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
   - Navigate to `/CNVizard/resources/annotsv_format.txt`.
   - Modify the standard column names as needed.

3. **Add/Remove a panel-list:**
   - Navigate to `/CNVizard/resources/candidate_lists/`.
   - Add or remove new panel `.txt` files.

4. **Modify OMIM List:**
   - Navigate to `/CNVizard/resources/omim.txt`.
   - Modify the `.txt` file as needed.

5. **Change Reference List:**
   - Navigate to `/CNVizard/resources/references/`.
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

## Application Structure

```bash
CNVizard/
├── cnvizard/
│   ├── app.py                     # Main application script
│   ├── exporter.py                # Handles exporting data
│   ├── helpers.py                 # Helper functions for filtering
│   ├── __init__.py                # Package initializer
│   ├── plotter.py                 # Plotting functions
│   ├── reference_processing.py    # Functions for processing reference files
│   ├── styler.py                  # Styling functions for dataframes
│   └── visualizer.py              # Main visualizer class
├── CNVizard.egg-info/
├── CNVizard.png                   # Application logo
├── default.env                    # Example environment file
├── LICENSE                        # License file
├── MANIFEST.in                    # Manifest file for packaging
├── README.md                      # This README file
├── requirements.txt               # Required dependencies
├── resources/
│   ├── annotsv_format.txt         # AnnotSV format file
│   ├── candidate_lists/           # Candidate gene lists
│   │   ├── cardiac_rythm_disorders.txt
│   │   ├── ...                    # Other candidate lists
│   ├── omim.txt                   # OMIM annotation file
│   └── references/                # Reference files
│       ├── bintest_ref.parquet
│       ├── ...                    # Other reference files
└── setup.py                       # Setup script for packaging
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.