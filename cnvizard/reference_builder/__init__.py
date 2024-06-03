"""
CNVizard Reference Builder Package

This package contains modules for creating and merging reference files used by CNVizard.
"""

# Import functions and classes from the main modules
from .create_reference_files import main as create_reference_files
from .merge_reference_files import main as merge_reference_files

# Import utilities
from .data_processing import (
    prepare_cnv_table,
    explode_cnv_table,
    _get_call_counts,
    _get_frequency,
)
