"""
File containing pandas stylers used for the CNVizard.
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import pandas as pd

def make_pretty(styler):
    """
    Styler function used to highlight and format the filtered .cnr/bintest DataFrame.

    Args:
        styler (pd.io.formats.style.Styler): Styler object to format.

    Returns:
        pd.io.formats.style.Styler: Formatted Styler object.
    """
    styler.format(precision=2)
    styler.format(subset="OMIMG", precision=0)
    styler.format(subset=["start", "end"], thousands=',', decimal='.')
    styler.highlight_between(subset='weight', color='green', axis=0, left=0.9, right=1)
    styler.highlight_between(subset='weight', color='yellow', axis=0, left=0.8, right=0.9)
    styler.highlight_between(subset='weight', color='red', axis=0, left=0, right=0.8)
    styler.highlight_between(subset='depth', color='red', axis=0, right=0)
    styler.highlight_between(subset='log2', axis=0, right=-0.65)
    return styler

def mark_log2(value):
    if value <= float(0.65):
        return "background-color: pink;"
