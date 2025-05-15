"""
Core functionality for the samovaR Python wrapper.
"""

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Initialize R package
try:
    samovaR = importr("samovaR")
except Exception as e:
    raise ImportError(
        "Failed to import samovaR R package. Please ensure it is installed in R."
    ) from e

# Enable pandas conversion
pandas2ri.activate()

def __getattr__(name):
    """Forward attribute access to the R package."""
    try:
        return getattr(samovaR, name)
    except AttributeError:
        raise AttributeError(f"'{name}' not found in samovaR package") 