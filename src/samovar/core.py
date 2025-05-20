"""
Core functionality for the samovaR Python wrapper.
"""

import json
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

def load_config():
    """Load configuration from config.json."""
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config.json')
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return {"r_lib_path": None, "r_version": None}

# Load configuration
config = load_config()

# Set R library path if specified
if config["r_lib_path"]:
    robjects.r(f'.libPaths("{config["r_lib_path"]}")')

# Initialize R package
try:
    samovaR = importr("samovaR")
except Exception as e:
    raise ImportError(
        f"Failed to import samovaR R package. Please ensure it is installed in R. "
        f"Using R library path: {config['r_lib_path']}"
    ) from e

# Enable pandas conversion
pandas2ri.activate()

def __getattr__(name):
    """Forward attribute access to the R package."""
    try:
        return getattr(samovaR, name)
    except AttributeError:
        raise AttributeError(f"'{name}' not found in samovaR package") 