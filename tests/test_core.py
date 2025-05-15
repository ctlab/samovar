"""
Tests for the core functionality of the samovaR Python wrapper.
"""

import pytest
from samovaR import core

def test_samovaR_import():
    """Test that the R package can be imported."""
    assert core.samovaR is not None

def test_pandas_conversion():
    """Test that pandas conversion is enabled."""
    import pandas as pd
    import numpy as np
    
    # Create a simple pandas DataFrame
    df = pd.DataFrame({
        'A': [1, 2, 3],
        'B': ['a', 'b', 'c']
    })
    
    # Convert to R and back
    r_df = core.pandas2ri.py2rpy(df)
    py_df = core.pandas2ri.rpy2py(r_df)
    
    # Check that the data is preserved
    assert isinstance(py_df, pd.DataFrame)
    assert py_df.equals(df) 