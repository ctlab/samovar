import os
import pandas as pd
from pathlib import Path
from samovar.table2iss import parse_annotation_table
import pytest

def test_parse_annotation_table():
    """Test reading taxonomy table and counting occurrences."""
    # Create output directory
    output_dir = "tests_outs/test_parse_annotation_table"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test table with the example data
    table_path = os.path.join(output_dir, "test_table.csv")
    df = pd.DataFrame({
        "seqID": ["a1", "a2", "a3"],
        "taxID_k1_0": ["0", "562", "0"],
        "taxID_k2_1": ["0", "562", "562"]
    })
    df.to_csv(table_path, index=False)
    
    # Read table
    result = parse_annotation_table(table_path)
    
    # Save output for inspection
    result.to_csv(os.path.join(output_dir, "result.csv"), index=False)
    
    # Check results
    assert isinstance(result, pd.DataFrame)
    assert set(result.columns) == {'taxid', 'N_k1', 'N_k2'}
    assert len(result) == 2  # 0 and 562
    
    # Check counts
    taxid_0 = result[result['taxid'] == '0'].iloc[0]
    assert taxid_0['N_k1'] == 2
    assert taxid_0['N_k2'] == 1 
    
    taxid_562 = result[result['taxid'] == '562'].iloc[0]
    assert taxid_562['N_k1'] == 1
    assert taxid_562['N_k2'] == 2 

def test_parse_real_annotation_table():
    """Test parsing the real annotation table and check output dimensions."""
    # Create output directory
    output_dir = "tests_outs/test_parse_real_annotation_table"
    os.makedirs(output_dir, exist_ok=True)
    
    # Path to real annotation file
    table_path = "tests/data/annotation.csv"
    
    # Parse table
    result = parse_annotation_table(table_path)
    
    # Save output for manual inspection
    output_path = os.path.join(output_dir, "result.csv")
    result.to_csv(output_path, index=False)
    
    # Check basic properties
    assert isinstance(result, pd.DataFrame)
    assert 'taxid' in result.columns
    assert all(col.startswith('N_') for col in result.columns if col != 'taxid')
    assert len(result) > 0  # Should have at least one taxid 