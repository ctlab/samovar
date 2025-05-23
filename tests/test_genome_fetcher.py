"""
Tests for genome fetching and taxonomy parsing functionality
"""

import os
import pytest
import pandas as pd
from pathlib import Path
from samovar.genome_fetcher import fetch_genome

@pytest.fixture
def test_output_dir():
    """Create and return test output directory"""
    output_dir = Path("tests_outs/test_genome_fetcher")
    output_dir.mkdir(parents=True, exist_ok=True)
    return str(output_dir)

@pytest.fixture
def test_annotation_table(test_output_dir):
    """Create a test taxonomy file"""
    data = {
        'species': ['E. coli', 'S. aureus', 'B. subtilis'],
        'taxid': ['562', '1280', '1423'],
        'other_taxid': ['562', '1280', '1423'],
        'non_taxid': ['A', 'B', 'C']
    }
    df = pd.DataFrame(data)
    file_path = os.path.join(test_output_dir, 'test_taxonomy.tsv')
    df.to_csv(file_path, sep='\t', index=False)
    return file_path

def test_fetch_genome_existing(test_output_dir):
    """Test fetching genome for a known taxid"""
    # Use E. coli as a test case (taxid: 562)
    email = "test@example.com"
    result = fetch_genome("562", test_output_dir, email, reference_only=False)  # Try without reference filter
    
    assert result is not None
    assert os.path.exists(result)
    assert result.endswith("562.fa")

def test_fetch_genome_invalid_taxid(test_output_dir):
    """Test fetching genome for an invalid taxid"""
    email = "test@example.com"
    result = fetch_genome("999999999", test_output_dir, email)
    assert result is None

def test_fetch_genome_already_exists(test_output_dir):
    """Test that existing genome is not re-downloaded"""
    # First download
    email = "test@example.com"
    result1 = fetch_genome("562", test_output_dir, email, reference_only=False)
    
    if result1 is None:
        pytest.skip("Could not download initial genome for testing")
    
    # Get file modification time
    mtime1 = os.path.getmtime(result1)
    
    # Try to download again
    result2 = fetch_genome("562", test_output_dir, email, reference_only=False)
    mtime2 = os.path.getmtime(result2)
    
    assert result1 == result2
    assert mtime1 == mtime2  # File should not have been modified