"""
Tests for genome fetching and taxonomy parsing functionality
"""

import os
import pytest
import pandas as pd
from pathlib import Path
from samovar.genome_fetcher import fetch_genome, parse_taxonomy_table

@pytest.fixture
def test_output_dir():
    """Create and return test output directory"""
    output_dir = Path("tests_outs/test_genome_fetcher")
    output_dir.mkdir(parents=True, exist_ok=True)
    return str(output_dir)

@pytest.fixture
def test_taxonomy_file(test_output_dir):
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

def test_parse_taxonomy_table(test_taxonomy_file, test_output_dir):
    """Test parsing taxonomy table and downloading genomes"""
    email = "test@example.com"
    results = parse_taxonomy_table(test_taxonomy_file, test_output_dir, email, reference_only=False)
    
    # Should have downloaded at least one genome
    assert len(results) > 0
    
    # Check that all files exist
    for result in results:
        assert os.path.exists(result)
        assert result.endswith('.fa')

def test_parse_taxonomy_table_no_taxid_columns(test_output_dir):
    """Test parsing taxonomy table with no taxid columns"""
    # Create a file with no taxid columns
    data = {'col1': ['A', 'B', 'C'], 'col2': ['D', 'E', 'F']}
    df = pd.DataFrame(data)
    file_path = os.path.join(test_output_dir, 'no_taxid.tsv')
    df.to_csv(file_path, sep='\t', index=False)
    
    email = "test@example.com"
    results = parse_taxonomy_table(file_path, test_output_dir, email)
    assert len(results) == 0

def test_parse_taxonomy_table_invalid_file(test_output_dir):
    """Test parsing non-existent taxonomy file"""
    email = "test@example.com"
    results = parse_taxonomy_table("nonexistent.tsv", test_output_dir, email)
    assert len(results) == 0 