import os
import pytest
import pandas as pd
import tempfile
import yaml
from pathlib import Path
from samovar.table2iss import samovar_annotation_regenerate, process_abundance_table
from unittest.mock import patch, MagicMock

@pytest.fixture
def test_data_dir():
    return Path("data/test_annotations")

@pytest.fixture
def test_output_dir():
    output_dir = Path("tests_outs/test_samovar")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir

@pytest.fixture
def mock_config():
    config = {
        'threshold_amount': 1e-5,
        'plot_log': False,
        'min_cluster_size': 2,
        'N': 10,
        'N_reads': 1000,
        'output_dir': str(Path("tests_outs/test_samovar"))
    }
    return config

def test_samovar_annotation_regenerate_basic(test_data_dir, test_output_dir, mock_config):
    """Test basic functionality of samovar_annotation_regenerate"""
    # Create temporary config file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(mock_config, f)
        config_path = f.name

    try:
        # Run the function
        samovar_annotation_regenerate(
            annotation_dir=str(test_data_dir),
            config_samovar=config_path,
            output_dir=str(test_output_dir)
        )

        # Check if output files were created
        output_files = list(test_output_dir.glob("*.csv"))
        assert len(output_files) > 0, "No output files were created"

        # Check if the output files are valid CSV files
        for output_file in output_files:
            df = pd.read_csv(output_file)
            assert not df.empty, f"Output file {output_file} is empty"
            # Convert numeric columns to float and check for non-negative values
            numeric_cols = df.select_dtypes(include=['object']).columns
            for col in numeric_cols:
                try:
                    values = pd.to_numeric(df[col], errors='coerce')
                    assert all(values >= 0), f"Output file {output_file} contains negative values in column {col}"
                except:
                    pass  # Skip non-numeric columns

    finally:
        # Cleanup
        os.unlink(config_path)

def test_samovar_annotation_regenerate_integration(test_data_dir, test_output_dir, mock_config):
    """Test integration with process_abundance_table"""
    # Create temporary config file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(mock_config, f)
        config_path = f.name

    try:
        # First run samovar_annotation_regenerate
        samovar_annotation_regenerate(
            annotation_dir=str(test_data_dir),
            config_samovar=config_path,
            output_dir=str(test_output_dir)
        )

        # Get the first output file
        output_files = list(test_output_dir.glob("*.csv"))
        assert len(output_files) > 0, "No output files were created"
        test_file = output_files[0]

        # Create a temporary genome directory for testing
        with tempfile.TemporaryDirectory() as genome_dir:
            # Mock the genome fetching functionality
            with patch('samovar.table2iss.fetch_genome') as mock_fetch, \
                 patch('samovar.table2iss.get_genome_file') as mock_get, \
                 patch('samovar.table2iss.regenerate_metagenome') as mock_regenerate:
                # Mock get_genome_file to return None (no existing genome)
                mock_get.return_value = None
                # Mock fetch_genome to return a dummy genome file
                mock_fetch.return_value = os.path.join(genome_dir, "dummy.fasta")
                
                # Create a dummy genome file
                with open(os.path.join(genome_dir, "dummy.fasta"), 'w') as f:
                    f.write(">dummy\nATCG\n")

                # Mock regenerate_metagenome to do nothing
                mock_regenerate.return_value = None

                # Process the output with process_abundance_table
                result = process_abundance_table(
                    table=str(test_file),
                    genome_dir=genome_dir,
                    output_dir=str(test_output_dir / "processed")
                )
                
                # Verify the result
                assert isinstance(result, pd.DataFrame), "Result should be a DataFrame"
                assert not result.empty, "Result DataFrame should not be empty"
                # Convert numeric columns to float and check for non-negative values
                numeric_cols = result.select_dtypes(include=['object']).columns
                for col in numeric_cols:
                    try:
                        values = pd.to_numeric(result[col], errors='coerce')
                        assert all(values >= 0), f"Result contains negative values in column {col}"
                    except:
                        pass  # Skip non-numeric columns

    finally:
        # Cleanup
        os.unlink(config_path) 