"""
Tests for genome fetching and taxonomy parsing functionality
"""

import gzip
import os
import pytest
import pandas as pd
from pathlib import Path
from samovar.genome_fetcher import bundled_test_genome, fetch_genome


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


def test_bundled_test_genome_has_shipped_taxids():
    assert bundled_test_genome("562") is not None
    assert bundled_test_genome("9606") is not None
    assert bundled_test_genome("999999999") is None


def test_fetch_genome_reuses_gzipped_processed(tmp_path):
    dest = tmp_path / "562-processed.fasta.gz"
    with gzip.open(dest, "wt") as handle:
        handle.write(">x\nACGT\n")
    result = fetch_genome("562", str(tmp_path), "test@example.com")
    assert result == str(dest)


def test_fetch_genome_uses_bundled_test_data(tmp_path):
    result = fetch_genome("562", str(tmp_path), "test@example.com")
    assert result is not None
    assert Path(result).is_file()
    assert "562" in Path(result).name


def test_fetch_genome_invalid_taxid(test_output_dir):
    """Test fetching genome for an invalid taxid"""
    email = "test@example.com"
    try:
        result = fetch_genome("999999999", test_output_dir, email)
    except Exception as exc:
        pytest.skip(f"NCBI unavailable: {exc}")
    assert result is None


def test_fetch_genome_already_exists(test_output_dir):
    """Test that existing genome is not re-downloaded"""
    email = "test@example.com"
    result1 = fetch_genome("562", test_output_dir, email, reference_only=False)

    if result1 is None:
        pytest.skip("Could not obtain genome 562 from bundled test data or NCBI")

    mtime1 = os.path.getmtime(result1)
    result2 = fetch_genome("562", test_output_dir, email, reference_only=False)
    mtime2 = os.path.getmtime(result2)

    assert result1 == result2
    assert mtime1 == mtime2


def test_generate_random_taxids(test_output_dir):
    """Test generating random taxids"""
    from samovar.genome_fetcher import generate_random_taxids
    from Bio import Entrez

    Entrez.email = "test@example.com"
    try:
        taxids = generate_random_taxids()
    except Exception as exc:
        pytest.skip(f"NCBI unavailable: {exc}")
    if not taxids:
        pytest.skip("NCBI returned no taxids (API unavailable or empty search)")
    assert len(taxids) == 10
    assert all(isinstance(taxid, str) for taxid in taxids)
    assert len(set(taxids)) == len(taxids)

    try:
        taxids = generate_random_taxids(group="Archaea", N=5)
    except Exception as exc:
        pytest.skip(f"NCBI unavailable: {exc}")
    if not taxids:
        pytest.skip("NCBI returned no Archaea taxids (API unavailable or empty search)")
    assert len(taxids) == 5
    assert all(isinstance(taxid, str) for taxid in taxids)
    assert len(set(taxids)) == len(taxids)


def test_generate_random_taxids_no_email():
    """Test that function raises error when Entrez.email is not set"""
    from samovar.genome_fetcher import generate_random_taxids
    from Bio import Entrez

    if hasattr(Entrez, 'email'):
        delattr(Entrez, 'email')

    with pytest.raises(ValueError, match="Entrez.email must be set"):
        generate_random_taxids()


def test_generate_random_taxids_invalid_group(test_output_dir):
    """Test generating taxids for invalid group"""
    from samovar.genome_fetcher import generate_random_taxids
    from Bio import Entrez

    Entrez.email = "test@example.com"
    try:
        taxids = generate_random_taxids(group="InvalidGroup123")
    except Exception as exc:
        pytest.skip(f"NCBI unavailable: {exc}")
    assert len(taxids) == 0
