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


def test_fetch_genome_reuses_library_not_bundled(tmp_path, monkeypatch):
    """Truncated test_genomes must not satisfy fetch_genome for a missing taxid."""
    from unittest.mock import patch

    cfg = tmp_path / "config.json"
    cfg.write_text("{}\n")
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "1")
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)

    lib = tmp_path / "lib"
    lib.mkdir()
    processed = lib / "562-processed.fasta.gz"
    with gzip.open(processed, "wt") as handle:
        handle.write(">x\nACGT\n")
    monkeypatch.setenv("SAMOVAR_GENOME_DIRS", str(lib))

    result = fetch_genome("562", str(tmp_path / "run"), "test@example.com")
    assert result is not None
    assert Path(result).is_file()
    assert "562" in Path(result).name

    empty = tmp_path / "empty_run"
    empty.mkdir()
    with patch("samovar.genome_fetcher._assembly_ftp_path", return_value=None):
        missing = fetch_genome("999999001", str(empty), "test@example.com")
    assert missing is None


def test_fetch_genome_does_not_use_bundled_test_data(tmp_path, monkeypatch):
    from unittest.mock import patch

    cfg = tmp_path / "config.json"
    cfg.write_text("{}\n")
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "0")
    with patch("samovar.genome_fetcher._assembly_ftp_path", return_value=None):
        result = fetch_genome("562", str(tmp_path / "run"), "test@example.com")
    assert result is None


def test_fetch_genome_invalid_taxid(test_output_dir):
    """Test fetching genome for an invalid taxid"""
    email = "test@example.com"
    try:
        result = fetch_genome("999999999", test_output_dir, email)
    except Exception as exc:
        pytest.skip(f"NCBI unavailable: {exc}")
    assert result is None


def test_fetch_genome_already_exists(tmp_path):
    """Test that existing processed genome is not re-downloaded"""
    dest = tmp_path / "562-processed.fasta.gz"
    with gzip.open(dest, "wt") as handle:
        handle.write(">x\nACGT\n")
    result1 = fetch_genome("562", str(tmp_path), "test@example.com")
    mtime1 = os.path.getmtime(result1)
    result2 = fetch_genome("562", str(tmp_path), "test@example.com")
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
