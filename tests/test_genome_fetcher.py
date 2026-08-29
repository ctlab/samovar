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


def test_parse_max_genome_mb_and_skip_list():
    from samovar.genome_fetcher import (
        UNLIMITED_GENOME_MB,
        genome_id_in_skip_list,
        parse_genome_skip_list,
        parse_max_genome_mb,
    )

    assert parse_max_genome_mb(None, default_from_env=False) == UNLIMITED_GENOME_MB
    assert parse_max_genome_mb("inf", default_from_env=False) == UNLIMITED_GENOME_MB
    assert parse_max_genome_mb(0, default_from_env=False) == UNLIMITED_GENOME_MB
    assert parse_max_genome_mb(100, default_from_env=False) == 100.0
    skip = parse_genome_skip_list("9606,GCF_009914755.1", default_from_env=False)
    assert genome_id_in_skip_list("9606", skip)
    assert genome_id_in_skip_list("GCF_009914755.1", skip)
    assert not genome_id_in_skip_list("562", skip)


def test_fetch_genome_skip_list_blocks_download(tmp_path, monkeypatch):
    from unittest.mock import patch

    monkeypatch.delenv("SAMOVAR_GENOME_SKIP_LIST", raising=False)
    monkeypatch.delenv("SAMOVAR_MAX_GENOME_MB", raising=False)
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    (tmp_path / "config.json").write_text("{}\n")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "0")
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)

    with patch("samovar.genome_fetcher._download_url") as download:
        with patch("samovar.genome_fetcher._assembly_ftp_path") as ftp:
            result = fetch_genome(
                "9606",
                str(tmp_path / "run"),
                "test@example.com",
                genome_skip_list="9606",
            )
    assert result is None
    download.assert_not_called()
    ftp.assert_not_called()


def test_fetch_genome_skip_list_allows_existing(tmp_path, monkeypatch):
    monkeypatch.delenv("SAMOVAR_GENOME_SKIP_LIST", raising=False)
    dest = tmp_path / "9606-processed.fasta.gz"
    with gzip.open(dest, "wt") as handle:
        handle.write(">x\nACGT\n")
    result = fetch_genome(
        "9606",
        str(tmp_path),
        "test@example.com",
        genome_skip_list="9606",
        max_genome_mb=1.0,
    )
    assert result == str(dest)


def test_fetch_genome_max_mb_blocks_new_download(tmp_path, monkeypatch):
    from unittest.mock import patch

    monkeypatch.delenv("SAMOVAR_MAX_GENOME_MB", raising=False)
    monkeypatch.delenv("SAMOVAR_GENOME_SKIP_LIST", raising=False)
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    (tmp_path / "config.json").write_text("{}\n")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "0")
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)

    with patch(
        "samovar.genome_fetcher._assembly_ftp_path",
        return_value="https://example.invalid/GCF_fake",
    ):
        with patch("samovar.genome_fetcher._remote_fasta_size_mb", return_value=900.0):
            with patch("samovar.genome_fetcher._download_url") as download:
                result = fetch_genome(
                    "562",
                    str(tmp_path / "run"),
                    "test@example.com",
                    max_genome_mb=100,
                )
    assert result is None
    download.assert_not_called()


def test_fetch_genome_default_max_mb_is_unlimited(tmp_path, monkeypatch):
    from unittest.mock import patch

    monkeypatch.delenv("SAMOVAR_MAX_GENOME_MB", raising=False)
    monkeypatch.delenv("SAMOVAR_GENOME_SKIP_LIST", raising=False)
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    (tmp_path / "config.json").write_text("{}\n")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "0")
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)

    def fake_download(url, dest, timeout=45):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"\x1f\x8b")

    with patch(
        "samovar.genome_fetcher._assembly_ftp_path",
        return_value="https://example.invalid/GCF_fake",
    ):
        with patch("samovar.genome_fetcher._remote_fasta_size_mb") as head:
            with patch("samovar.genome_fetcher._download_url", side_effect=fake_download) as download:
                with patch(
                    "samovar.genome_fetcher._assembly_doc_for_taxid",
                    return_value=None,
                ):
                    with patch(
                        "samovar.genome_fetcher.preprocess_fasta",
                        side_effect=lambda **kw: Path(kw["output_file"]).write_text(">x\nACGT\n"),
                    ):
                        fetch_genome("1280", str(tmp_path / "run"), "test@example.com")
    download.assert_called()
    head.assert_not_called()


def test_fetch_genome_downloads_into_raw(tmp_path, monkeypatch):
    from unittest.mock import patch

    from samovar.genome_fetcher import fetch_genome

    monkeypatch.delenv("SAMOVAR_MAX_GENOME_MB", raising=False)
    monkeypatch.delenv("SAMOVAR_GENOME_SKIP_LIST", raising=False)
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    (tmp_path / "config.json").write_text("{}\n")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "0")
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)

    dests = []

    def fake_download(url, dest, timeout=45):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"\x1f\x8b")
        dests.append(Path(dest))

    with patch(
        "samovar.genome_fetcher._assembly_ftp_path",
        return_value="https://example.invalid/GCF_fake",
    ):
        with patch("samovar.genome_fetcher._download_url", side_effect=fake_download):
            with patch(
                "samovar.genome_fetcher._assembly_doc_for_taxid",
                return_value=None,
            ):
                with patch(
                    "samovar.genome_fetcher.preprocess_fasta",
                    side_effect=lambda **kw: Path(kw["output_file"]).write_text(">x\nACGT\n"),
                ):
                    fetch_genome("1280", str(tmp_path / "run"), "test@example.com")
    assert dests
    assert dests[0].parent.name == "raw"
    assert dests[0].name == "1280.fna.gz"


def test_cached_ncbi_gzip_size_mismatch_warns(tmp_path, monkeypatch, caplog):
    import logging
    from unittest.mock import patch

    from samovar.genome_cache import genome_download_dir
    from samovar.genome_fetcher import fetch_genome

    monkeypatch.delenv("SAMOVAR_MAX_GENOME_MB", raising=False)
    monkeypatch.delenv("SAMOVAR_GENOME_SKIP_LIST", raising=False)
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    (tmp_path / "config.json").write_text("{}\n")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "0")
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)

    raw = genome_download_dir() / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    cached = raw / "1280.fna.gz"
    cached.write_bytes(b"\x1f\x8b")

    caplog.set_level(logging.WARNING)
    with patch(
        "samovar.genome_fetcher._assembly_ftp_path",
        return_value="https://example.invalid/GCF_fake",
    ):
        with patch(
            "samovar.genome_fetcher._remote_fasta_size_bytes",
            return_value=5_000_000,
        ):
            with patch("samovar.genome_fetcher._download_url") as download:
                with patch(
                    "samovar.genome_fetcher.preprocess_fasta",
                    side_effect=lambda **kw: Path(kw["output_file"]).write_text(">x\nACGT\n"),
                ):
                    fetch_genome("1280", str(tmp_path / "run"), "test@example.com")
    download.assert_not_called()
    assert "NCBI reports 5000000 bytes" in caplog.text


def test_gzip_sizes_disagree_tolerance():
    from samovar.genome_fetcher import gzip_sizes_disagree

    assert not gzip_sizes_disagree(1_000_000, 1_000_000)
    assert not gzip_sizes_disagree(1_000, 1_050)
    assert gzip_sizes_disagree(2, 5_000_000)


def test_raw_parent_is_sibling_of_processed(tmp_path):
    from samovar.genome_fetcher import _raw_parent_for_dest

    processed = tmp_path / ".genomes" / "processed"
    processed.mkdir(parents=True)
    raw = _raw_parent_for_dest(processed)
    assert raw == tmp_path / ".genomes" / "raw"
    assert raw.is_dir()


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
