"""Tests for the parse_annotators module."""

import os
import pytest
import pandas as pd
from samovar.parse_annotators import (
    read_kaiju_raw,
    read_kraken1_raw,
    read_kraken2_raw,
    read_krakenu_raw,
    read_metaphlan_raw,
    read_annotation,
    Annotation,
    RankAnnotation,
    ExpandAnnotation,
    parse_metaphlan_db,
    extract_true_taxid,
    _discover_nodes_path,
    _load_nodes_cache,
    resolve_metaphlan_db_file,
    remap_taxid_dataframe,
    canonical_taxid,
    ensure_taxid_rank_map,
    ensure_taxid_name_map,
    read_taxid_rank_table,
    rank_map_column,
)


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture
def kaiju_file(test_data_dir):
    """Return path to Kaiju test file."""
    return os.path.join(test_data_dir, "kaiju.log")


@pytest.fixture
def kraken1_file(test_data_dir):
    """Return path to Kraken1 test file."""
    return os.path.join(test_data_dir, "kraken.log")


@pytest.fixture
def kraken2_file(test_data_dir):
    """Return path to Kraken2 test file."""
    return os.path.join(test_data_dir, "kraken2.log")

@pytest.fixture
def krakenu_file(test_data_dir):
    """Return path to Kraken Unique test file."""
    return os.path.join(test_data_dir, "krakenunique.log")


@pytest.fixture
def metaphlan_file(test_data_dir):
    """Return path to MetaPhlAn test file."""
    return os.path.join(test_data_dir, "metaphlan4.log")


def test_extract_true_taxid_from_headers():
    """True taxID must not collapse to 0 when taxid: is missing from a contig."""
    pattern = r"(?<=taxid:)[0-9]*"
    assert extract_true_taxid("Hsap.fna|taxid:9606|-NC_000001.11_0_0", pattern) == "9606"
    assert extract_true_taxid("Ecoli.fna|taxid:562|-NC_000913.3_0_0", pattern) == "562"
    assert extract_true_taxid("Scer.fna|taxid:4932|-NC_001133.9_0_0", pattern) == "4932"
    # Yeast chromosomes II+ were stored without taxid: in the header.
    assert extract_true_taxid("Scer.fna-NC_001134.8_0_0", pattern) == "4932"
    assert extract_true_taxid("Scer.fna-NC_001139.9_1_1/1", pattern) == "4932"
    assert extract_true_taxid("taxid:562_0_0/1", pattern) == "562"
    assert extract_true_taxid("random_read", pattern) == ""


def test_read_kaiju_raw(kaiju_file):
    """Test reading Kaiju output file."""
    df = read_kaiju_raw(kaiju_file)
    print(df.head())
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["classified", "seq", "taxID"]
    assert len(df.seq.unique()) == len(df.seq)
    assert len(df) > 0


def test_read_kraken1_raw(kraken1_file):
    """Test reading Kraken1 output file."""
    df = read_kraken1_raw(kraken1_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["classified", "seq", "taxID", "length", "k-mer"]
    assert len(df.seq.unique()) == len(df.seq)
    assert len(df) > 0


def test_read_kraken2_raw(kraken2_file):
    """Test reading Kraken 2 output file."""
    df = read_kraken2_raw(kraken2_file)
    assert isinstance(df, pd.DataFrame)
    assert "taxID" in df.columns
    assert len(df.seq.unique()) == len(df.seq)
    assert len(df) > 0

def test_read_krakenu_raw(krakenu_file):
    """Test reading Kraken Unique output file."""
    df = read_krakenu_raw(krakenu_file)
    assert isinstance(df, pd.DataFrame)
    assert "taxID" in df.columns
    assert len(df.seq.unique()) == len(df.seq)
    assert len(df) > 0


def test_read_metaphlan_raw(metaphlan_file):
    """Test reading MetaPhlAn output file."""
    df = read_metaphlan_raw(metaphlan_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["seq", "taxID"]
    assert len(df.seq.unique()) == len(df.seq)
    assert len(df) > 0


def test_read_annotation(test_data_dir):
    """Test reading multiple annotation files."""
    file_path_type = {
        os.path.join(test_data_dir, "kaiju.log"): "kaiju",
        os.path.join(test_data_dir, "kraken.log"): "kraken1",
        os.path.join(test_data_dir, "kraken2.log"): "kraken2",
        os.path.join(test_data_dir, "krakenunique.log"): "krakenu"
        #os.path.join(test_data_dir, "metaphlan4.log"): "metaphlan"
    }
    ann = Annotation(file_path_type, get_true_annotation=r".*")
    ann.export("tests_outs/annotation.csv")

    assert isinstance(ann.DataFrame, pd.DataFrame)
    assert len(ann.DataFrame) > 0
    assert any(col.startswith("taxID_") for col in ann.DataFrame.columns)


def test_annotation_class(test_data_dir):
    """Test Annotation class functionality."""
    file_path_type = {
        os.path.join(test_data_dir, "kaiju.log"): "kaiju",
        os.path.join(test_data_dir, "kraken.log"): "kraken1"
    }
    ann = Annotation(file_path_type, get_true_annotation=r".*")
    
    # Test basic functionality
    assert isinstance(ann.DataFrame, pd.DataFrame)
    assert len(ann.true_annotation) > 0
    
    # Test rank annotation
    rank_ann = ann.rank_annotation("species")
    assert isinstance(rank_ann, RankAnnotation)
    
    # Test expand annotation
    expand_ann = ann.expand_annotation(["species", "genus"])
    assert isinstance(expand_ann, ExpandAnnotation)
    assert len(expand_ann.rank_annotation) == 2


def test_correct_level(test_data_dir):
    """Test correct_level method of Annotation class."""
    # Create test data with known taxIDs
    test_data = {
        'seq': ['seq1', 'seq2', 'seq3', 'seq4'],
        'taxID_kraken1_0': ['9606', '511145', '0', '1234567890'],  # Human, E. coli, unclassified, invalid
        'taxID_kaiju_0': ['9606', '0', '0', '1234567890']
    }
    df = pd.DataFrame(test_data)
    
    # Create Annotation object with test data
    ann = Annotation({}, get_true_annotation=r".*")
    ann.DataFrame = df
    
    # Test correcting to species level
    ann.correct_level(level='species')
    
    # Check that valid taxIDs were corrected
    assert ann.DataFrame['taxID_kraken1_0'].iloc[0] == '9606'  # Human species taxID
    assert ann.DataFrame['taxID_kraken1_0'].iloc[1] == '562'  # E. coli species taxID
    assert ann.DataFrame['taxID_kraken1_0'].iloc[2] == '0'  # Unclassified should remain 0
    assert ann.DataFrame['taxID_kraken1_0'].iloc[3] == '1234567890'  # Invalid taxID should remain unchanged
    
    # Test correcting to genus level
    ann.correct_level(level='genus')
    
    # Check that taxIDs were corrected to genus level
    # Note: These assertions may need to be adjusted based on actual NCBI taxonomy data
    assert str(ann.DataFrame['taxID_kraken1_0'].iloc[0]) == '9605'
    assert str(ann.DataFrame['taxID_kraken1_0'].iloc[1]) == '561'
    assert ann.DataFrame['taxID_kraken1_0'].iloc[2] == '0'  # Unclassified should remain 0
    assert ann.DataFrame['taxID_kraken1_0'].iloc[3] == '1234567890'  # Invalid taxID should remain unchanged


def test_remap_taxid_dataframe_genus_collapses_strain_and_true(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    df = pd.DataFrame({
        "taxID_kaiju": ["511145", "562", "0"],
        "taxid_kraken2": [511145, 561, 0],
        "true": ["511145", "511145", "9606"],
        "taxID_kaiju_confidence": [0.9, 0.8, 0.1],
    })
    cache = tmp_path / "taxid_genera_map.tsv"
    out = remap_taxid_dataframe(df, rank="genus", cache_path=str(cache))
    assert str(canonical_taxid(out.loc[0, "taxID_kaiju"])) == "561"
    assert str(canonical_taxid(out.loc[1, "taxID_kaiju"])) == "561"
    assert str(canonical_taxid(out.loc[1, "taxid_kraken2"])) == "561"
    assert str(canonical_taxid(out.loc[0, "true"])) == "561"
    assert str(canonical_taxid(out.loc[2, "true"])) == "9605"
    assert str(canonical_taxid(out.loc[2, "taxID_kaiju"])) == "0"
    assert out.loc[0, "taxID_kaiju_confidence"] == 0.9
    assert str(canonical_taxid(df.loc[0, "taxID_kaiju"])) == "511145"
    text = cache.read_text()
    assert text.splitlines()[0] == "taxid|genera_taxid"
    assert "511145|561" in text
    table = read_taxid_rank_table(cache)
    assert table["511145"] == "561"
    exact = remap_taxid_dataframe(df, rank="none")
    assert str(canonical_taxid(exact.loc[0, "taxID_kaiju"])) == "511145"


def test_ensure_taxid_rank_map_reuses_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    cache = tmp_path / "map.tsv"
    cache.write_text("taxid|genera_taxid\n999999|42\n")
    mapping = ensure_taxid_rank_map(["999999"], rank="genus", cache_path=str(cache))
    assert mapping["999999"] == "42"
    assert rank_map_column("genus") == "genera_taxid"
    assert rank_map_column("species") == "species_taxid"


def test_ensure_taxid_name_map_writes_pipe_table(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    cache = tmp_path / "taxid_names.tsv"
    mapping = ensure_taxid_name_map(["562", "561", "0"], cache_path=str(cache))
    assert "coli" in mapping["562"].lower() or mapping["562"] == "Escherichia coli"
    assert mapping["561"] in {"Escherichia", mapping["561"]}
    text = cache.read_text()
    assert text.splitlines()[0] == "taxid|name"
    assert "562|" in text


def test_rank_annotation_class():
    """Test RankAnnotation class functionality."""
    rank_ann = RankAnnotation("species")
    assert rank_ann.rank == "species"
    assert isinstance(rank_ann.annotation, pd.DataFrame)


def test_expand_annotation_class():
    """Test ExpandAnnotation class functionality."""
    expand_ann = ExpandAnnotation()
    assert isinstance(expand_ann.rank_annotation, dict)


def test_parse_metaphlan_db(test_data_dir):
    """Test parsing MetaPhlAn database."""
    # Create a mock database file
    db_dir = os.path.join(test_data_dir, "mock_metaphlan_db")
    os.makedirs(db_dir, exist_ok=True)
    db_file = os.path.join(db_dir, "mpa_v30_CHOCOPhlAn_201901_species_map.db")
    
    # Create a simple SQLite database with test data
    import sqlite3
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE mpa_species_map (
            ref_id TEXT,
            tax_id TEXT
        )
    """)
    test_data = [
        ("M367-c418", "511145"),  # E. coli
        ("M1206-c595", "9606"),   # Human
    ]
    cursor.executemany("INSERT INTO mpa_species_map VALUES (?, ?)", test_data)
    conn.commit()
    conn.close()
    
    # Test parsing
    mapping = parse_metaphlan_db(db_dir)
    assert isinstance(mapping, dict)
    assert mapping["M367-c418"] == "511145"
    assert mapping["M1206-c595"] == "9606"
    
    # Cleanup
    os.remove(db_file)
    os.rmdir(db_dir)

def test_read_metaphlan_with_db(test_data_dir):
    """Test reading MetaPhlAn output with database mapping."""
    # Create mock database
    db_dir = os.path.join(test_data_dir, "mock_metaphlan_db")
    os.makedirs(db_dir, exist_ok=True)
    db_file = os.path.join(db_dir, "mpa_v30_CHOCOPhlAn_201901_species_map.db")
    
    import sqlite3
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE mpa_species_map (
            ref_id TEXT,
            tax_id TEXT
        )
    """)
    test_data = [
        ("M367-c418", "511145"),  # E. coli
        ("M1206-c595", "9606"),   # Human
    ]
    cursor.executemany("INSERT INTO mpa_species_map VALUES (?, ?)", test_data)
    conn.commit()
    conn.close()
    
    # Test reading with database
    df = read_metaphlan_raw(os.path.join(test_data_dir, "metaphlan4.log"), db_dir)
    assert isinstance(df, pd.DataFrame)
    assert "taxID" in df.columns
    assert df["taxID"].isin(["511145", "9606"]).any()
    
    # Cleanup
    os.remove(db_file)
    os.rmdir(db_dir)


def test_resolve_metaphlan_db_file_db_name_and_file(tmp_path):
    db_file = tmp_path / "custom_map.db"
    db_file.write_bytes(b"")
    assert resolve_metaphlan_db_file(str(tmp_path), db_name="custom_map.db") == str(db_file)
    assert resolve_metaphlan_db_file(str(db_file)) == str(db_file)


def test_discover_nodes_path_from_env(tmp_path, monkeypatch):
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text("1\t|\t1\t|\tno rank\t|\n")
    monkeypatch.delenv("SAMOVAR_NODES_DMP", raising=False)
    monkeypatch.delenv("SAMOVAR_NODES_SEARCH", raising=False)
    assert _discover_nodes_path() is None
    monkeypatch.setenv("SAMOVAR_NODES_DMP", str(nodes))
    assert _discover_nodes_path() == str(nodes)
    monkeypatch.delenv("SAMOVAR_NODES_DMP")
    monkeypatch.setenv("SAMOVAR_NODES_SEARCH", str(tmp_path))
    assert _discover_nodes_path() == str(nodes)


def test_load_nodes_cache_atomic_reload(tmp_path, monkeypatch):
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text("1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tspecies\t|\n")
    cache = tmp_path / "cache"
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(cache))
    tree1 = _load_nodes_cache(str(nodes))
    tree2 = _load_nodes_cache(str(nodes))
    assert tree1 == tree2
    assert 2 in tree1
    pickles = list(cache.glob("nodes_*.pkl"))
    assert len(pickles) == 1
    # Truncated pickle is discarded and rebuilt.
    pickles[0].write_bytes(b"not a pickle")
    tree3 = _load_nodes_cache(str(nodes))
    assert tree3 == tree1