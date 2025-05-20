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
    ExpandAnnotation
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
def metaphlan_file(test_data_dir):
    """Return path to MetaPhlAn test file."""
    return os.path.join(test_data_dir, "metaphlan4.log")


def test_read_kaiju_raw(kaiju_file):
    """Test reading Kaiju output file."""
    df = read_kaiju_raw(kaiju_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["classified", "seq", "score", "taxID", "N"]
    assert len(df) > 0


def test_read_kraken1_raw(kraken1_file):
    """Test reading Kraken1 output file."""
    df = read_kraken1_raw(kraken1_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["classified", "seq", "taxID", "length", "k-mer"]
    assert len(df) > 0


def test_read_kraken2_raw(kraken2_file):
    """Test reading Kraken2 output file."""
    df = read_kraken2_raw(kraken2_file)
    assert isinstance(df, pd.DataFrame)
    assert "taxID" in df.columns
    assert len(df) > 0


def test_read_metaphlan_raw(metaphlan_file):
    """Test reading MetaPhlAn output file."""
    df = read_metaphlan_raw(metaphlan_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["seq", "taxID"]
    assert len(df) > 0


def test_read_annotation(test_data_dir):
    """Test reading multiple annotation files."""
    file_path_type = {
        os.path.join(test_data_dir, "kaiju.log"): "kaiju",
        os.path.join(test_data_dir, "kraken.log"): "kraken1"
    }
    df = read_annotation(file_path_type)
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert any(col.startswith("taxID_") for col in df.columns)


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


def test_rank_annotation_class():
    """Test RankAnnotation class functionality."""
    rank_ann = RankAnnotation("species")
    assert rank_ann.rank == "species"
    assert isinstance(rank_ann.annotation, pd.DataFrame)


def test_expand_annotation_class():
    """Test ExpandAnnotation class functionality."""
    expand_ann = ExpandAnnotation()
    assert isinstance(expand_ann.rank_annotation, dict) 