"""Shannon KS scoring for table_reads_generator selection."""

import numpy as np
import pandas as pd
import pytest

from samovar.annotation_io import read_annotation_dir
from samovar.table_scorers import (
    canonicalize_table_scorer,
    ks_shannon,
    pick_best_table_method,
    score_generated_tables,
    shannon_entropy,
    shannon_vector_from_annotation,
)


@pytest.fixture
def toy_annotation_dir(tmp_path):
    ann = tmp_path / "annotations"
    ann.mkdir()
    pd.DataFrame(
        {
            "seq": ["a", "b", "c"],
            "taxID_kaiju_0": [562, 562, 9606],
            "taxID_kraken2_0": [562, 9606, 9606],
            "true": [562, 562, 9606],
        }
    ).to_csv(ann / "1.annotation.csv", index=False)
    pd.DataFrame(
        {
            "seq": ["a", "b"],
            "taxID_kaiju_0": [562, 9606],
            "taxID_kraken2_0": [9606, 9606],
            "true": [562, 9606],
        }
    ).to_csv(ann / "2.annotation.csv", index=False)
    return ann
from samovar.table_scorers import (
    canonicalize_table_scorer,
    ks_shannon,
    pick_best_table_method,
    score_generated_tables,
    shannon_entropy,
    shannon_vector_from_annotation,
)


def test_shannon_entropy_uniform_two_taxa():
    value = shannon_entropy([1.0, 1.0])
    assert abs(value - np.log(2)) < 1e-9


def test_ks_identical_vectors_is_zero():
    vec = np.array([0.1, 0.4, 0.7])
    row = ks_shannon(vec, vec)
    assert row["ks_statistic"] == 0.0
    assert row["pvalue"] == 1.0
    assert row["rank_value"] == 0.0


def test_pick_best_prefers_lower_ks():
    winner = pick_best_table_method(
        [
            {"mode": "bootstrap", "rank_value": 0.4, "pvalue": 0.9, "ok": True},
            {"mode": "direct", "rank_value": 0.0, "pvalue": 1.0, "ok": True},
            {"mode": "glm", "rank_value": 1.0, "pvalue": 0.0, "ok": False},
        ]
    )
    assert winner == "direct"


def test_score_generated_tables_on_toy(toy_annotation_dir):
    from samovar.regenerate import regenerate_annotation_tables

    data = read_annotation_dir(toy_annotation_dir)
    obs = shannon_vector_from_annotation(data)
    assert obs.size >= 2
    tables = regenerate_annotation_tables(
        toy_annotation_dir, toy_annotation_dir.parent / "direct", {"regeneration_mode": "direct"}
    )
    score = score_generated_tables(data, tables, "shannon_ks")
    assert score["ks_statistic"] == 0.0
    assert canonicalize_table_scorer("alpha-ks") == "shannon_ks"
