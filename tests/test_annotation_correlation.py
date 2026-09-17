"""Spearman tax / feat / true correlation helper."""

import json

import pandas as pd

from samovar.annotation_correlation import (
    annotation_stage_allows_correlation,
    plot_spearman_annotation_correlation,
    spearman_annotation_correlation,
)


def test_spearman_includes_tax_feat_true():
    df = pd.DataFrame(
        {
            "seq": [f"r{i}" for i in range(12)],
            "taxID_kaiju_0": [562, 562, 9606, 9606] * 3,
            "taxID_kraken2_1": [562, 9606, 9606, 562] * 3,
            "feat_len_0": [10, 20, 30, 40] * 3,
            "true": [562, 562, 9606, 9606] * 3,
        }
    )
    corr = spearman_annotation_correlation(df)
    names = set(corr.index.astype(str))
    assert "taxID_kaiju_0" in names
    assert "taxID_kraken2_1" in names
    assert "feat_len_0" in names
    assert "true" in names
    assert corr.loc["true", "true"] == 1.0
    assert corr.loc["taxID_kaiju_0", "true"] > 0.9
    assert pd.notna(corr.loc["feat_len_0", "true"])
    assert corr.loc["feat_len_0", "true"] != 0.0


def test_constant_feat_vs_true_is_zero_not_nan():
    df = pd.DataFrame(
        {
            "taxID_kaiju_0": [562, 562, 9606, 9606] * 3,
            "feat_kraken2_1_length": [126] * 12,
            "true": [562, 562, 9606, 9606] * 3,
        }
    )
    corr = spearman_annotation_correlation(df)
    assert "feat_kraken2_1_length" in set(corr.index.astype(str))
    rho = float(corr.loc["feat_kraken2_1_length", "true"])
    assert rho == 0.0
    assert pd.notna(rho)


def test_top_20_features_sorted_by_importance():
    n = 25
    data = {"taxID_kaiju_0": [562, 9606] * 6, "true": [562, 9606] * 6}
    for i in range(n):
        data[f"feat_kmer2_0_{i:02d}"] = list(range(12))
    df = pd.DataFrame(data)
    importance = pd.Series(
        {f"feat_kmer2_0_{i:02d}": float(i) for i in range(n)}
        | {"taxid_kaiju": 100.0}
    )
    corr = spearman_annotation_correlation(df, importance=importance)
    names = [str(c) for c in corr.columns]
    feats = [c for c in names if c.startswith("feat_")]
    assert "taxID_kaiju_0" in names
    assert "true" in names
    assert names[-1] == "true"
    assert len(feats) == 20
    assert "feat_kmer2_0_00" not in feats
    assert "feat_kmer2_0_24" in feats
    assert names[0] == "taxID_kaiju_0"
    feat_order = [int(c.rsplit("_", 1)[-1]) for c in feats]
    assert feat_order == sorted(feat_order, reverse=True)


def test_writes_multiqc_and_optional_plots(tmp_path):
    df = pd.DataFrame(
        {
            "taxID_kaiju_0": [562] * 6 + [9606] * 6,
            "taxID_kraken2_1": [562] * 5 + [9606] * 7,
            "true": [562] * 6 + [9606] * 6,
        }
    )
    dest = tmp_path / "regenerated_annotations_plots"
    matrix = plot_spearman_annotation_correlation(df, dest)
    assert matrix is not None
    payload = json.loads((dest / "spearman_correlation_mqc.json").read_text())
    assert payload["plot_type"] == "heatmap"
    assert payload["pconfig"]["min"] == -1
    assert payload["pconfig"]["max"] == 1
    assert (dest / "spearman_correlation.tsv").is_file()


def test_stage_gate():
    assert annotation_stage_allows_correlation(
        "run/regenerated_annotations", "run/regenerated_annotations_plots"
    )
    assert annotation_stage_allows_correlation(
        None, "run/reprofiled_annotations_plots"
    )
    assert not annotation_stage_allows_correlation(
        "run/initial_annotations", "run/initial_annotations_plots"
    )
