"""Accuracy purity, F1 purity, and annotator quality barplots."""

import pandas as pd
import pytest

from samovar.scores import (
    accuracy_purity,
    annotation_scores,
    f1_purity,
    majority_vote,
    purity_by_taxon,
    score_annotators,
)
from samovar.viz_annotation import viz_annotation


def test_accuracy_purity_majority_within_true_taxid():
    true = ["A", "A", "A", "B", "B", "B"]
    pred = ["A", "A", "A", "A", "B", "B"]
    assert accuracy_purity(true, pred) == pytest.approx(5 / 6)
    table = purity_by_taxon(true, pred)
    by = table.set_index("true")
    assert by.loc["A", "majority_pred"] == "A"
    assert by.loc["A", "recall"] == pytest.approx(1.0)
    assert by.loc["A", "precision"] == pytest.approx(0.75)
    assert by.loc["B", "recall"] == pytest.approx(2 / 3)
    assert by.loc["B", "precision"] == pytest.approx(1.0)


def test_current_f1_equals_accuracy_and_macro_is_standard():
    from samovar.scores import format_f1_caption, format_r2_caption, standard_classification_metrics

    true = ["A", "A", "A", "B", "B", "B"]
    pred = ["A", "A", "A", "A", "B", "B"]
    stats = annotation_scores(true, pred)
    assert stats["f1"] == pytest.approx(5 / 6)
    assert stats["accuracy"] == pytest.approx(5 / 6)
    assert stats["f1_micro"] == pytest.approx(5 / 6)
    f1_a = 2 * 0.75 * 1.0 / 1.75
    f1_b = 2 * 1.0 * (2 / 3) / (1.0 + 2 / 3)
    assert stats["f1_macro"] == pytest.approx(0.5 * f1_a + 0.5 * f1_b)
    std = standard_classification_metrics(true, pred)
    caption = format_f1_caption(stats["f1"], std)
    assert "TP/N" in caption
    assert "F1-macro" in caption
    assert "1 −" in format_r2_caption(0.5) or "SS" in format_r2_caption(0.5)


def test_f1_purity_penalises_external_calls_into_majority():
    true = ["A", "A", "A", "B", "B", "B"]
    pred = ["A", "A", "A", "A", "B", "B"]
    f1_a = 2 * 1.0 * 0.75 / (1.0 + 0.75)
    f1_b = 2 * (2 / 3) * 1.0 / ((2 / 3) + 1.0)
    assert f1_purity(true, pred) == pytest.approx(0.5 * f1_a + 0.5 * f1_b)
    assert f1_purity(true, pred) <= accuracy_purity(true, pred) + 1e-12


def test_consistent_wrong_label_is_pure_but_not_accurate():
    true = ["A", "A", "B", "B"]
    pred = ["X", "X", "Y", "Y"]
    assert accuracy_purity(true, pred) == pytest.approx(1.0)
    assert f1_purity(true, pred) == pytest.approx(1.0)
    assert annotation_scores(true, pred)["f1"] == pytest.approx(0.0)


def test_shared_majority_label_lowers_f1_purity():
    true = ["A", "A", "B", "B"]
    pred = ["X", "X", "X", "X"]
    assert accuracy_purity(true, pred) == pytest.approx(1.0)
    assert f1_purity(true, pred) == pytest.approx(2 / 3)


def test_unclassified_majority_is_zero_purity():
    true = ["A", "A", "A"]
    pred = ["0", "0", "0"]
    assert accuracy_purity(true, pred) == pytest.approx(0.0)
    assert f1_purity(true, pred) == pytest.approx(0.0)


def test_majority_vote_prefers_classified():
    df = pd.DataFrame(
        {
            "kaiju": ["562", "562", "0"],
            "kraken2": ["562", "9606", "9606"],
        }
    )
    votes = majority_vote(df, ["kaiju", "kraken2"])
    assert list(votes) == ["562", "562", "9606"]


def test_score_annotators_adds_consensus_and_keeps_samovar():
    df = pd.DataFrame(
        {
            "kaiju": ["562", "562", "9606", "9606"],
            "kraken2": ["562", "9606", "9606", "9606"],
            "SAMOVAR": ["562", "562", "9606", "9606"],
            "true": ["562", "562", "9606", "9606"],
        }
    )
    table = score_annotators(df, ["kaiju", "kraken2", "SAMOVAR"])
    names = list(table["annotator"].astype(str))
    assert names[-1] == "SAMOVAR"
    assert names.count("SAMOVAR") == 1
    assert "consensus" not in names
    assert "kaiju" in names and "kraken2" in names
    assert (table["accuracy_purity"] >= 0).all()
    assert (table["f1_purity"] >= 0).all()


def test_score_annotators_names_majority_vote_samovar():
    df = pd.DataFrame(
        {
            "kaiju": ["562", "562", "9606", "9606"],
            "kraken2": ["562", "9606", "9606", "9606"],
            "true": ["562", "562", "9606", "9606"],
        }
    )
    table = score_annotators(df, ["kaiju", "kraken2"])
    names = list(table["annotator"].astype(str))
    assert names[-1] == "SAMOVAR"
    assert "consensus" not in names


def test_score_annotators_without_ground_truth():
    df = pd.DataFrame(
        {
            "kaiju": ["562", "562", "9606", "0"],
            "kraken2": ["562", "9606", "9606", "9606"],
        }
    )
    table = score_annotators(df, ["kaiju", "kraken2"])
    names = list(table["annotator"].astype(str))
    assert "kaiju" in names and "kraken2" in names
    assert names[-1] == "SAMOVAR"
    assert table["f1"].isna().all()
    assert table["r2"].isna().all()
    assert (table["n_reads"] == 4).all()


def test_viz_annotation_writes_scores_barplot(tmp_path):
    df = pd.DataFrame(
        {
            "seq": [f"r{i}" for i in range(20)],
            "taxID_kaiju_0": [562] * 10 + [9606] * 8 + [0, 0],
            "taxID_kraken2_1": [562] * 8 + [9606] * 10 + [0, 0],
            "taxid_SAMOVAR": [562] * 10 + [9606] * 10,
            "true": [562] * 10 + [9606] * 10,
        }
    )
    out = tmp_path / "plots"
    results = viz_annotation(
        df,
        type=("scores",),
        show_top=0,
        output_dir=str(out),
        plot=False,
        rank="none",
        use_names=False,
    )
    assert "scores" in results
    assert not results["scores"].empty
    assert (out / "scores.png").is_file()
    assert (out / "quality_scores.csv").is_file()
    assert (out / "purity_by_taxon.csv").is_file()
    csv = pd.read_csv(out / "quality_scores.csv")
    assert "accuracy_purity" in csv.columns
    assert "f1_purity" in csv.columns
    assert "accuracy" in csv.columns
    assert "f1_macro" in csv.columns
    assert "tnr" in csv.columns
    assert "completeness" in csv.columns
    assert (out / "opal_scores.png").is_file()
    import matplotlib.image as mpimg

    for name in ("scores.png", "opal_scores.png"):
        arr = mpimg.imread(out / name)
        assert arr.shape[1] >= 1400, f"{name} width {arr.shape[1]} is too low-res"
        assert arr.shape[0] >= 700, f"{name} height {arr.shape[0]} is too low-res"
    assert "f1" in csv.columns
    assert "r2" in csv.columns
    assert "SAMOVAR" in set(csv["annotator"].astype(str))
    assert "consensus" not in set(csv["annotator"].astype(str))
