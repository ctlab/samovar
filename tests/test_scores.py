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


def test_majority_vote_tie_keeps_the_first_column():
    df = pd.DataFrame({"a": ["562", "0"], "b": ["9606", "9606"], "c": ["562", "0"]})
    assert list(majority_vote(df, ["a", "b", "c"])) == ["562", "9606"]
    tied = pd.DataFrame({"a": ["111"], "b": ["222"]})
    assert list(majority_vote(tied, ["a", "b"])) == ["111"]
    assert list(majority_vote(tied, ["b", "a"])) == ["222"]


def test_default_scores_skip_opal_and_split_read_types_once():
    df = pd.DataFrame(
        {
            "kaiju": ["562", "562", "9606", "9606"],
            "kraken2": ["562", "9606", "9606", "9606"],
            "true": ["562", "562", "9606", "9606"],
            "read_type": ["illumina", "illumina", "ont", "ont"],
        }
    )
    full = score_annotators(df.drop(columns=["read_type"]), ["kaiju", "kraken2"], opal=True)
    plain = score_annotators(df, ["kaiju", "kraken2"])
    overall = plain[plain["read_type"] == "all"]
    assert overall["completeness"].isna().all()
    assert overall["accuracy"].tolist() == full["accuracy"].tolist()
    assert set(plain["read_type"]) == {"all", "illumina", "ont"}
    sampled = score_annotators(df.drop(columns=["read_type"]), ["kaiju"], metrics_sample=2, seed=1)
    assert int(sampled["n_sampled"].iloc[0]) == 2


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


@pytest.mark.parametrize("opal", [False, True])
def test_viz_annotation_writes_scores_barplot(tmp_path, monkeypatch, opal):
    calls = []
    monkeypatch.setattr("samovar.opal.maybe_run_opal", lambda *args, **kwargs: calls.append(kwargs))
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
        opal=opal,
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
    assert (out / "opal_scores.png").is_file() == opal
    assert (out / "opal_scores.html").is_file() == opal
    assert bool(calls) == opal
    assert csv["completeness"].notna().all() if opal else csv["completeness"].isna().all()
    import matplotlib.image as mpimg

    for name in (("scores.png", "opal_scores.png") if opal else ("scores.png",)):
        arr = mpimg.imread(out / name)
        assert arr.shape[1] >= 1400, f"{name} width {arr.shape[1]} is too low-res"
        assert arr.shape[0] >= 700, f"{name} height {arr.shape[0]} is too low-res"
    assert "f1" in csv.columns
    assert "r2" in csv.columns
    assert "SAMOVAR" in set(csv["annotator"].astype(str))
    assert "consensus" not in set(csv["annotator"].astype(str))


def test_classification_marginals_match_sklearn_and_opal_reference():
    import numpy as np
    from sklearn.metrics import precision_recall_fscore_support, multilabel_confusion_matrix
    from samovar.scores import standard_classification_metrics
    from samovar.opal import confusion_rates

    rng = np.random.default_rng(42)
    y = rng.choice(["A", "B", "C", "0"], 1003)
    p = rng.choice(["A", "B", "D", "0"], 1003)
    result = standard_classification_metrics(y, p)
    for average, keys in [("macro", ("precision_macro", "recall_macro", "f1_macro")),
                          ("micro", (None, None, "f1_micro")),
                          ("weighted", (None, None, "f1_weighted"))]:
        values = precision_recall_fscore_support(y, p, average=average, zero_division=0)
        for key, value in zip(keys, values):
            if key:
                assert result[key] == pytest.approx(value, abs=1e-12)
    assert result["accuracy"] == pytest.approx(float((y == p).mean()), abs=1e-12)
    categories = ["A", "B", "C", "D", "0", "unused"]
    categorical = standard_classification_metrics(
        pd.Series(pd.Categorical(y, categories=categories)),
        pd.Series(pd.Categorical(p, categories=categories)),
    )
    assert categorical == pytest.approx(result, abs=1e-12)
    matrices = multilabel_confusion_matrix(y, p, labels=["A", "B", "C"])
    tn, fp, fn, tp = (matrices[:, i, j] for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)])
    rates = confusion_rates(y, p)
    for key, a, b in [("tnr", tn, tn + fp), ("tpr", tp, tp + fn),
                       ("fpr", fp, fp + tn), ("fnr", fn, fn + tp)]:
        expected = np.divide(a, b, out=np.zeros_like(a, dtype=float), where=b > 0).mean()
        assert rates[key] == pytest.approx(expected, abs=1e-12)


def test_selected_metrics_skip_unused_computations(monkeypatch):
    import samovar.scores as scores
    def forbidden(*args, **kwargs):
        raise AssertionError("disabled metric was computed")
    monkeypatch.setattr(scores, "annotation_scores", forbidden)
    monkeypatch.setattr(scores, "_r2_abundance", forbidden)
    frame = pd.DataFrame({"kaiju": ["562", "0"], "true": ["562", "562"]})
    result = scores.score_annotators(frame, ["kaiju"], metrics=["accuracy"])
    assert result.loc[0, "accuracy"] == 0.5
    assert result["r2"].isna().all() and result["f1_purity"].isna().all()
    assert result.loc[0, "n_reads"] == 2
