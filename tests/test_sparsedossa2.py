"""SparseDOSSA2 table generator / CV scorer wiring (R is mocked)."""

from pathlib import Path

import pandas as pd
import pytest

from samovar.regenerate import resolve_regeneration_mode
from samovar.sparsedossa2 import (
    canonicalize_sparsedossa2_mode,
    remap_simulated_features,
    simulate_count_matrix,
    template_for_mode,
)
from samovar.table_regenerators import apply_extra_flags, get_table_regenerator
from samovar.table_scorers import canonicalize_table_scorer, score_generated_tables


def test_sparsedossa2_mode_aliases():
    assert resolve_regeneration_mode("sd2") == ("builtin", "sparsedossa2-fit")
    assert resolve_regeneration_mode("sparsedossa2-stool") == (
        "builtin",
        "sparsedossa2-stool",
    )
    assert canonicalize_sparsedossa2_mode("sparsedossa2_vaginal") == "sparsedossa2-vaginal"
    assert template_for_mode("sparsedossa2-ibd") == "IBD"
    assert template_for_mode("sparsedossa2-fit") is None
    regen = get_table_regenerator("sparsedossa2-stool")
    assert regen.mode == "sparsedossa2-stool"


def test_apply_extra_flags_sparsedossa2_parallel():
    cfg = apply_extra_flags(
        {
            "extra_flags": "--template Stool --workers 4 --cv-folds 2 --maxit 8 --verbose --fit"
        }
    )
    assert cfg["sparsedossa2_template"] == "Stool"
    assert cfg["sparsedossa2_workers"] == 4
    assert cfg["cv_folds"] == 2
    assert cfg["verbose"] is True
    assert cfg["sparsedossa2_fit"] is True
    assert cfg["maxit"] == 8


def test_remap_simulated_features_rank_and_taxid():
    sim = pd.DataFrame({"S1": [10.0, 1.0], "S2": [8.0, 0.0]}, index=["Feature1", "Feature2"])
    out = remap_simulated_features(sim, ["562", "9606"])
    assert list(out.index) == ["562", "9606"]
    assert out.loc["562", "S1"] == 10.0
    keep = pd.DataFrame({"S1": [3.0, 4.0]}, index=["562", "9606"])
    kept = remap_simulated_features(keep, ["562", "9606"])
    assert kept.loc["562", "S1"] == 3.0


def test_simulate_count_matrix_mocked(monkeypatch, tmp_path):
    observed = pd.DataFrame({"a": [10.0, 1.0], "b": [8.0, 2.0]}, index=["562", "9606"])

    def fake_run(argv, timeout=None, **_kwargs):
        out = Path(argv[argv.index("--output") + 1])
        pd.DataFrame({"Sample1": [9.0, 2.0], "Sample2": [7.0, 1.0]}, index=["562", "9606"]).to_csv(
            out
        )
        assert "--workers" in argv
        return ""

    monkeypatch.setattr("samovar.sparsedossa2.run_sparsedossa2_r", fake_run)
    sim = simulate_count_matrix(
        observed, mode="sparsedossa2-fit", n_sample=2, config={"seed": 1, "cores": 2}
    )
    assert list(sim.index) == ["562", "9606"]
    assert sim.shape[1] == 2


def test_sparsedossa2_cv_scorer_mocked(monkeypatch):
    def fake_fitcv(matrix, config=None):
        return {"cv_goodness_of_fit": 12.5, "K": 2, "n_sample": matrix.shape[1]}

    monkeypatch.setattr("samovar.sparsedossa2.fitcv_score_matrix", fake_fitcv)
    tables = {
        "kaiju": pd.DataFrame({"taxid": ["562", "9606"], "N_1": [10, 1], "N_2": [8, 2]})
    }
    row = score_generated_tables(pd.DataFrame(), tables, "sparsedossa2-cv")
    assert row["scorer"] == "sparsedossa2_cv"
    assert row["cv_goodness_of_fit"] == 12.5
    assert row["rank_value"] == pytest.approx(-12.5)
    assert canonicalize_table_scorer("sd2_cv") == "sparsedossa2_cv"


def test_rank_methods_per_annotator_independent_cv(monkeypatch):
    from samovar.table_scorers import rank_methods_per_annotator

    def fake_jobs(items, config=None):
        _ = config
        out = {}
        for key, matrix in items:
            out[key] = {
                "cv_goodness_of_fit": float(matrix.to_numpy().sum()),
                "K": 2,
                "n_sample": matrix.shape[1],
            }
        return out

    monkeypatch.setattr("samovar.sparsedossa2.fitcv_score_jobs", fake_jobs)
    kaiju_good = pd.DataFrame({"taxid": ["562"], "N_1": [20], "N_2": [20]})
    kaiju_bad = pd.DataFrame({"taxid": ["562"], "N_1": [1], "N_2": [1]})
    kraken_good = pd.DataFrame({"taxid": ["9606"], "N_1": [1], "N_2": [1]})
    kraken_bad = pd.DataFrame({"taxid": ["9606"], "N_1": [30], "N_2": [30]})
    tables = {
        "good": {"kaiju": kaiju_good, "kraken2": kraken_good},
        "bad": {"kaiju": kaiju_bad, "kraken2": kraken_bad},
    }
    ranked = rank_methods_per_annotator(
        pd.DataFrame(), tables, "sparsedossa2_cv", modes=["good", "bad"]
    )
    assert ranked["winner_by_annotator"]["kaiju"] == "good"
    assert ranked["winner_by_annotator"]["kraken2"] == "bad"
    assert ranked["winner"] == "mixed"
    assert "cv_matrix" in ranked["by_annotator"]["kaiju"]
    assert "score_matrix" in ranked["by_annotator"]["kaiju"]
    assert ranked["metric_name"] == "cv_goodness_of_fit"
    row = next(r for r in ranked["candidates"] if r["annotator"] == "kaiju" and r["mode"] == "good")
    assert row["metric_name"] == "cv_goodness_of_fit"
    assert "pvalue" in row
    assert "n_observed" in row
    assert "n_generated" in row


def test_sd2_regenerate_from_phyloseq_abundance(monkeypatch, tmp_path):
    from samovar.parse_annotators import Annotation
    from samovar.sparsedossa2 import regenerate

    phyloseq = pd.DataFrame(
        {"taxid": ["562", "9606", "4932"], "CL3": [10, 1, 0], "CC1": [8, 2, 4], "SV1": [1, 3, 5]}
    )

    def fake_sim(matrix, mode="sparsedossa2-fit", n_sample=None, config=None):
        assert list(matrix.index) == ["562", "9606", "4932"]
        assert matrix.shape[1] == 3
        n = int(n_sample or 2)
        out = pd.DataFrame(index=matrix.index)
        for i in range(n):
            out[f"s{i}"] = matrix.iloc[:, i % matrix.shape[1]].to_numpy()
        return out

    monkeypatch.setattr("samovar.sparsedossa2.simulate_count_matrix", fake_sim)
    tables = regenerate(
        Annotation.from_abundance_tables({"table": phyloseq}),
        None,
        {"regeneration_mode": "sparsedossa2-fit", "N": 2, "N_reads": 40},
    )
    table = tables["table"]
    assert list(table.columns)[0] == "taxid"
    assert [c for c in table.columns if str(c).startswith("N_")]
    assert set(table["taxid"].astype(str)) <= {"562", "9606", "4932"}


def test_sd2_cv_accepts_phyloseq_columns(monkeypatch):
    phyloseq = pd.DataFrame({"taxid": ["562", "9606"], "CL3": [10, 1], "CC1": [8, 2]})

    def fake_fitcv(matrix, config=None):
        assert matrix.shape == (2, 2)
        assert list(matrix.index) == ["562", "9606"]
        return {"cv_goodness_of_fit": 3.0, "n_sample": matrix.shape[1]}

    monkeypatch.setattr("samovar.sparsedossa2.fitcv_score_matrix", fake_fitcv)
    row = score_generated_tables(pd.DataFrame(), {"gp": phyloseq}, "sparsedossa2-cv")
    assert row["cv_goodness_of_fit"] == 3.0
    assert row["n_generated"] == 2


def test_write_feature_matrix_keeps_taxa_rows_when_fewer_than_samples(tmp_path):
    from samovar.sparsedossa2 import _write_feature_matrix

    mat = pd.DataFrame(
        {"s1": [1, 2], "s2": [3, 4], "s3": [5, 6], "s4": [7, 8], "s5": [9, 9]},
        index=["562", "9606"],
    )
    dest = tmp_path / "obs.csv"
    _write_feature_matrix(mat, dest)
    loaded = pd.read_csv(dest, index_col=0)
    assert loaded.shape == (2, 5)
    assert list(loaded.index.astype(str)) == ["562", "9606"]

