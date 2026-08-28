"""Shannon KS scoring for table_reads_generator selection."""

import numpy as np
import pandas as pd
import pytest

from pathlib import Path

from samovar.annotation_io import read_annotation_dir
from samovar.table_scorers import (
    canonicalize_table_scorer,
    ks_bray,
    ks_shannon,
    pick_best_table_method,
    rank_methods_per_annotator,
    score_generated_tables,
    shannon_entropy,
    shannon_vector_from_annotation,
)

GLOBALPATTERNS = Path("tests/data/globalpatterns_abundance.csv")
BRAY_PLUGIN = Path("tests/data/bray_ks_table_scorer.py")


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


def test_shannon_entropy_uniform_two_taxa():
    value = shannon_entropy([1.0, 1.0])
    assert abs(value - np.log(2)) < 1e-9


def test_ks_identical_vectors_is_zero():
    vec = np.array([0.1, 0.4, 0.7])
    row = ks_shannon(vec, vec)
    assert row["ks_statistic"] == 0.0
    assert row["pvalue"] == 1.0
    assert row["rank_value"] == 0.0
    assert row["metric"] == 0.0
    assert row["metric_name"] == "ks_statistic"
    assert "n_observed" in row and "n_generated" in row


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
    assert score["metric_name"] == "ks_statistic"
    assert canonicalize_table_scorer("alpha-ks") == "shannon_ks"


def test_rank_methods_shannon_has_score_matrix(toy_annotation_dir, tmp_path):
    from samovar.regenerate import regenerate_annotation_tables
    from samovar.table_scorers import CANDIDATE_KEYS, rank_methods_per_annotator

    data = read_annotation_dir(toy_annotation_dir)
    direct = regenerate_annotation_tables(
        toy_annotation_dir, tmp_path / "direct", {"regeneration_mode": "direct", "N": 2, "seed": 1}
    )
    bootstrap = regenerate_annotation_tables(
        toy_annotation_dir,
        tmp_path / "boot",
        {"regeneration_mode": "bootstrap", "N": 2, "seed": 1},
    )
    ranked = rank_methods_per_annotator(
        data,
        {"direct": direct, "bootstrap": bootstrap},
        "shannon_ks",
        modes=["direct", "bootstrap"],
    )
    assert ranked["metric_name"] == "ks_statistic"
    block = ranked["by_annotator"]["kaiju"]
    assert "score_matrix" in block
    assert "direct" in block["score_matrix"]["direct"]
    row = next(r for r in ranked["candidates"] if r["mode"] == "direct" and r["annotator"] == "kaiju")
    for key in CANDIDATE_KEYS:
        assert key in row
    assert row["ks_statistic"] == 0.0


def test_ks_bray_identical_tables_is_zero():
    table = pd.DataFrame(
        {
            "taxid": [562, 9606],
            "N_a": [10, 2],
            "N_b": [8, 4],
            "N_c": [1, 9],
        }
    )
    row = ks_bray(table, table.copy())
    assert row["scorer"] == "bray_ks"
    assert row["ks_statistic"] == 0.0
    assert row["rank_value"] == 0.0
    assert canonicalize_table_scorer("bray") == "bray_ks"


def test_annotation_scoring_is_not_table_scoring(tmp_path, monkeypatch):
    from samovar.paths import write_config
    from samovar.tools_import import import_tool

    script = tmp_path / "viz_score.py"
    script.write_text("def score(inputs, output_dir, config):\n    return None\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="viz_score",
        tool_type="scoring",
        exec_path=str(script),
        also_repo_build=False,
    )
    with pytest.raises(ValueError, match="table-scoring"):
        canonicalize_table_scorer("viz_score")


def test_globalpatterns_three_generators_three_scorers(tmp_path, monkeypatch):
    from samovar.abundance import n_sample_columns
    from samovar.paths import write_config
    from samovar.regenerate import regenerate_annotation_tables
    from samovar.table_scorers import CANDIDATE_KEYS
    from samovar.tools_import import import_tool

    if not GLOBALPATTERNS.is_file():
        pytest.skip("tests/data/globalpatterns_abundance.csv missing")
    gp = pd.read_csv(GLOBALPATTERNS)
    sample_cols = [c for c in gp.columns if c != "taxid"][:6]
    gp = gp[["taxid", *sample_cols]].head(12)
    src = tmp_path / "gp.csv"
    gp.to_csv(src, index=False)
    tables_by_mode = {}
    for mode in ("direct", "bootstrap", "glm"):
        tables_by_mode[mode] = regenerate_annotation_tables(
            src,
            tmp_path / mode,
            {"regeneration_mode": mode, "N": 6, "seed": 1, "rescale_abundance": True},
        )
        assert n_sample_columns(next(iter(tables_by_mode[mode].values())))
    annotator = next(iter(tables_by_mode["direct"]))

    tables_by_mode = {}
    for mode in ("direct", "bootstrap", "glm"):
        tables_by_mode[mode] = regenerate_annotation_tables(
            src,
            tmp_path / mode,
            {"regeneration_mode": mode, "N": 6, "seed": 1, "rescale_abundance": True},
        )
        assert n_sample_columns(next(iter(tables_by_mode[mode].values())))

    shannon = rank_methods_per_annotator(
        gp, tables_by_mode, "shannon_ks", modes=["direct", "bootstrap", "glm"]
    )
    direct_sh = next(
        r
        for r in shannon["candidates"]
        if r["mode"] == "direct" and r["annotator"] == annotator
    )
    for key in CANDIDATE_KEYS:
        assert key in direct_sh
    assert direct_sh["ks_statistic"] == 0.0
    assert shannon["winner_by_annotator"][annotator] == "direct"

    bray = rank_methods_per_annotator(
        gp, tables_by_mode, "bray_ks", modes=["direct", "bootstrap", "glm"]
    )
    direct_br = next(
        r for r in bray["candidates"] if r["mode"] == "direct" and r["annotator"] == annotator
    )
    assert direct_br["ks_statistic"] == 0.0
    assert bray["scorer"] == "bray_ks"

    def fake_jobs(items, *, config=None):
        out = {}
        for key, matrix in items:
            out[str(key)] = {"cv_goodness_of_fit": float(matrix.to_numpy().sum())}
        return out

    monkeypatch.setattr("samovar.sparsedossa2.fitcv_score_jobs", fake_jobs)
    sd2 = rank_methods_per_annotator(
        gp, tables_by_mode, "sparsedossa2_cv", modes=["direct", "bootstrap", "glm"]
    )
    assert sd2["scorer"] == "sparsedossa2_cv"
    assert sd2["by_annotator"][annotator]["candidates"]
    assert sd2["winner_by_annotator"][annotator] in {"direct", "bootstrap", "glm"}

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    if not BRAY_PLUGIN.is_file():
        pytest.skip("tests/data/bray_ks_table_scorer.py missing")
    import_tool(
        name="bray_ks_plugin",
        tool_type="table-scoring",
        exec_path=str(BRAY_PLUGIN.resolve()),
        also_repo_build=False,
    )
    assert canonicalize_table_scorer("bray_ks_plugin") == "bray_ks_plugin"
    custom = rank_methods_per_annotator(
        gp, tables_by_mode, "bray_ks_plugin", modes=["direct", "bootstrap", "glm"]
    )
    assert custom["scorer"] == "bray_ks_plugin"
    custom_direct = next(
        r
        for r in custom["candidates"]
        if r["mode"] == "direct" and r["annotator"] == annotator
    )
    assert custom_direct["ks_statistic"] == 0.0
    assert custom["winner_by_annotator"][annotator] == bray["winner_by_annotator"][annotator]
