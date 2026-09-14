"""Sample-level abundance quality scoring (contract, prepare, pipeline stages)."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest
import yaml

from samovar.abundance import materialize_observed_abundance, regenerated_abundance_dir, write_abundance_dir
from samovar.config import PipelineConfig, setup_pipeline
from samovar.paths import repo_root, write_config
from samovar.regenerate import stage_regenerate_tables
from samovar.sample_scorers import (
    bray_curtis_distance,
    canonicalize_sample_scorer,
    main as sample_scorers_main,
    parse_sample_score_tokens,
    resolve_sample_scorer,
    sample_qc_configured,
    score_samples,
    score_samples_bray_curtis,
    stage_score_sample_qc,
)
from samovar.table_scorers import stage_score_regenerated_tables
from samovar.tools_import import import_tool


REPO = repo_root()
SAMPLE_QC_ROOT = REPO / "samovar" / "sampleQC"


def _tiny(cols=None, extra_refs=None):
    data = {"taxid": ["562", "2886930", "9606"], "N_s1": [10, 4, 2], "N_s2": [8, 3, 1]}
    if extra_refs:
        data.update(extra_refs)
    if cols:
        data = {"taxid": data["taxid"], **cols}
    return pd.DataFrame(data)


def _publish_sample_qc(src: Path) -> None:
    SAMPLE_QC_ROOT.mkdir(parents=True, exist_ok=True)
    for phase in ("full", "final"):
        dest = SAMPLE_QC_ROOT / phase
        dest.mkdir(parents=True, exist_ok=True)
        incoming = src / "sampleQC" / phase
        if not incoming.is_dir():
            continue
        for path in incoming.rglob("*.csv"):
            target = dest / path.relative_to(incoming)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)


def test_bray_curtis_higher_is_better():
    ref = _tiny(extra_refs={"N_s3": [9, 4, 2], "N_s4": [7, 3, 1]})
    close = _tiny(cols={"N_x": [10, 4, 2]})
    far = _tiny(cols={"N_x": [1, 20, 30]})
    q_close = score_samples_bray_curtis(close, ref)["quality"].iloc[0]
    q_far = score_samples_bray_curtis(far, ref)["quality"].iloc[0]
    assert q_close > q_far
    assert 0.0 <= q_far <= 1.0
    assert 0.0 <= q_close <= 1.0


def test_nearest_three_mean_distance():
    ref = pd.DataFrame(
        {
            "taxid": ["a", "b"],
            "N_r1": [1.0, 0.0],
            "N_r2": [0.0, 1.0],
            "N_r3": [0.5, 0.5],
            "N_r4": [0.9, 0.1],
        }
    )
    gen = pd.DataFrame({"taxid": ["a", "b"], "N_g": [1.0, 0.0]})
    frame = score_samples_bray_curtis(gen, ref)
    vec = gen.set_index("taxid")["N_g"].to_numpy(dtype=float)
    dists = sorted(
        bray_curtis_distance(vec, ref.set_index("taxid")[c].to_numpy(dtype=float))
        for c in ["N_r1", "N_r2", "N_r3", "N_r4"]
    )
    expected = 1.0 - sum(dists[:3]) / 3.0
    assert frame["n_neighbors"].iloc[0] == 3
    assert frame["quality"].iloc[0] == pytest.approx(expected)
    assert score_samples(gen, ref, scorer="bray_curtis")["quality"].iloc[0] == pytest.approx(
        expected
    )


def test_fewer_than_three_reference_samples():
    ref = _tiny(cols={"N_only": [10, 4, 2], "N_two": [8, 3, 1]})
    gen = _tiny(cols={"N_g": [10, 4, 2]})
    frame = score_samples_bray_curtis(gen, ref)
    assert int(frame["n_neighbors"].iloc[0]) == 2
    assert frame["ok"].iloc[0]


def test_zero_and_invalid_samples():
    ref = _tiny()
    gen = pd.DataFrame(
        {
            "taxid": ["562", "2886930", "9606"],
            "N_zero": [0, 0, 0],
            "N_ok": [10, 4, 2],
        }
    )
    frame = score_samples_bray_curtis(gen, ref).set_index("sample")
    assert frame.loc["zero", "ok"]
    assert frame.loc["ok", "quality"] >= frame.loc["zero", "quality"]
    empty_ref = score_samples_bray_curtis(gen, pd.DataFrame({"taxid": ["562"]}))
    assert list(empty_ref["n_neighbors"]) == [0, 0]
    assert list(empty_ref["quality"]) == [0.0, 0.0]


def test_identical_sample_scores_one():
    table = _tiny(cols={"N_only": [10, 4, 2]})
    frame = score_samples_bray_curtis(table, table)
    assert frame["quality"].iloc[0] == pytest.approx(1.0)
    assert int(frame["n_neighbors"].iloc[0]) == 1


def test_independent_invocation_cli(tmp_path):
    gen = tmp_path / "gen.csv"
    ref = tmp_path / "ref.csv"
    out = tmp_path / "scores.csv"
    _tiny().to_csv(ref, index=False)
    _tiny(cols={"N_g": [9, 4, 2]}).to_csv(gen, index=False)
    rc = sample_scorers_main(
        ["score", "--generated", str(gen), "--reference", str(ref), "-o", str(out)]
    )
    assert rc == 0
    frame = pd.read_csv(out)
    assert "quality" in frame.columns
    assert len(frame) == 1


def test_resolve_precedence():
    assert resolve_sample_scorer("kaiju", global_name="bray_curtis") == "bray_curtis"
    assert (
        resolve_sample_scorer(
            "kaiju",
            global_name="bray_curtis",
            by_annotator={"kaiju": "bray-curtis"},
        )
        == "bray_curtis"
    )
    assert canonicalize_sample_scorer("nearest3_bray") == "bray_curtis"
    g, by = parse_sample_score_tokens(["bray_curtis", "kaiju:bray_curtis"])[:2]
    assert g == "bray_curtis"
    assert by["kaiju"] == "bray_curtis"
    assert sample_qc_configured("", {}) is False
    assert sample_qc_configured("bray_curtis", {}) is True


def test_prepare_absent_scorer_unchanged(tmp_path):
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(tmp_path / "reads"),
        output_dir=str(tmp_path / "out"),
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=[["kaiju /tmp/kaiju"]],
    )
    (tmp_path / "reads").mkdir()
    result = setup_pipeline(args)
    script = Path(result["pipeline"]).read_text()
    a2iss = yaml.safe_load(Path(result["configs"]["annotation2iss"]).read_text())
    assert "sample_scorers" not in script
    assert "score_sample_qc_full" not in script
    assert "score_sample_qc_final" not in script
    assert not (a2iss.get("sample_score") or "").strip()
    assert not a2iss.get("sample_score_by_annotator")


def test_prepare_global_and_annotator_specific(tmp_path):
    (tmp_path / "reads").mkdir()
    yaml_path = tmp_path / "in.yaml"
    yaml_path.write_text(
        yaml.dump(
            {
                "input_dir": str(tmp_path / "reads"),
                "annotators": [],
                "sample_score": "bray_curtis",
                "sample_score_by_annotator": {"kaiju": "bray_curtis"},
            }
        )
    )
    global_cfg = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=None,
            input_dir=str(tmp_path / "reads"),
            output_dir=str(tmp_path / "global"),
            sample_score=["bray_curtis"],
            kraken2=None,
            kaiju=None,
        )
    )
    assert global_cfg.sample_score == "bray_curtis"
    script = Path(global_cfg.generate_pipeline(str(tmp_path / "global"))).read_text()
    assert "score_sample_qc_full" in script
    assert "score_sample_qc_final" in script
    assert "--phase full" in script
    mixed = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=str(yaml_path),
            input_dir=None,
            output_dir=str(tmp_path / "mixed"),
            sample_score=["kraken2:bray_curtis"],
            kraken2=None,
            kaiju=None,
        )
    )
    assert mixed.sample_score == "bray_curtis"
    assert mixed.sample_score_by_annotator["kaiju"] == "bray_curtis"
    assert mixed.sample_score_by_annotator["kraken2"] == "bray_curtis"
    off = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=str(yaml_path),
            input_dir=None,
            output_dir=str(tmp_path / "off"),
            sample_score=["none"],
            kraken2=None,
            kaiju=None,
        )
    )
    assert not sample_qc_configured(off.sample_score, off.sample_score_by_annotator)
    specific_only = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=None,
            input_dir=str(tmp_path / "reads"),
            output_dir=str(tmp_path / "ann"),
            sample_score=["kaiju:bray_curtis"],
            kraken2=None,
            kaiju=None,
        )
    )
    assert not specific_only.sample_score
    assert specific_only.sample_score_by_annotator["kaiju"] == "bray_curtis"
    assert (
        resolve_sample_scorer(
            "kraken2",
            global_name=specific_only.sample_score,
            by_annotator=specific_only.sample_score_by_annotator,
        )
        == ""
    )
    assert (
        resolve_sample_scorer(
            "kaiju",
            global_name=specific_only.sample_score,
            by_annotator=specific_only.sample_score_by_annotator,
        )
        == "bray_curtis"
    )


def test_tools_import_sample_scoring(tmp_path, monkeypatch):
    from samovar.main_config import normalize_tool_group, parse_tool_entry

    assert normalize_tool_group("sample-score") == "sample_scoring"
    script = Path("tests/tools/dummy_sample_scorer.py").resolve()
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="dummy_sample",
        tool_type="sample-qc",
        exec_path=str(script),
        also_repo_build=False,
    )
    assert spec[3] == "sample_scoring"
    raw = json.loads(cfg.read_text())["tools"]["dummy_sample"]
    assert parse_tool_entry(raw, "dummy_sample")[3] == "sample_scoring"
    gen = _tiny()
    frame = score_samples(gen, gen, scorer="dummy_sample")
    assert (frame["quality"] == 1.0).all()


def test_pipeline_full_and_final_stages(tmp_path):
    src = tmp_path / "gp.csv"
    pd.DataFrame({"taxid": [562, 9606], "CL3": [10, 1], "CC1": [8, 2]}).to_csv(src, index=False)
    materialize_observed_abundance(tmp_path)
    tables = stage_regenerate_tables(
        tmp_path, {"regeneration_mode": "direct", "table_reads_generator": "direct"}
    )
    assert tables
    cfg = {
        "sample_score": "bray_curtis",
        "sample_score_by_annotator": {},
    }
    full = stage_score_sample_qc(tmp_path, cfg, phase="full")
    assert full["enabled"]
    assert (tmp_path / "sampleQC" / "full" / "all.csv").is_file()
    ranked = stage_score_regenerated_tables(tmp_path, {"table_score": "shannon_ks"})
    assert ranked
    final = stage_score_sample_qc(tmp_path, cfg, phase="final")
    assert final["enabled"]
    assert (tmp_path / "sampleQC" / "final" / "all.csv").is_file()
    scores = pd.read_csv(tmp_path / "sampleQC" / "final" / "all.csv")
    assert "sample" in scores.columns
    assert "quality" in scores.columns
    assert scores["sample"].nunique() >= 1
    disabled = stage_score_sample_qc(tmp_path, {}, phase="full")
    assert disabled["enabled"] is False
    _publish_sample_qc(tmp_path)
    assert (SAMPLE_QC_ROOT / "full" / "all.csv").is_file()
    assert (SAMPLE_QC_ROOT / "final" / "all.csv").is_file()


def test_annotator_specific_stage_skips_unmapped(tmp_path):
    obs = tmp_path / "initial_abundance"
    regen = regenerated_abundance_dir(tmp_path)
    obs.mkdir(parents=True)
    write_abundance_dir(
        obs,
        {
            "kaiju": _tiny(),
            "kraken2": _tiny(cols={"N_s1": [1, 2, 3], "N_s2": [3, 2, 1]}),
        },
    )
    write_abundance_dir(
        regen,
        {
            "kaiju": _tiny(cols={"N_s1": [10, 4, 2], "N_s2": [8, 3, 1]}),
            "kraken2": _tiny(cols={"N_s1": [1, 2, 3], "N_s2": [3, 2, 1]}),
        },
    )
    result = stage_score_sample_qc(
        tmp_path,
        {"sample_score": "", "sample_score_by_annotator": {"kaiju": "bray_curtis"}},
        phase="final",
    )
    written = [Path(p).name for p in result["written"]]
    assert "kaiju.csv" in written
    assert "kraken2.csv" not in written
    dest = SAMPLE_QC_ROOT / "final"
    dest.mkdir(parents=True, exist_ok=True)
    src = tmp_path / "sampleQC" / "final" / "kaiju.csv"
    if src.is_file():
        shutil.copy2(src, dest / "kaiju.csv")
