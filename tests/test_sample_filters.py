"""Sample-level abundance filtering (baselines, prepare, pipeline stages)."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import pandas as pd
import pytest

from samovar.abundance import (
    materialize_observed_abundance,
    n_sample_columns,
    regenerated_abundance_dir,
)
from samovar.config import PipelineConfig
from samovar.regenerate import stage_regenerate_tables
from samovar.sample_filters import (
    KEEP_ONE_WARNING,
    canonicalize_sample_filter,
    filter_samples,
    parse_sample_filter_stage,
    stage_filter_sample_qc,
    unsupported_param_message,
    validate_filter_params_at_prepare,
    validate_pipeline_sample_filters,
)
from samovar.table_scorers import load_tables_by_mode_from_run, stage_score_regenerated_tables


def _table(n_samples=8, n_taxa=4):
    data = {"taxid": [str(560 + i) for i in range(n_taxa)]}
    for i in range(n_samples):
        data[f"N_s{i+1}"] = [float((i + 1) * (j + 1)) for j in range(n_taxa)]
    return pd.DataFrame(data)


def test_canonicalize_filters():
    assert canonicalize_sample_filter("top_n") == "top"
    assert canonicalize_sample_filter("median-sd") == "median_sd"
    assert parse_sample_filter_stage("before_compare") == "candidates"
    assert parse_sample_filter_stage("both") == "both"


def test_top_n_keeps_highest_quality():
    table = _table(6)
    scores = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(1, 7)],
            "quality": [0.1, 0.9, 0.2, 0.8, 0.3, 0.7],
        }
    )
    out = filter_samples(table, scores, {"n": 3}, filt="top")
    kept = [c[2:] for c in n_sample_columns(out)]
    assert kept == ["s2", "s4", "s6"]


def test_top_frac_alone():
    table = _table(10)
    scores = pd.DataFrame(
        {"sample": [f"s{i}" for i in range(1, 11)], "quality": list(range(10))}
    )
    out = filter_samples(table, scores, {"frac": 0.4}, filt="top")
    assert len(n_sample_columns(out)) == 4


def test_top_n_and_frac_together():
    table = _table(10)
    scores = pd.DataFrame(
        {"sample": [f"s{i}" for i in range(1, 11)], "quality": list(range(10))}
    )
    out = filter_samples(table, scores, {"n": 8, "frac": 0.3}, filt="top")
    assert len(n_sample_columns(out)) == 3


def test_median_sd_keeps_band():
    table = _table(7)
    scores = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(1, 8)],
            "quality": [0.0, 0.49, 0.5, 0.51, 0.5, 0.5, 1.0],
        }
    )
    out = filter_samples(table, scores, {"sd": 1.0}, filt="median_sd")
    kept = n_sample_columns(out)
    assert "N_s1" not in kept
    assert "N_s7" not in kept
    assert len(kept) >= 2


def test_prepare_rejects_unsupported_param():
    msg = unsupported_param_message("median_sd", ["n"])
    assert msg.startswith("Not supported by the exact abundance table filtering function.")
    assert "--sample-filter-sd" in msg
    with pytest.raises(ValueError, match="Not supported by the exact abundance table filtering function"):
        validate_filter_params_at_prepare("median_sd", {"n": 3})
    with pytest.raises(ValueError, match="Not supported by the exact abundance table filtering function"):
        validate_filter_params_at_prepare("top", {"sd": 2})


def test_prepare_rejects_bad_n_and_frac():
    with pytest.raises(ValueError, match="top N need to keep more than 1 sample"):
        validate_filter_params_at_prepare("top", {"n": 1})
    with pytest.raises(ValueError, match="top N need to be smaller than the amount of the generated samples"):
        validate_filter_params_at_prepare("top", {"n": 8}, n_generated=8)
    with pytest.raises(ValueError, match="top % need to be more than 0 and less than 1"):
        validate_filter_params_at_prepare("top", {"frac": 1.0})
    with pytest.raises(ValueError, match="top % need to be more than 0 and less than 1"):
        validate_filter_params_at_prepare("top", {"frac": 0.0})
    with pytest.raises(ValueError, match="top % need to be more than 0 and less than 1"):
        validate_filter_params_at_prepare("top", {"frac": 50})


def test_prepare_requires_sample_score():
    with pytest.raises(ValueError, match="regenerated sample quality function"):
        validate_pipeline_sample_filters(sample_filter="top", params={"n": 3})


def test_prepare_cli_rejects_unsupported_combo():
    args = argparse.Namespace(
        input_config=None,
        input_dir="/tmp/reads",
        output_dir="/tmp/out",
        sample_score=["bray_curtis"],
        sample_filter=["median_sd"],
        sample_filter_n=["4"],
        kraken2=None,
        kaiju=None,
    )
    with pytest.raises(ValueError, match="Not supported by the exact abundance table filtering function"):
        PipelineConfig.from_args(args)


def test_runtime_top_n_bounds():
    table = _table(4)
    scores = pd.DataFrame({"sample": ["s1", "s2", "s3", "s4"], "quality": [1, 2, 3, 4]})
    with pytest.raises(ValueError, match="top N need to keep more than 1 sample"):
        filter_samples(table, scores, {"n": 1}, filt="top")
    with pytest.raises(ValueError, match="top N need to be smaller than the amount of the generated samples"):
        filter_samples(table, scores, {"n": 4}, filt="top")


def test_keep_one_warning_when_frac_too_small():
    table = _table(5)
    scores = pd.DataFrame({"sample": [f"s{i}" for i in range(1, 6)], "quality": [5, 4, 3, 2, 1]})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = filter_samples(table, scores, {"frac": 0.1}, filt="top")
    assert len(n_sample_columns(out)) == 1
    assert any(KEEP_ONE_WARNING in str(item.message) for item in caught)


def test_median_sd_all_equal_keeps_all():
    table = _table(5)
    scores = pd.DataFrame({"sample": [f"s{i}" for i in range(1, 6)], "quality": [0.5] * 5})
    out = filter_samples(table, scores, {"sd": 2.0}, filt="median_sd")
    assert len(n_sample_columns(out)) == 5


def _pipeline_text(cfg, dest: Path) -> str:
    dest.mkdir(parents=True, exist_ok=True)
    return Path(cfg.generate_pipeline(str(dest))).read_text()


def test_prepare_absent_filter_unchanged(tmp_path):
    (tmp_path / "reads").mkdir()
    cfg = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=None,
            input_dir=str(tmp_path / "reads"),
            output_dir=str(tmp_path / "out"),
            sample_score=["bray_curtis"],
            kraken2=[["kraken2 /tmp/k2"]],
            kaiju=None,
        )
    )
    script = _pipeline_text(cfg, tmp_path / "out")
    assert "sample_filters" not in script
    assert "filter_sample_qc_full" not in script
    assert "filter_sample_qc_final" not in script


def test_prepare_filter_stages(tmp_path):
    (tmp_path / "reads").mkdir()
    final_only = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=None,
            input_dir=str(tmp_path / "reads"),
            output_dir=str(tmp_path / "final"),
            sample_score=["bray_curtis"],
            sample_filter=["top"],
            sample_filter_n=["3"],
            sample_filter_stage="final",
            regeneration_n=8,
            table_reads_generator=["bootstrap"],
            kraken2=None,
            kaiju=None,
        )
    )
    text = _pipeline_text(final_only, tmp_path / "final")
    assert "filter_sample_qc_final" in text
    assert "filter_sample_qc_full" not in text
    assert "sample_filters stage" in text
    before = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=None,
            input_dir=str(tmp_path / "reads"),
            output_dir=str(tmp_path / "before"),
            sample_score=["bray_curtis"],
            sample_filter=["top"],
            sample_filter_n=["3"],
            sample_filter_stage="candidates",
            regeneration_n=8,
            table_reads_generator=["bootstrap", "vae"],
            kraken2=None,
            kaiju=None,
        )
    )
    text_b = _pipeline_text(before, tmp_path / "before")
    assert "filter_sample_qc_full" in text_b
    assert "filter_sample_qc_final" not in text_b
    both = PipelineConfig.from_args(
        argparse.Namespace(
            input_config=None,
            input_dir=str(tmp_path / "reads"),
            output_dir=str(tmp_path / "both"),
            sample_score=["bootstrap:bray_curtis"],
            sample_filter=["bootstrap:top", "vae:median_sd"],
            sample_filter_n=["bootstrap:3"],
            sample_filter_sd=["vae:2"],
            sample_filter_stage="both",
            regeneration_n=8,
            table_reads_generator=["bootstrap", "vae"],
            kraken2=None,
            kaiju=None,
        )
    )
    assert both.sample_score_by_method["bootstrap"] == "bray_curtis"
    assert both.sample_filter_by_method["bootstrap"] == "top"
    assert both.sample_filter_by_method["vae"] == "median_sd"
    text_both = _pipeline_text(both, tmp_path / "both")
    assert "filter_sample_qc_full" in text_both
    assert "filter_sample_qc_final" in text_both


def _seed_abundance(tmp_path, n_samples=8):
    src = tmp_path / "obs.csv"
    frame = pd.DataFrame(
        {
            "taxid": ["562", "9606", "28901", "1280"],
            **{f"CL{i}": [10 + i, 3, 2, 1] for i in range(n_samples)},
        }
    )
    frame.to_csv(src, index=False)
    materialize_observed_abundance(tmp_path)
    return tmp_path


def test_pipeline_filter_final_bootstrap_vae_glm(tmp_path):
    _seed_abundance(tmp_path, n_samples=8)
    cfg = {
        "table_reads_generators": ["bootstrap", "vae", "glm"],
        "table_reads_generator": ["bootstrap", "vae", "glm"],
        "N": 8,
        "N_reads": 200,
        "seed": 1,
        "sample_score": "bray_curtis",
        "sample_filter": "top",
        "n": 3,
        "sample_filter_stage": "final",
        "table_score": "shannon_ks",
    }
    tables = stage_regenerate_tables(tmp_path, cfg)
    assert tables
    modes = load_tables_by_mode_from_run(tmp_path)
    assert set(modes) >= {"bootstrap", "vae", "glm"}
    assert "camisim-table" not in modes
    ranked = stage_score_regenerated_tables(tmp_path, cfg)
    assert ranked
    before_n = len(n_sample_columns(next(iter(ranked["tables"].values()))))
    result = stage_filter_sample_qc(tmp_path, cfg, phase="final")
    assert result["enabled"]
    after = list(regenerated_abundance_dir(tmp_path).glob("*.csv"))
    assert after
    kept = n_sample_columns(pd.read_csv(after[0]))
    assert len(kept) == 3
    assert before_n >= 3
    assert (tmp_path / "sampleFilter" / "final" / "report.csv").is_file()


def test_pipeline_filter_before_compare(tmp_path):
    _seed_abundance(tmp_path, n_samples=8)
    cfg = {
        "table_reads_generators": ["direct", "bootstrap", "vae"],
        "N": 8,
        "N_reads": 200,
        "seed": 2,
        "sample_score": "bray_curtis",
        "sample_filter": "top",
        "n": 4,
        "sample_filter_stage": "candidates",
        "table_score": "shannon_ks",
    }
    stage_regenerate_tables(tmp_path, cfg)
    full = stage_filter_sample_qc(tmp_path, cfg, phase="full")
    assert full["enabled"]
    modes = load_tables_by_mode_from_run(tmp_path)
    boot = next(iter(modes["bootstrap"].values()))
    assert len(n_sample_columns(boot)) == 4
    ranked = stage_score_regenerated_tables(tmp_path, cfg)
    assert ranked.get("tables")


def test_pipeline_filter_both_modes(tmp_path):
    _seed_abundance(tmp_path, n_samples=8)
    cfg = {
        "table_reads_generators": ["bootstrap", "glm"],
        "N": 8,
        "N_reads": 150,
        "seed": 3,
        "sample_score": "bray_curtis",
        "sample_filter": "median_sd",
        "sd": 2.0,
        "sample_filter_stage": "both",
        "table_score": "shannon_ks",
    }
    stage_regenerate_tables(tmp_path, cfg)
    assert stage_filter_sample_qc(tmp_path, cfg, phase="full")["enabled"]
    stage_score_regenerated_tables(tmp_path, cfg)
    assert stage_filter_sample_qc(tmp_path, cfg, phase="final")["enabled"]
    assert (tmp_path / "sampleFilter" / "full" / "report.csv").is_file()
    assert (tmp_path / "sampleFilter" / "final" / "report.csv").is_file()


def test_independent_filter_cli(tmp_path):
    from samovar.sample_filters import main as filter_main

    table = tmp_path / "t.csv"
    scores = tmp_path / "s.csv"
    out = tmp_path / "o.csv"
    _table(6).to_csv(table, index=False)
    pd.DataFrame(
        {"sample": [f"s{i}" for i in range(1, 7)], "quality": [1, 2, 3, 4, 5, 6]}
    ).to_csv(scores, index=False)
    rc = filter_main(
        [
            "filter",
            "--table",
            str(table),
            "--scores",
            str(scores),
            "-o",
            str(out),
            "--filter",
            "top",
            "--n",
            "3",
        ]
    )
    assert rc == 0
    assert len(n_sample_columns(pd.read_csv(out))) == 3
