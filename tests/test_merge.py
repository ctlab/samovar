"""samovar merge unions completed runs; sample names must be distinct."""

from __future__ import annotations

import os
import subprocess
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest
import yaml

from samovar.config import setup_pipeline
from samovar.exec_control import CHECKPOINT_STEPS, listed_done, mark_done
from samovar.merge import (
    SAMPLE_NAME_ERROR,
    MergeError,
    main as merge_main,
    merge_runs,
)

REPO = Path(__file__).resolve().parents[1]


def _ns(**kwargs):
    base = {
        "input_config": None,
        "input_dir": None,
        "output_dir": None,
        "add_annotator": False,
        "kraken2": None,
        "kaiju": None,
        "dummy": None,
        "cores": 1,
        "no_multiqc": True,
        "run_multiqc": False,
        "startpoint": None,
        "table_reads_generator": None,
        "use_test_genomes": False,
    }
    base.update(kwargs)
    return Namespace(**base)


def _ann_csv(sample: str, taxid: int = 562) -> str:
    return (
        "seq,taxID_dummy,length,sample,true\n"
        f"{sample}_r0,{taxid},12,{sample},{taxid}\n"
        f"{sample}_r1,{taxid},12,{sample},{taxid}\n"
    )


def _fastq(dest: Path, sample: str) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    rec = f"@{sample}|taxid:562\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
    (dest / f"{sample}_R1.fastq").write_text(rec)
    (dest / f"{sample}_R2.fastq").write_text(rec)


def _seed_run(root: Path, sample: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    reads = root / "reads"
    _fastq(reads, sample)
    out = root / "run"
    setup_pipeline(
        _ns(
            input_dir=str(reads),
            output_dir=str(out),
            dummy=[["dummy"]],
            run_multiqc=False,
        )
    )
    _fastq(out / "initial", sample)
    _fastq(out / "initial_trimmed", sample)
    reports = out / "initial_reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / f"{sample}.out").write_text("ok\n")
    ann = out / "initial_annotations"
    ann.mkdir(parents=True, exist_ok=True)
    (ann / f"{sample}.annotation.csv").write_text(_ann_csv(sample))
    (ann / "combined_annotation_table.csv").write_text(_ann_csv(sample))
    abund = out / "initial_abundance"
    abund.mkdir(parents=True, exist_ok=True)
    (abund / "dummy.csv").write_text(f"taxid,N_{sample}\n562,2\n")
    regen_ab = out / "regenerated" / ".regenerated_abundance"
    regen_ab.mkdir(parents=True, exist_ok=True)
    (regen_ab / "dummy.csv").write_text(f"taxid,N_{sample}\n562,2\n")
    _fastq(out / "regenerated", f"{sample}_k1")
    regen_ann = out / "regenerated_annotations"
    regen_ann.mkdir(parents=True, exist_ok=True)
    (regen_ann / f"{sample}_k1.annotation.csv").write_text(_ann_csv(f"{sample}_k1"))
    (regen_ann / "combined_annotation_table.csv").write_text(_ann_csv(f"{sample}_k1"))
    for step in CHECKPOINT_STEPS:
        mark_done(out, step)
    return out


def test_cli_help_lists_merge():
    cfg = REPO / "build" / "config.json"
    if not cfg.exists():
        cfg.parent.mkdir(parents=True, exist_ok=True)
        cfg.write_text('{"python_path": "python3", "r_path": "R", "r_lib_path": "."}\n')
    result = subprocess.run(
        ["bash", str(REPO / "bin" / "samovar"), "help"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "merge" in result.stdout
    assert "--mode initial" in result.stdout
    assert "--mode regenerated" in result.stdout


def test_merge_rejects_duplicate_sample_names(tmp_path):
    a = _seed_run(tmp_path / "a", "same")
    b = _seed_run(tmp_path / "b", "same")
    with pytest.raises(MergeError, match="Sample names need to be distinct"):
        merge_runs([a, b], tmp_path / "out", "initial")
    rc = merge_main(
        ["--mode", "initial", "--output_dir", str(tmp_path / "cli"), str(a), str(b)]
    )
    assert rc == 1


def test_merge_initial_unions_annotations_redoes_downstream(tmp_path):
    a = _seed_run(tmp_path / "a", "s1")
    b = _seed_run(tmp_path / "b", "s2")
    dest = tmp_path / "merged_initial"
    meta = merge_runs([a, b], dest, "initial")
    assert meta["mode"] == "initial"
    assert meta["next_step"] == "viz_initial"
    assert set(meta["samples"]) >= {"s1", "s2"}
    assert (dest / "initial" / "s1_R1.fastq").is_file()
    assert (dest / "initial" / "s2_R1.fastq").is_file()
    assert (dest / "initial_annotations" / "s1.annotation.csv").is_file()
    assert (dest / "initial_annotations" / "s2.annotation.csv").is_file()
    combined = pd.read_csv(dest / "initial_annotations" / "combined_annotation_table.csv")
    assert set(combined["sample"]) == {"s1", "s2"}
    assert not (dest / "regenerated_annotations").exists()
    assert not (dest / "initial_abundance").exists()
    done = set(listed_done(dest))
    assert "combine_initial" in done
    assert "viz_initial" not in done
    assert "abundance_tables" not in done
    assert "reprofile" not in done
    assert (dest / ".log" / "merge.yaml").is_file()
    assert (dest / ".log" / "samovar.sh").is_file()


def test_merge_regenerated_unions_both_combined_tables(tmp_path):
    a = _seed_run(tmp_path / "a", "s1")
    b = _seed_run(tmp_path / "b", "s2")
    dest = tmp_path / "merged_regen"
    meta = merge_runs([a, b], dest, "regenerated")
    assert meta["mode"] == "regenerated"
    assert meta["next_step"] == "viz_regenerated"
    init = pd.read_csv(dest / "initial_annotations" / "combined_annotation_table.csv")
    regen = pd.read_csv(
        dest / "regenerated_annotations" / "combined_annotation_table.csv"
    )
    assert set(init["sample"]) == {"s1", "s2"}
    assert {"s1_k1", "s2_k1"} <= set(regen["sample"])
    abund = pd.read_csv(dest / "regenerated" / ".regenerated_abundance" / "dummy.csv")
    assert "N_s1" in abund.columns and "N_s2" in abund.columns
    observed = pd.read_csv(dest / "initial_abundance" / "dummy.csv")
    assert "N_s1" in observed.columns and "N_s2" in observed.columns
    done = set(listed_done(dest))
    assert "combine_regenerated" in done
    assert "viz_regenerated" not in done
    assert "reprofile" not in done


def test_merge_regenerated_skips_shared_table_score_plots(tmp_path):
    a = _seed_run(tmp_path / "a", "s1")
    b = _seed_run(tmp_path / "b", "s2")
    for run in (a, b):
        plots = run / "regenerated" / ".regenerated_abundance" / "table_score_plots"
        plots.mkdir(parents=True)
        (plots / "TableScore_kaiju.png").write_bytes(b"\x89PNG\r\n")
        (plots / "TableScore_quality_scores_mqc.json").write_text("{}\n")
    dest = tmp_path / "merged_regen"
    merge_runs([a, b], dest, "regenerated")
    assert not (dest / "regenerated" / ".regenerated_abundance" / "table_score_plots").exists()
    abund = pd.read_csv(dest / "regenerated" / ".regenerated_abundance" / "dummy.csv")
    assert "N_s1" in abund.columns and "N_s2" in abund.columns


def test_merge_regenerated_alias_and_cli(tmp_path):
    a = _seed_run(tmp_path / "a", "left")
    b = _seed_run(tmp_path / "b", "right")
    dest = tmp_path / "cli_out"
    rc = merge_main(
        [
            "--mode",
            "regen",
            "--output_dir",
            str(dest),
            str(a),
            str(b),
        ]
    )
    assert rc == 0
    data = yaml.safe_load((dest / ".log" / "merge.yaml").read_text())
    assert data["mode"] == "regenerated"


def test_merge_needs_two_runs(tmp_path):
    a = _seed_run(tmp_path / "a", "s1")
    with pytest.raises(MergeError, match="two or more"):
        merge_runs([a], tmp_path / "out", "initial")


def _exec_env():
    env = os.environ.copy()
    env["PATH"] = str(REPO / "bin") + os.pathsep + env.get("PATH", "")
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return env


def test_merge_initial_exec_rebuilds_abundance(tmp_path):
    a = _seed_run(tmp_path / "a", "s1")
    b = _seed_run(tmp_path / "b", "s2")
    dest = tmp_path / "merged_initial"
    merge_runs([a, b], dest, "initial")
    proc = subprocess.run(
        [
            "bash",
            str(REPO / "bin" / "samovar"),
            "exec",
            "--output_dir",
            str(dest),
            "--startpoint",
            "abundance_tables",
            "--endpoint",
            "abundance_tables",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        env=_exec_env(),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    dummy = dest / "initial_abundance" / "dummy.csv"
    assert dummy.is_file()
    table = pd.read_csv(dummy)
    cols = list(table.columns)
    assert any("s1" in str(c) for c in cols)
    assert any("s2" in str(c) for c in cols)


def test_merge_regenerated_exec_starts_at_reprofile(tmp_path):
    a = _seed_run(tmp_path / "a", "s1")
    b = _seed_run(tmp_path / "b", "s2")
    dest = tmp_path / "merged_regen"
    merge_runs([a, b], dest, "regenerated")
    proc = subprocess.run(
        [
            "bash",
            str(REPO / "bin" / "samovar"),
            "exec",
            "--output_dir",
            str(dest),
            "--startpoint",
            "reprofile",
            "--endpoint",
            "reprofile",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        env=_exec_env(),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (dest / "reprofiled_annotations").is_dir()
    assert list((dest / "reprofiled_annotations").glob("*.csv")) or (
        dest / "reprofiled_annotations" / "trained_model.joblib"
    ).is_file()
