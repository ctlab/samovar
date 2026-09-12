"""CLI and integration tests for ``samovar apply``.

Fixtures are built from production APIs and ``tests/tools`` only — not from
``samovar/samovar_apply`` examples.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

from samovar.apply import (
    ApplyError,
    IncompatibleInputError,
    MissingModelError,
    MissingPipelineError,
    apply_pipeline,
    configured_steps,
    load_pipeline_state,
    main as apply_main,
)
from samovar.paths import write_config
from samovar.reprofilers import run_reprofiler
from samovar.tools_import import import_tool

REPO = Path(__file__).resolve().parents[1]
TOOLS = Path(__file__).resolve().parent / "tools"
ADAPTING = TOOLS / "adapting_reprofiler.py"
COUNT = TOOLS / "count_annotations.py"


def _fastq_pair(dest: Path, n: int = 8) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    r1 = []
    r2 = []
    for i in range(n):
        taxid = 562 if i % 2 == 0 else 9606
        rid = f"r{i}|taxid:{taxid}"
        rec = f"@{rid}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
        r1.append(rec)
        r2.append(rec)
    (dest / "sample_R1.fastq").write_text("".join(r1), encoding="utf-8")
    (dest / "sample_R2.fastq").write_text("".join(r2), encoding="utf-8")


def _training_tables():
    n = 16
    tax = [9606, 9606, 562, 562] * (n // 4)
    regenerated = pd.DataFrame(
        {
            "seq": [f"t{i}" for i in range(n)],
            "taxid_dummy": tax,
            "length": [12] * n,
            "true": tax,
        }
    )
    initial = {
        "sample.annotation": pd.DataFrame(
            {
                "seq": ["a", "b", "c", "d"],
                "taxid_dummy": [9606, 562, 9606, 562],
                "length": [12, 12, 12, 12],
            }
        )
    }
    ground = {"dummy": pd.DataFrame({"taxid": [9606, 562], "N_1": [8, 8]})}
    return regenerated, ground, initial


def _dump(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def _source_pipeline(tmp_path: Path, monkeypatch, *, with_model: bool = True) -> Path:
    cfg = tmp_path / "install.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(REPO), "tools": {}}, also_repo_build=False)
    import_tool(
        name="adaptwrap",
        tool_type="ml",
        exec_path=str(ADAPTING),
        also_repo_build=False,
    )
    import_tool(
        name="apply_counts",
        tool_type="scoring",
        exec_path=str(COUNT),
        inputs="*annotations",
        also_repo_build=False,
    )
    source = tmp_path / "source_run"
    configs = source / ".log" / "configs"
    configs.mkdir(parents=True)
    _dump(
        configs / "config_init.yaml",
        {
            "r1_dir": str(source / "initial_trimmed"),
            "r2_dir": str(source / "initial_trimmed"),
            "output_dir": str(source / "initial_reports"),
            "run_config": [
                {
                    "run_name": "dummy",
                    "type": "dummy",
                    "cmd": "dummy",
                    "db_path": ".",
                    "threads": 1,
                }
            ],
        },
    )
    _dump(configs / "config_qc.yaml", {"output_dir": str(source), "qc": ""})
    _dump(
        configs / "config_scoring.yaml",
        {
            "output_dir": str(source),
            "scoring_tools": ["apply_counts"],
        },
    )
    _dump(
        configs / "config_export.yaml",
        {
            "output_dir": str(source),
            "export": "identity",
            "export_formats": ["abundance"],
        },
    )
    _dump(
        configs / "config_reprofiling.yaml",
        {
            "output_dir": str(source / "reprofiled_annotations"),
            "initial_dir": str(source / "initial_annotations"),
            "regenerated_path": str(
                source / "regenerated_annotations" / "combined_annotation_table.csv"
            ),
            "ground_truth_dir": str(source / "regenerated" / ".regenerated_abundance"),
            "reprofiler": "adaptwrap",
            "seed": 0,
        },
    )
    _dump(
        configs / "config_annotation2iss.yaml",
        {
            "annotation_dir": str(source / "initial_annotations"),
            "observed_abundance_dir": str(source / "initial_abundance"),
            "abundance_dir": str(source / "regenerated" / ".regenerated_abundance"),
            "output_dir": str(source / "regenerated"),
            "regeneration_mode": "direct",
            "table_reads_generator": "direct",
            "seed": 0,
            "cores": 1,
        },
    )
    (source / ".log" / "window.env").write_text(
        "export SAMOVAR_START=setup_reads\nexport SAMOVAR_END=viz_reprofiled\n",
        encoding="utf-8",
    )
    regenerated, ground, initial = _training_tables()
    regen_dir = source / "regenerated_annotations"
    regen_dir.mkdir(parents=True)
    regenerated.to_csv(regen_dir / "combined_annotation_table.csv", index=False)
    gt = source / "regenerated" / ".regenerated_abundance"
    gt.mkdir(parents=True)
    ground["dummy"].to_csv(gt / "dummy.csv", index=False)
    if with_model:
        run_reprofiler(
            "linear",
            regenerated=regenerated,
            ground_truth=ground,
            initial=initial,
            output_dir=source / "reprofiled_annotations",
            config={"seed": 0},
        )
    return source


def test_cli_help_lists_apply():
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
    assert "apply" in result.stdout
    assert "--full" in result.stdout


def test_apply_module_help():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [os.environ.get("PYTHON_PATH") or "python3", "-m", "samovar.apply", "--help"],
        cwd=REPO,
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode == 0
    assert "--full" in proc.stdout
    assert "--pipeline" in proc.stdout


def test_apply_missing_pipeline(tmp_path):
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    with pytest.raises(MissingPipelineError):
        apply_pipeline(reads, tmp_path / "missing", tmp_path / "out")
    rc = apply_main(
        [
            "--input_dir",
            str(reads),
            "--pipeline",
            str(tmp_path / "missing"),
            "--output_dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 1


def test_apply_refuses_inplace_overwrite(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    with pytest.raises(ApplyError, match="must not be the source"):
        apply_pipeline(reads, source, source, full=False)


def test_apply_incompatible_input(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(IncompatibleInputError):
        apply_pipeline(empty, source, tmp_path / "out")


def test_apply_missing_model(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch, with_model=False)
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    with pytest.raises(MissingModelError):
        apply_pipeline(reads, source, tmp_path / "out", full=False)


def test_apply_full_missing_regenerated(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    import shutil

    shutil.rmtree(source / "regenerated_annotations")
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    with pytest.raises(ApplyError):
        apply_pipeline(reads, source, tmp_path / "out", full=True)


def test_configured_steps_normal_skips_table_regen(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    state = load_pipeline_state(source)
    normal = configured_steps(state, full=False)
    full = configured_steps(state, full=True)
    assert "annotate_initial" in normal
    assert "reprofile" in normal
    assert "viz_reprofiled" in normal
    assert "regenerate_tables" not in normal
    assert "regenerate_tables" in full
    assert "regenerate_reads" not in full
    assert "regenerate_reads" not in normal


def test_normal_apply_integration_and_contracts(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    out = tmp_path / "apply_out"
    prov = apply_pipeline(reads, source, out, full=False, cores=1)
    assert prov["retrained"] is False
    assert prov["mode"] == "normal"
    profiles = list((out / "reprofiled_annotations").glob("*_reprofiled.csv"))
    assert profiles
    table = pd.read_csv(profiles[0])
    assert "taxid_SAMOVAR" in table.columns
    assert (out / "initial_annotations").is_dir()
    assert any((out / "initial_annotations").glob("*.csv"))
    assert (out / "reprofiled_annotations" / "trained_model.joblib").is_file()
    assert not (out / "reprofiled_annotations" / "adapted.txt").exists()
    counts = out / "apply_counts_scores" / "stage.txt"
    assert counts.is_file()
    stages = (out / "apply_counts_scores" / "stage.txt").read_text(encoding="utf-8")
    assert "viz_reprofiled" in stages or "viz_initial" in stages
    exports = out / "exports"
    assert exports.is_dir()
    apply_yaml = yaml.safe_load((out / ".log" / "apply.yaml").read_text())
    assert apply_yaml["pipeline"] == str(source.resolve())
    assert apply_yaml["retrained"] is False


@pytest.mark.optional
def test_apply_repeatability(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    a = tmp_path / "apply_a"
    b = tmp_path / "apply_b"
    apply_pipeline(reads, source, a, full=False, cores=1)
    apply_pipeline(reads, source, b, full=False, cores=1)
    pa = sorted((a / "reprofiled_annotations").glob("*_reprofiled.csv"))
    pb = sorted((b / "reprofiled_annotations").glob("*_reprofiled.csv"))
    assert pa and pb
    da = pd.read_csv(pa[0])
    db = pd.read_csv(pb[0])
    assert list(da["taxid_SAMOVAR"]) == list(db["taxid_SAMOVAR"])


@pytest.mark.optional
def test_full_apply_retrains_and_regenerates(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    out = tmp_path / "apply_full"
    prov = apply_pipeline(reads, source, out, full=True, cores=1)
    assert prov["retrained"] is True
    assert prov["mode"] == "full"
    assert (out / "reprofiled_annotations" / "adapted.txt").read_text().strip() == "retrained"
    assert (out / "regenerated" / ".regenerated_abundance").is_dir()
    assert any((out / "regenerated" / ".regenerated_abundance").glob("*.csv"))
    assert (out / "initial_abundance").is_dir()
    profiles = list((out / "reprofiled_annotations").glob("*_reprofiled.csv"))
    assert profiles
    assert "taxid_SAMOVAR" in pd.read_csv(profiles[0]).columns
    counts = out / "apply_counts_scores"
    assert counts.is_dir()


def test_normal_vs_full_semantics(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    normal = tmp_path / "normal"
    full = tmp_path / "full"
    apply_pipeline(reads, source, normal, full=False, cores=1)
    apply_pipeline(reads, source, full, full=True, cores=1)
    assert not (normal / "reprofiled_annotations" / "adapted.txt").exists()
    assert (full / "reprofiled_annotations" / "adapted.txt").is_file()
    assert not (normal / "regenerated" / ".regenerated_abundance").exists() or not any(
        (normal / "regenerated" / ".regenerated_abundance").glob("*.csv")
    )
    assert any((full / "regenerated" / ".regenerated_abundance").glob("*.csv"))
    for dest in (normal, full):
        table = pd.read_csv(
            next((dest / "reprofiled_annotations").glob("*_reprofiled.csv"))
        )
        assert "taxid_SAMOVAR" in table.columns
        assert (dest / "apply_counts_scores").is_dir()


def test_cli_apply_success(tmp_path, monkeypatch):
    source = _source_pipeline(tmp_path, monkeypatch)
    reads = tmp_path / "reads"
    _fastq_pair(reads)
    out = tmp_path / "cli_out"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = str(REPO / "bin") + os.pathsep + env.get("PATH", "")
    env["SAMOVAR_CONFIG"] = os.environ["SAMOVAR_CONFIG"]
    proc = subprocess.run(
        [
            "bash",
            str(REPO / "bin" / "samovar"),
            "apply",
            "--input_dir",
            str(reads),
            "--pipeline",
            str(source),
            "--output_dir",
            str(out),
            "--cores",
            "1",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert list((out / "reprofiled_annotations").glob("*_reprofiled.csv"))
