"""Checkpoints and optional .tmp cleanup for samovar exec."""

import os
import stat
import subprocess
from pathlib import Path

from samovar.exec_control import (
    CHECKPOINT_STEPS,
    cleanup_tmp,
    clear_checkpoints,
    listed_done,
    mark_done,
    should_skip,
)
from samovar.exec_control import main as exec_control_main


REPO = Path(__file__).resolve().parents[1]


def test_checkpoint_skip_mark_clear_and_redo(tmp_path, monkeypatch):
    assert "regenerate_reads" in CHECKPOINT_STEPS
    assert not should_skip(tmp_path, "annotate_initial")
    marker = mark_done(tmp_path, "annotate_initial")
    assert marker.is_file()
    assert should_skip(tmp_path, "annotate_initial")
    assert listed_done(tmp_path) == ["annotate_initial"]
    assert not should_skip(tmp_path, "annotate_initial", redo=True)
    monkeypatch.setenv("SAMOVAR_REDO", "1")
    assert not should_skip(tmp_path, "annotate_initial")
    monkeypatch.delenv("SAMOVAR_REDO")
    assert should_skip(tmp_path, "annotate_initial")
    clear_checkpoints(tmp_path, "annotate_initial")
    assert not should_skip(tmp_path, "annotate_initial")
    mark_done(tmp_path, "setup_reads")
    mark_done(tmp_path, "regenerate_reads")
    removed = clear_checkpoints(tmp_path)
    assert {p.name for p in removed} == {"setup_reads.done", "regenerate_reads.done"}
    assert listed_done(tmp_path) == []


def test_exec_control_cli_skip_mark_list(tmp_path):
    assert exec_control_main(["skip", str(tmp_path), "setup_reads"]) == 1
    assert exec_control_main(["mark", str(tmp_path), "setup_reads"]) == 0
    assert exec_control_main(["skip", str(tmp_path), "setup_reads"]) == 0
    assert exec_control_main(["list", str(tmp_path)]) == 0


def test_cleanup_tmp_removes_scratch_keeps_log(tmp_path):
    (tmp_path / ".log" / "checkpoints").mkdir(parents=True)
    (tmp_path / ".log" / "checkpoints" / "setup_reads.done").write_text("keep\n")
    (tmp_path / ".generate").mkdir()
    (tmp_path / ".generate" / "generate.sh").write_text("echo hi\n")
    (tmp_path / ".tmp" / "nested").mkdir(parents=True)
    (tmp_path / ".tmp" / "nested" / "x").write_text("x")
    (tmp_path / ".iss_full").mkdir()
    (tmp_path / ".iss_full" / "pool.fasta").write_text(">a\nA\n")
    (tmp_path / "initial_reports" / "tmp_sample.out").mkdir(parents=True)
    leftover = tmp_path / "regenerated" / "s.iss.tmp.foo"
    leftover.parent.mkdir()
    leftover.write_text("tmp")
    (tmp_path / "reprofiled_annotations").mkdir()
    (tmp_path / "reprofiled_annotations" / "keep.csv").write_text("ok\n")

    removed = cleanup_tmp(tmp_path)
    assert (tmp_path / ".log" / "checkpoints" / "setup_reads.done").is_file()
    assert (tmp_path / ".generate" / "generate.sh").is_file()
    assert (tmp_path / "reprofiled_annotations" / "keep.csv").is_file()
    assert (tmp_path / "initial_reports" / "tmp_sample.out").is_dir()
    assert not (tmp_path / ".tmp").exists()
    assert not (tmp_path / ".iss_full").exists()
    assert not leftover.exists()
    assert any(".tmp" in path for path in removed)


def test_mark_done_is_atomic(tmp_path):
    marker = mark_done(tmp_path, "setup_reads")
    assert marker.is_file()
    leftovers = list((tmp_path / ".log" / "checkpoints").glob(".setup_reads.done.*"))
    assert leftovers == []


def test_exec_redo_clears_checkpoints(tmp_path):
    log = tmp_path / ".log"
    ckpt = log / "checkpoints"
    ckpt.mkdir(parents=True)
    (ckpt / "setup_reads.done").write_text("old\n")
    script = log / "samovar.sh"
    script.write_text("#!/bin/bash\necho ran\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    result = subprocess.run(
        [
            "bash",
            str(REPO / "bin" / "samovar_exec"),
            "--output_dir",
            str(tmp_path),
            "--redo",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, "PYTHONPATH": str(REPO / "src")},
    )
    assert "Clearing pipeline checkpoints" in result.stdout
    assert not (ckpt / "setup_reads.done").exists()
    assert "Pipeline execution completed" in result.stdout


def test_exec_requires_prepare_after_generate(tmp_path):
    gen = tmp_path / ".generate"
    gen.mkdir()
    (gen / "generate.sh").write_text("#!/bin/bash\necho generated\n")
    (gen / "generate.sh").chmod(0o755)
    (tmp_path / "initial").mkdir()
    (tmp_path / "initial" / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    result = subprocess.run(
        [
            "bash",
            str(REPO / "bin" / "samovar_exec"),
            "--output_dir",
            str(tmp_path),
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(REPO / "src")},
    )
    assert result.returncode == 1
    assert "samovar.sh is missing" in result.stdout
    assert "samovar prepare" in result.stdout
    assert "Skipping generate" in result.stdout


def test_exec_skips_generate_when_reads_exist(tmp_path):
    gen = tmp_path / ".generate"
    gen.mkdir()
    marker = tmp_path / "generated.flag"
    (gen / "generate.sh").write_text(f"#!/bin/bash\ntouch '{marker}'\n")
    (gen / "generate.sh").chmod(0o755)
    (tmp_path / "initial").mkdir()
    (tmp_path / "initial" / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    log = tmp_path / ".log"
    log.mkdir()
    (log / "samovar.sh").write_text("#!/bin/bash\necho ran-main\n")
    (log / "samovar.sh").chmod(0o755)
    result = subprocess.run(
        [
            "bash",
            str(REPO / "bin" / "samovar_exec"),
            "--output_dir",
            str(tmp_path),
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, "PYTHONPATH": str(REPO / "src")},
    )
    assert "Skipping generate" in result.stdout
    assert "ran-main" in result.stdout
    assert not marker.exists()


def test_exec_cleanup_tmp_flag(tmp_path):
    log = tmp_path / ".log"
    log.mkdir()
    (log / "samovar.sh").write_text("#!/bin/bash\necho ran\n")
    scratch = tmp_path / ".tmp" / "junk"
    scratch.mkdir(parents=True)
    (scratch / "a").write_text("x")
    result = subprocess.run(
        [
            "bash",
            str(REPO / "bin" / "samovar_exec"),
            "--output_dir",
            str(tmp_path),
            "--cleanup-tmp",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, "PYTHONPATH": str(REPO / "src")},
    )
    assert "Cleaning temporary directories" in result.stdout
    assert not (tmp_path / ".tmp").exists()
    assert log.is_dir()
