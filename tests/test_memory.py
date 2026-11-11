"""Memory budget, streaming feature rows, and chunked reprofile prediction."""

import os

import numpy as np
import pandas as pd

from samovar.gc import main as gc_main
from samovar.memory import (
    job_mem_mb,
    parallel_jobs,
    parse_mem,
    predict_n_jobs,
    proba_chunk_rows,
    read_csv_memmap,
)
from samovar.qc import gc_fraction
from samovar.reprofiling import predict_taxid, train_models, preprocess_data

import pytest


def test_allocation_caps_memory_and_cores(monkeypatch):
    from samovar.memory import available_cores, memory_limit_mb

    monkeypatch.setattr("samovar.memory.available_mem_mb", lambda: 16_000)
    monkeypatch.setattr("samovar.memory._cgroup_text", lambda name: None)
    monkeypatch.delenv("SAMOVAR_MEM_MB", raising=False)
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "8192")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    assert memory_limit_mb() == 8192
    assert memory_limit_mb(1000) == 1000
    assert memory_limit_mb(20_000) == 8192
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.setenv("SLURM_MEM_PER_CPU", "1024")
    assert memory_limit_mb() == 2048
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    monkeypatch.setattr("samovar.memory._cgroup_text", lambda name: "4294967296" if name == "memory.max" else None)
    assert memory_limit_mb() == 4096
    monkeypatch.setattr("samovar.memory._cgroup_text", lambda name: "max" if name == "memory.max" else None)
    assert memory_limit_mb() == int(16_000 * 0.85)

    monkeypatch.setattr("samovar.memory.os.cpu_count", lambda: 8)
    monkeypatch.setattr("samovar.memory.os.sched_getaffinity", lambda _pid: set(range(8)))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "4")
    assert available_cores() == 2
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.delenv("SLURM_CPUS_ON_NODE", raising=False)

    def cgroup(name):
        if name == "cpu.max":
            return "200000 100000"
        if name == "cpuset.cpus":
            return "0-1,4"
        return None

    monkeypatch.setattr("samovar.memory._cgroup_text", cgroup)
    assert available_cores() == 2


def test_parse_mem_matches_slurm():
    assert parse_mem("16000") == 16000
    assert parse_mem("16G") == 16 * 1024
    assert parse_mem("16GB") == 16 * 1024
    assert parse_mem("1.5G") == 1536
    assert parse_mem("2T") == 2 * 1024 * 1024
    assert parse_mem("512K") == 1
    assert parse_mem(32) == 32


def test_prepare_keeps_requested_cores_for_exec(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    (tmp_path / "reads").mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(tmp_path / "reads"),
            "output_dir": str(tmp_path / "out"),
            "threads": None,
            "cores": 4,
            "mem_mb": None,
            "custom_cli_flags": {},
        },
    )()
    config = PipelineConfig.from_args(args)
    assert config.cores == 4
    script = open(config.generate_pipeline(str(tmp_path / "out")), encoding="utf-8").read()
    assert 'export SAMOVAR_CORES="4"' in script
    assert 'if [ "$DETECTED_CORES" -lt "$SAMOVAR_CORES" ]; then' in script
    init = (tmp_path / "out" / ".log" / "configs" / "config_init.yaml").read_text(encoding="utf-8")
    assert "cores: 4" in init


def test_prepare_bakes_mem_into_the_pipeline(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    (tmp_path / "reads").mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(tmp_path / "reads"),
            "output_dir": str(tmp_path / "out"),
            "threads": None,
            "cores": 4,
            "mem_mb": parse_mem("4G"),
            "custom_cli_flags": {},
        },
    )()
    config = PipelineConfig.from_args(args)
    assert config.mem_mb == 4096
    script = open(config.generate_pipeline(str(tmp_path / "out")), encoding="utf-8").read()
    assert 'PREPARED_MEM="4096"' in script
    assert "export SAMOVAR_MEM_MB" in script


def test_job_mem_shares_sum_to_the_ceiling():
    limit = 10_000
    cores = 10
    shares = [job_mem_mb(1, cores, limit) for _ in range(cores)]
    assert sum(shares) == limit
    assert job_mem_mb(4, cores, limit) == 4_000


def test_parallel_jobs_shrink_when_the_matrix_does_not_fit():
    assert parallel_jobs(16, nbytes=1000, mem_mb=64_000) == 16
    huge = 8 * 1024 * 1024 * 1024
    assert parallel_jobs(16, nbytes=huge, mem_mb=512) == 1


def test_predict_threads_are_capped_and_restored():
    assert predict_n_jobs(32) == 2
    assert predict_n_jobs(1) == 1
    rng = np.random.default_rng(0)
    frame = pd.DataFrame(
        {
            "taxid_a": rng.integers(1, 4, 30),
            "taxid_b": rng.integers(1, 4, 30),
            "length": 100,
            "true": rng.integers(0, 3, 30),
        }
    )
    model, _models, _metrics, cols = train_models(
        preprocess_data(frame), test_size=0.25, n_jobs=8, methods="random_forest"
    )
    assert model.n_jobs == 8
    seen = []
    real = model.predict_proba

    def _spy(block):
        seen.append(model.n_jobs)
        return real(block)

    model.predict_proba = _spy
    os.environ["SAMOVAR_PREDICT_CHUNK"] = "4"
    try:
        out = predict_taxid(frame, model_path=model, feature_cols=cols)
    finally:
        os.environ.pop("SAMOVAR_PREDICT_CHUNK", None)
    assert seen
    assert max(seen) <= 2
    assert model.n_jobs == 8
    assert "taxid_SAMOVAR" in out.columns
    assert proba_chunk_rows(3, 3, mem_mb=64) >= 1


def test_gc_streams_paired_reads(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    r1.write_text("@r0\nAAAA\n+\nIIII\n")
    r2.write_text("@r0\nGGGG\n+\nIIII\n")
    out = tmp_path / "gc.out"
    assert gc_main(["-i", str(r1), "-I", str(r2), "-o", str(out)]) == 0
    text = out.read_text().strip().splitlines()
    assert text[0] == "seq\tGC"
    assert text[1].startswith("r0\t")
    assert abs(float(text[1].split("\t")[1]) - gc_fraction("AAAAGGGG")) < 1e-6


@pytest.mark.optional
def test_numeric_csv_uses_memmap(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_MEMMAP_MIN_BYTES", "0")
    src = tmp_path / "ann.csv"
    pd.DataFrame({"seq": ["a", "b"], "length": [10, 12], "taxid_kaiju": [562, 9606]}).to_csv(
        src, index=False
    )
    frame = read_csv_memmap(src, tmp_path / "store")
    assert list(frame["length"]) == [10, 12]
    assert list(frame["seq"]) == ["a", "b"]
    assert any((tmp_path / "store").glob("*.f64"))


def test_reprofile_writes_one_csv_per_sample(tmp_path):
    from samovar.reprofilers import run_reprofiler

    regenerated = pd.DataFrame(
        {
            "taxid_kaiju": [562, 562, 9606, 9606, 562, 9606, 562, 9606],
            "length": [100] * 8,
            "true": [562, 562, 9606, 9606, 562, 9606, 562, 9606],
        }
    )
    trained = tmp_path / "trained"
    run_reprofiler(
        "random_forest",
        regenerated=regenerated,
        initial={
            "s0": pd.DataFrame({"taxid_kaiju": [562], "length": [100]}),
        },
        output_dir=trained,
        config={"seed": 0},
    )
    model = trained / "trained_model.joblib"
    initial = tmp_path / "initial"
    initial.mkdir()
    for name, taxid in (("s0.annotation", 562), ("s1.annotation", 9606)):
        pd.DataFrame({"taxid_kaiju": [taxid], "length": [90]}).to_csv(
            initial / f"{name}.csv", index=False
        )
    dest = tmp_path / "applied"
    from samovar.reprofilers import apply_saved_reprofiler

    apply_saved_reprofiler(model_path=model, initial_dir=initial, output_dir=dest)
    written = sorted(p.name for p in dest.glob("*_reprofiled.csv"))
    assert written == ["s0_reprofiled.csv", "s1_reprofiled.csv"]
    for path in dest.glob("*_reprofiled.csv"):
        table = pd.read_csv(path)
        assert list(table["taxid_SAMOVAR"])
        assert len(table) == 1
