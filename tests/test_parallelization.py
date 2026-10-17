"""Scheduler caps, one Kaiju index load, and per-stage assembly threads."""

import threading
import time
from argparse import Namespace
from pathlib import Path

from samovar.annotators_wrapper import kaiju_multi_script, use_kaiju_multi
from samovar.assembly_profiling import assembly_stage_threads, run_megahit, run_prodigal
from samovar.config import AnnotatorConfig, PipelineConfig, _apply_translated_cli_flags
from samovar.memory import annotator_mem_mb, job_mem_mb
from samovar.tool_spec import annotator_priority, scheduler_threads


def test_heavy_databases_take_the_whole_memory_ceiling():
    limit = 32_000
    cores = 32
    assert annotator_mem_mb("kraken2", 16, cores, limit) == limit
    assert annotator_mem_mb("krakenuniq", 16, cores, limit) == limit
    assert annotator_mem_mb("metax", 4, cores, limit) == limit
    kaiju = annotator_mem_mb("kaiju", 8, cores, limit)
    assert kaiju == job_mem_mb(8, cores, limit)
    assert kaiju < limit
    assert annotator_priority("kaiju") > annotator_priority("kraken2")
    assert annotator_priority("kraken2") > annotator_priority("metax")
    assert annotator_priority("metax") > annotator_priority("gc")


def test_scheduler_reserves_only_what_the_program_uses():
    assert scheduler_threads("gc", 200) == 1
    assert scheduler_threads("kmer_encoder", 200) == 8
    assert scheduler_threads("kmer2", 200) == 8
    assert scheduler_threads("metax", 200) == 4
    assert scheduler_threads("kraken2", 200) == 16
    assert scheduler_threads("kaiju", 200) == 32
    assert scheduler_threads("assembly", 200) == 16
    assert scheduler_threads("custom_tool", 200) == 8
    assert scheduler_threads("gc", 1) == 1
    assert scheduler_threads("kraken2", 4) == 4


def test_assembly_stages_do_not_share_the_machine():
    assert assembly_stage_threads("megahit", 200) == 16
    assert assembly_stage_threads("flye", 200) == 16
    assert assembly_stage_threads("prodigal", 200) == 8
    assert assembly_stage_threads("anvio", 200) == 8
    assert assembly_stage_threads("gtdbtk", 200) == 8
    assert assembly_stage_threads("checkm2", 200) == 4
    assert assembly_stage_threads("identity", 3) == 3
    assert assembly_stage_threads("megahit", 2) == 2


def test_kaiju_multi_expands_one_fmi_glob(tmp_path):
    fmi = tmp_path / "kaiju.fmi"
    fmi.write_bytes(b"x")
    script = kaiju_multi_script(
        binary="kaiju-multi",
        nodes=str(tmp_path / "nodes.dmp"),
        database=str(tmp_path / "*.fmi"),
        r1=["/r/a_R1.fastq"],
        r2=[],
        outputs=["/o/a.out"],
        threads=2,
    )
    assert str(fmi) in script
    assert "*.fmi" not in script


def test_kaiju_multi_starts_at_two_samples():
    assert use_kaiju_multi(1) is False
    assert use_kaiju_multi(2) is True


def test_kaiju_multi_script_is_one_process():
    script = kaiju_multi_script(
        binary="kaiju-multi",
        nodes="/db/nodes.dmp",
        database="/db/kaiju.fmi",
        r1=["/r/s1_R1.fastq", "/r/s2_R1.fastq"],
        r2=["/r/s1_R2.fastq", "/r/s2_R2.fastq"],
        outputs=["/o/s1.out", "/o/s2.out"],
        threads=32,
        empty=["/o/empty.out"],
    )
    assert "kaiju-multi" in script
    assert script.count("kaiju-multi") == 1
    assert "-z 32" in script
    assert "/r/s1_R1.fastq,/r/s2_R1.fastq" in script
    assert "touch /o/empty.out" in script


def test_metax_extra_does_not_repeat_threads():
    cfg = PipelineConfig(
        cores=200,
        annotators=[
            AnnotatorConfig(
                run_name="metax",
                type="metax",
                db_path="/db",
                cmd="custom.sh",
                extra="",
            )
        ],
    )
    _apply_translated_cli_flags(cfg, Namespace(custom_cli_flags={}))
    extra = cfg.annotators[0].extra or ""
    assert "-t" not in extra
    assert "200" not in extra


def test_megahit_command_uses_the_stage_cap(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        "samovar.assembly_profiling.which_tool", lambda *args, **kwargs: "megahit"
    )
    monkeypatch.setattr(
        "samovar.assembly_profiling._run", lambda cmd, **kwargs: calls.append(list(cmd))
    )
    r1 = tmp_path / "sample_R1.fastq"
    r1.write_text("@r0\nACGT\n+\nIIII\n", encoding="utf-8")
    dest = tmp_path / "contigs.fa"
    run_megahit(str(r1), "", str(dest), threads=200)
    cmd = calls[0]
    assert cmd[cmd.index("-t") + 1] == "16"


def test_prodigal_runs_genomes_together(tmp_path, monkeypatch):
    depth = {"now": 0, "max": 0}
    lock = threading.Lock()

    def fake_run(cmd, **kwargs):
        del cmd, kwargs
        with lock:
            depth["now"] += 1
            depth["max"] = max(depth["max"], depth["now"])
        time.sleep(0.05)
        with lock:
            depth["now"] -= 1

    monkeypatch.setattr("samovar.assembly_profiling.which_tool", lambda *a, **k: "prodigal")
    monkeypatch.setattr("samovar.assembly_profiling._run", fake_run)
    src = tmp_path / "mags"
    src.mkdir()
    for name in ("a", "b", "c"):
        (src / f"{name}.fa").write_text(f">{name}\nACGT\n", encoding="utf-8")
    run_prodigal(str(src), str(tmp_path / "genes"), threads=4)
    assert depth["max"] >= 2
    assert not list(Path(src).glob("*.faa"))
