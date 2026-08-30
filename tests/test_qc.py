"""Read QC contract: identity default, GC filter, stage/postfix selection."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import yaml

from samovar.config import PipelineConfig, setup_pipeline
from samovar.paths import write_config
from samovar.qc import (
    canonical_qc_name,
    count_fastq_records,
    default_postfix_for_qc,
    gc_fraction,
    parse_qc_postfix,
    qc_name_for_sample,
    require_known_qc,
    trim_directory,
    trim_stage,
)
from samovar.seqio import iter_fastq_records
from samovar.tools_import import import_tool, main as import_main


REPO = Path(__file__).resolve().parents[1]
GC_FILTER = REPO / "tests" / "tools" / "gc_filter.py"
TOY_R1 = REPO / "tests" / "data" / "reads" / "1_full_R1.fastq"
TOY_R2 = REPO / "tests" / "data" / "reads" / "1_full_R2.fastq"


def _fastq(path: Path, records: list) -> None:
    lines = []
    for i, seq in enumerate(records):
        lines.extend([f"@r{i}", seq, "+", "I" * len(seq)])
    path.write_text("\n".join(lines) + "\n")


def test_gc_filter_drops_out_of_range(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="gc_filter",
        tool_type="QC",
        exec_path=str(GC_FILTER),
        also_repo_build=False,
    )
    src = tmp_path / "src"
    dest = tmp_path / "dest"
    src.mkdir()
    _fastq(src / "1_full_R1.fastq", ["AAAA", "GGGG", "ACGT"])
    _fastq(src / "1_full_R2.fastq", ["AAAA", "CCCC", "TGCA"])
    trim_directory(
        src,
        dest,
        stage_qc="gc_filter",
        config={"min_gc": 0.4, "max_gc": 0.6},
    )
    kept = [rec[1].strip() for rec in iter_fastq_records(dest / "1_full_R1.fastq")]
    assert kept == ["ACGT"]
    assert count_fastq_records(dest / "1_full_R2.fastq") == 1


def test_identity_when_qc_absent(tmp_path):
    src = tmp_path / "src"
    dest = tmp_path / "dest"
    src.mkdir()
    _fastq(src / "s_R1.fastq", ["ACGT"])
    _fastq(src / "s_R2.fastq", ["TGCA"])
    trim_directory(src, dest, stage_qc="")
    assert (dest / "s_R1.fastq").is_file()
    assert count_fastq_records(dest / "s_R1.fastq") == 1


def test_different_qc_initial_vs_generated(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="gc_filter",
        tool_type="qc",
        exec_path=str(GC_FILTER),
        also_repo_build=False,
    )
    run = tmp_path / "run"
    (run / "initial").mkdir(parents=True)
    (run / "regenerated").mkdir()
    _fastq(run / "initial" / "1_full_R1.fastq", ["GGGG", "ACGT"])
    _fastq(run / "initial" / "1_full_R2.fastq", ["CCCC", "TGCA"])
    _fastq(run / "regenerated" / "1_k_R1.fastq", ["GGGG", "ACGT"])
    _fastq(run / "regenerated" / "1_k_R2.fastq", ["CCCC", "TGCA"])
    qc_cfg = {
        "qc_initial": "gc_filter",
        "qc_generated": "",
        "min_gc": 0.4,
        "max_gc": 0.6,
    }
    trim_stage(run, "initial", qc_cfg)
    trim_stage(run, "generated", qc_cfg)
    assert count_fastq_records(run / "initial_trimmed" / "1_full_R1.fastq") == 1
    assert count_fastq_records(run / "regenerated_trimmed" / "1_k_R1.fastq") == 2


def test_hybrid_postfix_qc(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="gc_filter",
        tool_type="QC",
        exec_path=str(GC_FILTER),
        also_repo_build=False,
    )
    src = tmp_path / "src"
    dest = tmp_path / "dest"
    src.mkdir()
    for tech in ("illumina", "ont"):
        _fastq(src / f"1_{tech}_R1.fastq", ["GGGG", "ACGT"])
        _fastq(src / f"1_{tech}_R2.fastq", ["CCCC", "TGCA"])
    _fastq(src / "1_full_R1.fastq", ["GGGG", "ACGT"])
    _fastq(src / "1_full_R2.fastq", ["CCCC", "TGCA"])
    mapping = parse_qc_postfix(["illumina:gc_filter", "ont:"])
    assert mapping["illumina"] == "gc_filter"
    assert mapping["ont"] == ""
    assert qc_name_for_sample("1_illumina", stage_qc="", postfix_map=mapping) == "gc_filter"
    trim_directory(
        src,
        dest,
        stage_qc="",
        postfix_map=mapping,
        config={"min_gc": 0.4, "max_gc": 0.6},
    )
    assert count_fastq_records(dest / "1_illumina_R1.fastq") == 1
    assert count_fastq_records(dest / "1_ont_R1.fastq") == 2
    assert count_fastq_records(dest / "1_full_R1.fastq") == 2


def test_prepare_same_and_split_qc(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="gc_filter",
        tool_type="QC",
        exec_path=str(GC_FILTER),
        also_repo_build=False,
    )
    reads = tmp_path / "reads"
    reads.mkdir()
    _fastq(reads / "1_full_R1.fastq", ["ACGT"])
    _fastq(reads / "1_full_R2.fastq", ["ACGT"])
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(reads),
        output_dir=str(tmp_path / "out_same"),
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=None,
        qc="gc_filter",
        qc_initial=None,
        qc_generated=None,
        qc_postfix=None,
        tool_flags=[["qc", "--min-gc 0.3 --max-gc 0.7"]],
    )
    same = PipelineConfig.from_args(args)
    assert same.qc_initial == "gc_filter"
    assert same.qc_generated == "gc_filter"
    assert "--min-gc" in (same.qc_flags or "")
    text = Path(setup_pipeline(args)["pipeline"]).read_text()
    assert "samovar.qc trim" in text
    assert "qc_initial" in text
    yaml_qc = yaml.safe_load(
        (tmp_path / "out_same" / ".log" / "configs" / "config_qc.yaml").read_text()
    )
    assert yaml_qc["qc_initial"] == "gc_filter"
    init_yaml = yaml.safe_load(
        (tmp_path / "out_same" / ".log" / "configs" / "config_init.yaml").read_text()
    )
    assert init_yaml["r1_dir"].endswith("initial_trimmed")

    args2 = argparse.Namespace(
        input_config=None,
        input_dir=str(reads),
        output_dir=str(tmp_path / "out_split"),
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=None,
        qc=None,
        qc_initial="gc_filter",
        qc_generated="none",
        qc_postfix=["illumina:gc_filter"],
        tool_flags=None,
    )
    split = PipelineConfig.from_args(args2)
    assert split.qc_initial == "gc_filter"
    assert split.qc_generated == ""
    assert split.qc_postfix["illumina"] == "gc_filter"


def test_fastp_aliases_and_postfix(tmp_path, monkeypatch):
    """Illumina/BGI names resolve to imported ``fastp``; postfix map is filled."""
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="fastp",
        tool_type="QC",
        exec_path=str(GC_FILTER),
        lazy_install="conda install -y bioconda::fastp",
        also_repo_build=False,
    )
    assert canonical_qc_name("illumina") == "fastp"
    assert canonical_qc_name("bgi") == "fastp"
    assert canonical_qc_name("mgi") == "fastp"
    assert require_known_qc("illumina") == "fastp"
    assert require_known_qc("bgi") == "fastp"
    assert default_postfix_for_qc("fastp")["illumina"] == "fastp"
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(tmp_path),
        output_dir=str(tmp_path / "out"),
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=None,
        qc="illumina",
        qc_initial=None,
        qc_generated=None,
        qc_postfix=None,
        tool_flags=None,
    )
    pipe = PipelineConfig.from_args(args)
    assert pipe.qc_initial == "fastp"
    assert pipe.qc_generated == "fastp"
    assert pipe.qc_postfix["illumina"] == "fastp"
    assert pipe.qc_postfix["bgi"] == "fastp"
    assert pipe.qc_postfix["mgi"] == "fastp"


def test_import_qc_pytest(tmp_path, monkeypatch):
    from samovar.paths import update_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    rc = import_main(
        [
            "-n",
            "gc_filter",
            "--type",
            "QC",
            "--exec-path",
            str(GC_FILTER),
            "--pytest",
        ]
    )
    assert rc == 0
    tools = json.loads(cfg.read_text())["tools"]
    assert any(str(k).startswith("gc_filter") for k in tools)


def test_toy_reads_gc_filter(tmp_path, monkeypatch):
    """Bundled toy FASTQ: kept reads stay inside the GC window."""
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="gc_filter",
        tool_type="QC",
        exec_path=str(GC_FILTER),
        also_repo_build=False,
    )
    src = tmp_path / "toy"
    dest = tmp_path / "trimmed"
    src.mkdir()
    shutil.copy2(TOY_R1, src / TOY_R1.name)
    if TOY_R2.is_file():
        shutil.copy2(TOY_R2, src / TOY_R2.name)
    before = count_fastq_records(src / TOY_R1.name)
    assert before > 10
    trim_directory(
        src,
        dest,
        stage_qc="gc_filter",
        config={"min_gc": 0.35, "max_gc": 0.65},
    )
    after = count_fastq_records(dest / TOY_R1.name)
    assert after <= before
    for _h, seq, _p, _q in iter_fastq_records(dest / TOY_R1.name):
        assert 0.35 <= gc_fraction(seq) <= 0.65
    assert require_known_qc("gc_filter") == "gc_filter"
