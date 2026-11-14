"""Repeat package manifest: layout, seeds, and refusal before overwrite."""

import argparse
import json
import os
from pathlib import Path

import pytest


def test_repeat_seed_is_stable_and_distinct():
    from samovar.repeats import repeat_seed

    assert repeat_seed(42, 1, "repeat") == repeat_seed(42, 1, "repeat")
    assert repeat_seed(42, 1, "repeat") != repeat_seed(42, 2, "repeat")
    assert repeat_seed(42, 1, "repeat") != repeat_seed(42, 1, "sample:A")
    assert repeat_seed(42, 1, "repeat") != repeat_seed(7, 1, "repeat")


def test_parse_repeat_rejects_bad_counts():
    from samovar.repeats import parse_repeat

    assert parse_repeat(None) is None
    assert parse_repeat(1) == 1
    assert parse_repeat("3") == 3
    for bad in (0, -2, 1.5, True, "1.5", "no"):
        with pytest.raises(ValueError):
            parse_repeat(bad)


def test_package_roundtrip_and_resume(tmp_path):
    from samovar.repeats import load_manifest, publish_package, repeat_seed

    root = tmp_path / "a batch"
    path = publish_package(
        root, count=3, base_seed=42, kind="reads", identity="reads:/data/in"
    )
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["schema"] == 1
    assert [row["out_dir"] for row in saved["repeats"]] == ["out_0001", "out_0002", "out_0003"]
    assert saved["repeats"][0]["seed"] == repeat_seed(42, 1, "repeat")
    assert (root / "out_0001").is_dir()
    before = path.read_bytes()
    old = os.getcwd()
    try:
        os.chdir("/")
        publish_package(root, count=3, base_seed=42, kind="reads", identity="reads:/data/in")
    finally:
        os.chdir(old)
    assert path.read_bytes() == before
    assert load_manifest(root)["repeats"][2]["repeat_id"] == "0003"


def test_atomic_replace_keeps_the_previous_manifest(tmp_path, monkeypatch):
    from samovar.repeats import _atomic_json

    path = tmp_path / ".log" / "repeats.json"
    _atomic_json(path, {"schema": 1, "repeat": 1})
    before = path.read_text(encoding="utf-8")

    def fail(src, dst):
        raise OSError("replace failed")

    monkeypatch.setattr("samovar.repeats.os.replace", fail)
    with pytest.raises(OSError):
        _atomic_json(path, {"schema": 1, "repeat": 9})
    assert path.read_text(encoding="utf-8") == before
    assert '"repeat": 1' in before


def test_conflict_does_not_add_children(tmp_path):
    from samovar.repeats import publish_package

    root = tmp_path / "pkg"
    publish_package(root, count=1, base_seed=42, kind="reads", identity="reads:/in")
    marker = root / "out_0001" / "keep.txt"
    marker.write_text("stay", encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        publish_package(root, count=3, base_seed=42, kind="reads", identity="reads:/in")
    assert marker.read_text(encoding="utf-8") == "stay"
    assert not (root / "out_0002").exists()


def test_refuse_nonempty_and_nested(tmp_path):
    from samovar.repeats import publish_package

    busy = tmp_path / "busy"
    busy.mkdir()
    (busy / "samovar.sh").write_text("x", encoding="utf-8")
    with pytest.raises(ValueError, match="not empty"):
        publish_package(busy, count=1, base_seed=42, kind="reads", identity="reads:/in")
    assert not (busy / ".log" / "repeats.json").exists()

    root = tmp_path / "pkg"
    publish_package(root, count=1, base_seed=42, kind="generated", identity="g")
    with pytest.raises(ValueError, match="nest"):
        publish_package(root / "out_0001", count=2, base_seed=42, kind="reads", identity="reads:/in")
    from samovar.config import PipelineConfig
    from samovar.repeats import prepare_package

    with pytest.raises(ValueError, match="nest"):
        prepare_package(PipelineConfig(output_dir=str(root / "out_0001"), input_dir=str(tmp_path / "reads")))
    assert list((root / "out_0001").iterdir()) == []


def test_prepare_resumes_generated_package(tmp_path):
    from samovar.config import PipelineConfig
    from samovar.repeats import prepare_package, publish_package

    root = tmp_path / "generated"
    publish_package(root, count=3, base_seed=42, kind="generated", identity="genomes|iss")
    (root / "initial").mkdir()
    config = PipelineConfig(
        output_dir=str(root),
        input_dir=str(root / "initial"),
        repeat=3,
        regeneration_seed=42,
    )
    prepare_package(config)
    assert not (root / "out_0001" / "out_0001").exists()
    with pytest.raises(ValueError, match="does not match"):
        prepare_package(
            PipelineConfig(
                output_dir=str(root),
                input_dir=str(root / "initial"),
                repeat=2,
                regeneration_seed=42,
            )
        )
    assert not (root / "out_0004").exists()


def test_prepare_cli_overrides_yaml_and_skips_root_pipeline(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig, setup_pipeline
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    reads = tmp_path / "reads"
    reads.mkdir()
    spec = tmp_path / "in.yaml"
    spec.write_text(
        f"repeat: 3\nseed: 42\ninput_dir: {reads}\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        input_config=str(spec),
        input_dir=str(reads),
        output_dir=str(tmp_path / "out"),
        threads=None,
        cores=1,
        mem_mb=None,
        custom_cli_flags={},
        repeat=1,
        add_annotator=False,
    )
    result = setup_pipeline(args)
    manifest = json.loads(Path(result["repeats"]).read_text(encoding="utf-8"))
    assert manifest["repeat"] == 1
    assert (tmp_path / "out" / "out_0001").is_dir()
    assert not (tmp_path / "out" / "out_0003").exists()
    assert not (tmp_path / "out" / ".log" / "samovar.sh").exists()
    plain = tmp_path / "plain.yaml"
    plain.write_text(f"input_dir: {reads}\n", encoding="utf-8")
    single = PipelineConfig.from_args(
        argparse.Namespace(**{**vars(args), "repeat": None, "input_config": str(plain), "output_dir": str(tmp_path / "single")})
    )
    assert single.repeat is None


def test_generate_repeat_package_and_single_seed(tmp_path):
    import argparse

    from samovar.iss_config import setup_iss_test

    single = tmp_path / "single"
    setup_iss_test(
        argparse.Namespace(
            genome_dir=str(tmp_path / "genomes"),
            output_dir=str(single),
            host_genome="",
            n_samples=1,
            total_reads=100,
            host_fraction="0",
            seed=123,
            model="hiseq",
        )
    )
    assert not (single / ".log" / "repeats.json").exists()
    text = (single / ".generate" / "configs" / "iss_config.yaml").read_text(encoding="utf-8")
    assert "seed: 123" in text

    batch = tmp_path / "batch dir"
    setup_iss_test(
        argparse.Namespace(
            genome_dir=str(tmp_path / "genomes"),
            output_dir=str(batch),
            host_genome="",
            n_samples=1,
            total_reads=100,
            host_fraction="0",
            seed=123,
            model="hiseq",
            repeat=1,
        )
    )
    manifest = json.loads((batch / ".log" / "repeats.json").read_text(encoding="utf-8"))
    assert manifest["repeat"] == 1
    assert manifest["base_seed"] == 123
    assert (batch / "out_0001").is_dir()
    assert (batch / ".generate" / "generate.sh").is_file()
    with pytest.raises(SystemExit, match="does not match"):
        setup_iss_test(
            argparse.Namespace(
                genome_dir=str(tmp_path / "genomes"),
                output_dir=str(batch),
                host_genome="",
                n_samples=1,
                total_reads=100,
                host_fraction="0",
                seed=123,
                model="hiseq",
                repeat=3,
            )
        )
    assert not (batch / "out_0003").exists()


def _generate_args(tmp_path, output, **extra):
    import argparse

    fields = dict(
        genome_dir=str(tmp_path / "genomes"),
        output_dir=str(output),
        host_genome="",
        n_samples=2,
        total_reads=20,
        host_fraction="0",
        seed=42,
        model="hiseq",
        repeat=2,
    )
    fields.update(extra)
    return argparse.Namespace(**fields)


def test_fixed_abundance_is_shared_and_seeds_differ(tmp_path):
    import yaml

    from samovar.iss_config import setup_iss_test
    from samovar.repeats import repeat_seed

    table = tmp_path / "abundance.csv"
    original = "taxid,N_A,N_B\n562,10,0\n9606,0,4\n"
    table.write_text(original, encoding="utf-8")
    batch = tmp_path / "batch"
    setup_iss_test(_generate_args(tmp_path, batch, abundance_table=str(table)))
    assert table.read_text(encoding="utf-8") == original
    seeds = []
    for name in ("out_0001", "out_0002"):
        cfg = yaml.safe_load((batch / name / ".generate" / "configs" / "iss_config.yaml").read_text())
        assert Path(cfg["abundance_table"]).resolve() == table.resolve()
        assert cfg["sample_names"] == ["A", "B"]
        assert cfg["output_dir"].endswith(f"{name}/initial")
        seeds.append(cfg["seed"])
        script = (batch / name / ".generate" / "generate.sh").read_text()
        assert "workflow/iss_test/Snakefile" in script
    assert seeds == [repeat_seed(42, 1, "reads"), repeat_seed(42, 2, "reads")]
    root_script = (batch / ".generate" / "generate.sh").read_text()
    assert "out_0001" in root_script and "out_0002" in root_script
    assert "mark" in root_script


def test_community_repeats_use_the_community_seed(tmp_path):
    import yaml

    from samovar.iss_config import setup_iss_test
    from samovar.repeats import repeat_seed

    batch = tmp_path / "drawn"
    setup_iss_test(_generate_args(tmp_path, batch))
    cfg = yaml.safe_load((batch / "out_0002" / ".generate" / "configs" / "iss_config.yaml").read_text())
    assert cfg["seed"] == repeat_seed(42, 2, "community")
    assert "abundance_table" not in cfg
    again = tmp_path / "drawn-again"
    setup_iss_test(_generate_args(tmp_path, again))
    cfg2 = yaml.safe_load((again / "out_0002" / ".generate" / "configs" / "iss_config.yaml").read_text())
    assert cfg2["seed"] == cfg["seed"]


def test_generator_failure_is_recorded_without_fake_reads(tmp_path):
    import subprocess

    from samovar.repeats import publish_package, write_repeat_orchestrator

    root = tmp_path / "pkg"
    publish_package(root, count=2, base_seed=42, kind="generated", identity="g")
    rows = [
        {"repeat_id": "0001", "out_dir": "out_0001"},
        {"repeat_id": "0002", "out_dir": "out_0002"},
    ]
    for row, code in zip(rows, (0, 1)):
        script = root / row["out_dir"] / ".generate" / "generate.sh"
        script.parent.mkdir(parents=True, exist_ok=True)
        script.write_text(f"#!/bin/bash\nexit {code}\n", encoding="utf-8")
        script.chmod(0o755)
    write_repeat_orchestrator(root, rows)
    result = subprocess.run(["bash", str(root / ".generate" / "generate.sh")], check=False, text=True)
    manifest = json.loads((root / ".log" / "repeats.json").read_text(encoding="utf-8"))
    states = {row["repeat_id"]: row["state"] for row in manifest["repeats"]}
    assert result.returncode != 0
    assert states == {"0001": "generated", "0002": "error"}
    assert manifest["repeats"][1]["error"] == "generate failed"
    assert list((root / "out_0002").rglob("*.fastq")) == []


def _prepared_args(tmp_path, output, **extra):
    import argparse

    fields = dict(
        input_config=None,
        input_dir=None,
        output_dir=str(output),
        threads=None,
        cores=1,
        mem_mb=None,
        custom_cli_flags={},
        repeat=None,
        add_annotator=False,
        mode=None,
    )
    fields.update(extra)
    return argparse.Namespace(**fields)


def test_prepare_inherits_manifest_and_rejects_a_different_n(tmp_path, monkeypatch):
    from samovar.config import setup_pipeline
    from samovar.paths import write_config
    from samovar.repeats import publish_package

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    root = tmp_path / "pkg"
    publish_package(root, count=2, base_seed=42, kind="generated", identity="g")
    setup_pipeline(_prepared_args(tmp_path, root, mode="benchmark"))
    assert (root / "out_0001" / ".log" / "samovar.sh").is_file()
    assert (root / "out_0002" / ".log" / "samovar.sh").is_file()
    assert "viz_initial" in (root / "out_0001" / ".log" / "window.env").read_text()
    full = tmp_path / "full"
    publish_package(full, count=1, base_seed=42, kind="generated", identity="g")
    setup_pipeline(_prepared_args(tmp_path, full, mode="full"))
    assert "viz_reprofiled" in (full / "out_0001" / ".log" / "window.env").read_text()
    with pytest.raises(SystemExit):
        setup_pipeline(_prepared_args(tmp_path, root, repeat=9))
    assert not (root / "out_0009").exists()


def test_existing_fastq_cannot_be_repeated(tmp_path, monkeypatch):
    from samovar.config import setup_pipeline
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    reads = tmp_path / "reads"
    reads.mkdir()
    out = tmp_path / "out"
    with pytest.raises(SystemExit):
        setup_pipeline(_prepared_args(tmp_path, out, input_dir=str(reads), repeat=3))
    assert not (out / "out_0001").exists()


def test_exec_continues_after_one_failure_and_resume_skips_successes(tmp_path):
    from samovar.repeats import exec_package, publish_package

    root = tmp_path / "pkg"
    publish_package(root, count=3, base_seed=42, kind="generated", identity="g")
    before = (root / ".log" / "repeats.json").read_bytes()
    log = root / "ran.txt"
    for name, code in (("out_0001", 0), ("out_0002", 1), ("out_0003", 0)):
        script = root / name / ".log" / "samovar.sh"
        script.parent.mkdir(parents=True, exist_ok=True)
        script.write_text(
            f"#!/bin/bash\necho {name} >> {log}\nexit {code}\n",
            encoding="utf-8",
        )
        script.chmod(0o755)
        (root / name / "keep.txt").write_text(name, encoding="utf-8")
    assert exec_package(root) == 1
    summary = json.loads((root / "summary" / "status.json").read_text())
    assert summary["status"] == "incomplete"
    assert summary["n_expected"] == 3
    assert summary["n_success"] == 2
    assert summary["errors"][0]["repeat_id"] == "0002"
    assert "mean" not in summary
    assert log.read_text().split() == ["out_0001", "out_0002", "out_0003"]
    script = root / "out_0002" / ".log" / "samovar.sh"
    script.write_text(f"#!/bin/bash\necho out_0002 >> {log}\nexit 0\n", encoding="utf-8")
    assert exec_package(root) == 0
    assert log.read_text().split() == ["out_0001", "out_0002", "out_0003", "out_0002"]
    assert (root / ".log" / "repeats.json").read_bytes() != before
    again = json.loads((root / ".log" / "repeats.json").read_text())
    assert [row["seed"] for row in again["repeats"]] == [
        json.loads(before)["repeats"][i]["seed"] for i in range(3)
    ]
    assert (root / "out_0003" / "keep.txt").read_text() == "out_0003"
    summary = json.loads((root / "summary" / "status.json").read_text())
    assert summary["status"] == "complete"
    assert summary["n_success"] == 3


def test_status_with_no_success_is_no_data(tmp_path):
    from samovar.repeats import publish_package, set_repeat_state, write_status_summary

    root = tmp_path / "pkg"
    publish_package(root, count=1, base_seed=7, kind="generated", identity="g")
    set_repeat_state(root, "0001", "failed", "generate failed", "pipeline")
    payload = json.loads(write_status_summary(root).read_text())
    assert payload == {
        "errors": [{"error": "generate failed", "repeat_id": "0001"}],
        "n_expected": 1,
        "n_success": 0,
        "status": "no data",
    }


def test_ont_repeat_keeps_the_ont_contract(tmp_path):
    import yaml

    from samovar.iss_config import setup_iss_test
    from samovar.repeats import repeat_seed

    genomes = Path("data/test_genomes/meta")
    batch = tmp_path / "ont"
    setup_iss_test(
        _generate_args(
            tmp_path,
            batch,
            repeat=1,
            genome_dir=str(genomes.resolve()),
            reads_generator="ont",
            simulator="camisim",
        )
    )
    cfg = yaml.safe_load((batch / "out_0001" / ".generate" / "configs" / "camisim.yaml").read_text())
    assert str(cfg.get("camisim_mode") or cfg.get("mode")) == "ont"
    assert cfg["seed"] == repeat_seed(42, 1, "community")


def _profile(taxids, columns):
    import pandas as pd

    frame = {"taxid": taxids}
    frame.update(columns)
    return pd.DataFrame(frame)


def test_long_tables_keep_successful_repeats(tmp_path):
    import csv
    import math

    from samovar.profiles import write_profiles
    from samovar.repeats import publish_package, set_repeat_state, write_long_tables

    root = tmp_path / "pkg"
    publish_package(root, count=4, base_seed=42, kind="generated", identity="generated:long")
    species = dict(kind="abundance", unit="counts", basis="reads", taxonomy="ncbi",
                   rank="species", normalization="all")
    genus = dict(kind="abundance", unit="counts", basis="reads", taxonomy="ncbi",
                 rank="genus", normalization="classified")
    presence = dict(kind="presence", unit="binary", basis="presence", taxonomy="ncbi",
                    rank="species", normalization="classified")

    stored = {}

    def put(repeat, phase, tool, table, meta):
        stored.setdefault((repeat, phase), {})[tool] = (table, meta)

    for repeat, a_562, a_4932 in (("0001", 10, 5), ("0002", 7, 1), ("0003", 8, None)):
        columns = {"N_A": [a_562, a_4932] if a_4932 is not None else [a_562]}
        taxids = ["562", "4932"] if a_4932 is not None else ["562"]
        if repeat == "0001":
            taxids = ["562", "4932", "0"]
            columns = {"N_A": [10, 5, 1], "N_B": [4, 6, 2]}
        elif repeat == "0002":
            columns["N_B"] = [3, 9]
        put(repeat, "initial", "kraken", _profile(taxids, columns), species)
        sylph_cols = {"N_A": [1 if repeat != "0003" else 3]}
        if repeat != "0003":
            sylph_cols["N_B"] = [1]
        put(repeat, "initial", "sylph", _profile(["2"], sylph_cols), genus)
        set_repeat_state(root, repeat, "succeeded", stage="pipeline")
    put("0001", "initial", "bracken", _profile(["562"], {"N_A": [4]}), species)
    put("0001", "initial", "present", _profile(["562"], {"N_A": [1]}), presence)
    put("0001", "regenerated", "kraken", _profile(["562"], {"N_A": [11]}), species)
    put("0004", "initial", "kraken", _profile(["999"], {"N_A": [50], "N_B": [50]}), species)
    for (repeat, phase), tools in stored.items():
        write_profiles(
            root / f"out_{repeat}" / f"{phase}_profiles",
            {name: table for name, (table, _meta) in tools.items()},
            {name: meta for name, (_table, meta) in tools.items()},
        )
    set_repeat_state(root, "0004", "failed", "exec failed", "pipeline")
    plots = root / "out_0001" / "initial_annotations_plots"
    plots.mkdir()
    (plots / "profile_quality_scores.csv").write_text(
        "annotator,reason\nkraken,no profile gold standard available\n", encoding="utf-8"
    )
    runtime = root / "out_0001" / "initial_reports" / "runtime" / "samples" / "kraken" / "classify"
    runtime.mkdir(parents=True)
    for sample, seconds, rss in (("A", 1.5, 10), ("B", 2.5, 20)):
        (runtime / f"{sample}.tsv").write_text(f"s\tmax_rss\n{seconds}\t{rss}\n", encoding="utf-8")
    (root / "out_0001" / "runtime" / "stages").mkdir(parents=True)
    (root / "out_0001" / "runtime" / "stages" / "viz.json").write_text(
        '{"stage": "viz_initial", "elapsed_s": 4}\n', encoding="utf-8"
    )
    marker = root / "out_0001" / ".log" / "checkpoints" / "keep"
    marker.parent.mkdir(parents=True)
    marker.write_text("stay", encoding="utf-8")
    before = marker.read_bytes()

    written = write_long_tables(root)
    first = written["profiles"].read_bytes()
    assert write_long_tables(root)["profiles"].read_bytes() == first
    assert marker.read_bytes() == before
    assert not any(row["taxon"] == "999" for row in csv.DictReader(first.decode().splitlines()))

    with written["profiles"].open(encoding="utf-8", newline="") as handle:
        profiles = list(csv.DictReader(handle))
    kraken = [row for row in profiles if row["run_name"] == "kraken" and row["phase"] == "initial"]
    assert {(row["repeat_id"], row["sample_id"]) for row in kraken} == {
        ("0001", "A"), ("0001", "B"), ("0002", "A"), ("0002", "B"), ("0003", "A"),
    }
    filled = next(row for row in kraken if row["repeat_id"] == "0003" and row["taxon"] == "4932")
    assert float(filled["value"]) == 0
    assert not any(row["repeat_id"] == "0003" and row["sample_id"] == "B" and row["run_name"] == "kraken" for row in profiles)
    assert not any(row["run_name"] == "bracken" and row["repeat_id"] != "0001" for row in profiles)
    unclassified = next(row for row in kraken if row["repeat_id"] == "0001" and row["sample_id"] == "A" and row["taxon"] == "0")
    assert float(unclassified["value"]) == 1
    assert not any(row["repeat_id"] == "0003" and row["taxon"] == "0" for row in kraken)
    sylph = [row for row in profiles if row["run_name"] == "sylph"]
    assert {row["rank"] for row in sylph} == {"genus"}
    assert {row["rank"] for row in kraken} == {"species"}
    assert sum(1 for row in sylph if row["sample_id"] == "A") == 3
    present = [row for row in profiles if row["run_name"] == "present"]
    assert [(row["taxon"], row["value"]) for row in present] == [("562", "1.0")]
    regenerated = [row for row in profiles if row["phase"] == "regenerated"]
    assert [(row["repeat_id"], row["sample_id"], float(row["value"])) for row in regenerated] == [("0001", "A", 11)]
    assert {float(row["value"]) for row in kraken if row["taxon"] == "562" and row["sample_id"] == "A"} == {10, 7, 8}

    with written["results"].open(encoding="utf-8", newline="") as handle:
        results = list(csv.DictReader(handle))
    absent = [row for row in results if row["reason"] == "no profile gold standard available"]
    assert len(absent) == 1 and absent[0]["value"] == "" and absent[0]["metric"] == ""
    assert absent[0]["rank"] == "species" and absent[0]["normalization"] == "all"
    rss = [float(row["value"]) for row in results if row["metric"] == "max_rss_mb"]
    assert rss == [10.0, 20.0]
    assert [row["reason"] for row in results if row["metric"] == "stage_elapsed_s"] == ["viz_initial"]
    with written["statistics"].open(encoding="utf-8", newline="") as handle:
        stats = list(csv.DictReader(handle))
    taxon_562 = next(
        row for row in stats
        if row["run_name"] == "kraken" and row["sample_id"] == "A"
        and row["phase"] == "initial" and row["taxon"] == "562"
    )
    assert math.isclose(float(taxon_562["mean"]), (10 + 7 + 8) / 3)
    assert (taxon_562["n_expected"], taxon_562["n_success"], taxon_562["n_valid"]) == ("4", "3", "3")
    assert not any(row["taxon"] == "999" for row in stats)


def test_statistics_uses_whole_repeats_and_known_moments():
    import math

    import numpy as np

    from samovar.repeats import repeat_seed, statistics_rows, summary_seed

    def row(repeat, taxon, value, **extra):
        base = {
            "repeat_id": repeat, "sample_id": "A", "phase": "initial", "run_name": "kraken",
            "kind": "abundance", "taxonomy": "ncbi", "rank": "species", "unit": "counts",
            "basis": "reads", "normalization": "all", "taxon": taxon, "value": value,
        }
        base.update(extra)
        return base

    rows = [
        row("0003", "562", 3),
        row("0001", "562", 1),
        row("0002", "562", 2),
        row("0001", "one", 7),
        row("0001", "bad", "inf"),
        row("0002", "bad", "nan"),
        row("0003", "bad", ""),
        row("0001", "flat", 5),
        row("0002", "flat", 5),
        row("0003", "flat", 5),
        row("0001", "pair", 10),
        row("0002", "pair", 14),
        row("0001", "split", 1, sample_id="A"),
        row("0001", "split", 100, sample_id="B"),
        row("0001", "dup", 4),
        row("0001", "dup", 4),
        row("0002", "dup", 6),
    ]
    kwargs = dict(
        n_expected=4, n_success=3, base_seed=42,
        input_kind="generated", fixed_abundance=True,
    )
    first = {item["taxon"] + item["sample_id"]: item for item in statistics_rows(rows, [], **kwargs)}
    second = {item["taxon"] + item["sample_id"]: item for item in statistics_rows(list(reversed(rows)), [], **kwargs)}
    assert first["562A"]["ci_low"] == second["562A"]["ci_low"]
    assert first["562A"]["ci_high"] == second["562A"]["ci_high"]

    got = first["562A"]
    assert float(got["mean"]) == 2
    assert float(got["sd"]) == 1
    assert got["n_valid"] == 3
    assert "n_valid<5" in got["limitation"]
    assert "incomplete" in got["limitation"]
    assert "fixed abundance" in got["limitation"]
    key = (
        "A", "initial", "kraken", "abundance", "ncbi", "species", "counts",
        "reads", "all", "", "562", "",
    )
    seed = repeat_seed(summary_seed(42), 1, "\x1f".join(key))
    values = np.array([1.0, 2.0, 3.0])
    drawn = np.random.default_rng(seed).integers(0, 3, size=(10000, 3))
    low, high = np.quantile(values[drawn].mean(axis=1), [0.025, 0.975])
    assert got["ci_low"] == format(float(low), ".12g")
    assert got["ci_high"] == format(float(high), ".12g")
    assert got["n_resamples"] == 10000
    assert got["bootstrap_method"] == "percentile"

    assert float(first["oneA"]["mean"]) == 7
    assert first["oneA"]["sd"] == "" and first["oneA"]["ci_low"] == ""
    assert first["oneA"]["reason"] == "n_valid=1"
    assert first["badA"]["n_valid"] == 0 and first["badA"]["mean"] == ""
    assert first["badA"]["reason"] == "no-data"
    assert float(first["flatA"]["mean"]) == 5 and float(first["flatA"]["sd"]) == 0
    assert first["flatA"]["ci_low"] == first["flatA"]["ci_high"]
    assert "zero-width" in first["flatA"]["limitation"]
    pair_sd = math.sqrt(((10 - 12) ** 2 + (14 - 12) ** 2) / 1)
    assert float(first["pairA"]["mean"]) == 12
    assert math.isclose(float(first["pairA"]["sd"]), pair_sd)
    assert float(first["splitA"]["mean"]) == 1 and float(first["splitB"]["mean"]) == 100
    assert first["dupA"]["n_valid"] == 2 and float(first["dupA"]["mean"]) == 5

    copies = statistics_rows(
        [row(f"{i:04d}", "562", 9) for i in range(1, 4)],
        [],
        n_expected=3, n_success=3, base_seed=42, input_kind="profiles",
    )
    assert copies[0]["n_valid"] == 0 and copies[0]["mean"] == ""
    assert "not an independent replicate" in copies[0]["reason"]
    plain = statistics_rows(
        rows[:3], [], n_expected=3, n_success=3, base_seed=42, input_kind="generated",
    )
    assert "population interval" in plain[0]["limitation"]
    assert "fixed abundance" not in plain[0]["limitation"]

    stages = statistics_rows([], [
        {"repeat_id": "0001", "metric": "stage_elapsed_s", "value": 4, "reason": "viz_initial"},
        {"repeat_id": "0001", "metric": "stage_elapsed_s", "value": 9, "reason": "qc_initial"},
    ], n_expected=1, n_success=1, base_seed=1, input_kind="generated")
    assert {item["run_name"]: float(item["mean"]) for item in stages} == {
        "viz_initial": 4, "qc_initial": 9,
    }


def test_public_cli_repeat_package_preserves_single_run_semantics(tmp_path, monkeypatch):
    import subprocess
    import sys

    import yaml

    from samovar.repeats import repeat_seed

    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "root": str(Path(__file__).resolve().parents[1]),
        "compilers": {"python": sys.executable},
        "tools": {},
        "custom-flags": [],
    }), encoding="utf-8")
    monkeypatch.setenv("SAMOVAR_CONFIG", str(config))
    samovar = Path(__file__).resolve().parents[1] / "bin" / "samovar"
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    abundance = tmp_path / "abundance.csv"
    original = "taxid,N_A,N_B\n562,10,0\n9606,0,4\n"
    abundance.write_text(original, encoding="utf-8")

    def generate(dest, *extra):
        result = subprocess.run(
            [str(samovar), "generate", "--genome_dir", str(genomes), "--host_genome", "",
             "--n_samples", "2", "--total_reads", "20", "--host_fraction", "0",
             "--seed", "42", "--model", "hiseq", "--autoimport", "0",
             "--directory", str(dest), *extra],
            check=False, text=True, capture_output=True,
        )
        assert result.returncode == 0, result.stderr
        return result

    single = tmp_path / "single"
    package = tmp_path / "batch dir"
    generate(single, "--abundance", str(abundance))
    generate(package, "--abundance", str(abundance), "--repeat", "2")
    assert abundance.read_text(encoding="utf-8") == original
    assert not (single / ".log" / "repeats.json").exists()
    single_cfg = yaml.safe_load((single / ".generate" / "configs" / "iss_config.yaml").read_text())
    child_cfg = yaml.safe_load(
        (package / "out_0001" / ".generate" / "configs" / "iss_config.yaml").read_text()
    )
    assert single_cfg["seed"] == 42
    assert child_cfg["seed"] == repeat_seed(42, 1, "reads")
    assert child_cfg["sample_names"] == ["A", "B"]
    assert Path(child_cfg["abundance_table"]).resolve() == abundance.resolve()
    assert "out_0001/initial" in child_cfg["output_dir"].replace("\\", "/")

    drawn = tmp_path / "drawn"
    generate(drawn, "--repeat", "1")
    drawn_cfg = yaml.safe_load((drawn / "out_0001" / ".generate" / "configs" / "iss_config.yaml").read_text())
    assert drawn_cfg["seed"] == repeat_seed(42, 1, "community")
    assert "abundance_table" not in drawn_cfg

    def prepare(dest, mode):
        result = subprocess.run(
            [str(samovar), "prepare", "--directory", str(dest), "--mode", mode,
             "--cores", "1", "--autoimport", "0"],
            check=False, text=True, capture_output=True,
        )
        assert result.returncode == 0, result.stderr

    prepare(package, "benchmark")
    assert not (package / ".log" / "samovar.sh").exists()
    for name in ("out_0001", "out_0002"):
        window = (package / name / ".log" / "window.env").read_text()
        script = (package / name / ".log" / "samovar.sh").read_text()
        assert "viz_initial" in window
        assert "viz_reprofiled" not in window
        assert "trained_model.joblib" not in script
        assert name in script
    other = (package / "out_0001" / ".log" / "samovar.sh").read_text()
    assert "out_0002" not in other

    prepare(drawn, "full")
    full_window = (drawn / "out_0001" / ".log" / "window.env").read_text()
    assert "viz_reprofiled" in full_window

    refused = subprocess.run(
        [str(samovar), "prepare", "--directory", str(package), "--repeat", "9", "--autoimport", "0"],
        check=False, text=True, capture_output=True,
    )
    assert refused.returncode != 0
    assert not (package / "out_0009").exists()
