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
