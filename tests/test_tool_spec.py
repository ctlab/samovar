import argparse

from samovar.main_config import disk_payload, migrate_legacy, parse_tool_entry
from samovar.tool_spec import (
    apply_translated_flags,
    parse_flags_translate,
    parse_known_with_custom,
    parse_tool_record,
    record_to_spec,
)


def test_parse_flags_translate_cli_string():
    mapping = parse_flags_translate("--threads:--threads --cores:-t")
    assert mapping["--threads"] == "--threads"
    assert mapping["--cores"] == "-t"


def test_record_roundtrip_object_and_legacy_list():
    rec = parse_tool_record(
        ["", "bash", "/usr/bin/kraken2", "annotator"],
        "kraken2:2.1.3",
    )
    assert rec["exec"]["path"] == "/usr/bin/kraken2"
    assert rec["type"] == "annotator"
    assert rec["lazy-install"]
    assert rec["flags-translate"]["--threads"] == "--threads"
    spec = record_to_spec(rec, "kraken2")
    assert parse_tool_entry(spec, "kraken2")[2].endswith("kraken2")


def test_apply_translated_flags_kaiju():
    extra = apply_translated_flags(
        "",
        name="kaiju",
        canonical={"--threads": 8, "--cores": 8},
    )
    assert "-z 8" in extra


def test_parse_known_with_custom_accepts_threads():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir")
    args, leftover = parse_known_with_custom(
        parser, ["--output_dir", "out", "--threads", "4"]
    )
    assert leftover == []
    assert int(args.threads) == 4
    assert args.custom_cli_flags.get("--threads") == "4"


def test_disk_payload_versions_tool_keys():
    cfg = migrate_legacy(
        {
            "root": "/opt/samovar",
            "tools": {
                "iss": {
                    "exec": {"env": "", "parser": "bash", "path": "/usr/bin/iss"},
                    "type": "reads_generator",
                    "lazy-install": "pip install insilicoseq",
                    "flags": "",
                    "flags-translate": {"--threads": "--cpus"},
                }
            },
        }
    )
    payload = disk_payload(cfg)
    keys = list(payload["tools"])
    assert any(k == "iss" or str(k).startswith("iss:") for k in keys)
    rec = next(payload["tools"][k] for k in keys if k == "iss" or str(k).startswith("iss:"))
    assert rec["flags-translate"]["--threads"] == "--cpus"
    assert "custom-flags" in payload


def test_prepare_threads_reach_iss_and_snakemake(tmp_path, monkeypatch):
    from pathlib import Path

    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config(
        {"root": str(tmp_path), "tools": {}, "custom-flags": ["--threads", "--cores"]},
        also_repo_build=False,
    )
    (tmp_path / "reads").mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(tmp_path / "reads"),
            "output_dir": str(tmp_path / "out"),
            "threads": 6,
            "cores": 1,
            "custom_cli_flags": {"--threads": "6"},
        },
    )()
    config = PipelineConfig.from_args(args)
    assert config.cores == 6
    paths = config.generate_configs(str(tmp_path / "out"))
    a2iss = Path(paths["annotation2iss"]).read_text()
    assert "cores: 6" in a2iss
    script = Path(config.generate_pipeline(str(tmp_path / "out"))).read_text()
    assert "--cores 6" in script
