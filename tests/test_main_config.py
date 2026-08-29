from pathlib import Path

from samovar.main_config import (
    as_install_config,
    build_install_config,
    disk_payload,
    extra_genome_dirs,
    format_install_report,
    iter_tools,
    migrate_legacy,
    parse_tool_entry,
    sync_by_keys,
    tool_inputs,
    tool_path,
)
from samovar.version import get_version


def test_migrate_drops_duplicate_flat_keys():
    raw = {
        "version": "0.10.11",
        "root": "/opt/samovar",
        "python_path": "/opt/conda/bin/python",
        "iss_path": "/opt/conda/bin/iss",
        "ncbi_email": "a@b.c",
        "test_genomes": "/opt/samovar/data/test_genomes",
        "genomes": "/scratch/genomes",
        "processed_genomes": "/scratch/genomes",
        "genome_dirs": ["/scratch/genomes"],
        "tools": {
            "python": "/opt/conda/bin/python",
            "iss": "/opt/conda/bin/iss",
            "kraken2": "/opt/conda/bin/kraken2",
        },
        "tool_envs": {},
    }
    cfg = migrate_legacy(raw)
    payload = disk_payload(cfg)
    assert "python_path" not in payload
    assert "iss_path" not in payload
    assert "tool_envs" not in payload
    assert payload["compilers"]["python"] == "/opt/conda/bin/python"
    assert payload["API"]["ncbi_email"] == "a@b.c"
    assert payload["genomes"]["raw"]["default"] == "/scratch/genomes"
    kraken_key = next(
        k for k in payload["tools"] if k == "kraken2" or str(k).startswith("kraken2:")
    )
    rec = payload["tools"][kraken_key]
    assert isinstance(rec, dict)
    assert rec["type"] == "annotator"
    assert rec.get("lazy-install")
    assert "--threads" in (rec.get("flags-translate") or {})
    kraken = parse_tool_entry(rec, "kraken2")
    assert kraken[2].endswith("kraken2")
    assert kraken[3] == "annotator"
    assert payload.get("custom-flags")


def test_parse_tool_entry_keeps_empty_flags_when_inputs_set():
    spec = parse_tool_entry(
        ["", "bash", "/opt/counts.py", "scoring", "", "*annotations"],
        "counts",
    )
    assert spec == ["", "bash", "/opt/counts.py", "scoring", "", "*annotations"]
    assert tool_inputs(spec, "counts") == "*annotations"
    table = parse_tool_entry(
        {
            "env": "",
            "workflow": "bash",
            "path": "/opt/table.py",
            "group": "scoring",
            "inputs": "combined_annotation_table.csv",
        },
        "table_score",
    )
    assert table[5] == "combined_annotation_table.csv"
    assert table[4] == ""
    four = parse_tool_entry(["", "bash", "/usr/bin/iss", "reads_generator"], "iss")
    assert len(four) == 4


def test_install_config_legacy_get():
    cfg = as_install_config(
        {
            "compilers": {"python": "/usr/bin/python3"},
            "API": {"ncbi_email": "x@y.z"},
            "tools": {"iss": ["", "bash", "/usr/bin/iss", "reads_generator"]},
        }
    )
    assert cfg.get("python_path") == "/usr/bin/python3"
    assert cfg.get("ncbi_email") == "x@y.z"
    assert cfg.get("iss_path") == "/usr/bin/iss"


def test_set_tool_via_legacy_assignment():
    cfg = as_install_config({})
    cfg["opal_path"] = "/usr/bin/opal.py"
    assert tool_path(iter_tools(cfg)["opal.py"], "opal.py") == "/usr/bin/opal.py"


def test_build_install_config_nested(tmp_path):
    py = tmp_path / "python"
    py.write_text("#!/bin/sh\n")
    py.chmod(0o755)
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(py),
        version=get_version(),
        discovered_tools={"iss": str(tmp_path / "iss")},
        ncbi_email="n@e.c",
    )
    payload = disk_payload(cfg)
    assert set(payload) <= {
        "version",
        "root",
        "compilers",
        "API",
        "genomes",
        "databases",
        "workflows",
        "tools",
        "custom-flags",
    }
    assert payload["version"] == get_version()
    assert payload["compilers"]["python"] == str(py)
    assert "python_path" not in payload
    assert extra_genome_dirs(cfg) == extra_genome_dirs(payload)


def test_sync_by_keys_fills_missing_and_keeps_existing():
    template = {
        "compilers": {"python": "", "cpp_libs": []},
        "genomes": {"samovar_database": "/new", "data": {}},
        "databases": {},
    }
    existing = {
        "compilers": {"python": "/usr/bin/python3"},
        "genomes": {
            "samovar_database": "/old/store",
            "data": {"562": ["", "test", "meta/562.fna"]},
        },
        "databases": {"kaiju": [["phage_test", "/db", ""]]},
    }
    merged = sync_by_keys(existing, template)
    assert merged["compilers"]["python"] == "/usr/bin/python3"
    assert merged["compilers"]["cpp_libs"] == []
    assert merged["genomes"]["samovar_database"] == "/old/store"
    assert merged["databases"]["kaiju"][0][0] == "phage_test"


def test_install_migrates_legacy_genome_data_and_keeps_catalog(tmp_path):
    py = tmp_path / "python"
    py.write_text("#!/bin/sh\n")
    py.chmod(0o755)
    existing = {
        "root": str(tmp_path),
        "compilers": {"python": str(py)},
        "genomes": {
            "samovar_database": str(tmp_path / "store"),
            "data": {
                "GCF_000819615.1": [
                    "GCF_000819615.1",
                    "samovar_database",
                    "GCF_000819615.1.fa.gz",
                    "10847",
                ]
            },
        },
        "databases": {"kraken2": [["phage_test", "/k", ""]]},
    }
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(py),
        version=get_version(),
        existing=existing,
        samovar_database=str(tmp_path / "store"),
    )
    payload = disk_payload(cfg)
    rec = payload["genomes"]["data"]["10847"]
    assert rec == [
        "10847",
        "GCF_000819615.1",
        "samovar_database",
        "GCF_000819615.1.fa.gz",
    ]
    assert "GCF_000819615.1" not in payload["genomes"]["data"]
    assert payload["databases"]["kraken2"][0][0] == "phage_test"


def test_home_genome_dirs_are_dropped(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    scratch = tmp_path / "scratch" / "genomes"
    scratch.mkdir(parents=True)
    home_cache = home / ".cache" / "samovar" / "genomes"
    home_cache.mkdir(parents=True)
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(tmp_path / "python"),
        version=get_version(),
        extra_genome_dirs=[str(home_cache), str(scratch)],
        genomes_default=str(home_cache),
    )
    payload = disk_payload(cfg)
    raw_paths = list(payload["genomes"]["raw"].values())
    assert str(scratch) in raw_paths
    assert not any(str(home) in p for p in raw_paths)
    assert str(home_cache) not in extra_genome_dirs(payload)


def test_home_symlink_to_scratch_is_not_stored(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    scratch = tmp_path / "scratch" / "genomes"
    scratch.mkdir(parents=True)
    cache = home / ".cache" / "samovar"
    cache.parent.mkdir(parents=True)
    cache.symlink_to(scratch.parent)
    linked = cache / "genomes"
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(tmp_path / "python"),
        version=get_version(),
        extra_genome_dirs=[str(linked), str(scratch)],
        genomes_default=str(scratch),
    )
    payload = disk_payload(cfg)
    raw_paths = list(payload["genomes"]["raw"].values())
    assert raw_paths == [str(scratch)]
    assert not any(str(home) in p for p in raw_paths)


def test_format_install_report_first_install(tmp_path):
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(tmp_path / "python"),
        version=get_version(),
        ncbi_email="n@e.c",
        samovar_database=str(tmp_path / "genomes"),
    )
    text = format_install_report(
        payload=disk_payload(cfg),
        previous=None,
        previous_path=str(tmp_path / "config.json"),
    )
    assert "Previous installation config: not found" in text
    assert "writing" in text
    assert "New / updated options this install:" in text
    assert "genomes.samovar_database" in text
    assert "Options aligned from the previous install:" not in text


def test_format_install_report_keeps_previous_options(tmp_path):
    store = str(tmp_path / "library")
    previous = {
        "version": "0.0.1",
        "root": str(tmp_path),
        "genomes": {"samovar_database": store},
        "databases": {"kaiju": {"index": str(tmp_path / "kaiju.fmi")}},
        "tools": {"kraken2": str(tmp_path / "bin" / "kraken2")},
    }
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(tmp_path / "python"),
        version=get_version(),
        existing=previous,
        samovar_database=store,
        ncbi_email="new@e.c",
    )
    text = format_install_report(
        payload=disk_payload(cfg),
        previous=previous,
        previous_path=str(tmp_path / "config.json"),
    )
    assert "Previous installation config: found" in text
    assert "Options aligned from the previous install:" in text
    assert "databases.kaiju" in text
    assert "genomes.samovar_database" in text
    assert "New / updated options this install:" in text
    assert "API.ncbi_email" in text
    assert "new@e.c" in text
