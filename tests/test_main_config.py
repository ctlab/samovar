from pathlib import Path

from samovar.main_config import (
    as_install_config,
    build_install_config,
    disk_payload,
    extra_genome_dirs,
    iter_tools,
    migrate_legacy,
    parse_tool_entry,
    tool_path,
)


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
    kraken = parse_tool_entry(payload["tools"]["kraken2"], "kraken2")
    assert kraken[2].endswith("kraken2")
    assert kraken[3] == "annotator"


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
        version="0.10.20",
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
    }
    assert payload["compilers"]["python"] == str(py)
    assert "python_path" not in payload
    assert extra_genome_dirs(cfg) == extra_genome_dirs(payload)


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
        version="0.10.20",
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
        version="0.10.20",
        extra_genome_dirs=[str(linked), str(scratch)],
        genomes_default=str(scratch),
    )
    payload = disk_payload(cfg)
    raw_paths = list(payload["genomes"]["raw"].values())
    assert raw_paths == [str(scratch)]
    assert not any(str(home) in p for p in raw_paths)
