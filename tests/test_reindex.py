"""``samovar reindex`` moves processed FASTAs into samovar_database and updates config."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

import pytest

from samovar.genome_index import (
    index_processed_file,
    iter_processed_fastas,
    processed_filename,
    reindex,
    resolve_indexed_path,
    samovar_database_dir,
)
from samovar.paths import load_config, write_install_config_pointer


def _isolate(tmp_path, monkeypatch):
    root = tmp_path / "install"
    root.mkdir()
    cfg = tmp_path / "cfg" / "config.json"
    cfg.parent.mkdir()
    store = tmp_path / "samovar_db"
    payload = {
        "version": "0.0.0",
        "root": str(root),
        "compilers": {"python": "/usr/bin/python3"},
        "API": {},
        "genomes": {
            "samovar_database": str(store),
            "test": [],
            "raw": {},
            "processed": {"samovar_database": str(store / "processed")},
            "data": {},
        },
        "databases": {},
        "workflows": {},
        "tools": {},
    }
    cfg.write_text(json.dumps(payload) + "\n")
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    monkeypatch.setenv("SAMOVAR_ROOT", str(root))
    monkeypatch.setenv("SAMOVAR_DATABASE", str(store))
    write_install_config_pointer(cfg, root=root)
    return cfg, store


def _write_fa(path: Path, seq: str = "ACGT") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz" or path.name.endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(f">g\n{seq}\n")
    else:
        path.write_text(f">g\n{seq}\n")
    return path


def test_reindex_moves_processed_and_updates_config(tmp_path, monkeypatch):
    cfg, store = _isolate(tmp_path, monkeypatch)
    src = tmp_path / "run" / ".genomes" / "processed"
    plain = _write_fa(src / "GCF_000819615.1.fa", "AAAA")
    result = reindex([tmp_path / "run"], dest=store / "processed")
    assert result["count"] == 1
    dest = store / "processed" / "GCF_000819615.1.fa.gz"
    assert dest.is_file()
    assert not plain.exists()
    assert resolve_indexed_path("GCF_000819615.1") == dest
    loaded = load_config()
    rec = loaded["genomes"]["data"]["GCF_000819615.1"]
    assert rec[1] == "GCF_000819615.1"
    assert rec[2] == "samovar_database"
    assert rec[3] == "GCF_000819615.1.fa.gz"


def test_reindex_errors_when_processed_empty(tmp_path, monkeypatch):
    _isolate(tmp_path, monkeypatch)
    empty = tmp_path / "empty"
    empty.mkdir()
    (empty / "processed").mkdir()
    with pytest.raises(RuntimeError, match="nothing in processed"):
        reindex([empty])


def test_reindex_no_args_moves_indexed_off_dest(tmp_path, monkeypatch):
    cfg, store = _isolate(tmp_path, monkeypatch)
    stray = tmp_path / "other" / "processed"
    gz = _write_fa(stray / "GCF_000840245.1.fa.gz", "CCCC")
    index_processed_file(gz, accession="GCF_000840245.1", move_to=stray)
    assert resolve_indexed_path("GCF_000840245.1") == stray / "GCF_000840245.1.fa.gz"
    result = reindex(None, dest=store / "processed")
    assert result["count"] == 1
    dest = store / "processed" / "GCF_000840245.1.fa.gz"
    assert dest.is_file()
    assert not gz.exists()
    assert resolve_indexed_path("GCF_000840245.1") == dest


def test_reindex_no_args_errors_without_processed(tmp_path, monkeypatch):
    _isolate(tmp_path, monkeypatch)
    with pytest.raises(RuntimeError, match="nothing in processed"):
        reindex(None)


def test_iter_processed_fastas_nested(tmp_path):
    folder = tmp_path / "run" / ".genomes" / "processed"
    _write_fa(folder / "GCF_000836945.1.fa.gz")
    found = iter_processed_fastas(tmp_path / "run")
    assert [p.name for p in found] == ["GCF_000836945.1.fa.gz"]
    assert processed_filename("GCF_000836945.1") == "GCF_000836945.1.fa.gz"


def test_samovar_reindex_cli(tmp_path, monkeypatch):
    cfg, store = _isolate(tmp_path, monkeypatch)
    src = tmp_path / "run" / "processed"
    _write_fa(src / "GCF_000844825.1.fa.gz")
    from samovar.genome_index import main

    assert main([str(tmp_path / "run"), "--dest", str(store / "processed")]) == 0
    assert (store / "processed" / "GCF_000844825.1.fa.gz").is_file()
    pointer = tmp_path / "install" / "build" / "config_path"
    assert pointer.read_text().strip() == str(cfg.resolve())
    assert samovar_database_dir() == store


def test_reindex_skips_already_stored_hardlinks(tmp_path, monkeypatch):
    cfg, store = _isolate(tmp_path, monkeypatch)
    dest = store / "processed"
    dest.mkdir(parents=True)
    stored = _write_fa(dest / "GCF_000819615.1.fa.gz", "AAAA")
    index_processed_file(stored, accession="GCF_000819615.1", move_to=dest)
    run = tmp_path / "run" / ".genomes" / "processed"
    run.mkdir(parents=True)
    link = run / "GCF_000819615.1.fa.gz"
    try:
        os.link(stored, link)
    except OSError:
        import shutil

        shutil.copy2(stored, link)
    fresh = _write_fa(run / "GCF_000867865.1.fa.gz", "TTTT")
    result = reindex([tmp_path / "run"], dest=dest)
    assert result["count"] == 1
    assert (dest / "GCF_000867865.1.fa.gz").is_file()
    assert stored.is_file()


def test_materialize_reindex_modes(tmp_path, monkeypatch):
    from unittest.mock import patch

    from samovar.genome_fetcher import materialize_accessions
    from samovar.genome_index import genome_data_map, resolve_indexed_path

    _isolate(tmp_path, monkeypatch)
    store = tmp_path / "samovar_db"
    acc = "GCF_000819615.1"

    def fake_fetch(accession, dest_dir, email, **kwargs):
        dest = Path(dest_dir) / f"{accession}.fa.gz"
        dest.parent.mkdir(parents=True, exist_ok=True)
        _write_fa(dest, "ACGT")
        if kwargs.get("index"):
            from samovar.genome_index import index_processed_file

            index_processed_file(dest, accession=accession, move_to=Path(dest_dir))
        return str(dest)

    out1 = tmp_path / "run1"
    with patch("samovar.genome_fetcher.fetch_assembly_processed", side_effect=fake_fetch):
        paths = materialize_accessions(
            [acc], output_dir=str(out1), email="t@e", reindex_mode=1
        )
    assert paths
    assert (store / "processed" / f"{acc}.fa.gz").is_file()
    assert (out1 / ".genomes" / "processed" / f"{acc}.fa.gz").is_file()
    assert resolve_indexed_path(acc) is not None

    out0 = tmp_path / "run0"
    acc2 = "GCF_000867865.1"
    with patch("samovar.genome_fetcher.fetch_assembly_processed", side_effect=fake_fetch):
        materialize_accessions(
            [acc, acc2], output_dir=str(out0), email="t@e", reindex_mode=0
        )
    assert acc2 not in genome_data_map()
    assert (out0 / ".genomes" / "processed" / f"{acc}.fa.gz").is_file()
    assert (out0 / ".genomes" / "processed" / f"{acc2}.fa.gz").is_file()
    # Indexed genome was not re-fetched into the store as a second copy name.
    assert resolve_indexed_path(acc).name == f"{acc}.fa.gz"


def test_expand_named_database(tmp_path, monkeypatch):
    from samovar.config import _expand_indexed_database
    from samovar.genome_index import register_database

    _isolate(tmp_path, monkeypatch)
    db = tmp_path / "kaiju_db"
    db.mkdir()
    register_database("kaiju", "phage_test", str(db), "--threads 2")
    path, extra = _expand_indexed_database("kaiju", "phage_test", None)
    assert path == str(db)
    assert extra == "--threads 2"
    path2, extra2 = _expand_indexed_database("kaiju", str(db), "x")
    assert path2 == str(db)
    assert extra2 == "x"


def test_normalize_genome_data_legacy_and_current():
    from samovar.genome_index import normalize_genome_data

    legacy_acc = normalize_genome_data(
        {"GCF_000819615.1": ["GCF_000819615.1", "samovar_database", "GCF_000819615.1.fa.gz"]}
    )
    assert legacy_acc["GCF_000819615.1"] == [
        "",
        "GCF_000819615.1",
        "samovar_database",
        "GCF_000819615.1.fa.gz",
    ]
    with_tax = normalize_genome_data(
        {
            "GCF_000819615.1": [
                "GCF_000819615.1",
                "samovar_database",
                "GCF_000819615.1.fa.gz",
                "10847",
            ]
        }
    )
    assert with_tax["10847"] == [
        "10847",
        "GCF_000819615.1",
        "samovar_database",
        "GCF_000819615.1.fa.gz",
    ]
    stub = normalize_genome_data({"562": ["", "test", "meta/562.fna"]})
    assert stub["562"] == ["562", "562", "test", "meta/562.fna"]
    current = normalize_genome_data(
        {"10847": ["10847", "GCF_000819615.1", "samovar_database", "GCF_000819615.1.fa.gz"]}
    )
    assert current["10847"][1] == "GCF_000819615.1"
    assert current["10847"][3] == "GCF_000819615.1.fa.gz"
