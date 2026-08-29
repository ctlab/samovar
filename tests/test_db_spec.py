import json

from samovar.config import PipelineConfig, _expand_indexed_database
from samovar.db_spec import (
    databases_to_rows,
    iter_database_records,
    lazy_download_for,
    parse_database_record,
)
from samovar.genome_index import register_database
from samovar.paths import write_config, update_config
from samovar.tools_import import import_database, main as import_main


def test_parse_legacy_list_and_object():
    rec = parse_database_record(["phage_test", "/db", "--threads 2"], tool="kaiju")
    assert rec["name"] == "phage_test"
    assert rec["path"] == "/db"
    assert rec["flags"] == "--threads 2"
    obj = parse_database_record(
        {"path": "/k", "flags": "--memory-mapping", "lazy-download": "curl x"},
        tool="kraken2",
        name="standard_8GB:2025oct",
    )
    assert obj["name"] == "standard_8GB"
    assert obj["_version"] == "2025oct"
    assert "curl" in obj["lazy-download"]


def test_iter_and_disk_roundtrip():
    cfg = {
        "databases": {
            "kaiju": [["phage_test", "/db", ""]],
            "kraken2": {
                "standard_8GB:2025oct": {
                    "path": "/k",
                    "flags": "",
                    "lazy-download": "curl k2",
                }
            },
        }
    }
    rows = databases_to_rows(cfg)
    assert rows["kaiju"][0][0] == "phage_test"
    assert rows["kraken2"][0][0] in {"standard_8GB:2025oct", "standard_8GB"}
    grouped = iter_database_records(cfg)
    assert grouped["kraken2"]["standard_8GB:2025oct"]["path"] == "/k"


def test_set_database_and_expand(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg_path))
    monkeypatch.setattr(
        "samovar.paths.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    write_config({"root": str(tmp_path), "tools": {}, "databases": {}}, also_repo_build=False)
    db = tmp_path / "kaiju_db"
    db.mkdir()
    register_database("kaiju", "phage_test", str(db), "--threads 2")
    stored = json.loads(cfg_path.read_text())["databases"]
    assert stored["kaiju"]["phage_test"]["path"] == str(db)
    assert stored["kaiju"]["phage_test"]["flags"] == "--threads 2"
    path, extra = _expand_indexed_database("kaiju", "phage_test", None)
    assert path == str(db)
    assert extra == "--threads 2"


def test_prepare_from_args_expands_indexed_db(tmp_path, monkeypatch):
    import argparse

    cfg_path = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg_path))
    monkeypatch.setattr(
        "samovar.paths.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    write_config({"root": str(tmp_path), "tools": {}, "databases": {}}, also_repo_build=False)
    db = tmp_path / "k2"
    db.mkdir()
    register_database("kraken2", "standard_8GB", str(db), "--memory-mapping", version="2025oct")
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(tmp_path),
        output_dir=str(tmp_path / "out"),
        kraken2=None,
        kaiju=None,
        dummy=None,
    )
    setattr(args, "cmd_kraken2-test", [["kraken2 standard_8GB"]])
    config = PipelineConfig.from_args(args)
    k2 = next(a for a in config.annotators if a.type == "kraken2")
    assert k2.db_path == str(db.resolve()) or k2.db_path == str(db)
    assert "--memory-mapping" in (k2.extra or "")


def test_import_database_cli(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    write_config({"root": str(tmp_path), "tools": {}, "databases": {}}, also_repo_build=False)
    db = tmp_path / "k2"
    db.mkdir()
    (db / "hash.k2d").write_text("x")
    rec = import_database(
        name="standard_8GB",
        tool="kraken2",
        exec_path=str(db),
        version="2025oct",
        flags="--memory-mapping",
        also_repo_build=False,
    )
    assert rec["name"] == "standard_8GB"
    tools_cfg = json.loads(cfg.read_text())
    recs = tools_cfg["databases"]["kraken2"]["standard_8GB:2025oct"]
    assert recs["path"] == str(db.resolve())
    assert recs["flags"] == "--memory-mapping"
    assert recs["lazy-download"]
    rc = import_main(
        [
            "-n",
            "kraken2:virus",
            "--type",
            "database",
            "--exec-path",
            str(db),
            "--tool-version",
            "2025oct",
        ]
    )
    assert rc == 0
    payload = json.loads(cfg.read_text())["databases"]["kraken2"]
    assert "virus:2025oct" in payload


def test_metaphlan_flags_injected(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    monkeypatch.setattr(
        "samovar.paths.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    write_config({"root": str(tmp_path), "databases": {}}, also_repo_build=False)
    db = tmp_path / "mpa"
    db.mkdir()
    (db / "mpa.pkl").write_text("x")
    register_database("metaphlan", "toy", str(db), "")
    path, extra = _expand_indexed_database("metaphlan", "toy", None)
    assert path == str(db)
    assert "--bowtie2db" in extra
    assert str(db) in extra


def test_lazy_download_known_url():
    recipe = lazy_download_for("kraken2", "standard_8GB", "2025oct")
    assert "k2_standard_08gb_20251015" in recipe
    assert "curl" in recipe
    assert "tar -xzf" in recipe


def test_lazy_download_plain_tar_and_versioned_qiime():
    mpa = lazy_download_for("metaphlan", "jan25", "2025jan")
    assert "mpa_vJan25" in mpa
    assert "tar -xf " in mpa
    assert "tar -xzf" not in mpa
    qza = lazy_download_for("qiime", "silva-138-99", "2024.10")
    assert "2024.10/common/silva-138-99-nb-classifier.qza" in qza
    assert "tar " not in qza


def test_two_versions_kept(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg_path))
    monkeypatch.setattr(
        "samovar.paths.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    write_config({"root": str(tmp_path), "tools": {}, "databases": {}}, also_repo_build=False)
    a = tmp_path / "silva_a"
    b = tmp_path / "silva_b"
    a.write_text("a")
    b.write_text("b")
    register_database("qiime", "silva-138-99", str(a), "", version="2024.2")
    register_database("qiime", "silva-138-99", str(b), "", version="2024.10")
    stored = json.loads(cfg_path.read_text())["databases"]["qiime"]
    assert stored["silva-138-99:2024.2"]["path"] == str(a)
    assert stored["silva-138-99:2024.10"]["path"] == str(b)
    assert stored["silva-138-99:2024.10"]["url"].endswith("2024.10/common/silva-138-99-nb-classifier.qza")
