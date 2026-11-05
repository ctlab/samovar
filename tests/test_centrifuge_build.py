"""samovar build --type centrifuge wraps centrifuge-build."""

from pathlib import Path
from unittest.mock import patch

import yaml

from samovar.build_database import build_database_from_config

import pytest

# Centrifuge-specific database adapter.
pytestmark = pytest.mark.optional


def test_centrifuge_build_writes_library_and_calls_centrifuge_build(tmp_path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "562.fna").write_text(">ecoli extra\nATGCATGC\n", encoding="utf-8")
    (genomes / "4932.faa").write_text(">prot\nMKTLLL\n", encoding="utf-8")
    cfg = tmp_path / "db.yaml"
    cfg.write_text(yaml.dump({"input_dir": [str(genomes)]}), encoding="utf-8")
    db = tmp_path / "centrifuge_db"
    seen = {}

    def fake_tax(db_path="kraken_db", taxonomy_path=None):
        seen["taxonomy"] = db_path
        Path(db_path).mkdir(parents=True, exist_ok=True)
        (Path(db_path) / "nodes.dmp").write_text("1\t|\t1\t|\n", encoding="utf-8")
        (Path(db_path) / "names.dmp").write_text("1\t|\troot\t|\n", encoding="utf-8")

    def fake_run(cmd, check=True):
        seen["cmd"] = cmd
        return None

    with patch("samovar.build_database.get_taxonomy_db", side_effect=fake_tax), patch(
        "samovar.build_database.run_command", side_effect=fake_run
    ):
        build_database_from_config(
            str(cfg), db_type="centrifuge", db_path=str(db), threads=3
        )

    library = (db / "library.fna").read_text(encoding="utf-8")
    assert library.startswith(">ecoli\nATGCATGC")
    assert "MKT" not in library
    assert (db / "seqid2taxid.map").read_text(encoding="utf-8") == "ecoli\t562\n"
    assert seen["taxonomy"].endswith("taxonomy")
    cmd = seen["cmd"]
    assert cmd[0] == "centrifuge-build"
    assert cmd[cmd.index("-p") + 1] == "3"
    assert cmd[cmd.index("--conversion-table") + 1] == str(db / "seqid2taxid.map")
    assert cmd[cmd.index("--taxonomy-tree") + 1] == str(db / "taxonomy" / "nodes.dmp")
    assert cmd[cmd.index("--name-table") + 1] == str(db / "taxonomy" / "names.dmp")
    assert cmd[-2] == str(db / "library.fna")
    assert cmd[-1] == str(db / "centrifuge")


def test_centrifuge_build_rejects_proteome_only(tmp_path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "4932.faa").write_text(">prot\nMKTLLL\n", encoding="utf-8")
    cfg = tmp_path / "db.yaml"
    cfg.write_text(yaml.dump({"input_dir": [str(genomes)]}), encoding="utf-8")
    try:
        build_database_from_config(
            str(cfg), db_type="centrifuge", db_path=str(tmp_path / "db")
        )
    except ValueError as exc:
        assert "nucleotide" in str(exc)
    else:
        raise AssertionError("proteome-only centrifuge build should fail")
