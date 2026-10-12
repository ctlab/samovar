"""KrakenUniq library layout: uncompressed FASTA in library/added, dumps in taxonomy/."""

import gzip
import os
import stat
from pathlib import Path
from unittest.mock import patch

import yaml

from samovar.build_database import (
    KrakenUniqLibrarySkip,
    add_database_krakenunique,
    build_database_from_config,
    build_database_krakenunique,
)


def _nucleotide(path: Path, seq: str = "ATGCATGCATGC") -> None:
    path.write_text(f">seq\n{seq}\n", encoding="utf-8")


def test_krakenunique_writes_added_fasta_and_skips_proteome(tmp_path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    _nucleotide(genomes / "562.fna")
    gz_path = genomes / "9606.fna.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
        handle.write(">human\nAAAAAAAAAAAA\n")
    (genomes / "4932.faa").write_text(">prot\nMKTLLL\n", encoding="utf-8")
    cfg = tmp_path / "db.yaml"
    cfg.write_text(yaml.dump({"input_dir": [str(genomes)]}), encoding="utf-8")
    db = tmp_path / "krakenuniq_db"
    seen = {}

    def fake_tax(db_path="kraken_db", taxonomy_path=None):
        seen["taxonomy"] = db_path

    def fake_build(db_path="krakenuniq_db", threads=1):
        seen["build"] = db_path
        seen["threads"] = threads

    with patch("samovar.build_database.get_taxonomy_db", side_effect=fake_tax), patch(
        "samovar.build_database.build_database_krakenunique", side_effect=fake_build
    ):
        build_database_from_config(str(cfg), db_type="krakenunique", db_path=str(db), threads=3)

    added = db / "library" / "added"
    genome = (added / "562.fa").read_text(encoding="utf-8")
    assert genome.startswith(">")
    assert "ATGC" in genome
    assert not (added / "562.fa").read_bytes().startswith(b"\x1f\x8b")
    human = (added / "9606.fa").read_text(encoding="utf-8")
    assert "AAAAAAAAAAAA" in human
    assert not (added / "4932.fa").exists()
    mapped = (added / "562.map").read_text(encoding="utf-8").strip().split("\t")
    assert mapped[1] == "562"
    assert len(mapped) == 3
    missing = (db / "missing-taxa.tsv").read_text(encoding="utf-8")
    assert "4932\tproteome\t" in missing
    assert not (db / "genomes.map").exists()
    assert not (db / "562.fna").exists()
    assert not (db / "library-files.txt").exists()
    assert seen["taxonomy"].endswith(os.path.join("taxonomy")) or seen["taxonomy"].endswith("taxonomy")
    assert seen["build"] == str(db)
    assert seen["threads"] == 3


def test_krakenunique_skip_does_not_abort_when_one_genome_is_missing(tmp_path):
    db = tmp_path / "db"
    present = tmp_path / "562.fna"
    _nucleotide(present)
    add_database_krakenunique(str(present), "562", "562.fna", str(db))
    try:
        add_database_krakenunique(str(tmp_path / "gone.fna"), "999", "gone.fna", str(db))
    except KrakenUniqLibrarySkip as exc:
        assert exc.reason == "missing"
        assert exc.taxid == "999"
    else:
        raise AssertionError("missing genome should be reported, not indexed")
    assert (db / "library" / "added" / "562.fa").is_file()
    assert not (db / "library" / "added" / "999.fa").exists()


def test_krakenunique_build_drops_stale_library_list(tmp_path, monkeypatch):
    db = tmp_path / "db"
    added = db / "library" / "added"
    added.mkdir(parents=True)
    (added / "562.fa").write_text(">a\nATGC\n", encoding="utf-8")
    tax = db / "taxonomy"
    tax.mkdir()
    (tax / "nodes.dmp").write_text("1\t|\t1\t|\n", encoding="utf-8")
    (tax / "names.dmp").write_text("1\t|\troot\t|\n", encoding="utf-8")
    stale = db / "library-files.txt"
    stale.write_text("", encoding="utf-8")
    fake_bin = tmp_path / "krakenuniq-build"
    fake_bin.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    fake_bin.chmod(fake_bin.stat().st_mode | stat.S_IEXEC)
    jelly = tmp_path / "jellyfish"
    jelly.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    jelly.chmod(jelly.stat().st_mode | stat.S_IEXEC)

    def resolve(name, tool_key=None):
        token = str(name)
        if "jellyfish" in token:
            return str(jelly)
        if "krakenuniq-build" in token:
            return str(fake_bin)
        return token

    recorded = {}

    def fake_run(cmd, check=True, text=True, env=None):
        recorded["cmd"] = list(cmd)
        class Result:
            returncode = 0
        return Result()

    monkeypatch.setattr("samovar.paths.resolve_executable", resolve)
    monkeypatch.setattr("samovar.build_database.subprocess.run", fake_run)
    build_database_krakenunique(str(db), threads=2)
    assert not stale.exists()
    assert recorded["cmd"][0] == str(fake_bin)
    assert "--db" in recorded["cmd"]
    assert str(db) in recorded["cmd"]
    assert "--taxids-for-genomes" in recorded["cmd"]
