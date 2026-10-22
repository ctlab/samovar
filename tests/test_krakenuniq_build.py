"""KrakenUniq library layout: uncompressed FASTA in library/added, dumps in taxonomy/."""

import gzip
import os
import stat
from pathlib import Path
from unittest.mock import patch

import yaml
import pytest

from samovar.build_database import (
    KrakenUniqLibrarySkip,
    add_database_krakenunique,
    build_database_from_config,
    build_database_krakenunique,
)


def _nucleotide(path: Path, seq: str = "ATGCATGCATGC") -> None:
    path.write_text(f">seq\n{seq}\n", encoding="utf-8")


@pytest.mark.parametrize("minimizer_len", [None, 11])
def test_krakenunique_writes_added_fasta_and_skips_proteome(tmp_path, minimizer_len):
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

    def fake_build(db_path="krakenuniq_db", threads=1, minimizer_len=None):
        seen["build"] = db_path
        seen["threads"] = threads
        seen["minimizer_len"] = minimizer_len

    with patch("samovar.build_database.get_taxonomy_db", side_effect=fake_tax), patch(
        "samovar.build_database.build_database_krakenunique", side_effect=fake_build
    ):
        build_database_from_config(str(cfg), db_type="krakenunique", db_path=str(db), threads=3, minimizer_len=minimizer_len)

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
    assert seen["minimizer_len"] == minimizer_len


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


@pytest.mark.parametrize("library_bytes,override,expected", [
    (16, None, 12),
    (512 * 1024**2, None, 13),
    (2 * 1024**3, None, 14),
    (8 * 1024**3, None, 15),
    (16, 10, 10),
    (16, 15, 15),
    (16, 0, None),
    (16, 31, None),
])
def test_krakenunique_build_drops_stale_list_and_sizes_index(
    tmp_path, monkeypatch, library_bytes, override, expected,
):
    db = tmp_path / "db"
    added = db / "library" / "added"
    added.mkdir(parents=True)
    (added / "562.fa").write_text(">a\nATGC\n", encoding="utf-8")
    with (added / "562.fa").open("r+b") as handle:
        handle.truncate(library_bytes)  # Sparse files; no large fixture allocation.
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
    if expected is None:
        with pytest.raises(ValueError, match="between 1 and 30"):
            build_database_krakenunique(str(db), threads=2, minimizer_len=override)
        assert "cmd" not in recorded
        assert stale.exists()
        return
    build_database_krakenunique(str(db), threads=2, minimizer_len=override)
    assert not stale.exists()
    assert recorded["cmd"][0] == str(fake_bin)
    assert "--db" in recorded["cmd"]
    assert "--rebuild" in recorded["cmd"]
    assert str(db) in recorded["cmd"]
    assert "--taxids-for-genomes" not in recorded["cmd"]
    assert "--taxids-for-sequences" not in recorded["cmd"]
    args = recorded["cmd"]
    assert args[args.index("--minimizer-len") + 1] == str(expected)


@pytest.mark.parametrize('db_type', ['kraken', 'krakenunique'])
def test_build_cli_forwards_minimizer_without_changing_annotation_flags(monkeypatch, db_type):
    import runpy
    import sys

    entrypoint = Path(__file__).resolve().parents[1] / 'workflow/database_prep/samovar_build_database.py'
    monkeypatch.setattr(sys, 'argv', [str(entrypoint), '--type', db_type,
                                    '--config_path', 'genomes.yaml', '--db_path', '/tmp/db',
                                    '--minimizer-len', '11', '--threads', '2', '--index', 'test',
                                    '--flags=--preload'])
    with patch('samovar.build_database.build_database_from_config') as build, patch(
        'samovar.genome_index.register_database'
    ) as register:
        runpy.run_path(str(entrypoint), run_name='__main__')
    assert build.call_args.kwargs['minimizer_len'] == 11
    assert build.call_args.kwargs['threads'] == 2
    assert build.call_args.kwargs['db_type'] == 'krakenunique'
    register.assert_called_once_with('krakenunique', 'test', '/tmp/db', '--preload')


@pytest.mark.parametrize('length', ['0', '31', 'not-an-int'])
def test_build_cli_rejects_invalid_minimizer(monkeypatch, length):
    import runpy
    import sys

    entrypoint = Path(__file__).resolve().parents[1] / 'workflow/database_prep/samovar_build_database.py'
    monkeypatch.setattr(sys, 'argv', [str(entrypoint), '--type', 'krakenunique',
                                    '--config_path', 'unused.yaml', '--db_path', '/tmp/db',
                                    '--minimizer-len', length])
    with patch('samovar.build_database.build_database_from_config') as build:
        with pytest.raises(SystemExit) as exc:
            runpy.run_path(str(entrypoint), run_name='__main__')
    assert exc.value.code == 2
    build.assert_not_called()
