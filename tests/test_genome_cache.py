"""Bundled test_genomes are never an NCBI library. Index/reuse lives in test_reindex.py."""

from __future__ import annotations

from samovar.genome_cache import (
    is_bundled_test_genomes_path,
    run_genome_view,
    seed_run_genomes,
)
from samovar.paths import test_genomes_dir as bundled_test_genomes_dir


def test_refuses_bundled_test_genomes_as_library():
    tg = bundled_test_genomes_dir()
    assert is_bundled_test_genomes_path(tg)
    assert is_bundled_test_genomes_path(tg / "meta")
    assert is_bundled_test_genomes_path(tg / "host" / "9606.fna")
    assert not is_bundled_test_genomes_path(tg.parent)


def test_seed_creates_the_hidden_view_only_when_a_genome_is_linked(tmp_path):
    empty = tmp_path / "empty-run"
    stats = seed_run_genomes(empty / ".genomes", generate_output_dir=empty, reuse=False)
    assert stats["linked"] == 0 and stats["copied"] == 0
    assert not (empty / ".genomes").exists()
    assert not (empty / "genomes").exists()

    source = tmp_path / "library"
    source.mkdir()
    (source / "562.fna").write_text(">562\nACGT\n", encoding="utf-8")
    (source / "562.fna.sample").write_text("no\n", encoding="utf-8")
    generate = tmp_path / "generated"
    cfg = generate / ".generate" / "configs"
    cfg.mkdir(parents=True)
    (cfg / "iss_config.yaml").write_text(
        "genome_dir: {path}\n".format(path=source), encoding="utf-8"
    )
    dest = generate / ".genomes"
    stats = seed_run_genomes(dest, generate_output_dir=generate, reuse=False)
    linked = dest / "562.fna"
    assert stats["linked"] == 1 and stats["copied"] == 0
    assert linked.is_symlink()
    assert linked.resolve() == (source / "562.fna").resolve()


def test_run_genome_view_keeps_an_old_directory_and_does_not_create_one(tmp_path):
    fresh = tmp_path / "fresh"
    assert run_genome_view(fresh) == fresh / ".genomes"
    assert not (fresh / ".genomes").exists()
    assert not (fresh / "genomes").exists()

    old = tmp_path / "old"
    (old / "genomes").mkdir(parents=True)
    (old / "genomes" / "562.fna").write_text(">562\nAC\n", encoding="utf-8")
    assert run_genome_view(old) == old / "genomes"

    blank = tmp_path / "blank"
    (blank / "genomes").mkdir(parents=True)
    assert run_genome_view(blank) == blank / ".genomes"
