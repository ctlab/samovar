"""Tests for NCBI genome cache reuse (never silent-fallback to truncated test_genomes)."""

from __future__ import annotations

from pathlib import Path

from samovar.genome_cache import (
    dest_name_for_taxid,
    find_library_genome,
    is_bundled_test_genomes_path,
    place_genome,
    register_genome_dir,
    seed_run_genomes,
)
from samovar.paths import test_genomes_dir as bundled_test_genomes_dir


def _isolate_config(tmp_path, monkeypatch):
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    cache = tmp_path / "cache"
    cache.mkdir()
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg / "config.json"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache))
    monkeypatch.delenv("SAMOVAR_GENOMES", raising=False)
    monkeypatch.delenv("SAMOVAR_PROCESSED_GENOMES", raising=False)
    monkeypatch.delenv("SAMOVAR_GENOME_DIRS", raising=False)
    monkeypatch.delenv("SAMOVAR_ALLOW_TEST_GENOMES", raising=False)
    monkeypatch.setenv("SAMOVAR_REUSE_GENOMES", "1")
    (cfg / "config.json").write_text("{}\n")
    return cfg / "config.json"


def test_refuses_bundled_test_genomes_as_library():
    tg = bundled_test_genomes_dir()
    assert is_bundled_test_genomes_path(tg)
    assert is_bundled_test_genomes_path(tg / "meta")
    assert is_bundled_test_genomes_path(tg / "host" / "9606.fna")
    assert not is_bundled_test_genomes_path(tg.parent)


def test_register_genome_dir_rejects_test_genomes(tmp_path, monkeypatch, caplog):
    _isolate_config(tmp_path, monkeypatch)
    assert register_genome_dir(bundled_test_genomes_dir() / "meta") is None
    lib = tmp_path / "ncbi_lib"
    lib.mkdir()
    (lib / "562.fna").write_text(">x\nACGT\n")
    registered = register_genome_dir(lib)
    assert registered is not None
    from samovar.paths import load_config

    dirs = load_config().get("genome_dirs") or []
    assert any(str(lib.resolve()) in str(d) for d in dirs)


def test_place_genome_prefers_symlink(tmp_path):
    src = tmp_path / "src" / "562.fna"
    src.parent.mkdir()
    src.write_text(">x\nACGT\n")
    dest = tmp_path / "dest"
    placed = place_genome(src, dest, "562")
    assert placed.is_symlink()
    assert placed.read_text() == ">x\nACGT\n"
    assert placed.name == "562.fna"


def test_dest_name_strips_library_prefix():
    assert dest_name_for_taxid("562", "Bacteria_562-processed.fasta") == "562-processed.fasta"
    assert dest_name_for_taxid("562", "562.fna.gz") == "562.fna.gz"


def test_find_library_skips_test_genomes(tmp_path, monkeypatch):
    _isolate_config(tmp_path, monkeypatch)
    monkeypatch.setenv("SAMOVAR_GENOME_DIRS", "")
    assert find_library_genome("562", extra=[], include_test_genomes=False) is None
    lib = tmp_path / "lib"
    lib.mkdir()
    (lib / "562-processed.fasta").write_text(">x\nACGT\n")
    found = find_library_genome("562", extra=[lib], include_test_genomes=False)
    assert found is not None
    assert found.name == "562-processed.fasta"


def test_seed_from_generate_sources_not_registered(tmp_path, monkeypatch):
    _isolate_config(tmp_path, monkeypatch)
    run = tmp_path / "run"
    gen_cfg = run / ".generate" / "configs"
    gen_cfg.mkdir(parents=True)
    meta = bundled_test_genomes_dir() / "meta"
    host = bundled_test_genomes_dir() / "host" / "9606.fna"
    (gen_cfg / "iss_config.yaml").write_text(
        f"genome_dir: {meta}\nhost_genome: {host}\n",
        encoding="utf-8",
    )
    dest = run / "genomes"
    stats = seed_run_genomes(dest, reuse=True, generate_output_dir=run)
    assert stats["linked"] or stats["copied"]
    assert any(p.name.startswith("562") for p in dest.iterdir())
    from samovar.paths import load_config

    dirs = [str(Path(d).resolve()) for d in (load_config().get("genome_dirs") or [])]
    assert not any("test_genomes" in d for d in dirs)


def test_seed_from_camisim_generate_yaml(tmp_path, monkeypatch):
    _isolate_config(tmp_path, monkeypatch)
    run = tmp_path / "run"
    gen_cfg = run / ".generate" / "configs"
    gen_cfg.mkdir(parents=True)
    meta = bundled_test_genomes_dir() / "meta"
    host = bundled_test_genomes_dir() / "host" / "9606.fna"
    (gen_cfg / "camisim.yaml").write_text(
        f"simulator: camisim\ngenome_dir: {meta}\nhost_genome: {host}\n",
        encoding="utf-8",
    )
    dest = run / "genomes"
    stats = seed_run_genomes(dest, reuse=True, generate_output_dir=run)
    assert stats["linked"] or stats["copied"]
    assert any(p.name.startswith("562") for p in dest.iterdir())
    assert any(p.name.startswith("9606") for p in dest.iterdir())


def test_seed_reuses_prefixed_processed(tmp_path, monkeypatch):
    _isolate_config(tmp_path, monkeypatch)
    src = tmp_path / "realistic"
    src.mkdir()
    (src / "Bacteria_1664297-processed.fasta").write_text(">x\nACGT\n")
    dest = tmp_path / "out_genomes"
    stats = seed_run_genomes(dest, reuse=True, extra_dirs=[src])
    assert stats["linked"] or stats["copied"]
    placed = dest / "1664297-processed.fasta"
    assert placed.exists()
    assert placed.is_symlink()
