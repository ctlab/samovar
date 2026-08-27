"""``--reannotation-level`` catalog lookup (up and down the NCBI ranks)."""

import gzip
from pathlib import Path

from samovar.genome_resolve import (
    catalog_hits_at_rank,
    fetch_taxid_for_level,
    normalize_reannotation_level,
    pick_catalog_genome,
    ranks_for_level,
)


def _fake_exact(taxid, rank):
    lineage = {
        511145: {"species": 562, "genus": 561, "family": 543},
        562: {"species": 562, "genus": 561, "family": 543},
        561: {"genus": 561, "family": 543},
        543: {"family": 543},
        10847: {"species": 10847},
        2886930: {"species": 2886930},
    }
    mapped = lineage.get(int(taxid), {}).get(rank)
    return str(mapped) if mapped else None


def _cfg(tmp_path: Path, taxid="562"):
    store = tmp_path / "processed"
    store.mkdir()
    fasta = store / "GCF_test.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">x\nACGT\n")
    return {
        "genomes": {
            "samovar_database": str(tmp_path),
            "processed": {"samovar_database": str(store)},
            "data": {
                taxid: [taxid, "GCF_test", "samovar_database", "GCF_test.fa.gz"],
            },
        }
    }


def test_normalize_reannotation_level():
    assert normalize_reannotation_level("t") == "taxid"
    assert normalize_reannotation_level("species") == "species"
    assert normalize_reannotation_level("g") == "genus"
    assert normalize_reannotation_level("f") == "family"
    assert normalize_reannotation_level("a") == "any"
    assert ranks_for_level("taxid") == ()
    assert "genus" in ranks_for_level("g")
    assert "phylum" in ranks_for_level("any")


def test_pick_exact_taxid(tmp_path, monkeypatch):
    monkeypatch.setattr("samovar.genome_resolve._resolve_taxid_by_rank_exact", _fake_exact)
    monkeypatch.setattr("samovar.genome_resolve.canonical_ncbi_taxid", lambda x, cfg=None: str(x))
    cfg = _cfg(tmp_path)
    hit = pick_catalog_genome("562", "taxid", cfg)
    assert hit is not None
    assert hit[0] == "562"
    assert hit[2] == "taxid"
    assert pick_catalog_genome("561", "taxid", cfg) is None


def test_species_up_from_strain(tmp_path, monkeypatch):
    monkeypatch.setattr("samovar.genome_resolve._resolve_taxid_by_rank_exact", _fake_exact)
    monkeypatch.setattr("samovar.genome_resolve.canonical_ncbi_taxid", lambda x, cfg=None: str(x))
    cfg = _cfg(tmp_path)
    hit = pick_catalog_genome("511145", "species", cfg)
    assert hit is not None
    assert hit[0] == "562"
    assert hit[2] == "species"


def test_genus_down_from_genus_taxid(tmp_path, monkeypatch):
    monkeypatch.setattr("samovar.genome_resolve._resolve_taxid_by_rank_exact", _fake_exact)
    monkeypatch.setattr("samovar.genome_resolve.canonical_ncbi_taxid", lambda x, cfg=None: str(x))
    cfg = _cfg(tmp_path)
    hits = catalog_hits_at_rank("561", "genus", cfg)
    assert hits and hits[0][0] == "562"
    hit = pick_catalog_genome("561", "genus", cfg)
    assert hit is not None and hit[2] == "genus"


def test_family_down(tmp_path, monkeypatch):
    monkeypatch.setattr("samovar.genome_resolve._resolve_taxid_by_rank_exact", _fake_exact)
    monkeypatch.setattr("samovar.genome_resolve.canonical_ncbi_taxid", lambda x, cfg=None: str(x))
    cfg = _cfg(tmp_path)
    hit = pick_catalog_genome("543", "family", cfg)
    assert hit is not None and hit[0] == "562"


def test_any_walks_to_genus(tmp_path, monkeypatch):
    monkeypatch.setattr("samovar.genome_resolve._resolve_taxid_by_rank_exact", _fake_exact)
    monkeypatch.setattr("samovar.genome_resolve.canonical_ncbi_taxid", lambda x, cfg=None: str(x))
    cfg = _cfg(tmp_path)
    hit = pick_catalog_genome("561", "any", cfg)
    assert hit is not None and hit[2] == "genus"


def test_fetch_taxid_for_level_progresses():
    assert fetch_taxid_for_level("562", "taxid", []) == "562"
    assert fetch_taxid_for_level("562", "taxid", ["562"]) is None


def test_resolve_genome_file_uses_catalog(tmp_path, monkeypatch):
    from samovar.genome_resolve import resolve_genome_file

    monkeypatch.setattr("samovar.genome_resolve._resolve_taxid_by_rank_exact", _fake_exact)
    monkeypatch.setattr("samovar.genome_resolve.canonical_ncbi_taxid", lambda x, cfg=None: str(x))
    cfg = _cfg(tmp_path)
    empty = tmp_path / "run"
    empty.mkdir()
    path = resolve_genome_file(
        "511145",
        str(empty),
        "t@e",
        level="species",
        cfg=cfg,
    )
    assert path and Path(path).is_file()


def test_catalog_processed_updates_config(tmp_path, monkeypatch):
    from samovar.genome_index import catalog_processed_genome, genome_data_map
    from samovar.paths import write_config

    store = tmp_path / "genomes"
    cfg_path = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg_path))
    write_config(
        {
            "root": str(tmp_path),
            "genomes": {"samovar_database": str(store), "processed": {}, "data": {}},
        },
        also_repo_build=False,
    )
    src = tmp_path / "run" / "562-processed.fasta.gz"
    src.parent.mkdir()
    with gzip.open(src, "wt") as handle:
        handle.write(">x\nACGT\n")
    placed = catalog_processed_genome(src, taxid="562", keep_src=True)
    assert placed.is_file()
    data = genome_data_map()
    assert "562" in data
    assert src.is_file()
