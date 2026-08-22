"""Toy-only Kraken2 / Kaiju taxon gaps (example, not production)."""

from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from samovar.build_database import build_database_from_config
from samovar.example_db import (
    ESCHERICHIA_TAXIDS,
    EXAMPLE_OMIT_WARNING,
    PHIX_TAXIDS,
    filter_example_omit,
)


def test_filter_example_omit_kraken2_drops_escherichia():
    mapping = {
        "562.fna": "562",
        "4932.fna": "4932",
        "2886930.fna": "2886930",
        "9606.fna": "9606",
    }
    kept, omitted = filter_example_omit(mapping, "kraken2")
    assert "562" in omitted.values()
    assert set(omitted.values()) <= ESCHERICHIA_TAXIDS
    assert "4932" in kept.values()
    assert "2886930" in kept.values()


def test_filter_example_omit_kaiju_drops_phix():
    mapping = {
        "562.fna": "562",
        "4932.fna": "4932",
        "2886930.fna": "2886930",
        "9606.fna": "9606",
    }
    kept, omitted = filter_example_omit(mapping, "kaiju")
    assert "2886930" in omitted.values()
    assert set(omitted.values()) <= PHIX_TAXIDS
    assert "562" in kept.values()


def _tiny_genome_dir(tmp_path: Path) -> Path:
    genomes = tmp_path / "test_genomes" / "meta"
    genomes.mkdir(parents=True)
    (genomes / "562.fna").write_text(">e\n" + "ATGGAATTCGGT" * 8 + "TAA\n")
    (genomes / "4932.fna").write_text(">y\n" + "ATGGAATTCGGT" * 8 + "TAA\n")
    (genomes / "2886930.fna").write_text(">p\n" + "ATGGAATTCGGT" * 8 + "TAA\n")
    return genomes


def test_build_kraken2_example_omit_skips_escherichia(tmp_path):
    genomes = _tiny_genome_dir(tmp_path)
    cfg = tmp_path / "db.yaml"
    cfg.write_text(yaml.dump({"input_dir": [str(genomes)], "output_dir": str(tmp_path / "prep")}))
    db_path = tmp_path / "kraken2_db"
    added = []

    def fake_add(input_file, taxid, db_path="kraken_db"):
        added.append(str(taxid))

    with patch("samovar.build_database.add_database_kraken2", side_effect=fake_add), patch(
        "samovar.build_database.get_taxonomy_db"
    ), patch("samovar.build_database.build_database_kraken2"), pytest.warns(UserWarning, match="EXAMPLE"):
        build_database_from_config(
            str(cfg), db_type="kraken2", db_path=str(db_path), example_omit=True
        )
    assert "562" not in added
    assert "4932" in added
    assert "2886930" in added
    assert "Escherichia" in EXAMPLE_OMIT_WARNING


def test_build_kaiju_example_omit_skips_phix(tmp_path):
    genomes = _tiny_genome_dir(tmp_path)
    cfg = tmp_path / "db.yaml"
    cfg.write_text(yaml.dump({"input_dir": [str(genomes)], "output_dir": str(tmp_path / "prep")}))
    db_path = tmp_path / "kaiju_db"
    added = []

    def fake_add(input_file, taxid, db_path="kaiju_db", protein=False, fetch_missing=True):
        added.append(str(taxid))
        Path(db_path).mkdir(parents=True, exist_ok=True)
        lib = Path(db_path) / "library.faa"
        lib.write_text(f">{taxid}\nMEEF\n")

    with patch("samovar.build_database.add_database_kaiju", side_effect=fake_add), patch(
        "samovar.build_database.get_taxonomy_db"
    ), patch("samovar.build_database.build_database_kaiju"), pytest.warns(UserWarning, match="Phi X"):
        build_database_from_config(
            str(cfg), db_type="kaiju", db_path=str(db_path), example_omit=True
        )
    assert "2886930" not in added
    assert "562" in added
    assert "4932" in added


def test_prepare_custom_kraken2_kaiju_db_paths(tmp_path):
    from samovar.config import setup_pipeline

    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    k2 = tmp_path / "kraken2_db"
    kj = tmp_path / "kaiju_db"
    k2.mkdir()
    kj.mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(reads),
            "output_dir": str(tmp_path / "out"),
            "kraken2": [[f"kraken2 {k2}"]],
            "kaiju": [[f"kaiju {kj}"]],
        },
    )()
    result = setup_pipeline(args)
    init = Path(result["configs"]["init_annotator"]).read_text()
    assert "kraken2" in init and "kaiju" in init
    assert str(k2.resolve()) in init
    assert str(kj.resolve()) in init
