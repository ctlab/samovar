"""Phage OTU + scoring example (bootstrap + bray_ks, not direct / shannon_ks)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd

from samovar.abundance import (
    convert_to_abundance_dir,
    n_sample_columns,
    normalize_abundance_table,
)
from samovar.regenerate import regenerate_annotation_tables, stage_regenerate_tables
from samovar.table2iss import process_annotation_tables
from samovar.table_scorers import canonicalize_table_scorer, score_generated_tables, stage_score_regenerated_tables

PHAGE_TAXIDS = {"562", "2886930", "9606"}
ANN = Path("data/test_annotations")
GENOMES = Path("data/test_genomes")


def _phage_otu():
    return pd.DataFrame(
        {
            "taxid": ["562", "2886930", "9606"],
            "N_1": [40, 12, 8],
            "N_2": [35, 18, 5],
            "N_3": [50, 10, 6],
            "N_4": [28, 20, 9],
        }
    )


def test_annotation_dir_converts_without_taxid_star_in_abundance(tmp_path):
    dest = tmp_path / "abund"
    tables = convert_to_abundance_dir(ANN, dest, {"regeneration_mode": "direct"})
    assert tables
    for table in tables.values():
        out = normalize_abundance_table(table)
        assert list(out.columns)[0] == "taxid"
        assert n_sample_columns(out)
        assert not any(str(c).startswith("taxID_") for c in out.columns)


def test_phage_otu_bootstrap_bray_ks_not_direct(tmp_path):
    src = tmp_path / "initial_abundance"
    src.mkdir()
    _phage_otu().to_csv(src / "phage.csv", index=False)
    gen = tmp_path / "regenerated"
    tables = regenerate_annotation_tables(
        src,
        gen,
        {
            "regeneration_mode": "bootstrap",
            "N": 4,
            "N_reads": 80,
            "seed": 1,
            "rescale_abundance": True,
        },
        select_best=False,
    )
    assert tables
    table = next(iter(tables.values()))
    assert set(table["taxid"].astype(str)) & {"562", "2886930"}
    assert canonicalize_table_scorer("bray_ks") == "bray_ks"
    row = score_generated_tables(
        _phage_otu(),
        tables,
        scorer="bray_ks",
    )
    assert row.get("scorer") == "bray_ks"
    assert "ks_statistic" in row or "rank_value" in row


def test_scoring_example_stage_bray_ks(tmp_path):
    root = tmp_path / "run"
    (root / "initial_abundance").mkdir(parents=True)
    _phage_otu().to_csv(root / "initial_abundance" / "phage.csv", index=False)
    written = stage_regenerate_tables(
        root,
        {
            "table_reads_generator": "bootstrap",
            "table_score": "bray_ks",
            "N": 4,
            "N_reads": 60,
            "seed": 2,
        },
    )
    assert written
    payload = stage_score_regenerated_tables(
        root,
        {
            "table_score": "bray_ks",
            "table_reads_generator": "bootstrap",
        },
    )
    assert payload.get("scorer") == "bray_ks"
    sel = root / "regenerated" / ".regenerated_abundance" / "table_selection.json"
    assert sel.is_file()


def test_abundance2iss_from_phage_otu_without_annotation_csv(tmp_path):
    abund = tmp_path / "abund"
    abund.mkdir()
    _phage_otu().to_csv(abund / "phage.csv", index=False)
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for taxid in PHAGE_TAXIDS:
        src = GENOMES / "meta" / f"{taxid}.fna"
        if not src.is_file() and taxid == "9606":
            src = GENOMES / "host" / "9606.fna"
        if src.is_file():
            (genome_dir / src.name).write_text(src.read_text())
        else:
            (genome_dir / f"{taxid}.fa").write_text(f">{taxid}\n{'ATGC' * 80}\n")
    out = tmp_path / "regenerated"
    with patch("samovar.reads_generators.simulate_from_sample_tables") as mock_sim:
        mock_sim.return_value = None
        process_annotation_tables(
            table_paths=[],
            genome_dir=str(genome_dir),
            output_dir=str(out),
            email="test@samovar.com",
            abundance_dir=str(abund),
            regeneration_config={"regeneration_mode": "direct", "abundance_dir": str(abund)},
        )
    assert mock_sim.called
    sample_tables = mock_sim.call_args[0][0]
    assert sample_tables
    merged = next(iter(sample_tables.values()))
    assert "taxid" in merged.columns
    assert not any(str(c).startswith("taxID_") for c in merged.columns)
