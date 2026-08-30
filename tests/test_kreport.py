from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from samovar.kreport import (
    KReportTaxonomy,
    counts_to_kreport,
    counts_to_mpa,
    dump_kreport,
    rollup_clade_counts,
)
from samovar.parse_annotators import Annotation


def write_mini_taxdump(root: Path) -> Path:
    """Tiny NCBI-like dump: root → Bacteria → Escherichia → E. coli → strain."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "2\t|\t1\t|\tsuperkingdom\t|\n"
        "561\t|\t2\t|\tgenus\t|\n"
        "562\t|\t561\t|\tspecies\t|\n"
        "83333\t|\t562\t|\tno rank\t|\n",
        encoding="utf-8",
    )
    (root / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
        "561\t|\tEscherichia\t|\t\t|\tscientific name\t|\n"
        "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n"
        "83333\t|\tEscherichia coli K-12\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    return root


def _parse_kreport(text: str):
    rows = []
    for line in text.strip().splitlines():
        pct, clade, direct, rank, taxid, name = line.split("\t")
        rows.append(
            {
                "pct": float(pct),
                "clade": int(clade),
                "direct": int(direct),
                "rank": rank,
                "taxid": int(taxid),
                "name": name,
            }
        )
    return rows


def test_rank_codes_and_clade_rollup(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    tax = KReportTaxonomy.from_taxdump(dump)
    assert tax.rank_code(0) == "U"
    assert tax.rank_code(1) == "R"
    assert tax.rank_code(2) == "D"
    assert tax.rank_code(561) == "G"
    assert tax.rank_code(562) == "S"
    assert tax.rank_code(83333) == "S1"
    direct = {562: 2, 0: 1}
    clade = rollup_clade_counts(direct, tax)
    assert clade[562] == 2
    assert clade[561] == 2
    assert clade[2] == 2
    assert clade[1] == 2
    assert 0 not in clade


def test_counts_to_kreport_layout(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    tax = KReportTaxonomy.from_taxdump(write_mini_taxdump(tmp_path / "taxdump"))
    text = counts_to_kreport({0: 1, 562: 2}, tax)
    rows = _parse_kreport(text)
    by_id = {row["taxid"]: row for row in rows}
    assert rows[0]["taxid"] == 0
    assert rows[0]["rank"] == "U"
    assert rows[0]["clade"] == 1 and rows[0]["direct"] == 1
    assert by_id[1]["rank"] == "R" and by_id[1]["direct"] == 0 and by_id[1]["clade"] == 2
    assert by_id[2]["rank"] == "D" and by_id[2]["name"].endswith("Bacteria")
    assert by_id[2]["name"].startswith("  ")
    assert by_id[562]["rank"] == "S" and by_id[562]["direct"] == 2
    assert by_id[562]["name"].strip() == "Escherichia coli"
    assert abs(by_id[0]["pct"] - 33.33) < 0.01
    assert abs(by_id[562]["pct"] - 66.67) < 0.01


def test_mpa_style_lineage_and_clade_counts(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    tax = KReportTaxonomy.from_taxdump(write_mini_taxdump(tmp_path / "taxdump"))
    text = counts_to_mpa({562: 2, 0: 1}, tax)
    lines = dict(line.split("\t") for line in text.strip().splitlines())
    assert lines["d__Bacteria"] == "2"
    assert lines["d__Bacteria|g__Escherichia"] == "2"
    assert lines["d__Bacteria|g__Escherichia|s__Escherichia_coli"] == "2"
    assert "unclassified" not in text
    assert "root" not in text


def test_annotation_to_kraken2_reuses_abundance(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    obj = Annotation.from_long_table(
        pd.DataFrame(
            {
                "seq": ["a", "b", "c"],
                "taxID_kaiju_0": ["562", "562", "0"],
                "sample": ["1", "1", "1"],
            }
        )
    )
    tables = obj.to_abundance()
    n562 = int(tables["kaiju"][tables["kaiju"]["taxid"].astype(str) == "562"]["N_1"].iloc[0])
    reports = obj.to_kraken2(taxdump=dump)
    rows = _parse_kreport(reports["kaiju"]["1"])
    eco = next(row for row in rows if row["taxid"] == 562)
    assert eco["direct"] == n562 == 2
    dest = tmp_path / "out.kreport"
    written = obj.to_kraken2(dest, taxdump=dump)
    assert Path(written).read_text().startswith(" ")
    assert "Escherichia coli" in Path(written).read_text()


def test_missing_taxdump_message(monkeypatch):
    monkeypatch.delenv("SAMOVAR_TAXDUMP", raising=False)
    monkeypatch.delenv("SAMOVAR_NODES_DMP", raising=False)
    monkeypatch.setattr("samovar.kreport.nodes_dmp", lambda cfg=None: None)
    monkeypatch.setattr("samovar.kreport.names_dmp", lambda cfg=None: None)
    with pytest.raises(FileNotFoundError, match="taxdump"):
        KReportTaxonomy.from_taxdump()


def test_dump_kreport_directory_layout(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    tables = {
        "kaiju": pd.DataFrame({"taxid": [562, 0], "N_1": [2, 1], "N_2": [1, 0]}),
        "kraken2": pd.DataFrame({"taxid": [562], "N_1": [3]}),
    }
    dest = dump_kreport(tables, tmp_path / "reports", taxdump=dump)
    assert (dest / "kaiju" / "1.kreport").is_file()
    assert (dest / "kaiju" / "2.kreport").is_file()
    assert (dest / "kraken2" / "1.kreport").is_file()
    mpa_dest = dump_kreport(tables, tmp_path / "mpa", taxdump=dump, mpa=True)
    assert (mpa_dest / "kaiju" / "1.mpa").is_file()
    assert "d__Bacteria" in (mpa_dest / "kaiju" / "1.mpa").read_text()
