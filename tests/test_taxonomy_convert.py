from __future__ import annotations

from pathlib import Path

import pandas as pd

from samovar.annotation_convert import convert_annotation, main as convert_main
from samovar.kreport import KReportTaxonomy, counts_to_kreport
from samovar.parse_annotators import Annotation
from samovar.taxonomy import normalize_taxonomy
from tests.test_cami_profile import assert_rfc, parse_profile


def write_mini_gtdb_dmp(root: Path) -> Path:
    """GTDB-like dump: domain Bacteria → Escherichia → E. coli."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "2\t|\t1\t|\tdomain\t|\n"
        "561\t|\t2\t|\tgenus\t|\n"
        "562\t|\t561\t|\tspecies\t|\n",
        encoding="utf-8",
    )
    (root / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
        "561\t|\tEscherichia\t|\t\t|\tscientific name\t|\n"
        "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    return root


def write_mini_gtdb_tsv(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    table = root / "bac120_taxonomy.tsv"
    table.write_text(
        "RS_GCF_000005845.2\t"
        "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
        "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;"
        "s__Escherichia coli\n",
        encoding="utf-8",
    )
    return root


def test_normalize_taxonomy():
    assert normalize_taxonomy("") == "ncbi"
    assert normalize_taxonomy("NCBI") == "ncbi"
    assert normalize_taxonomy("gtdb-tk") == "gtdb"
    try:
        normalize_taxonomy("silva")
        assert False
    except ValueError:
        pass


def test_gtdb_dmp_cami_taxonomy_id(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_gtdb_dmp(tmp_path / "gtdb")
    obj = Annotation.from_long_table(
        pd.DataFrame(
            {
                "seq": ["a", "b"],
                "taxID_kaiju_0": ["562", "562"],
                "sample": ["1", "1"],
            }
        )
    )
    reports = obj.to_cami(taxdump=dump, taxonomy="gtdb")
    text = reports["kaiju"]["1"]
    assert_rfc(text)
    header = parse_profile(text)[0][0]
    assert header["taxonomyid"] == "gtdb"
    assert header["ranks"].startswith("domain|")
    by_id = {r["TAXID"]: r for r in parse_profile(text)[0][2]}
    assert by_id["2"]["RANK"] == "domain"
    assert by_id["562"]["TAXPATH"] == "2|||||561|562"


def test_gtdb_tsv_resolves_names(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_gtdb_tsv(tmp_path / "gtdb")
    tax = KReportTaxonomy.from_taxdump(dump, taxonomy="gtdb")
    assert tax.system == "gtdb"
    species = tax.resolve("s__Escherichia coli")
    assert species > 1
    assert tax.rank_of(species) == "species"
    assert tax.resolve("RS_GCF_000005845.2") == species
    text = counts_to_kreport({species: 3}, tax)
    assert "Escherichia coli" in text
    assert "\tD\t" in text


def test_convert_cli_taxonomy_gtdb(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_gtdb_dmp(tmp_path / "gtdb")
    table = tmp_path / "kaiju.csv"
    pd.DataFrame({"taxid": [562, 0], "N_1": [4, 1]}).to_csv(table, index=False)
    dest = tmp_path / "out.profile"
    rc = convert_main(
        [
            "-i",
            str(table),
            "-o",
            str(dest),
            "--from",
            "abundance",
            "--to",
            "cami",
            "--taxdump",
            str(dump),
            "--taxonomy",
            "gtdb",
        ]
    )
    assert rc == 0
    text = dest.read_text()
    assert_rfc(text)
    assert "@TaxonomyID:gtdb" in text
    assert "@Ranks:domain|" in text


def test_convert_default_is_ncbi(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    from tests.test_kreport import write_mini_taxdump

    dump = write_mini_taxdump(tmp_path / "ncbi")
    table = tmp_path / "kaiju.csv"
    pd.DataFrame({"taxid": [562], "N_1": [2]}).to_csv(table, index=False)
    dest = tmp_path / "out.profile"
    convert_annotation(
        table, dest, source_format="abundance", dest_format="cami", taxdump=str(dump)
    )
    text = Path(dest).read_text()
    assert "@TaxonomyID:ncbi" in text
    assert "@Ranks:superkingdom|" in text


def test_annotation_to_gtdb(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_gtdb_dmp(tmp_path / "gtdb")
    obj = Annotation.from_abundance_tables(
        {"kaiju": pd.DataFrame({"taxid": [562], "N_1": [2]})}
    )
    dest = tmp_path / "gtdb.profile"
    written = obj.to_GTDB(dest, taxdump=dump)
    text = Path(written).read_text()
    assert "@TaxonomyID:gtdb" in text


def test_gtdb_metadata_maps_ncbi_taxid(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = tmp_path / "gtdb"
    dump.mkdir()
    (dump / "taxonomy.tsv").write_text(
        "accession\tgtdb_taxonomy\tncbi_taxid\n"
        "GCF_000005845.2\t"
        "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
        "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;"
        "s__Escherichia coli\t562\n",
        encoding="utf-8",
    )
    tax = KReportTaxonomy.from_taxdump(dump, taxonomy="gtdb")
    leaf = tax.resolve("562")
    assert leaf == tax.resolve("s__Escherichia coli")
    assert tax.rank_of(leaf) == "species"
