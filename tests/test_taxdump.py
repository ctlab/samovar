from pathlib import Path

from samovar.taxdump import (
    TAXDUMP_DMP_FILES,
    find_dmp,
    link_taxdump_into,
    load_scientific_names,
    parse_names_dmp,
    taxdump_dir,
)
from samovar.build_database import get_taxonomy_db
from samovar.main_config import build_install_config, disk_payload
from samovar.version import get_version


def test_link_taxdump_into_symlinks_dmp_only(tmp_path):
    src = tmp_path / "taxdump"
    src.mkdir()
    (src / "nodes.dmp").write_text("1\t|\t1\t|\n")
    (src / "names.dmp").write_text("1\t|\troot\t|\n")
    (src / "prot.accession2taxid.gz").write_text("huge")
    dest = tmp_path / "kaiju_db"
    link_taxdump_into(dest, src)
    assert (dest / "nodes.dmp").is_symlink()
    assert (dest / "nodes.dmp").resolve() == (src / "nodes.dmp").resolve()
    assert (dest / "names.dmp").is_symlink()
    assert not (dest / "prot.accession2taxid.gz").exists()


def test_get_taxonomy_db_uses_source_dir(tmp_path, monkeypatch):
    src = tmp_path / "shared"
    src.mkdir()
    (src / "nodes.dmp").write_text("1\t|\t1\t|\n")
    (src / "names.dmp").write_text("1\t|\troot\t|\n")
    db = tmp_path / "db"
    get_taxonomy_db(db_path=str(db), taxonomy_path=str(src))
    assert (db / "nodes.dmp").is_symlink()


def test_build_install_config_records_taxdump(tmp_path):
    dump = tmp_path / "mytax"
    dump.mkdir()
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(tmp_path / "python"),
        version=get_version(),
        samovar_database=str(tmp_path / "genomes"),
        taxdump=str(dump),
    )
    payload = disk_payload(cfg)
    assert payload["genomes"]["taxdump"] == str(dump)


def test_taxdump_dir_defaults_under_samovar_database(tmp_path, monkeypatch):
    monkeypatch.delenv("SAMOVAR_TAXDUMP", raising=False)
    store = tmp_path / "genomes"
    cfg = {
        "genomes": {"samovar_database": str(store), "taxdump": ""},
    }
    assert taxdump_dir(cfg) == store / "taxdump"


def test_taxdump_dir_falls_back_to_indexed_database(tmp_path, monkeypatch):
    monkeypatch.delenv("SAMOVAR_TAXDUMP", raising=False)
    dump = tmp_path / "ncbi_taxdump"
    dump.mkdir()
    (dump / "nodes.dmp").write_text("1\t|\n")
    cfg = {
        "genomes": {"samovar_database": str(tmp_path / "store"), "taxdump": ""},
        "databases": {
            "taxdump": {
                "ncbi": {"path": str(dump), "type": "database"},
            }
        },
    }
    assert taxdump_dir(cfg) == dump


def test_find_dmp_nested_taxonomy(tmp_path):
    nested = tmp_path / "taxonomy"
    nested.mkdir()
    (nested / "nodes.dmp").write_text("1\t|\n")
    assert find_dmp("nodes.dmp", tmp_path) == nested / "nodes.dmp"
    assert "nodes.dmp" in TAXDUMP_DMP_FILES


def test_parse_and_cache_scientific_names(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    names = tmp_path / "names.dmp"
    names.write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "1\t|\tall\t|\t\t|\tsynonym\t|\n"
        "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    parsed = parse_names_dmp(names)
    assert parsed[1] == "root"
    assert parsed[562] == "Escherichia coli"
    assert 1 in load_scientific_names(names)
    assert load_scientific_names(names)[562] == "Escherichia coli"
