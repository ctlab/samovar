from pathlib import Path

from samovar.taxdump import (
    TAXDUMP_DMP_FILES,
    attach_taxdump_catalog,
    ensure_taxdump,
    find_dmp,
    link_taxdump_into,
    load_scientific_names,
    metax_dmp_dir,
    parse_names_dmp,
    taxdump_dir,
    usable_taxdump_dir,
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


def test_usable_taxdump_prefers_configured_dump_over_download(tmp_path, monkeypatch):
    monkeypatch.delenv("SAMOVAR_TAXDUMP", raising=False)
    monkeypatch.delenv("SAMOVAR_TAXDUMP_VERSION", raising=False)
    real = tmp_path / "ncbi"
    real.mkdir()
    nodes = real / "nodes.dmp"
    nodes.write_text("1\t|\t1\t|\n")
    empty = tmp_path / "empty-taxdump"
    cfg = {
        "genomes": {
            "samovar_database": str(tmp_path / "store"),
            "taxdump": str(empty),
        },
        "databases": {
            "taxdump": {
                "ncbi:2024-01-01": {
                    "path": str(real),
                    "version": "2024-01-01",
                    "type": "database",
                    "lazy-download": "curl should not run",
                }
            }
        },
    }

    def boom(*_args, **_kwargs):
        raise AssertionError("taxdump download was attempted")

    monkeypatch.setattr("subprocess.run", boom)
    assert usable_taxdump_dir(cfg, version="2024-01-01") == real
    assert usable_taxdump_dir(cfg, version="1999-01-01") is None
    assert ensure_taxdump(cfg).resolve() == real.resolve()
    assert metax_dmp_dir(cfg) == str(real)


def test_disk_payload_records_taxdump_version_without_download_recipe(tmp_path):
    dump = tmp_path / "taxdump"
    dump.mkdir()
    nodes = dump / "nodes.dmp"
    nodes.write_text("1\t|\t1\t|\n")
    (dump / "taxdump.version").write_text("2024-08-14\n")
    cfg = build_install_config(
        root=str(tmp_path),
        python_path=str(tmp_path / "python"),
        version=get_version(),
        samovar_database=str(tmp_path / "genomes"),
        taxdump=str(dump),
    )
    cfg.setdefault("databases", {})["taxdump"] = {
        "ncbi": {
            "path": str(dump),
            "type": "database",
            "lazy-download": "curl -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        }
    }
    payload = disk_payload(cfg)
    assert payload["genomes"]["taxdump_version"] == "2024-08-14"
    recorded = payload["databases"]["taxdump"]["ncbi:2024-08-14"]
    assert recorded["path"] == str(dump)
    assert recorded["version"] == "2024-08-14"
    assert recorded["lazy-download"] == ""
    assert "ncbi" not in payload["databases"]["taxdump"]
    attached = attach_taxdump_catalog({"genomes": {"taxdump": str(dump)}, "databases": {}})
    assert attached["genomes"]["taxdump_version"] == "2024-08-14"


def test_attach_keeps_distinct_taxdump_versions(tmp_path):
    older = tmp_path / "old"
    newer = tmp_path / "new"
    older.mkdir()
    newer.mkdir()
    (older / "nodes.dmp").write_text("1\t|\n")
    (newer / "nodes.dmp").write_text("1\t|\n")
    (older / "taxdump.version").write_text("2020-01-01\n")
    (newer / "taxdump.version").write_text("2024-01-01\n")
    payload = attach_taxdump_catalog(
        {
            "genomes": {"taxdump": str(older)},
            "databases": {
                "taxdump": {
                    "ncbi:2024-01-01": {
                        "path": str(newer),
                        "version": "2024-01-01",
                        "type": "database",
                        "lazy-download": "curl -L https://example.invalid/taxdump.tar.gz",
                    }
                }
            },
        }
    )
    recorded = payload["databases"]["taxdump"]
    assert recorded["ncbi:2020-01-01"]["path"] == str(older)
    assert recorded["ncbi:2024-01-01"]["path"] == str(newer)
    assert recorded["ncbi:2024-01-01"]["lazy-download"] == ""
    assert payload["genomes"]["taxdump_version"] == "2020-01-01"


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
