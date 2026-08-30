from __future__ import annotations

from pathlib import Path

import pandas as pd

from samovar.abundance import n_sample_columns
from samovar.annotation_convert import convert_annotation, main as convert_main
from samovar.parse_annotators import Annotation
from samovar.paths import write_config
from samovar.tools_import import import_tool

from tests.test_kreport import write_mini_taxdump


def _ann_dir(tmp_path: Path) -> Path:
    ann = tmp_path / "ann"
    ann.mkdir()
    pd.DataFrame(
        {
            "seq": ["a", "b", "c"],
            "taxID_kaiju_0": [562, 562, 9606],
            "taxID_kraken2_0": [562, 9606, 9606],
        }
    ).to_csv(ann / "1.annotation.csv", index=False)
    return ann


def test_annotation_to_abundance_method():
    df = pd.DataFrame(
        {
            "seq": ["a", "b", "c"],
            "taxID_kaiju_0": ["562", "562", "9606"],
            "sample": ["1", "1", "1"],
        }
    )
    obj = Annotation.from_long_table(df)
    tables = obj.to_abundance()
    assert "kaiju" in tables
    kaiju = tables["kaiju"]
    assert "taxid" in kaiju.columns
    assert n_sample_columns(kaiju)
    assert int(kaiju[kaiju["taxid"].astype(str) == "562"]["N_1"].iloc[0]) == 2
    scaled = obj.to_abundance(n_reads=100)
    n1 = n_sample_columns(scaled["kaiju"])[0]
    assert int(scaled["kaiju"][n1].sum()) == 100
    long = obj.to("annotation")
    assert any(str(c).startswith("taxID") for c in long.columns)


def test_convert_cli_annotation_to_abundance(tmp_path):
    ann = _ann_dir(tmp_path)
    dest = tmp_path / "abund"
    rc = convert_main(["-i", str(ann), "-o", str(dest), "--to", "abundance"])
    assert rc == 0
    csvs = list(dest.glob("*.csv"))
    assert csvs
    kaiju = pd.read_csv(next(p for p in csvs if "kaiju" in p.name))
    assert "taxid" in kaiju.columns
    assert n_sample_columns(kaiju)


def test_convert_custom_converter(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    script = tmp_path / "toyconv.py"
    script.write_text(
        "from pathlib import Path\n"
        "from samovar.parse_annotators import Annotation\n"
        "\n"
        "def dump(annotation, dest, config):\n"
        "    dest = Path(dest)\n"
        "    dest.parent.mkdir(parents=True, exist_ok=True)\n"
        "    tables = annotation.to_abundance()\n"
        "    dest.write_text(','.join(tables) + '\\n')\n"
        "    return dest\n"
        "\n"
        "def load(path, config):\n"
        "    return Annotation.from_abundance_tables({})\n"
    )
    import_tool(
        name="toyfmt",
        tool_type="annotation-converter",
        exec_path=str(script),
        also_repo_build=False,
    )
    ann = _ann_dir(tmp_path)
    dest = tmp_path / "out.toy"
    written = convert_annotation(ann, dest, dest_format="toyfmt")
    assert Path(written).read_text().strip()
    assert "kaiju" in Path(written).read_text()


def test_convert_cli_annotation_to_kraken2(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    ann = _ann_dir(tmp_path)
    dest = tmp_path / "kreports"
    rc = convert_main(
        ["-i", str(ann), "-o", str(dest), "--to", "kraken2", "--taxdump", str(dump)]
    )
    assert rc == 0
    reports = list(dest.rglob("*.kreport"))
    assert reports
    text = reports[0].read_text()
    parts = text.strip().splitlines()[0].split("\t")
    assert len(parts) == 6
    assert any("Escherichia" in p.read_text() for p in reports)


def test_convert_to_kraken2_mpa(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    ann = _ann_dir(tmp_path)
    dest = tmp_path / "mpa"
    written = convert_annotation(
        ann, dest, dest_format="kraken2-mpa", taxdump=str(dump)
    )
    mpa_files = list(Path(written).rglob("*.mpa"))
    assert mpa_files
    assert "d__Bacteria" in mpa_files[0].read_text()
    assert "s__Escherichia_coli" in mpa_files[0].read_text()


def test_convert_cli_annotation_to_cami(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    ann = _ann_dir(tmp_path)
    dest = tmp_path / "profiles"
    rc = convert_main(
        ["-i", str(ann), "-o", str(dest), "--to", "cami", "--taxdump", str(dump)]
    )
    assert rc == 0
    from tests.test_cami_profile import assert_rfc

    profiles = list(dest.glob("*.profile"))
    assert profiles
    for path in profiles:
        text = path.read_text()
        assert_rfc(text)
        assert "@Version:0.10.0" in text
        assert "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE" in text
        assert "2|||||561|562" in text


def test_convert_abundance_to_cami(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    table = tmp_path / "kaiju.csv"
    pd.DataFrame({"taxid": [562, 0], "N_1": [4, 1]}).to_csv(table, index=False)
    dest = tmp_path / "sample.profile"
    written = convert_annotation(
        table,
        dest,
        source_format="abundance",
        dest_format="bioboxes",
        taxdump=str(dump),
    )
    text = Path(written).read_text()
    from tests.test_cami_profile import assert_rfc

    assert_rfc(text)
    assert "@SampleID:1" in text
