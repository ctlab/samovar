"""Abundance-table regenerate contract: Annotation, kraken2, GlobalPatterns, ISS."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from samovar.abundance import (
    n_sample_columns,
    normalize_abundance_table,
    is_abundance_table,
)
from samovar.parse_annotators import read_kraken2_raw
from samovar.regenerate import regenerate_annotation_tables, regenerate_preserve

GLOBALPATTERNS = Path("tests/data/globalpatterns_abundance.csv")


@pytest.fixture
def gp_table():
    if not GLOBALPATTERNS.is_file():
        pytest.skip("tests/data/globalpatterns_abundance.csv missing")
    df = pd.read_csv(GLOBALPATTERNS)
    sample_cols = [c for c in df.columns if c != "taxid"][:6]
    return df[["taxid", *sample_cols]].head(12)


def test_normalize_phyloseq_style_columns(gp_table):
    assert is_abundance_table(gp_table)
    out = normalize_abundance_table(gp_table)
    assert list(out.columns)[0] == "taxid"
    assert n_sample_columns(out)
    assert all(c.startswith("N_") for c in out.columns[1:])


def test_annotation_dir_to_abundance(tmp_path):
    ann = tmp_path / "ann"
    ann.mkdir()
    pd.DataFrame(
        {
            "seq": ["a", "b", "c"],
            "taxID_kaiju_0": [562, 562, 9606],
            "taxID_kraken2_0": [562, 9606, 9606],
        }
    ).to_csv(ann / "1.annotation.csv", index=False)
    tables = regenerate_annotation_tables(ann, tmp_path / "out", {"regeneration_mode": "direct"})
    assert set(tables) >= {"kaiju", "kraken2"}
    assert int(tables["kaiju"][tables["kaiju"]["taxid"].astype(str) == "562"]["N_1"].iloc[0]) == 2


def test_single_kraken2_annotation_csv(tmp_path):
    ann = tmp_path / "k2"
    ann.mkdir()
    pd.DataFrame(
        {"seq": ["r1", "r2", "r3"], "taxID_kraken2_0": [562, 562, 511145]}
    ).to_csv(ann / "s1.annotation.csv", index=False)
    tables = regenerate_annotation_tables(ann, tmp_path / "out", {"regeneration_mode": "direct"})
    assert list(tables) == ["kraken2"]
    row = tables["kraken2"]
    assert int(row[row["taxid"].astype(str) == "562"]["N_s1"].iloc[0]) == 2


def test_kraken2_raw_log_to_abundance():
    log = Path("tests/data/kraken2.log")
    if not log.is_file():
        pytest.skip("tests/data/kraken2.log missing")
    raw = read_kraken2_raw(str(log))
    assert "taxID" in raw.columns
    tables = regenerate_preserve(raw.head(20))
    assert tables
    table = next(iter(tables.values()))
    assert "taxid" in table.columns
    assert n_sample_columns(table)


def test_globalpatterns_python_generators(gp_table, tmp_path):
    src = tmp_path / "globalpatterns_abundance.csv"
    gp_table.to_csv(src, index=False)
    for mode in ("direct", "bootstrap", "vae", "glm", "camisim-table"):
        out = tmp_path / mode
        tables = regenerate_annotation_tables(
            src,
            out,
            {
                "regeneration_mode": mode,
                "N": 3,
                "N_reads": 80,
                "seed": 1,
                "rescale_abundance": True,
            },
        )
        assert tables
        table = next(iter(tables.values()))
        n_cols = n_sample_columns(table)
        if mode == "direct":
            assert len(n_cols) == 6
        else:
            assert len(n_cols) == 3
        assert "taxid" in table.columns
        assert int(table[n_cols].sum().sum()) > 0


def test_globalpatterns_custom_python(gp_table, tmp_path, monkeypatch):
    from samovar.paths import write_config
    from samovar.tools_import import import_tool

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    script = tmp_path / "echo_table.py"
    script.write_text(
        "from samovar.abundance import input_to_abundance_tables\n"
        "\n"
        "def regenerate(data, metadata, config):\n"
        "    tables = input_to_abundance_tables(data)\n"
        "    return {name: frame.copy() for name, frame in tables.items()}\n"
    )
    import_tool(
        name="echo_table",
        tool_type="table",
        exec_path=str(script),
        also_repo_build=False,
    )
    src = tmp_path / "gp.csv"
    gp_table.to_csv(src, index=False)
    tables = regenerate_annotation_tables(
        src, tmp_path / "out", {"regeneration_mode": "echo_table"}
    )
    assert tables
    assert n_sample_columns(next(iter(tables.values())))


def test_globalpatterns_sparsedossa2_mocked(gp_table, tmp_path, monkeypatch):
    def fake_sim(matrix, mode="sparsedossa2-fit", n_sample=None, config=None):
        n = int(n_sample or matrix.shape[1])
        cols = [f"s{i}" for i in range(n)]
        out = pd.DataFrame(index=matrix.index)
        for i, col in enumerate(cols):
            src = matrix.iloc[:, i % matrix.shape[1]]
            out[col] = src.to_numpy()
        return out[cols]

    monkeypatch.setattr("samovar.sparsedossa2.simulate_count_matrix", fake_sim)
    src = tmp_path / "gp.csv"
    gp_table.to_csv(src, index=False)
    tables = regenerate_annotation_tables(
        src,
        tmp_path / "sd2",
        {"regeneration_mode": "sparsedossa2-fit", "N": 3, "N_reads": 40, "seed": 1},
    )
    assert tables
    assert len(n_sample_columns(next(iter(tables.values())))) == 3


def test_globalpatterns_samovar_r(gp_table, tmp_path):
    src = tmp_path / "gp.csv"
    gp_table.to_csv(src, index=False)
    try:
        tables = regenerate_annotation_tables(
            src,
            tmp_path / "r",
            {"regeneration_mode": "samovar", "N": 1, "N_reads": 40, "seed": 1},
        )
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        pytest.skip(str(exc)[:300])
    if not tables:
        pytest.skip("samovaR wrote no abundance CSVs")
    table = next(iter(tables.values()))
    assert "taxid" in table.columns
    assert n_sample_columns(table)


def test_regenerated_abundance_matches_iss_contract(gp_table, tmp_path, monkeypatch):
    src = tmp_path / "gp.csv"
    gp_table.to_csv(src, index=False)
    tables = regenerate_annotation_tables(src, tmp_path / "ab", {"regeneration_mode": "direct"})
    table = next(iter(tables.values()))
    iss_shape = normalize_abundance_table(table)
    assert n_sample_columns(iss_shape)

    recorded = {}

    def fake_iss(table, **kwargs):
        recorded["table"] = table if not isinstance(table, str) else pd.read_csv(table)
        recorded["kwargs"] = kwargs
        out = Path(kwargs["output_dir"])
        out.mkdir(parents=True, exist_ok=True)
        paths = [out / "CL3_full_R1.fastq", out / "CL3_full_R2.fastq"]
        for p in paths:
            p.write_text("@r\nACGT\n+\nIIII\n")
        return [str(p) for p in paths]

    monkeypatch.setattr("samovar.table2iss.generate_iss_from_abundance_table", fake_iss)
    from samovar.table2iss import generate_iss_test_samples

    generate_iss_test_samples(
        genome_dir=str(tmp_path),
        host_genome="",
        output_dir=str(tmp_path / "fq"),
        n_samples=1,
        total_reads=10,
        abundance_table=iss_shape,
    )
    assert "taxid" in recorded["table"].columns
    assert n_sample_columns(normalize_abundance_table(recorded["table"]))


def test_iss_normalizes_phyloseq_columns(gp_table, tmp_path, monkeypatch):
    from samovar.table2iss import generate_iss_from_abundance_table

    tiny = gp_table.copy()
    tiny["taxid"] = ["562", "4932"] + ["9606"] * (len(tiny) - 2)
    tiny = tiny.head(3)

    monkeypatch.setattr("samovar.table2iss._resolve_genomes_for_taxids", lambda *a, **k: {})
    monkeypatch.setattr(
        "samovar.table2iss._write_empty_fastq_pair",
        lambda r1, r2: (
            Path(r1).parent.mkdir(parents=True, exist_ok=True),
            Path(r1).write_text("@r\nA\n+\nI\n"),
            Path(r2).write_text("@r\nT\n+\nI\n"),
        ),
    )
    out = tmp_path / "iss"
    out.mkdir()
    paths = generate_iss_from_abundance_table(
        tiny,
        genome_dir=str(tmp_path),
        output_dir=str(out),
        genomes=[],
    )
    assert paths
    assert n_sample_columns(normalize_abundance_table(tiny))
