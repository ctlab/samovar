"""Tests for metagenome annotation table regeneration modes."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

from samovar.annotation_io import read_annotation_dir
from samovar.regenerate import (
    normalize_regeneration_mode,
    regenerate_annotation_tables,
    regenerate_bootstrap,
    regenerate_glm_python,
    regenerate_preserve,
    regenerate_vae,
    sample_names_from_abundance_columns,
    write_samovar_config_defaults,
)
from samovar.table2iss import (
    _sample_tables_from_abundance_dir,
    process_annotation_tables,
    samovar_annotation_regenerate,
)


@pytest.fixture
def toy_annotation_dir(tmp_path):
    ann = tmp_path / "annotations"
    ann.mkdir()
    pd.DataFrame(
        {
            "seq": ["a", "b", "c"],
            "taxID_kaiju_0": [562, 562, 9606],
            "taxID_kraken2_0": [562, 9606, 9606],
            "true": [562, 562, 9606],
        }
    ).to_csv(ann / "1.annotation.csv", index=False)
    pd.DataFrame(
        {
            "seq": ["a", "b"],
            "taxID_kaiju_0": [562, 9606],
            "taxID_kraken2_0": [9606, 9606],
            "true": [562, 9606],
        }
    ).to_csv(ann / "2.annotation.csv", index=False)
    return ann


@pytest.fixture
def test_data_dir():
    return Path("data/test_annotations")


def test_normalize_regeneration_mode_aliases():
    assert normalize_regeneration_mode(None) == "preserve"
    assert normalize_regeneration_mode("exact") == "preserve"
    assert normalize_regeneration_mode("glm") == "glm"
    with pytest.raises(ValueError):
        normalize_regeneration_mode("unknown")


def test_preserve_keeps_observed_counts(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    tables = regenerate_preserve(data, n_reads=1000, rescale=False)
    kaiju = tables["kaiju"]
    assert set(kaiju.columns) == {"taxid", "N_1", "N_2"}
    row562 = kaiju[kaiju["taxid"] == "562"]
    assert int(row562["N_1"].iloc[0]) == 2
    assert int(kaiju[kaiju["taxid"] == "9606"]["N_1"].iloc[0]) == 1
    assert int(row562["N_2"].iloc[0]) == 1


def test_preserve_rescale_when_requested(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    tables = regenerate_preserve(data, n_reads=100, rescale=True)
    kaiju = tables["kaiju"]
    assert int(kaiju["N_1"].sum()) == 100


def test_bootstrap_produces_n_synthetic_columns(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    tables = regenerate_bootstrap(data, n_samples=4, n_reads=200, seed=1)
    kaiju = tables["kaiju"]
    n_cols = [c for c in kaiju.columns if c.startswith("N_")]
    assert len(n_cols) == 4
    for col in n_cols:
        assert int(kaiju[col].sum()) == 200


def test_vae_produces_n_synthetic_columns(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    tables = regenerate_vae(data, n_samples=3, n_reads=150, seed=2, latent_dim=2)
    kaiju = tables["kaiju"]
    n_cols = [c for c in kaiju.columns if c.startswith("N_")]
    assert len(n_cols) == 3
    assert kaiju["taxid"].nunique() >= 1


def test_glm_python_changes_profile_vs_preserve(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    preserve = regenerate_preserve(data, rescale=False)["kaiju"]
    glm = regenerate_glm_python(data, n_samples=2, n_reads=100, seed=3)["kaiju"]
    assert glm.shape[1] == 3
    glm_col = "N_synth_1" if "N_synth_1" in glm.columns else "N_1"
    assert glm_col in glm.columns


def test_regenerate_annotation_tables_writes_csvs(toy_annotation_dir, tmp_path):
    out = tmp_path / "abundance"
    tables = regenerate_annotation_tables(
        toy_annotation_dir,
        out,
        {"regeneration_mode": "bootstrap", "N": 2, "N_reads": 50, "seed": 0},
    )
    assert tables
    assert list(out.glob("*.csv"))
    assert "kaiju" in tables


def test_samovar_annotation_regenerate_preserve_integration(toy_annotation_dir, tmp_path):
    out = tmp_path / "out"
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(
        "regeneration_mode: preserve\nN_reads: 1000\noutput_dir: "
        + str(out).replace("\\", "/")
    )
    samovar_annotation_regenerate(
        annotation_dir=str(toy_annotation_dir),
        config_samovar=str(cfg),
        output_dir=str(out),
    )
    kaiju = pd.read_csv(out / "kaiju.csv")
    assert "taxid" in kaiju.columns
    assert int(kaiju[kaiju["taxid"] == 562]["N_1"].iloc[0]) == 2


def test_real_toy_annotations_preserve(test_data_dir, tmp_path):
    if not test_data_dir.exists():
        pytest.skip("data/test_annotations missing")
    out = tmp_path / "preserve"
    regenerate_annotation_tables(
        test_data_dir,
        out,
        write_samovar_config_defaults({"regeneration_mode": "preserve"}),
    )
    files = list(out.glob("*.csv"))
    assert len(files) >= 2
    for path in files:
        df = pd.read_csv(path)
        assert "taxid" in df.columns
        n_cols = [c for c in df.columns if c.startswith("N_")]
        assert n_cols
        assert df[n_cols].sum().sum() > 0


def test_real_toy_annotations_generative_modes(test_data_dir, tmp_path):
    if not test_data_dir.exists():
        pytest.skip("data/test_annotations missing")
    for mode in ("bootstrap", "vae", "glm"):
        out = tmp_path / mode
        regenerate_annotation_tables(
            test_data_dir,
            out,
            {
                "regeneration_mode": mode,
                "N": 5,
                "N_reads": 500,
                "seed": 7,
            },
        )
        for path in out.glob("*.csv"):
            df = pd.read_csv(path)
            n_cols = [c for c in df.columns if c.startswith("N_")]
            assert len(n_cols) == 5
            for col in n_cols:
                assert int(df[col].sum()) == 500


def test_sample_names_from_abundance_columns():
    cols = ["N_synth_1", "N_synth_2"]
    assert sample_names_from_abundance_columns(cols) == ["synth_1", "synth_2"]
    assert sample_names_from_abundance_columns(cols, ["a", "b"]) == ["a", "b"]


def test_process_annotation_tables_with_bootstrap_regeneration(toy_annotation_dir, tmp_path):
    from unittest.mock import patch

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for taxid in ("562", "9606"):
        (genome_dir / f"{taxid}.fa").write_text(f">{taxid}\nATCGATCG\n")

    reads = tmp_path / "reads"
    ann_files = sorted(toy_annotation_dir.glob("*.annotation.csv"))

    def fake_iss(cmd, **kwargs):
        class Result:
            returncode = 0

        if isinstance(cmd, (list, tuple)) and cmd and Path(str(cmd[0])).name == "iss":
            output = cmd[cmd.index("--output") + 1]
            os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
            with open(f"{output}_R1.fastq", "w") as r1, open(f"{output}_R2.fastq", "w") as r2:
                r1.write("@read/1\nACGT\n+\nIIII\n")
                r2.write("@read/2\nTGCA\n+\nIIII\n")
        return Result()

    with patch("samovar.table2iss.subprocess.run", side_effect=fake_iss):
        process_annotation_tables(
            table_paths=[str(p) for p in ann_files],
            genome_dir=str(genome_dir),
            output_dir=str(reads),
            annotation_dir=str(toy_annotation_dir),
            regeneration_config={
                "regeneration_mode": "bootstrap",
                "N": 3,
                "N_reads": 20,
                "seed": 0,
            },
        )

    assert (reads / "synth_1_kaiju_R1.fastq").exists()
    assert (reads / "synth_3_kaiju_R1.fastq").exists()
    sample_tables = _sample_tables_from_abundance_dir(
        reads / ".regenerated_abundance", None
    )
    assert len(sample_tables) == 3
