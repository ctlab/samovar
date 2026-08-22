"""Tests for metagenome annotation table regeneration modes."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

from samovar.annotation_io import read_annotation_dir
from samovar.regenerate import (
    is_direct_mode,
    normalize_regeneration_mode,
    regenerate_annotation_tables,
    regenerate_bootstrap,
    regenerate_glm_python,
    regenerate_preserve,
    regenerate_vae,
    sample_names_from_abundance_columns,
    synthetic_sample_names,
    write_samovar_config_defaults,
    _correlation_matrix,
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


def _n_matrix(df):
    out = df.copy()
    out["taxid"] = out["taxid"].astype(str)
    cols = sorted(c for c in out.columns if str(c).startswith("N_"))
    return out.set_index("taxid")[cols].sort_index()


def _fake_iss(cmd, **kwargs):
    class Result:
        returncode = 0

    if isinstance(cmd, (list, tuple)) and cmd and Path(str(cmd[0])).name == "iss":
        output = cmd[cmd.index("--output") + 1]
        os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
        with open(f"{output}_R1.fastq", "w") as r1, open(f"{output}_R2.fastq", "w") as r2:
            r1.write("@read/1\nACGT\n+\nIIII\n")
            r2.write("@read/2\nTGCA\n+\nIIII\n")
    return Result()


def test_normalize_regeneration_mode_aliases():
    assert normalize_regeneration_mode(None) == "direct"
    assert normalize_regeneration_mode("preserve") == "direct"
    assert normalize_regeneration_mode("exact") == "direct"
    assert is_direct_mode("preserve")
    assert normalize_regeneration_mode("glm") == "glm"
    assert normalize_regeneration_mode("samovar") == "samovar"
    assert normalize_regeneration_mode("boil") == "samovar"
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
    tables = regenerate_bootstrap(
        data, n_samples=4, n_reads=200, seed=1, rescale=True
    )
    kaiju = tables["kaiju"]
    n_cols = [c for c in kaiju.columns if c.startswith("N_")]
    assert len(n_cols) == 4
    for col in n_cols:
        assert int(kaiju[col].sum()) == 200


def test_bootstrap_does_not_rescale_by_default(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    tables = regenerate_bootstrap(
        data, n_samples=4, n_reads=200, seed=1, rescale=False
    )
    kaiju = tables["kaiju"]
    n_cols = [c for c in kaiju.columns if c.startswith("N_")]
    assert len(n_cols) == 4
    for col in n_cols:
        assert int(kaiju[col].sum()) != 200


def test_vae_produces_n_synthetic_columns(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    tables = regenerate_vae(data, n_samples=3, n_reads=150, seed=2, latent_dim=2)
    kaiju = tables["kaiju"]
    n_cols = [c for c in kaiju.columns if c.startswith("N_")]
    assert len(n_cols) == 3
    assert kaiju["taxid"].nunique() >= 1


def test_glm_python_changes_profile_vs_direct(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    direct = regenerate_preserve(data, rescale=False)["kaiju"]
    glm = regenerate_glm_python(data, n_samples=2, n_reads=100, seed=3)["kaiju"]
    assert glm.shape[1] == 3
    assert "N_1" in glm.columns
    assert "N_2" in glm.columns
    assert not _n_matrix(direct).equals(_n_matrix(glm))


def test_generative_modes_do_not_force_n_reads(toy_annotation_dir, tmp_path):
    for mode in ("bootstrap", "vae", "glm"):
        tables = regenerate_annotation_tables(
            toy_annotation_dir,
            tmp_path / mode,
            {
                "regeneration_mode": mode,
                "N": 2,
                "N_reads": 999,
                "seed": 0,
                "rescale_abundance": False,
            },
        )
        for table in tables.values():
            n_cols = [c for c in table.columns if str(c).startswith("N_")]
            for col in n_cols:
                assert int(table[col].sum()) != 999


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


def test_regenerate_annotation_tables_samovar_is_not_python(toy_annotation_dir, tmp_path):
    with pytest.raises(ValueError, match="optional R regenerator"):
        regenerate_annotation_tables(
            toy_annotation_dir,
            tmp_path / "out",
            {"regeneration_mode": "samovar"},
        )


def test_samovar_annotation_regenerate_direct_integration(toy_annotation_dir, tmp_path):
    out = tmp_path / "out"
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(
        "regeneration_mode: direct\nN_reads: 1000\noutput_dir: "
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


def test_real_toy_annotations_direct(test_data_dir, tmp_path):
    if not test_data_dir.exists():
        pytest.skip("data/test_annotations missing")
    out = tmp_path / "direct"
    regenerate_annotation_tables(
        test_data_dir,
        out,
        write_samovar_config_defaults({"regeneration_mode": "direct"}),
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
                "rescale_abundance": True,
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


def test_synthetic_sample_names_stable_with_n():
    assert synthetic_sample_names(["1", "2"], None) == ["1", "2"]
    assert synthetic_sample_names(["1", "2"], 2) == ["1", "2"]
    assert synthetic_sample_names(["1", "2"], 1) == ["1"]
    assert synthetic_sample_names(["1", "2"], 5) == ["1", "2", "1_r2", "2_r2", "1_r3"]


def test_glm_corrcoef_single_taxon_is_finite():
    data = pd.DataFrame(
        {
            "sample": ["a", "a", "b"],
            "taxID_kaiju_0": [562, 562, 562],
        }
    )
    tables = regenerate_glm_python(data, n_samples=3, n_reads=20, seed=0)
    kaiju = tables["kaiju"]
    n_cols = [c for c in kaiju.columns if c.startswith("N_")]
    assert n_cols == ["N_a", "N_b", "N_a_r2"]
    assert not kaiju[n_cols].isna().any().any()


def test_correlation_matrix_guards_nan():
    import numpy as np

    constant = np.ones((2, 3))
    corr = _correlation_matrix(np.log1p(constant))
    assert corr.shape == (2, 2)
    assert np.isfinite(corr).all()
    assert corr[0, 0] == 1.0
    single = _correlation_matrix(np.log1p(np.array([[1.0, 2.0, 3.0]])))
    assert single.shape == (1, 1)
    assert single[0, 0] == 1.0


@pytest.mark.parametrize(
    "mode,n_samples,expected_samples",
    [
        ("direct", None, ["1", "2"]),
        ("bootstrap", 3, ["1", "2", "1_r2"]),
        ("vae", 3, ["1", "2", "1_r2"]),
        ("glm", 3, ["1", "2", "1_r2"]),
    ],
)
def test_process_annotation_tables_modes_compatible_with_iss(
    toy_annotation_dir, tmp_path, mode, n_samples, expected_samples
):
    from unittest.mock import patch

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for taxid in ("562", "9606"):
        (genome_dir / f"{taxid}.fa").write_text(f">{taxid}\nATCGATCG\n")

    reads = tmp_path / "reads"
    ann_files = sorted(toy_annotation_dir.glob("*.annotation.csv"))
    regen_cfg = {
        "regeneration_mode": mode,
        "N_reads": 40,
        "seed": 0,
        "rescale_abundance": False,
    }
    if n_samples is not None:
        regen_cfg["N"] = n_samples

    direct_tables = regenerate_annotation_tables(
        toy_annotation_dir,
        tmp_path / "direct_tables",
        {"regeneration_mode": "direct"},
    )
    if mode != "direct":
        gen_tables = regenerate_annotation_tables(
            toy_annotation_dir, tmp_path / "gen_tables", regen_cfg
        )
        differed = False
        for name in direct_tables:
            if not _n_matrix(direct_tables[name]).equals(_n_matrix(gen_tables[name])):
                differed = True
            n_cols = [c for c in gen_tables[name].columns if str(c).startswith("N_")]
            assert "taxid" in gen_tables[name].columns
            assert n_cols
        assert differed, f"{mode} tables should differ from direct"

    with patch("samovar.table2iss.subprocess.run", side_effect=_fake_iss):
        process_annotation_tables(
            table_paths=[str(p) for p in ann_files],
            genome_dir=str(genome_dir),
            output_dir=str(reads),
            annotation_dir=str(toy_annotation_dir),
            regeneration_config=regen_cfg,
        )

    for sample in expected_samples:
        for annotator in ("kaiju", "kraken2"):
            assert (reads / f"{sample}_{annotator}_R1.fastq").exists()
            assert (reads / f"{sample}_{annotator}_R2.fastq").exists()

    if mode == "direct":
        assert not (reads / ".regenerated_abundance").exists()
    else:
        sample_tables = _sample_tables_from_abundance_dir(
            reads / ".regenerated_abundance", None
        )
        assert sorted(sample_tables) == sorted(expected_samples)


def test_process_annotation_tables_samovar_mode_feeds_iss(
    toy_annotation_dir, tmp_path, monkeypatch
):
    from unittest.mock import patch

    script = tmp_path / "annotation_regenerate.R"
    script.write_text("# optional R regenerator\n")
    monkeypatch.setattr("samovar.table2iss.annotation_regenerate_r", lambda: script)

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for taxid in ("562", "9606"):
        (genome_dir / f"{taxid}.fa").write_text(f">{taxid}\nATCGATCG\n")

    reads = tmp_path / "reads"
    ann_files = sorted(toy_annotation_dir.glob("*.annotation.csv"))

    def fake_run(cmd, **kwargs):
        if isinstance(cmd, (list, tuple)) and cmd and Path(str(cmd[0])).name == "iss":
            return _fake_iss(cmd, **kwargs)
        output_dir = cmd[cmd.index("--output_dir") + 1]
        os.makedirs(output_dir, exist_ok=True)
        pd.DataFrame(
            {"taxid": [562, 9606], "N_1": [12, 8], "N_2": [3, 17]}
        ).to_csv(Path(output_dir) / "kaiju.csv", index=False)
        pd.DataFrame(
            {"taxid": [562, 9606], "N_1": [1, 19], "N_2": [10, 10]}
        ).to_csv(Path(output_dir) / "kraken2.csv", index=False)

        class Result:
            returncode = 0

        return Result()

    with patch("samovar.table2iss._resolve_r_executable", return_value=("R", None)), \
         patch("samovar.table2iss.subprocess.run", side_effect=fake_run):
        process_annotation_tables(
            table_paths=[str(p) for p in ann_files],
            genome_dir=str(genome_dir),
            output_dir=str(reads),
            annotation_dir=str(toy_annotation_dir),
            regeneration_config={"regeneration_mode": "samovar", "N": 2},
        )

    for sample in ("1", "2"):
        for annotator in ("kaiju", "kraken2"):
            assert (reads / f"{sample}_{annotator}_R1.fastq").exists()
    sample_tables = _sample_tables_from_abundance_dir(
        reads / ".regenerated_abundance", None
    )
    assert sorted(sample_tables) == ["1", "2"]
    kaiju = pd.read_csv(reads / ".regenerated_abundance" / "kaiju.csv")
    direct = regenerate_annotation_tables(
        toy_annotation_dir,
        tmp_path / "direct_tables",
        {"regeneration_mode": "direct"},
    )["kaiju"]
    assert not _n_matrix(direct).equals(_n_matrix(kaiju))
