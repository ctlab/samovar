"""Tests for metagenome annotation table regeneration modes."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
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


def _cosine(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def test_bootstrap_mimics_real_samples_with_error(test_data_dir):
    if not test_data_dir.exists():
        pytest.skip("data/test_annotations missing")
    data = read_annotation_dir(test_data_dir)
    n_obs = max(int(data["sample"].nunique()) if "sample" in data.columns else 1, 1)
    direct = regenerate_preserve(data, rescale=False)
    boot = regenerate_bootstrap(
        data, n_samples=n_obs, seed=3, rescale=False, error_scale=0.08
    )
    differed = False
    similar = False
    for name in direct:
        d = _n_matrix(direct[name])
        b = _n_matrix(boot[name])
        shared = [c for c in d.columns if c in b.columns]
        assert shared
        if not d[shared].equals(b[shared]):
            differed = True
        for col in shared:
            idx = d.index.union(b.index)
            u = d[col].reindex(idx).fillna(0)
            v = b[col].reindex(idx).fillna(0)
            if _cosine(u, v) >= 0.75:
                similar = True
    assert differed
    assert similar


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


def test_glm_uses_cluster_walk_not_copy(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    glm = regenerate_glm_python(
        data, n_samples=2, seed=11, min_cluster_size=1, rescale=False
    )["kaiju"]
    n_cols = [c for c in glm.columns if str(c).startswith("N_")]
    assert n_cols
    assert (glm[n_cols].fillna(0) >= 0).all().all()
    assert glm["taxid"].nunique() >= 1


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


def test_python_modes_are_self_contained(toy_annotation_dir, tmp_path, monkeypatch):
    monkeypatch.delenv("SAMOVAR_R_REGENERATE", raising=False)
    monkeypatch.delenv("SAMOVAR_ANNOTATION_REGENERATE_R", raising=False)
    monkeypatch.setattr("samovar.table2iss.annotation_regenerate_r", lambda: None)
    for mode in ("direct", "bootstrap", "vae", "glm"):
        out = tmp_path / mode
        tables = regenerate_annotation_tables(
            toy_annotation_dir,
            out,
            {"regeneration_mode": mode, "N": 2, "seed": 1},
        )
        assert tables
        for table in tables.values():
            assert "taxid" in table.columns
            n_cols = [c for c in table.columns if str(c).startswith("N_")]
            assert n_cols
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

    with patch("samovar.table2iss.subprocess.run", side_effect=_fake_iss) as mocked:
        process_annotation_tables(
            table_paths=[str(p) for p in ann_files],
            genome_dir=str(genome_dir),
            output_dir=str(reads),
            annotation_dir=str(toy_annotation_dir),
            regeneration_config=regen_cfg,
            seed=regen_cfg["seed"],
        )
    iss_cmds = [
        list(c.args[0])
        for c in mocked.call_args_list
        if c.args and isinstance(c.args[0], (list, tuple)) and c.args[0] and Path(str(c.args[0][0])).name == "iss"
    ]
    assert iss_cmds
    for cmd in iss_cmds:
        assert "--seed" in cmd

    for sample in expected_samples:
        for annotator in ("kaiju", "kraken2"):
            assert (reads / f"{sample}_{annotator}_R1.fastq").exists()
            assert (reads / f"{sample}_{annotator}_R2.fastq").exists()

    if mode == "direct":
        sample_tables = _sample_tables_from_abundance_dir(
            reads / ".regenerated_abundance", None
        )
        assert sorted(sample_tables) == sorted(expected_samples)
        direct_from_disk = pd.read_csv(reads / ".regenerated_abundance" / "kaiju.csv")
        assert _n_matrix(direct_tables["kaiju"]).equals(_n_matrix(direct_from_disk))
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


def test_bootstrap_seed_reproducible_but_not_constant():
    rng = np.random.default_rng(0)
    taxa = [562, 9606, 28901, 1423]
    rows = []
    for sample in ("1", "2", "3"):
        for _ in range(40):
            rows.append(
                {
                    "seq": f"{sample}_{len(rows)}",
                    "sample": sample,
                    "taxID_kaiju_0": int(rng.choice(taxa, p=[0.5, 0.3, 0.15, 0.05])),
                }
            )
    data = pd.DataFrame(rows)
    a = regenerate_bootstrap(data, n_samples=4, seed=11, rescale=False, error_scale=0.2)
    b = regenerate_bootstrap(data, n_samples=4, seed=11, rescale=False, error_scale=0.2)
    c = regenerate_bootstrap(data, n_samples=4, seed=99, rescale=False, error_scale=0.2)
    ma, mb, mc = _n_matrix(a["kaiju"]), _n_matrix(b["kaiju"]), _n_matrix(c["kaiju"])
    assert ma.equals(mb)
    assert not ma.equals(mc)
    profiles = {tuple(ma[col].tolist()) for col in ma.columns}
    assert len(profiles) > 1


def test_generative_modes_differ_from_each_other_and_direct(toy_annotation_dir):
    data = read_annotation_dir(toy_annotation_dir)
    tables = {
        "direct": regenerate_preserve(data, rescale=False)["kaiju"],
        "bootstrap": regenerate_bootstrap(data, n_samples=3, seed=5, rescale=False)["kaiju"],
        "vae": regenerate_vae(data, n_samples=3, seed=5, rescale=False)["kaiju"],
        "glm": regenerate_glm_python(data, n_samples=3, seed=5, rescale=False)["kaiju"],
    }
    matrices = {k: _n_matrix(v) for k, v in tables.items()}
    keys = list(matrices)
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            assert not matrices[a].equals(matrices[b]), f"{a} == {b}"


def test_iss_readcounts_match_regenerated_glm_tables(
    toy_annotation_dir, tmp_path
):
    from unittest.mock import patch

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for taxid in ("562", "9606"):
        (genome_dir / f"{taxid}.fa").write_text(f">{taxid}\nATCGATCG\n")
    reads = tmp_path / "reads"
    captured = {}

    def fake(cmd, **kwargs):
        if isinstance(cmd, (list, tuple)) and "--readcount_file" in cmd:
            path = Path(cmd[cmd.index("--readcount_file") + 1])
            captured[cmd[cmd.index("--output") + 1]] = path.read_text()
        return _fake_iss(cmd, **kwargs)

    regen_cfg = {"regeneration_mode": "glm", "N": 2, "seed": 8, "rescale_abundance": False}
    with patch("samovar.table2iss.subprocess.run", side_effect=fake):
        process_annotation_tables(
            table_paths=[str(p) for p in sorted(toy_annotation_dir.glob("*.annotation.csv"))],
            genome_dir=str(genome_dir),
            output_dir=str(reads),
            annotation_dir=str(toy_annotation_dir),
            regeneration_config=regen_cfg,
            seed=8,
        )
    glm = pd.read_csv(reads / ".regenerated_abundance" / "kaiju.csv")
    n_cols = [c for c in glm.columns if str(c).startswith("N_")]
    expected_total = int(glm[n_cols].sum().sum())
    kaiju_key = next(k for k in captured if "kaiju" in k)
    iss_total = sum(int(line.split()[-1]) for line in captured[kaiju_key].splitlines() if line.strip())
    # Abundance rows count pairs; ISS's readcount file counts individual reads.
    assert iss_total == expected_total * 2
    assert expected_total > 0

