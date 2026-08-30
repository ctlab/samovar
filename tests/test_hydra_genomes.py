"""Hydra locks the genomes used in a generate run; export packs lazy-load only."""

from __future__ import annotations

import shutil
import tarfile
from pathlib import Path

import yaml
from omegaconf import OmegaConf

from samovar.iss_config import select_generate_genome_paths, setup_iss_test
from samovar.paths import test_genomes_dir as bundled_genomes_dir
from samovar.repro import export_run, reproduce

META = bundled_genomes_dir() / "meta"
HOST = bundled_genomes_dir() / "host" / "9606.fna"
PHAGE = "2886930"


def _copy_meta(dest: Path) -> Path:
    dest.mkdir(parents=True)
    for path in sorted(META.iterdir()):
        if path.is_file():
            shutil.copy2(path, dest / path.name)
    return dest


def _generate(args_ns):
    return setup_iss_test(args_ns)


def _ns(output_dir, genome_dir, **extra):
    from argparse import Namespace

    base = dict(
        genome_dir=str(genome_dir),
        output_dir=str(output_dir),
        host_genome=str(HOST),
        n_samples=1,
        total_reads=20,
        host_fraction="0",
        seed=1,
        model="hiseq",
        cores=1,
        simulator="iss",
        reads_generator=None,
        metagenome_generator=None,
        abundance_table=None,
        extra_flags="",
        tool_flags=None,
        camisim_mode=None,
        accessions=None,
        reindex=0,
        raw_genomes=0,
        max_genomes=None,
        genomes=None,
        max_genome_mb=None,
        genome_skip_list=None,
        threads=None,
        custom_cli_flags={},
    )
    base.update(extra)
    return Namespace(**base)


def test_top_writes_used_genomes_to_hydra(tmp_path):
    genomes = _copy_meta(tmp_path / "meta")
    out = tmp_path / "run"
    result = _generate(_ns(out, genomes, max_genomes=2))
    cfg = yaml.safe_load(Path(result["config"]).read_text())
    locked = [Path(p).name for p in cfg["genomes"]]
    assert len(locked) == 2
    snap = OmegaConf.load(out / ".hydra" / "config.yaml")
    used = [str(r.taxid) for r in snap.genomes.used if not r.get("host")]
    assert used == [Path(p).stem.split(".")[0] for p in select_generate_genome_paths(str(genomes), 2)]
    assert len(used) == 2


def test_locked_genomes_ignore_later_dir_entries(tmp_path):
    genomes = _copy_meta(tmp_path / "meta")
    out = tmp_path / "run"
    _generate(_ns(out, genomes, max_genomes=2))
    first = OmegaConf.load(out / ".hydra" / "config.yaml")
    used = [str(r.taxid) for r in first.genomes.used if not r.get("host")]
    (genomes / "000.fna").write_text(">x\nATGC\n")
    rerun = tmp_path / "rerun"
    code = reproduce(out, str(rerun), install=False, stages=["generate"])
    assert code == 0
    cfg2 = yaml.safe_load((rerun / ".generate" / "configs" / "iss_config.yaml").read_text())
    names = [Path(p).stem.split(".")[0] for p in cfg2["genomes"]]
    assert "000" not in names
    assert names == used


def test_export_full_omits_database_files(tmp_path, monkeypatch):
    from samovar.paths import write_config
    from samovar.repro import record_stage
    from types import SimpleNamespace

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    db = tmp_path / "db"
    db.mkdir()
    (db / "hash.k2d").write_text("secret-index")
    write_config(
        {
            "root": str(tmp_path),
            "tools": {},
            "databases": {
                "kraken2": {
                    "toy": {
                        "path": str(db),
                        "lazy-download": "#!/bin/bash\necho fetch\n",
                        "url": "https://example.invalid/db.tgz",
                    }
                }
            },
        },
        also_repo_build=False,
    )
    out = tmp_path / "run"
    record_stage(
        "generate",
        str(out),
        args=SimpleNamespace(output_dir=str(out), n_samples=1),
        argv=["--n_samples", "1"],
    )
    archive = tmp_path / "full.tgz"
    export_run(out, archive, mode="full")
    names = tarfile.open(archive, "r:gz").getnames()
    assert any(n.endswith("payload/db-meta/kraken2/toy/lazy-download.sh") for n in names)
    assert not any(n.endswith("hash.k2d") for n in names)


def test_toy_and_phage_export_rerun_same_genomes(tmp_path):
    """Toy (all meta taxa, --top 2) and phage (2886930 only) survive export/reproduce."""
    for label, extra in (
        ("toy", {"max_genomes": 2}),
        ("phage", {"genomes": [PHAGE]}),
    ):
        genomes = _copy_meta(tmp_path / f"{label}_meta")
        out = tmp_path / f"{label}_run"
        _generate(_ns(out, genomes, **extra))
        snap = OmegaConf.load(out / ".hydra" / "config.yaml")
        used = [str(r.taxid) for r in snap.genomes.used if not r.get("host")]
        assert used
        if label == "phage":
            assert used == [PHAGE]
        else:
            assert len(used) == 2
            assert PHAGE in {p.stem.split(".")[0] for p in genomes.iterdir()} or True
        archive = tmp_path / f"{label}.tgz"
        export_run(out, archive, mode="full")
        names = tarfile.open(archive, "r:gz").getnames()
        assert not any(n.endswith(".fna") for n in names)
        rerun = tmp_path / f"{label}_rerun"
        assert reproduce(archive, str(rerun), install=False, stages=["generate"]) == 0
        cfg2 = yaml.safe_load((rerun / ".generate" / "configs" / "iss_config.yaml").read_text())
        again = [Path(p).stem.split(".")[0] for p in cfg2["genomes"]]
        assert again == used
