"""samovar branch copies a run up to --start-point; prepare cannot change annotators."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest
import yaml

from samovar.branch import (
    copy_run_until,
    main as branch_main,
    paths_produced_before,
    reject_annotator_change,
)
from samovar.config import setup_pipeline
from samovar.exec_control import mark_done


def _ns(**kwargs):
    base = {
        "input_config": None,
        "input_dir": None,
        "output_dir": None,
        "add_annotator": False,
        "kraken2": None,
        "kaiju": None,
        "dummy": None,
        "cores": 1,
        "no_multiqc": True,
        "run_multiqc": False,
        "startpoint": None,
        "table_reads_generator": None,
        "use_test_genomes": False,
    }
    base.update(kwargs)
    return Namespace(**base)


def _seed_run(root: Path, *, dummy: bool = True) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    reads = root / "reads"
    reads.mkdir()
    (reads / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    (reads / "1_full_R2.fastq").write_text("@r\nT\n+\nI\n")
    out = root / "run"
    setup_pipeline(
        _ns(
            input_dir=str(reads),
            output_dir=str(out),
            dummy=[["dummy"]] if dummy else None,
            run_multiqc=False,
        )
    )
    for step in (
        "setup_reads",
        "qc_initial",
        "annotate_initial",
        "combine_initial",
        "viz_initial",
        "abundance_tables",
        "regenerate_tables",
        "seed_genomes",
    ):
        mark_done(out, step)
    (out / "initial").mkdir(exist_ok=True)
    (out / "initial" / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    (out / "initial" / "1_full_R2.fastq").write_text("@r\nT\n+\nI\n")
    (out / "initial_trimmed").mkdir(exist_ok=True)
    (out / "initial_reports").mkdir(exist_ok=True)
    (out / "initial_reports" / "1.out").write_text("ok\n")
    (out / "initial_annotations").mkdir(exist_ok=True)
    (out / "initial_annotations" / "1.csv").write_text("taxid,n\n1,1\n")
    abund = out / "initial_abundance"
    abund.mkdir(exist_ok=True)
    (abund / "dummy.csv").write_text("taxid,N_1_full\n1,10\n")
    regen = out / "regenerated" / ".regenerated_abundance"
    regen.mkdir(parents=True)
    (regen / "dummy.csv").write_text("taxid,N_1_full\n1,10\n")
    (out / "regenerated" / "1_full_dummy_R1.fastq").write_text("@r\nA\n+\nI\n")
    (out / "regenerated_annotations").mkdir(exist_ok=True)
    (out / "regenerated_annotations" / "1.csv").write_text("taxid,n\n1,1\n")
    return out


def test_paths_produced_before_skips_later_fastq():
    rels = paths_produced_before("regenerate_tables")
    assert "initial_abundance" in rels
    assert "initial" in rels
    assert "regenerated" not in rels
    assert not any(r.startswith("regenerated") for r in rels)


def test_copy_run_until_keeps_early_drops_later(tmp_path):
    src = _seed_run(tmp_path / "src")
    dest = tmp_path / "branch"
    assert branch_main(
        ["--start-point", "regenerate_tables", str(src), str(dest)]
    ) == 0
    assert (dest / "initial_abundance" / "dummy.csv").is_file()
    assert (dest / "initial" / "1_full_R1.fastq").is_file()
    assert not (dest / "regenerated" / "1_full_dummy_R1.fastq").exists()
    assert not (dest / "regenerated" / ".regenerated_abundance" / "dummy.csv").exists()
    assert (dest / ".log" / "checkpoints" / "abundance_tables.done").is_file()
    assert not (dest / ".log" / "checkpoints" / "regenerate_tables.done").exists()
    assert not (dest / ".log" / "checkpoints" / "seed_genomes.done").exists()
    meta = yaml.safe_load((dest / ".log" / "branch.yaml").read_text())
    assert meta["startpoint"] == "regenerate_tables"
    assert "dummy" in meta["annotators"]


def test_branch_prepare_merges_hydra_rejects_new_annotator(tmp_path):
    src = _seed_run(tmp_path / "src")
    dest = tmp_path / "branch"
    copy_run_until(src, dest, "regenerate_tables")
    original_abund = (dest / "initial_abundance" / "dummy.csv").read_text()

    second = setup_pipeline(
        _ns(
            output_dir=str(dest),
            table_reads_generator=["bootstrap"],
            run_multiqc=False,
        )
    )
    assert Path(second["pipeline"]).name == "samovar.sh"
    window = (dest / ".log" / "window.env").read_text()
    assert "SAMOVAR_START=regenerate_tables" in window
    assert (dest / "initial_abundance" / "dummy.csv").read_text() == original_abund
    iss = (dest / ".log" / "configs" / "config_annotation2iss.yaml").read_text()
    assert "bootstrap" in iss

    with pytest.raises(ValueError, match="cannot change initial annotators"):
        reject_annotator_change(
            dest,
            _ns(output_dir=str(dest), **{"cmd_kraken2-test": [["kraken2 toy"]]}),
        )

    with pytest.raises(SystemExit):
        setup_pipeline(
            _ns(
                output_dir=str(dest),
                **{"cmd_kraken2-test": [["kraken2 toy"]]},
            )
        )
