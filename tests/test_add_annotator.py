"""prepare --add-annotator writes samovar_vN.sh and keeps original Hydra params."""

from __future__ import annotations

import os
import stat
import subprocess
from argparse import Namespace
from pathlib import Path

import yaml

from samovar.add_annotator import (
    active_pipeline_script,
    annotator_fastqs_complete,
    existing_run_names,
    invalidate_add_annotator,
    next_pipeline_version,
)
from samovar.config import setup_pipeline
from samovar.exec_control import listed_done, mark_done

REPO = Path(__file__).resolve().parents[1]


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
    }
    base.update(kwargs)
    return Namespace(**base)


def test_add_annotator_writes_versioned_script(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    (reads / "1_full_R2.fastq").write_text("@r\nT\n+\nI\n")
    out = tmp_path / "run"
    first = setup_pipeline(
        _ns(
            input_dir=str(reads),
            output_dir=str(out),
            dummy=[["dummy"]],
        )
    )
    original = Path(first["pipeline"])
    assert original.name == "samovar.sh"
    original_text = original.read_text()
    init = yaml.safe_load((out / ".log" / "configs" / "config_init.yaml").read_text())
    assert {r["run_name"] for r in init["run_config"]} == {"dummy"}

    mark_done(out, "setup_reads")
    mark_done(out, "qc_initial")
    mark_done(out, "annotate_initial")
    mark_done(out, "seed_genomes")
    (out / "initial_abundance").mkdir(parents=True)
    (out / "initial_abundance" / "dummy.csv").write_text("taxid,N_1\n1,1\n")
    (out / "regenerated").mkdir(parents=True)
    (out / "regenerated" / ".process_annotations.done").write_text("ok\n")

    second = setup_pipeline(
        _ns(
            output_dir=str(out),
            add_annotator=True,
            **{"cmd_kraken2-test": [["kraken2 /tmp/kraken2_db"]]},
        )
    )
    versioned = Path(second["pipeline"])
    assert versioned.name == "samovar_v1.sh"
    assert versioned.is_file()
    assert original.read_text() == original_text
    assert next_pipeline_version(out) == 2
    assert active_pipeline_script(out).name == "samovar_v1.sh"

    init2 = yaml.safe_load((out / ".log" / "configs" / "config_init.yaml").read_text())
    names = {r["run_name"] for r in init2["run_config"]}
    assert names == {"dummy", "kraken2-test"}
    done = set(listed_done(out))
    assert "setup_reads" in done
    assert "qc_initial" in done
    assert "seed_genomes" in done
    assert "annotate_initial" not in done
    assert not (out / "regenerated" / ".process_annotations.done").is_file()
    assert not (out / "initial_abundance" / "dummy.csv").is_file()


def test_add_annotator_rejects_duplicate(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    (reads / "1_full_R2.fastq").write_text("@r\nT\n+\nI\n")
    out = tmp_path / "run"
    setup_pipeline(
        _ns(input_dir=str(reads), output_dir=str(out), dummy=[["dummy"]])
    )
    import pytest

    with pytest.raises(SystemExit):
        setup_pipeline(
            _ns(
                output_dir=str(out),
                add_annotator=True,
                dummy=[["dummy"]],
            )
        )


def test_abundance_dir_does_not_use_annotator_stems_as_samples(tmp_path):
    from samovar.table2iss import _sample_tables_from_abundance_dir

    folder = tmp_path / "abund"
    folder.mkdir()
    (folder / "kaiju.csv").write_text("taxid,N_1_full,N_2_full\n1,10,20\n")
    (folder / "kraken2.csv").write_text("taxid,N_1_full,N_2_full\n1,3,4\n")
    tables = _sample_tables_from_abundance_dir(
        folder, ["kaiju", "kraken2"]
    )
    assert sorted(tables) == ["1_full", "2_full"]
    assert "N_kaiju" in tables["1_full"].columns
    assert "N_kraken2" in tables["1_full"].columns


def test_annotator_fastqs_complete(tmp_path):
    r1 = tmp_path / "s1_kaiju_R1.fastq"
    r1.write_text("@a\nACGT\n+\nIIII\n")
    (tmp_path / "s1_kaiju_R2.fastq").write_text("@a\nACGT\n+\nIIII\n")
    assert annotator_fastqs_complete(str(tmp_path), ["s1"], "kaiju")
    assert not annotator_fastqs_complete(str(tmp_path), ["s1", "s2"], "kaiju")


def test_exec_runs_versioned_pipeline(tmp_path):
    log = tmp_path / ".log"
    log.mkdir()
    v1 = log / "samovar_v1.sh"
    marker = tmp_path / "ran_v1"
    v1.write_text(f"#!/bin/bash\necho ok > '{marker}'\n")
    v1.chmod(v1.stat().st_mode | stat.S_IEXEC)
    (log / "samovar.sh").write_text("#!/bin/bash\nexit 1\n")
    (log / "active_pipeline").write_text("samovar_v1.sh\n")
    result = subprocess.run(
        ["bash", str(REPO / "bin" / "samovar_exec"), "--output_dir", str(tmp_path)],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, "PYTHONPATH": str(REPO / "src")},
    )
    assert marker.is_file()
    assert "samovar_v1.sh" in result.stdout
