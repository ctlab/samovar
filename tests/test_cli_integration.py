"""Fast integration tests for prepare / build / exec wiring (no real annotators)."""

import os
import stat
import subprocess
from pathlib import Path
from unittest.mock import patch

import yaml

from samovar.build_database import build_database_from_config
from samovar.config import PipelineConfig, setup_pipeline


REPO = Path(__file__).resolve().parents[1]


def _ensure_config_json():
    cfg = REPO / "build" / "config.json"
    if cfg.exists():
        return
    cfg.parent.mkdir(parents=True, exist_ok=True)
    example = REPO / "build" / "config.json.example"
    if example.exists():
        cfg.write_text(example.read_text())
    else:
        cfg.write_text(
            '{"python_path": "python3", "r_path": "R", "r_lib_path": "."}\n'
        )


def test_cli_help_lists_prepare_build_exec():
    _ensure_config_json()
    result = subprocess.run(
        ["bash", str(REPO / "bin" / "samovar"), "help"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    )
    out = result.stdout
    assert "prepare" in out
    assert "build" in out
    assert "exec" in out
    assert "--redo" in out
    assert "--cleanup-tmp" in out


def test_prepare_writes_pipeline_and_configs(tmp_path):
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(tmp_path / "reads"),
            "output_dir": str(tmp_path / "out"),
            "kraken2": [["kraken2 /tmp/kraken_db"]],
            "kaiju": [["kaiju /tmp/kaiju_db/kaiju_db.fmi"]],
        },
    )()
    (tmp_path / "reads").mkdir()
    result = setup_pipeline(args)
    pipeline = Path(result["pipeline"])
    assert pipeline.exists()
    text = pipeline.read_text()
    assert "snakemake -s " in text
    assert "workflow/annotators/Snakefile" in text
    assert "workflow/combine_annotation_tables.py" in text
    assert "workflow/ML.py" in text
    assert str(tmp_path / "reads") in text
    init_cfg = yaml.safe_load(Path(result["configs"]["init_annotator"]).read_text())
    types = {run["type"] for run in init_cfg["run_config"]}
    assert types == {"kraken2", "kaiju"}


def test_build_from_config_invokes_kaiju_indexers(tmp_path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "562.fna").write_text(">seq\n" + "ATGGAATTCGGT" * 8 + "TAA\n")
    cfg = tmp_path / "db.yaml"
    cfg.write_text(
        yaml.dump({"input_dir": [str(genomes)], "output_dir": str(tmp_path / "prep")})
    )
    db_path = tmp_path / "kaiju_db"

    def fake_run(cmd, check=True, text=True):
        class Result:
            returncode = 0

        (db_path / "nodes.dmp").write_text("1\t|\t1\t|\n")
        (db_path / "kaiju_db.fmi").write_text("fake")
        return Result()

    with patch("samovar.build_database.run_command", side_effect=fake_run), patch(
        "samovar.build_database.get_taxonomy_db"
    ) as get_tax, patch(
        "samovar.build_database.fetch_proteome", return_value=None
    ), patch(
        "samovar.build_database.fetch_gff", return_value=None
    ):
        get_tax.side_effect = lambda db_path="kraken_db", taxonomy_path=None: (
            Path(db_path).mkdir(parents=True, exist_ok=True)
        )
        build_database_from_config(str(cfg), db_type="kaiju", db_path=str(db_path))

    library = db_path / "library.faa"
    library_gz = db_path / "library.faa.gz"
    if library_gz.exists() and not library.exists():
        import gzip as _gzip
        content = _gzip.open(library_gz, "rt").read()
    else:
        assert library.exists()
        content = library.read_text()
    assert "_562" in content or content.startswith(">562")
    assert "*" not in content.replace(">\n", "")


def test_exec_runs_generated_samovar_sh(tmp_path):
    log = tmp_path / ".log"
    log.mkdir()
    script = log / "samovar.sh"
    marker = tmp_path / "ran"
    script.write_text(f"#!/bin/bash\necho ok > '{marker}'\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    result = subprocess.run(
        ["bash", str(REPO / "bin" / "samovar_exec"), "--output_dir", str(tmp_path)],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    )
    assert marker.exists()
    assert "Pipeline execution completed" in result.stdout


def test_prepare_to_exec_mocked_end_to_end(tmp_path):
    """prepare generates a script; exec runs it. External tools are not required."""
    out_dir = tmp_path / "run"
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    (reads / "1_full_R2.fastq").write_text("@r\nT\n+\nI\n")
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(reads),
            "output_dir": str(out_dir),
            "kraken2": None,
            "kaiju": [["kaiju /tmp/kaiju_db.fmi"]],
        },
    )()
    result = setup_pipeline(args)
    pipeline = Path(result["pipeline"])
    # Replace the generated pipeline with a stub that still consumes the same path.
    pipeline.write_text("#!/bin/bash\nset -e\necho prepared-and-execed\n")
    proc = subprocess.run(
        ["bash", str(REPO / "bin" / "samovar_exec"), "--output_dir", str(out_dir)],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "prepared-and-execed" in proc.stdout or "Pipeline execution completed" in proc.stdout
    assert "Running main SamovaR pipeline" in proc.stdout
