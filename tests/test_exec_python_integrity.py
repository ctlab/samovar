"""Integrity: generated exec pipeline uses Python viz, not R."""

import argparse
from pathlib import Path

from samovar.config import setup_pipeline


def test_generated_exec_pipeline_is_python_only(tmp_path):
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(tmp_path / "reads"),
        output_dir=str(tmp_path / "out"),
        cores=1,
        max_genomes=5,
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=[["kaiju /tmp/kaiju"]],
    )
    (tmp_path / "reads").mkdir()
    result = setup_pipeline(args)
    script = Path(result["pipeline"]).read_text()
    assert "compare_annotations.py" in script
    assert "compare_annotations.R" not in script
    assert "R_PATH" not in script
    assert "workflow/ML.py" in script
    assert "snakemake -s " in script
    assert "workflow/annotators/Snakefile" in script
