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
    # PATH restore uses _SAMOVAR_CALLER_PATH; that name contains the letters R_PATH=.
    assert "export R_PATH=" not in script
    assert "workflow/ML.py" in script
    assert '"$SNAKEMAKE" -s ' in script
    assert '--directory "$out_dir"' in script
    assert "workflow/annotators/Snakefile" in script
    assert "samovar.exec_control" in script
    for step in ("setup_reads", "abundance_tables", "regenerate_tables", "score_regenerated_tables", "regenerate_reads", "reprofile"):
        assert f"ckpt_skip {step}" in script
    assert "visualization failed; continuing" not in script
    assert '[ -f "$CKPT/$1.done" ]' in script
