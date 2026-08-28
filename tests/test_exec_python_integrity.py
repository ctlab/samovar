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
    assert "export R_PATH=" not in script
    assert "R_PATH=" not in script.replace("SAMOVAR_PATH", "")
    assert "workflow/ML.py" in script
    assert "snakemake -s " in script
    assert "workflow/annotators/Snakefile" in script
    assert "samovar.exec_control" in script
    for step in ("setup_reads", "abundance_tables", "regenerate_tables", "score_regenerated_tables", "regenerate_reads", "reprofile"):
        assert f"ckpt_skip {step}" in script
    assert "visualization failed; continuing" not in script
    assert '[ -f "$CKPT/$1.done" ]' in script


def test_workflow_scripts_are_python_or_cpp():
    """Former R workflow steps must be Python or the C++ combiner."""
    from samovar.paths import repo_root

    root = repo_root()
    workflow = root / "workflow"
    r_files = list(workflow.rglob("*.R")) + list(workflow.rglob("*.Rmd"))
    assert r_files == [], f"R workflow leftovers: {r_files}"
    assert (workflow / "compare_annotations.py").is_file()
    assert (workflow / "annotation_regenerate.py").is_file()
    assert (workflow / "combine_annotation_tables.py").is_file()
    assert (workflow / "remap_taxids.py").is_file()
    assert (workflow / "ML.py").is_file()
    cpp = root / "src" / "cpp" / "combine_annotations.cpp"
    assert cpp.is_file()
    combine = (workflow / "combine_annotation_tables.py").read_text(encoding="utf-8")
    assert "samovar_combine_annotations" in combine or "combine_with_cpp" in combine


def test_ensemble_pipeline_is_self_contained_without_r(tmp_path):
    """prepare/exec configs must not require R sources from this tree."""
    from samovar.paths import repo_root

    root = repo_root()
    assert not (root / "DESCRIPTION").exists()
    assert not (root / "R").exists()
    assert not (root / "workflow" / "compare_annotations.R").exists()
    assert (root / "workflow" / "compare_annotations.py").is_file()
    assert (root / "workflow" / "ML.py").is_file()
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
    assert ".R" not in script or "compare_annotations.R" not in script
