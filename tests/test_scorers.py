import json
from pathlib import Path

from samovar.config import PipelineConfig
from samovar.main_config import parse_tool_entry
from samovar.paths import write_config
from samovar.scorers import (
    DEFAULT_SCORING_INPUTS,
    expand_scoring_inputs,
    flags_apply_to_scorer,
    run_custom_scorers,
    run_scorer,
)
from samovar.tools_import import import_tool

COUNT_WRAPPER = Path(__file__).resolve().parent / "tools" / "count_annotations.py"


def test_expand_scoring_inputs_annotations_and_table(tmp_path):
    initial = tmp_path / "initial_annotations"
    regen = tmp_path / "regenerated_annotations"
    extra = tmp_path / "other_output"
    initial.mkdir()
    regen.mkdir()
    extra.mkdir()
    (initial / "combined_annotation_table.csv").write_text("taxid,n\n1,2\n")
    (regen / "a.annotation.csv").write_text("x\n")
    dirs = expand_scoring_inputs(tmp_path, "*annotations")
    assert {p.name for p in dirs} == {"initial_annotations", "regenerated_annotations"}
    tables = expand_scoring_inputs(tmp_path, "annotation_table")
    assert tables == [initial / "combined_annotation_table.csv"]
    both = expand_scoring_inputs(
        tmp_path, "*annotations,combined_annotation_table.csv"
    )
    names = {p.name for p in both}
    assert "initial_annotations" in names
    assert "combined_annotation_table.csv" in names


def test_custom_python_scorer_and_flags(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="counts",
        tool_type="scoring",
        exec_path=str(COUNT_WRAPPER),
        flags="--min-files 0",
        inputs="*annotations",
        also_repo_build=False,
    )
    ann = tmp_path / "out" / "initial_annotations"
    ann.mkdir(parents=True)
    (ann / "s.annotation.csv").write_text("a\tb\n")
    dest = run_scorer(
        "counts",
        run_dir=tmp_path / "out",
        config={"stage": "viz_initial", "scoring_flags": "--tag x"},
    )
    report = (dest / "file_counts.tsv").read_text()
    assert "initial_annotations" in report
    assert (dest / "stage.txt").read_text().strip() == "viz_initial"
    assert flags_apply_to_scorer("scoring", "counts")
    assert flags_apply_to_scorer("viz", "counts")
    assert flags_apply_to_scorer("counts", "counts")
    assert not flags_apply_to_scorer("annotator", "counts")


def test_run_custom_scorers_skips_unmatched_glob(tmp_path, monkeypatch):
    script = tmp_path / "noop.py"
    script.write_text(
        "def score(inputs, output_dir, config):\n"
        "    output_dir.mkdir(parents=True, exist_ok=True)\n"
        "    (output_dir / 'ran.txt').write_text('ok')\n"
    )
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="noop",
        tool_type="scoring",
        exec_path=str(script),
        inputs="missing_glob_zzz",
        also_repo_build=False,
    )
    written = run_custom_scorers(tmp_path / "empty")
    assert written == []


def test_cli_scorer_fallback(tmp_path, monkeypatch):
    script = tmp_path / "cli_score.py"
    script.write_text(
        "#!/usr/bin/env python3\n"
        "import argparse\n"
        "from pathlib import Path\n"
        "\n"
        "def main():\n"
        "    p = argparse.ArgumentParser()\n"
        "    p.add_argument('-i', action='append', dest='inputs', default=[])\n"
        "    p.add_argument('-o', dest='out')\n"
        "    p.add_argument('--mark', default='')\n"
        "    args = p.parse_args()\n"
        "    dest = Path(args.out)\n"
        "    dest.mkdir(parents=True, exist_ok=True)\n"
        "    dest.joinpath('cli.txt').write_text(args.mark + '|' + ','.join(args.inputs))\n"
        "\n"
        "if __name__ == '__main__':\n"
        "    main()\n"
    )
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="cli_score",
        tool_type="scoring",
        exec_path=str(script),
        flags="--mark hello",
        also_repo_build=False,
    )
    folder = tmp_path / "run" / "initial_annotations"
    folder.mkdir(parents=True)
    dest = run_scorer("cli_score", run_dir=tmp_path / "run")
    text = (dest / "cli.txt").read_text()
    assert "hello" in text
    assert "initial_annotations" in text


def test_prepare_scoring_flags_and_pipeline(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="counts",
        tool_type="scoring",
        exec_path=str(COUNT_WRAPPER),
        flags="--min-files 0",
        also_repo_build=False,
    )
    (tmp_path / "reads").mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(tmp_path / "reads"),
            "output_dir": str(tmp_path / "out"),
            "tool_flags": [
                ["scoring", "--tag group"],
                ["counts", "--min-files 1"],
            ],
            "scoring": None,
        },
    )()
    config = PipelineConfig.from_args(args)
    assert "--tag group" in (config.scoring_flags or "")
    assert "--min-files 1" in (config.scoring_tool_flags or {}).get("counts", "")
    paths = config.generate_configs(str(tmp_path / "out"))
    scoring_yaml = Path(paths["scoring"]).read_text()
    assert "--tag group" in scoring_yaml
    assert "counts" in scoring_yaml
    script = Path(config.generate_pipeline(str(tmp_path / "out"))).read_text()
    assert "samovar.scorers run" in script
    assert "viz_initial" in script
    raw = json.loads(cfg.read_text())["tools"]["counts"]
    spec = parse_tool_entry(raw, "counts")
    assert spec[5] == DEFAULT_SCORING_INPUTS
