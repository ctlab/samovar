"""Algorithmic methods.md from prepare configs (no LLM prose)."""

from pathlib import Path

import pytest
import yaml

from samovar.config import AnnotatorConfig, PipelineConfig
from samovar.exec_control import CHECKPOINT_STEPS
from samovar.methods import (
    collect_stage_rows,
    join_en,
    latex_escape,
    load_template,
    main,
    methods_template_dir,
    parse_used_citation_keys,
    render_methods_md,
    render_methods_tex,
    render_methods_tsv,
    write_methods,
)
from samovar.tool_contracts import CONTRACTS


def _dump_yaml(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data), encoding="utf-8")


def _toy_run(tmp_path: Path, *, start="setup_reads", end="viz_reprofiled") -> Path:
    configs = tmp_path / ".log" / "configs"
    configs.mkdir(parents=True)
    (tmp_path / ".log" / "window.env").write_text(
        f'export SAMOVAR_START={start}\nexport SAMOVAR_END={end}\n',
        encoding="utf-8",
    )
    _dump_yaml(
        configs / "config_init.yaml",
        {
            "run_config": [
                {
                    "run_name": "k2",
                    "type": "kraken2",
                    "cmd": "kraken2",
                    "db_path": "/db/std8",
                    "db_name": "std8",
                },
                {
                    "run_name": "kj",
                    "type": "kaiju",
                    "cmd": "kaiju",
                    "db_path": "/db/nr",
                    "db_name": "nr",
                },
            ]
        },
    )
    _dump_yaml(
        configs / "config_annotation2iss.yaml",
        {
            "reads_generator": "iss",
            "table_reads_generators": ["direct"],
            "regeneration_mode": "direct",
            "table_score": "",
            "sample_score": "",
            "sample_filter": "",
        },
    )
    _dump_yaml(configs / "config_qc.yaml", {"qc": "", "qc_initial": "off", "qc_generated": "off"})
    _dump_yaml(configs / "config_scoring.yaml", {"scoring_tools": []})
    _dump_yaml(
        configs / "config_reprofiling.yaml",
        {"reprofiler": "ensemble", "feature_importance": "off"},
    )
    _dump_yaml(configs / "config_export.yaml", {"export": "off", "export_formats": []})
    (tmp_path / "used_citations.bib").write_text(
        "% --- kraken2 ---\n"
        "@article{Wood_2019, title={Improved metagenomic analysis with Kraken 2}, year={2019}}\n\n"
        "% --- kaiju ---\n"
        "@article{Menzel_2016, title={Fast and sensitive taxonomic classification}, year={2016}}\n\n"
        "% --- iss ---\n"
        "@article{Gourle_2019, title={Simulating Illumina metagenomic data}, year={2019}}\n",
        encoding="utf-8",
    )
    return tmp_path


def test_join_en_and_latex_escape():
    assert join_en([]) == ""
    assert join_en(["a"]) == "a"
    assert join_en(["a", "b"]) == "a and b"
    assert join_en(["a", "b", "c"]) == "a, b, and c"
    assert r"\_" in latex_escape("foo_bar")
    assert r"\&" in latex_escape("a & b")


def test_parse_used_citation_keys_by_tool_section():
    text = (
        "% --- kraken2 ---\n"
        "@article{Wood_2019, title={x}}\n"
        "% --- kaiju ---\n"
        "@article{Menzel_2016, title={y}}\n"
        "@article{Extra_2016, title={z}}\n"
    )
    mapping = parse_used_citation_keys(text)
    assert mapping["kraken2"] == ["Wood_2019"]
    assert mapping["kaiju"] == ["Menzel_2016", "Extra_2016"]


def test_every_contract_has_a_methods_sentence_template():
    dest = methods_template_dir()
    missing = []
    for name in CONTRACTS:
        tmpl = load_template(name, dest)
        if not tmpl.get("md") or not tmpl.get("tex"):
            missing.append(name)
    assert not missing, f"cite/methods lacks md/tex for {missing}"
    annot = load_template("annotator")
    assert "{items}" in annot["md"]
    assert "{cite_tex}" in annot["tex"]


def test_collect_stage_rows_fills_annotator_ensemble_sentence(tmp_path):
    out = _toy_run(tmp_path)
    rows = collect_stage_rows(out)
    annot = [r for r in rows if r["stage"] == "annotate_initial" and r["contract"] == "annotator"]
    assert len(annot) == 1
    sentence = annot[0]["md"]
    assert "kraken2 on the std8 database" in sentence
    assert "kaiju on the nr database" in sentence
    assert "were used" in sentence
    assert "@Wood_2019" in sentence
    assert "@Menzel_2016" in sentence
    assert annot[0]["keys"] == ["Wood_2019", "Menzel_2016"]
    qc_rows = [r for r in rows if r["contract"] == "qc" and r["md"]]
    assert qc_rows == []
    reads = next(r for r in rows if r["stage"] == "setup_reads" and r["contract"] == "reads_generator")
    assert "iss" in reads["md"]
    assert "@Gourle_2019" in reads["md"]
    workflow = [r for r in rows if r["contract"] == "workflow"]
    assert workflow
    stages = {r["stage"] for r in rows}
    assert stages <= set(CHECKPOINT_STEPS)
    assert "setup_reads" in stages
    assert "reprofile" in stages


def test_window_env_limits_stages(tmp_path):
    out = _toy_run(tmp_path, start="annotate_initial", end="combine_initial")
    rows = collect_stage_rows(out)
    stages = [r["stage"] for r in rows]
    assert stages == ["annotate_initial", "combine_initial"]
    assert all(r["stage"] != "setup_reads" for r in rows)


def test_render_md_tex_tsv_structure(tmp_path):
    rows = collect_stage_rows(_toy_run(tmp_path))
    md = render_methods_md(rows)
    assert md.startswith("# Methods\n")
    assert "## Stages" in md
    assert "### `annotate_initial` (annotator)" in md
    assert "**Citation keys:**" in md
    tex = render_methods_tex(rows)
    assert r"\section{Methods}" in tex
    assert r"\cite{Wood_2019,Menzel_2016}" in tex
    tsv = render_methods_tsv(rows)
    lines = tsv.strip().splitlines()
    assert lines[0] == "contract\tmethods\tcitation"
    annot_line = next(line for line in lines if line.startswith("annotator\t"))
    cols = annot_line.split("\t")
    assert len(cols) == 3
    assert "kraken2 on the std8 database" in cols[1]
    assert "Wood_2019" in cols[2]
    assert "kaiju" in cols[1]


def test_write_methods_and_cli(tmp_path):
    out = _toy_run(tmp_path)
    paths = write_methods(out)
    for name in ("methods.md", "samovar-method.md", "methods.tex", "methods.tsv", "methods.json"):
        assert Path(paths[name]).is_file()
        assert (out / ".log" / "configs" / name).is_file()
    assert (out / "methods.md").read_text(encoding="utf-8") == (
        out / "samovar-method.md"
    ).read_text(encoding="utf-8")
    assert main(["--output_dir", str(out)]) == 0
    assert main(["--output_dir", str(tmp_path / "missing")]) == 1


def test_generate_configs_writes_methods(tmp_path):
    cfg = PipelineConfig()
    cfg.output_dir = str(tmp_path)
    cfg.annotators = [
        AnnotatorConfig(
            run_name="k2", type="kraken2", db_path="/db/std8", cmd="kraken2", db_name="std8"
        ),
        AnnotatorConfig(
            run_name="kj", type="kaiju", db_path="/db/nr", cmd="kaiju", db_name="nr"
        ),
    ]
    cfg.reads_generator = "iss"
    cfg.regeneration_modes = ["direct"]
    cfg.regeneration_mode = "direct"
    cfg.run_multiqc = False
    cfg.export_corrector = "off"
    cfg.scoring_tools = []
    cfg.feature_importance = "off"
    cfg.qc_initial = ""
    cfg.startpoint = "setup_reads"
    cfg.endpoint = "annotate_initial"
    with pytest.warns(UserWarning, match="No BibTeX citation"):
        configs = cfg.generate_configs(str(tmp_path))
    md = Path(configs["methods"])
    assert md.is_file()
    text = md.read_text(encoding="utf-8")
    assert "kraken2 on the std8 database" in text
    assert "kaiju on the nr database" in text
    assert (tmp_path / "methods.tsv").is_file()
    tsv = (tmp_path / "methods.tsv").read_text(encoding="utf-8")
    assert tsv.splitlines()[0] == "contract\tmethods\tcitation"
    window = Path(configs["window"]).read_text(encoding="utf-8")
    assert "SAMOVAR_END=annotate_initial" in window
    rows = collect_stage_rows(tmp_path)
    assert [r["stage"] for r in rows if r["contract"] != "workflow"] == [
        "setup_reads",
        "qc_initial",
        "annotate_initial",
    ]
