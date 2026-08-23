"""Stage JSON summaries and MultiQC custom-content staging."""

import json
from pathlib import Path

from samovar.exec_control import mark_done
from samovar.stage_report import bundle_multiqc, write_stage_report


def test_mark_done_writes_samovar_json(tmp_path):
    (tmp_path / "initial").mkdir()
    (tmp_path / "initial" / "1_full_R1.fastq").write_text("@r\nA\n+\nI\n")
    mark_done(tmp_path, "setup_reads")
    report = tmp_path / ".log" / "multiqc" / "setup_reads.samovar.json"
    assert report.is_file()
    text = report.read_text()
    assert '"stage": "setup_reads"' in text
    assert "Read setup" in text
    assert "samovar_version" in text


def test_bundle_multiqc_includes_plots_and_scores(tmp_path):
    plots = tmp_path / "initial_annotations_plots"
    plots.mkdir()
    (plots / "scores.png").write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 16)
    (plots / "F1_kaiju.html").write_text(
        '<html><body><div id="vis"></div>'
        "<script>var spec = {\"mark\": \"bar\"}; var embedOpt = {};</script>"
        "</body></html>"
    )
    (plots / "F1_kaiju_mqc.json").write_text(
        json.dumps(
            {
                "plot_type": "heatmap",
                "section_name": "F1 heatmap — kaiju",
                "data": [[1, 0], [0, 1]],
                "xcats": ["562", "9606"],
                "ycats": ["562", "9606"],
            }
        )
    )
    (plots / "quality_scores.csv").write_text(
        "annotator,accuracy,f1,completeness\nkaiju,0.5,0.4,0.3\nkraken2,0.6,0.5,0.2\nconsensus,0.55,0.45,0.4\n"
    )
    write_stage_report(tmp_path, "viz_initial")
    staged = bundle_multiqc(tmp_path)
    assert (staged / "00_run_options_mqc.html").is_file()
    html = (staged / "00_run_options_mqc.html").read_text()
    assert "SAMOVAR" in html or "Options used" in html
    assert not (staged / "00_SamovaR_pipeline_mqc.html").exists()
    assert not list(staged.glob("*altair*"))
    jsons = list(staged.glob("*_mqc.json"))
    assert jsons
    scores = json.loads(next(p for p in jsons if "quality_scores" in p.name).read_text())
    assert "kaiju" in scores["data"] and "kraken2" in scores["data"]
    assert "SAMOVAR" in scores["data"]
    assert "consensus" not in scores["data"]
    bars = json.loads(next(p for p in jsons if "score_bars" in p.name).read_text())
    assert bars["pconfig"]["stacking"] == "group"
    assert bars["pconfig"]["cpswitch"] is False
    assert bars["pconfig"]["sort_samples"] is False
    assert "Accuracy" in bars["data"]
    assert "kaiju" in bars["data"]["Accuracy"]
    cfg = (staged / "multiqc_config.yaml").read_text()
    assert "intro_text: SAMOVAR report" in cfg
    assert (staged / "99_conclusion_mqc.html").is_file()
    assert "Raw" in next(p for p in jsons if "F1_kaiju" in p.name).read_text()


def test_write_heatmap_mqc_is_native_multiqc(tmp_path):
    import pandas as pd

    from samovar.stage_report import write_heatmap_mqc

    mat = pd.DataFrame([[1, 2], [3, 4]], index=["a", "b"], columns=["x", "y"])
    path = tmp_path / "initial_annotations_plots" / "F1_kaiju_mqc.json"
    write_heatmap_mqc(
        mat,
        path,
        section_name="F1 heatmap — kaiju",
        description="exportable",
        xlab="pred",
        ylab="true",
    )
    payload = json.loads(path.read_text())
    assert payload["plot_type"] == "heatmap"
    assert payload["data"] == [[1, 2], [3, 4]]
    assert payload["xcats"] == ["x", "y"]
    assert payload["id"].startswith("samovar_raw_")
    assert payload["section_id"] == payload["id"]
    assert payload["pconfig"]["cluster_rows"] is False
    assert payload["pconfig"]["xcats_samples"] is False


def test_run_multiqc_invokes_real_cli(tmp_path, monkeypatch):
    import subprocess

    from samovar.stage_report import run_multiqc

    recorded = []

    def fake_call(cmd, *args, **kwargs):
        recorded.append(cmd)
        return 0

    monkeypatch.setattr("samovar.paths.discover_multiqc", lambda: "/usr/bin/multiqc")
    monkeypatch.setattr(subprocess, "call", fake_call)
    write_stage_report(tmp_path, "setup_reads")
    assert run_multiqc(tmp_path) == 0
    assert recorded
    cmd = recorded[0]
    assert cmd[0] == "/usr/bin/multiqc"
    assert "--interactive" in cmd
    assert "--export" not in cmd
    assert run_multiqc(tmp_path, extra_args=["--", "--export"]) == 0
    assert "--export" in recorded[-1]


def test_bundle_keeps_stage_plot_ids_unique(tmp_path):
    import pandas as pd

    from samovar.stage_report import write_heatmap_mqc

    mat = pd.DataFrame([[1.0, 0.0], [0.0, 2.0]], index=["562", "9606"], columns=["562", "9606"])
    for folder in ("initial_annotations_plots", "regenerated_annotations_plots"):
        path = tmp_path / folder / "F1_kaiju_mqc.json"
        write_heatmap_mqc(mat, path, section_name="F1 heatmap — kaiju", description="d", xlab="pred", ylab="true")
    write_stage_report(tmp_path, "viz_initial")
    write_stage_report(tmp_path, "viz_regenerated")
    staged = bundle_multiqc(tmp_path)
    payloads = [json.loads(p.read_text()) for p in staged.glob("*F1_kaiju_mqc.json")]
    ids = [p["id"] for p in payloads]
    assert len(ids) == 2
    assert len(set(ids)) == 2
    assert any(i.startswith("samovar_raw_") for i in ids)
    assert any(i.startswith("samovar_regenerated_") for i in ids)


def test_multiqc_render_is_adequate(tmp_path):
    import os
    import subprocess

    import pandas as pd
    import pytest

    from samovar.paths import discover_multiqc
    from samovar.stage_report import bundle_multiqc
    from samovar.viz_annotation import viz_annotation

    exe = discover_multiqc()
    if not exe:
        pytest.skip("MultiQC CLI is not installed")

    df = pd.DataFrame(
        {
            "seq": [f"r{i}" for i in range(20)],
            "taxID_kaiju_0": [562] * 10 + [9606] * 8 + [0, 0],
            "taxID_kraken2_1": [562] * 8 + [9606] * 10 + [0, 0],
            "true": [562] * 10 + [9606] * 10,
        }
    )
    for folder in ("initial_annotations_plots", "regenerated_annotations_plots"):
        viz_annotation(
            df,
            type=("f1", "R2", "cv", "scores"),
            show_top=0,
            output_dir=str(tmp_path / folder),
            plot=False,
            rank="none",
            use_names=False,
        )
    write_stage_report(tmp_path, "viz_initial")
    write_stage_report(tmp_path, "viz_regenerated")
    staged = bundle_multiqc(tmp_path)
    out = tmp_path / "multiqc"
    proc = subprocess.run(
        [exe, str(staged), "-o", str(out), "-f", "-c", str(staged / "multiqc_config.yaml"), "--interactive"],
        capture_output=True,
        text=True,
        env={**os.environ, "MPLBACKEND": "Agg"},
        timeout=120,
    )
    log = proc.stdout + "\n" + proc.stderr
    assert proc.returncode == 0, log
    for needle in ("Error parsing", "Error creating", "not recognised"):
        assert needle not in log, log
    htmls = list(out.glob("*multiqc_report.html"))
    assert htmls, f"no report in {list(out.iterdir())}\n{log}"
    html_path = htmls[0]
    html = html_path.read_text(encoding="utf-8", errors="replace")
    assert "SamovaR" in html or "SAMOVAR" in html
    assert "Raw" in html
    assert "Regenerated" in html
    assert "vegaEmbed" not in html
    data_dirs = list(out.glob("*multiqc_report_data"))
    assert data_dirs
    plot_data = json.loads((data_dirs[0] / "multiqc_data.json").read_text())["report_plot_data"]
    keys = " ".join(plot_data)
    assert "samovar_raw_F1_" in keys
    assert "samovar_regenerated_F1_" in keys
    assert "samovar_raw_R2_" in keys
    assert "samovar_regenerated_R2_" in keys
    assert "samovar_raw_score_bars" in keys
    assert html_path.stat().st_size > 50_000


def test_bundle_hides_hybrid_platform_sections(tmp_path):
    plots = tmp_path / "initial_annotations_plots"
    plots.mkdir()
    (plots / "F1_kaiju_mqc.json").write_text(
        json.dumps({"plot_type": "heatmap", "section_name": "F1 heatmap — kaiju", "data": [[1]], "xcats": ["a"], "ycats": ["a"]})
    )
    (plots / "F1_kaiju.illumina_mqc.json").write_text(
        json.dumps({"plot_type": "heatmap", "section_name": "F1 heatmap — kaiju", "data": [[1]], "xcats": ["a"], "ycats": ["a"]})
    )
    (plots / "quality_scores.csv").write_text(
        "annotator,accuracy,f1,read_type\n"
        "kaiju,0.5,0.4,all\nkraken2,0.6,0.5,all\n"
        "kaiju,0.2,0.1,illumina\nkraken2,0.3,0.2,illumina\n"
        "kaiju,0.9,0.8,ont\nkraken2,0.1,0.1,ont\n"
    )
    staged = bundle_multiqc(tmp_path)
    payloads = [json.loads(p.read_text()) for p in staged.glob("*_mqc.json")]
    parents = {p.get("parent_id") for p in payloads}
    assert "samovar_raw" in parents
    assert "samovar_raw_illumina" in parents
    assert "samovar_raw_ont" in parents
    css = (staged / "samovar_multiqc.css").read_text()
    assert "samovar_raw_illumina" in css
    assert "samovar-show-platforms" in css
    html = (staged / "00_run_options_mqc.html").read_text()
    assert "samovar-toggle-platforms" in html


def test_bundle_unknown_truth_omits_nan_metrics(tmp_path):
    plots = tmp_path / "initial_annotations_plots"
    plots.mkdir()
    (plots / "quality_scores.csv").write_text(
        "annotator,n_reads,n_taxa,accuracy,f1,r2\nkaiju,10,2,,,\nkraken2,10,3,,,\nSAMOVAR,10,2,,,\n"
    )
    (plots / "CV_kaiju_vs_kraken2_mqc.json").write_text(
        json.dumps({"plot_type": "heatmap", "section_name": "Cross-validation", "data": [[1]], "xcats": ["a"], "ycats": ["a"]})
    )
    staged = bundle_multiqc(tmp_path)
    tables = [json.loads(p.read_text()) for p in staged.glob("*quality_scores_mqc.json")]
    assert tables
    data = tables[0]["data"]
    assert "kaiju" in data
    assert "f1" not in data["kaiju"]
    assert "n_reads" in data["kaiju"]
    conclusion = (staged / "99_conclusion_mqc.html").read_text()
    assert "SAMOVAR" in conclusion
    bars = list(staged.glob("*score_bars_mqc.json"))
    if bars:
        payload = json.loads(bars[0].read_text())
        assert "F1 (current)" not in payload.get("data", {})

