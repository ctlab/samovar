"""Stage JSON summaries and MultiQC custom-content staging."""

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
    (plots / "quality_scores.csv").write_text(
        "annotator,accuracy,f1,completeness\nkaiju,0.5,0.4,0.3\nkraken2,0.6,0.5,0.2\n"
    )
    write_stage_report(tmp_path, "viz_initial")
    staged = bundle_multiqc(tmp_path)
    assert (staged / "00_SamovaR_pipeline_mqc.html").is_file()
    html = (staged / "00_SamovaR_pipeline_mqc.html").read_text()
    assert "Initial quality plots" in html
    assert not list(staged.glob("*_mqc.png"))
    jsons = list(staged.glob("*_mqc.json"))
    htmls = [p for p in staged.glob("*_mqc.html") if "F1_kaiju" in p.name]
    assert jsons
    assert htmls
    chart = htmls[0].read_text()
    assert "vegaEmbed" in chart
    assert "F1 heatmap" in chart
    scores = next(p for p in jsons if "quality_scores" in p.name).read_text()
    assert "kaiju" in scores and "kraken2" in scores
