"""Runtime is a measured elapsed time, not a quality score and not a filled-in zero."""

import json

from samovar.config import _checkpoint_block
from samovar.multiqc_report import effective_settings
from samovar.runtime import load_sample_runtimes, mean_per_run, record_stage
from samovar.stage_report import bundle_multiqc


def _benchmark(path, seconds, rss=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if rss is None:
        path.write_text(f"s\th:m:s\n{seconds}\t0:00:01\n", encoding="utf-8")
    else:
        path.write_text(
            f"s\th:m:s\tmax_rss\n{seconds}\t0:00:01\t{rss}\n",
            encoding="utf-8",
        )


def test_checkpoint_records_only_when_the_stage_runs():
    text = _checkpoint_block("viz_initial", "echo hi")
    assert 'record-stage "$out_dir" viz_initial' in text
    assert text.index("record-stage") < text.index("ckpt_finish viz_initial")
    skipped = text.split("else", 1)[0]
    assert "record-stage" not in skipped


def test_per_sample_mean_is_not_the_stage_wall_clock(tmp_path):
    reports = tmp_path / "initial_reports" / "runtime" / "samples" / "kraken2_std" / "kraken2"
    _benchmark(reports / "a.tsv", 3, rss=100)
    _benchmark(reports / "b.tsv", 4, rss=300)
    broken = reports / "c.tsv"
    broken.write_text("max_rss\n10\n", encoding="utf-8")
    record_stage(str(tmp_path), "annotate_initial", "10", "15")
    rows = load_sample_runtimes(tmp_path)
    assert [row["sample"] for row in rows] == ["a", "b"]
    means = mean_per_run(rows)[0]
    assert abs(means["mean_elapsed_s"] - 3.5) < 1e-9
    assert abs(means["mean_max_rss_mb"] - 200) < 1e-9
    staged = bundle_multiqc(tmp_path, effective_settings([{"sections": ["runtime"]}]))
    box = json.loads((staged / "runtime_box_mqc.json").read_text())
    assert box["plot_type"] == "box"
    assert box["data"]["initial: kraken2_std"] == [3.0, 4.0]
    memory = json.loads((staged / "runtime_memory_box_mqc.json").read_text())
    assert memory["plot_type"] == "box"
    assert memory["data"]["initial: kraken2_std"] == [100.0, 300.0]
    summary = json.loads((staged / "runtime_summary_mqc.json").read_text())
    by_kind = {row["kind"]: row for row in summary["data"].values()}
    assert by_kind["per_sample_mean"]["elapsed_s"] == 3.5
    assert by_kind["per_sample_mean"]["max_rss_mb"] == 200.0
    assert by_kind["absolute"]["elapsed_s"] == 5.0
    assert by_kind["absolute"]["max_rss_mb"] == ""
    hidden = bundle_multiqc(tmp_path, effective_settings([{"sections": ["profiles"], "exclude": ["runtime"]}]))
    assert not list(hidden.glob("runtime*"))
