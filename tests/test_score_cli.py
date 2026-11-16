"""Public samovar score: annotators and profiles are separate invocations."""

import pandas as pd
import pytest

from samovar.viz_annotation import score_annotators_main, score_main, score_profiles_main


def test_score_annotators_writes_the_table(tmp_path):
    ann = tmp_path / "ann"
    ann.mkdir()
    pd.DataFrame({
        "seq": ["r1", "r2"],
        "taxID_kaiju_0": [562, 9606],
        "true": [562, 9606],
    }).to_csv(ann / "s1.csv", index=False)
    plots = tmp_path / "plots"
    assert score_main([
        "annotators", "--annotation-dir", str(ann), "--output-dir", str(plots), "--rank", "none",
    ]) == 0
    table = pd.read_csv(plots / "quality_scores.csv")
    assert "accuracy" in table.columns
    assert set(table["annotator"].astype(str)) == {"kaiju"}


def test_score_annotators_skips_a_missing_directory(tmp_path, capsys):
    assert score_annotators_main(["--annotation-dir", str(tmp_path / "missing")]) == 0
    assert "skipping visualization" in capsys.readouterr().err


def test_score_annotators_rejects_a_file(tmp_path):
    path = tmp_path / "not-a-dir"
    path.write_text("x", encoding="utf-8")
    with pytest.raises(SystemExit) as exc:
        score_annotators_main(["--annotation-dir", str(path)])
    assert exc.value.code != 0


def test_score_profiles_on_an_empty_stage_exits_zero(tmp_path):
    assert score_profiles_main(["--output-dir", str(tmp_path), "--stage", "initial"]) == 0
    assert not (tmp_path / "initial_annotations_plots" / "profile_quality_scores.csv").exists()


def test_score_profiles_writes_the_profile_table(tmp_path):
    reads = tmp_path / "initial"
    reads.mkdir()
    (reads / "A.fastq").write_text("@r1|taxid:562\nACGT\n+\nIIII\n", encoding="utf-8")
    ann = tmp_path / "initial_annotations"
    ann.mkdir()
    pd.DataFrame({"seq": ["r1"], "taxID_kraken2_0": [562]}).to_csv(ann / "A.csv", index=False)
    assert score_main(["profiles", "--output-dir", str(tmp_path), "--stage", "initial"]) == 0
    scores = pd.read_csv(tmp_path / "initial_annotations_plots" / "profile_quality_scores.csv")
    assert "kraken2" in set(scores["annotator"].astype(str))
    assert "A" in set(scores["sample"].astype(str))
