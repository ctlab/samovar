"""OPAL-style profiling metrics and optional CAMI profile export."""

from pathlib import Path

import pandas as pd
import pytest

from samovar.opal import (
    confusion_rates,
    opal_command,
    opal_executable,
    opal_rank_range,
    opal_style_metrics,
    series_to_counts,
    write_cami_profile,
)
from samovar.scores import score_annotators
from samovar.viz_annotation import viz_annotation


def test_confusion_rates_tn_fn():
    true = ["A", "A", "A", "B", "B", "B"]
    pred = ["A", "A", "A", "A", "B", "B"]
    rates = confusion_rates(true, pred)
    # Class A: TP=3 FN=0 FP=1 TN=2 → TNR=2/3 FNR=0
    # Class B: TP=2 FN=1 FP=0 TN=3 → TNR=1 FNR=1/3
    assert rates["tnr"] == pytest.approx((2 / 3 + 1) / 2)
    assert rates["fnr"] == pytest.approx((0 + 1 / 3) / 2)
    assert rates["tn_pct"] == pytest.approx(100 * rates["tnr"])
    assert rates["fn_pct"] == pytest.approx(100 * rates["fnr"])


def test_opal_style_presence_and_l1():
    true = ["A", "A", "B", "B"]
    pred = ["A", "A", "A", "A"]
    stats = opal_style_metrics(true, pred)
    assert stats["opal_tp"] == 1
    assert stats["opal_fp"] == 0
    assert stats["opal_fn"] == 1
    assert stats["completeness"] == pytest.approx(0.5)
    assert stats["opal_purity"] == pytest.approx(1.0)
    assert stats["jaccard"] == pytest.approx(0.5)
    assert stats["l1_norm"] == pytest.approx(1.0)
    assert 0 <= stats["bray_curtis"] <= 1


def test_write_cami_profile(tmp_path):
    path = write_cami_profile({"562": 3, "9606": 1, "0": 2}, tmp_path / "gold.profile", "s1", "genus")
    text = path.read_text()
    assert "@SampleID:s1" in text
    assert "@Version:0.10.0" in text
    assert "562" in text and "9606" in text
    assert "\t0\t" not in text
    assert "PERCENTAGE" in text


def test_write_cami_profile_keeps_phage_without_genus(tmp_path):
    path = write_cami_profile({"2886930": 5, "9606": 1, "0": 2}, tmp_path / "gold.profile", "s1", "genus")
    text = path.read_text()
    assert "2886930" in text
    assert "9606" in text
    assert "\t0\t" not in text


def test_opal_rank_range_includes_species():
    assert opal_rank_range("genus") == "genus,species"
    assert opal_rank_range("species") == "species,species"


def test_series_to_counts_skips_unclassified():
    counts = series_to_counts(["562", "0", "562", "other"])
    assert counts == {"562": 2.0}


def test_opal_executable_optional():
    exe = opal_executable()
    assert exe is None or Path(exe).exists()


def test_opal_command_uses_python_for_script(tmp_path, monkeypatch):
    script = tmp_path / "opal.py"
    script.write_text("#!/usr/bin/env python3\n")
    monkeypatch.setattr("samovar.opal.opal_executable", lambda: str(script))
    cmd = opal_command()
    assert cmd is not None
    assert cmd[-1] == str(script)
    assert len(cmd) >= 1


def test_score_table_includes_opal_and_tn(tmp_path):
    df = pd.DataFrame(
        {
            "kaiju": ["562", "562", "9606", "9606"],
            "kraken2": ["562", "9606", "9606", "9606"],
            "true": ["562", "562", "9606", "9606"],
        }
    )
    table = score_annotators(df, ["kaiju", "kraken2"])
    for col in ("tnr", "fnr", "tn_pct", "fn_pct", "completeness", "opal_purity", "l1_norm"):
        assert col in table.columns
        assert table[col].notna().all()


def test_viz_writes_opal_scores_plot(tmp_path):
    df = pd.DataFrame(
        {
            "seq": [f"r{i}" for i in range(12)],
            "taxID_kaiju_0": [562] * 6 + [9606] * 6,
            "taxID_kraken2_1": [562] * 4 + [9606] * 8,
            "true": [562] * 6 + [9606] * 6,
        }
    )
    out = tmp_path / "plots"
    viz_annotation(
        df,
        type=("scores",),
        show_top=0,
        output_dir=str(out),
        plot=False,
        rank="none",
        use_names=False,
    )
    csv = pd.read_csv(out / "quality_scores.csv")
    assert "tn_pct" in csv.columns
    assert "completeness" in csv.columns
    assert (out / "scores.png").is_file()
    assert (out / "opal_scores.png").is_file()
