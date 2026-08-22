"""Tests for Python visualization and annotation I/O."""

import pandas as pd
from pathlib import Path

from samovar.annotation_io import annotation_to_abundance, read_annotation_dir
from samovar.viz_annotation import (
    _clip_taxon_label,
    compare_annotations,
    fpc_taxon_order,
    viz_annotation,
)


def test_read_annotation_dir_skips_combined(tmp_path):
    pd.DataFrame({"seq": ["r1"], "taxID_kaiju_0": [562], "true": [562]}).to_csv(
        tmp_path / "1.annotation.csv", index=False
    )
    pd.DataFrame({"seq": ["r1"], "taxID_kaiju_0": [1]}).to_csv(
        tmp_path / "combined_annotation_table.csv", index=False
    )
    df = read_annotation_dir(tmp_path)
    assert len(df) == 1
    assert df["sample"].iloc[0] == "1"
    assert "taxID_kaiju" in df.columns


def test_annotation_to_abundance_counts(tmp_path):
    df = pd.DataFrame(
        {
            "taxID_kaiju": [562, 562, 9606],
            "sample": ["1", "1", "1"],
        }
    )
    tables = annotation_to_abundance(df, n_reads=100)
    assert "taxID_kaiju" in tables
    out = tables["taxID_kaiju"]
    assert set(out["taxid"].astype(str)) == {"562", "9606"}
    assert int(out["N_1"].sum()) == 100


def test_clip_taxon_label():
    assert _clip_taxon_label("Escherichia") == "Escherichia"
    long = "Candidatus Blochmanniella extra"
    clipped = _clip_taxon_label(long, 25)
    assert len(clipped) == 25
    assert clipped.endswith("…")


def test_fpc_taxon_order_stable():
    freq = pd.DataFrame(
        {"x": ["1", "2", "1"], "y": ["1", "2", "2"], "Freq": [10, 5, 1]}
    )
    order = fpc_taxon_order(freq, "x", "y")
    assert set(order) == {"1", "2"}
    assert "0" not in order


def test_viz_annotation_writes_png(tmp_path):
    df = pd.DataFrame(
        {
            "seq": [f"r{i}" for i in range(20)],
            "taxID_kaiju_0": [562] * 10 + [9606] * 8 + [0, 0],
            "taxID_kraken2_1": [562] * 8 + [9606] * 10 + [0, 0],
            "true": [562] * 10 + [9606] * 10,
        }
    )
    out = tmp_path / "plots"
    results = viz_annotation(
        df,
        type=("f1", "R2", "cv"),
        show_top=0,
        output_dir=str(out),
        plot=False,
        rank="none",
        use_names=False,
    )
    assert "F1" in results
    assert "R2" in results
    assert "CV" in results
    pngs = list(out.glob("*.png"))
    assert pngs, "expected publication PNGs from cnsplots/matplotlib"
    assert any(p.name.startswith("F1_") for p in pngs)
    assert any(p.name.startswith("R2_") for p in pngs)
    assert any(p.name.startswith("CV_") for p in pngs)

    import matplotlib.image as mpimg

    f1_png = next(p for p in pngs if p.name.startswith("F1_"))
    img = mpimg.imread(f1_png)
    assert img.shape[0] > 50 and img.shape[1] > 50


def test_require_viz_backend():
    from samovar.viz_annotation import require_viz_backend

    backend = require_viz_backend()
    assert backend in {"cnsplots", "altair"}


def test_compare_annotations_cli_fails_on_missing_dir(tmp_path):
    import os
    import subprocess
    import sys

    from samovar.paths import repo_root

    script = repo_root() / "workflow" / "compare_annotations.py"
    env = os.environ.copy()
    src = str(repo_root() / "src")
    env["PYTHONPATH"] = src + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    proc = subprocess.run(
        [sys.executable, str(script), "--annotation_dir", str(tmp_path / "missing")],
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode != 0


def test_compare_annotations_fails_on_empty_dir(tmp_path):
    import pytest

    with pytest.raises(FileNotFoundError, match="per-sample"):
        compare_annotations(annotation_dir=str(tmp_path))


def test_compare_annotations_fails_on_combined_only(tmp_path):
    import pytest

    pd.DataFrame({"seq": ["r1"]}).to_csv(
        tmp_path / "combined_annotation_table.csv", index=False
    )
    with pytest.raises(FileNotFoundError, match="per-sample"):
        compare_annotations(annotation_dir=str(tmp_path))


def test_compare_annotations_cli_fails_on_empty_dir(tmp_path):
    import os
    import subprocess
    import sys

    from samovar.paths import repo_root

    script = repo_root() / "workflow" / "compare_annotations.py"
    env = os.environ.copy()
    src = str(repo_root() / "src")
    env["PYTHONPATH"] = src + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    proc = subprocess.run(
        [sys.executable, str(script), "--annotation_dir", str(tmp_path)],
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode != 0
    assert "per-sample" in (proc.stderr + proc.stdout)


def test_compare_annotations_cli_integrity(tmp_path):
    src = Path("data/test_annotations")
    out = tmp_path / "plots"
    csv = tmp_path / "combined.csv"
    data = compare_annotations(
        annotation_dir=str(src),
        output_dir=str(out),
        csv_file=str(csv),
        show_top=0,
        types=("f1", "R2", "cv"),
        rank="none",
    )
    assert not data.empty
    assert csv.exists()
    assert list(out.glob("F1_*.png"))
    # Combined table on disk is not rewritten by viz
    orig = pd.read_csv(src / "1.annotation.csv")
    assert "taxID_kraken2_0" in orig.columns or "taxID_kraken2" in list(map(str, orig.columns))
