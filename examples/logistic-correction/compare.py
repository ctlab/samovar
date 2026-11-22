#!/usr/bin/env python3
"""Score raw vs logistic-corrected abundances against ``true`` taxIDs.

Used by ``examples/logistic-correction/pipeline.sh``. Works on any completed
SAMOVAR run that has ``initial_annotations`` and ``regenerated_annotations``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from samovar.abundance import (
    input_to_abundance_tables,
    n_sample_columns,
    normalize_abundance_table,
)
from samovar.abundance_correctors import (
    apply_taxon_efficiency,
    as_annotation,
    fit_taxon_efficiency,
    maybe_attach_samovar_reference,
)
from samovar.parse_annotators import Annotation


def _profile_metrics(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """Mean canonical profile metrics over shared samples."""
    from samovar.profiles import metadata
    from samovar.scores import profile_scores

    pred = normalize_abundance_table(a)
    truth = normalize_abundance_table(b)
    cols = [col for col in n_sample_columns(pred) if col in n_sample_columns(truth)]
    if not cols:
        return {"l1_norm": float("nan"), "r2": float("nan")}
    scores = profile_scores(
        pred[["taxid", *cols]],
        truth[["taxid", *cols]],
        metadata(),
        metadata(),
    )
    return {name: float(scores[name].mean()) for name in ("l1_norm", "r2")}


def _true_table(annotation: Annotation) -> Optional[pd.DataFrame]:
    frame = annotation.DataFrame
    if frame is None or frame.empty or "true" not in frame.columns:
        return None
    work = frame.copy()
    if "sample" not in work.columns:
        work["sample"] = "1"
    piece = work[["true", "sample"]].copy()
    piece["true"] = piece["true"].astype(str)
    mat = piece.groupby(["true", "sample"], dropna=False).size().unstack(fill_value=0)
    out = mat.reset_index().rename(columns={"true": "taxid"})
    out.columns = ["taxid"] + [f"N_{c}" for c in out.columns[1:]]
    return normalize_abundance_table(out)


def _row(name, table, corrected, truth):
    work = normalize_abundance_table(table)
    cols = n_sample_columns(work)
    row = {"annotator": name, "raw_reads": float(work[cols].sum().sum()) if cols else 0.0}
    if truth is not None:
        raw = _profile_metrics(table, truth)
        row["L1_raw"] = raw["l1_norm"]
        row["R2_raw"] = raw["r2"]
        if name in corrected:
            logistic = _profile_metrics(corrected[name], truth)
            row["L1_logistic"] = logistic["l1_norm"]
            row["R2_logistic"] = logistic["r2"]
            row["delta_L1"] = row["L1_raw"] - row["L1_logistic"]
            row["delta_R2"] = row["R2_logistic"] - row["R2_raw"]
    return row


def compare_run(run_dir: Path, dest: Path) -> pd.DataFrame:
    initial = as_annotation(run_dir / "initial_annotations")
    regenerated = as_annotation(run_dir / "regenerated_annotations")
    truth = _true_table(initial)
    if truth is None:
        truth = _true_table(regenerated)
    model = run_dir / "reprofiled_annotations" / "trained_model.joblib"
    ref = maybe_attach_samovar_reference(regenerated, model if model.is_file() else None)
    rates = fit_taxon_efficiency(ref, C=1.0, min_efficiency=0.05)
    raw = input_to_abundance_tables(initial)
    corrected = apply_taxon_efficiency(raw, rates, min_efficiency=0.05)
    rows = [_row(name, table, corrected, truth) for name, table in raw.items()]
    repro = run_dir / "reprofiled_annotations"
    if repro.is_dir():
        ensemble = as_annotation(repro)
        frame = ensemble.DataFrame
        if frame is not None and not frame.empty and "sample" in frame.columns:
            frame = frame.copy()
            frame["sample"] = (
                frame["sample"].astype(str).str.replace(r"_reprofiled$", "", regex=True)
            )
            ensemble = Annotation.from_long_table(frame)
        ens_tables = input_to_abundance_tables(ensemble)
        ens_corr = apply_taxon_efficiency(ens_tables, rates, min_efficiency=0.05)
        for name, table in ens_tables.items():
            if str(name).lower() == "samovar":
                rows.append(_row(name, table, ens_corr, truth))
                raw[name] = table
                if name in ens_corr:
                    corrected[name] = ens_corr[name]
    summary = pd.DataFrame(rows)
    dest.mkdir(parents=True, exist_ok=True)
    summary.to_csv(dest / "comparison.csv", index=False)
    _plot(
        summary,
        dest / "comparison_L1.png",
        "L1_raw",
        "L1_logistic",
        "L1 (relative abundance vs true)",
        "Abundance L1 vs true (classified taxids)",
    )
    _plot(
        summary,
        dest / "comparison_R2.png",
        "R2_raw",
        "R2_logistic",
        "R² (relative abundance vs true)",
        "Abundance R² vs true (classified taxids)",
    )
    return summary


def _plot(
    summary: pd.DataFrame,
    dest: Path,
    raw_col: str,
    logi_col: str,
    ylabel: str,
    title: str,
) -> None:
    if summary.empty or raw_col not in summary.columns:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = summary["annotator"].astype(str).tolist()
    raw = summary[raw_col].tolist()
    logi = summary[logi_col].tolist() if logi_col in summary.columns else [0] * len(raw)
    x = range(len(names))
    width = 0.35
    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(names)), 4.2))
    ax.bar([i - width / 2 for i in x], raw, width, label="raw counts", color="#fdae6b")
    ax.bar([i + width / 2 for i in x], logi, width, label="logistic correction", color="#3182bd")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=30, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, help="SAMOVAR output directory")
    parser.add_argument("-o", "--output", default="", help="Figures directory (default: RUN/figures)")
    args = parser.parse_args(argv)
    run = Path(args.run)
    dest = Path(args.output) if args.output else run / "figures"
    summary = compare_run(run, dest)
    print(summary.to_string(index=False))
    print(f"Wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
