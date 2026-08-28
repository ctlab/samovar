#!/usr/bin/env python3
"""Linear (logistic) reprofiler for annotation tables with known truth.

Import::

    samovar tools import -n linear \\
        --exec-path examples/reprofiling/linear_classifier.py --type ml \\
        --flags "--max-iter 500"

Then ``samovar prepare --reprofiler linear``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
from sklearn.linear_model import LogisticRegression

from samovar.reprofiling import preprocess_data, save_model
from samovar.reprofilers import ReprofileResult, profile_tables_with_model


def _priors_from_ground_truth(tables: Dict[str, pd.DataFrame]) -> Optional[Dict[int, float]]:
    counts: Dict[int, float] = {}
    total = 0.0
    for table in tables.values():
        if table is None or table.empty or "taxid" not in table.columns:
            continue
        n_cols = [c for c in table.columns if str(c).startswith("N_")]
        if not n_cols:
            continue
        for _, row in table.iterrows():
            try:
                taxid = int(row["taxid"])
            except (TypeError, ValueError):
                continue
            n = float(pd.to_numeric(row[n_cols], errors="coerce").fillna(0).sum())
            counts[taxid] = counts.get(taxid, 0.0) + n
            total += n
    if total <= 0 or len(counts) < 2:
        return None
    return {k: v / total for k, v in counts.items()}


def _parse_flags(argv) -> dict:
    extra = list(argv or [])
    out = {"C": 1.0, "max_iter": 500, "use_priors": False}
    i = 0
    while i < len(extra):
        tok = extra[i]
        if tok in {"--C", "-C"} and i + 1 < len(extra):
            out["C"] = float(extra[i + 1])
            i += 2
            continue
        if tok in {"--max-iter", "--max_iter"} and i + 1 < len(extra):
            out["max_iter"] = int(extra[i + 1])
            i += 2
            continue
        if tok in {"--use-priors", "--use_priors"}:
            out["use_priors"] = True
            i += 1
            continue
        i += 1
    return out


def reprofile(regenerated, ground_truth, initial, config) -> ReprofileResult:
    cfg = dict(config or {})
    flags = _parse_flags(cfg.get("extra_argv") or [])
    seed = int(cfg.get("seed") or 42)
    features_df = cfg.get("features_df")
    df = preprocess_data(regenerated)
    if df.empty or "true" not in df.columns:
        return ReprofileResult(tables={})
    feature_cols = sorted([c for c in df.columns if c not in {"seq", "true"}])
    X = df[feature_cols]
    y = df["true"]
    kwargs = {
        "C": flags["C"],
        "max_iter": flags["max_iter"],
        "random_state": seed,
        "solver": "lbfgs",
    }
    if flags["use_priors"]:
        priors = _priors_from_ground_truth(ground_truth or {})
        if priors:
            kwargs["class_weight"] = "balanced"
    model = LogisticRegression(**kwargs)
    model.fit(X, y)
    tables = profile_tables_with_model(
        initial,
        model,
        features_df=features_df,
        classify_unclassified=bool(cfg.get("classify_unclassified")),
    )
    dest = Path(cfg.get("output_dir") or ".")
    dest.mkdir(parents=True, exist_ok=True)
    save_model(model, str(dest / "trained_model.joblib"))
    return ReprofileResult(
        tables=tables,
        model=model,
        models={"LogisticRegression": model},
        feature_cols=feature_cols,
    )


def _cli(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regenerated", required=True)
    parser.add_argument("--ground-truth", dest="ground_truth", required=True)
    parser.add_argument("--initial", required=True)
    parser.add_argument("-o", dest="output_dir", required=True)
    parser.add_argument("--C", type=float, default=1.0)
    parser.add_argument("--max-iter", dest="max_iter", type=int, default=500)
    parser.add_argument("--use-priors", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    args, extra = parser.parse_known_args(argv)
    from samovar.reprofilers import load_csv_tables, load_regenerated_table

    config = {
        "output_dir": args.output_dir,
        "seed": args.seed,
        "extra_argv": [
            "--C",
            str(args.C),
            "--max-iter",
            str(args.max_iter),
            *(["--use-priors"] if args.use_priors else []),
            *extra,
        ],
    }
    result = reprofile(
        load_regenerated_table(args.regenerated),
        load_csv_tables(args.ground_truth, skip_prefixes=()),
        load_csv_tables(args.initial),
        config,
    )
    from samovar.reprofilers import write_reprofile_result

    write_reprofile_result(result, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
