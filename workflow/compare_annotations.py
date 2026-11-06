#!/usr/bin/env python3
"""Python replacement for workflow/compare_annotations.R."""

from __future__ import annotations

import argparse
from pathlib import Path
import yaml

from samovar.viz_annotation import compare_annotations, normalize_viz_config


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Visualize SamovaR annotation tables")
    parser.add_argument("--annotation_dir", required=True)
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--csv", dest="csv_file", default=None)
    parser.add_argument("--show_top", type=int, default=0)
    parser.add_argument("--rank", default="genus")
    parser.add_argument("--split", action="store_true")
    parser.add_argument("--config", help="Stage YAML containing a viz mapping")
    parser.add_argument("--threads", "--cores", type=int, default=None,
                        help="Native Spearman threads; scalar metrics use vectorized counts")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--viz-metrics", action="append", default=None,
                        help="Comma-separated score column names; all or none")
    parser.add_argument(
        "--full",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Also write F1, R2, and CV plots. Default is the scores table and bar plot.",
    )
    parser.add_argument(
        "--opal",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Fill OPAL profile columns and write opal_scores.png. Default leaves them empty.",
    )
    parser.add_argument(
        "--viz-metrics-sample",
        type=int,
        default=None,
        help="Score a seeded subset of this many rows. 0 scores every row.",
    )
    parser.add_argument(
        "--spearman-full",
        action="store_true",
        help="Spearman on the concatenated per-read table (default: 10000 reads per sample, in C++)",
    )
    args = parser.parse_args(argv)
    config = {}
    if args.config:
        with open(args.config) as handle:
            config = yaml.safe_load(handle) or {}
    try:
        settings = normalize_viz_config(config.get("viz"))
    except ValueError as exc:
        parser.error(str(exc))
    for key, value in (("full", args.full), ("opal", args.opal),
                       ("metrics_sample", args.viz_metrics_sample), ("metrics", args.viz_metrics)):
        if value is not None:
            settings[key] = value
    try:
        viz = normalize_viz_config(settings)
    except ValueError as exc:
        parser.error(str(exc))

    if args.output_dir:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    from samovar.profile_pipeline import score
    from samovar.annotation_io import read_annotation_dir
    root = Path(args.annotation_dir).parent
    stage = Path(args.annotation_dir).name.removesuffix("_annotations")
    has_profiles = score(root, stage, {**config, "viz": viz})
    if has_profiles and read_annotation_dir(args.annotation_dir).empty:
        return 0

    compare_annotations(
        annotation_dir=args.annotation_dir,
        output_dir=args.output_dir,
        csv_file=args.csv_file,
        show_top=args.show_top,
        types=("f1", "R2", "cv", "scores") if viz["full"] else ("scores",),
        rank=args.rank,
        split=args.split,
        spearman_full=True if args.spearman_full else None,
        opal=viz["opal"],
        metrics_sample=viz["metrics_sample"],
        metrics=viz["metrics"],
        seed=args.seed if args.seed is not None else int(config.get("seed", 42)),
        threads=args.threads if args.threads is not None else int(config.get("cores", 1)),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
