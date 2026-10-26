#!/usr/bin/env python3
"""Python replacement for workflow/compare_annotations.R."""

from __future__ import annotations

import argparse
from pathlib import Path

from samovar.viz_annotation import compare_annotations


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Visualize SamovaR annotation tables")
    parser.add_argument("--annotation_dir", required=True)
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--csv", dest="csv_file", default=None)
    parser.add_argument("--show_top", type=int, default=0)
    parser.add_argument("--rank", default="genus")
    parser.add_argument("--split", action="store_true")
    parser.add_argument(
        "--full",
        action="store_true",
        help="Also write F1, R2, and CV plots. Default is the scores table and bar plot.",
    )
    parser.add_argument(
        "--opal",
        action="store_true",
        help="Fill OPAL profile columns and write opal_scores.png. Default leaves them empty.",
    )
    parser.add_argument(
        "--viz-metrics-sample",
        type=int,
        default=0,
        help="Score a seeded subset of this many rows. 0 scores every row.",
    )
    parser.add_argument(
        "--spearman-full",
        action="store_true",
        help="Spearman on the concatenated per-read table (default: 10000 reads per sample, in C++)",
    )
    args = parser.parse_args(argv)

    if args.output_dir:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    compare_annotations(
        annotation_dir=args.annotation_dir,
        output_dir=args.output_dir,
        csv_file=args.csv_file,
        show_top=args.show_top,
        types=("f1", "R2", "cv", "scores") if args.full else ("scores",),
        rank=args.rank,
        split=args.split,
        spearman_full=True if args.spearman_full else None,
        opal=args.opal,
        metrics_sample=args.viz_metrics_sample,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
