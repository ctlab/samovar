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
    args = parser.parse_args(argv)

    if args.output_dir:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    compare_annotations(
        annotation_dir=args.annotation_dir,
        output_dir=args.output_dir,
        csv_file=args.csv_file,
        show_top=args.show_top,
        types=("f1", "R2", "cv"),
        rank=args.rank,
        split=args.split,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
