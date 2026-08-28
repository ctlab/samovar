#!/usr/bin/env python3
"""Linear (logistic) reprofiler — thin wrapper around the builtin.

Prefer ``samovar prepare --reprofiler linear``. Import this file only if you
need a custom name/flags overlay::

    samovar tools import -n linear \\
        --exec-path examples/reprofiling/linear_classifier.py --type ml \\
        --flags "--max-iter 500"
"""

from __future__ import annotations

import argparse
from samovar.reprofilers import (
    load_csv_tables,
    load_regenerated_table,
    reprofile_linear,
    write_reprofile_result,
)


def reprofile(regenerated, ground_truth, initial, config):
    return reprofile_linear(regenerated, ground_truth, initial, config)


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
    write_reprofile_result(result, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
