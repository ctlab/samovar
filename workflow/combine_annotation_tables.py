"""Combine annotator *.out files into per-sample CSVs.

The join is a C++ external sort-merge (bin/samovar_combine_annotations) so
Kraken k-mer columns are never loaded. Species-rank remapping is a streaming
pass over the resulting CSV (unique taxIDs only go through taxonomy lookup).
"""

import argparse
import sys
from pathlib import Path

from samovar.combine_tables import apply_species_level_csv, combine_with_cpp


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Sort-merge annotator outputs into per-sample annotation CSVs."
    )
    parser.add_argument("--input_dir", "-i", type=str, required=True)
    parser.add_argument("--output_dir", "-o", type=str, required=True)
    parser.add_argument(
        "--true_annotation",
        "-t",
        type=str,
        required=False,
        default="(?<=taxid:)[0-9]+",
        help="Kept for CLI compatibility; true taxIDs are taken from taxid:<digits> in seq IDs.",
    )
    parser.add_argument("--split_sample_name", "-s", type=int, required=False, default=1)
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=500000,
        help="External-sort run size (rows) for the C++ combiner.",
    )
    parser.add_argument(
        "--skip-species",
        action="store_true",
        help="Do not remap taxIDs to species after the merge.",
    )
    args = parser.parse_args(argv)

    combine_with_cpp(args.input_dir, args.output_dir, args.split_sample_name, args.chunk_rows)

    if not args.skip_species:
        for csv_path in sorted(Path(args.output_dir).glob("*.annotation.csv")):
            apply_species_level_csv(str(csv_path), level="species")
            print(f"Species-mapped {csv_path.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
