"""Combine annotator *.out files into per-sample CSVs.

The join is a C++ external sort-merge (bin/samovar_combine_annotations) so
Kraken k-mer columns are never loaded.

TaxIDs are kept at the annotators' original resolution. Collapsing to genus
(or another rank) happens only at visualization time via a cached
``taxid|genera_taxid`` translation table — rewriting this table would break
genome fetch (NCBI assemblies are strain/species IDs, not genus IDs).
"""

import argparse
import sys

from samovar.combine_tables import combine_with_cpp


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
    parser.add_argument(
        "--split_sample_name",
        "-s",
        type=int,
        default=1,
        help="How many underscore-separated tokens form the sample name (regenerated ISS uses 2).",
    )
    parser.add_argument(
        "--truth-table",
        "--ground-truth-table",
        dest="truth_table",
        default=None,
        help="CAMI-style or seq\\ttaxid table for the true column (initial annotation).",
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=500000,
        help="External-sort run size (rows) for the C++ combiner.",
    )
    parser.add_argument(
        "--rank",
        type=str,
        default=None,
        help="Ignored. Rank collapse is applied only during visualization.",
    )
    parser.add_argument(
        "--skip-rank",
        action="store_true",
        help="No-op (tables are never rank-collapsed after merge).",
    )
    parser.add_argument(
        "--skip-species",
        action="store_true",
        help="Deprecated no-op; species/genus remapping is not applied to combined tables.",
    )
    args = parser.parse_args(argv)

    combine_with_cpp(
        args.input_dir,
        args.output_dir,
        args.split_sample_name,
        args.chunk_rows,
        truth_table=args.truth_table,
    )
    from samovar.annotation_qc import EmptyAnnotatorsError, filter_classified_abundance_tables
    from samovar.abundance import input_to_abundance_tables, load_table_input

    try:
        loaded = load_table_input(args.output_dir)
        tables = input_to_abundance_tables(loaded)
        if tables:
            filter_classified_abundance_tables(
                tables,
                reports_dir=args.input_dir,
                fatal_if_none=False,
            )
    except EmptyAnnotatorsError:
        raise
    except Exception:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
