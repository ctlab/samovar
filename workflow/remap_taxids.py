"""Build or apply a taxid→rank translation table (default genus) via ete3.

Does not rewrite annotation tables unless --output is given. The combined
annotation CSV must stay at original taxIDs for genome fetch.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from samovar.parse_annotators import (
    DEFAULT_TAX_RANK,
    default_rank_map_path,
    ensure_taxid_rank_map,
    remap_taxid_dataframe,
    taxid_value_columns,
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Cache taxid|{rank}_taxid translations; optionally remap a copy."
    )
    parser.add_argument(
        "--input",
        "-i",
        help="Annotation CSV used only to collect taxIDs (not rewritten).",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Optional remapped copy. Omit to update the cache only.",
    )
    parser.add_argument(
        "--cache",
        default=None,
        help="Translation table path (default: user cache taxid_genera_taxid.tsv).",
    )
    parser.add_argument(
        "--rank",
        default=DEFAULT_TAX_RANK,
        help="NCBI rank (genus, species, ...). 'none' skips lookup.",
    )
    args = parser.parse_args(argv)

    cache = args.cache or str(default_rank_map_path(args.rank))
    taxids = set()
    df = None
    if args.input:
        df = pd.read_csv(args.input)
        for col in taxid_value_columns(df.columns):
            taxids.update(df[col].dropna().unique())
    ensure_taxid_rank_map(taxids, args.rank, cache)
    if args.output:
        if df is None:
            raise SystemExit("--output requires --input")
        dest = Path(args.output)
        dest.parent.mkdir(parents=True, exist_ok=True)
        remap_taxid_dataframe(df, rank=args.rank, cache_path=cache).to_csv(dest, index=False)
    print(f"Rank map ({args.rank}): {cache}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
