#!/usr/bin/env python3
"""Dummy annotator for contract tests (``--type annotator``).

CLI: -i R1 -I R2 -d DB -o OUT -t threads
Output: two-column TSV seq, taxID. parse_output() reads the same shape.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_output(file_path: str) -> pd.DataFrame:
    path = Path(file_path)
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["seq", "taxID"])
    df = pd.read_table(path, header=None)
    df = df.iloc[:, :2]
    df.columns = ["seq", "taxID"]
    return df


def _iter_fastq_ids(path: str):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return
    with p.open() as handle:
        for i, line in enumerate(handle):
            if i % 4 == 0:
                yield line[1:].split()[0].strip()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-I", dest="r2", default="")
    parser.add_argument("-d", dest="db", default="")
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for seq in _iter_fastq_ids(args.i):
        rows.append(f"{seq}\t562")
    dest.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
