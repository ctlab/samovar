#!/usr/bin/env python3
"""Dummy MAG QC: MAG dir → completeness/contamination TSV."""

from __future__ import annotations

import argparse
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", required=True)
    parser.add_argument("-d", dest="db", default="")
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    mag_dir = Path(args.c)
    rows = ["mag_id\tcompleteness\tcontamination\tscore"]
    if mag_dir.is_dir():
        for path in sorted(mag_dir.iterdir()):
            if path.suffix.lower() in {".fa", ".fna", ".fasta"}:
                rows.append(f"{path.stem}\t100.0\t0.0\t1.0")
    if len(rows) == 1:
        rows.append("mag1\t100.0\t0.0\t1.0")
    dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
