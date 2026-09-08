#!/usr/bin/env python3
"""Dummy taxon quantifier: join MAG counts with MAG taxonomy."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", required=True)
    parser.add_argument("-x", required=True)
    parser.add_argument("-o", required=True)
    args = parser.parse_args(argv)
    tax = {}
    tax_path = Path(args.x)
    if tax_path.is_file():
        for line in tax_path.read_text(encoding="utf-8").splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 2:
                tax[parts[0]] = parts[1]
    counts = defaultdict(int)
    abund = Path(args.a)
    if abund.is_file():
        lines = abund.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t") if lines else []
        n_idx = 1
        for i, name in enumerate(header):
            if name == "N" or name.startswith("N_"):
                n_idx = i
                break
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) <= n_idx:
                continue
            mag = parts[0]
            n = int(float(parts[n_idx]))
            counts[tax.get(mag, "0")] += n
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not counts:
        counts["562"] = 10
    out = ["taxid\tN_1"]
    for taxid, n in sorted(counts.items()):
        out.append(f"{taxid}\t{n}")
    dest.write_text("\n".join(out) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
