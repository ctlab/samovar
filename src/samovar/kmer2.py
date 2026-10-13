#!/usr/bin/env python3
"""Dinucleotide (k=2) feature extractor.

Reads FASTQ via the annotator CLI (``-i``/``-I``/``-o``/``-t``) and writes a
headered TSV: ``seq`` plus ACGT 2-mer counts (``AA``, ``AC``, …, ``TT``).
Counting runs in ``bin/samovar_kmer`` (OpenMP over reads, ``-t`` threads).
``--kmer-max N`` keeps the N k-mers with the largest across-read variance
(at most 16; ties keep the earlier AA..TT column). Those are the columns
that land in the combined annotation table.
This is a Feature sub-contract of Annotation: no taxID column.
"""

from __future__ import annotations

import argparse
import gzip
import os
import subprocess
from pathlib import Path
from typing import Iterator, List, Optional, TextIO, Tuple


def _open_text(path: str) -> TextIO:
    gzipped = str(path).endswith(".gz")
    if not gzipped:
        try:
            with open(path, "rb") as raw:
                gzipped = raw.read(2) == b"\x1f\x8b"
        except OSError:
            gzipped = False
    if gzipped:
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_fastq_records(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(read_id, sequence)`` from a FASTQ file."""
    p = Path(path)
    if not path or not p.is_file() or p.stat().st_size == 0:
        return
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            handle.readline()
            handle.readline()
            stripped = header.strip()
            if not stripped:
                continue
            if stripped.startswith("@"):
                stripped = stripped[1:]
            tokens = stripped.split()
            if not tokens:
                continue
            read_id = tokens[0].replace("/1", "").replace("/2", "")
            yield read_id, seq.strip().upper()


def parse_output(file_path: str):
    """Contract helper: headered seq + k-mer columns."""
    from samovar.parse_annotators import read_custom_raw

    return read_custom_raw(file_path)


def kmer_binary() -> Path:
    from samovar.paths import repo_root

    return repo_root() / "bin" / "samovar_kmer"


def ensure_kmer_binary() -> Path:
    """Build ``bin/samovar_kmer`` when the source is newer than the binary."""
    from samovar.paths import cxx_compiler, repo_root

    binary = kmer_binary()
    root = repo_root()
    source = root / "src" / "cpp" / "kmer_profile.cpp"
    makefile = root / "src" / "cpp" / "Makefile"
    if not source.is_file():
        raise FileNotFoundError(f"C++ k-mer profiler source not found: {source}")
    newest = source.stat().st_mtime
    if makefile.is_file():
        newest = max(newest, makefile.stat().st_mtime)
    needs_build = (not binary.is_file()) or newest > binary.stat().st_mtime
    if needs_build:
        cxx = cxx_compiler()
        if not cxx:
            raise FileNotFoundError(
                "No C++ compiler found (g++, c++, or clang++). Install a compiler or set CXX."
            )
        env = os.environ.copy()
        env["CXX"] = cxx
        subprocess.check_call(["make", "-C", str(makefile.parent), "kmer"], env=env)
    if not os.access(binary, os.X_OK):
        raise RuntimeError(f"C++ k-mer profiler is not executable: {binary}")
    return binary


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Count ACGT dinucleotides (k=2) per read for the Feature contract."
    )
    parser.add_argument("-i", "-1", dest="r1", required=True, help="R1 FASTQ")
    parser.add_argument("-I", "-2", dest="r2", default="", help="R2 FASTQ (optional)")
    parser.add_argument("-d", dest="db", default="", help="Unused (annotator CLI)")
    parser.add_argument("-o", "--output", dest="o", required=True, help="Output TSV")
    parser.add_argument("-t", dest="threads", default="1", help="OpenMP threads for the C++ counter")
    parser.add_argument(
        "--kmer-max",
        dest="kmer_max",
        type=int,
        default=None,
        help="Max distinct k-mers written (highest across-read variance; at most 16)",
    )
    args = parser.parse_args(argv)
    if args.kmer_max is not None and args.kmer_max < 0:
        parser.error("--kmer-max must be >= 0")
    threads = int(args.threads)
    if threads < 1:
        parser.error("-t must be >= 1")
    binary = ensure_kmer_binary()
    cmd = [str(binary), "-i", args.r1, "-o", args.o, "-t", str(threads)]
    if args.r2:
        cmd.extend(["-I", args.r2])
    if args.db:
        cmd.extend(["-d", args.db])
    if args.kmer_max is not None:
        cmd.extend(["--kmer-max", str(args.kmer_max)])
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
