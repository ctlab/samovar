#!/usr/bin/env python3
"""Compatibility wrapper: constant-taxID annotator lives in ``samovar.baselines``."""

from samovar.baselines.constant_taxid import classify_fastq, iter_fastq_ids, main, parse_output

__all__ = ["classify_fastq", "iter_fastq_ids", "main", "parse_output"]


if __name__ == "__main__":
    raise SystemExit(main())
