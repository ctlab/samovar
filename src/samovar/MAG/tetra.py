"""Tetranucleotide profiles and a normal-error homology test.

Each 4-mer frequency is treated as normal with binomial sampling variance.
The variance is floored at an effective length so a long contig is not
judged more strictly than a typical assembly fragment. Two profiles match
when the root-mean-square z-score of their difference stays within that error.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np

K = 4
N_KMERS = 4**K
# Uncertainty of a profile estimated from about 20 kb, not from the whole contig.
EFFECTIVE_N = 20_000
DEFAULT_THRESHOLD = 1.6
_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def read_fasta_sequence(path: Path) -> str:
    opener = gzip.open if str(path).endswith(".gz") else open
    chunks = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith(">"):
                continue
            chunks.append(line.strip())
    raw = "".join(chunks).upper()
    return "".join(ch for ch in raw if ch in _BASE)


def kmer_counts(sequence: str) -> Tuple[np.ndarray, int]:
    counts = np.zeros(N_KMERS, dtype=np.float64)
    encoded = 0
    filled = 0
    n = 0
    for base in sequence:
        code = _BASE.get(base)
        if code is None:
            encoded = 0
            filled = 0
            continue
        encoded = ((encoded << 2) | code) & (N_KMERS - 1)
        filled += 1
        if filled >= K:
            counts[encoded] += 1.0
            n += 1
    return counts, n


def profile_rms_z(
    left: np.ndarray,
    n_left: int,
    right: np.ndarray,
    n_right: int,
    *,
    effective_n: int = EFFECTIVE_N,
) -> float:
    """RMS z-score of the frequency difference under a diagonal normal model."""
    n_left = int(n_left)
    n_right = int(n_right)
    if n_left < K or n_right < K:
        return 0.0
    freq_left = np.asarray(left, dtype=np.float64) / float(n_left)
    freq_right = np.asarray(right, dtype=np.float64) / float(n_right)
    scale_left = max(min(n_left, int(effective_n)), 1)
    scale_right = max(min(n_right, int(effective_n)), 1)
    var = freq_left * (1.0 - freq_left) / scale_left + freq_right * (1.0 - freq_right) / scale_right
    z = (freq_left - freq_right) / np.sqrt(var + 1e-18)
    return float(np.sqrt(np.mean(z * z)))


def profiles_match(rms: float, threshold: float = DEFAULT_THRESHOLD) -> bool:
    return float(rms) <= float(threshold)


def load_fasta_profiles(paths: Iterable[Path]) -> list:
    """Return dicts with path, counts, n, and length for each FASTA."""
    loaded = []
    for path in paths:
        sequence = read_fasta_sequence(Path(path))
        counts, n = kmer_counts(sequence)
        loaded.append(
            {
                "path": Path(path),
                "name": Path(path).stem,
                "counts": counts,
                "n": n,
                "length": len(sequence),
            }
        )
    return loaded
