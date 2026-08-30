"""GC-content FASTQ filter (``samovar tools import --type QC``).

Drops a read (and its mate) when GC fraction of ACGT bases is outside
``config['min_gc']``–``config['max_gc']`` (defaults 0 and 1). Extra CLI:
``--min-gc`` / ``--max-gc``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from samovar.qc import gc_fraction, iter_fastq_records, write_fastq_records


def _keep(seq: str, min_gc: float, max_gc: float) -> bool:
    value = gc_fraction(seq)
    return min_gc <= value <= max_gc


def trim(
    r1: str,
    r2: Optional[str],
    dest_r1: str,
    dest_r2: Optional[str],
    config: Dict[str, Any],
) -> List[str]:
    cfg = dict(config or {})
    min_gc = float(cfg.get("min_gc", 0.0))
    max_gc = float(cfg.get("max_gc", 1.0))
    kept_r1: List[Tuple[str, str, str, str]] = []
    kept_r2: List[Tuple[str, str, str, str]] = []
    r1_recs = list(iter_fastq_records(Path(r1)))
    r2_path = Path(r2) if r2 else None
    r2_recs = list(iter_fastq_records(r2_path)) if r2_path is not None and r2_path.is_file() else []
    if r2_recs and len(r2_recs) != len(r1_recs):
        n = min(len(r1_recs), len(r2_recs))
        r1_recs = r1_recs[:n]
        r2_recs = r2_recs[:n]
    if r2_recs:
        for a, b in zip(r1_recs, r2_recs):
            if _keep(a[1], min_gc, max_gc) and _keep(b[1], min_gc, max_gc):
                kept_r1.append(a)
                kept_r2.append(b)
    else:
        for rec in r1_recs:
            if _keep(rec[1], min_gc, max_gc):
                kept_r1.append(rec)
    write_fastq_records(Path(dest_r1), kept_r1)
    out = [str(dest_r1)]
    if dest_r2:
        write_fastq_records(Path(dest_r2), kept_r2 if r2_recs else [])
        out.append(str(dest_r2))
    return out
