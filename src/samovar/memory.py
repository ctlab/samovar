"""Memory and thread budgets shared by annotators, reprofile, and apply.

``--cores`` is the machine CPU budget. A stage may use only the jobs whose
working set fits in ``SAMOVAR_MEM_MB`` (or most of MemAvailable). Prediction
does not inherit the training thread count: its probability matrix is sized
per chunk, and its thread ceiling is separate.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Optional, Union

import numpy as np
import pandas as pd

# Prediction parallelizes over trees on one chunk. A higher count multiplies
# temporary buffers without shrinking the reads × classes matrix.
PREDICT_THREAD_CAP = 2
# Leave this fraction of the detected RAM for the OS and file cache.
_MEM_FRACTION = 0.85
_MEMMAP_DEFAULT = 8 * 1024 * 1024


def _meminfo_kb(key: str) -> Optional[int]:
    try:
        with open("/proc/meminfo", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(key):
                    return int(line.split()[1])
    except OSError:
        return None
    return None


def available_mem_mb() -> int:
    """RAM this process may plan for, before the headroom fraction."""
    kb = _meminfo_kb("MemAvailable:")
    if kb is None:
        kb = _meminfo_kb("MemTotal:")
    if not kb:
        return 8192
    return max(1, kb // 1024)


_MEM_UNITS = {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 * 1024}
_MEM_TEXT = re.compile(r"^(\d+(?:\.\d+)?)([KMGT])?B?$")


def parse_mem(value: Union[str, int, float, None]) -> int:
    """SLURM ``--mem`` value as MiB.

    A bare number is MiB. ``K``, ``M``, ``G``, and ``T`` scale it; a trailing
    ``B`` (``16GB``) is accepted.
    """
    if isinstance(value, bool) or value is None:
        raise argparse.ArgumentTypeError("memory value is required (SLURM --mem, e.g. 16G)")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        mb = int(math.ceil(float(value)))
        if mb < 1:
            raise argparse.ArgumentTypeError("memory must be at least 1 MiB")
        return mb
    text = str(value).strip().upper().replace(" ", "")
    match = _MEM_TEXT.fullmatch(text)
    if not match:
        raise argparse.ArgumentTypeError(
            f"invalid memory {value!r}; use SLURM --mem form: 16000, 16G, 512M"
        )
    number = float(match.group(1))
    scale = _MEM_UNITS[match.group(2) or "M"]
    mb = int(math.ceil(number * scale))
    if mb < 1:
        raise argparse.ArgumentTypeError("memory must be at least 1 MiB")
    return mb


def add_mem_argument(parser: argparse.ArgumentParser, *, default=None) -> None:
    """``--mem`` next to ``--cores``, with ``--mem-mb`` as an alias."""
    parser.add_argument(
        "--mem",
        "--mem-mb",
        "--mem_mb",
        dest="mem_mb",
        type=parse_mem,
        default=default,
        metavar="MEM",
        help=(
            "Maximum memory for this run, as in SLURM --mem. "
            "A bare number is MiB; K, M, G, and T are accepted (16G). "
            "Default: 85%% of available RAM."
        ),
    )


def memory_limit_mb(requested: Optional[int] = None) -> int:
    """Ceiling in MiB. Explicit value, else ``SAMOVAR_MEM_MB``, else 85% of RAM."""
    if requested is not None:
        try:
            return max(1, int(requested))
        except (TypeError, ValueError):
            pass
    raw = os.environ.get("SAMOVAR_MEM_MB", "").strip()
    if raw:
        try:
            return max(1, int(raw))
        except ValueError:
            pass
    return max(1, int(available_mem_mb() * _MEM_FRACTION))


def rss_mb() -> int:
    """Resident set of this process, or 0 when /proc is unavailable."""
    try:
        with open("/proc/self/status", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    return max(0, int(line.split()[1]) // 1024)
    except OSError:
        return 0
    return 0


def remaining_mb(limit_mb: Optional[int] = None) -> int:
    """MiB still under the ceiling. At least 1 so a chunk is never empty."""
    limit = memory_limit_mb(limit_mb)
    return max(1, limit - rss_mb())


def job_mem_mb(threads: int, cores: int, limit_mb: Optional[int] = None) -> int:
    """Snakemake ``mem_mb`` share. Thread shares of one stage sum to the ceiling."""
    limit = memory_limit_mb(limit_mb)
    width = max(1, int(cores or 1))
    slots = max(1, min(int(threads or 1), width))
    return max(1, min(limit, (limit * slots) // width))


def parallel_jobs(
    requested: int,
    *,
    nbytes: int = 0,
    mem_mb: Optional[int] = None,
    ceiling: Optional[int] = None,
) -> int:
    """Jobs that fit in memory.

    Each worker is assumed to be able to hold one copy of ``nbytes``. When the
    table is small this returns ``requested`` (use the CPUs). When it is large
    the count drops so the copies stay under the ceiling.
    """
    try:
        jobs = int(requested if requested is not None else 1)
    except (TypeError, ValueError):
        jobs = 1
    if jobs < 1:
        jobs = 1
    if ceiling is not None:
        jobs = min(jobs, max(1, int(ceiling)))
    if nbytes > 0:
        copy_mb = max(1.0, float(nbytes) / (1024.0 * 1024.0))
        affordable = max(1, int((memory_limit_mb(mem_mb) * 0.7) / copy_mb))
        jobs = min(jobs, affordable)
    return max(1, jobs)


def predict_n_jobs(requested: int = 1) -> int:
    """Thread count for ``predict`` / ``predict_proba``. Not the training width."""
    return parallel_jobs(requested, ceiling=PREDICT_THREAD_CAP)


def proba_chunk_rows(
    n_classes: int,
    n_features: int,
    mem_mb: Optional[int] = None,
) -> int:
    """Rows whose feature block and probability matrix fit in a quarter of the budget."""
    per = max(1, (max(1, int(n_classes)) + max(1, int(n_features))) * 8)
    budget = max(8, memory_limit_mb(mem_mb) // 4) * 1024 * 1024
    rows = max(1, budget // per)
    override = os.environ.get("SAMOVAR_PREDICT_CHUNK", "").strip()
    if override:
        try:
            return max(1, min(rows, int(override)))
        except ValueError:
            pass
    return int(rows)


def memmap_min_bytes() -> int:
    raw = os.environ.get("SAMOVAR_MEMMAP_MIN_BYTES", "").strip()
    if raw:
        try:
            return max(0, int(raw))
        except ValueError:
            pass
    return _MEMMAP_DEFAULT


@contextmanager
def predict_threads(model, requested: int = 1) -> Iterator[int]:
    """Cap ``n_jobs`` for the duration of a predict call, then restore it."""
    jobs = predict_n_jobs(requested)
    saved = []
    targets = [model, getattr(model, "estimator", None)]
    for obj in targets:
        if obj is None or not hasattr(obj, "n_jobs"):
            continue
        saved.append((obj, obj.n_jobs))
        obj.n_jobs = jobs
    try:
        yield jobs
    finally:
        for obj, old in saved:
            obj.n_jobs = old


def read_csv_memmap(path, dest_dir=None) -> pd.DataFrame:
    """Read a CSV in chunks. Numeric columns stay in ``numpy.memmap`` files.

    Small files (under ``memmap_min_bytes``) use ``pandas.read_csv``.
    """
    file = Path(path)
    if file.stat().st_size < memmap_min_bytes():
        return pd.read_csv(file)
    store = Path(dest_dir) if dest_dir else file.parent / f".{file.stem}.memmap"
    store.mkdir(parents=True, exist_ok=True)
    numeric: list[str] = []
    handles = []
    object_parts: dict[str, list] = {}
    rows = 0
    try:
        for chunk in pd.read_csv(file, chunksize=50_000):
            if not numeric and not object_parts:
                for col in chunk.columns:
                    if pd.api.types.is_numeric_dtype(chunk[col]):
                        numeric.append(str(col))
                    else:
                        object_parts[str(col)] = []
                handles = [
                    open(store / f"{index}.f64", "wb") for index in range(len(numeric))
                ]
            for name, handle in zip(numeric, handles):
                values = pd.to_numeric(chunk[name], errors="coerce").fillna(0)
                handle.write(np.ascontiguousarray(values, dtype=np.float64).tobytes())
            for name, parts in object_parts.items():
                if name in chunk.columns:
                    parts.append(chunk[name].astype(str).to_numpy())
            rows += len(chunk)
    finally:
        for handle in handles:
            handle.close()
    if rows == 0:
        return pd.read_csv(file)
    data = {}
    for index, name in enumerate(numeric):
        data[name] = np.memmap(
            store / f"{index}.f64", dtype=np.float64, mode="r", shape=(rows,)
        )
    for name, parts in object_parts.items():
        data[name] = np.concatenate(parts) if parts else np.array([], dtype=object)
    ordered = list(pd.read_csv(file, nrows=0).columns)
    return pd.DataFrame({col: data[col] for col in ordered if col in data}, copy=False)
