"""Intentional toy-database gaps used only in tests and example pipelines.

This is **not** a production recommendation. Custom Kraken2 / Kaiju indexes
built from ``data/test_genomes`` can drop one taxon each so missing-taxon
behaviour is visible in the toy MultiQC report.
"""

from __future__ import annotations

import logging
import os
import warnings
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

# Escherichia coli / genus Escherichia (NCBI).
ESCHERICHIA_TAXIDS = frozenset({"561", "562", "511145"})
# Escherichia phage phiX174 (and the historic phiX taxid).
PHIX_TAXIDS = frozenset({"2886930", "10847"})

EXAMPLE_OMIT = {
    "kraken2": {
        "label": "Escherichia",
        "taxids": ESCHERICHIA_TAXIDS,
    },
    "kaiju": {
        "label": "Phage Phi X",
        "taxids": PHIX_TAXIDS,
    },
}

EXAMPLE_OMIT_WARNING = (
    "WARNING: this is just an EXAMPLE (test / toy databases only), not a "
    "recommended production setting. The Kraken2 toy index omits Escherichia "
    "(taxids 562, 561, 511145). The Kaiju toy index omits Phage Phi X "
    "(taxids 2886930, 10847). Real metagenome indexes should keep every taxon."
)


def example_omit_requested(explicit: Optional[bool] = None) -> bool:
    if explicit is True:
        return True
    if explicit is False:
        return False
    flag = os.environ.get("SAMOVAR_EXAMPLE_OMIT", "").strip().lower()
    return flag in {"1", "true", "yes", "on"}


def input_dirs_look_like_test_genomes(directories: Sequence[str]) -> bool:
    for raw in directories:
        text = str(raw).replace("\\", "/").lower()
        if "test_genomes" in text:
            return True
        try:
            if (Path(raw).resolve() / "meta").is_dir() and (Path(raw).name == "test_genomes"):
                return True
        except OSError:
            continue
    return False


def should_apply_example_omit(
    directories: Sequence[str],
    explicit: Optional[bool] = None,
) -> bool:
    if example_omit_requested(explicit):
        return True
    return input_dirs_look_like_test_genomes(directories)


def warn_example_omit(db_type: str, omitted: Iterable[str]) -> None:
    skipped = ", ".join(sorted({str(t) for t in omitted})) or "(none matched)"
    spec = EXAMPLE_OMIT.get(db_type, {})
    label = spec.get("label", db_type)
    message = (
        f"{EXAMPLE_OMIT_WARNING} "
        f"Building {db_type}: omitting {label} (matched taxids: {skipped})."
    )
    warnings.warn(message, UserWarning, stacklevel=2)
    logger.warning(message)


def filter_example_omit(
    file_taxid_map: Mapping[str, str],
    db_type: str,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return (kept, omitted) path→taxid maps for a toy index."""
    spec = EXAMPLE_OMIT.get(str(db_type).lower())
    if not spec:
        return dict(file_taxid_map), {}
    drop = spec["taxids"]
    kept: Dict[str, str] = {}
    omitted: Dict[str, str] = {}
    for path, taxid in file_taxid_map.items():
        token = str(taxid).split(".")[0]
        if token in drop:
            omitted[path] = str(taxid)
        else:
            kept[path] = str(taxid)
    return kept, omitted
