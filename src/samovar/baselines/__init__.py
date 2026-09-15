"""Built-in baseline tools for every ``samovar tools import`` contract.

These are first-class implementations (identity / constant / passthrough),
not test stubs. Contract pytest and ``--pytest`` on import run against them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

_DIR = Path(__file__).resolve().parent

# group → baseline script (package-relative)
BASELINE_TOOLS: Dict[str, Path] = {
    "annotator": _DIR / "constant_taxid.py",
    "table_reads_generator": _DIR / "identity_table.py",
    "sample_scoring": _DIR / "constant_sample_score.py",
    "sample_filtering": _DIR / "passthrough_sample_filter.py",
    "scoring": _DIR / "count_files.py",
    "reads_generator": _DIR / "identity_reads.py",
    "metagenome_generator": _DIR / "identity_reads.py",
    "annotation_converter": _DIR / "identity_converter.py",
    "export": _DIR / "identity_export.py",
    "qc": _DIR / "gc_filter.py",
    "assembler": _DIR / "identity_assembler.py",
    "gene_caller": _DIR / "translate_orfs.py",
    "binner": _DIR / "identity_binner.py",
    "binner_qc": _DIR / "passthrough_binner_qc.py",
    "binner_combine": _DIR / "identity_binner_combine.py",
    "mag_taxonomy": _DIR / "constant_mag_taxonomy.py",
    "aligner": _DIR / "identity_aligner.py",
    "mag_quantifier": _DIR / "constant_mag_quantifier.py",
    "taxon_quantifier": _DIR / "identity_taxon_quantifier.py",
    "read_assigner": _DIR / "identity_read_assigner.py",
}

# Names that resolve to the constant-taxID annotator (legacy aliases kept).
CONSTANT_TAXID_NAMES = frozenset(
    {
        "constant_taxid",
        "constant9606",
        "constant",
        "dummy",
        "dummy9606",
        "random",
    }
)


def baseline_path(group: str) -> Path:
    path = BASELINE_TOOLS.get(group)
    if path is None:
        raise KeyError(f"no built-in baseline for contract group {group!r}")
    return path
