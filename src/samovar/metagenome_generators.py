"""Combined community+read simulators (``metagenome_generator``).

Unlike ``reads_generator`` (abundance → FASTQ), a metagenome generator owns
both community design and sequencing. Built-ins are CAMISIM modes. Custom
tools are imported with ``samovar tools import --type meta``.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from samovar.main_config import (
    flags_target_matches,
    imported_flags_for_names,
    iter_tools,
    merge_flag_strings,
    parse_tool_entry,
    tool_flags,
    tool_path,
)
from samovar.paths import load_config
from samovar.reads_generators import (
    CamisimReadsGenerator,
    CustomReadsGenerator,
    ReadsGenerator,
    extra_flags_argv,
)
from samovar.table_regenerators import extra_flags_argv as _split_flags

METAGENOME_GENERATOR_GROUP = "metagenome_generator"

_BUILTIN_ALIASES = {
    "camisim": "camisim",
    "cami": "camisim",
    "hybrid": "hybrid",
    "nanosim": "nanosim",
    "nanosim3": "nanosim",
    "ont": "nanosim",
    "nanopore": "nanosim",
    "simulator.py": "nanosim",
}


class MissingMetagenomeGeneratorError(ValueError):
    """``tools.<name>`` is missing or is not a ``metagenome_generator``."""


def resolve_metagenome_generator(name: Optional[str]) -> Tuple[str, str]:
    key = str(name or "").strip()
    if not key:
        raise MissingMetagenomeGeneratorError(
            "Empty metagenome_generator name. Import a tool with "
            "`samovar tools import -n NAME --type meta` or use camisim|hybrid|nanosim."
        )
    low = key.lower().replace("-", "_")
    if low in _BUILTIN_ALIASES:
        return "builtin", _BUILTIN_ALIASES[low]
    return "custom", key


def require_known_metagenome_generator(name: Optional[str]) -> str:
    kind, canon = resolve_metagenome_generator(name)
    if kind == "custom":
        lookup_metagenome_generator(canon)
    return canon


def flags_apply_to_metagenome_generator(target: str, name: Optional[str]) -> bool:
    try:
        kind, canon = resolve_metagenome_generator(name)
    except MissingMetagenomeGeneratorError:
        return False
    names = [canon, name]
    if kind == "builtin":
        if canon == "nanosim":
            names.extend(["nanosim3", "ont", "simulator.py"])
        if canon == "camisim":
            names.extend(["cami", "hybrid"])
        if canon == "hybrid":
            names.append("camisim")
    return flags_target_matches(
        target,
        *names,
        groups=("metagenome_generator", "metagenome", "meta"),
    )


def lookup_metagenome_generator(name: str) -> list:
    key = str(name or "").strip()
    if not key:
        raise MissingMetagenomeGeneratorError(
            "Empty metagenome_generator name. Import a tool with "
            "`samovar tools import -n NAME --type meta`."
        )
    tools = iter_tools(load_config())
    spec = tools.get(key)
    matched = key
    if spec is None:
        low = key.lower()
        for stored, value in tools.items():
            if stored.lower() == low:
                spec = value
                matched = stored
                break
    if spec is None:
        raise MissingMetagenomeGeneratorError(
            f"metagenome_generator {key!r} is not in the main install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type meta` "
            "before generate / prepare."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != METAGENOME_GENERATOR_GROUP:
        raise MissingMetagenomeGeneratorError(
            f"tools.{matched} has group {group!r}, expected "
            f"{METAGENOME_GENERATOR_GROUP!r}. Re-import with --type meta."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingMetagenomeGeneratorError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return parsed


def attach_metagenome_flags(name: Optional[str], config: Dict[str, Any]) -> Dict[str, Any]:
    cfg = dict(config or {})
    kind, canon = resolve_metagenome_generator(name or cfg.get("metagenome_generator"))
    tools = iter_tools(load_config())
    lookup_names = [canon, name, cfg.get("metagenome_generator")]
    if kind == "builtin" and canon == "nanosim":
        lookup_names.extend(["nanosim3", "ont", "simulator.py"])
    imported = imported_flags_for_names(tools, *[str(n) for n in lookup_names if n])
    if kind == "custom":
        spec = lookup_metagenome_generator(canon)
        imported = merge_flag_strings(tool_flags(spec, canon), imported)
    existing = merge_flag_strings(
        cfg.get("extra_flags"),
        cfg.get("metagenome_generator_flags"),
    )
    if imported and imported in existing:
        cfg["extra_flags"] = existing
    else:
        cfg["extra_flags"] = merge_flag_strings(imported, existing)
    cfg["extra_argv"] = extra_flags_argv(cfg.get("extra_flags"))
    cfg["metagenome_generator"] = canon
    return cfg


def get_metagenome_generator(name: Optional[str]) -> ReadsGenerator:
    kind, canon = resolve_metagenome_generator(name)
    if kind == "builtin":
        return CamisimReadsGenerator(canon)
    lookup_metagenome_generator(canon)
    return CustomReadsGenerator(canon)


def constant_abundance_frame(
    table: pd.DataFrame,
    n_reads: Optional[int] = None,
) -> pd.DataFrame:
    """Equal counts per taxid (combined community design for constant+ISS)."""
    if table is None or table.empty or "taxid" not in table.columns:
        return table
    out = table.copy()
    n_cols = [c for c in out.columns if c != "taxid" and "n" in str(c).lower()]
    if not n_cols:
        return out
    n_taxa = max(int(out["taxid"].nunique()), 1)
    target = int(n_reads) if n_reads not in (None, "", 0) else None
    if target is None:
        total = int(out[n_cols].sum().sum())
        target = total if total > 0 else n_taxa
    per = max(1, int(target) // n_taxa)
    for col in n_cols:
        out[col] = per
    return out


def parse_constant_flags(argv: Optional[List[str]]) -> Dict[str, Any]:
    """``--n-reads`` / ``--N_reads`` / ``--model`` from extra_argv."""
    args = list(argv or [])
    out: Dict[str, Any] = {}
    i = 0
    leftover: List[str] = []
    while i < len(args):
        tok = args[i]
        if tok in {"--n-reads", "--N_reads", "--total-reads"} and i + 1 < len(args):
            out["n_reads"] = int(args[i + 1])
            i += 2
            continue
        if tok.startswith("--n-reads=") or tok.startswith("--N_reads="):
            out["n_reads"] = int(tok.split("=", 1)[1])
            i += 1
            continue
        if tok == "--model" and i + 1 < len(args):
            out["model"] = args[i + 1]
            i += 2
            continue
        leftover.append(tok)
        i += 1
    out["iss_extra"] = leftover
    _ = _split_flags
    return out
