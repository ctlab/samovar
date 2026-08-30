"""Canonical ``~/.config/samovar/config.json`` schema (no duplicated paths).

On disk the file is nested: ``compilers``, ``API``, ``genomes``, ``databases``,
``workflows``, ``tools``. Readers still accept the old flat keys
(``python_path``, ``iss_path``, ``tools.<name> = "/bin/..."``, …).
``write_config`` always stores the nested form.
"""

from __future__ import annotations

import os
import shutil
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

CANONICAL_TOP = (
    "version",
    "root",
    "compilers",
    "API",
    "genomes",
    "databases",
    "workflows",
    "tools",
    "custom-flags",
)

LEGACY_TOP = (
    "python_path",
    "r_path",
    "r_lib_path",
    "iss_path",
    "opal_path",
    "multiqc_path",
    "camisim_path",
    "nextflow_path",
    "nanosim_path",
    "art_path",
    "wgsim_path",
    "ncbi_email",
    "test_genomes",
    "genomes",
    "processed_genomes",
    "genome_dirs",
    "path",
    "extra_path",
    "tool_envs",
    "annotation_regenerate_r",
)

TOOL_GROUPS = (
    "runtime",
    "compiler",
    "annotator",
    "reads_generator",
    "metagenome_generator",
    "table_reads_generator",
    "table_scoring",
    "scoring",
    "reprofiler",
    "annotation_converter",
    "qc",
    "workflow",
)

TOOL_GROUP_ALIASES: Dict[str, str] = {
    "a": "annotator",
    "ann": "annotator",
    "annotator": "annotator",
    "annotators": "annotator",
    "runtime": "runtime",
    "compiler": "compiler",
    "reads": "reads_generator",
    "read": "reads_generator",
    "reads_generator": "reads_generator",
    "iss": "reads_generator",
    "meta": "metagenome_generator",
    "metagenome": "metagenome_generator",
    "metagenome_generator": "metagenome_generator",
    "table": "table_reads_generator",
    "table_reads": "table_reads_generator",
    "table_reads_generator": "table_reads_generator",
    "table_scoring": "table_scoring",
    "table_score": "table_scoring",
    "tablescoring": "table_scoring",
    "score": "scoring",
    "scoring": "scoring",
    "viz": "scoring",
    "visualisation": "scoring",
    "visualization": "scoring",
    "visualizations": "scoring",
    "visualisations": "scoring",
    "plots": "scoring",
    "ml": "reprofiler",
    "reprofile": "reprofiler",
    "reprofiler": "reprofiler",
    "reprofiling": "reprofiler",
    "converter": "annotation_converter",
    "convert": "annotation_converter",
    "annotation_converter": "annotation_converter",
    "annotation-converter": "annotation_converter",
    "qc": "qc",
    "quality": "qc",
    "trim": "qc",
    "trimming": "qc",
    "read_qc": "qc",
    "read-qc": "qc",
    "workflow": "workflow",
    "wf": "workflow",
}

# Scoring / viz tools: glob under the run output directory (6th tools.* slot).
DEFAULT_SCORING_INPUTS = "*annotations"

TOOL_GROUP_BY_NAME: Dict[str, str] = {
    "bash": "runtime",
    "python": "runtime",
    "python3": "runtime",
    "R": "runtime",
    "Rscript": "runtime",
    "g++": "compiler",
    "c++": "compiler",
    "clang++": "compiler",
    "kraken2": "annotator",
    "kaiju": "annotator",
    "kraken": "annotator",
    "krakenuniq": "annotator",
    "metaphlan": "annotator",
    "metaphlan4": "annotator",
    "centrifuge": "annotator",
    "iss": "reads_generator",
    "art": "reads_generator",
    "art_illumina": "reads_generator",
    "wgsim": "reads_generator",
    "seqtk": "reads_generator",
    "samtools": "reads_generator",
    "camisim": "metagenome_generator",
    "nanosim": "metagenome_generator",
    "nanosim3": "metagenome_generator",
    "simulator.py": "metagenome_generator",
    "sparsedossa2-fit": "table_reads_generator",
    "sparsedossa2-stool": "table_reads_generator",
    "sparsedossa2-vaginal": "table_reads_generator",
    "sparsedossa2-ibd": "table_reads_generator",
    "sparsedossa2-cv": "table_scoring",
    "opal": "scoring",
    "opal.py": "scoring",
    "multiqc": "scoring",
    "fastp": "qc",
    "cutadapt": "qc",
    "trimmomatic": "qc",
    "chopper": "qc",
    "nanofilt": "qc",
    "random_forest": "reprofiler",
    "adaboost": "reprofiler",
    "snakemake": "workflow",
    "nextflow": "workflow",
}

DEFAULT_WORKFLOWS: Dict[str, List[str]] = {
    "snakemake": [
        "annotators",
        "annotation2iss",
        "read_processing",
        "database_prep",
        "iss_test",
    ],
    "nextflow": ["camisim"],
    "conda": [],
}

# Legacy scalar key → (tool name in tools{}, compiler key, or special)
LEGACY_TOOL_KEYS = {
    "python_path": "python",
    "r_path": "R",
    "iss_path": "iss",
    "opal_path": "opal.py",
    "multiqc_path": "multiqc",
    "camisim_path": "camisim",
    "nextflow_path": "nextflow",
    "nanosim_path": "nanosim",
    "art_path": "art_illumina",
    "wgsim_path": "wgsim",
}

PathLike = Union[str, os.PathLike]


def _raw_get(cfg: Dict[str, Any], key: str, default: Any = None) -> Any:
    """Bypass InstallConfig.legacy_view to avoid recursion."""
    if isinstance(cfg, dict):
        return dict.get(cfg, key, default)
    return default


def _as_list(value: Any) -> List[str]:
    if value is None or value is False:
        return []
    if isinstance(value, (list, tuple)):
        out: List[str] = []
        for item in value:
            out.extend(_as_list(item))
        return out
    text = str(value).strip()
    if not text or text.startswith("_"):
        return []
    return [text]


def _split_dirs(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, dict):
        return [str(v).strip() for v in value.values() if str(v).strip()]
    if isinstance(value, (list, tuple)):
        items: List[str] = []
        for part in value:
            items.extend(_split_dirs(part))
        return items
    text = str(value).strip()
    if not text or text.startswith("_"):
        return []
    return [
        piece.strip()
        for piece in text.replace(";", ":").split(":")
        if piece.strip() and piece.strip() not in {"$PATH", "${PATH}"}
    ]


def tool_group_for(name: str) -> str:
    return TOOL_GROUP_BY_NAME.get(name, TOOL_GROUP_BY_NAME.get(Path(name).name, "runtime"))


def normalize_tool_group(value: str) -> str:
    """Map CLI ``--type`` / aliases onto ``TOOL_GROUPS``."""
    key = str(value or "").strip().lower().replace("-", "_")
    if key in TOOL_GROUP_ALIASES:
        return TOOL_GROUP_ALIASES[key]
    if key in TOOL_GROUPS:
        return key
    raise ValueError(
        f"Unknown tool type {value!r}. Use one of: {', '.join(TOOL_GROUPS)}"
    )


def parse_tool_entry(value: Any, name: str = "") -> List[str]:
    """Normalize a tools.* value to ``[env, workflow, path, group]``.

    On disk the preferred form is an object with ``exec`` / ``type`` /
    ``lazy-install`` / ``flags-translate``. List rows and path strings still parse.
    """
    from samovar.tool_spec import record_to_spec

    return record_to_spec(value, name)


def _trim_tool_spec(spec: List[str]) -> List[str]:
    """Keep 4, 5 (flags), or 6 (flags + inputs) elements; drop empty tails."""
    out = [str(x) if x is not None else "" for x in spec]
    while len(out) < 4:
        out.append("")
    if len(out) > 6:
        out = out[:6]
    flags = str(out[4]).strip() if len(out) > 4 else ""
    inputs = str(out[5]).strip() if len(out) > 5 else ""
    if inputs:
        return [out[0], out[1], out[2], out[3], flags, inputs]
    if flags:
        return [out[0], out[1], out[2], out[3], flags]
    return out[:4]


def tool_flags(entry: Any, name: str = "") -> str:
    spec = parse_tool_entry(entry, name)
    if len(spec) > 4:
        return str(spec[4]).strip()
    return ""


def tool_inputs(entry: Any, name: str = "") -> str:
    """6th tools.* slot (input glob). Scoring defaults to ``*annotations``."""
    spec = parse_tool_entry(entry, name)
    if len(spec) > 5 and str(spec[5]).strip():
        return str(spec[5]).strip()
    group = str(spec[3] or "").strip() if len(spec) > 3 else ""
    if group == "scoring":
        return DEFAULT_SCORING_INPUTS
    return ""


def merge_flag_strings(*parts: Optional[str]) -> str:
    chunks = [str(p).strip() for p in parts if p is not None and str(p).strip()]
    return " ".join(chunks)


def normalize_flag_target(value: str) -> str:
    return str(value or "").strip().lower().replace("-", "_")


def flags_target_matches(
    target: str,
    *names: str,
    groups: Optional[Sequence[str]] = None,
) -> bool:
    """True if ``--flags TARGET …`` should attach to a tool with these names."""
    token = normalize_flag_target(target)
    if not token:
        return False
    aliases = {normalize_flag_target(name) for name in names if name}
    if groups:
        aliases.update(normalize_flag_target(group) for group in groups)
    return token in aliases


def imported_flags_for_names(tools: Dict[str, List[str]], *names: str) -> str:
    """Return ``tools.<name>[4]`` for the first matching key (case-insensitive)."""
    from samovar.tool_spec import bare_tool_name

    for name in names:
        key = str(name or "").strip()
        if not key:
            continue
        if key in tools:
            return tool_flags(tools[key], key)
        want = bare_tool_name(key).lower()
        for stored, spec in tools.items():
            if stored.lower() == key.lower() or bare_tool_name(stored).lower() == want:
                return tool_flags(spec, stored)
    return ""


def tool_path(entry: Any, name: str = "") -> str:
    spec = parse_tool_entry(entry, name)
    path = spec[2]
    env = spec[0]
    if not path:
        return ""
    candidate = Path(path).expanduser()
    if env.lower() == "conda" or (candidate.is_dir() and not candidate.is_file()):
        binary = name if name not in {"nanosim", "nanosim3"} else "simulator.py"
        if name == "art":
            binary = "art_illumina"
        nested = candidate / "bin" / binary
        try:
            if nested.is_file():
                return str(nested)
        except OSError:
            pass
        nested2 = candidate / binary
        try:
            if nested2.is_file():
                return str(nested2)
        except OSError:
            pass
    return path


def tool_env_prefix(entry: Any, name: str = "") -> str:
    spec = parse_tool_entry(entry, name)
    env, _workflow, path, _group = spec[:4]
    if not path:
        return ""
    candidate = Path(path).expanduser()
    try:
        is_dir = candidate.is_dir()
        is_file = candidate.is_file()
    except OSError:
        is_dir = False
        is_file = False
    if env.lower() == "conda" and is_dir:
        return str(candidate)
    if is_file and candidate.parent.name == "bin":
        return str(candidate.parent.parent)
    return ""


def iter_tools(cfg: Dict[str, Any]) -> Dict[str, List[str]]:
    raw = _raw_get(cfg, "tools") if isinstance(_raw_get(cfg, "tools"), dict) else {}
    out: Dict[str, List[str]] = {}
    from samovar.tool_spec import bare_tool_name, parse_tool_record, record_to_spec

    by_bare: Dict[str, str] = {}
    for name, value in raw.items():
        if str(name).startswith("_"):
            continue
        spec = record_to_spec(value, str(name))
        out[str(name)] = spec
        bare = bare_tool_name(name)
        if bare and bare not in out:
            out[bare] = spec
            by_bare[bare] = str(name)
        elif bare:
            # Prefer versioned key as the canonical bare mapping.
            if ":" in str(name) and ":" not in by_bare.get(bare, ""):
                out[bare] = spec
                by_bare[bare] = str(name)
    return out


def iter_tool_records(cfg: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    raw = _raw_get(cfg, "tools") if isinstance(_raw_get(cfg, "tools"), dict) else {}
    from samovar.tool_spec import parse_tool_record

    out: Dict[str, Dict[str, Any]] = {}
    for name, value in raw.items():
        if str(name).startswith("_"):
            continue
        out[str(name)] = parse_tool_record(value, str(name))
    return out


def lookup_tool_record(cfg: Dict[str, Any], name: str) -> Optional[Dict[str, Any]]:
    from samovar.tool_spec import bare_tool_name

    records = iter_tool_records(cfg)
    if name in records:
        return records[name]
    want = bare_tool_name(name).lower()
    for key, rec in records.items():
        if bare_tool_name(key).lower() == want:
            return rec
    return None


def set_tool(
    cfg: Dict[str, Any],
    name: str,
    *,
    path: str = "",
    env: str = "",
    workflow: str = "",
    group: str = "",
    flags: Optional[str] = None,
    inputs: Optional[str] = None,
    lazy_install: Optional[str] = None,
    flags_translate: Any = None,
    version: Optional[str] = None,
) -> Dict[str, Any]:
    from samovar.tool_spec import (
        bare_tool_name,
        join_tool_key,
        parse_flags_translate,
        parse_tool_record,
        probe_tool_version,
        split_tool_key,
    )

    tools = dict(_raw_get(cfg, "tools") or {})
    bare, key_ver = split_tool_key(name)
    existing_key = name if name in tools else ""
    if not existing_key:
        for stored in list(tools):
            if bare_tool_name(stored) == bare:
                existing_key = stored
                break
    previous = parse_tool_record(tools.get(existing_key or name), bare or name)
    exe = dict(previous.get("exec") or {})
    if env:
        exe["env"] = env
    if workflow:
        exe["parser"] = workflow
    elif env and str(exe.get("parser") or "") in {"", "bash"}:
        exe["parser"] = env
    if path:
        exe["path"] = path
    if not exe.get("parser"):
        exe["parser"] = exe.get("env") or "bash"
    rec = {
        "exec": {
            "env": str(exe.get("env") or ""),
            "parser": str(exe.get("parser") or "bash"),
            "path": str(exe.get("path") or ""),
        },
        "type": group or previous.get("type") or tool_group_for(bare or name),
        "lazy-install": (
            str(lazy_install).strip()
            if lazy_install is not None
            else str(previous.get("lazy-install") or "")
        ),
        "flags": str(flags).strip() if flags is not None else str(previous.get("flags") or ""),
        "flags-translate": (
            parse_flags_translate(flags_translate)
            if flags_translate is not None
            else dict(previous.get("flags-translate") or {})
        ),
    }
    in_val = inputs if inputs is not None else previous.get("inputs")
    if in_val:
        rec["inputs"] = str(in_val).strip()
    if not rec["lazy-install"]:
        from samovar.tool_spec import lazy_install_for

        rec["lazy-install"] = lazy_install_for(bare or name, version or key_ver)
    if not rec["flags-translate"]:
        from samovar.tool_spec import DEFAULT_FLAGS_TRANSLATE

        rec["flags-translate"] = dict(DEFAULT_FLAGS_TRANSLATE.get(bare) or {})
    ver = str(version or key_ver or previous.get("_version") or "").strip()
    if not ver and rec["exec"]["path"] and bare not in {"samovar"}:
        ver = probe_tool_version(rec["exec"]["path"], bare)
    disk_key = join_tool_key(bare or name, ver)
    for stored in list(tools):
        if bare_tool_name(stored) == bare:
            tools.pop(stored, None)
    tools[disk_key] = rec
    dict.__setitem__(cfg, "tools", tools)
    return cfg


def compilers_of(cfg: Dict[str, Any]) -> Dict[str, Any]:
    block = _raw_get(cfg, "compilers")
    return dict(block) if isinstance(block, dict) else {}


def api_of(cfg: Dict[str, Any]) -> Dict[str, Any]:
    block = _raw_get(cfg, "API")
    if isinstance(block, dict):
        return dict(block)
    alt = _raw_get(cfg, "api")
    return dict(alt) if isinstance(alt, dict) else {}


def genomes_block(cfg: Dict[str, Any]) -> Dict[str, Any]:
    block = _raw_get(cfg, "genomes")
    if isinstance(block, dict):
        if any(k in block for k in ("test", "raw", "processed", "data")):
            return dict(block)
        return dict(block)
    return {}


def _looks_like_path_map_only(block: Dict[str, Any]) -> bool:
    return False


def compiler_python(cfg: Dict[str, Any]) -> str:
    compilers = compilers_of(cfg)
    return str(compilers.get("python") or _raw_get(cfg, "python_path") or "").strip()


def compiler_r(cfg: Dict[str, Any]) -> str:
    compilers = compilers_of(cfg)
    return str(compilers.get("R") or _raw_get(cfg, "r_path") or "").strip()


def compiler_r_libs(cfg: Dict[str, Any]) -> List[str]:
    compilers = compilers_of(cfg)
    libs = compilers.get("R_libs")
    if libs:
        return _split_dirs(libs)
    return _split_dirs(_raw_get(cfg, "r_lib_path"))


def compiler_python_libs(cfg: Dict[str, Any]) -> List[str]:
    compilers = compilers_of(cfg)
    libs = compilers.get("python_libs")
    if libs:
        return _split_dirs(libs)
    return _split_dirs(_raw_get(cfg, "path") or _raw_get(cfg, "extra_path"))


def compiler_bash(cfg: Dict[str, Any]) -> str:
    compilers = compilers_of(cfg)
    return str(compilers.get("bash") or "").strip()


def compiler_cpp(cfg: Dict[str, Any]) -> str:
    compilers = compilers_of(cfg)
    return str(compilers.get("cpp") or "").strip()


def ncbi_email_from_cfg(cfg: Dict[str, Any]) -> str:
    api = api_of(cfg)
    return str(api.get("ncbi_email") or _raw_get(cfg, "ncbi_email") or "").strip()


def test_genome_dirs_from_cfg(cfg: Dict[str, Any]) -> List[str]:
    block = genomes_block(cfg)
    if block:
        return _split_dirs(block.get("test"))
    return _split_dirs(_raw_get(cfg, "test_genomes"))


def is_home_path(path: str, home: Optional[Path] = None) -> bool:
    """True if ``path`` is under ``$HOME`` before or after resolving symlinks.

    HPC home quotas cannot hold NCBI caches. A symlink from ``~/.cache`` onto
    scratch still must not be stored as a ``$HOME`` path.
    """
    text = str(path or "").strip()
    if not text:
        return False
    home_lex = (home or Path.home()).expanduser()
    candidate = Path(text).expanduser()
    home_s = str(home_lex)
    cand_s = str(candidate)
    if cand_s == home_s or cand_s.startswith(home_s + os.sep):
        return True
    try:
        return candidate.resolve().is_relative_to(home_lex.resolve())
    except (OSError, ValueError):
        return False


def drop_home_paths(mapping: Dict[str, str], home: Optional[Path] = None) -> Dict[str, str]:
    cleaned = {
        fid: path for fid, path in mapping.items() if path and not is_home_path(path, home)
    }
    return unique_folder_map(cleaned)


def unique_folder_map(mapping: Dict[str, str]) -> Dict[str, str]:
    """Keep one id per resolved directory; prefer ``default``."""
    ordered: List[Tuple[str, str]] = []
    if "default" in mapping:
        ordered.append(("default", mapping["default"]))
    for fid, path in mapping.items():
        if fid != "default":
            ordered.append((fid, path))
    seen: set = set()
    out: Dict[str, str] = {}
    for fid, path in ordered:
        try:
            key = str(Path(path).expanduser().resolve())
        except OSError:
            key = path
        if key in seen:
            continue
        seen.add(key)
        out[fid] = path
    return out


def folder_map(value: Any) -> Dict[str, str]:
    """``raw`` / ``processed``: ``{folder_id: path}``."""
    if not value:
        return {}
    if isinstance(value, dict):
        # {"folder": "folder_id", "path": "..."} single record
        if "folder" in value and ("path" in value or "dir" in value):
            fid = str(value.get("folder") or "default").strip() or "default"
            path = str(value.get("path") or value.get("dir") or "").strip()
            return {fid: path} if path else {}
        out: Dict[str, str] = {}
        for key, val in value.items():
            if str(key).startswith("_"):
                continue
            if isinstance(val, dict):
                path = str(val.get("path") or val.get("dir") or val.get("folder") or "").strip()
            else:
                path = str(val or "").strip()
            if path:
                out[str(key)] = path
        return out
    if isinstance(value, str) and value.strip():
        return {"default": value.strip()}
    if isinstance(value, (list, tuple)):
        out = {}
        for i, item in enumerate(value):
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                out[str(item[0])] = str(item[1])
            elif str(item).strip():
                out[str(i) if i else "default"] = str(item).strip()
        return out
    return {}


def raw_genome_dirs(cfg: Dict[str, Any]) -> List[str]:
    block = genomes_block(cfg)
    if block:
        mapped = drop_home_paths(folder_map(block.get("raw")))
        dirs = list(mapped.values())
        if dirs:
            return dirs
    raw = _raw_get(cfg, "genomes")
    if isinstance(raw, str) and raw.strip() and not is_home_path(raw):
        return [raw.strip()]
    return []


def processed_genome_dirs(cfg: Dict[str, Any]) -> List[str]:
    block = genomes_block(cfg)
    if block:
        mapped = drop_home_paths(folder_map(block.get("processed")))
        dirs = list(mapped.values())
        if dirs:
            return dirs
    proc = _raw_get(cfg, "processed_genomes")
    if isinstance(proc, str) and proc.strip() and not is_home_path(proc):
        return [proc.strip()]
    return raw_genome_dirs(cfg)


def extra_genome_dirs(cfg: Dict[str, Any]) -> List[str]:
    """All library folders: raw map values plus legacy genome_dirs."""
    dirs: List[str] = []
    seen = set()

    def add(raw: str) -> None:
        text = str(raw or "").strip()
        if not text or text in seen or is_home_path(text):
            return
        seen.add(text)
        dirs.append(text)

    block = genomes_block(cfg)
    if block:
        for path in folder_map(block.get("raw")).values():
            add(path)
        for path in folder_map(block.get("processed")).values():
            add(path)
        data = block.get("data") if isinstance(block.get("data"), dict) else {}
        folders = folder_map(block.get("raw"))
        folders.update(folder_map(block.get("processed")))
        for _tax, rec in data.items():
            if not isinstance(rec, (list, tuple)) or len(rec) < 2:
                continue
            # New schema: [species, genome_id, database, file]; legacy: [acc, folder, file]
            fid = str(rec[2] if len(rec) >= 4 else rec[1])
            if fid in folders:
                add(folders[fid])
    for item in _split_dirs(_raw_get(cfg, "genome_dirs")):
        add(item)
    return dirs


def first_dir(values: Sequence[str]) -> str:
    for item in values:
        if str(item).strip():
            return str(item).strip()
    return ""


def databases_of(cfg: Dict[str, Any]) -> Dict[str, List[List[str]]]:
    """Legacy view: ``{tool: [[name, path, flags], …]}``.

    Preferred on-disk form is ``databases.<tool>.<name:version>`` objects
    (``path`` / ``flags`` / ``lazy-download``). List rows still parse.
    """
    from samovar.db_spec import databases_to_rows

    return databases_to_rows(cfg)


def set_database(
    cfg: Dict[str, Any],
    tool: str,
    name: str,
    *,
    path: str = "",
    flags: Optional[str] = None,
    lazy_download: Optional[str] = None,
    version: Optional[str] = None,
    url: Optional[str] = None,
) -> Dict[str, Any]:
    """Upsert ``databases.<tool>.<name:version>`` and return the record."""
    from samovar.db_spec import (
        databases_for_disk,
        iter_database_records,
        join_db_key,
        lazy_download_for,
        parse_database_record,
        split_db_key,
    )

    tool_name = str(tool or "").strip()
    bare, key_ver = split_db_key(name)
    if not tool_name:
        raise ValueError("database annotator --tool is required")
    if not bare:
        raise ValueError("database --name is required")
    grouped = iter_database_records(cfg)
    tool_block = dict(grouped.get(tool_name) or {})
    disk_key = join_db_key(bare, str(version or key_ver or "").strip())
    previous = tool_block.get(disk_key) or {}
    if not previous:
        for stored, rec in tool_block.items():
            if stored == name or (rec.get("name") == bare and not rec.get("_version") and not (version or key_ver)):
                previous = rec
                break
    rec = parse_database_record(previous, tool=tool_name, name=bare)
    rec["tool"] = tool_name
    rec["name"] = bare
    if path:
        rec["path"] = str(path).strip()
    if flags is not None:
        rec["flags"] = str(flags).strip()
    if lazy_download is not None:
        rec["lazy-download"] = str(lazy_download).strip()
    if url is not None:
        rec["url"] = str(url).strip()
    ver = str(version or key_ver or rec.get("_version") or "").strip()
    rec["_version"] = ver
    if not rec.get("lazy-download"):
        rec["lazy-download"] = lazy_download_for(tool_name, bare, ver, rec.get("url") or "")
    disk_key = join_db_key(bare, ver)
    tool_block[disk_key] = rec
    grouped[tool_name] = tool_block
    cfg["databases"] = databases_for_disk({"databases": {
        t: {k: v for k, v in block.items()}
        for t, block in grouped.items()
    }})
    return rec


def scan_test_genome_index(test_root: Path) -> Dict[str, List[str]]:
    """``{taxid}_test`` → 4-field records for bundled ISS stubs (never NCBI taxids)."""
    from samovar.genome_index import as_test_taxid, TEST_FOLDER_ID

    data: Dict[str, List[str]] = {}
    if not test_root.is_dir():
        return data
    for path in sorted(test_root.rglob("*")):
        if not path.is_file():
            continue
        name = path.name
        if name.startswith(".") or name == "test.fa":
            continue
        lower = name.lower()
        if lower.endswith(".faa") or lower.endswith(".faa.gz"):
            continue
        if not any(
            lower.endswith(ext)
            for ext in (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
        ):
            continue
        stem = name.split(".")[0]
        taxid = stem.split("-")[0]
        if taxid.endswith("_test"):
            taxid = taxid[: -len("_test")]
        if not taxid.isdigit():
            continue
        try:
            rel_parent = path.parent.relative_to(test_root).as_posix()
        except ValueError:
            rel_parent = path.parent.name
        rel = name if rel_parent in {".", ""} else f"{rel_parent}/{name}"
        key = as_test_taxid(taxid)
        if key in data and name.endswith(".gz"):
            continue
        data[key] = [key, key, TEST_FOLDER_ID, rel]
    return data


def empty_canonical(*, version: str = "", root: str = "") -> Dict[str, Any]:
    return {
        "version": version,
        "root": root,
        "compilers": {
            "bash": "",
            "python": "",
            "python_libs": [],
            "R": "",
            "R_libs": [],
            "cpp": "",
            "cpp_libs": [],
        },
        "API": {"ncbi_email": ""},
        "genomes": {
            "samovar_database": str(Path(root) / "genomes") if root else "",
            "taxdump": "",
            "gtdb": "",
            "test": [],
            "raw": {},
            "processed": {},
            "data": {},
        },
        "databases": {},
        "workflows": {k: list(v) for k, v in DEFAULT_WORKFLOWS.items()},
        "tools": {},
        "custom-flags": ["--threads", "--cores"],
    }


def _merge_folder_map(existing: Any, extra: Dict[str, str]) -> Dict[str, str]:
    mapped = folder_map(existing)
    mapped.update({k: v for k, v in extra.items() if v})
    return mapped


def apply_legacy_updates(cfg: Dict[str, Any], updates: Dict[str, Any]) -> Dict[str, Any]:
    """Write either canonical or legacy keys into ``cfg`` (mutates, returns cfg)."""
    for key, value in updates.items():
        if key in CANONICAL_TOP and key not in {"tools", "genomes", "custom-flags", "databases"}:
            cfg[key] = value
            continue
        if key == "compilers" and isinstance(value, dict):
            block = compilers_of(cfg)
            block.update(value)
            cfg["compilers"] = block
            continue
        if key == "API" and isinstance(value, dict):
            block = api_of(cfg)
            block.update(value)
            cfg["API"] = block
            continue
        if key == "samovar_database":
            block = genomes_block(cfg) or empty_canonical()["genomes"]
            block["samovar_database"] = str(value or "").strip()
            proc = folder_map(block.get("processed"))
            if block["samovar_database"]:
                proc["samovar_database"] = str(Path(block["samovar_database"]).expanduser() / "processed")
                proc.setdefault("default", proc["samovar_database"])
            block["processed"] = proc
            cfg["genomes"] = block
            continue
        if key == "genomes":
            if isinstance(value, dict) and any(
                k in value
                for k in ("test", "raw", "processed", "data", "samovar_database", "taxdump", "gtdb")
            ):
                block = genomes_block(cfg) or empty_canonical()["genomes"]
                for sub, val in value.items():
                    block[sub] = val
                cfg["genomes"] = block
            elif isinstance(value, str):
                block = genomes_block(cfg) or empty_canonical()["genomes"]
                raw = folder_map(block.get("raw"))
                if value.strip():
                    raw["default"] = value.strip()
                elif "default" in raw:
                    raw.pop("default", None)
                block["raw"] = raw
                cfg["genomes"] = block
            continue
        if key == "processed_genomes":
            block = genomes_block(cfg) or empty_canonical()["genomes"]
            proc = folder_map(block.get("processed"))
            text = str(value or "").strip()
            if text:
                proc["default"] = text
            else:
                proc.pop("default", None)
            block["processed"] = proc
            cfg["genomes"] = block
            continue
        if key == "genome_dirs":
            block = genomes_block(cfg) or empty_canonical()["genomes"]
            raw = folder_map(block.get("raw"))
            # Replace extra library folders; keep default/scratch if present
            reserved = {k: v for k, v in raw.items() if k in {"default", "scratch", "run"}}
            raw = dict(reserved)
            for i, item in enumerate(_split_dirs(value)):
                fid = Path(item).name or f"lib{i}"
                n = fid
                k = 2
                while n in raw and raw[n] != item:
                    n = f"{fid}_{k}"
                    k += 1
                raw[n] = item
            block["raw"] = raw
            cfg["genomes"] = block
            continue
        if key == "test_genomes":
            block = genomes_block(cfg) or empty_canonical()["genomes"]
            block["test"] = _split_dirs(value)
            cfg["genomes"] = block
            continue
        if key in {"taxdump", "taxdump_path", "ncbi_taxdump"}:
            block = genomes_block(cfg) or empty_canonical()["genomes"]
            block["taxdump"] = str(value or "").strip()
            cfg["genomes"] = block
            continue
        if key in {"gtdb", "gtdb_taxdump", "taxdump_gtdb"}:
            block = genomes_block(cfg) or empty_canonical()["genomes"]
            block["gtdb"] = str(value or "").strip()
            cfg["genomes"] = block
            continue
        if key == "ncbi_email":
            api = api_of(cfg)
            api["ncbi_email"] = str(value or "")
            cfg["API"] = api
            continue
        if key == "annotation_regenerate_r":
            set_tool(
                cfg,
                "annotation_regenerate.R",
                path=str(value or ""),
                env="",
                workflow="R",
                group="table_reads_generator",
            )
            continue
        if key in LEGACY_TOOL_KEYS:
            set_tool(cfg, LEGACY_TOOL_KEYS[key], path=str(value or ""))
            if key == "python_path":
                compilers = compilers_of(cfg)
                compilers["python"] = str(value or "")
                cfg["compilers"] = compilers
            if key == "r_path":
                compilers = compilers_of(cfg)
                compilers["R"] = str(value or "")
                cfg["compilers"] = compilers
            continue
        if key == "r_lib_path":
            compilers = compilers_of(cfg)
            compilers["R_libs"] = _split_dirs(value)
            cfg["compilers"] = compilers
            continue
        if key in {"path", "extra_path"}:
            compilers = compilers_of(cfg)
            compilers["python_libs"] = _split_dirs(value)
            cfg["compilers"] = compilers
            continue
        if key == "custom-flags" or key == "custom_flags":
            from samovar.tool_spec import custom_flags_from_cfg

            cfg["custom-flags"] = custom_flags_from_cfg({"custom-flags": value})
            continue
        if key == "databases":
            if isinstance(value, dict):
                from samovar.db_spec import iter_database_records

                incoming = iter_database_records({"databases": value})
                for tool, grouped in incoming.items():
                    for key_name, rec in grouped.items():
                        set_database(
                            cfg,
                            rec.get("tool") or tool,
                            rec.get("name") or key_name,
                            path=str(rec.get("path") or ""),
                            flags=rec.get("flags"),
                            lazy_download=rec.get("lazy-download"),
                            version=str(rec.get("_version") or rec.get("version") or ""),
                            url=rec.get("url"),
                        )
            continue
        if key == "tools":
            if isinstance(value, dict):
                for name, entry in value.items():
                    if str(name).startswith("_"):
                        continue
                    from samovar.tool_spec import parse_tool_record

                    rec = parse_tool_record(entry, str(name))
                    exe = rec.get("exec") or {}
                    set_tool(
                        cfg,
                        str(name),
                        env=str(exe.get("env") or ""),
                        workflow=str(exe.get("parser") or ""),
                        path=str(exe.get("path") or ""),
                        group=str(rec.get("type") or ""),
                        flags=str(rec.get("flags") or ""),
                        inputs=rec.get("inputs"),
                        lazy_install=rec.get("lazy-install"),
                        flags_translate=rec.get("flags-translate"),
                        version=str(rec.get("_version") or ""),
                    )
            continue
        if key == "tool_envs":
            if isinstance(value, dict):
                for name, prefix in value.items():
                    if str(name).startswith("_") or not prefix:
                        continue
                    existing_rec = lookup_tool_record(cfg, str(name))
                    existing = parse_tool_entry(existing_rec, str(name)) if existing_rec else ["", "", "", ""]
                    path = existing[2] or str(prefix)
                    set_tool(
                        cfg,
                        str(name),
                        env="conda",
                        workflow=str(name),
                        path=path if Path(str(existing[2] or "")).is_file() else str(prefix),
                        group=existing[3],
                    )
            continue
        if key == "databases" and isinstance(value, dict):
            cfg["databases"] = value
            continue
        if key == "workflows" and isinstance(value, dict):
            cfg["workflows"] = value
            continue
        cfg[key] = value
    return cfg


def migrate_legacy(raw: Dict[str, Any]) -> Dict[str, Any]:
    """Lift a flat 0.10.x config (or mixed) into the nested schema."""
    cfg = empty_canonical(
        version=str(raw.get("version") or ""),
        root=str(raw.get("root") or ""),
    )
    # Start from nested pieces already present
    if isinstance(raw.get("compilers"), dict):
        cfg["compilers"].update({k: v for k, v in raw["compilers"].items() if not str(k).startswith("_")})
    if isinstance(raw.get("API"), dict):
        cfg["API"].update(raw["API"])
    elif isinstance(raw.get("api"), dict):
        cfg["API"].update(raw["api"])
    if isinstance(raw.get("genomes"), dict) and any(
        k in raw["genomes"]
        for k in ("test", "raw", "processed", "data", "samovar_database", "taxdump", "gtdb")
    ):
        for key in ("test", "raw", "processed", "data", "samovar_database", "taxdump", "gtdb"):
            if key in raw["genomes"]:
                cfg["genomes"][key] = raw["genomes"][key]
    if isinstance(raw.get("databases"), dict):
        cfg["databases"] = dict(raw["databases"])
    if isinstance(raw.get("workflows"), dict):
        wf = dict(DEFAULT_WORKFLOWS)
        wf.update(raw["workflows"])
        cfg["workflows"] = wf
    if isinstance(raw.get("tools"), dict):
        from samovar.tool_spec import parse_tool_record

        for name, entry in raw["tools"].items():
            if str(name).startswith("_"):
                continue
            cfg["tools"][str(name)] = parse_tool_record(entry, str(name))
    if raw.get("custom-flags") or raw.get("custom_flags"):
        from samovar.tool_spec import custom_flags_from_cfg

        cfg["custom-flags"] = custom_flags_from_cfg(raw)

    # Overlay leftover flat keys
    leftovers = {k: v for k, v in raw.items() if k not in CANONICAL_TOP}
    apply_legacy_updates(cfg, leftovers)

    # Flat genomes string when genomes was not a nested dict
    if isinstance(raw.get("genomes"), str) and raw["genomes"].strip():
        apply_legacy_updates(cfg, {"genomes": raw["genomes"]})
    return cfg


def sync_by_keys(existing: Any, template: Any) -> Any:
    """Merge ``existing`` onto ``template`` by key (install-time config sync).

    Nested dicts recurse. Values already set in ``existing`` win, including
    extra keys that the template does not know about. Empty existing scalars
    and empty collections fall back to the template.
    """
    if isinstance(template, dict):
        src = existing if isinstance(existing, dict) else {}
        out: Dict[str, Any] = {}
        for key, tval in template.items():
            if str(key).startswith("_"):
                continue
            if key in src:
                out[key] = sync_by_keys(src[key], tval)
            else:
                out[key] = deepcopy(tval)
        for key, eval_ in src.items():
            if str(key).startswith("_"):
                continue
            if key not in out:
                out[key] = deepcopy(eval_)
        return out
    if existing is None:
        return deepcopy(template)
    if existing == "" and template not in (None, ""):
        return deepcopy(template)
    if existing in ([], {}) and template not in (None, "", [], {}):
        return deepcopy(template)
    return existing


def to_canonical(data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    raw = dict(data or {})
    nested = migrate_legacy(raw)
    nested = sync_by_keys(
        nested,
        empty_canonical(
            version=str(nested.get("version") or ""),
            root=str(nested.get("root") or ""),
        ),
    )
    genomes = nested.get("genomes") if isinstance(nested.get("genomes"), dict) else {}
    from samovar.genome_index import normalize_genome_data

    genomes["data"] = normalize_genome_data(genomes.get("data"))
    nested["genomes"] = genomes
    return nested


def legacy_view(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """In-memory aliases so ``cfg.get('python_path')`` still works."""
    tools = iter_tools(cfg)
    view: Dict[str, Any] = {}
    view["python_path"] = compiler_python(cfg)
    view["r_path"] = compiler_r(cfg)
    view["r_lib_path"] = ":".join(compiler_r_libs(cfg))
    view["iss_path"] = tool_path(tools.get("iss"), "iss")
    view["opal_path"] = tool_path(tools.get("opal.py") or tools.get("opal"), "opal.py")
    view["multiqc_path"] = tool_path(tools.get("multiqc"), "multiqc")
    view["camisim_path"] = tool_path(tools.get("camisim"), "camisim")
    view["nextflow_path"] = tool_path(tools.get("nextflow"), "nextflow")
    view["nanosim_path"] = tool_path(
        tools.get("nanosim") or tools.get("simulator.py") or tools.get("nanosim3"),
        "nanosim",
    )
    view["art_path"] = tool_path(tools.get("art_illumina") or tools.get("art"), "art_illumina")
    view["wgsim_path"] = tool_path(tools.get("wgsim"), "wgsim")
    view["ncbi_email"] = ncbi_email_from_cfg(cfg)
    tests = test_genome_dirs_from_cfg(cfg)
    view["test_genomes"] = tests[0] if tests else ""
    view["genomes"] = first_dir(raw_genome_dirs(cfg))
    view["processed_genomes"] = first_dir(processed_genome_dirs(cfg))
    view["samovar_database"] = str(
        (genomes_block(cfg) or {}).get("samovar_database") or ""
    )
    if not view["samovar_database"] and view["processed_genomes"]:
        view["samovar_database"] = view["processed_genomes"]
    view["taxdump"] = str((genomes_block(cfg) or {}).get("taxdump") or "")
    view["gtdb"] = str((genomes_block(cfg) or {}).get("gtdb") or "")
    view["genome_dirs"] = extra_genome_dirs(cfg)
    view["path"] = compiler_python_libs(cfg)
    regen = tools.get("annotation_regenerate.R")
    view["annotation_regenerate_r"] = tool_path(regen, "annotation_regenerate.R") if regen else ""
    view["tools"] = {name: spec[2] for name, spec in tools.items() if spec[2]}
    envs: Dict[str, str] = {}
    for name, spec in tools.items():
        prefix = tool_env_prefix(spec, name)
        if prefix:
            envs[name] = prefix
    view["tool_envs"] = envs
    return view


class InstallConfig(dict):
    """Dict that stores the nested schema and answers legacy ``.get()`` keys."""

    def get(self, key, default=None):  # type: ignore[override]
        if dict.__contains__(self, key):
            return dict.get(self, key, default)
        view = legacy_view(self)
        if key in view:
            return view[key]
        return default

    def __getitem__(self, key):  # type: ignore[override]
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            view = legacy_view(self)
            if key in view:
                return view[key]
            raise

    def __setitem__(self, key, value):  # type: ignore[override]
        if key == "tools":
            apply_legacy_updates(self, {"tools": value})
            return
        if key == "genomes" and not (
            isinstance(value, dict)
            and any(k in value for k in ("test", "raw", "processed", "data", "samovar_database"))
        ):
            apply_legacy_updates(self, {"genomes": value})
            return
        if key in CANONICAL_TOP:
            dict.__setitem__(self, key, value)
            return
        apply_legacy_updates(self, {key: value})

    def __contains__(self, key):  # type: ignore[override]
        if dict.__contains__(self, key):
            return True
        return key in legacy_view(self)


def as_install_config(data: Optional[Dict[str, Any]]) -> InstallConfig:
    canonical = to_canonical(data)
    return InstallConfig(canonical)


def infer_tool_env(path: str, conda_prefix: str = "") -> Tuple[str, str]:
    """Return ``(env, workflow_name)`` for a discovered binary or prefix."""
    if not path:
        return "", "bash"
    candidate = Path(path).expanduser()
    try:
        is_dir = candidate.is_dir()
        is_file = candidate.is_file()
    except OSError:
        is_dir = False
        is_file = False
    prefix = ""
    if is_file and candidate.parent.name == "bin":
        prefix = str(candidate.parent.parent)
    elif is_dir:
        prefix = str(candidate)
        if candidate.name == "bin":
            prefix = str(candidate.parent)
    if prefix and conda_prefix and os.path.normpath(prefix) == os.path.normpath(conda_prefix):
        return "", "bash"
    if prefix and (Path(prefix) / "conda-meta").is_dir():
        return "conda", Path(prefix).name
    if is_dir and not is_file:
        return "conda", candidate.name
    return "", "bash"


def build_install_config(
    *,
    root: str,
    python_path: str,
    version: str,
    existing: Optional[Dict[str, Any]] = None,
    discovered_tools: Optional[Dict[str, str]] = None,
    ncbi_email: str = "",
    genomes_default: str = "",
    processed_default: str = "",
    extra_genome_dirs: Optional[Sequence[str]] = None,
    extra_path: Optional[Sequence[str]] = None,
    samovar_database: str = "",
    taxdump: str = "",
    bash: str = "",
    cxx: str = "",
    r_path: str = "",
    r_libs: Optional[Sequence[str]] = None,
    conda_prefix: str = "",
    conda_sidecars: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Assemble the nested install config from discovery + previous file."""
    cfg = as_install_config(existing)
    cfg["version"] = version
    cfg["root"] = root
    compilers = compilers_of(cfg)
    compilers["python"] = python_path or compilers.get("python") or ""
    compilers["bash"] = bash or compilers.get("bash") or shutil.which("bash") or ""
    compilers["cpp"] = cxx or compilers.get("cpp") or shutil.which("g++") or shutil.which("c++") or shutil.which("clang++") or ""
    compilers["R"] = r_path or compilers.get("R") or shutil.which("R") or ""
    if r_libs is not None:
        compilers["R_libs"] = list(r_libs)
    elif not compilers.get("R_libs"):
        compilers["R_libs"] = []
    if extra_path is not None:
        compilers["python_libs"] = [str(x) for x in extra_path if str(x).strip()]
    elif not compilers.get("python_libs"):
        compilers["python_libs"] = []
    if not compilers.get("cpp_libs"):
        compilers["cpp_libs"] = []
    cfg["compilers"] = compilers

    api = api_of(cfg)
    if ncbi_email:
        api["ncbi_email"] = ncbi_email
    elif not api.get("ncbi_email"):
        api["ncbi_email"] = ""
    cfg["API"] = api

    genomes = genomes_block(cfg) or empty_canonical()["genomes"]
    test_root = str(Path(root) / "data" / "test_genomes")
    genomes["test"] = [test_root]
    db = str(samovar_database or genomes.get("samovar_database") or "").strip()
    if not db:
        db = str(Path(root) / "genomes")
    if db and not is_home_path(db):
        genomes["samovar_database"] = db
        proc_store = str(Path(db).expanduser() / "processed")
    else:
        genomes["samovar_database"] = str(Path(root) / "genomes")
        proc_store = str(Path(root) / "genomes" / "processed")
    dump = str(taxdump or genomes.get("taxdump") or "").strip()
    if not dump:
        dump = str(Path(genomes["samovar_database"]).expanduser() / "taxdump")
    if dump and not is_home_path(dump):
        genomes["taxdump"] = dump
    else:
        genomes["taxdump"] = str(Path(genomes["samovar_database"]).expanduser() / "taxdump")
    raw = drop_home_paths(folder_map(genomes.get("raw")))
    if genomes_default and not is_home_path(genomes_default):
        raw["default"] = genomes_default
    for i, item in enumerate(extra_genome_dirs or []):
        text = str(item).strip()
        if not text or is_home_path(text):
            continue
        fid = Path(text).name or f"lib{i}"
        n = fid
        k = 2
        while n in raw and raw[n] != text:
            n = f"{fid}_{k}"
            k += 1
        raw[n] = text
    genomes["raw"] = raw
    proc = drop_home_paths(folder_map(genomes.get("processed")))
    proc["samovar_database"] = proc_store
    if processed_default and not is_home_path(processed_default):
        proc["default"] = processed_default
    else:
        proc.setdefault("default", proc_store)
    genomes["processed"] = proc
    from samovar.genome_index import normalize_genome_data

    migrated = normalize_genome_data(genomes.get("data"))
    for taxid, rec in scan_test_genome_index(Path(test_root)).items():
        migrated.setdefault(taxid, rec)
    genomes["data"] = migrated
    cfg["genomes"] = genomes

    if not isinstance(cfg.get("databases"), dict):
        cfg["databases"] = {}

    wf = dict(DEFAULT_WORKFLOWS)
    if isinstance(cfg.get("workflows"), dict):
        wf.update(cfg["workflows"])
    sidecars = list(conda_sidecars or [])
    if sidecars:
        conda_wf = list(wf.get("conda") or [])
        for name in sidecars:
            if name not in conda_wf:
                conda_wf.append(name)
        wf["conda"] = conda_wf
    cfg["workflows"] = wf

    discovered = dict(discovered_tools or {})
    # Keep previously configured tools, overlay discovery when missing path
    for name, path in discovered.items():
        existing_spec = parse_tool_entry(
            lookup_tool_record(cfg, name) or (cfg.get("tools") or {}).get(name),
            name,
        )
        if existing_spec[2]:
            continue
        env, workflow = infer_tool_env(path, conda_prefix=conda_prefix)
        set_tool(cfg, name, path=path, env=env, workflow=workflow, group=tool_group_for(name))

    # Compilers also as tools
    if compilers.get("python"):
        set_tool(cfg, "python", path=str(compilers["python"]), group="runtime")
        py3 = str(Path(compilers["python"]).parent / "python3")
        if Path(py3).is_file():
            set_tool(cfg, "python3", path=py3, group="runtime")
    if compilers.get("R"):
        set_tool(cfg, "R", path=str(compilers["R"]), workflow="R", group="runtime")
        rscript = shutil.which("Rscript") or str(Path(compilers["R"]).parent / "Rscript")
        if Path(rscript).is_file() or shutil.which("Rscript"):
            set_tool(
                cfg,
                "Rscript",
                path=str(Path(rscript).resolve()) if Path(rscript).is_file() else (shutil.which("Rscript") or rscript),
                workflow="R",
                group="runtime",
            )
    if compilers.get("bash"):
        set_tool(cfg, "bash", path=str(compilers["bash"]), workflow="bash", group="runtime")
    if compilers.get("cpp"):
        set_tool(cfg, "g++", path=str(compilers["cpp"]), workflow="bash", group="compiler")

    from samovar.version import get_version as _pkg_ver

    set_tool(
        cfg,
        "samovar",
        path=str(compilers.get("python") or ""),
        workflow="bash",
        group="runtime",
        lazy_install="pip install -e .",
        version=_pkg_ver(),
    )

    if not cfg.get("custom-flags"):
        from samovar.tool_spec import DEFAULT_CUSTOM_FLAGS

        cfg["custom-flags"] = list(DEFAULT_CUSTOM_FLAGS)

    return dict(cfg)


def disk_payload(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Nested keys only — what is written to config.json."""
    canonical = to_canonical(cfg)
    payload = {key: canonical.get(key) for key in CANONICAL_TOP}
    genomes = payload.get("genomes") if isinstance(payload.get("genomes"), dict) else {}
    if genomes:
        genomes["raw"] = drop_home_paths(folder_map(genomes.get("raw")))
        genomes["processed"] = drop_home_paths(folder_map(genomes.get("processed")))
        payload["genomes"] = genomes
    payload["tools"] = _tools_for_disk(canonical)
    payload["databases"] = _databases_for_disk(canonical)
    from samovar.tool_spec import custom_flags_from_cfg

    payload["custom-flags"] = custom_flags_from_cfg(canonical)
    return payload


def _tools_for_disk(cfg: Dict[str, Any]) -> Dict[str, Any]:
    from samovar.tool_spec import (
        DEFAULT_FLAGS_TRANSLATE,
        bare_tool_name,
        join_tool_key,
        lazy_install_for,
        parse_tool_record,
        probe_tool_version,
        split_tool_key,
    )
    from samovar.version import get_version

    raw = _raw_get(cfg, "tools") if isinstance(_raw_get(cfg, "tools"), dict) else {}
    out: Dict[str, Any] = {}
    for name, value in raw.items():
        if str(name).startswith("_"):
            continue
        rec = parse_tool_record(value, str(name))
        exe = dict(rec.get("exec") or {})
        bare, key_ver = split_tool_key(name)
        ver = str(rec.pop("_version", None) or key_ver or "").strip()
        if bare == "samovar" and not ver:
            ver = get_version()
        elif not ver and exe.get("path") and bare not in {"samovar"}:
            ver = probe_tool_version(str(exe.get("path")), bare)
        if not rec.get("lazy-install"):
            rec["lazy-install"] = lazy_install_for(bare, ver)
        if not rec.get("flags-translate"):
            rec["flags-translate"] = dict(DEFAULT_FLAGS_TRANSLATE.get(bare) or {})
        rec["exec"] = {
            "env": str(exe.get("env") or ""),
            "parser": str(exe.get("parser") or "bash"),
            "path": str(exe.get("path") or ""),
        }
        rec["type"] = rec.get("type") or tool_group_for(bare)
        rec["flags"] = str(rec.get("flags") or "")
        key = join_tool_key(bare, ver)
        out[key] = rec
    return out


def _databases_for_disk(cfg: Dict[str, Any]) -> Dict[str, Any]:
    from samovar.db_spec import databases_for_disk

    return databases_for_disk(cfg)


def _fmt_install_value(value: Any) -> str:
    if value is None or value is False:
        return ""
    if isinstance(value, dict):
        parts: List[str] = []
        for key in ("path", "index", "dir", "flags", "db", "name"):
            item = value.get(key)
            if item not in (None, "", [], {}):
                parts.append(f"{key}={item}")
        if parts:
            return ", ".join(parts)
        leftover = [
            f"{key}={item}"
            for key, item in value.items()
            if not str(key).startswith("_") and item not in (None, "", [], {})
        ]
        return ", ".join(leftover[:6])
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item).strip() for item in value if str(item).strip())
    return str(value).strip()


def install_option_rows(cfg: Optional[Dict[str, Any]]) -> List[Tuple[str, str]]:
    """User-facing install options (paths, databases, tools) — not the full JSON."""
    if not cfg:
        return []
    payload = disk_payload(cfg)
    rows: List[Tuple[str, str]] = []

    def add(label: str, value: Any) -> None:
        text = _fmt_install_value(value)
        if text:
            rows.append((label, text))

    add("version", payload.get("version"))
    add("root", payload.get("root"))
    compilers = payload.get("compilers") if isinstance(payload.get("compilers"), dict) else {}
    for key in ("python", "bash", "cpp", "R"):
        add(f"compilers.{key}", compilers.get(key))
    add("compilers.python_libs", compilers.get("python_libs"))
    add("compilers.R_libs", compilers.get("R_libs"))
    api = payload.get("API") if isinstance(payload.get("API"), dict) else {}
    add("API.ncbi_email", api.get("ncbi_email"))
    genomes = payload.get("genomes") if isinstance(payload.get("genomes"), dict) else {}
    add("genomes.samovar_database", genomes.get("samovar_database"))
    add("genomes.taxdump", genomes.get("taxdump"))
    add("genomes.gtdb", genomes.get("gtdb"))
    for name, path in folder_map(genomes.get("raw")).items():
        add(f"genomes.raw.{name}", path)
    for name, path in folder_map(genomes.get("processed")).items():
        add(f"genomes.processed.{name}", path)
    data = genomes.get("data")
    if isinstance(data, dict) and data:
        add("genomes.data", f"{len(data)} catalog entries")
    databases = payload.get("databases") if isinstance(payload.get("databases"), dict) else {}
    from samovar.db_spec import looks_like_db_map

    for name, spec in databases.items():
        if isinstance(spec, dict) and looks_like_db_map(spec):
            for dbname, rec in spec.items():
                add(f"databases.{name}.{dbname}", rec)
        else:
            add(f"databases.{name}", spec if spec not in (None, "", {}, []) else "(named, empty)")
    skip_tools = {"python", "python3", "g++", "c++", "clang++", "R", "Rscript", "bash"}
    for name, spec in iter_tools(payload).items():
        if name in skip_tools:
            continue
        path = tool_path(spec, name)
        if path:
            add(f"tools.{name}", path)
    workflows = payload.get("workflows") if isinstance(payload.get("workflows"), dict) else {}
    add("workflows.conda", workflows.get("conda"))
    return rows


def format_install_report(
    *,
    payload: Dict[str, Any],
    previous: Optional[Dict[str, Any]] = None,
    previous_path: str = "",
) -> str:
    """Human-readable previous-vs-new install options for ``install.sh``."""
    had_previous = isinstance(previous, dict) and bool(previous)
    new_rows = install_option_rows(payload)
    old_map = dict(install_option_rows(previous if had_previous else None))
    aligned: List[Tuple[str, str]] = []
    installed: List[Tuple[str, str]] = []
    for label, value in new_rows:
        old = old_map.get(label)
        if had_previous and old == value:
            aligned.append((label, value))
        else:
            installed.append((label, value))

    lines: List[str] = []
    if had_previous:
        where = previous_path.strip() or "config.json"
        lines.append(f"Previous installation config: found ({where})")
    else:
        where = previous_path.strip()
        if where:
            lines.append(f"Previous installation config: not found (writing {where})")
        else:
            lines.append("Previous installation config: not found")

    def _block(title: str, rows: List[Tuple[str, str]], empty: str) -> None:
        lines.append("")
        lines.append(title)
        if not rows:
            lines.append(f"  {empty}")
            return
        width = max(len(label) for label, _ in rows)
        for label, value in rows:
            lines.append(f"  {label:<{width}}  {value}")

    if had_previous:
        _block(
            "Options aligned from the previous install:",
            aligned,
            "(none — previous values were empty or all refreshed this run)",
        )
    _block(
        "New / updated options this install:",
        installed,
        "(none — this install only reused previous values)",
    )
    return "\n".join(lines)
