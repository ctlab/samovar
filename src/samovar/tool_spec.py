"""Canonical tools.* records: exec, type, lazy-install, flags-translate.

On disk a tool is an object keyed by ``name:version``::

    "kraken2:2.1.3": {
      "exec": {"env": "", "parser": "bash", "path": "/usr/bin/kraken2"},
      "type": "annotator",
      "lazy-install": "conda install -y bioconda::kraken2=2.1.3",
      "flags": "",
      "flags-translate": {"--threads": "--threads"}
    }

Legacy list rows ``[env, parser, path, group, flags?, inputs?]`` still parse.
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

DEFAULT_CUSTOM_FLAGS = ["--threads", "--cores"]

# Versions this tree was last exercised against (install.sh warns on mismatch).
TESTED_TOOL_VERSIONS: Dict[str, str] = {
    "samovar": "",  # filled from package version at write time
    "python": "3.12",
    "kraken2": "2.1.3",
    "kaiju": "1.10.1",
    "iss": "2.0.0",
    "snakemake": "8.0.0",
    "multiqc": "1.21",
}

DEFAULT_LAZY_INSTALL: Dict[str, str] = {
    "kraken2": "conda install -y bioconda::kraken2",
    "kaiju": "conda install -y bioconda::kaiju",
    "kraken": "conda install -y bioconda::kraken",
    "krakenuniq": "conda install -y bioconda::krakenuniq",
    "metaphlan": "conda install -y bioconda::metaphlan",
    "metaphlan4": "conda install -y bioconda::metaphlan",
    "centrifuge": "conda install -y bioconda::centrifuge",
    "iss": "pip install 'insilicoseq>2.0.0'",
    "snakemake": "pip install snakemake",
    "multiqc": "pip install 'multiqc>=1.21'",
    "nextflow": "conda install -y bioconda::nextflow",
    "wgsim": "conda install -y bioconda::wgsim",
    "art_illumina": "conda install -y bioconda::art",
    "seqtk": "conda install -y bioconda::seqtk",
    "samtools": "conda install -y bioconda::samtools",
    "opal.py": "pip install opal",
    "camisim": "pip install CAMISIM",
    "nanosim": "conda install -y bioconda::nanosim",
    "nanosim3": "conda install -y bioconda::nanosim",
    "fastp": "conda install -y bioconda::fastp",
    "cutadapt": "conda install -y bioconda::cutadapt",
    "trimmomatic": "conda install -y bioconda::trimmomatic",
    "chopper": "conda install -y bioconda::chopper",
    "nanofilt": "conda install -y bioconda::nanofilt",
}

DEFAULT_FLAGS_TRANSLATE: Dict[str, Dict[str, str]] = {
    "kraken2": {"--threads": "--threads", "--cores": "--threads"},
    "kaiju": {"--threads": "-z", "--cores": "-z"},
    "kraken": {"--threads": "--threads", "--cores": "--threads"},
    "krakenuniq": {"--threads": "--threads", "--cores": "--threads"},
    "metaphlan": {"--threads": "--nproc", "--cores": "--nproc"},
    "metaphlan4": {"--threads": "--nproc", "--cores": "--nproc"},
    "centrifuge": {"--threads": "--threads", "--cores": "--threads"},
    "iss": {"--threads": "--cpus", "--cores": "--cpus"},
    "snakemake": {"--threads": "--cores", "--cores": "--cores"},
    "art_illumina": {"--threads": "--p", "--cores": "--p"},
    "wgsim": {"--threads": "--threads", "--cores": "--threads"},
    "custom": {"--threads": "-t", "--cores": "-t"},
    "fastp": {"--threads": "--thread", "--cores": "--thread"},
    "cutadapt": {"--threads": "--cores", "--cores": "--cores"},
    "trimmomatic": {"--threads": "-threads", "--cores": "-threads"},
    "chopper": {"--threads": "--threads", "--cores": "--threads"},
}

_VERSION_TOKEN = re.compile(r"(\d+\.\d+(?:\.\d+)?)")


def split_tool_key(key: str) -> Tuple[str, str]:
    text = str(key or "").strip()
    if not text:
        return "", ""
    if ":" in text:
        name, _, ver = text.partition(":")
        return name.strip(), ver.strip()
    return text, ""


def join_tool_key(name: str, version: str = "") -> str:
    bare = str(name or "").strip()
    ver = str(version or "").strip()
    if ver:
        return f"{bare}:{ver}"
    return bare


def bare_tool_name(key: str) -> str:
    return split_tool_key(key)[0]


def parse_flags_translate(raw: Any) -> Dict[str, str]:
    """``--flags-translate \"--threads:--threads --cores:-t\"`` or a mapping."""
    out: Dict[str, str] = {}
    if raw is None or raw is False:
        return out
    if isinstance(raw, dict):
        for src, dest in raw.items():
            s, d = str(src).strip(), str(dest).strip()
            if s and d:
                out[_norm_flag(s)] = d if d.startswith("-") else d
        return out
    if isinstance(raw, (list, tuple)):
        for item in raw:
            if isinstance(item, dict):
                out.update(parse_flags_translate(item))
            elif isinstance(item, (list, tuple)) and len(item) >= 2:
                out[_norm_flag(item[0])] = str(item[1]).strip()
            else:
                out.update(parse_flags_translate(str(item)))
        return out
    text = str(raw).strip()
    if not text:
        return out
    for tok in text.split():
        if ":" not in tok:
            continue
        src, _, dest = tok.partition(":")
        src, dest = src.strip(), dest.strip()
        if src and dest:
            out[_norm_flag(src)] = dest
    return out


def _norm_flag(flag: str) -> str:
    text = str(flag or "").strip()
    if not text:
        return ""
    if not text.startswith("-"):
        text = f"--{text}"
    return text


def flags_translate_for(name: str, record: Optional[Mapping[str, Any]] = None) -> Dict[str, str]:
    bare = bare_tool_name(name)
    mapping = dict(DEFAULT_FLAGS_TRANSLATE.get(bare) or {})
    if record:
        mapping.update(parse_flags_translate(record.get("flags-translate") or record.get("flags_translate")))
    return mapping


def lazy_install_for(name: str, version: str = "") -> str:
    bare = bare_tool_name(name)
    cmd = DEFAULT_LAZY_INSTALL.get(bare, "")
    if cmd and version and "conda install" in cmd and "=" not in cmd.split()[-1]:
        pkg = cmd.rsplit(None, 1)[-1]
        if "::" in pkg and "=" not in pkg:
            return cmd.rsplit(None, 1)[0] + f" {pkg}={version}"
    return cmd


def parse_tool_record(value: Any, name: str = "") -> Dict[str, Any]:
    """Return the on-disk object form (always a dict)."""
    from samovar.main_config import tool_group_for

    bare, key_ver = split_tool_key(name)
    group = tool_group_for(bare or name)
    env, parser, path, group_out, flags, inputs = "", "bash", "", group, "", ""
    lazy = ""
    translate: Dict[str, str] = {}
    version = key_ver

    if isinstance(value, dict) and isinstance(value.get("exec"), dict):
        exe = value.get("exec") or {}
        env = str(exe.get("env") or "").strip()
        parser = str(exe.get("parser") or exe.get("workflow") or "").strip()
        path = str(exe.get("path") or "").strip()
        group_out = str(value.get("type") or value.get("group") or group).strip() or group
        flags = str(value.get("flags") or "").strip()
        inputs = str(value.get("inputs") or value.get("glob") or "").strip()
        lazy = str(value.get("lazy-install") or value.get("lazy_install") or "").strip()
        translate = parse_flags_translate(value.get("flags-translate") or value.get("flags_translate"))
        version = str(value.get("version") or version).strip()
    elif isinstance(value, dict):
        env = str(value.get("env") or "").strip()
        parser = str(value.get("workflow") or value.get("parser") or "").strip()
        path = str(value.get("path") or "").strip()
        group_out = str(value.get("group") or value.get("type") or group).strip() or group
        flags = str(value.get("flags") or value.get("extra") or "").strip()
        inputs = str(
            value.get("inputs") or value.get("input") or value.get("glob") or ""
        ).strip()
        lazy = str(value.get("lazy-install") or value.get("lazy_install") or "").strip()
        translate = parse_flags_translate(value.get("flags-translate") or value.get("flags_translate"))
        version = str(value.get("version") or version).strip()
        if isinstance(value.get("exec"), str):
            path = path or str(value.get("exec")).strip()
    elif isinstance(value, (list, tuple)):
        parts = [str(x).strip() if x is not None else "" for x in value]
        while len(parts) < 4:
            parts.append("")
        env, parser, path, group_out = parts[0], parts[1], parts[2], parts[3] or group
        if len(parts) > 4:
            flags = parts[4]
        if len(parts) > 5:
            inputs = parts[5]
        if not path and len(parts) == 1:
            path = parts[0]
    elif value not in (None, False, ""):
        path = str(value).strip()

    if not parser:
        parser = env if env else "bash"
    if not group_out:
        group_out = group
    if not translate:
        translate = dict(DEFAULT_FLAGS_TRANSLATE.get(bare) or {})
    if not lazy:
        lazy = lazy_install_for(bare, version)

    rec: Dict[str, Any] = {
        "exec": {"env": env, "parser": parser, "path": path},
        "type": group_out,
        "lazy-install": lazy,
        "flags": flags,
        "flags-translate": translate,
    }
    if inputs:
        rec["inputs"] = inputs
    if version:
        rec["_version"] = version
    return rec


def _trim_spec(spec: List[str]) -> List[str]:
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


def record_to_spec(record: Mapping[str, Any], name: str = "") -> List[str]:
    """List form used by existing ``parse_tool_entry`` callers."""
    rec = parse_tool_record(value=record, name=name)
    exe = rec.get("exec") or {}
    return _trim_spec(
        [
            str(exe.get("env") or ""),
            str(exe.get("parser") or "bash"),
            str(exe.get("path") or ""),
            str(rec.get("type") or ""),
            str(rec.get("flags") or ""),
            str(rec.get("inputs") or ""),
        ]
    )


def probe_tool_version(path: str, name: str = "") -> str:
    exe = str(path or "").strip()
    if not exe:
        return ""
    candidate = Path(exe).expanduser()
    cmd = str(candidate) if candidate.is_file() else exe
    for args in ([cmd, "--version"], [cmd, "-v"], [cmd, "version"]):
        try:
            proc = subprocess.run(
                args,
                capture_output=True,
                timeout=8,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        blob = (proc.stdout or b"") + b"\n" + (proc.stderr or b"")
        if isinstance(blob, bytes):
            blob = blob.decode("utf-8", errors="replace")
        match = _VERSION_TOKEN.search(blob)
        if match:
            return match.group(1)
    bare = bare_tool_name(name)
    if bare in {"python", "python3"}:
        try:
            proc = subprocess.run(
                [cmd, "-c", "import sys; print(sys.version.split()[0])"],
                capture_output=True,
                text=True,
                timeout=8,
            )
            ver = (proc.stdout or "").strip()
            if ver:
                return ver
        except (OSError, subprocess.TimeoutExpired):
            pass
    return ""


def custom_flags_from_cfg(cfg: Optional[Mapping[str, Any]] = None) -> List[str]:
    raw = []
    if cfg is not None:
        raw = cfg.get("custom-flags") or cfg.get("custom_flags") or []
    if not raw:
        return list(DEFAULT_CUSTOM_FLAGS)
    out: List[str] = []
    for item in raw if isinstance(raw, (list, tuple)) else [raw]:
        flag = _norm_flag(item)
        if flag and flag not in out:
            out.append(flag)
    for flag in DEFAULT_CUSTOM_FLAGS:
        if flag not in out:
            out.append(flag)
    return out


def translated_argv(canonical_flag: str, value: Any, mapping: Mapping[str, str]) -> str:
    dest = mapping.get(_norm_flag(canonical_flag)) or mapping.get(str(canonical_flag).strip())
    if not dest:
        return ""
    if value is True or value in (None, ""):
        return str(dest)
    if value is False:
        return ""
    return f"{dest} {value}".strip()


def threads_cli_token(name: str, threads: Any, record: Optional[Mapping[str, Any]] = None) -> str:
    mapping = flags_translate_for(name, record)
    return translated_argv("--threads", threads, mapping)


def apply_translated_flags(
    extra: str,
    *,
    name: str,
    record: Optional[Mapping[str, Any]] = None,
    canonical: Optional[Mapping[str, Any]] = None,
) -> str:
    """Append mapped CLI flags (e.g. ``--threads 8`` → ``-z 8`` for kaiju)."""
    mapping = flags_translate_for(name, record)
    parts = [str(extra or "").strip()]
    for src, value in (canonical or {}).items():
        if value in (None, False, ""):
            continue
        token = translated_argv(src, value, mapping)
        if not token:
            continue
        dest = token.split()[0]
        if dest and dest in " ".join(parts):
            continue
        parts.append(token)
    return " ".join(p for p in parts if p)


def parse_unknown_custom_flags(
    unknown: Sequence[str],
    custom_flags: Sequence[str],
) -> Tuple[Dict[str, str], List[str]]:
    """Split leftover argv into custom-flag values vs still-unknown tokens."""
    wanted = {_norm_flag(f) for f in custom_flags if f}
    values: Dict[str, str] = {}
    leftover: List[str] = []
    i = 0
    tokens = list(unknown)
    while i < len(tokens):
        tok = tokens[i]
        key = tok
        attached = ""
        if tok.startswith("--") and "=" in tok:
            key, _, attached = tok.partition("=")
        flag = _norm_flag(key) if str(tok).startswith("-") else ""
        if flag in wanted:
            if attached:
                values[flag] = attached
                i += 1
                continue
            if i + 1 < len(tokens) and not str(tokens[i + 1]).startswith("-"):
                values[flag] = tokens[i + 1]
                i += 2
                continue
            values[flag] = "true"
            i += 1
            continue
        leftover.append(tok)
        i += 1
    return values, leftover


def version_mismatch_rows(cfg: Optional[Mapping[str, Any]] = None) -> List[Tuple[str, str, str, str]]:
    """Return (name, found, tested_or_previous, lazy-install) for mismatched tools."""
    from samovar.paths import load_config

    data = cfg or load_config()
    raw = data.get("tools") if isinstance(data.get("tools"), dict) else {}
    rows: List[Tuple[str, str, str, str]] = []
    for key, value in raw.items():
        rec = parse_tool_record(value, str(key))
        bare, prev = split_tool_key(key)
        found = probe_tool_version(str((rec.get("exec") or {}).get("path") or ""), bare)
        tested = TESTED_TOOL_VERSIONS.get(bare) or prev
        if bare == "samovar":
            from samovar.version import get_version

            found = get_version()
            tested = get_version()
        if not found or not tested:
            continue
        if found.split(".")[:2] != tested.split(".")[:2] and found != tested:
            rows.append((bare, found, tested, str(rec.get("lazy-install") or "")))
    return rows


def load_custom_flags(cfg: Optional[Mapping[str, Any]] = None) -> List[str]:
    if cfg is not None:
        return custom_flags_from_cfg(cfg)
    try:
        from samovar.paths import load_config

        return custom_flags_from_cfg(load_config())
    except Exception:
        return list(DEFAULT_CUSTOM_FLAGS)


def parse_known_with_custom(parser, argv: Optional[Sequence[str]] = None):
    """``parse_known_args`` then accept ``custom-flags`` from the install config."""
    args, unknown = parser.parse_known_args(None if argv is None else list(argv))
    custom = load_custom_flags()
    values, leftover = parse_unknown_custom_flags(unknown, custom)
    setattr(args, "custom_cli_flags", values)
    if values.get("--threads") and not getattr(args, "threads", None):
        try:
            args.threads = int(values["--threads"])
        except (TypeError, ValueError):
            args.threads = values["--threads"]
    if values.get("--cores") and getattr(args, "cores", None) in (None, 1):
        try:
            args.cores = int(values["--cores"])
        except (TypeError, ValueError):
            pass
    return args, leftover
