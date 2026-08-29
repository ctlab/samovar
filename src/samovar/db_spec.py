"""Canonical ``databases.<tool>.<name:version>`` records.

On disk a database is an object nested under its annotator, keyed like tools::

    "databases": {
      "kraken2": {
        "standard_8GB:2025oct": {
          "path": "/path/to/db",
          "flags": "",
          "lazy-download": "curl -L … && tar -xzf …",
          "url": "https://…"
        }
      }
    }

Legacy ``[[name, path, flags], …]`` rows and a single ``{name, path, flags}``
dict still parse. ``flags`` are native CLI tokens for the annotator that uses
this DB (merged into ``AnnotatorConfig.extra`` at prepare).
"""

from __future__ import annotations

import shlex
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from samovar.tool_spec import join_tool_key, split_tool_key

SCHEMA_KEYS = {
    "path",
    "database_path",
    "name",
    "database_name",
    "flags",
    "database_flags",
    "lazy-download",
    "lazy_download",
    "lazy-install",
    "version",
    "url",
    "type",
    "tool",
}

# (tool, name) or (tool, name, version) → official archive used to rebuild the index.
DEFAULT_LAZY_DOWNLOAD: Dict[Tuple[str, ...], str] = {
    (
        "kraken2",
        "standard_8GB",
    ): "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20251015.tar.gz",
    (
        "kraken2",
        "virus",
    ): "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz",
    (
        "kraken2",
        "pracken",
    ): "https://genome-idx.s3.amazonaws.com/kraken/k2_NCBI_reference_20251007.tar.gz",
    (
        "kaiju",
        "fungi",
    ): "https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_fungi_2024-08-16.tgz",
    (
        "kaiju",
        "refseq",
    ): "https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_2024-08-14.tgz",
    (
        "centrifuge",
        "p_compressed+h+v",
    ): "https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz",
    (
        "kraken",
        "minikraken_4GB",
    ): "https://ccb.jhu.edu/data/minikraken/minikraken_20171013_4GB.tgz",
    (
        "kraken",
        "minikraken_8GB",
    ): "https://ccb.jhu.edu/data/minikraken/minikraken_20171019_8GB.tgz",
    (
        "krakenuniq",
        "microbial",
    ): "https://genome-idx.s3.amazonaws.com/kraken/kuniq_microbialdb_minus_kdb.20230808.tgz",
    (
        "metaphlan",
        "toy",
    ): "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.tar",
    (
        "metaphlan",
        "jan25",
    ): "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503.tar",
    (
        "qiime",
        "silva-138-99",
    ): "https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza",
    (
        "qiime",
        "silva-138-99",
        "2024.2",
    ): "https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza",
    (
        "qiime",
        "silva-138-99",
        "2024.10",
    ): "https://data.qiime2.org/2024.10/common/silva-138-99-nb-classifier.qza",
    (
        "qiime",
        "unite-dynamic",
    ): "https://github.com/colinbrislawn/unite-train/releases/download/v10.0-v04.04.2024-qiime2-2024.2/unite_ver10_dynamic_04.04.2024-Q2-2024.2.qza",
    (
        "qiime",
        "unite-99",
    ): "https://github.com/colinbrislawn/unite-train/releases/download/v10.0-v04.04.2024-qiime2-2024.2/unite_ver10_99_04.04.2024-Q2-2024.2.qza",
    (
        "taxdump",
        "ncbi",
    ): "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
}


def split_db_key(key: str) -> Tuple[str, str]:
    """Inner key ``name:version`` → (name, version)."""
    return split_tool_key(key)


def join_db_key(name: str, version: str = "") -> str:
    return join_tool_key(name, version)


def parse_tool_and_name(raw: str, tool: str = "") -> Tuple[str, str]:
    """``kraken2:standard_8GB`` or ``-n standard_8GB --tool kraken2``."""
    text = str(raw or "").strip()
    hinted = str(tool or "").strip()
    if not text:
        return hinted, ""
    if hinted:
        return hinted, text
    if ":" in text:
        left, _, right = text.partition(":")
        if right:
            return left.strip(), right.strip()
    return "", text


def _archive_basename(url: str) -> str:
    path = str(url or "").strip().split("?", 1)[0].rstrip("/")
    name = path.rsplit("/", 1)[-1] if path else "archive"
    return name or "archive"


def _tar_extract_flag(url: str) -> str:
    lower = str(url or "").lower()
    if lower.endswith(".tar") and not lower.endswith((".tar.gz", ".tar.bz2", ".tar.xz")):
        return "-xf"
    if lower.endswith((".tar.bz2", ".tbz2", ".tbz")):
        return "-xjf"
    if lower.endswith((".tar.xz", ".txz")):
        return "-xJf"
    return "-xzf"


def curl_tarball_recipe(url: str, *, extract: bool = True) -> str:
    quoted = shlex.quote(str(url).strip())
    archive = shlex.quote(_archive_basename(url))
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        'DEST="${PREFIX:-.}"',
        'mkdir -p "$DEST"',
        f'ARCHIVE="$DEST/{archive}"',
        f'curl -L --fail --retry 20 --retry-delay 30 -C - -o "$ARCHIVE" {quoted}',
    ]
    if extract:
        flag = _tar_extract_flag(url)
        lines.append(f'tar {flag} "$ARCHIVE" -C "$DEST"')
    return "\n".join(lines) + "\n"


def curl_file_recipe(url: str) -> str:
    quoted = shlex.quote(str(url).strip())
    fname = shlex.quote(_archive_basename(url))
    return (
        "#!/bin/bash\n"
        "set -euo pipefail\n"
        'DEST="${PREFIX:-.}"\n'
        "mkdir -p \"$DEST\"\n"
        f'curl -L --fail --retry 20 --retry-delay 30 -C - -o "$DEST/{fname}" {quoted}\n'
    )


def official_url_for(tool: str, name: str, version: str = "", url: str = "") -> str:
    href = str(url or "").strip()
    if href:
        return href
    bare = split_db_key(name)[0]
    tool_name = str(tool).strip()
    ver = str(version or "").strip()
    if ver:
        href = DEFAULT_LAZY_DOWNLOAD.get((tool_name, bare, ver), "")
    if not href:
        href = DEFAULT_LAZY_DOWNLOAD.get((tool_name, bare), "") or DEFAULT_LAZY_DOWNLOAD.get(
            (tool_name, str(name).strip()), ""
        )
    return href


def lazy_download_for(tool: str, name: str, version: str = "", url: str = "") -> str:
    href = official_url_for(tool, name, version, url)
    if not href:
        return ""
    lower = href.lower()
    if lower.endswith((".qza", ".qza.gz", ".fa.gz", ".fasta.gz", ".fna.gz")):
        return curl_file_recipe(href)
    return curl_tarball_recipe(href)


def _as_str(value: Any) -> str:
    return "" if value is None else str(value).strip()


def parse_database_record(
    value: Any,
    *,
    tool: str = "",
    name: str = "",
) -> Dict[str, Any]:
    """Normalize a databases.* value to an object record."""
    bare, ver = split_db_key(name)
    rec: Dict[str, Any] = {
        "tool": str(tool or "").strip(),
        "name": bare,
        "path": "",
        "flags": "",
        "lazy-download": "",
        "url": "",
        "type": "database",
        "_version": ver,
    }
    if value is None or value is False:
        return rec
    if isinstance(value, dict):
        rec["tool"] = _as_str(value.get("tool")) or rec["tool"]
        rec["name"] = (
            _as_str(value.get("name") or value.get("database_name")) or rec["name"]
        )
        rec["path"] = _as_str(value.get("path") or value.get("database_path"))
        rec["flags"] = _as_str(value.get("flags") or value.get("database_flags"))
        rec["lazy-download"] = _as_str(
            value.get("lazy-download")
            or value.get("lazy_download")
            or value.get("lazy-install")
        )
        rec["url"] = _as_str(value.get("url"))
        rec["type"] = _as_str(value.get("type")) or "database"
        rec["_version"] = _as_str(value.get("version") or value.get("_version")) or ver
        if not rec["name"]:
            rec["name"] = bare
        return rec
    if isinstance(value, (list, tuple)):
        parts = [_as_str(x) for x in value]
        while len(parts) < 3:
            parts.append("")
        rec["name"] = parts[0] or rec["name"]
        rec["path"] = parts[1]
        rec["flags"] = parts[2]
        return rec
    rec["path"] = _as_str(value)
    return rec


def record_to_row(rec: Mapping[str, Any]) -> List[str]:
    name = str(rec.get("name") or "").strip()
    ver = str(rec.get("_version") or rec.get("version") or "").strip()
    token = join_db_key(name, ver) if ver and ":" not in name else name
    return [token, str(rec.get("path") or "").strip(), str(rec.get("flags") or "").strip()]


def annotator_flags_for_db(tool: str, path: str, stored_flags: str = "") -> str:
    """Native CLI flags the annotator should see, including tool-specific DB path."""
    extra = str(stored_flags or "").strip()
    dest = str(path or "").strip()
    if not dest:
        return extra
    if tool in {"metaphlan", "metaphlan4"} and "--bowtie2db" not in extra:
        extra = f"--bowtie2db {shlex.quote(dest)} {extra}".strip()
    return extra


def looks_like_db_map(value: Mapping[str, Any]) -> bool:
    """True if ``value`` is ``{name: record, …}`` rather than a single record."""
    if not value:
        return False
    keys = set(value)
    if keys <= SCHEMA_KEYS:
        return False
    nested = [
        k
        for k, item in value.items()
        if isinstance(item, (dict, list, tuple)) and k not in SCHEMA_KEYS
    ]
    return bool(nested)


def iter_database_records(cfg: Optional[Mapping[str, Any]]) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """``{tool: {name:version → record}}`` (full objects)."""
    raw = (cfg or {}).get("databases")
    if not isinstance(raw, dict):
        return {}
    out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for tool, value in raw.items():
        tool_name = str(tool).strip()
        if not tool_name or tool_name.startswith("_"):
            continue
        grouped: Dict[str, Dict[str, Any]] = {}
        if isinstance(value, dict) and looks_like_db_map(value):
            for inner, item in value.items():
                if str(inner).startswith("_"):
                    continue
                rec = parse_database_record(item, tool=tool_name, name=str(inner))
                rec["tool"] = rec.get("tool") or tool_name
                key = join_db_key(rec["name"] or str(inner), rec.get("_version") or "")
                if not rec["name"]:
                    rec["name"] = split_db_key(str(inner))[0]
                grouped[key] = rec
        elif isinstance(value, dict):
            rec = parse_database_record(value, tool=tool_name, name=str(value.get("name") or tool_name))
            rec["tool"] = rec.get("tool") or tool_name
            key = join_db_key(rec["name"] or tool_name, rec.get("_version") or "")
            grouped[key] = rec
        elif isinstance(value, (list, tuple)):
            rows: Iterable[Any]
            if value and isinstance(value[0], (list, tuple, dict)):
                rows = value
            else:
                rows = [value]
            for item in rows:
                rec = parse_database_record(item, tool=tool_name)
                rec["tool"] = rec.get("tool") or tool_name
                if not rec["name"] and rec["path"]:
                    rec["name"] = Path(rec["path"]).name
                key = join_db_key(rec["name"], rec.get("_version") or "")
                if key:
                    grouped[key] = rec
        elif value:
            rec = parse_database_record(value, tool=tool_name)
            rec["tool"] = tool_name
            key = join_db_key(rec["name"] or tool_name, rec.get("_version") or "")
            grouped[key] = rec
        if grouped:
            out[tool_name] = grouped
    return out


def databases_to_rows(cfg: Optional[Mapping[str, Any]]) -> Dict[str, List[List[str]]]:
    """Legacy ``{tool: [[name, path, flags], …]}`` view."""
    out: Dict[str, List[List[str]]] = {}
    for tool, grouped in iter_database_records(cfg).items():
        rows = []
        for rec in grouped.values():
            row = record_to_row(rec)
            if any(row):
                rows.append(row)
        if rows:
            out[tool] = rows
    return out


def lookup_database_record(
    cfg: Optional[Mapping[str, Any]],
    tool: str,
    name: str,
) -> Optional[Dict[str, Any]]:
    grouped = iter_database_records(cfg).get(str(tool).strip()) or {}
    token = str(name or "").strip()
    if not token:
        return None
    if token in grouped:
        return grouped[token]
    matches: List[Dict[str, Any]] = []
    for key, rec in grouped.items():
        bare, ver = split_db_key(key)
        if key == token or rec.get("name") == token or bare == token:
            matches.append(rec)
        elif join_db_key(bare, ver) == token:
            matches.append(rec)
    if not matches:
        return None
    live = [rec for rec in matches if Path(str(rec.get("path") or "")).expanduser().exists()]
    return (live or matches)[-1]


def databases_for_disk(cfg: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    """Object form nested under each annotator."""
    out: Dict[str, Any] = {}
    for tool, grouped in iter_database_records(cfg).items():
        inner: Dict[str, Any] = {}
        for key, rec in grouped.items():
            version = str(rec.get("_version") or "").strip()
            url = official_url_for(
                tool, rec.get("name") or key, version, str(rec.get("url") or "")
            )
            disk = {
                "path": str(rec.get("path") or ""),
                "flags": str(rec.get("flags") or ""),
                "lazy-download": str(rec.get("lazy-download") or ""),
                "type": "database",
            }
            if url:
                disk["url"] = url
            if version:
                disk["version"] = version
            if not disk["lazy-download"]:
                disk["lazy-download"] = lazy_download_for(
                    tool, rec.get("name") or key, version, url
                )
            inner[join_db_key(str(rec.get("name") or key), version)] = disk
        if inner:
            out[tool] = inner
    return out
