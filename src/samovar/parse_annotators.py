"""Module for parsing and processing taxonomic annotation outputs from various tools.

This module provides functionality to read and process taxonomic annotation outputs
from tools like Kraken, Kaiju, and MetaPhlAn. It standardizes the output format
to allow for easy comparison and analysis of results from different tools.
"""

import pandas as pd
import os
import re
import sys
import hashlib
import tempfile
from typing import Dict, List, Optional
import sqlite3
import pickle
from pathlib import Path

from samovar.annotators_wrapper import get_annotator_instance
from samovar.taxonomy_engine import NCBITaxonomyParser

# `ete3` imports its webplugin unconditionally, and the webplugin imports the
# stdlib `cgi` module. Python 3.13 removed `cgi`, so importing ete3 fails
# unless we provide a tiny compatibility shim.
try:
    from ete3 import NCBITaxa
except ModuleNotFoundError as exc:
    if str(exc).endswith(": cgi") or getattr(exc, "name", None) == "cgi":
        # Minimal `cgi` module shim: ete3 only needs `cgi.parse_qs` at import
        # time (runtime usage is irrelevant for our CLI/library usage).
        import types
        import urllib.parse

        cgi_stub = types.ModuleType("cgi")
        cgi_stub.parse_qs = urllib.parse.parse_qs
        sys.modules["cgi"] = cgi_stub
        from ete3 import NCBITaxa
    else:
        raise

# Common ISS / test-genome prefixes when headers omit `taxid:`.
_TRUE_TAXID_PREFIXES = {
    "scer.fna": "4932",
    "ecoli.fna": "562",
    "hsap.fna": "9606",
    "phix.fna": "2886930",
}


def extract_true_taxid(seq_id: str, pattern: Optional[str] = None) -> str:
    """Extract a true taxID from a read/sequence identifier.

    Prefers an explicit `taxid:<digits>` token. The default lookbehind pattern
    `(?<=taxid:)[0-9]*` can match an empty string, and some yeast contigs in
    the test genomes were stored as `Scer.fna-NC_...` without `taxid:`.
    """
    seq_id = "" if seq_id is None else str(seq_id)
    if pattern:
        match = re.search(pattern, seq_id)
        if match:
            token = match.group(0)
            if token:
                return token

    match = re.search(r"taxid:([0-9]+)", seq_id, flags=re.IGNORECASE)
    if match:
        return match.group(1)

    lower = seq_id.lower()
    for prefix, taxid in _TRUE_TAXID_PREFIXES.items():
        if lower.startswith(prefix) or f"|{prefix}" in lower:
            return taxid
    return ""


def extract_read_type(seq_id: str) -> str:
    """``read_type:<token>`` from a read/sequence identifier (hybrid CAMISIM)."""
    text = "" if seq_id is None else str(seq_id)
    match = re.search(r"read_type:([A-Za-z0-9_+\-]+)", text, flags=re.IGNORECASE)
    return match.group(1).lower() if match else ""


ncbi = None
_taxonomy_parent_rank = None

# Metrics and plots default to genus: annotators often return strain/species/genus
# taxIDs that only agree after collapsing to a shared rank via NCBI lineage.
DEFAULT_TAX_RANK = "genus"

_RANK_ALIASES = {
    "sp": "species",
    "species": "species",
    "g": "genus",
    "genus": "genus",
    "genera": "genus",
    "fam": "family",
    "family": "family",
    "order": "order",
    "class": "class",
    "phylum": "phylum",
    "kingdom": "kingdom",
    "superkingdom": "superkingdom",
}

# Minimal built-in fallback for common test taxids if no taxonomy DB is available.
_FALLBACK_LINEAGE = {
    9606: {"species": 9606, "genus": 9605},
    9605: {"genus": 9605},
    511145: {"species": 562, "genus": 561},
    562: {"species": 562, "genus": 561},
    561: {"genus": 561},
    4932: {"species": 4932, "genus": 4930},
    4930: {"genus": 4930},
}


def normalize_tax_rank(rank: Optional[str]) -> Optional[str]:
    """Map CLI/user rank names to NCBI rank strings. None means keep exact taxIDs."""
    if rank is None:
        return None
    key = str(rank).strip().lower()
    if key in {"", "none", "exact", "taxid", "raw", "false", "off"}:
        return None
    return _RANK_ALIASES.get(key, key)


def canonical_taxid(value) -> str:
    """Normalize a table cell to a taxID string (unclassified -> 0)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "0"
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "na", "<na>"}:
        return "0"
    if text.endswith(".0") and text[:-2].lstrip("-").isdigit():
        text = text[:-2]
    match = re.search(r"(\d+)", text)
    return match.group(1) if match else text


def taxid_value_columns(columns) -> List[str]:
    """Annotator taxID columns plus ``true``; skip confidence fields."""
    out = []
    for col in columns:
        name = str(col).lower()
        if "confidence" in name:
            continue
        if name == "true" or name.startswith("taxid"):
            out.append(col)
    return out


def build_taxid_rank_map(taxids, rank: str) -> Dict[str, str]:
    """Map unique taxID strings to the ancestor (or self) at ``rank``."""
    rank_name = normalize_tax_rank(rank)
    mapping: Dict[str, str] = {}
    if rank_name is None:
        for taxid in taxids:
            key = canonical_taxid(taxid)
            mapping[str(taxid)] = key
            mapping[key] = key
        return mapping
    for taxid in taxids:
        raw = str(taxid) if taxid is not None else "0"
        key = canonical_taxid(taxid)
        if key == "0":
            mapping[raw] = "0"
            mapping[key] = "0"
            continue
        ranked = _resolve_taxid_by_rank(key, rank_name)
        mapping[raw] = ranked if ranked is not None else key
        mapping[key] = mapping[raw]
    return mapping


def rank_map_column(rank: str) -> str:
    """Column name for a rank translation table (genus -> ``genera_taxid``)."""
    resolved = normalize_tax_rank(rank) or str(rank).strip().lower() or "genus"
    if resolved == "genus":
        return "genera_taxid"
    return f"{resolved}_taxid"


def default_rank_map_path(rank: str) -> Path:
    return _cache_dir() / f"taxid_{rank_map_column(rank)}.tsv"


def read_taxid_rank_table(path) -> Dict[str, str]:
    """Load a ``taxid|{rank}_taxid`` translation table (pipe, tab, or csv)."""
    table = Path(path)
    if not table.exists():
        return {}
    mapping: Dict[str, str] = {}
    with table.open(encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "|" in line:
                parts = [p.strip() for p in line.split("|")]
            elif "\t" in line:
                parts = [p.strip() for p in line.split("\t")]
            else:
                parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                continue
            if i == 0 and parts[0].lower() in {"taxid", "tax_id"}:
                continue
            src = canonical_taxid(parts[0])
            dst = canonical_taxid(parts[1])
            mapping[src] = dst
            mapping[str(parts[0]).strip()] = dst
    return mapping


def write_taxid_rank_table(path, mapping: Dict[str, str], rank: str) -> None:
    """Write ``taxid|{rank}_taxid`` (pipe-separated) for reuse across viz runs."""
    table = Path(path)
    table.parent.mkdir(parents=True, exist_ok=True)
    col = rank_map_column(rank)
    rows: Dict[str, str] = {}
    for key, value in mapping.items():
        src = canonical_taxid(key)
        if src == "0":
            continue
        rows[src] = canonical_taxid(value)
    tmp = table.with_suffix(table.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write(f"taxid|{col}\n")
        for src in sorted(rows, key=lambda x: (len(x), x)):
            handle.write(f"{src}|{rows[src]}\n")
    tmp.replace(table)


def ensure_taxid_rank_map(taxids, rank: str, cache_path=None) -> Dict[str, str]:
    """Resolve missing taxIDs to ``rank`` and persist a translation table."""
    rank_name = normalize_tax_rank(rank) or rank
    paths = []
    if cache_path:
        paths.append(Path(cache_path))
    try:
        global_path = default_rank_map_path(rank_name)
        if not cache_path or Path(cache_path).resolve() != global_path.resolve():
            paths.append(global_path)
    except Exception:
        pass

    mapping: Dict[str, str] = {}
    for path in paths:
        mapping.update(read_taxid_rank_table(path))

    added = False
    for taxid in taxids:
        key = canonical_taxid(taxid)
        raw = str(taxid) if taxid is not None else "0"
        if key == "0":
            mapping[key] = "0"
            mapping[raw] = "0"
            continue
        if key in mapping:
            mapping[raw] = mapping[key]
            continue
        ranked = _resolve_taxid_by_rank(key, rank_name)
        mapping[key] = ranked if ranked is not None else key
        mapping[raw] = mapping[key]
        added = True

    requested = {canonical_taxid(t) for t in taxids}
    local_path = Path(cache_path).resolve() if cache_path else None
    for path in paths:
        if local_path and path.resolve() == local_path:
            subset = {
                key: value
                for key, value in mapping.items()
                if canonical_taxid(key) in requested
            }
            write_taxid_rank_table(path, subset, rank_name)
        elif added or not path.exists():
            write_taxid_rank_table(path, mapping, rank_name)
    return mapping


def default_name_map_path() -> Path:
    return _cache_dir() / "taxid_name_map.tsv"


_FALLBACK_NAMES = {
    9606: "Homo sapiens",
    9605: "Homo",
    562: "Escherichia coli",
    561: "Escherichia",
    511145: "Escherichia coli str. K-12 substr. MG1655",
    4932: "Saccharomyces cerevisiae",
    4930: "Saccharomyces",
    2886930: "Escherichia phage phiX174",
}


def ensure_taxid_name_map(taxids, cache_path=None) -> Dict[str, str]:
    """Resolve scientific names via ete3 and cache ``taxid|name`` rows."""
    paths = []
    if cache_path:
        paths.append(Path(cache_path))
    try:
        global_path = default_name_map_path()
        if not cache_path or Path(cache_path).resolve() != global_path.resolve():
            paths.append(global_path)
    except Exception:
        pass

    mapping: Dict[str, str] = {}
    for path in paths:
        mapping.update(_read_taxid_name_table(path))

    added = False
    missing = []
    for taxid in taxids:
        key = canonical_taxid(taxid)
        if key == "0":
            mapping[key] = "0"
            continue
        if key in mapping and mapping[key] not in {"", "0"}:
            continue
        missing.append(key)

    if missing:
        translated: Dict[int, str] = {}
        try:
            ids = [int(t) for t in missing]
            translated = _get_ncbi().get_taxid_translator(ids) or {}
        except Exception:
            translated = {}
        for key in missing:
            name = None
            try:
                name = translated.get(int(key))
            except (TypeError, ValueError):
                name = None
            if not name:
                try:
                    name = _FALLBACK_NAMES.get(int(key))
                except (TypeError, ValueError):
                    name = None
            mapping[key] = name if name else key
            added = True

    requested = {canonical_taxid(t) for t in taxids if canonical_taxid(t) != "0"}
    local_path = Path(cache_path).resolve() if cache_path else None
    for path in paths:
        if local_path and path.resolve() == local_path:
            subset = {k: v for k, v in mapping.items() if canonical_taxid(k) in requested}
            _write_taxid_name_table(path, subset)
        elif added or not path.exists():
            _write_taxid_name_table(path, mapping)
    return mapping


def _read_taxid_name_table(path) -> Dict[str, str]:
    table = Path(path)
    if not table.exists():
        return {}
    mapping: Dict[str, str] = {}
    with table.open(encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "|" not in line:
                continue
            src, _, name = line.partition("|")
            src = src.strip()
            name = name.strip()
            if i == 0 and src.lower() in {"taxid", "tax_id"}:
                continue
            if not src or not name:
                continue
            mapping[canonical_taxid(src)] = name
    return mapping


def _write_taxid_name_table(path, mapping: Dict[str, str]) -> None:
    table = Path(path)
    table.parent.mkdir(parents=True, exist_ok=True)
    rows = {}
    for key, value in mapping.items():
        src = canonical_taxid(key)
        if src == "0" or not value:
            continue
        rows[src] = str(value).replace("|", "/")
    tmp = table.with_suffix(table.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write("taxid|name\n")
        for src in sorted(rows, key=lambda x: (len(x), x)):
            handle.write(f"{src}|{rows[src]}\n")
    tmp.replace(table)


def remap_taxid_dataframe(
    df: pd.DataFrame,
    rank: str = DEFAULT_TAX_RANK,
    cache_path: Optional[str] = None,
) -> pd.DataFrame:
    """Rewrite taxID columns in a *copy* using a cached taxid→rank table.

    Annotation CSVs on disk are left unchanged. ``cache_path`` is a
    ``taxid|genera_taxid`` (or other rank) table written during visualization.
    """
    rank_name = normalize_tax_rank(rank)
    out = df.copy()
    if rank_name is None:
        return out
    cols = taxid_value_columns(out.columns)
    if not cols:
        return out
    unique = set()
    for col in cols:
        unique.update(out[col].dropna().unique())
    mapping = ensure_taxid_rank_map(unique, rank_name, cache_path)
    for col in cols:
        out[col] = out[col].map(
            lambda v, m=mapping: m.get(str(v), m.get(canonical_taxid(v), canonical_taxid(v)))
            if pd.notna(v)
            else v
        )
    return out


def _get_ncbi() -> "NCBITaxa":
    global ncbi
    if ncbi is None:
        ncbi = NCBITaxa()
    return ncbi


def _cache_dir() -> Path:
    """Cache directory for taxonomy maps (XDG or SAMOVAR_CACHE_DIR)."""
    override = os.environ.get("SAMOVAR_CACHE_DIR")
    if override:
        path = Path(override)
    else:
        xdg = os.environ.get("XDG_CACHE_HOME")
        path = Path(xdg) / "samovar" if xdg else Path.home() / ".cache" / "samovar"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _load_nodes_cache(nodes_path: str) -> Dict[int, tuple]:
    """Load parent/rank map from nodes.dmp with an atomic filesystem cache.

    Concurrent job-array writers reload after parsing and replace the pickle
    via a temp file so a partial write cannot be observed.
    """
    cache_dir = _cache_dir()
    stat = os.stat(nodes_path)
    digest = hashlib.sha256(
        f"{os.path.abspath(nodes_path)}:{int(stat.st_mtime)}:{stat.st_size}".encode()
    ).hexdigest()[:16]
    cache_file = cache_dir / f"nodes_{digest}.pkl"

    def _try_load():
        try:
            with open(cache_file, "rb") as fh:
                return pickle.load(fh)
        except FileNotFoundError:
            return None
        except (OSError, pickle.UnpicklingError, EOFError, AttributeError):
            try:
                cache_file.unlink()
            except OSError:
                pass
            return None

    loaded = _try_load()
    if loaded is not None:
        return loaded

    parser = NCBITaxonomyParser(nodes_path)
    tree = parser.tree

    loaded = _try_load()
    if loaded is not None:
        return loaded

    cache_dir.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=cache_file.name + ".", suffix=".tmp", dir=str(cache_dir)
    )
    try:
        with os.fdopen(fd, "wb") as fh:
            pickle.dump(tree, fh, protocol=pickle.HIGHEST_PROTOCOL)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_name, cache_file)
    except OSError:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        loaded = _try_load()
        if loaded is not None:
            return loaded
        return tree
    return tree


def _discover_nodes_path() -> Optional[str]:
    """Find nodes.dmp from the environment, not from machine-specific paths.

    Lookup order:
      1. ``SAMOVAR_NODES_DMP`` — exact file
      2. ``SAMOVAR_NODES_SEARCH`` — colon-separated files or directories to scan
    """
    env_path = os.environ.get("SAMOVAR_NODES_DMP")
    if env_path and os.path.exists(env_path):
        return env_path

    extra = os.environ.get("SAMOVAR_NODES_SEARCH", "")
    for raw in extra.split(":"):
        raw = raw.strip()
        if not raw:
            continue
        if os.path.isfile(raw) and os.path.basename(raw) == "nodes.dmp":
            return raw
        if os.path.isdir(raw):
            for candidate in (
                os.path.join(raw, "nodes.dmp"),
                os.path.join(raw, "taxonomy", "nodes.dmp"),
            ):
                if os.path.exists(candidate):
                    return candidate
    return None


def _get_taxid_by_rank_from_nodes(taxid: int, rank_name: str) -> Optional[str]:
    """Resolve rank using local nodes.dmp parser."""
    global _taxonomy_parent_rank
    if _taxonomy_parent_rank is None:
        nodes_path = _discover_nodes_path()
        if nodes_path is None:
            return None
        _taxonomy_parent_rank = _load_nodes_cache(nodes_path)

    tree = _taxonomy_parent_rank
    current = taxid
    visited = set()
    while current in tree and current not in visited:
        visited.add(current)
        parent, rank = tree[current]
        if rank == rank_name:
            return str(current)
        if parent == current:
            break
        current = parent
    return None


def _resolve_taxid_by_rank(taxid: str, rank_name: str) -> Optional[str]:
    """Resolve target rank with fast and robust fallback strategy."""
    rank_name = normalize_tax_rank(rank_name) or rank_name
    taxid = canonical_taxid(taxid)
    if taxid == "0":
        return "0"

    try:
        taxid_int = int(taxid)
    except (TypeError, ValueError):
        return None

    # 1) Try ete3-backed resolver.
    try:
        lineage = _get_ncbi().get_lineage(taxid_int)
        ranks = _get_ncbi().get_rank(lineage)
        for tid, rank in ranks.items():
            if rank == rank_name:
                return str(tid)
    except Exception:
        pass

    # 2) Try local nodes.dmp cache.
    local_hit = _get_taxid_by_rank_from_nodes(taxid_int, rank_name)
    if local_hit is not None:
        return local_hit

    # 3) Built-in tiny fallback for deterministic unit tests/common IDs.
    fallback = _FALLBACK_LINEAGE.get(taxid_int, {}).get(rank_name)
    if fallback is not None:
        return str(fallback)

    return None


def resolve_metaphlan_db_file(
    db_path: str, db_name: Optional[str] = None
) -> str:
    """Resolve a MetaPhlAn SQLite mapping file from a path or directory.

    ``db_path`` may be the ``.db`` file itself. Otherwise ``db_name`` (config)
    or the sole ``*.db`` under the directory is used.
    """
    if os.path.isfile(db_path):
        return db_path
    if db_name:
        db_file = os.path.join(db_path, db_name)
        if not os.path.exists(db_file):
            raise FileNotFoundError(f"MetaPhlAn database not found at {db_file}")
        return db_file
    if not os.path.isdir(db_path):
        raise FileNotFoundError(f"MetaPhlAn database not found at {db_path}")
    matches = sorted(Path(db_path).glob("*.db"))
    if not matches:
        raise FileNotFoundError(
            f"No *.db MetaPhlAn mapping file under {db_path}; "
            "pass db_path as the file or set db_name in the annotator config"
        )
    if len(matches) > 1:
        raise FileNotFoundError(
            f"Multiple *.db files under {db_path}; set db_name in the annotator config"
        )
    return str(matches[0])


def parse_metaphlan_db(
    db_path: str, db_name: Optional[str] = None
) -> Dict[str, str]:
    """Parse MetaPhlAn database to map reference IDs to NCBI taxIDs.
    
    Args:
        db_path: Path to MetaPhlAn database directory or SQLite file
        db_name: Optional filename inside ``db_path``
        
    Returns:
        Dictionary mapping MetaPhlAn reference IDs to NCBI taxIDs
    """
    db_file = resolve_metaphlan_db_file(db_path, db_name=db_name)
        
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # Query to get mapping between reference IDs and taxIDs
    cursor.execute("SELECT ref_id, tax_id FROM mpa_species_map")
    mapping = {row[0]: row[1] for row in cursor.fetchall()}
    
    conn.close()
    return mapping


def _read_table_or_empty(file_path: str, columns: List[str]) -> pd.DataFrame:
    """Read a whitespace table, returning an empty frame for empty files."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return pd.DataFrame(columns=columns)
    try:
        df = pd.read_table(file_path, header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns)
    if df.empty:
        return pd.DataFrame(columns=columns)
    df.columns = columns[: df.shape[1]]
    return df


def read_kaiju_raw(file_path: str, db_path: Optional[str] = None) -> pd.DataFrame:
    """Read raw Kaiju output file.
    
    Args:
        file_path: Path to Kaiju output file
        db_path: Optional path to Kaiju database file
        
    Returns:
        DataFrame with columns: classified, seq, score, taxID, N
    """
    df = _read_table_or_empty(file_path, ["classified", "seq", "taxID"])
    return df


def read_kraken1_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken1 output file.
    
    Args:
        file_path: Path to Kraken1 output file
        
    Returns:
        DataFrame with columns: classified, seq, taxID, length, k-mer
    """
    df = _read_table_or_empty(
        file_path, ["classified", "seq", "taxID", "length", "k-mer"]
    )
    return df


def read_kraken2_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken2 output file.
    
    Args:
        file_path: Path to Kraken2 output file
        
    Returns:
        DataFrame with columns: classified, seq, taxa, length, k-mer, taxID
    """
    df = _read_table_or_empty(
        file_path, ["classified", "seq", "taxa", "length", "k-mer"]
    )
    if df.empty:
        df["taxID"] = []
        return df
    df["length"] = [re.sub(r"\|.*", "", i) for i in df["length"]]
    df["taxID"] = [re.search(r'(?<=taxid )[0-9]*', i).group(0) for i in df["taxa"]]
    df["taxa"] = [re.sub(r' \(.*', "", i) for i in df["taxa"]]
    return df


def read_krakenu_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken-Unique output file.
    
    Args:
        file_path: Path to Kraken-Unique output file
        
    Returns:
        DataFrame with columns: classified, seq, taxID, length, k-mer
    """
    df = _read_table_or_empty(
        file_path, ["classified", "seq", "taxID", "length", "k-mer"]
    )
    return df


def read_metaphlan_raw(file_path: str, db_path: Optional[str] = None) -> pd.DataFrame:
    """Read raw MetaPhlAn output file.
    
    Args:
        file_path: Path to MetaPhlAn output file
        db_path: Optional path to MetaPhlAn database directory
        
    Returns:
        DataFrame with columns: seq, taxID
    """
    df = _read_table_or_empty(file_path, ["seq", "taxID"])
    if df.empty:
        return df
    
    if db_path is not None:
        # Parse database and map reference IDs to taxIDs
        db_mapping = parse_metaphlan_db(db_path)
        
        # Extract reference ID from the taxID field
        df["ref_id"] = df["taxID"].apply(lambda x: re.search(r'M\d+-c\d+', x).group(0) if re.search(r'M\d+-c\d+', x) else None)
        
        # Map reference IDs to taxIDs
        df["taxID"] = df["ref_id"].map(db_mapping)
        df = df.drop("ref_id", axis=1)
    
    return df


def read_custom_raw(file_path: str) -> pd.DataFrame:
    """Read custom annotator output (seq, taxID)."""
    return _read_table_or_empty(file_path, ["seq", "taxID"])


# Dictionary mapping tool names to their respective read functions
READ_FUNCTIONS = {
    "kraken": read_kraken1_raw,
    "kraken1": read_kraken1_raw,
    "kraken2": read_kraken2_raw,
    "krakenunique": read_krakenu_raw,
    "krakenuniq": read_krakenu_raw,
    "krakenu": read_krakenu_raw,
    "metaphlan": read_metaphlan_raw,
    "mpa": read_metaphlan_raw,
    "mp4": read_metaphlan_raw,
    "kaiju": read_kaiju_raw,
    "dummy": read_custom_raw,
    "dummy9606": read_custom_raw,
    "constant9606": read_custom_raw,
    "constant": read_custom_raw,
    "random": read_custom_raw,
    "custom": read_custom_raw,
    "centrifuge": read_custom_raw,
    "metauto": read_custom_raw,
    "assembly_hybrid": read_custom_raw,
}


def read_tool_output(file_path: str, tool_type: str) -> pd.DataFrame:
    """Read an annotator file via native parsers or the OOP wrapper."""
    reader = READ_FUNCTIONS.get(tool_type)
    if reader is not None:
        return reader(file_path)
    annotator = get_annotator_instance(tool_type, run_config={}, config={})
    return annotator.parse_output(file_path)


def read_annotation(file_path_type: Dict[str, str], trimmed: bool = True) -> pd.DataFrame:
    """Read and combine multiple annotation files.
    
    Args:
        file_path_type: Dictionary mapping file paths to their tool types
        trimmed: Whether to only keep taxID columns
        
    Returns:
        Combined DataFrame with annotations from all tools
    """
    res = pd.DataFrame()
    for path, tool_type in file_path_type.items():
        df = read_tool_output(path, tool_type).set_index("seq")
        if trimmed:
            df = df[["taxID"]]
        tool_name = tool_type.split("_")[-1].split(".")[0]
        df.columns = [f"{col}_{tool_name}" if col != "seq" else col for col in df.columns]
        res = pd.concat([res, df], axis=1)
    return res


class Annotation:
    """Class for handling taxonomic annotations from multiple tools.
    
    This class provides functionality to read, process, and analyze taxonomic
    annotations from various tools. It supports operations like expanding
    taxonomic ranks, comparing annotations, and exporting results.
    """
    
    def __init__(self, file_path_type: Dict[str, str], get_true_annotation: Optional[str] = None):
        """Initialize Annotation object.
        
        Args:
            file_path_type: Dictionary mapping file paths to their tool types
            get_true_annotation: Regex pattern to extract true annotation from sequence IDs
        """
        self.id = 0
        self.DataFrame = pd.DataFrame()
        
        # Read and combine all annotation files
        for path, tool_type in file_path_type.items():
            df = read_tool_output(path, tool_type).set_index("seq").astype({"taxID": 'string'})
            df.columns = [f"{col}_{tool_type}_{self.id}" for col in df.columns]
            
            # Check for duplicate indices
            if df.index.duplicated().any():
                df = df[~df.index.duplicated(keep='first')]
            
            try:
                self.DataFrame = pd.concat([self.DataFrame, df], axis=1)
            except:
                raise ValueError(f"Error concatenating {path}")
            self.id += 1

        self.DataFrame = self.DataFrame.fillna("0")
        self.DataFrame["read_type"] = [
            extract_read_type(seq_id) for seq_id in self.DataFrame.index
        ]

        # Extract true annotations if pattern provided
        self.true_annotation = []
        if get_true_annotation is not None:
            self.true_annotation = [
                extract_true_taxid(seq_id, get_true_annotation)
                for seq_id in self.DataFrame.index
            ]
            print("True annotations extracted")

        self.DataFrame["read_type"] = [
            extract_read_type(seq_id) for seq_id in self.DataFrame.index
        ]

        # Get unique annotations and ranks
        set_columns = []
        for name, column in self.tr().items():
            set_columns += list(set(column))

        self.annotation_list = self.list2set(set_columns)
        self.true_annotation_list = self.list2set(self.true_annotation)
        self.rank_list = self.list2set([*self.annotation_list, *self.true_annotation_list])

    def true_annotation_unique(self) -> set:
        """Get unique true annotations."""
        return set(self.true_annotation)

    def true_annotation_rename(self, change_dict: Dict[str, str]):
        """Rename true annotations using provided mapping.
        
        Args:
            change_dict: Dictionary mapping old names to new names
        """
        self.true_annotation = [change_dict.get(TA) if TA in change_dict else "" for TA in self.true_annotation]

    def rank_annotation(self, rank: str = DEFAULT_TAX_RANK) -> 'RankAnnotation':
        """Get annotations at specified taxonomic rank.
        
        Args:
            rank: Taxonomic rank to get annotations for
            
        Returns:
            RankAnnotation object with annotations at specified rank
        """
        rank_list = [self.rank(j, rank) for j in self.rank_list]
        rank_dict = dict(zip(self.rank_list, rank_list))
        return RankAnnotation(rank).make(self.full(), rank_dict)

    def expand_annotation(self, rank: List[str] = None) -> 'ExpandAnnotation':
        """Expand annotations to multiple taxonomic ranks.
        
        Args:
            rank: List of taxonomic ranks to expand to
            
        Returns:
            ExpandAnnotation object with annotations at multiple ranks
        """
        if rank is None:
            rank = [DEFAULT_TAX_RANK]
        full_rank_annotation = ExpandAnnotation()
        for i in rank:
            full_rank_annotation.add(self.rank_annotation(i))
        return full_rank_annotation

    def correct_annotations(self, rank: str = DEFAULT_TAX_RANK) -> pd.DataFrame:
        """Get correct annotations at specified rank.
        
        Args:
            rank: Taxonomic rank to check
            
        Returns:
            DataFrame with counts of correct annotations
        """
        return pd.DataFrame(self.rank_annotation(rank).correct_annotation().annotation.value_counts())

    def full(self) -> pd.DataFrame:
        """Get full annotation DataFrame including true annotations."""
        tmp = self.tr()
        tmp["true"] = self.true_annotation
        return tmp

    def correct_level(self, level: str = DEFAULT_TAX_RANK) -> None:
        """Correct annotations to specified taxonomic level.
        
        Args:
            level: Taxonomic level to correct to (default genus). Use
                ``none``/``exact`` to skip remapping.
            
        This method:
        1. Gets unique taxIDs from the DataFrame (and ``true`` / true_annotation)
        2. Maps them to the specified taxonomic level via ete3 NCBITaxa
        3. Replaces original taxIDs with their corrected versions
        """
        self.DataFrame = remap_taxid_dataframe(self.DataFrame, rank=level)
        if self.true_annotation:
            mapping = build_taxid_rank_map(self.true_annotation, level)
            self.true_annotation = [
                mapping.get(str(ta), canonical_taxid(ta)) for ta in self.true_annotation
            ]
            self.true_annotation_list = self.list2set(self.true_annotation)

    def tr(self) -> pd.DataFrame:
        """Get DataFrame with only taxID columns."""
        return self.DataFrame.copy().filter(regex="taxID.*")

    @staticmethod
    def list2set(a: List) -> List[str]:
        """Convert list to set of strings."""
        return list(set([str(i) for i in a]))

    @staticmethod
    def list2rank(a: List, at_rank: str) -> List:
        """Convert list of taxIDs to nodes at specified rank.
        
        Args:
            a: List of taxIDs
            at_rank: Target taxonomic rank
            
        Returns:
            List of taxIDs at specified rank
        """
        b = Annotation.list2set(a)
        return [Annotation.rank(i, at_rank) for i in b]

    @staticmethod
    def rank(j: str, i: str) -> Optional[str]:
        """Get taxID at specified rank.
        
        Args:
            j: Input taxID
            i: Target rank
            
        Returns:
            taxID at specified rank or None if not found
        """
        if j == "0":
            return "0"
        try:
            return _resolve_taxid_by_rank(j, i)
        except Exception:
            return None

    def export(self, file: Optional[str] = None) -> pd.DataFrame:
        """Export annotations to CSV file.
        
        Args:
            file: Optional path to save CSV file
            
        Returns:
            DataFrame with exported annotations
        """
        df_return = self.DataFrame.loc[:, [col for col in self.DataFrame if col.startswith('taxID')]]
        lencol = [col for col in self.DataFrame if col.startswith('len')]
        if len(lencol) > 0:
            df_return['length'] = self.DataFrame.loc[:, lencol[0]].to_list()
        df_return["true"] = self.true_annotation
        if file is not None:
            df_return.to_csv(file)
        return df_return


class RankAnnotation:
    """Class for handling annotations at a specific taxonomic rank."""
    
    def __init__(self, rank: str):
        """Initialize RankAnnotation object.
        
        Args:
            rank: Taxonomic rank
        """
        self.rank = rank
        self.annotation = pd.DataFrame()

    def add(self, name: str, annotation: List):
        """Add annotation column.
        
        Args:
            name: Column name
            annotation: List of annotations
        """
        self.annotation[name] = annotation

    def make(self, annotation: pd.DataFrame, rank_dict: Dict[str, str]) -> 'RankAnnotation':
        """Create RankAnnotation from DataFrame.
        
        Args:
            annotation: Input DataFrame
            rank_dict: Dictionary mapping taxIDs to ranks
            
        Returns:
            RankAnnotation object
        """
        for name, column in annotation.items():
            self.add(str(name), [rank_dict.get(TA) for TA in column])
        self.reindex(annotation.index)
        return self

    def reindex(self, index):
        """Set DataFrame index."""
        self.annotation.index = index

    def y(self) -> pd.Series:
        """Get true annotations."""
        return self.annotation["true"]

    def x(self) -> pd.DataFrame:
        """Get predicted annotations."""
        return self.annotation.copy().drop("true", axis=1)

    def correct_annotation(self) -> 'RankAnnotation':
        """Get correct annotations."""
        tmp = RankAnnotation(self.rank)
        for name, column in self.x().items():
            tmp.add(str(name), pd.DataFrame(column == self.y()))
        tmp.reindex(self.annotation.index)
        return tmp


class ExpandAnnotation:
    """Class for handling annotations at multiple taxonomic ranks."""
    
    def __init__(self):
        """Initialize ExpandAnnotation object."""
        self.rank_annotation = {}

    def add(self, rank_annotation: RankAnnotation):
        """Add RankAnnotation object.
        
        Args:
            rank_annotation: RankAnnotation object to add
        """
        self.rank_annotation[rank_annotation.rank] = rank_annotation

    def get(self, rank: str) -> RankAnnotation:
        """Get RankAnnotation for specified rank.
        
        Args:
            rank: Taxonomic rank
            
        Returns:
            RankAnnotation object
        """
        return self.rank_annotation[rank]

def match_annotation(annotation_name:str) -> Optional[str]:
    """Match annotation name to tool name.

    Recognizes native `*.<tool>.out` files and custom outputs
    `*.custom_<tool>.out` from CustomAnnotator / ConstantTaxidAnnotator.
    """
    annotation_name = os.path.basename(annotation_name)
    if not annotation_name.endswith(".out"):
        print("Skipping", annotation_name)
        return None

    arg = annotation_name.split(".")[-2]
    if arg.startswith("custom_"):
        return arg[len("custom_"):]

    aliases = {
        "kraken": "kraken",
        "kraken1": "kraken1",
        "kraken2": "kraken2",
        "krakenuniq": "krakenuniq",
        "krakenu": "krakenu",
        "metaphlan": "metaphlan",
        "metaphlan4": "metaphlan",
        "mpa": "metaphlan",
        "mp4": "metaphlan",
        "kaiju": "kaiju",
        "dummy": "dummy",
        "dummy9606": "dummy9606",
        "constant9606": "constant9606",
        "constant": "constant",
        "random": "random",
        "centrifuge": "centrifuge",
        "metauto": "metauto",
        "assembly_hybrid": "assembly_hybrid",
    }
    if arg in aliases:
        return aliases[arg]
    # Unknown .out files are treated as custom tools (two-column seq/taxID).
    return arg