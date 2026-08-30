"""NCBI vs GTDB taxonomy for ``samovar convert``.

Annotation tables still carry taxon tokens (NCBI taxids, GTDB names, or
genome accessions). Convert walks the matching dump:

* **NCBI** — ``nodes.dmp`` / ``names.dmp`` (``genomes.taxdump`` / ``SAMOVAR_TAXDUMP``)
* **GTDB** — gtdb-taxdump (same dmp files, rank ``domain``) or official
  ``bac120_taxonomy.tsv`` / ``ar53_taxonomy.tsv`` (and optional ``taxid.map``)

CAMI ``@TaxonomyID`` and ``@Ranks`` follow the chosen system.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

PathLike = Union[str, Path]
TaxTree = Dict[int, Tuple[int, str]]

TAXONOMY_NCBI = "ncbi"
TAXONOMY_GTDB = "gtdb"
TAXONOMY_ALIASES = {
    "ncbi": TAXONOMY_NCBI,
    "ncbitax": TAXONOMY_NCBI,
    "ncbitaxonomy": TAXONOMY_NCBI,
    "taxdump": TAXONOMY_NCBI,
    "gtdb": TAXONOMY_GTDB,
    "gtdb_tk": TAXONOMY_GTDB,
    "gtdbtk": TAXONOMY_GTDB,
    "gtdbtaxonomy": TAXONOMY_GTDB,
    "gtdb_taxdump": TAXONOMY_GTDB,
}

NCBI_CAMI_RANKS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
GTDB_CAMI_RANKS = (
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
GTDB_RANK_PREFIX = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
GTDB_TSV_NAMES = (
    "bac120_taxonomy.tsv",
    "ar53_taxonomy.tsv",
    "ar122_taxonomy.tsv",
    "taxonomy.tsv",
)


class UnknownTaxonomyError(ValueError):
    """``--taxonomy`` is not NCBI or GTDB."""


def normalize_taxonomy(name: Any) -> str:
    key = str(name or "").strip().lower().replace("-", "_").replace(" ", "")
    if not key:
        return TAXONOMY_NCBI
    system = TAXONOMY_ALIASES.get(key, key)
    if system not in {TAXONOMY_NCBI, TAXONOMY_GTDB}:
        raise UnknownTaxonomyError(
            f"unknown taxonomy {name!r}; use ncbi or gtdb"
        )
    return system


def cami_ranks_for(system: str) -> Tuple[str, ...]:
    return GTDB_CAMI_RANKS if normalize_taxonomy(system) == TAXONOMY_GTDB else NCBI_CAMI_RANKS


def gtdb_dir(cfg: Optional[Dict[str, Any]] = None) -> Path:
    """Install-level GTDB dump directory (dmp files or taxonomy TSV)."""
    env = (
        os.environ.get("SAMOVAR_GTDB", "").strip()
        or os.environ.get("SAMOVAR_GTDB_TAXDUMP", "").strip()
    )
    if env:
        return Path(env).expanduser()
    if cfg is None:
        try:
            from samovar.paths import load_config

            cfg = load_config()
        except Exception:
            cfg = {}
    from samovar.genome_index import samovar_database_dir
    from samovar.main_config import genomes_block

    block = genomes_block(cfg or {}) if cfg else {}
    if isinstance(block, dict):
        for key in ("gtdb", "gtdb_taxdump", "taxdump_gtdb"):
            raw = str(block.get(key) or "").strip()
            if raw:
                return Path(raw).expanduser()
    try:
        from samovar.db_spec import lookup_database_record

        rec = lookup_database_record(cfg, "taxdump", "gtdb")
        if rec and rec.get("path"):
            return Path(str(rec["path"])).expanduser()
    except Exception:
        pass
    return samovar_database_dir(cfg) / "gtdb"


def taxonomy_root(system: str, taxdump: Optional[PathLike] = None) -> Path:
    if taxdump is not None:
        return Path(taxdump).expanduser()
    if normalize_taxonomy(system) == TAXONOMY_GTDB:
        return gtdb_dir()
    from samovar.taxdump import taxdump_dir

    return taxdump_dir()


def gtdb_taxonomy_tables(root: PathLike) -> List[Path]:
    base = Path(root).expanduser()
    found: List[Path] = []
    seen = set()

    def _add(path: Path) -> None:
        try:
            resolved = path.resolve()
        except OSError:
            resolved = path
        if resolved in seen or not path.is_file():
            return
        seen.add(resolved)
        found.append(path)

    if base.is_file():
        _add(base)
        return found
    for name in GTDB_TSV_NAMES:
        _add(base / name)
        _add(base / "taxonomy" / name)
    for pattern in ("*taxonomy*.tsv", "*_taxonomy.tsv"):
        for path in sorted(base.glob(pattern)):
            _add(path)
        nested = base / "taxonomy"
        if nested.is_dir():
            for path in sorted(nested.glob(pattern.replace("*taxonomy*.tsv", "*.tsv"))):
                _add(path)
            for path in sorted(nested.glob("*.tsv")):
                if "metadata" in path.name.lower():
                    continue
                _add(path)
    return found


def parse_taxid_map(path: PathLike) -> Dict[str, int]:
    """Accession (or token) → integer taxid (gtdb-taxdump ``taxid.map``)."""
    aliases: Dict[str, int] = {}
    with open(path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            parts = text.split("\t") if "\t" in text else text.split()
            if len(parts) < 2:
                continue
            token, raw_id = parts[0].strip(), parts[1].strip()
            try:
                taxid = int(float(raw_id))
            except (TypeError, ValueError):
                continue
            if not token:
                continue
            aliases[token] = taxid
            aliases[token.lower()] = taxid
            stripped = token
            for prefix in ("RS_", "GB_", "RS ", "GB "):
                if stripped.upper().startswith(prefix):
                    stripped = stripped[len(prefix) :]
                    aliases[stripped] = taxid
                    aliases[stripped.lower()] = taxid
    return aliases


def _stable_node_id(rank: str, name: str, used: set) -> int:
    digest = hashlib.sha1(f"{rank}|{name.lower()}".encode("utf-8")).hexdigest()
    n = int(digest[:8], 16) & 0x7FFFFFFF
    if n < 2:
        n = 2
    while n in used:
        n = n + 1 if n < 0x7FFFFFFF else 2
    used.add(n)
    return n


def _split_gtdb_lineage(lineage: str) -> List[Tuple[str, str, str]]:
    slots: List[Tuple[str, str, str]] = []
    for piece in str(lineage or "").replace("|", ";").split(";"):
        token = piece.strip()
        if not token or token in {".", "n/a", "NA"}:
            continue
        rank = ""
        name = token
        if "__" in token:
            prefix, _, rest = token.partition("__")
            rank = GTDB_RANK_PREFIX.get(prefix.lower(), "")
            name = rest.strip() or token
        if not name:
            continue
        slots.append((rank or "no rank", name, token))
    return slots


def parse_gtdb_taxonomy_files(
    paths: Sequence[PathLike],
) -> Tuple[TaxTree, Dict[int, str], Dict[str, int]]:
    """Build parent/rank tree + name/accession aliases from GTDB taxonomy TSV."""
    tree: TaxTree = {1: (1, "no rank")}
    names: Dict[int, str] = {1: "root"}
    aliases: Dict[str, int] = {"root": 1, "1": 1}
    used = {1}

    def _node(rank: str, name: str, token: str) -> int:
        key = f"{rank}|{name.lower()}"
        existing = aliases.get(f"node:{key}")
        if existing:
            taxid = existing
        else:
            taxid = _stable_node_id(rank, name, used)
            aliases[f"node:{key}"] = taxid
            names[taxid] = name
        aliases[token] = taxid
        aliases[token.lower()] = taxid
        aliases[name.lower()] = taxid
        aliases[str(taxid)] = taxid
        return taxid

    for path in paths:
        for accession, lineage, ncbi_taxid in _iter_gtdb_rows(path):
            parent = 1
            leaf = 1
            for rank, name, token in _split_gtdb_lineage(lineage):
                taxid = _node(rank, name, token)
                if taxid not in tree:
                    tree[taxid] = (parent, rank)
                parent = taxid
                leaf = taxid
            if accession:
                aliases[accession] = leaf
                aliases[accession.lower()] = leaf
            if ncbi_taxid:
                aliases[str(ncbi_taxid)] = leaf
    return tree, names, aliases


def _iter_gtdb_rows(path: PathLike) -> Iterable[Tuple[str, str, Optional[int]]]:
    file_path = Path(path)
    with open(file_path, encoding="utf-8", errors="ignore", newline="") as handle:
        sample = handle.readline()
        if not sample:
            return
        body = sample + handle.read()
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
    except csv.Error:
        dialect = csv.excel_tab
    rows = csv.reader(body.splitlines(), dialect)
    header = next(rows, None)
    if header is None:
        return
    lower = [c.strip().lower() for c in header]
    if "gtdb_taxonomy" in lower:
        idx = {name: i for i, name in enumerate(lower)}
        tax_i = idx["gtdb_taxonomy"]
        acc_i = idx.get("accession", idx.get("user_genome", idx.get("genome", 0)))
        ncbi_i = idx.get("ncbi_taxid", idx.get("ncbi_taxonomy_id"))
        for row in rows:
            if len(row) <= tax_i:
                continue
            lineage = row[tax_i].strip()
            accession = row[acc_i].strip() if acc_i < len(row) else ""
            ncbi: Optional[int] = None
            if ncbi_i is not None and ncbi_i < len(row) and row[ncbi_i].strip():
                try:
                    ncbi = int(float(row[ncbi_i].strip()))
                except (TypeError, ValueError):
                    ncbi = None
            yield accession, lineage, ncbi
        return
    lineage0 = header[1].strip() if len(header) > 1 else ""
    if "d__" in lineage0 or "p__" in lineage0:
        yield header[0].strip(), header[1].strip(), None
    for row in rows:
        if len(row) < 2:
            continue
        yield row[0].strip(), row[1].strip(), None


def cami_rank_for(system: str, rank: str) -> Optional[str]:
    """Map a dump rank string onto a CAMI ``RANKS`` token."""
    token = str(rank or "").strip().lower()
    ranks = cami_ranks_for(system)
    if normalize_taxonomy(system) == TAXONOMY_GTDB:
        if token in {"superkingdom", "domain"}:
            return "domain"
    elif token in {"superkingdom", "domain"}:
        return "superkingdom"
    if token in ranks:
        return token
    return None
