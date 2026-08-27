"""Resolve a table taxid to an indexed (or NCBI) genome at a rank.

``--reannotation-level``:

- ``taxid`` / ``t``: exact catalog key, else NCBI for that taxid.
- ``species`` / ``s``: taxid, then species (up from strain, down from species).
- ``genus`` / ``g``, ``family`` / ``f``, ``any`` / ``a``: walk ranks coarser
  until the catalog has ≥1 hit at that rank; otherwise NCBI under that id.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from samovar.genome_index import (
    genome_data_map,
    is_test_taxid,
    resolve_indexed_path,
)
from samovar.parse_annotators import (
    _resolve_taxid_by_rank_exact,
    canonical_taxid,
)
from samovar.seqio import find_genome_file
from samovar.taxdump import canonical_ncbi_taxid

logger = logging.getLogger(__name__)

REANNOTATION_ALIASES = {
    "t": "taxid",
    "taxid": "taxid",
    "exact": "taxid",
    "raw": "taxid",
    "s": "species",
    "sp": "species",
    "species": "species",
    "g": "genus",
    "genus": "genus",
    "f": "family",
    "fam": "family",
    "family": "family",
    "a": "any",
    "any": "any",
    "all": "any",
}

LEVEL_RANKS: Dict[str, Tuple[str, ...]] = {
    "taxid": (),
    "species": ("species",),
    "genus": ("species", "genus"),
    "family": ("species", "genus", "family"),
    "any": ("species", "genus", "family", "order", "class", "phylum"),
}


def normalize_reannotation_level(value: Optional[str]) -> str:
    key = str(value or "taxid").strip().lower()
    if key in {"", "none", "default"}:
        return "taxid"
    if key not in REANNOTATION_ALIASES:
        raise ValueError(
            f"Unknown reannotation-level={value!r}. "
            "Use taxid/t, species/s, genus/g, family/f, or any/a."
        )
    return REANNOTATION_ALIASES[key]


def ranks_for_level(level: Optional[str]) -> Tuple[str, ...]:
    return LEVEL_RANKS[normalize_reannotation_level(level)]


def _canon(taxid: str) -> str:
    key = canonical_taxid(taxid)
    if key == "0":
        return "0"
    try:
        return canonical_ncbi_taxid(key)
    except Exception:
        return key


def rank_bucket(taxid: str, rank: str) -> str:
    """Ancestor at ``rank``, or the taxid itself when that rank is missing (downward search)."""
    key = _canon(taxid)
    if key in {"", "0"}:
        return key
    exact = _resolve_taxid_by_rank_exact(key, rank)
    if exact and exact not in {"", "0"}:
        return _canon(exact)
    return key


def _catalog_rows(cfg: Optional[Dict[str, Any]] = None) -> List[Tuple[str, List[str]]]:
    rows = []
    for key, rec in genome_data_map(cfg).items():
        if is_test_taxid(key):
            continue
        rows.append((str(key), list(rec)))
    return rows


def catalog_hits_at_rank(
    query: str,
    rank: str,
    cfg: Optional[Dict[str, Any]] = None,
) -> List[Tuple[str, Path]]:
    """Indexed genomes whose taxon at ``rank`` matches ``query`` (up or down)."""
    want = rank_bucket(query, rank)
    if want in {"", "0"}:
        return []
    hits: List[Tuple[str, Path]] = []
    for key, rec in _catalog_rows(cfg):
        species = rec[0] if rec else ""
        bucket = rank_bucket(key, rank)
        if species and _canon(species) != _canon(key):
            alt = rank_bucket(species, rank)
            if alt == want or _canon(species) == want:
                bucket = alt if alt == want else bucket
        if bucket != want and _canon(key) != want and _canon(species) != want:
            continue
        path = resolve_indexed_path(key, cfg)
        if path is None or not path.is_file():
            continue
        hits.append((key, path))
    hits.sort(key=lambda item: (item[0] != _canon(query), item[0]))
    return hits


def catalog_hits_for_taxid(
    query: str,
    cfg: Optional[Dict[str, Any]] = None,
) -> List[Tuple[str, Path]]:
    key = _canon(query)
    path = resolve_indexed_path(key, cfg) or resolve_indexed_path(query, cfg)
    if path is not None and path.is_file():
        return [(key, path)]
    return []


def pick_catalog_genome(
    query: str,
    level: Optional[str] = "taxid",
    cfg: Optional[Dict[str, Any]] = None,
) -> Optional[Tuple[str, Path, str]]:
    """Return ``(matched_taxid, path, how)`` from the install catalog, or None."""
    query = str(query).split(".")[0]
    if canonical_taxid(query) == "0":
        return None
    exact = catalog_hits_for_taxid(query, cfg)
    if exact:
        return exact[0][0], exact[0][1], "taxid"
    for rank in ranks_for_level(level):
        hits = catalog_hits_at_rank(query, rank, cfg)
        if hits:
            return hits[0][0], hits[0][1], rank
    return None


def fetch_taxid_for_level(query: str, level: Optional[str], tried: Iterable[str]) -> Optional[str]:
    """Next NCBI taxid to download after catalog misses (exact, then rank ancestors)."""
    seen = {canonical_taxid(x) for x in tried}
    query = str(query).split(".")[0]
    cand = _canon(query)
    if cand not in seen and cand not in {"", "0"}:
        return cand
    for rank in ranks_for_level(level):
        bucket = rank_bucket(query, rank)
        if bucket not in seen and bucket not in {"", "0"}:
            return bucket
    return None


def resolve_genome_file(
    taxid: str,
    genome_dir: str,
    email: str,
    *,
    level: Optional[str] = "taxid",
    reference_only: bool = True,
    gzip_genomes: bool = True,
    cfg: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Local run dir, then catalog at ``level``, then NCBI (catalogued immediately)."""
    from samovar.genome_fetcher import fetch_genome
    from samovar.seqio import is_gzip_path
    from samovar.seqio import gzip_file, gunzip_file

    taxid = str(taxid).split(".")[0]
    local = find_genome_file(genome_dir, taxid)
    if local:
        return local

    hit = pick_catalog_genome(taxid, level, cfg)
    if hit is not None:
        _key, path, how = hit
        logger.info(
            "Reannotation %s: taxid %s -> catalog %s (%s) %s",
            normalize_reannotation_level(level),
            taxid,
            _key,
            how,
            path,
        )
        return str(path)

    tried: List[str] = []
    while True:
        fetch_id = fetch_taxid_for_level(taxid, level, tried)
        if not fetch_id:
            return None
        tried.append(fetch_id)
        try:
            genome_file = fetch_genome(
                fetch_id,
                genome_dir,
                email,
                reference_only=reference_only,
                gzip_genomes=gzip_genomes,
                index=True,
            )
        except Exception as exc:
            logger.warning("Failed to fetch genome for taxid %s: %s", fetch_id, exc)
            genome_file = None
        if genome_file:
            try:
                if gzip_genomes and not is_gzip_path(genome_file):
                    genome_file = str(gzip_file(genome_file))
                elif not gzip_genomes and is_gzip_path(genome_file):
                    genome_file = str(gunzip_file(genome_file, remove_source=True))
            except OSError as exc:
                logger.warning("Could not apply gzip_genomes=%s to %s: %s", gzip_genomes, genome_file, exc)
            return str(genome_file)
        if normalize_reannotation_level(level) == "taxid":
            return None
    return None
