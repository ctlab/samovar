"""CAMI / bioboxes taxonomic profiling format (RFC 0.10.0).

https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd

Direct counts come from abundance tables (``Annotation.to_abundance()``).
When NCBI taxdump is available, clade percentages and ``TAXPATH`` / ``TAXPATHSN``
are filled by walking ``nodes.dmp`` / ``names.dmp``. With ``--taxonomy gtdb``
the same walk uses GTDB (``@TaxonomyID:gtdb``, rank ``domain``).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from samovar.abundance import n_sample_columns
from samovar.kreport import (
    KReportTaxonomy,
    direct_counts_from_series,
    rollup_clade_counts,
    sample_from_n_column,
)
from samovar.taxonomy import (
    NCBI_CAMI_RANKS,
    cami_rank_for,
    cami_ranks_for,
    normalize_taxonomy,
)

PathLike = Union[str, Path]
ProfileMap = Dict[str, Dict[str, str]]

CAMI_FORMAT_VERSION = "0.10.0"
CAMI_RANKS = NCBI_CAMI_RANKS
CAMI_SUFFIXES = {".profile", ".cami", ".txt"}
_SAMPLE_ID_RE = re.compile(r"[A-Za-z0-9._]+")


def cami_rank_name(ncbi_rank: str, system: str = "ncbi") -> Optional[str]:
    """Map a dump rank string onto a CAMI ``RANKS`` token, or None."""
    return cami_rank_for(system, ncbi_rank)


def cami_sample_id(raw: Any) -> str:
    text = "".join(
        ch if ch.isalnum() or ch in "._" else "_" for ch in str(raw or "")
    )
    text = text.strip("._") or "sample"
    if not _SAMPLE_ID_RE.fullmatch(text):
        text = "sample"
    return text


def _cami_field(text: str) -> str:
    """Keep characters allowed in a CAMI output field."""
    return "".join(
        ch if ch.isalnum() or ch in " ,.;()_-" else " " for ch in str(text or "")
    ).strip()


def _fmt_pct(value: float) -> str:
    pct = max(0.0, min(100.0, float(value)))
    return f"{pct:.6f}"


def _lineage_slots(taxid: int, taxonomy: KReportTaxonomy) -> Dict[str, int]:
    slots: Dict[str, int] = {}
    system = getattr(taxonomy, "system", "ncbi")
    for node in taxonomy.ancestors_inclusive(taxid):
        mapped = cami_rank_name(taxonomy.rank_of(node), system)
        if mapped:
            slots[mapped] = node
    return slots


def _taxpath_pair(
    taxid: int,
    taxonomy: KReportTaxonomy,
    *,
    stop_rank: Optional[str],
    extra_leaf: bool,
) -> Tuple[str, str]:
    ranks = cami_ranks_for(getattr(taxonomy, "system", "ncbi"))
    slots = _lineage_slots(taxid, taxonomy)
    ids: List[str] = []
    names: List[str] = []
    if extra_leaf or not stop_rank:
        end = len(ranks)
    else:
        end = ranks.index(stop_rank) + 1
    for rank in ranks[:end]:
        node = slots.get(rank)
        ids.append(str(node) if node else "")
        names.append(_cami_field(taxonomy.name_of(node)) if node else "")
    while ids and ids[-1] == "":
        ids.pop()
        names.pop()
    if extra_leaf:
        ids.append(str(taxid))
        names.append(_cami_field(taxonomy.name_of(taxid)))
    if not ids:
        ids = [str(taxid)]
        names = [_cami_field(taxonomy.name_of(taxid))]
    return "|".join(ids), "|".join(names)


def _header(
    sample_id: str,
    *,
    with_names: bool = True,
    taxonomy: str = "ncbi",
) -> List[str]:
    system = normalize_taxonomy(taxonomy)
    ranks = cami_ranks_for(system)
    cols = (
        "TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"
        if with_names
        else "TAXID\tRANK\tTAXPATH\tPERCENTAGE"
    )
    return [
        f"@SampleID:{cami_sample_id(sample_id)}",
        f"@Version:{CAMI_FORMAT_VERSION}",
        f"@Ranks:{'|'.join(ranks)}",
        f"@TaxonomyID:{system}",
        f"@@{cols}",
    ]


def counts_to_cami_profile(
    direct: Mapping[int, int],
    taxonomy: KReportTaxonomy,
    sample_id: str = "1",
    *,
    total: Optional[float] = None,
) -> str:
    """One sample as a CAMI ``.profile`` (taxonomy-expanded TAXPATH)."""
    unclassified = int(direct.get(0, 0))
    classified = sum(n for tid, n in direct.items() if tid != 0 and n > 0)
    denom = (
        float(total)
        if total is not None and float(total) > 0
        else float(unclassified + classified)
    )
    if denom <= 0:
        denom = float(unclassified + classified)
    clade = rollup_clade_counts(direct, taxonomy)
    ranks = cami_ranks_for(getattr(taxonomy, "system", "ncbi"))
    rows: List[Tuple[int, str, int, float, str, str]] = []
    seen = set()
    for taxid, n_clade in clade.items():
        if taxid in {0, 1} or n_clade <= 0:
            continue
        mapped = cami_rank_name(taxonomy.rank_of(taxid), taxonomy.system)
        extra = mapped is None
        if mapped is None and int(direct.get(taxid, 0)) <= 0:
            continue
        if taxid in seen:
            continue
        seen.add(taxid)
        path, names = _taxpath_pair(
            taxid, taxonomy, stop_rank=mapped, extra_leaf=extra
        )
        rank = mapped or ""
        pct = 0.0 if denom <= 0 else 100.0 * float(n_clade) / denom
        order = ranks.index(mapped) if mapped else len(ranks)
        rows.append((order, rank, taxid, pct, path, names))
    rows.sort(key=lambda item: (item[0], -item[3], item[2]))
    lines = _header(sample_id, with_names=True, taxonomy=taxonomy.system)
    for _order, rank, taxid, pct, path, names in rows:
        lines.append(f"{taxid}\t{rank}\t{path}\t{names}\t{_fmt_pct(pct)}")
    return "\n".join(lines) + "\n"


def counts_to_cami_flat(
    counts: Mapping[Any, float],
    sample_id: str = "1",
    rank: str = "genus",
    *,
    total: Optional[float] = None,
    taxonomy: str = "ncbi",
) -> str:
    """RFC-shaped profile without taxdump (TAXPATH is the taxid itself)."""
    from samovar.kreport import _as_taxid
    from samovar.parse_annotators import taxid_ncbi_rank
    from samovar.viz_annotation import is_special_taxon, normalize_taxon_token

    system = normalize_taxonomy(taxonomy)
    ranks = cami_ranks_for(system)
    requested = rank if rank in ranks else "genus"
    classified = 0.0
    rows_raw: List[Tuple[str, float]] = []
    for taxid, n in counts.items():
        token = normalize_taxon_token(taxid)
        if is_special_taxon(token) or float(n) <= 0:
            continue
        if _as_taxid(token) == 0 and not str(token).strip():
            continue
        if str(token).strip().lower() in {"0", "unclassified"}:
            continue
        classified += float(n)
        rows_raw.append((token, float(n)))
    denom = float(total) if total is not None and float(total) > 0 else classified
    lines = _header(sample_id, with_names=False, taxonomy=system)
    for token, n in sorted(rows_raw, key=lambda kv: (-kv[1], kv[0])):
        actual = taxid_ncbi_rank(token) if str(token).replace(".", "", 1).isdigit() else None
        mapped = cami_rank_name(actual or "", system) if actual else None
        row_rank = mapped if mapped else ("" if actual else requested)
        pct = 0.0 if denom <= 0 else 100.0 * n / denom
        lines.append(f"{token}\t{row_rank}\t{token}\t{_fmt_pct(pct)}")
    return "\n".join(lines) + "\n"


def join_cami_sections(sections: Sequence[str]) -> str:
    """Concatenate sample sections (CAMI 0.10.0 multi-sample)."""
    parts = [str(s).strip("\n") for s in sections if str(s).strip()]
    return "\n\n".join(parts) + ("\n" if parts else "")


def try_taxonomy(
    taxdump: Optional[PathLike] = None,
    taxonomy: str = "ncbi",
) -> Optional[KReportTaxonomy]:
    try:
        return KReportTaxonomy.from_taxdump(taxdump, taxonomy=taxonomy)
    except FileNotFoundError:
        return None


def write_cami_profile(
    counts: Mapping[Any, float],
    dest: PathLike,
    sample_id: str = "1",
    rank: str = "genus",
    total: Optional[float] = None,
    taxdump: Optional[PathLike] = None,
    taxonomy: Optional[KReportTaxonomy] = None,
    taxonomy_id: str = "ncbi",
) -> Path:
    """Write a single-sample CAMI bioboxes ``.profile``."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    system = normalize_taxonomy(taxonomy_id)
    if taxonomy is not None:
        tax = taxonomy
    elif taxdump is not None:
        tax = try_taxonomy(taxdump, taxonomy=system)
    else:
        tax = None
    if tax is not None:
        direct = direct_counts_from_series(
            list(counts.keys()), counts.values(), tax
        )
        text = counts_to_cami_profile(direct, tax, sample_id, total=total)
    else:
        text = counts_to_cami_flat(
            counts, sample_id, rank, total=total, taxonomy=system
        )
    dest_path.write_text(text, encoding="utf-8")
    return dest_path


def abundance_tables_to_cami(
    tables: Mapping[str, Any],
    taxonomy: Optional[KReportTaxonomy],
) -> ProfileMap:
    """``{annotator: {sample: profile_text}}`` from ``taxid`` × ``N_*`` tables."""
    reports: ProfileMap = {}
    for annotator, table in tables.items():
        if table is None or getattr(table, "empty", True):
            continue
        cols = n_sample_columns(table)
        if not cols or "taxid" not in table.columns:
            continue
        sample_map: Dict[str, str] = {}
        for col in cols:
            total = float(table[col].sum())
            sid = sample_from_n_column(col)
            if taxonomy is not None:
                sample_map[sid] = counts_to_cami_profile(
                    direct_counts_from_series(
                        table["taxid"], table[col], taxonomy
                    ),
                    taxonomy,
                    sid,
                    total=total,
                )
            else:
                sample_map[sid] = counts_to_cami_flat(
                    dict(zip(table["taxid"], table[col])),
                    sid,
                    total=total,
                    taxonomy=getattr(taxonomy, "system", "ncbi"),
                )
        if sample_map:
            reports[str(annotator)] = sample_map
    return reports


def _safe_token(value: str) -> str:
    return (
        "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value))
        or "sample"
    )


def write_cami_profiles(dest: PathLike, reports: ProfileMap) -> Path:
    dest_path = Path(dest)
    n_annotators = len(reports)
    n_files = sum(len(samples) for samples in reports.values())
    if not reports or n_files == 0:
        raise ValueError("no CAMI profiles to write (empty abundance tables)")
    suffix = dest_path.suffix.lower()
    if suffix in CAMI_SUFFIXES and n_annotators == 1:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        samples = next(iter(reports.values()))
        dest_path.write_text(
            join_cami_sections(list(samples.values())), encoding="utf-8"
        )
        return dest_path
    if suffix in CAMI_SUFFIXES:
        dest_path = dest_path.with_suffix("")
    dest_path.mkdir(parents=True, exist_ok=True)
    for annotator, samples in reports.items():
        out = dest_path / f"{_safe_token(annotator)}.profile"
        out.write_text(join_cami_sections(list(samples.values())), encoding="utf-8")
    return dest_path


def dump_cami(
    tables: Mapping[str, Any],
    dest: PathLike,
    *,
    taxdump: Optional[PathLike] = None,
    taxonomy: str = "ncbi",
) -> Path:
    tax = try_taxonomy(taxdump, taxonomy=taxonomy)
    reports = abundance_tables_to_cami(tables, tax)
    return write_cami_profiles(dest, reports)
