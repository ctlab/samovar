"""Kraken 2 sample report (``--report`` / kreport) from per-taxid counts.

Kraken 2 writes this during classification; the same layout is produced by
KrakenTools ``make_kreport.py`` from a per-read table plus taxonomy. SamovaR
reuses ``Annotation.to_abundance()`` for the direct (leaf) counts, then rolls
them up the NCBI tree from taxdump (``nodes.dmp`` + ``names.dmp``).

Standard report columns (tab-separated):

1. Percentage of fragments in the clade rooted at this taxon
2. Number of fragments in that clade
3. Fragments assigned *directly* to this taxon
4. Rank code (``U``/``R``/``D``/``K``/``P``/``C``/``O``/``F``/``G``/``S``,
   or closest ancestor + distance, e.g. ``G2``)
5. NCBI taxID
6. Scientific name, indented by two spaces per tree level
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union

import pandas as pd

from samovar.abundance import n_sample_columns
from samovar.taxdump import find_dmp, load_scientific_names, names_dmp, nodes_dmp

PathLike = Union[str, Path]
TaxTree = Dict[int, Tuple[int, str]]
ReportMap = Dict[str, Dict[str, str]]

STANDARD_RANK_CODE = {
    "superkingdom": "D",
    "domain": "D",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
}
MPA_RANKS = ("D", "K", "P", "C", "O", "F", "G", "S")
KREPORT_SUFFIXES = {".kreport", ".report", ".kraken", ".kraken2", ".txt", ".tsv"}
MPA_SUFFIXES = {".mpa", ".txt", ".tsv"}
_UNCLASSIFIED = frozenset(
    {"0", "unclassified", "unclassified_root", "nan", "none", "na", ""}
)

_TAXDUMP_HINT = (
    "Kraken2-style reports need NCBI taxdump (nodes.dmp, names.dmp). "
    "Set genomes.taxdump, SAMOVAR_TAXDUMP, or pass --taxdump DIR."
)


class KReportTaxonomy:
    """Parent/rank tree plus scientific names used to expand taxIDs."""

    def __init__(
        self,
        tree: TaxTree,
        names: Optional[Mapping[int, str]] = None,
    ):
        self.tree = dict(tree)
        self.tree.setdefault(1, (1, "no rank"))
        self.names = {1: "root", **dict(names or {})}
        if 1 not in self.names:
            self.names[1] = "root"

    @classmethod
    def from_taxdump(cls, taxdump: Optional[PathLike] = None) -> "KReportTaxonomy":
        nodes_path, names_path = resolve_taxdump_files(taxdump)
        from samovar.parse_annotators import _load_nodes_cache

        tree = _load_nodes_cache(str(nodes_path))
        names: Dict[int, str] = {}
        if names_path is not None:
            names = load_scientific_names(names_path)
        return cls(tree, names)

    def parent_of(self, taxid: int) -> int:
        if taxid in self.tree:
            parent, _rank = self.tree[taxid]
            return parent
        return 1

    def rank_of(self, taxid: int) -> str:
        if taxid in self.tree:
            return self.tree[taxid][1]
        return "no rank"

    def name_of(self, taxid: int) -> str:
        if taxid == 0:
            return "unclassified"
        return self.names.get(taxid) or str(taxid)

    def rank_code(self, taxid: int) -> str:
        """Kraken 2 rank letter, including numbered intermediate ranks."""
        if taxid == 0:
            return "U"
        steps = 0
        current = taxid
        seen = set()
        while current not in seen:
            seen.add(current)
            if current == 1:
                return "R" if steps == 0 else f"R{steps}"
            rank = self.rank_of(current)
            code = STANDARD_RANK_CODE.get(rank)
            if code:
                return code if steps == 0 else f"{code}{steps}"
            parent = self.parent_of(current)
            if parent == current:
                return "-" if steps else "R"
            current = parent
            steps += 1
        return "-"

    def depth(self, taxid: int) -> int:
        if taxid in {0, 1}:
            return 0
        depth = 0
        current = taxid
        seen = set()
        while current not in seen and current != 1:
            seen.add(current)
            parent = self.parent_of(current)
            if parent == current:
                break
            current = parent
            depth += 1
        return depth

    def ancestors_inclusive(self, taxid: int) -> List[int]:
        """Walk taxid → root (inclusive), stopping at self-parent."""
        if taxid == 0:
            return [0]
        chain: List[int] = []
        current = taxid
        seen = set()
        while current not in seen:
            seen.add(current)
            chain.append(current)
            if current == 1:
                break
            parent = self.parent_of(current)
            if parent == current:
                if current != 1:
                    chain.append(1)
                break
            current = parent
        return chain

    def mpa_lineage(self, taxid: int) -> str:
        """``d__Bacteria|g__Escherichia|s__Escherichia_coli`` (standard ranks)."""
        parts: List[str] = []
        for node in reversed(self.ancestors_inclusive(taxid)):
            code = STANDARD_RANK_CODE.get(self.rank_of(node))
            if not code:
                continue
            name = self.name_of(node).replace(" ", "_")
            parts.append(f"{code.lower()}__{name}")
        return "|".join(parts)


def resolve_taxdump_files(
    taxdump: Optional[PathLike] = None,
) -> Tuple[Path, Optional[Path]]:
    """Locate ``nodes.dmp`` and (optional) ``names.dmp``."""
    if taxdump is not None:
        root = Path(taxdump).expanduser()
        if root.is_file() and root.name == "nodes.dmp":
            names = root.parent / "names.dmp"
            return root, names if names.is_file() else None
        if root.is_file() and root.name == "names.dmp":
            nodes = root.parent / "nodes.dmp"
            if not nodes.is_file():
                raise FileNotFoundError(_TAXDUMP_HINT)
            return nodes, root
        nodes = find_dmp("nodes.dmp", root)
        names = find_dmp("names.dmp", root)
        if nodes is None:
            raise FileNotFoundError(f"nodes.dmp not found under {root}. {_TAXDUMP_HINT}")
        return nodes, names
    nodes = nodes_dmp()
    names = names_dmp()
    if nodes is None:
        raise FileNotFoundError(_TAXDUMP_HINT)
    return nodes, names


def _as_taxid(value: Any) -> int:
    text = str(value).strip()
    if text.lower() in _UNCLASSIFIED:
        return 0
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return 0


def sample_from_n_column(column: str) -> str:
    text = str(column)
    return text[2:] if text.startswith("N_") else text


def direct_counts_from_series(
    taxids: Iterable[Any],
    values: Iterable[Any],
) -> Dict[int, int]:
    """Aggregate a taxa×count column into integer direct assignments."""
    counts: Dict[int, int] = defaultdict(int)
    for taxid, raw in zip(taxids, values):
        try:
            n = int(round(float(raw)))
        except (TypeError, ValueError):
            continue
        if n == 0:
            continue
        counts[_as_taxid(taxid)] += n
    return dict(counts)


def rollup_clade_counts(
    direct: Mapping[int, int],
    taxonomy: KReportTaxonomy,
) -> Dict[int, int]:
    """Clade size = direct + all descendants (walk each assigned taxid to root)."""
    clade: Dict[int, int] = defaultdict(int)
    for taxid, n in direct.items():
        if taxid == 0 or n <= 0:
            continue
        for node in taxonomy.ancestors_inclusive(taxid):
            clade[node] += n
    return dict(clade)


def _children_map(clade: Mapping[int, int], taxonomy: KReportTaxonomy) -> Dict[int, List[int]]:
    children: Dict[int, List[int]] = defaultdict(list)
    for taxid in clade:
        if taxid in {0, 1}:
            continue
        parent = taxonomy.parent_of(taxid)
        if parent == taxid:
            parent = 1
        children[parent].append(taxid)
    return children


def _dfs_preorder(
    clade: Mapping[int, int],
    children: Mapping[int, List[int]],
) -> List[int]:
    """Depth-first, children with larger clade counts first (Kraken 2 / make_kreport)."""
    if clade.get(1, 0) <= 0:
        return []
    order: List[int] = []
    stack: List[int] = [1]
    while stack:
        node = stack.pop(0)
        if clade.get(node, 0) <= 0:
            continue
        order.append(node)
        kids = [c for c in children.get(node, []) if clade.get(c, 0) > 0]
        kids.sort(key=lambda c: (clade.get(c, 0), c))
        for child in kids:
            stack.insert(0, child)
    return order


def _kreport_line(
    pct: float,
    clade: int,
    direct: int,
    rank: str,
    taxid: int,
    name: str,
    indent: int,
) -> str:
    return (
        f"{pct:6.2f}\t{int(clade)}\t{int(direct)}\t{rank}\t{int(taxid)}\t"
        f"{' ' * (indent * 2)}{name}"
    )


def counts_to_kreport(
    direct: Mapping[int, int],
    taxonomy: KReportTaxonomy,
    *,
    total: Optional[int] = None,
) -> str:
    """Render one sample as a Kraken 2 ``--report`` file (text)."""
    unclassified = int(direct.get(0, 0))
    classified = sum(n for tid, n in direct.items() if tid != 0 and n > 0)
    denom = int(total) if total is not None else unclassified + classified
    if denom <= 0:
        denom = unclassified + classified
    clade = rollup_clade_counts(direct, taxonomy)
    children = _children_map(clade, taxonomy)
    lines: List[str] = []
    if unclassified > 0:
        pct = 0.0 if denom <= 0 else 100.0 * unclassified / denom
        lines.append(
            _kreport_line(pct, unclassified, unclassified, "U", 0, "unclassified", 0)
        )
    for taxid in _dfs_preorder(clade, children):
        n_clade = int(clade.get(taxid, 0))
        n_direct = int(direct.get(taxid, 0))
        pct = 0.0 if denom <= 0 else 100.0 * n_clade / denom
        lines.append(
            _kreport_line(
                pct,
                n_clade,
                n_direct,
                taxonomy.rank_code(taxid),
                taxid,
                taxonomy.name_of(taxid),
                taxonomy.depth(taxid),
            )
        )
    return "\n".join(lines) + ("\n" if lines else "")


def counts_to_mpa(
    direct: Mapping[int, int],
    taxonomy: KReportTaxonomy,
) -> str:
    """Kraken 2 ``--use-mpa-style`` report (lineage + clade fragment count)."""
    clade = rollup_clade_counts(direct, taxonomy)
    children = _children_map(clade, taxonomy)
    lines: List[str] = []
    for taxid in _dfs_preorder(clade, children):
        code = taxonomy.rank_code(taxid)
        if not code or code[0] not in MPA_RANKS:
            continue
        if len(code) > 1 and code[1:].isdigit():
            continue
        lineage = taxonomy.mpa_lineage(taxid)
        if not lineage:
            continue
        lines.append(f"{lineage}\t{int(clade[taxid])}")
    return "\n".join(lines) + ("\n" if lines else "")


def abundance_tables_to_reports(
    tables: Mapping[str, pd.DataFrame],
    taxonomy: KReportTaxonomy,
    *,
    mpa: bool = False,
) -> ReportMap:
    """``{annotator: {sample: report_text}}`` from ``taxid`` × ``N_*`` tables."""
    reports: ReportMap = {}
    writer = counts_to_mpa if mpa else counts_to_kreport
    for annotator, table in tables.items():
        if table is None or getattr(table, "empty", True):
            continue
        cols = n_sample_columns(table)
        if not cols or "taxid" not in table.columns:
            continue
        sample_map: Dict[str, str] = {}
        for col in cols:
            direct = direct_counts_from_series(table["taxid"], table[col])
            sample_map[sample_from_n_column(col)] = writer(direct, taxonomy)
        if sample_map:
            reports[str(annotator)] = sample_map
    return reports


def _safe_token(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value)) or "sample"


def _report_suffix(mpa: bool) -> str:
    return ".mpa" if mpa else ".kreport"


def write_reports(
    dest: PathLike,
    reports: ReportMap,
    *,
    mpa: bool = False,
) -> Path:
    """Write kreport/MPA files. One sample+annotator and a file suffix → that file."""
    dest_path = Path(dest)
    n_files = sum(len(samples) for samples in reports.values())
    if not reports or n_files == 0:
        raise ValueError("no kraken2-style reports to write (empty abundance tables)")
    suffix = dest_path.suffix.lower()
    allowed = MPA_SUFFIXES if mpa else KREPORT_SUFFIXES
    if suffix in allowed and n_files == 1:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        annotator, samples = next(iter(reports.items()))
        dest_path.write_text(next(iter(samples.values())), encoding="utf-8")
        return dest_path
    if suffix in allowed:
        dest_path = dest_path.with_suffix("")
    dest_path.mkdir(parents=True, exist_ok=True)
    ext = _report_suffix(mpa)
    for annotator, samples in reports.items():
        folder = dest_path / _safe_token(annotator)
        folder.mkdir(parents=True, exist_ok=True)
        for sample, text in samples.items():
            (folder / f"{_safe_token(sample)}{ext}").write_text(text, encoding="utf-8")
    return dest_path


def dump_kreport(
    tables: Mapping[str, pd.DataFrame],
    dest: PathLike,
    *,
    taxdump: Optional[PathLike] = None,
    mpa: bool = False,
) -> Path:
    taxonomy = KReportTaxonomy.from_taxdump(taxdump)
    reports = abundance_tables_to_reports(tables, taxonomy, mpa=mpa)
    return write_reports(dest, reports, mpa=mpa)
