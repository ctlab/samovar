#!/usr/bin/env python3
"""Assign stable MAG ids by tetranucleotide homology.

Query FASTA stems are not ids. Profiles are clustered, then matched to a
catalog from the initial MAGs. A hit keeps that catalog id. A new cluster
receives the next ``M0001`` id, which is what combined annotation stores.
"""

from __future__ import annotations

import argparse
import fcntl
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np

from samovar.MAG.tetra import (
    DEFAULT_THRESHOLD,
    N_KMERS,
    load_fasta_profiles,
    profile_rms_z,
    profiles_match,
    read_fasta_sequence,
)

_FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def iter_mag_fastas(mag_dir: Path) -> List[Path]:
    src = Path(mag_dir)
    if src.is_file() and _is_fasta(src):
        return [src]
    if not src.is_dir():
        return []
    found = []
    for path in sorted(src.rglob("*")):
        if path.is_file() and _is_fasta(path) and not any(part.startswith(".") for part in path.relative_to(src).parts):
            found.append(path)
    return found


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suf) for suf in _FASTA_SUFFIXES)


def catalog_fasta_dir(catalog: Path) -> Path:
    return Path(catalog).parent / "fastas"


def _parse_counts(text: str) -> np.ndarray:
    values = np.zeros(N_KMERS, dtype=np.float64)
    parts = [p for p in str(text).split(",") if p != ""]
    for i, part in enumerate(parts[:N_KMERS]):
        values[i] = float(part)
    return values


def _format_counts(counts: np.ndarray) -> str:
    return ",".join(str(int(v)) for v in np.asarray(counts).tolist())


def load_catalog(path: Path) -> List[dict]:
    src = Path(path)
    if not src.is_file():
        return []
    rows = []
    lines = src.read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 3 or not parts[0]:
            continue
        rows.append(
            {
                "mag_id": parts[0],
                "n": int(float(parts[1] or 0)),
                "counts": _parse_counts(parts[2]),
                "source": parts[3] if len(parts) > 3 else "",
            }
        )
    return rows


def _write_catalog(path: Path, rows: Sequence[dict]) -> None:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    lines = ["mag_id\tn\tcounts\tsource"]
    for row in rows:
        lines.append(
            f"{row['mag_id']}\t{int(row['n'])}\t{_format_counts(row['counts'])}\t{row.get('source', '')}"
        )
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _next_ids(existing: Sequence[str], n_new: int) -> List[str]:
    used = 0
    for name in existing:
        text = str(name)
        if len(text) > 1 and text[0] == "M" and text[1:].isdigit():
            used = max(used, int(text[1:]))
    return [f"M{used + i + 1:04d}" for i in range(n_new)]


class _UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, i: int) -> int:
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def _cluster(records: Sequence[dict], threshold: float) -> List[List[int]]:
    uf = _UnionFind(len(records))
    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            rms = profile_rms_z(
                records[i]["counts"], records[i]["n"], records[j]["counts"], records[j]["n"]
            )
            if profiles_match(rms, threshold):
                uf.union(i, j)
    groups: Dict[int, List[int]] = {}
    for i in range(len(records)):
        groups.setdefault(uf.find(i), []).append(i)
    return list(groups.values())


def _representative(records: Sequence[dict], members: Sequence[int]) -> int:
    return max(members, key=lambda i: (records[i]["length"], records[i]["name"]))


def assign_ids(
    query_fastas: Sequence[Path],
    catalog_rows: Optional[Sequence[dict]] = None,
    *,
    threshold: float = DEFAULT_THRESHOLD,
) -> List[dict]:
    """Map each query FASTA to a stable mag id.

    Homologous queries share an id. A catalog hit reuses that id; otherwise
    the cluster gets a new ``M####`` id. Returned rows follow the query order.
    """
    records = load_fasta_profiles(query_fastas)
    if not records:
        return []
    if catalog_rows is None:
        catalog = []
    else:
        catalog = catalog_rows
    clusters = _cluster(records, threshold)
    ranked = sorted(
        clusters,
        key=lambda members: (
            int(np.argmax(records[_representative(records, members)]["counts"])),
            -records[_representative(records, members)]["length"],
            records[_representative(records, members)]["name"],
        ),
    )
    fresh = _next_ids([row["mag_id"] for row in catalog], len(ranked))
    by_index = {}
    new_rows = []
    fresh_i = 0
    for members in ranked:
        rep = _representative(records, members)
        best_id = ""
        best_rms = None
        for row in catalog:
            rms = profile_rms_z(records[rep]["counts"], records[rep]["n"], row["counts"], row["n"])
            if best_rms is None or rms < best_rms:
                best_rms = rms
                best_id = row["mag_id"]
        matched = best_id if best_rms is not None and profiles_match(best_rms, threshold) else ""
        if matched:
            mag_id = matched
            reference = matched
            rms_value = float(best_rms)
        else:
            mag_id = fresh[fresh_i]
            fresh_i += 1
            reference = ""
            rms_value = float(best_rms) if best_rms is not None else 0.0
            new_rows.append(
                {
                    "mag_id": mag_id,
                    "n": int(records[rep]["n"]),
                    "counts": records[rep]["counts"],
                    "source": records[rep]["name"],
                }
            )
        for i in members:
            pair = profile_rms_z(
                records[i]["counts"], records[i]["n"], records[rep]["counts"], records[rep]["n"]
            )
            by_index[i] = {
                "query": records[i]["name"],
                "path": records[i]["path"],
                "mag_id": mag_id,
                "reference": reference,
                "rms_z": pair if not reference else (rms_value if i == rep else pair),
                "match": 1 if reference or len(members) > 1 else 0,
            }
    catalog.extend(new_rows)
    # Stash updated catalog on the function result via attribute for callers
    # that pass the same list object: ``catalog_rows`` is mutated above.
    return [by_index[i] for i in range(len(records))]


def write_assignment_table(rows: Sequence[dict], dest: Path) -> None:
    lines = ["query\tmag_id\treference\trms_z\tmatch"]
    for row in rows:
        lines.append(
            f"{row['query']}\t{row['mag_id']}\t{row.get('reference', '')}\t"
            f"{float(row['rms_z']):.6g}\t{int(row['match'])}"
        )
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def stabilize_mags(
    query_dir: Path,
    catalog_path: Path,
    fasta_out: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    reference_dir: Optional[Path] = None,
    preserve_contigs: bool = False,
) -> List[dict]:
    """Rewrite ``query_dir`` into ``fasta_out`` using catalog ids.

    The catalog file is locked and extended with new clusters. Homologous
    FASTAs are concatenated under one ``<mag_id>.fa``.
    """
    catalog_path = Path(catalog_path)
    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = catalog_path.with_suffix(catalog_path.suffix + ".lock")
    with lock_path.open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        catalog = load_catalog(catalog_path)
        if reference_dir is not None and not catalog:
            for rec in load_fasta_profiles(iter_mag_fastas(Path(reference_dir))):
                catalog.append(
                    {
                        "mag_id": rec["name"],
                        "n": int(rec["n"]),
                        "counts": rec["counts"],
                        "source": rec["name"],
                    }
                )
        rows = assign_ids(iter_mag_fastas(query_dir), catalog, threshold=threshold)
        _write_catalog(catalog_path, catalog)
        fasta_out = Path(fasta_out)
        if fasta_out.exists():
            for old in fasta_out.glob("*"):
                if old.is_file():
                    old.unlink()
        fasta_out.mkdir(parents=True, exist_ok=True)
        grouped: Dict[str, List[Path]] = {}
        for row in rows:
            grouped.setdefault(row["mag_id"], []).append(Path(row["path"]))
        stored = catalog_fasta_dir(catalog_path)
        stored.mkdir(parents=True, exist_ok=True)
        for mag_id, paths in grouped.items():
            if preserve_contigs:
                from samovar.seqio import iter_seqio_fasta
                records = []
                for source_i, path in enumerate(paths):
                    for contig_i, record in enumerate(iter_seqio_fasta(path)):
                        sequence = str(record.seq)
                        wrapped = "\n".join(sequence[i:i+80] for i in range(0,len(sequence),80))
                        records.append(f">{mag_id}_c{source_i}_{contig_i}\n{wrapped}\n")
                body = "".join(records)
            else:
                sequence = "".join(read_fasta_sequence(path) for path in paths)
                wrapped = "\n".join(sequence[i : i + 80] for i in range(0, len(sequence), 80))
                body = f">{mag_id}\n{wrapped}\n"
            (fasta_out / f"{mag_id}.fa").write_text(body, encoding="utf-8")
            stored_fa = stored / f"{mag_id}.fa"
            if not stored_fa.is_file():
                stored_fa.write_text(body, encoding="utf-8")
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    return rows


def _reference_catalog(reference: str) -> List[dict]:
    if not reference:
        return []
    src = Path(reference)
    if src.is_file():
        return load_catalog(src)
    if src.is_dir():
        rows = []
        for rec in load_fasta_profiles(iter_mag_fastas(src)):
            rows.append(
                {
                    "mag_id": rec["name"],
                    "n": int(rec["n"]),
                    "counts": rec["counts"],
                    "source": rec["name"],
                }
            )
        return rows
    return []


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m samovar.MAG.homology")
    parser.add_argument("-c", required=True, help="Query MAG directory")
    parser.add_argument("-r", default="", help="Reference catalog TSV or MAG directory")
    parser.add_argument("-o", required=True, help="Assignment TSV")
    parser.add_argument("--fasta-out", default="", help="Directory of FASTAs renamed to stable ids")
    parser.add_argument("--catalog", default="", help="Catalog TSV to create or extend")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(list(argv) if argv is not None else None)
    query = Path(args.c)
    if args.fasta_out or args.catalog:
        catalog = Path(args.catalog or args.r or (Path(args.o).with_suffix(".catalog.tsv")))
        rows = stabilize_mags(
            query,
            catalog,
            Path(args.fasta_out or (query.parent / (query.name + "_stable"))),
            threshold=float(args.threshold),
            reference_dir=Path(args.r) if args.r and Path(args.r).is_dir() else None,
        )
    else:
        catalog_rows = _reference_catalog(args.r)
        rows = assign_ids(iter_mag_fastas(query), catalog_rows, threshold=float(args.threshold))
    write_assignment_table(rows, Path(args.o))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# A completed map-initial reference is immutable during read assignment.
def mapping_manifest(catalog: Path) -> Path:
    return Path(catalog).parent / 'mapping.manifest.json'


def _file_sha256(path: Path) -> str:
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def extract_gtdb_markers(work: Path, dest: Path) -> Path:
    """Join identify marker summaries/top-hits to GTDB-Tk's own Prodigal genes.

    Coordinates are zero-based half-open. Only HMMs in the GTDB marker
    summaries qualify; a generic Pfam/TIGRFAM hit is not marker evidence.
    """
    import csv
    from samovar.seqio import iter_seqio_fasta

    allowed = {}
    for summary in sorted(Path(work).rglob('*markers_summary.tsv')):
        with summary.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                mag = row['name']
                allowed.setdefault(mag, set()).update(
                    marker for key in ('list_unique_genes', 'list_multiple_genes',
                                       'list_multiple_unique_genes')
                    for marker in row.get(key, '').split(',') if marker)
    rows = set()
    for fasta in sorted(Path(work).rglob('*_protein.faa')):
        mag = fasta.name[:-len('_protein.faa')]
        genes = {}
        for gene in iter_seqio_fasta(fasta):
            fields = gene.description.split(' # ')
            if len(fields) < 4:
                raise ValueError(f'Missing GTDB gene coordinates: {gene.id}')
            genes[gene.id] = (gene.id.rsplit('_', 1)[0], int(fields[1])-1, int(fields[2]))
        for hits in sorted(fasta.parent.glob('*_tophit.tsv')):
            with hits.open() as handle:
                next(handle, None)
                for line in handle:
                    gene, matches = line.rstrip('\n').split('\t')
                    for match in matches.split(';'):
                        marker = match.split(',')[0]
                        if marker not in allowed.get(mag, set()):
                            continue
                        if gene not in genes:
                            raise ValueError(f'GTDB marker gene missing from protein FASTA: {gene}')
                        contig, start, end = genes[gene]
                        rows.add((mag,contig,start,end,gene,marker))
    found = {row[0] for row in rows}
    missing = [mag for mag, markers in allowed.items() if markers and mag not in found]
    if missing:
        raise ValueError(f'GTDB marker coordinates are missing for: {", ".join(sorted(missing))}')
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open('w') as handle:
        handle.write('mag_id\tcontig\tstart\tend\tgene_id\tmarker_id\n')
        for row in sorted(rows):
            handle.write('\t'.join(map(str,row)) + '\n')
    return dest


def finalize_mapping_catalog(catalog: Path, *, db: str = '', threads: int = 1,
                             taxonomy_tool: str = 'gtdbtk', samples=()) -> Path:
    import json
    from samovar.assembly_profiling import run_mag_taxonomy, _fasta_has_records, contig_to_mag_map

    catalog = Path(catalog)
    root = catalog.parent
    root.mkdir(parents=True, exist_ok=True)
    manifest = mapping_manifest(catalog)
    manifest.unlink(missing_ok=True)
    fastas = catalog_fasta_dir(catalog)
    fastas.mkdir(exist_ok=True)
    tax = root / 'taxonomy.tsv'
    markers = root / 'markers.tsv'
    if any(_fasta_has_records(p) for p in iter_mag_fastas(fastas)):
        run_mag_taxonomy(str(fastas), str(tax), db=db, threads=threads, name=taxonomy_tool,
                         require_markers=True)
        if taxonomy_tool.lower() in {'gtdbtk', 'gtdb', 'mag_taxonomy', ''}:
            extract_gtdb_markers(tax.with_suffix('.gtdbtk'), markers)
        else:
            # Custom taxonomy may supply this sidecar itself; never label
            # arbitrary gene-caller output as GTDB marker evidence.
            if not markers.is_file():
                markers.write_text('mag_id\tcontig\tstart\tend\tgene_id\tmarker_id\n')
    else:
        tax.write_text('mag_id\ttaxid\tlineage\n')
        markers.write_text('mag_id\tcontig\tstart\tend\tgene_id\tmarker_id\n')
    contigs = root / 'contigs.tsv'
    with contigs.open('w') as handle:
        handle.write('contig\tmag_id\n')
        for contig, mag in sorted(contig_to_mag_map(fastas).items()):
            handle.write(f'{contig}\t{mag}\n')
    files = [*iter_mag_fastas(fastas), tax, markers, contigs]
    if catalog.is_file():
        files.append(catalog)
    payload = {'schema_version':1, 'samples':sorted(samples), 'files':{
        str(path.relative_to(root)): _file_sha256(path) for path in files}}
    pending = manifest.with_suffix('.pending')
    pending.write_text(json.dumps(payload, indent=2) + '\n')
    pending.replace(manifest)
    return manifest


def validate_mapping_catalog(catalog: Path, *, samples=None) -> dict:
    import json

    catalog = Path(catalog)
    manifest = mapping_manifest(catalog)
    if not manifest.is_file():
        raise ValueError('Initial MAG catalog is incomplete or legacy; rerun initial annotation and downstream stages')
    payload = json.loads(manifest.read_text())
    if payload.get('schema_version') != 1:
        raise ValueError('Unsupported initial MAG catalog format; rerun initial annotation')
    if samples is not None and payload.get('samples') != sorted(samples):
        raise ValueError('Initial MAG cohort changed; prepare a fresh run and rebuild initial annotation')
    root = catalog.parent.resolve()
    expected = payload.get('files', {})
    for name in ('taxonomy.tsv', 'markers.tsv', 'contigs.tsv'):
        if name not in expected:
            raise ValueError(f'Incomplete initial MAG catalog: {name}')
    for name, checksum in expected.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file() or _file_sha256(path) != checksum:
            raise ValueError(f'Initial MAG reference changed or is corrupt: {name}')
    actual_fastas = {str(p.relative_to(root)) for p in iter_mag_fastas(root / 'fastas')}
    if actual_fastas != {name for name in expected if name.startswith('fastas/')}:
        raise ValueError('Initial MAG FASTA set changed; rerun initial annotation')
    return payload
