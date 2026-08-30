"""Optional CAMISIM backend for ``samovar generate``.

https://github.com/CAMI-challenge/CAMISIM

Modes
-----
``table``
    CAMISIM-style community design only, then InSilicoSeq (ISS) for reads.
``illumina`` / ``ont`` / ``wgsim``
    CAMISIM read simulation for that technology (no ISS) when the type exists.
``hybrid``
    Same community, one CAMISIM run per available technology; FASTQ headers
    get ``read_type:<tech>`` so annotation tables keep a read-type column.

CAMISIM 2 (Nextflow) is preferred; CAMISIM 1.3 ``metagenomesimulation.py`` is
the fallback. Community design always has a Python implementation so ``table``
mode works without Nextflow.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import random
import re
import shutil
import subprocess
import tarfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import yaml

from samovar.table_regenerators import extra_flags_argv
from samovar.regenerate import cap_generate_genome_rows, parse_max_genomes
from samovar.main_config import iter_tools, tool_path as tool_entry_path
from samovar.paths import (
    absolute_path,
    conda_prefix_for_executable,
    discover_art,
    discover_nanosim,
    load_config,
    python_path,
    repo_root,
    runtime_path_prefix,
    tool_env_prefix,
    user_cache_dir,
)
from samovar.seqio import (
    concat_fastq_files,
    is_fastq_name,
    iter_fastq_records,
    list_fasta_files,
    open_text,
    taxid_from_fasta_name,
    write_text_lines,
)

logger = logging.getLogger(__name__)

CAMISIM_REPO = "https://github.com/CAMI-challenge/CAMISIM"

# CAMISIM 2 ``params.type`` values and SamovaR mode aliases.
CAMISIM_TYPES = {
    "art": {
        "mode": "illumina",
        "read_type": "illumina",
        "paired": True,
        "label": "Illumina (ART)",
    },
    "nanosim3": {
        "mode": "ont",
        "read_type": "ont",
        "paired": False,
        "label": "Oxford Nanopore (NanoSim3)",
    },
    "wgsim": {
        "mode": "wgsim",
        "read_type": "wgsim",
        "paired": True,
        "label": "wgsim",
    },
}

MODE_TO_TYPE = {
    "illumina": "art",
    "art": "art",
    "hiseq": "art",
    "miseq": "art",
    "novaseq": "art",
    "ont": "nanosim3",
    "nanopore": "nanosim3",
    "nanosim": "nanosim3",
    "nanosim3": "nanosim3",
    "wgsim": "wgsim",
    "etc": "wgsim",
}

READ_TYPE_RE = re.compile(r"read_type:([A-Za-z0-9_+\-]+)", re.IGNORECASE)
TAXID_RE = re.compile(r"taxid:([0-9]+)", re.IGNORECASE)


def normalize_camisim_mode(raw: Optional[str]) -> str:
    text = str(raw or "table").strip().lower()
    if text in {"", "table", "community", "abundance", "design"}:
        return "table"
    if text in {"hybrid", "mix", "multi"}:
        return "hybrid"
    if text in MODE_TO_TYPE:
        return CAMISIM_TYPES[MODE_TO_TYPE[text]]["mode"]
    if text in {"reads", "read"}:
        return "illumina"
    raise ValueError(
        f"Unknown CAMISIM mode {raw!r}. Use table, illumina, ont, wgsim, or hybrid."
    )


def extract_read_type(seq_id: Optional[str]) -> str:
    """``read_type:<token>`` from a FASTQ / annotation sequence id."""
    text = "" if seq_id is None else str(seq_id)
    match = READ_TYPE_RE.search(text)
    return match.group(1).lower() if match else ""


def tag_read_id(header_or_id: str, read_type: str, taxid: str = "") -> str:
    """Append ``taxid:`` / ``read_type:`` tokens to a FASTQ id if missing."""
    raw = str(header_or_id or "").rstrip("\r\n")
    prefix = "@" if raw.startswith("@") else ""
    body = raw[1:] if prefix else raw
    if not body:
        body = "read"
    parts = body.split(None, 1)
    rid = parts[0]
    rest = (" " + parts[1]) if len(parts) > 1 else ""
    if taxid and not TAXID_RE.search(rid):
        rid = f"{rid}|taxid:{taxid}"
    token = str(read_type or "").strip().lower()
    if token and not READ_TYPE_RE.search(rid):
        rid = f"{rid}|read_type:{token}"
    return f"{prefix}{rid}{rest}"


def discover_camisim() -> Optional[str]:
    """Locate a CAMISIM checkout (Nextflow ``main.nf`` or 1.3 Python)."""
    cfg = load_config()
    candidates: List[str] = []

    def add(raw: Optional[str]) -> None:
        text = str(raw or "").strip()
        if text:
            candidates.append(text)

    add(os.environ.get("SAMOVAR_CAMISIM") or os.environ.get("SAMOVAR_CAMISIM_ROOT"))
    add(str(cfg.get("camisim_path") or cfg.get("camisim_root") or ""))
    add(tool_entry_path(iter_tools(cfg).get("camisim"), "camisim"))
    add(str(user_cache_dir() / "CAMISIM"))
    add(str(Path.home() / ".config" / "samovar" / "CAMISIM"))
    try:
        add(str(repo_root() / "third_party" / "CAMISIM"))
    except OSError:
        pass
    seen = set()
    for raw in candidates:
        if raw in seen:
            continue
        seen.add(raw)
        path = Path(raw).expanduser()
        try:
            if not path.exists():
                continue
        except OSError:
            continue
        if path.is_file():
            path = path.parent
        if _camisim_kind(path):
            try:
                return str(path.resolve())
            except OSError:
                return str(path)
    return None


def discover_nextflow() -> Optional[str]:
    cfg = load_config()
    nxt = tool_entry_path(iter_tools(cfg).get("nextflow"), "nextflow")
    for raw in (
        os.environ.get("SAMOVAR_NEXTFLOW"),
        cfg.get("nextflow_path"),
        nxt,
        shutil.which("nextflow"),
    ):
        text = str(raw or "").strip()
        if not text:
            continue
        path = Path(text).expanduser()
        try:
            if path.is_file():
                return str(path.resolve())
        except OSError:
            continue
        found = shutil.which(text)
        if found:
            return found
    return None


def _camisim_kind(root: Path) -> Optional[str]:
    try:
        if (root / "main.nf").is_file():
            return "nextflow"
        if (root / "metagenomesimulation.py").is_file():
            return "legacy"
        if (root / "scripts" / "metagenomesimulation.py").is_file():
            return "legacy"
    except OSError:
        return None
    return None


def available_camisim_types(root: Optional[str] = None) -> List[str]:
    """Types advertised by a CAMISIM checkout (defaults to the 2.0 set)."""
    path = Path(root) if root else (Path(discover_camisim()) if discover_camisim() else None)
    if path is None:
        return list(CAMISIM_TYPES)
    kind = _camisim_kind(path)
    if kind == "nextflow":
        found = []
        cfg_dir = path / "pipelines" / "metagenomic" / "config"
        for name in ("art", "nanosim3", "wgsim"):
            if (cfg_dir / f"{name}.config").is_file() or name in CAMISIM_TYPES:
                found.append(name)
        return found or list(CAMISIM_TYPES)
    # 1.3 ini documents art / nanosim / wgsim
    return list(CAMISIM_TYPES)


def types_for_mode(mode: str, root: Optional[str] = None) -> List[str]:
    mode = normalize_camisim_mode(mode)
    available = available_camisim_types(root)
    if mode == "table":
        return []
    if mode == "hybrid":
        # Prefer Illumina + ONT when both exist; otherwise every available type.
        prefer = [t for t in ("art", "nanosim3") if t in available]
        chosen = prefer if len(prefer) >= 2 else list(available)
        if len(chosen) < 2:
            raise RuntimeError(
                "hybrid CAMISIM mode needs at least two read simulators "
                f"(have {available or 'none'}). Install CAMISIM or pick illumina/ont/wgsim."
            )
        return chosen
    cam_type = MODE_TO_TYPE.get(mode)
    if cam_type not in available:
        raise RuntimeError(
            f"CAMISIM mode {mode!r} ({cam_type}) is not present in this CAMISIM "
            f"install (available: {', '.join(available) or 'none'})."
        )
    return [cam_type]


def collect_generate_genomes(
    genome_dir: str,
    host_genome: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Map input FASTA files to CAMISIM genome_ID / NCBI_ID rows."""
    rows: List[Dict[str, str]] = []
    seen = set()

    def add(path: Path, *, host: bool = False) -> None:
        try:
            if not path.is_file():
                return
            key = str(path.resolve())
        except OSError:
            key = str(path)
        if key in seen:
            return
        seen.add(key)
        taxid = taxid_from_fasta_name(path) or ("9606" if host else path.stem.split(".")[0])
        gid = str(taxid)
        suffix = 2
        existing = {r["genome_ID"] for r in rows}
        while gid in existing:
            gid = f"{taxid}_{suffix}"
            suffix += 1
        rows.append(
            {
                "genome_ID": gid,
                "path": str(path),
                "NCBI_ID": str(taxid) if str(taxid).isdigit() else "1",
                "OTU": gid,
                "novelty_category": "known_strain",
                "host": "1" if host else "0",
            }
        )

    for fasta in list_fasta_files(genome_dir, nucleotide=True, protein=False):
        add(fasta, host=False)
    if host_genome:
        add(Path(host_genome).expanduser(), host=True)
    if not rows:
        raise FileNotFoundError(f"No nucleotide FASTA files in {genome_dir}")
    return rows


def select_camisim_genome_rows(cfg: Dict[str, Any]) -> List[Dict[str, str]]:
    """Use locked ``cfg['genomes']`` when present; else cap by ``max_genomes``."""
    rows = collect_generate_genomes(cfg["genome_dir"], cfg.get("host_genome"))
    locked = [str(x) for x in (cfg.get("genomes") or []) if str(x).strip()]
    if locked:
        from samovar.iss_config import resolve_locked_genome_paths

        resolved = resolve_locked_genome_paths(cfg["genome_dir"], locked)
        want = set()
        for path in resolved:
            try:
                want.add(str(path.resolve()))
            except OSError:
                want.add(str(path))
        if not want:
            want = {str(Path(x).expanduser()) for x in locked}
        ordered: List[Dict[str, str]] = []
        seen = set()
        for row in rows:
            try:
                key = str(Path(row["path"]).resolve())
            except OSError:
                key = str(row["path"])
            if key in want and key not in seen:
                ordered.append(row)
                seen.add(key)
        if ordered:
            hosts = [r for r in rows if str(r.get("host") or "") == "1"]
            meta = [r for r in ordered if str(r.get("host") or "") != "1"]
            extra_hosts = []
            for host in hosts:
                try:
                    key = str(Path(host["path"]).resolve())
                except OSError:
                    key = str(host["path"])
                if key not in seen:
                    extra_hosts.append(host)
            return meta + extra_hosts
    return cap_generate_genome_rows(rows, cfg.get("max_genomes"))


def camisim_size_gbp(
    total_reads: int,
    read_length: int = 150,
    min_gbp: float = 1e-6,
    read_ends: int = 1,
) -> float:
    """Convert requested records/pairs to CAMISIM's per-sample Gbp size."""
    gbp = (
        max(int(total_reads), 1)
        * max(int(read_length), 1)
        * max(int(read_ends), 1)
        / 1e9
    )
    return max(round(gbp, 8), float(min_gbp))


def camisim_sizes_by_type(
    total_reads: int,
    types: Sequence[str],
    *,
    art_read_length: int = 150,
    nanosim_read_length: int = 4508,
    wgsim_read_length: int = 150,
) -> Dict[str, float]:
    """Gbp per simulator, splitting a hybrid run's requested records.

    ART's size is sequenced bases (two ends per pair), while CAMISIM's wgsim
    module divides size by one read length to obtain the number of pairs.
    NanoSim emits single-end reads whose mean length comes from the profile.
    """
    chosen = [str(t) for t in types]
    if not chosen:
        return {}
    total = max(int(total_reads), 1)
    base, remainder = divmod(total, len(chosen))
    counts = {
        cam_type: max(base + (1 if i < remainder else 0), 1)
        for i, cam_type in enumerate(chosen)
    }
    sizes: Dict[str, float] = {}
    for cam_type, count in counts.items():
        if cam_type == "art":
            sizes[cam_type] = camisim_size_gbp(
                count, art_read_length, min_gbp=1e-8, read_ends=2
            )
        elif cam_type == "nanosim3":
            sizes[cam_type] = camisim_size_gbp(
                count, nanosim_read_length, min_gbp=1e-8
            )
        else:
            sizes[cam_type] = camisim_size_gbp(
                count, wgsim_read_length, min_gbp=1e-8
            )
    return sizes


def camisim_taxdump(root: Optional[str]) -> str:
    """NCBI taxdump tarball: install ``genomes.taxdump``, else CAMISIM bundled dump."""
    try:
        from samovar.taxdump import taxdump_tarball

        archive = taxdump_tarball()
        if archive is not None:
            return str(archive.resolve())
    except Exception:
        pass
    if not root:
        return ""
    tools = Path(root) / "tools"
    preferred = tools / "ncbi-taxonomy_20170222.tar.gz"
    try:
        if preferred.is_file():
            return str(preferred.resolve())
        matches = sorted(tools.glob("ncbi-taxonomy*.tar.gz"))
        if matches:
            return str(matches[0].resolve())
    except OSError:
        return ""
    return ""


# Modern NCBI re-ids that are absent from CAMISIM's 2017 taxdump.
TAXID_ALIASES = {
    "2886930": "10847",  # Escherichia phage phiX174
}

_TAXDUMP_IDS_CACHE: Dict[str, Tuple[float, set]] = {}


def load_taxdump_taxids(taxdump: Optional[str]) -> set:
    """Taxids present in ``nodes.dmp`` (tarball or extracted taxdump directory)."""
    if not taxdump:
        return set()
    path = Path(taxdump)
    try:
        if path.is_dir():
            nodes = path / "nodes.dmp"
            if not nodes.is_file():
                nodes = path / "taxonomy" / "nodes.dmp"
            if not nodes.is_file():
                return set()
            mtime = nodes.stat().st_mtime
            key = str(nodes)
            cached = _TAXDUMP_IDS_CACHE.get(key)
            if cached and cached[0] == mtime:
                return cached[1]
            ids: set = set()
            with nodes.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    taxid = line.split("\t", 1)[0].strip()
                    if taxid.isdigit():
                        ids.add(taxid)
            _TAXDUMP_IDS_CACHE[key] = (mtime, ids)
            return ids
        if not path.is_file():
            return set()
        mtime = path.stat().st_mtime
    except OSError:
        return set()
    key = str(path)
    cached = _TAXDUMP_IDS_CACHE.get(key)
    if cached and cached[0] == mtime:
        return cached[1]
    ids = set()
    try:
        with tarfile.open(path) as tar:
            member = None
            for info in tar.getmembers():
                name = info.name.replace("\\", "/").rstrip("/")
                if name == "nodes.dmp" or name.endswith("/nodes.dmp"):
                    member = info
                    break
            if member is None:
                return set()
            handle = tar.extractfile(member)
            if handle is None:
                return set()
            for raw in handle:
                line = raw.decode("utf-8", errors="replace")
                taxid = line.split("\t", 1)[0].strip()
                if taxid.isdigit():
                    ids.add(taxid)
    except (OSError, tarfile.TarError):
        return set()
    _TAXDUMP_IDS_CACHE[key] = (mtime, ids)
    return ids


def map_ncbi_id_for_taxdump(taxid: str, dump_ids: Optional[set] = None) -> str:
    """Return a taxid CAMISIM's taxdump can resolve (or the original if unknown)."""
    tid = str(taxid or "").strip()
    if not tid:
        return "1"
    if dump_ids is None or not dump_ids:
        return tid
    if tid in dump_ids:
        return tid
    alias = TAXID_ALIASES.get(tid)
    if alias and alias in dump_ids:
        return alias
    logger.warning(
        "taxid %s is not in the CAMISIM NCBI taxdump; using taxid 1 in CAMISIM metadata",
        tid,
    )
    return "1"


def _deep_overlay(base: Dict[str, Any], overlay: Dict[str, Any]) -> Dict[str, Any]:
    """Merge ``overlay`` onto ``base``; nested dicts are updated, not replaced."""
    merged = dict(base)
    nested = ("distribution", "art", "nanosim3", "wgsim", "size_gbp_by_type")
    for key, value in overlay.items():
        if key in nested and isinstance(value, dict) and isinstance(merged.get(key), dict):
            nested_val = dict(merged[key])
            nested_val.update(value)
            merged[key] = nested_val
        elif key not in merged:
            merged[key] = value
    for key in nested:
        extra = overlay.get(key)
        if isinstance(extra, dict) and isinstance(merged.get(key), dict):
            nested_val = dict(merged[key])
            nested_val.update(extra)
            merged[key] = nested_val
    return merged


def design_communities(
    genome_ids: Sequence[str],
    n_samples: int,
    seed: int = 42,
    mode: str = "differential",
    log_mu: float = 1.0,
    log_sigma: float = 2.0,
    gauss_mu: float = 1.0,
    gauss_sigma: float = 1.0,
) -> List[Dict[str, float]]:
    """CAMISIM-like log-normal abundances (one dict per sample)."""
    ids = [str(g) for g in genome_ids]
    if not ids:
        raise ValueError("No genomes for community design")
    rng = np.random.default_rng(int(seed))
    n = len(ids)
    samples: List[Dict[str, float]] = []
    prev = None
    kind = str(mode or "differential").strip().lower()
    for _ in range(max(int(n_samples), 1)):
        if kind == "replicates" and prev is not None:
            noise = rng.normal(gauss_mu, max(gauss_sigma, 1e-9), size=n)
            vals = np.clip(prev + np.abs(noise) * 0.05, 1e-12, None)
        elif kind in {"timeseries_lognormal", "timeseries_normal"} and prev is not None:
            draw = rng.lognormal(log_mu, max(log_sigma, 1e-9), size=n)
            vals = (prev + draw) / 2.0
        else:
            vals = rng.lognormal(log_mu, max(log_sigma, 1e-9), size=n)
        vals = vals / float(vals.sum())
        prev = vals
        samples.append({gid: float(v) for gid, v in zip(ids, vals)})
    return samples


def add_host_abundance(
    samples: Sequence[Dict[str, float]],
    rows: Sequence[Dict[str, str]],
    host_fraction: Any,
    *,
    seed: int = 42,
) -> List[Dict[str, float]]:
    """Insert host genomes at the requested fraction for CAMISIM read modes."""
    host_ids = [str(row["genome_ID"]) for row in rows if row.get("host") == "1"]
    if not host_ids:
        return [dict(sample) for sample in samples]
    random_host = str(host_fraction).strip().upper() == "RANDOM"
    rng = random.Random(int(seed) + 7919)
    result: List[Dict[str, float]] = []
    for sample in samples:
        fraction = rng.random() if random_host else float(host_fraction)
        fraction = min(max(fraction, 0.0), 1.0)
        meta_total = sum(max(float(value), 0.0) for value in sample.values())
        scaled = {
            str(gid): (
                max(float(value), 0.0) / meta_total * (1.0 - fraction)
                if meta_total > 0
                else 0.0
            )
            for gid, value in sample.items()
        }
        per_host = fraction / len(host_ids)
        for host_id in host_ids:
            scaled[host_id] = per_host
        result.append(scaled)
    return result


def filter_wgsim_zero_read_taxa(
    samples: Sequence[Dict[str, float]],
    rows: Sequence[Dict[str, str]],
    *,
    size_gbp: float,
    read_length: int,
) -> List[Dict[str, float]]:
    """Drop taxa for which CAMISIM would invoke wgsim with ``-N 0``.

    CAMISIM length-weights the supplied abundances before calculating
    ``N = size * 1e9 * abundance / read_length``. Its downstream SAM converter
    crashes on wgsim's empty output, so remove only taxa that round to zero and
    renormalize the remaining community.
    """
    genome_sizes: Dict[str, int] = {}
    for row in rows:
        size = 0
        try:
            with open_text(row["path"]) as handle:
                for line in handle:
                    if not line.startswith(">"):
                        size += len(line.strip())
        except OSError:
            size = 0
        genome_sizes[str(row["genome_ID"])] = size

    total_pairs = float(size_gbp) * 1e9 / max(int(read_length), 1)
    filtered: List[Dict[str, float]] = []
    for sample_index, sample in enumerate(samples, start=1):
        weighted = {
            str(gid): max(float(abundance), 0.0) * genome_sizes.get(str(gid), 0)
            for gid, abundance in sample.items()
        }
        weighted_total = sum(weighted.values())
        if weighted_total <= 0:
            filtered.append(dict(sample))
            continue
        keep = {
            gid
            for gid, value in weighted.items()
            if total_pairs * value / weighted_total >= 0.5
        }
        if not keep:
            keep = {max(weighted, key=weighted.get)}
        omitted = sorted(set(map(str, sample)) - keep)
        if omitted:
            logger.warning(
                "wgsim sample %s: omitting zero-read genome IDs %s",
                sample_index,
                ", ".join(omitted),
            )
        kept_total = sum(
            max(float(value), 0.0)
            for gid, value in sample.items()
            if str(gid) in keep
        )
        filtered.append(
            {
                str(gid): max(float(value), 0.0) / kept_total
                for gid, value in sample.items()
                if str(gid) in keep and kept_total > 0
            }
        )
    return filtered


def write_distribution_files(
    dest_dir: Path,
    samples: Sequence[Dict[str, float]],
) -> List[Path]:
    """CAMISIM 2 expects TSV with no header; sample id is ``distribution_<id>``."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    for i, table in enumerate(samples):
        path = dest_dir / f"distribution_{i}.txt"
        lines = [f"{gid}\t{float(abund)}" for gid, abund in table.items()]
        path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
        paths.append(path)
    return paths


def read_distribution_files(directory: Path) -> List[Dict[str, float]]:
    files = sorted(directory.glob("distribution_*.txt"))
    if not files:
        files = sorted(directory.glob("distribution_*"))
    samples: List[Dict[str, float]] = []
    for path in files:
        table: Dict[str, float] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.lower().startswith("genome"):
                continue
            if "," in line:
                gid, abund = line.split(",", 1)
            elif "\t" in line:
                gid, abund = line.split("\t", 1)
            else:
                parts = line.split()
                if len(parts) < 2:
                    continue
                gid, abund = parts[0], parts[1]
            try:
                table[gid.strip()] = float(abund)
            except ValueError:
                continue
        if table:
            samples.append(table)
    return samples


def sanitize_fasta_for_camisim(src: Path, dest: Path, genome_id: str) -> Path:
    """Rewrite FASTA IDs so CAMISIM ``wgsim_to_sam.py`` can split on ``_``.

    wgsim names reads ``{seqid}_{start}_{end}_...``. Headers like
    ``Ecoli.fna|taxid:562|-NC_000913.3`` add extra underscores and crash the
    SAM converter (expected 6 fields, got 7).
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    token = re.sub(r"[^A-Za-z0-9]", "", str(genome_id)) or "g"
    n = 0

    def chunks() -> Iterable[str]:
        nonlocal n
        with open_text(src) as handle:
            for line in handle:
                if line.startswith(">"):
                    n += 1
                    yield f">{token}c{n}\n"
                else:
                    yield line

    write_text_lines(dest, chunks())
    return dest


def stage_genomes_for_camisim(
    rows: Sequence[Dict[str, str]],
    dest_dir: Path,
) -> List[Dict[str, str]]:
    staged_dir = dest_dir / "genomes"
    staged_dir.mkdir(parents=True, exist_ok=True)
    staged: List[Dict[str, str]] = []
    for row in rows:
        copy = dict(row)
        src = Path(row["path"])
        dest = staged_dir / f"{row['genome_ID']}.fna"
        sanitize_fasta_for_camisim(src, dest, row["genome_ID"])
        copy["path"] = str(dest)
        staged.append(copy)
    return staged


def write_genome_tables(
    dest_dir: Path,
    rows: Sequence[Dict[str, str]],
    taxdump: Optional[str] = None,
) -> Tuple[Path, Path]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    loc = dest_dir / "genome_locations.tsv"
    meta = dest_dir / "meta_data.tsv"
    loc.write_text(
        "".join(f"{r['genome_ID']}\t{r['path']}\n" for r in rows),
        encoding="utf-8",
    )
    dump_ids = load_taxdump_taxids(taxdump)
    meta_lines = ["genome_ID\tOTU\tNCBI_ID\tnovelty_category\n"]
    for r in rows:
        ncbi = map_ncbi_id_for_taxdump(str(r["NCBI_ID"]), dump_ids)
        meta_lines.append(
            f"{r['genome_ID']}\t{r['OTU']}\t{ncbi}\t{r['novelty_category']}\n"
        )
    meta.write_text("".join(meta_lines), encoding="utf-8")
    return loc, meta


def default_generate_config(
    *,
    genome_dir: str,
    output_dir: str,
    host_genome: str,
    n_samples: int = 10,
    total_reads: int = 2000,
    host_fraction: str = "RANDOM",
    seed: int = 42,
    model: str = "hiseq",
    cores: int = 1,
    mode: str = "table",
    max_genomes: Any = None,
) -> Dict[str, Any]:
    mode = normalize_camisim_mode(mode)
    root = discover_camisim()
    types = [] if mode == "table" else types_for_mode(mode, root)
    read_length = 150
    nanosim_read_length = 4508
    sizes = camisim_sizes_by_type(
        int(total_reads),
        types,
        art_read_length=read_length,
        nanosim_read_length=nanosim_read_length,
        wgsim_read_length=read_length,
    )
    return {
        "simulator": "camisim",
        "mode": mode,
        "camisim_types": types,
        "n_samples": int(n_samples),
        "total_reads": int(total_reads),
        "size_gbp": (
            next(iter(sizes.values()))
            if sizes
            else camisim_size_gbp(
                int(total_reads), read_length, min_gbp=1e-8, read_ends=2
            )
        ),
        "size_gbp_by_type": sizes,
        "host_fraction": host_fraction,
        "seed": int(seed),
        "cores": int(cores),
        "iss_model": model,
        "read_length": 150,
        "genome_dir": genome_dir,
        "host_genome": host_genome,
        "output_dir": str(Path(output_dir) / "initial"),
        "camisim_root": root or "",
        "nextflow": discover_nextflow() or "nextflow",
        "ncbi_taxdump_file": camisim_taxdump(root),
        "just_community_design": mode == "table",
        "gsa": False,
        "pooled_gsa": False,
        "anonymization": False,
        "distribution": {
            "mode": "differential",
            "log_mu": 1.0,
            "log_sigma": 2.0,
            "gauss_mu": 1.0,
            "gauss_sigma": 1.0,
        },
        "art": {
            "profile_read_length": 150,
            "fragment_size_mean": 270,
            "fragment_size_sd": 27,
        },
        "nanosim3": {
            "read_length": nanosim_read_length,
            "simulate_fastq_directly": True,
        },
        "wgsim": {
            "profile_read_length": 150,
            "fragment_size_mean": 270,
            "fragment_size_sd": 27,
            "base_error_rate": 0.0,
        },
        "max_genomes": parse_max_genomes(
            max_genomes, default_from_env=max_genomes is None
        ),
    }


def generate_dir(output_dir: str) -> Path:
    return Path(absolute_path(output_dir)) / ".generate"


def config_yaml_path(output_dir: str) -> Path:
    return generate_dir(output_dir) / "configs" / "camisim.yaml"


def write_camisim_configs(cfg: Dict[str, Any], output_dir: str) -> Dict[str, str]:
    """Write the editable YAML plus CAMISIM Nextflow / 1.3 files."""
    base = Path(absolute_path(output_dir))
    cfg_dir = generate_dir(str(base)) / "configs"
    cam_dir = cfg_dir / "camisim"
    cam_dir.mkdir(parents=True, exist_ok=True)
    rows = select_camisim_genome_rows(cfg)
    cfg = dict(cfg)
    if not cfg.get("camisim_root"):
        cfg["camisim_root"] = discover_camisim() or ""
    if not cfg.get("ncbi_taxdump_file"):
        cfg["ncbi_taxdump_file"] = camisim_taxdump(cfg.get("camisim_root"))
    if normalize_camisim_mode(cfg.get("mode")) != "table":
        rows = stage_genomes_for_camisim(rows, cam_dir)
    loc, meta = write_genome_tables(cam_dir, rows, taxdump=cfg.get("ncbi_taxdump_file"))
    cfg["genome_locations_file"] = str(loc)
    cfg["metadata_file"] = str(meta)
    cfg["n_genomes"] = len(rows)
    cfg["genomes"] = [
        r["path"] for r in rows if str(r.get("host") or "") != "1"
    ]
    yaml_path = cfg_dir / "camisim.yaml"
    yaml_path.write_text(
        "# SamovaR CAMISIM generate config. Edit and re-run .generate/generate.sh\n"
        + yaml.safe_dump(cfg, sort_keys=False),
        encoding="utf-8",
    )
    nf_path = write_nextflow_config(cfg, cam_dir, rows)
    ini_path = write_legacy_ini(cfg, cam_dir, rows)
    return {
        "yaml": str(yaml_path),
        "nextflow": str(nf_path),
        "legacy_ini": str(ini_path),
        "genome_locations": str(loc),
        "metadata": str(meta),
    }


def write_nextflow_config(
    cfg: Dict[str, Any],
    cam_dir: Path,
    rows: Sequence[Dict[str, str]],
) -> Path:
    root = str(cfg.get("camisim_root") or "")
    types = cfg.get("camisim_types") or []
    cam_type = types[0] if types else "art"
    just = bool(cfg.get("just_community_design"))
    n_gen = int(cfg.get("n_genomes") or len(rows) or 1)
    dist = cfg.get("distribution") or {}
    art = cfg.get("art") or {}
    nano = cfg.get("nanosim3") or {}
    wgsim = cfg.get("wgsim") or {}
    include_root = root.replace("\\", "/")
    # CAMISIM only honours just_community_design when distribution_files is empty.
    if just:
        dist_glob = ""
    else:
        dist_glob = str((cam_dir / "distributions" / "distribution_*.txt")).replace("\\", "/")
    taxdump = str(cfg.get("ncbi_taxdump_file") or camisim_taxdump(root) or "").replace("\\", "/")
    if cam_type == "wgsim":
        profile_len = int(wgsim.get("profile_read_length", 150))
        frag_mean = int(wgsim.get("fragment_size_mean", 270))
        frag_sd = int(wgsim.get("fragment_size_sd", 27))
    else:
        profile_len = int(art.get("profile_read_length", 150))
        frag_mean = int(art.get("fragment_size_mean", 270))
        frag_sd = int(art.get("fragment_size_sd", 27))
    art_profile = str(
        art.get("base_profile_name")
        or (f"{include_root}/tools/art_illumina-2.3.6/profiles/ART_MBARC-26_HiSeq_R" if include_root else "")
    ).replace("\\", "/")
    nano_profile = str(
        nano.get("base_profile_name")
        or (f"{include_root}/tools/nanosim_profile/nanosim323_r10_hc_lomanLab" if include_root else "")
    ).replace("\\", "/")
    conda_cache = str((user_cache_dir() / "camisim-conda").resolve()).replace("\\", "/")
    try:
        Path(conda_cache).mkdir(parents=True, exist_ok=True)
    except OSError:
        conda_cache = ""
    includes = []
    if include_root:
        includes.append(f'includeConfig "{include_root}/pipelines/shared/config/conda.config"')
        type_cfg = {
            "art": f"{include_root}/pipelines/metagenomic/config/art.config",
            "nanosim3": f"{include_root}/pipelines/metagenomic/config/nanosim.config",
            "wgsim": f"{include_root}/pipelines/metagenomic/config/wgsim.config",
        }.get(cam_type)
        if type_cfg:
            includes.append(f'includeConfig "{type_cfg}"')
    include_block = "\n".join(includes)
    size_map = cfg.get("size_gbp_by_type") or {}
    size = float(
        size_map.get(cam_type)
        or cfg.get("size_gbp")
        or camisim_sizes_by_type(
            int(cfg.get("total_reads") or 2000),
            [cam_type],
            art_read_length=int(art.get("profile_read_length", 150)),
            nanosim_read_length=int(nano.get("read_length", 4508)),
            wgsim_read_length=int(wgsim.get("profile_read_length", 150)),
        ).get(cam_type)
        or camisim_size_gbp(int(cfg.get("total_reads") or 2000), profile_len)
    )
    body = f"""// Generated by samovar generate --simulator camisim
// Exclusive Nextflow --config: include conda + simulator, then override params.
{include_block}

params {{
  outdir = "{Path(cfg['output_dir']).as_posix()}/camisim_raw"
  size = {size}
  type = "{cam_type}"
  number_of_samples = {int(cfg.get('n_samples') or 1)}
  gsa = {str(bool(cfg.get('gsa'))).lower()}
  pooled_gsa = {str(bool(cfg.get('pooled_gsa'))).lower()}
  anonymization = {str(bool(cfg.get('anonymization'))).lower()}
  seed = {int(cfg.get('seed') or 42)}
  biom_profile = ""
  genome_locations_file = "{Path(cfg['genome_locations_file']).as_posix()}"
  metadata_file = "{Path(cfg['metadata_file']).as_posix()}"
  ncbi_taxdump_file = "{taxdump}"
  max_strains_per_otu = 1
  min_strains_per_otu = 1
  distribution_files = "{dist_glob}"
  just_community_design = {str(just).lower()}
  mode = "{dist.get('mode', 'differential')}"
  log_mu = {float(dist.get('log_mu', 1))}
  log_sigma = {float(dist.get('log_sigma', 2))}
  gauss_mu = {float(dist.get('gauss_mu', 1))}
  gauss_sigma = {float(dist.get('gauss_sigma', 1))}
  genomes_total = {n_gen}
  genomes_real = {n_gen}
  id_to_gff_file = ""
  strain_simulation_template = ""
  verbose = "False"
  profile_read_length = {profile_len}
  fragment_size_mean = {frag_mean}
  fragment_size_sd = {frag_sd}
  base_error_rate = {float(wgsim.get('base_error_rate', 0.0))}
  create_cigar = false
  read_length = {int(nano.get('read_length', 4508))}
  simulate_fastq_directly = {str(bool(nano.get('simulate_fastq_directly'))).lower()}
  base_profile_name = "{art_profile if cam_type == 'art' else nano_profile}"
}}

process.time = '2h'
"""
    if conda_cache:
        body += f'\nconda.cacheDir = "{conda_cache}"\n'
    nano_env = (
        tool_env_prefix("nanosim")
        or tool_env_prefix("nanosim3")
        or conda_prefix_for_executable(discover_nanosim())
    )
    art_env = (
        tool_env_prefix("art")
        or tool_env_prefix("art_illumina")
        or conda_prefix_for_executable(discover_art())
    )
    overrides = []
    if nano_env:
        overrides.append(
            f"""
process {{
  withName: /simulate_reads_.*nanosim3/ {{
    conda = "{Path(nano_env).as_posix()}"
  }}
}}
"""
        )
    if art_env:
        overrides.append(
            f"""
process {{
  withName: 'simulate_reads_art' {{
    conda = "{Path(art_env).as_posix()}"
  }}
}}
"""
        )
    if overrides:
        body += "\n" + "\n".join(overrides)
    path = cam_dir / "nextflow.config"
    path.write_text(body, encoding="utf-8")
    return path


def write_legacy_ini(
    cfg: Dict[str, Any],
    cam_dir: Path,
    rows: Sequence[Dict[str, str]],
) -> Path:
    n_gen = int(cfg.get("n_genomes") or len(rows) or 1)
    dist = cfg.get("distribution") or {}
    art = cfg.get("art") or {}
    types = cfg.get("camisim_types") or []
    cam_type = types[0] if types else "art"
    taxdump = str(cfg.get("ncbi_taxdump_file") or camisim_taxdump(cfg.get("camisim_root")) or "")
    ini = f"""[Main]
seed={int(cfg.get('seed') or 42)}
phase=0
max_processors={int(cfg.get('cores') or 1)}
dataset_id=samovar
output_directory={Path(cfg['output_dir']) / 'camisim_raw'}
temp_directory={cam_dir / 'tmp'}
gsa={bool(cfg.get('gsa'))}
pooled_gsa={bool(cfg.get('pooled_gsa'))}
anonymous={bool(cfg.get('anonymization'))}
compress=0

[ReadSimulator]
readsim=art_illumina
error_profiles=
samtools=samtools
profile=mbarc
size={float(cfg.get('size_gbp') or 0.01)}
type={cam_type}
fragments_size_mean={int(art.get('fragment_size_mean', 270))}
fragment_size_standard_deviation={int(art.get('fragment_size_sd', 27))}

[CommunityDesign]
ncbi_taxdump={taxdump}
strain_simulation_template=
number_of_samples={int(cfg.get('n_samples') or 1)}

[community0]
metadata={cfg['metadata_file']}
id_to_genome_file={cfg['genome_locations_file']}
id_to_gff_file=
genomes_total={n_gen}
genomes_real={n_gen}
max_strains_per_otu=1
ratio=1
mode={dist.get('mode', 'differential')}
log_mu={float(dist.get('log_mu', 1))}
log_sigma={float(dist.get('log_sigma', 2))}
gauss_mu={float(dist.get('gauss_mu', 1))}
gauss_sigma={float(dist.get('gauss_sigma', 1))}
view=False
"""
    path = cam_dir / "mini_config.ini"
    path.write_text(ini, encoding="utf-8")
    return path


def write_generate_pipeline(output_dir: str, cfg: Dict[str, Any]) -> str:
    from samovar.paths import python_path as py_path
    from samovar.paths import shell_outdir_override_snippet, shell_source_install_env_snippet

    base = Path(absolute_path(output_dir))
    generate = generate_dir(str(base))
    generate.mkdir(parents=True, exist_ok=True)
    yaml_path = config_yaml_path(str(base))
    root = repo_root()
    py = py_path()
    tool_path = runtime_path_prefix()
    env_snippet = shell_source_install_env_snippet()
    pipeline = generate / "generate.sh"
    pipeline.write_text(
        f"""# Setup
set -e
export SAMOVAR_ROOT="{root}"
{env_snippet}export PATH="${{SAMOVAR_PATH:+$SAMOVAR_PATH:}}{tool_path}:{root}/bin:$PATH"
export PYTHONPATH="{root / 'src'}${{PYTHONPATH:+:$PYTHONPATH}}"
PYTHON_PATH="${{PYTHON_PATH:-{py}}}"
if [ -z "$PYTHON_PATH" ] || [ ! -x "$PYTHON_PATH" ]; then
  PYTHON_PATH="$(command -v python3 || command -v python || true)"
fi
PYTHON_PATH="${{PYTHON_PATH:-python3}}"

out_dir="{base}"
{shell_outdir_override_snippet()}
mkdir -p "$out_dir"

# CAMISIM community / reads (table mode continues with ISS)
"$PYTHON_PATH" -m samovar.camisim --config "{yaml_path}"
""",
        encoding="utf-8",
    )
    os.chmod(pipeline, 0o755)
    return str(pipeline)


def setup_camisim_generate(args: argparse.Namespace) -> Dict[str, str]:
    """Write CAMISIM configs + generate.sh from ``samovar generate`` args."""
    mode = normalize_camisim_mode(getattr(args, "camisim_mode", None) or "table")
    existing = getattr(args, "camisim_config", None)
    if existing:
        cfg = yaml.safe_load(Path(existing).read_text(encoding="utf-8")) or {}
    else:
        cfg = {}
    base = default_generate_config(
        genome_dir=absolute_path(args.genome_dir),
        output_dir=absolute_path(args.output_dir),
        host_genome=absolute_path(args.host_genome),
        n_samples=args.n_samples if args.n_samples is not None else cfg.get("n_samples", 10),
        total_reads=args.total_reads if args.total_reads is not None else cfg.get("total_reads", 2000),
        host_fraction=args.host_fraction if args.host_fraction is not None else cfg.get("host_fraction", "RANDOM"),
        seed=args.seed if args.seed is not None else cfg.get("seed", 42),
        model=args.model if args.model is not None else cfg.get("iss_model", "hiseq"),
        cores=getattr(args, "cores", None) if getattr(args, "cores", None) is not None else cfg.get("cores", 1),
        mode=mode,
        max_genomes=getattr(args, "max_genomes", None),
    )
    # CLI/defaults win for scalars; nested dicts from --camisim-config overlay.
    merged = _deep_overlay(base, cfg if isinstance(cfg, dict) else {})
    locked = getattr(args, "genomes", None)
    if locked:
        from samovar.iss_config import parse_genome_tokens, select_generate_genome_paths

        merged["genomes"] = [
            str(path)
            for path in select_generate_genome_paths(
                merged.get("genome_dir") or "",
                merged.get("max_genomes"),
                parse_genome_tokens(locked),
            )
        ]
    if not merged.get("camisim_root"):
        merged["camisim_root"] = discover_camisim() or ""
    if not merged.get("ncbi_taxdump_file"):
        merged["ncbi_taxdump_file"] = camisim_taxdump(merged.get("camisim_root"))
    if getattr(args, "size_gbp", None):
        merged["size_gbp"] = float(args.size_gbp)
        merged["size_gbp_by_type"] = {
            cam_type: float(args.size_gbp)
            for cam_type in (merged.get("camisim_types") or [])
        }
    paths = write_camisim_configs(merged, args.output_dir)
    pipeline = write_generate_pipeline(args.output_dir, merged)
    return {"config": paths["yaml"], "pipeline": pipeline, **paths}


def _load_run_config(path: Path) -> Dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"CAMISIM config must be a mapping: {path}")
    return data


def run_from_config(config_path: str) -> Dict[str, Any]:
    """Execute community design and/or CAMISIM/ISS from the editable YAML."""
    cfg_path = Path(config_path)
    cfg = _load_run_config(cfg_path)
    if not cfg.get("camisim_root"):
        cfg["camisim_root"] = discover_camisim() or ""
    if not cfg.get("ncbi_taxdump_file"):
        cfg["ncbi_taxdump_file"] = camisim_taxdump(cfg.get("camisim_root"))
    out_dir = Path(cfg.get("output_dir") or (cfg_path.parents[2] / "initial"))
    out_dir.mkdir(parents=True, exist_ok=True)
    work = cfg_path.parent / "camisim"
    work.mkdir(parents=True, exist_ok=True)
    rows = select_camisim_genome_rows(cfg)
    if normalize_camisim_mode(cfg.get("mode")) != "table":
        rows = stage_genomes_for_camisim(rows, work)
    loc, meta = write_genome_tables(work, rows, taxdump=cfg.get("ncbi_taxdump_file"))
    cfg["genome_locations_file"] = str(loc)
    cfg["metadata_file"] = str(meta)
    cfg["n_genomes"] = len(rows)
    mode = normalize_camisim_mode(cfg.get("mode"))
    dist_dir = work / "distributions"
    design_rows = [r for r in rows if r.get("host") != "1"] or list(rows)
    samples = design_communities(
        [r["genome_ID"] for r in design_rows],
        int(cfg.get("n_samples") or 1),
        seed=int(cfg.get("seed") or 42),
        mode=str((cfg.get("distribution") or {}).get("mode") or "differential"),
        log_mu=float((cfg.get("distribution") or {}).get("log_mu", 1)),
        log_sigma=float((cfg.get("distribution") or {}).get("log_sigma", 2)),
        gauss_mu=float((cfg.get("distribution") or {}).get("gauss_mu", 1)),
        gauss_sigma=float((cfg.get("distribution") or {}).get("gauss_sigma", 1)),
    )
    if mode != "table":
        samples = add_host_abundance(
            samples,
            rows,
            cfg.get("host_fraction", "RANDOM"),
            seed=int(cfg.get("seed") or 42),
        )
    if mode == "wgsim":
        wgsim_cfg = cfg.get("wgsim") or {}
        wgsim_size = float(
            (cfg.get("size_gbp_by_type") or {}).get("wgsim")
            or cfg.get("size_gbp")
        )
        samples = filter_wgsim_zero_read_taxa(
            samples,
            rows,
            size_gbp=wgsim_size,
            read_length=int(wgsim_cfg.get("profile_read_length", 150)),
        )
    write_distribution_files(dist_dir, samples)
    abundance_csv = work / "abundance_table.csv"
    _write_abundance_csv(abundance_csv, samples, rows)

    if mode == "table":
        # Native design only. CAMISIM Nextflow ignores just_community_design when
        # distribution_files is set, and would otherwise start a full ART run.
        outputs = _iss_from_distributions(cfg, rows, samples, out_dir)
        return {
            "mode": mode,
            "abundance": str(abundance_csv),
            "reads": outputs,
            "used_iss": True,
        }

    types = cfg.get("camisim_types") or types_for_mode(mode, cfg.get("camisim_root") or None)
    if not types:
        raise RuntimeError(f"No CAMISIM read types for mode {mode!r}")
    harvested: Dict[str, List[Path]] = {}
    used_camisim = False
    for cam_type in types:
        raw = work / f"raw_{cam_type}"
        ok = _run_camisim_reads(cfg, work, cam_type, raw)
        if not ok:
            log = raw / "camisim.nextflow.log"
            raise RuntimeError(
                f"CAMISIM {cam_type} simulation failed"
                + (f" (see {log})" if log.is_file() else "")
                + f". Need a CAMISIM checkout ({CAMISIM_REPO}) with Nextflow and "
                f"{CAMISIM_TYPES[cam_type]['label']}, or use --camisim-mode table."
            )
        used_camisim = True
        read_type = CAMISIM_TYPES[cam_type]["read_type"]
        harvested[read_type] = _harvest_camisim_reads(
            raw, out_dir, rows, read_type, int(cfg.get("n_samples") or 1)
        )
    if mode == "hybrid":
        outputs = _merge_hybrid_reads(harvested, out_dir, int(cfg.get("n_samples") or 1))
    else:
        read_type = next(iter(harvested))
        outputs = _promote_to_full(harvested[read_type], out_dir, int(cfg.get("n_samples") or 1), read_type)
    return {
        "mode": mode,
        "abundance": str(abundance_csv),
        "reads": outputs,
        "used_iss": False,
        "used_camisim": used_camisim,
        "read_types": list(harvested),
    }


def _write_abundance_csv(
    dest: Path,
    samples: Sequence[Dict[str, float]],
    rows: Sequence[Dict[str, str]],
) -> None:
    ids = [r["genome_ID"] for r in rows]
    header = ["taxid", "genome_ID"] + [f"N_{i+1}" for i in range(len(samples))]
    lines = [",".join(header)]
    for row in rows:
        gid = row["genome_ID"]
        vals = [str(samples[i].get(gid, 0.0) if i < len(samples) else 0.0) for i in range(len(samples))]
        lines.append(",".join([row["NCBI_ID"], gid, *vals]))
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _try_camisim_community(cfg: Dict[str, Any], work: Path) -> bool:
    """Best-effort CAMISIM community design; native tables already exist."""
    root = cfg.get("camisim_root") or discover_camisim()
    if not root:
        logger.info("CAMISIM not installed; using SamovaR log-normal community design")
        return False
    cfg = dict(cfg)
    cfg["camisim_root"] = root
    cfg["just_community_design"] = True
    try:
        return _run_camisim_reads(cfg, work, "art", work / "community_only", community_only=True)
    except Exception as exc:
        logger.warning("CAMISIM community design skipped (%s); using native abundances", exc)
        return False


def nextflow_run_cmd(
    nextflow_bin: str,
    camisim_root: str,
    config_path: Path,
    work_dir: Path,
) -> List[str]:
    """CAMISIM 2: ``nextflow run main.nf --config /abs/ours.config`` from the clone.

    Passing ``-c`` *and* ``--config`` to the same file skips ``metagenomic.config``
    twice; ``params.config`` already exclusive-includes that one file.
    """
    return [
        str(nextflow_bin),
        "run",
        "main.nf",
        "--config",
        str(Path(config_path).resolve()),
        "-work-dir",
        str(Path(work_dir).resolve()),
        "-resume",
    ]


def _run_camisim_reads(
    cfg: Dict[str, Any],
    work: Path,
    cam_type: str,
    dest: Path,
    community_only: bool = False,
) -> bool:
    root = cfg.get("camisim_root") or discover_camisim()
    if not root:
        return False
    dest.mkdir(parents=True, exist_ok=True)
    run_cfg = dict(cfg)
    run_cfg["camisim_root"] = root
    run_cfg["camisim_types"] = [cam_type]
    run_cfg["just_community_design"] = bool(community_only)
    run_cfg["output_dir"] = str(dest)
    if not run_cfg.get("ncbi_taxdump_file"):
        run_cfg["ncbi_taxdump_file"] = camisim_taxdump(root)
    rows = select_camisim_genome_rows(cfg)
    nf = write_nextflow_config(run_cfg, work, rows)
    kind = _camisim_kind(Path(root))
    log_path = dest / "camisim.nextflow.log"
    if kind == "nextflow":
        nxt = discover_nextflow() or run_cfg.get("nextflow") or "nextflow"
        cmd = nextflow_run_cmd(str(nxt), str(root), nf, work / f"nf_work_{cam_type}")
        extra = extra_flags_argv(
            run_cfg.get("extra_flags") or run_cfg.get("reads_generator_flags")
        )
        cmd.extend(extra)
        logger.info("Running CAMISIM: %s (cwd=%s)", " ".join(cmd), root)
        try:
            with open(log_path, "w", encoding="utf-8") as log:
                log.write("cwd: %s\ncmd: %s\n\n" % (root, " ".join(cmd)))
                log.flush()
                result = subprocess.run(
                    cmd,
                    cwd=str(root),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if result.returncode != 0:
                tail = _tail_text(log_path, 40)
                logger.warning(
                    "CAMISIM Nextflow failed (exit %s). Log: %s\n%s",
                    result.returncode,
                    log_path,
                    tail,
                )
                return False
            return True
        except OSError as exc:
            logger.warning("CAMISIM Nextflow failed: %s", exc)
            return False
    if kind == "legacy":
        script = Path(root) / "metagenomesimulation.py"
        if not script.is_file():
            script = Path(root) / "scripts" / "metagenomesimulation.py"
        ini = write_legacy_ini(run_cfg, work, rows)
        cmd = [python_path(), str(script), str(ini)]
        logger.info("Running CAMISIM 1.3: %s", " ".join(cmd))
        try:
            with open(log_path, "w", encoding="utf-8") as log:
                result = subprocess.run(
                    cmd,
                    cwd=str(root),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if result.returncode != 0:
                logger.warning("CAMISIM 1.3 failed (exit %s). Log: %s", result.returncode, log_path)
                return False
            return True
        except OSError as exc:
            logger.warning("CAMISIM 1.3 failed: %s", exc)
            return False
    return False


def _tail_text(path: Path, n_lines: int = 40) -> str:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-n_lines:])


def _iter_fastq_files(directory: Path) -> List[Path]:
    hits: List[Path] = []
    if not directory.is_dir():
        return hits
    for path in directory.rglob("*"):
        try:
            if path.is_file() and is_fastq_name(path.name):
                hits.append(path)
        except OSError:
            continue
    return hits


_POOLED_PAIRED = re.compile(
    r"^sample_(\d+)_0([12])\.(?:fq|fastq)(?:\.gz|\.bz2|\.xz)?$",
    re.IGNORECASE,
)
_POOLED_SINGLE = re.compile(
    r"^sample_(\d+)\.(?:fq|fastq)(?:\.gz|\.bz2|\.xz)?$",
    re.IGNORECASE,
)
_ANON = re.compile(r"anonymous_reads", re.IGNORECASE)


def _is_pooled_fastq(path: Path) -> bool:
    name = path.name
    return bool(_POOLED_PAIRED.match(name) or _POOLED_SINGLE.match(name) or _ANON.search(name))


def _sample_index_from_name(name: str, n_samples: int) -> Optional[int]:
    m = _POOLED_PAIRED.match(name) or _POOLED_SINGLE.match(name)
    if m:
        idx = int(m.group(1))
        if 0 <= idx < n_samples:
            return idx
    lower = name.lower()
    for i in range(n_samples):
        tokens = (f"sample_{i}_", f"sample_{i}.", f"sample{i}_")
        if any(tok in lower for tok in tokens):
            return i
    return None


def _sample_fastq_groups(raw_dir: Path, n_samples: int) -> List[List[Path]]:
    """Bucket CAMISIM FASTQs into n samples, preferring pooled files.

    CAMISIM writes per-genome ``sample{i}_{gid}1.fq.gz`` *and* pooled
    ``sample_{i}_01.fq.gz`` in the same ``reads/fastq`` folder. Concatenating
    both duplicates every read.
    """
    files = _iter_fastq_files(raw_dir)
    groups: List[List[Path]] = [[] for _ in range(max(n_samples, 1))]
    if not files:
        return groups
    pooled = [p for p in files if _is_pooled_fastq(p)]
    chosen = pooled if pooled else files
    for path in chosen:
        idx = _sample_index_from_name(path.name, n_samples)
        if idx is not None:
            groups[idx].append(path)
        elif n_samples == 1:
            groups[0].append(path)
    if not any(groups) and chosen:
        groups[0] = chosen
    return groups


def _mate_kind(path: Path, rows: Optional[Sequence[Dict[str, str]]] = None) -> str:
    name = path.name
    pooled = _POOLED_PAIRED.match(name)
    if pooled:
        return "R2" if pooled.group(2) == "2" else "R1"
    lower = name.lower()
    if re.search(r"(?:_r2|r2|_02)(?:\.(?:fq|fastq))", lower):
        return "R2"
    if rows:
        for gid in sorted((r["genome_ID"] for r in rows), key=len, reverse=True):
            if f"{gid}2.fq" in name or f"{gid}2.fastq" in name:
                return "R2"
            if f"{gid}1.fq" in name or f"{gid}1.fastq" in name:
                return "R1"
    if re.search(r"_2\.(?:fq|fastq)", lower) or lower.endswith("_2.fq") or lower.endswith("_2.fastq"):
        return "R2"
    return "R1"


def _taxid_for_filename(path: Path, rows: Sequence[Dict[str, str]]) -> str:
    by_id = {r["genome_ID"]: r["NCBI_ID"] for r in rows}
    name = path.name
    for gid, taxid in sorted(by_id.items(), key=lambda kv: len(kv[0]), reverse=True):
        if gid in name:
            return str(taxid)
    tax = taxid_from_fasta_name(path)
    return tax or ""


def _taxid_for_read_id(header_or_id: str, rows: Sequence[Dict[str, str]]) -> str:
    """CAMISIM pooled FASTQs omit genome_ID in the filename; wgsim IDs start with it."""
    raw = str(header_or_id or "")
    body = raw[1:] if raw.startswith("@") else raw
    rid = body.split(None, 1)[0] if body else ""
    for row in sorted(rows, key=lambda r: len(str(r["genome_ID"])), reverse=True):
        gid = str(row["genome_ID"])
        if rid == gid or rid.startswith(gid):
            return str(row["NCBI_ID"])
    return ""


def tag_fastq_file(
    src: Path,
    dest: Path,
    read_type: str,
    taxid: str = "",
    rows: Optional[Sequence[Dict[str, str]]] = None,
) -> None:
    src = Path(src)
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    same_file = src.resolve() == dest.resolve()
    write_to = dest.with_name(dest.name + ".tagging.tmp") if same_file else dest

    def chunks() -> Iterable[str]:
        for header, seq, plus, qual in iter_fastq_records(src):
            rec_tax = taxid or (_taxid_for_read_id(header, rows) if rows else "")
            yield tag_read_id(header, read_type, rec_tax) + "\n"
            yield seq
            yield plus
            yield qual

    write_text_lines(write_to, chunks())
    if same_file:
        write_to.replace(dest)


def _harvest_camisim_reads(
    raw_dir: Path,
    out_dir: Path,
    rows: Sequence[Dict[str, str]],
    read_type: str,
    n_samples: int,
) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    groups = _sample_fastq_groups(raw_dir, n_samples)
    written: List[Path] = []
    for i, files in enumerate(groups, start=1):
        r1_src = [p for p in files if _mate_kind(p, rows) == "R1"]
        r2_src = [p for p in files if _mate_kind(p, rows) == "R2"]
        r1_dest = out_dir / f"{i}_{read_type}_R1.fastq"
        r2_dest = out_dir / f"{i}_{read_type}_R2.fastq"
        _concat_tagged(r1_src, r1_dest, read_type, rows)
        if r2_src:
            _concat_tagged(r2_src, r2_dest, read_type, rows)
        else:
            _write_empty_fastq(r2_dest)
        written.extend([r1_dest, r2_dest])
    return written


def _concat_tagged(
    sources: Sequence[Path],
    dest: Path,
    read_type: str,
    rows: Sequence[Dict[str, str]],
) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not sources:
        _write_empty_fastq(dest)
        return
    tmp_dir = dest.parent / ".camisim_tag"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tagged = []
    for i, src in enumerate(sources):
        part = tmp_dir / f"{dest.stem}_{i}.fastq"
        tag_fastq_file(src, part, read_type, _taxid_for_filename(src, rows), rows=rows)
        tagged.append(part)
    concat_fastq_files(tagged, dest)
    for part in tagged:
        try:
            part.unlink()
        except OSError:
            pass


def _write_empty_fastq(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")


_TECH_FASTQ = re.compile(
    r"^(\d+)_(illumina|ont|wgsim|nanosim3|art)_(R[12])\.fastq(?:\.gz)?$",
    re.IGNORECASE,
)


def _remove_tech_fastq_copies(out_dir: Path) -> None:
    """Drop ``{n}_{tech}_R*.fastq`` after they have been copied to ``{n}_full_``."""
    if not out_dir.is_dir():
        return
    for path in list(out_dir.iterdir()):
        if path.is_file() and _TECH_FASTQ.match(path.name):
            try:
                path.unlink()
            except OSError:
                pass


def _merge_hybrid_reads(
    harvested: Dict[str, List[Path]],
    out_dir: Path,
    n_samples: int,
) -> List[str]:
    """Merge technologies without inventing mates for single-end ONT reads.

    A mixed paired-end/single-end library cannot be passed to classifiers as a
    paired FASTQ because R1 and R2 would contain different record counts. If a
    sample includes any single-end technology, flatten every sequence (both
    paired ends plus the long reads) into R1 and leave R2 empty. Annotator
    wrappers then select their single-end CLI form.
    """
    outputs: List[str] = []
    for i in range(1, n_samples + 1):
        by_type = {
            read_type: {
                mate: out_dir / f"{i}_{read_type}_{mate}.fastq"
                for mate in ("R1", "R2")
            }
            for read_type in harvested
        }
        has_single_end = any(
            paths["R1"].is_file()
            and paths["R1"].stat().st_size > 0
            and (
                not paths["R2"].is_file()
                or paths["R2"].stat().st_size == 0
            )
            for paths in by_type.values()
        )
        if has_single_end:
            r1_parts = []
            for paths in by_type.values():
                for mate in ("R1", "R2"):
                    path = paths[mate]
                    if path.is_file() and path.stat().st_size > 0:
                        r1_parts.append(path)
            r1_dest = out_dir / f"{i}_full_R1.fastq"
            r2_dest = out_dir / f"{i}_full_R2.fastq"
            concat_fastq_files(r1_parts, r1_dest)
            _write_empty_fastq(r2_dest)
            outputs.extend([str(r1_dest), str(r2_dest)])
            continue
        for mate in ("R1", "R2"):
            parts = [
                paths[mate]
                for paths in by_type.values()
                if paths[mate].is_file() and paths[mate].stat().st_size > 0
            ]
            dest = out_dir / f"{i}_full_{mate}.fastq"
            if parts:
                concat_fastq_files(parts, dest)
            else:
                _write_empty_fastq(dest)
            outputs.append(str(dest))
    _remove_tech_fastq_copies(out_dir)
    return outputs


def _promote_to_full(
    tagged: Sequence[Path],
    out_dir: Path,
    n_samples: int,
    read_type: str,
) -> List[str]:
    """Copy per-tech tagged FASTQs to the ``{sample}_full_R*.fastq`` names ISS uses."""
    outputs: List[str] = []
    for i in range(1, n_samples + 1):
        for mate in ("R1", "R2"):
            src = out_dir / f"{i}_{read_type}_{mate}.fastq"
            dest = out_dir / f"{i}_full_{mate}.fastq"
            if src.is_file():
                shutil.copy2(src, dest)
            else:
                _write_empty_fastq(dest)
            outputs.append(str(dest))
    _remove_tech_fastq_copies(out_dir)
    return outputs


def _iss_from_distributions(
    cfg: Dict[str, Any],
    rows: Sequence[Dict[str, str]],
    samples: Sequence[Dict[str, float]],
    out_dir: Path,
) -> List[str]:
    from samovar.table2iss import generate_reads_genome, generate_reads_metagenome, iss_cli_extra_flags

    out_dir.mkdir(parents=True, exist_ok=True)
    total_reads = int(cfg.get("total_reads") or 2000)
    host_fraction = cfg.get("host_fraction", "RANDOM")
    seed = int(cfg.get("seed") or 42)
    model = str(cfg.get("iss_model") or "hiseq")
    cpus = int(cfg.get("cores") or 1)
    extra = extra_flags_argv(cfg.get("extra_flags") or cfg.get("reads_generator_flags"))
    rng = random.Random(seed)
    host_row = next((r for r in rows if r.get("host") == "1"), None)
    if host_row is None and cfg.get("host_genome"):
        host_path = Path(cfg["host_genome"])
        host_row = {
            "path": str(host_path),
            "host": "1",
            "NCBI_ID": taxid_from_fasta_name(host_path) or "9606",
            "genome_ID": taxid_from_fasta_name(host_path) or "9606",
        }
    meta_rows = [r for r in rows if r.get("host") != "1"] or list(rows)
    outputs: List[str] = []
    with iss_cli_extra_flags(extra):
        for i, table in enumerate(samples, start=1):
            if str(host_fraction).upper() == "RANDOM":
                host_n = rng.randint(0, total_reads)
            else:
                host_n = int(float(host_fraction) * total_reads)
            meta_n = max(0, total_reads - host_n)
            sample_dir = out_dir / f".iss_s{i}"
            sample_dir.mkdir(parents=True, exist_ok=True)
            host_r1 = sample_dir / "host_R1.fastq"
            host_r2 = sample_dir / "host_R2.fastq"
            if host_n > 0 and host_row and Path(host_row["path"]).is_file():
                generate_reads_genome(
                    host_row["path"],
                    str(sample_dir / "host"),
                    host_n,
                    model=model,
                    cpus=cpus,
                    seed=seed + i,
                )
                host_tax = str(host_row.get("NCBI_ID") or "")
                if host_tax:
                    tag_fastq_file(host_r1, host_r1, "", host_tax)
                    if host_r2.is_file():
                        tag_fastq_file(host_r2, host_r2, "", host_tax)
            else:
                _write_empty_fastq(host_r1)
                _write_empty_fastq(host_r2)
            genome_files = [r["path"] for r in meta_rows]
            genome_ids = [r["NCBI_ID"] for r in meta_rows]
            weights = [float(table.get(r["genome_ID"], 0.0)) for r in meta_rows]
            weight_sum = sum(weights)
            if meta_n > 0 and genome_files:
                if weight_sum > 0:
                    amount = [max(0, int(round(w / weight_sum * meta_n))) for w in weights]
                else:
                    amount = [1] * len(genome_files)
                generate_reads_metagenome(
                    genome_files=genome_files,
                    output_dir=str(sample_dir),
                    amount=amount,
                    total_amount=meta_n,
                    sample_name="meta",
                    annotator_name="full",
                    model=model,
                    genome_ids=genome_ids,
                    cpus=cpus,
                    seed=seed + 100 + i,
                )
            else:
                _write_empty_fastq(sample_dir / "meta_full_R1.fastq")
                _write_empty_fastq(sample_dir / "meta_full_R2.fastq")
            for mate in ("R1", "R2"):
                dest = out_dir / f"{i}_full_{mate}.fastq"
                concat_fastq_files(
                    [sample_dir / f"host_{mate}.fastq", sample_dir / f"meta_full_{mate}.fastq"],
                    dest,
                )
                outputs.append(str(dest))
            shutil.rmtree(sample_dir, ignore_errors=True)
    return outputs


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SamovaR CAMISIM generate config")
    parser.add_argument("--config", required=True, help="Path to .generate/configs/camisim.yaml")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args(argv)
    result = run_from_config(args.config)
    print(json.dumps({k: v for k, v in result.items() if k != "reads"}, indent=2))
    print("reads:", len(result.get("reads") or []))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
