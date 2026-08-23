"""Reuse previously downloaded NCBI genomes instead of fetching them again.

Bundled ``data/test_genomes`` assemblies are truncated (ISS/CI stubs). They
are never treated as an NCBI library. Real metagenomes reuse only genomes
that SamovaR or the user already fetched from NCBI (or listed in
``genome_dirs``).

Config keys (written by ``./install.sh`` and updated on download / prepare):

* ``genomes`` — cache for NCBI FASTA downloads
* ``processed_genomes`` — ``{taxid}-processed.fasta.gz`` (defaults to ``genomes``)
* ``genome_dirs`` — extra directories that already hold full assemblies
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from samovar.paths import (
    genome_download_dir,
    load_config,
    processed_genomes_dir,
    repo_root,
    test_genomes_dir,
    update_config,
)
from samovar.seqio import (
    genome_lookup_extensions,
    is_fasta_name,
    taxid_from_fasta_name,
)

logger = logging.getLogger(__name__)

PathLike = Union[str, os.PathLike]


def _as_path(path: PathLike) -> Path:
    return Path(path).expanduser()


def _normalize_taxid(taxid: Union[str, int]) -> str:
    return str(taxid).split(".")[0]


def _split_dir_value(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        items: List[str] = []
        for part in value:
            items.extend(_split_dir_value(part))
        return items
    text = str(value).strip()
    if not text or text.startswith("_"):
        return []
    return [
        piece.strip()
        for piece in text.replace(";", ":").split(":")
        if piece.strip()
    ]


def reuse_enabled(override: Optional[bool] = None) -> bool:
    if override is not None:
        return bool(override)
    env = os.environ.get("SAMOVAR_REUSE_GENOMES", "").strip().lower()
    if env in {"0", "false", "no", "off"}:
        return False
    if env in {"1", "true", "yes", "on"}:
        return True
    cfg = load_config()
    if "reuse_genomes" in cfg:
        return bool(cfg.get("reuse_genomes"))
    return True


def allow_test_genomes(override: Optional[bool] = None) -> bool:
    if override is not None:
        return bool(override)
    env = os.environ.get("SAMOVAR_ALLOW_TEST_GENOMES", "").strip().lower()
    return env in {"1", "true", "yes", "on"}


def is_bundled_test_genomes_path(path: PathLike) -> bool:
    """True if ``path`` is (or is inside) the truncated ISS/CI test genomes tree."""
    try:
        resolved = _as_path(path).resolve()
    except OSError:
        resolved = _as_path(path)
    roots: List[Path] = []
    try:
        roots.append(test_genomes_dir().resolve())
    except OSError:
        roots.append(test_genomes_dir())
    try:
        roots.append((repo_root() / "data" / "test_genomes").resolve())
    except OSError:
        roots.append(repo_root() / "data" / "test_genomes")
    for root in roots:
        try:
            resolved.relative_to(root)
            return True
        except ValueError:
            continue
    return "test_genomes" in resolved.parts


def dest_name_for_taxid(taxid: Union[str, int], src: PathLike) -> str:
    """Canonical ``{taxid}.fna`` / ``{taxid}-processed.fasta.gz`` from any source name."""
    taxid = _normalize_taxid(taxid)
    src_p = _as_path(src)
    name = src_p.name
    idx = name.lower().rfind("-processed.")
    if idx >= 0:
        return f"{taxid}{name[idx:]}"
    return f"{taxid}{''.join(src_p.suffixes)}"


def place_genome(
    src: PathLike,
    dest_dir: PathLike,
    taxid: Union[str, int],
    *,
    dest_name: Optional[str] = None,
) -> Path:
    """Symlink ``src`` into ``dest_dir``; copy with a warning if linking fails."""
    src_p = _as_path(src)
    dest_p = _as_path(dest_dir)
    dest_p.mkdir(parents=True, exist_ok=True)
    name = dest_name or dest_name_for_taxid(taxid, src_p)
    target = dest_p / name
    try:
        src_resolved = src_p.resolve()
        dest_resolved = dest_p.resolve()
    except OSError:
        src_resolved = src_p
        dest_resolved = dest_p
    if src_resolved.parent == dest_resolved and src_p.name == name:
        return src_p if src_p.exists() else target
    if target.exists() or target.is_symlink():
        return target
    try:
        target.symlink_to(src_resolved)
        logger.info("Linked genome %s -> %s", target, src_resolved)
        return target
    except OSError as exc:
        logger.warning(
            "Could not symlink %s -> %s (%s); copying instead",
            src_resolved,
            target,
            exc,
        )
        shutil.copy2(src_resolved, target)
        return target


def _path_is_dir(path: Path) -> bool:
    try:
        return path.is_dir()
    except OSError:
        return False


def _path_is_file(path: Path) -> bool:
    try:
        return path.is_file()
    except OSError:
        return False


def _iter_taxid_files(directory: Path, taxid: str) -> List[Path]:
    hits: List[Path] = []
    seen = set()

    def add(path: Path) -> None:
        try:
            key = str(path.resolve())
        except OSError:
            key = str(path)
        if key in seen or not _path_is_file(path):
            return
        seen.add(key)
        hits.append(path)

    for ext in genome_lookup_extensions():
        add(directory / f"{taxid}{ext}")
    # Prefixed copies (e.g. Bacteria_562-processed.fasta from the realistic pipeline).
    for ext in genome_lookup_extensions():
        try:
            found = list(directory.glob(f"*_{taxid}{ext}"))
        except OSError:
            found = []
        for path in found:
            add(path)
    return hits


def genome_library_dirs(
    extra: Optional[Sequence[PathLike]] = None,
    *,
    include_test_genomes: bool = False,
    include_cache: bool = True,
) -> List[Path]:
    ordered: List[Path] = []
    seen = set()

    def add(raw: Optional[PathLike]) -> None:
        if raw is None:
            return
        path = _as_path(raw)
        if not str(path).strip():
            return
        try:
            exists = path.exists()
        except OSError:
            exists = False
        try:
            key = str(path.resolve()) if exists else str(path)
        except OSError:
            key = str(path)
        if key in seen:
            return
        if not include_test_genomes and is_bundled_test_genomes_path(path):
            return
        seen.add(key)
        ordered.append(path)

    for item in extra or []:
        add(item)
    env_extra = os.environ.get("SAMOVAR_GENOME_DIRS", "").strip()
    for item in _split_dir_value(env_extra):
        add(item)
    if include_cache:
        add(processed_genomes_dir())
        add(genome_download_dir())
    cfg = load_config()
    for item in _split_dir_value(cfg.get("genome_dirs")):
        add(item)
    return ordered


def find_library_genome(
    taxid: Union[str, int],
    extra: Optional[Sequence[PathLike]] = None,
    *,
    include_test_genomes: bool = False,
    processed_only: bool = False,
) -> Optional[Path]:
    taxid = _normalize_taxid(taxid)
    for directory in genome_library_dirs(
        extra, include_test_genomes=include_test_genomes
    ):
        if not _path_is_dir(directory):
            continue
        for candidate in _iter_taxid_files(directory, taxid):
            if processed_only and "-processed." not in candidate.name.lower():
                continue
            return candidate
    return None


def register_genome_dir(path: PathLike, *, force: bool = False) -> Optional[Path]:
    """Append ``path`` to config ``genome_dirs`` unless it is truncated test data."""
    directory = _as_path(path)
    try:
        directory = directory.resolve()
    except OSError:
        pass
    if not force and is_bundled_test_genomes_path(directory):
        logger.warning(
            "Refusing to register bundled test_genomes (truncated ISS stubs) as an "
            "NCBI library: %s",
            directory,
        )
        return None
    if not _path_is_dir(directory):
        logger.warning("Not registering missing genome directory: %s", directory)
        return None
    cfg = load_config()
    dirs = []
    seen = set()
    for item in _split_dir_value(cfg.get("genome_dirs")):
        try:
            key = str(Path(item).expanduser().resolve())
        except OSError:
            key = item
        if key in seen:
            continue
        seen.add(key)
        dirs.append(str(Path(item).expanduser()))
    key = str(directory)
    if key not in seen:
        dirs.append(key)
        update_config({"genome_dirs": dirs}, also_repo_build=False)
        logger.info("Registered genome library: %s", directory)
    return directory


def ensure_genome_config() -> Dict[str, Any]:
    """Write default ``genomes`` / ``processed_genomes`` paths into the install config."""
    genomes = str(genome_download_dir())
    processed = str(processed_genomes_dir())
    cfg = load_config()
    updates: Dict[str, Any] = {}
    if not str(cfg.get("genomes") or "").strip():
        updates["genomes"] = genomes
    if not str(cfg.get("processed_genomes") or "").strip():
        updates["processed_genomes"] = processed
    if "genome_dirs" not in cfg:
        updates["genome_dirs"] = list(_split_dir_value(cfg.get("genome_dirs")))
    if updates:
        cfg = update_config(updates, also_repo_build=False)
    return cfg


def _genome_dirs_from_generate_yaml(yaml_path: Path) -> List[Path]:
    if not _path_is_file(yaml_path):
        return []
    try:
        import yaml
    except ImportError:
        return []
    try:
        data = yaml.safe_load(yaml_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return []
    if not isinstance(data, dict):
        return []
    dirs: List[Path] = []
    seen = set()

    def add(raw: Optional[PathLike]) -> None:
        if not raw:
            return
        path = _as_path(raw)
        if _path_is_file(path):
            path = path.parent
        if not _path_is_dir(path):
            return
        try:
            key = str(path.resolve())
        except OSError:
            key = str(path)
        if key in seen:
            return
        seen.add(key)
        dirs.append(path)

    add(data.get("genome_dir"))
    add(data.get("host_genome"))
    return dirs


def generate_source_genome_dirs(output_dir: PathLike) -> List[Path]:
    """ISS/CAMISIM ``genome_dir`` / host parent from ``.generate/configs/*.yaml``.

    These are the assemblies used to *simulate* this run. They may be truncated
    test genomes; they are not registered as an NCBI library.
    """
    cfg_dir = _as_path(output_dir) / ".generate" / "configs"
    dirs: List[Path] = []
    seen = set()
    for name in ("iss_config.yaml", "camisim.yaml"):
        for path in _genome_dirs_from_generate_yaml(cfg_dir / name):
            try:
                key = str(path.resolve())
            except OSError:
                key = str(path)
            if key in seen:
                continue
            seen.add(key)
            dirs.append(path)
    return dirs


def seed_run_genomes(
    dest_dir: PathLike,
    *,
    reuse: Optional[bool] = None,
    extra_dirs: Optional[Sequence[PathLike]] = None,
    include_test_genomes: bool = False,
    generate_output_dir: Optional[PathLike] = None,
) -> Dict[str, Any]:
    """Populate ``$out_dir/genomes`` from generate sources and NCBI libraries."""
    dest = _as_path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    stats = {"linked": 0, "copied": 0, "skipped": 0, "sources": []}
    do_reuse = reuse_enabled(reuse)

    sources: List[Path] = []
    if generate_output_dir:
        sources.extend(generate_source_genome_dirs(generate_output_dir))
    if include_test_genomes or allow_test_genomes():
        logger.warning(
            "Seeding truncated bundled test_genomes into %s. These are ISS/CI "
            "stubs, not NCBI assemblies — do not use them for real metagenomes.",
            dest,
        )
        tg = test_genomes_dir()
        for sub in (tg, tg / "meta", tg / "host"):
            if _path_is_dir(sub):
                sources.append(sub)
    if do_reuse:
        # Explicit extra dirs (prepare --genome-dirs / check --src) are bulk-linked.
        # Config genome_dirs are *not* copied wholesale; fetch_genome looks them up
        # by taxid so a toy run does not inherit every cached assembly.
        for raw in extra_dirs or []:
            path = _as_path(raw)
            if _path_is_dir(path) and (
                include_test_genomes or not is_bundled_test_genomes_path(path)
            ):
                sources.append(path)

    seen_taxids = set()
    for src_dir in sources:
        if not _path_is_dir(src_dir):
            continue
        try:
            if src_dir.resolve() == dest.resolve():
                continue
        except OSError:
            pass
        stats["sources"].append(str(src_dir))
        try:
            children = sorted(src_dir.iterdir())
        except OSError:
            continue
        for path in children:
            try:
                is_file = path.is_file()
            except OSError:
                is_file = False
            if not is_file:
                continue
            if not is_fasta_name(path.name, protein=False, nucleotide=True):
                continue
            taxid = taxid_from_fasta_name(path)
            if taxid is None:
                # Prefixed processed files: Bacteria_562-processed.fasta
                name = path.name
                for ext in genome_lookup_extensions():
                    if name.endswith(ext):
                        stem = name[: -len(ext)]
                        if "_" in stem:
                            maybe = stem.rsplit("_", 1)[-1]
                            if maybe.isdigit():
                                taxid = maybe
                        break
            if taxid is None:
                continue
            key = f"{taxid}:{'-processed' if '-processed.' in path.name.lower() else 'raw'}"
            if key in seen_taxids:
                stats["skipped"] += 1
                continue
            before = dest / dest_name_for_taxid(taxid, path)
            existed = before.exists() or before.is_symlink()
            placed = place_genome(path, dest, taxid)
            seen_taxids.add(key)
            if existed:
                stats["skipped"] += 1
            elif placed.is_symlink():
                stats["linked"] += 1
            else:
                stats["copied"] += 1

    logger.info(
        "Seeded genomes into %s (linked=%s copied=%s skipped=%s)",
        dest,
        stats["linked"],
        stats["copied"],
        stats["skipped"],
    )
    return stats


def remember_prepare_genome_paths(
    extra_dirs: Optional[Sequence[PathLike]] = None,
) -> Dict[str, Any]:
    """Record cache locations (and extra libraries) in the install config."""
    cfg = ensure_genome_config()
    for directory in extra_dirs or []:
        register_genome_dir(directory)
    return cfg


def _cmd_seed(args: argparse.Namespace) -> int:
    extra = _split_dir_value(args.genome_dirs)
    stats = seed_run_genomes(
        args.dest,
        reuse=False if args.no_reuse else None,
        extra_dirs=extra,
        include_test_genomes=args.test_genomes,
        generate_output_dir=args.generate_dir,
    )
    print(
        f"seeded {args.dest}: linked={stats['linked']} "
        f"copied={stats['copied']} skipped={stats['skipped']}"
    )
    return 0


def _cmd_check(args: argparse.Namespace) -> int:
    src = _as_path(args.src)
    dest = _as_path(args.dest)
    extra = [src, *_split_dir_value(args.genome_dirs)]
    remember_prepare_genome_paths(extra_dirs=[src])
    stats = seed_run_genomes(
        dest,
        reuse=True,
        extra_dirs=extra,
        include_test_genomes=False,
        generate_output_dir=args.generate_dir,
    )
    missing: List[str] = []
    taxids = set()
    if src.is_dir():
        for path in src.iterdir():
            taxid = taxid_from_fasta_name(path)
            if taxid is None:
                continue
            taxids.add(taxid)
            found = _iter_taxid_files(dest, taxid)
            if not found:
                missing.append(taxid)
    print(
        f"reuse check: src={src} dest={dest} taxids={len(taxids)} "
        f"missing={len(missing)} linked={stats['linked']} copied={stats['copied']}"
    )
    if missing:
        print("missing taxids: " + ",".join(sorted(missing)[:20]))
        return 1
    n_links = sum(1 for p in dest.iterdir() if p.is_symlink()) if dest.is_dir() else 0
    print(f"symlinks in dest: {n_links}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="SamovaR genome cache helpers")
    sub = parser.add_subparsers(dest="cmd", required=True)

    seed = sub.add_parser("seed", help="Populate a run genomes/ directory")
    seed.add_argument("--dest", required=True)
    seed.add_argument("--generate-dir", default=None)
    seed.add_argument("--genome-dirs", default="")
    seed.add_argument("--test-genomes", action="store_true")
    seed.add_argument("--no-reuse", action="store_true", help="Skip config genome libraries")
    seed.set_defaults(func=_cmd_seed)

    check = sub.add_parser("check", help="Register src and seed dest; fail if taxids missing")
    check.add_argument("--src", required=True)
    check.add_argument("--dest", required=True)
    check.add_argument("--genome-dirs", default="")
    check.add_argument("--generate-dir", default=None)
    check.set_defaults(func=_cmd_check)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
