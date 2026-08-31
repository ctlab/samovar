"""Read simulators: builtins (ISS, ART, wgsim, NanoSim/CAMISIM) plus imported tools.

Analog of ``table_regenerators``. A generator takes either an abundance table
or generate-style genome dirs + host/community fractions and writes FASTQ
(``.fastq`` / ``.fastq.gz``) whose headers keep ``taxid:`` (true label) and,
for hybrid runs, ``read_type:``. Filenames are ``{sample}_full_R*`` on generate
and ``{sample}_{annotator}_R*`` on prepare; extra ids (``sequence_type``,
``illumina``, ``ont``, …) may appear between sample and mate.
"""

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

from samovar.main_config import (
    flags_target_matches,
    imported_flags_for_names,
    iter_tools,
    merge_flag_strings,
    parse_tool_entry,
    tool_flags,
    tool_path,
)
from samovar.paths import load_config
from samovar.seqio import concat_fastq_files, fastq_pair_paths, is_fastq_name
from samovar.table_regenerators import extra_flags_argv

READS_GENERATOR_GROUPS = ("reads_generator", "metagenome_generator")

# Canonical builtin names → generate backend.
_BUILTIN_ALIASES = {
    "iss": "iss",
    "insilicoseq": "iss",
    "insilico": "iss",
    "art": "art",
    "art_illumina": "art",
    "illumina": "art",
    "hiseq": "art",
    "miseq": "art",
    "novaseq": "art",
    "wgsim": "wgsim",
    "nanosim": "nanosim",
    "nanosim3": "nanosim",
    "ont": "nanosim",
    "nanopore": "nanosim",
    "simulator.py": "nanosim",
    "camisim": "camisim",
    "hybrid": "hybrid",
}

_CAMISIM_MODE = {
    "art": "illumina",
    "wgsim": "wgsim",
    "nanosim": "ont",
    "hybrid": "hybrid",
    "camisim": "table",
}

_READ_TYPE_FOR = {
    "art": "illumina",
    "wgsim": "wgsim",
    "nanosim": "ont",
    "iss": "",
    "camisim": "",
    "hybrid": "hybrid",
}

# Tokens allowed between sample id and ``_R1`` / ``_full``.
FILENAME_EXTRA_IDS = (
    "illumina",
    "bgi",
    "mgi",
    "ont",
    "wgsim",
    "nanosim3",
    "art",
    "sequence_type",
    "read_type",
)


logger = logging.getLogger(__name__)


class MissingReadsGeneratorError(ValueError):
    """``tools.<name>`` is missing or is not a reads/metagenome generator."""


def resolve_reads_generator(name: Optional[str]) -> Tuple[str, str]:
    """Return ``("builtin"|"custom", canonical_name)``."""
    key = str(name or "iss").strip() or "iss"
    low = key.lower().replace("-", "_")
    if low in _BUILTIN_ALIASES:
        return "builtin", _BUILTIN_ALIASES[low]
    return "custom", key


def require_known_reads_generator(name: Optional[str]) -> str:
    kind, canon = resolve_reads_generator(name)
    if kind == "custom":
        lookup_reads_generator(canon)
    return canon


def flags_apply_to_reads_generator(target: str, name: Optional[str]) -> bool:
    kind, canon = resolve_reads_generator(name)
    names = [canon, name, _READ_TYPE_FOR.get(canon, "")]
    if kind == "builtin":
        if canon == "art":
            names.extend(["art_illumina", "illumina"])
        if canon == "nanosim":
            names.extend(["nanosim3", "ont", "simulator.py"])
        if canon == "hybrid":
            names.extend(["camisim"])
    return flags_target_matches(
        target,
        *names,
        groups=("reads_generator", "reads", "read", "metagenome_generator", "meta"),
    )


def extra_ids_for_generator(name: Optional[str], config: Optional[Dict[str, Any]] = None) -> List[str]:
    cfg = dict(config or {})
    raw = cfg.get("filename_ids") or cfg.get("extra_ids") or cfg.get("read_types")
    if raw:
        if isinstance(raw, str):
            parts = [p.strip() for p in raw.replace(",", " ").split() if p.strip()]
        else:
            parts = [str(p).strip() for p in raw if str(p).strip()]
        return [p for p in parts if p and p.lower() != "full"]
    kind, canon = resolve_reads_generator(name)
    if kind != "builtin":
        return []
    if canon == "art":
        return ["illumina"]
    if canon == "wgsim":
        return ["wgsim"]
    if canon == "nanosim":
        return ["ont"]
    if canon == "hybrid":
        return ["illumina", "ont"]
    return []


def reads_output_stems(
    sample: str,
    *,
    annotator: Optional[str] = None,
    extra_ids: Sequence[str] = (),
    stage: str = "generate",
) -> List[str]:
    """FASTQ prefixes (no ``_R1``). Primary stem is last (``_full`` or ``_{annotator}``)."""
    sample = str(sample)
    extras = [str(x).strip() for x in extra_ids if str(x).strip() and str(x).strip().lower() != "full"]
    if str(stage).lower() == "regenerate":
        base = f"{sample}_{annotator}" if annotator else sample
        return [f"{base}_{eid}" for eid in extras] + [base]
    return [f"{sample}_{eid}" for eid in extras] + [f"{sample}_full"]


def attach_reads_flags(name: Optional[str], config: Dict[str, Any]) -> Dict[str, Any]:
    """Merge ``tools.*[4]`` with YAML / prepare ``--flags``."""
    cfg = dict(config or {})
    kind, canon = resolve_reads_generator(name or cfg.get("reads_generator"))
    tools = iter_tools(load_config())
    lookup_names = [canon, name, cfg.get("reads_generator")]
    if kind == "builtin":
        if canon == "art":
            lookup_names.extend(["art_illumina", "illumina"])
        if canon == "nanosim":
            lookup_names.extend(["nanosim3", "ont", "simulator.py"])
        if canon == "iss":
            lookup_names.append("iss")
    imported = imported_flags_for_names(tools, *[str(n) for n in lookup_names if n])
    if kind == "custom":
        try:
            spec = lookup_reads_generator(canon)
            imported = merge_flag_strings(tool_flags(spec, canon), imported)
        except MissingReadsGeneratorError:
            pass
    cfg["extra_flags"] = merge_flag_strings(
        imported,
        cfg.get("extra_flags"),
        cfg.get("reads_generator_flags"),
        cfg.get("metagenome_generator_flags"),
    )
    cfg["extra_argv"] = extra_flags_argv(cfg.get("extra_flags"))
    cfg["reads_generator"] = canon
    return cfg


def lookup_reads_generator(name: str) -> list:
    key = str(name or "").strip()
    if not key:
        raise MissingReadsGeneratorError(
            "Empty reads_generator name. Import a tool with "
            "`samovar tools import -n NAME --type reads`."
        )
    tools = iter_tools(load_config())
    spec = tools.get(key)
    matched = key
    if spec is None:
        low = key.lower()
        for stored, value in tools.items():
            if stored.lower() == low:
                spec = value
                matched = stored
                break
    if spec is None:
        raise MissingReadsGeneratorError(
            f"reads_generator {key!r} is not in the main install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type reads` "
            "before generate / prepare."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group not in READS_GENERATOR_GROUPS:
        raise MissingReadsGeneratorError(
            f"tools.{matched} has group {group!r}, expected one of "
            f"{READS_GENERATOR_GROUPS}. Re-import with --type reads or --type meta."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingReadsGeneratorError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return parsed


def camisim_mode_for_reads_generator(name: Optional[str], fallback: Optional[str] = None) -> Optional[str]:
    kind, canon = resolve_reads_generator(name)
    if kind != "builtin":
        return None
    if canon == "iss":
        return None
    if canon == "camisim":
        return str(fallback or "table")
    return _CAMISIM_MODE.get(canon)


class ReadsGenerator(ABC):
    """One abundance table or generate-style dirs in; FASTQ paths out."""

    @abstractmethod
    def generate(self, metadata: Optional[pd.DataFrame], config: Dict[str, Any]) -> List[str]:
        raise NotImplementedError


class IssReadsGenerator(ReadsGenerator):
    def generate(self, metadata, config):
        from samovar.table2iss import generate_iss_test_samples

        cfg = attach_reads_flags("iss", dict(config or {}))
        _ = metadata
        return generate_iss_test_samples(
            genome_dir=str(cfg.get("genome_dir") or ""),
            host_genome=str(cfg.get("host_genome") or ""),
            output_dir=str(cfg["output_dir"]),
            n_samples=int(cfg.get("n_samples") or 10),
            total_reads=int(cfg.get("total_reads") or cfg.get("N_reads") or 2000),
            host_fraction=cfg.get("host_fraction", "RANDOM"),
            seed=int(cfg.get("seed") or 42),
            model=str(cfg.get("model") or "hiseq"),
            genomes=cfg.get("genomes"),
            cpus=int(cfg.get("cores") or cfg.get("cpus") or 2),
            extra_flags=cfg.get("extra_argv") or extra_flags_argv(cfg.get("extra_flags")),
            abundance_table=cfg.get("abundance_table") or cfg.get("abundance"),
            max_genomes=cfg.get("max_genomes"),
        )


class CamisimReadsGenerator(ReadsGenerator):
    """ART / wgsim / NanoSim / hybrid / generic CAMISIM generate backend."""

    def __init__(self, name: str):
        self.name = name

    def generate(self, metadata, config):
        from samovar.camisim import run_from_config, setup_camisim_generate

        _ = metadata
        cfg = attach_reads_flags(self.name, dict(config or {}))
        mode = camisim_mode_for_reads_generator(self.name, cfg.get("camisim_mode"))
        args = type(
            "Args",
            (),
            {
                "genome_dir": cfg.get("genome_dir"),
                "output_dir": str(Path(cfg["output_dir"]).parent)
                if Path(cfg["output_dir"]).name == "initial"
                else cfg.get("output_root") or str(Path(cfg["output_dir"]).parent),
                "host_genome": cfg.get("host_genome") or "",
                "n_samples": int(cfg.get("n_samples") or 10),
                "total_reads": int(cfg.get("total_reads") or 2000),
                "host_fraction": cfg.get("host_fraction", "RANDOM"),
                "seed": int(cfg.get("seed") or 42),
                "model": str(cfg.get("model") or "hiseq"),
                "cores": int(cfg.get("cores") or 1),
                "camisim_mode": mode,
                "camisim_config": cfg.get("camisim_config"),
                "size_gbp": cfg.get("size_gbp"),
                "max_genomes": cfg.get("max_genomes"),
            },
        )()
        # Prefer writing YAML then running when this is invoked in-process.
        paths = setup_camisim_generate(args)
        yaml_path = paths.get("config") or paths.get("yaml")
        if yaml_path:
            data = Path(yaml_path).read_text(encoding="utf-8")
            extra = cfg.get("extra_flags") or ""
            if extra:
                import yaml

                mapped = yaml.safe_load(data) or {}
                mapped["extra_flags"] = extra
                mapped["reads_generator_flags"] = extra
                if cfg.get("max_genomes") is not None:
                    mapped["max_genomes"] = cfg.get("max_genomes")
                Path(yaml_path).write_text(yaml.dump(mapped), encoding="utf-8")
        if cfg.get("run", True) and yaml_path:
            result = run_from_config(yaml_path)
            return [str(p) for p in (result.get("reads") or [])]
        return [paths.get("pipeline") or ""]


class CustomReadsGenerator(ReadsGenerator):
    def __init__(self, name: str):
        self.name = str(name).strip()

    def generate(self, metadata, config):
        spec = lookup_reads_generator(self.name)
        path = Path(tool_path(spec, self.name)).expanduser()
        cfg = attach_reads_flags(self.name, dict(config or {}))
        py_fn = _try_python_callable(path, self.name)
        if py_fn is not None:
            out = py_fn(_job_spec(cfg), metadata, cfg)
            if not isinstance(out, (list, tuple)):
                raise TypeError(
                    f"{path} generate() must return list[str] of FASTQ paths, "
                    f"got {type(out).__name__}"
                )
            return [str(p) for p in out]
        return _run_cli_generator(path, metadata, cfg)


def _job_spec(cfg: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "abundance_table": cfg.get("abundance_table") or cfg.get("abundance"),
        "genome_dir": cfg.get("genome_dir"),
        "host_genome": cfg.get("host_genome"),
        "host_fraction": cfg.get("host_fraction"),
        "n_samples": cfg.get("n_samples"),
        "total_reads": cfg.get("total_reads") or cfg.get("N_reads"),
        "output_dir": cfg.get("output_dir"),
        "stage": cfg.get("stage") or "generate",
        "annotator": cfg.get("annotator_name") or cfg.get("annotator"),
        "extra_ids": extra_ids_for_generator(cfg.get("reads_generator"), cfg),
        "gzip_reads": bool(cfg.get("gzip_reads")),
        "max_genomes": cfg.get("max_genomes"),
    }


def _try_python_callable(path: Path, name: str) -> Optional[Callable]:
    if path.suffix.lower() != ".py" or not path.is_file():
        return None
    try:
        spec = importlib.util.spec_from_file_location(
            f"samovar_custom_reads_{name}", path
        )
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception:
        return None
    fn = getattr(module, "generate", None)
    if callable(fn):
        return fn
    cls = getattr(module, "ReadsGenerator", None)
    if cls is None:
        return None
    try:
        inst = cls()
    except Exception:
        return None
    run = getattr(inst, "generate", None) or getattr(inst, "run", None)
    return run if callable(run) else None


def _run_cli_generator(
    path: Path,
    metadata: Optional[pd.DataFrame],
    config: Dict[str, Any],
) -> List[str]:
    output_dir = config.get("output_dir")
    if not output_dir:
        raise ValueError(
            f"CLI reads_generator {path} needs output_dir in config "
            "(or define generate() in a Python module)."
        )
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cmd: List[str] = [str(path), "-o", str(out)]
    abundance = config.get("abundance_table") or config.get("abundance")
    if abundance:
        cmd.extend(["-i", str(abundance)])
    if config.get("genome_dir"):
        cmd.extend(["--genome-dir", str(config["genome_dir"])])
    if config.get("host_genome"):
        cmd.extend(["--host-genome", str(config["host_genome"])])
    if config.get("host_fraction") not in (None, ""):
        cmd.extend(["--host-fraction", str(config["host_fraction"])])
    if config.get("n_samples") not in (None, ""):
        cmd.extend(["--n-samples", str(int(config["n_samples"]))])
    total = config.get("total_reads") or config.get("N_reads")
    if total not in (None, ""):
        cmd.extend(["--total-reads", str(int(total))])
    if config.get("seed") is not None:
        cmd.extend(["--seed", str(config["seed"])])
    from samovar.regenerate import finite_max_genomes

    limit = finite_max_genomes(config.get("max_genomes"), default_from_env=False)
    if limit is not None:
        cmd.extend(["--max-genomes", str(limit)])
    if config.get("model"):
        cmd.extend(["--model", str(config["model"])])
    if config.get("gzip_reads"):
        cmd.append("--gzip-reads")
    stage = str(config.get("stage") or "generate")
    cmd.extend(["--stage", stage])
    annot = config.get("annotator_name") or config.get("annotator")
    if annot:
        cmd.extend(["--annotator", str(annot)])
    for eid in extra_ids_for_generator(config.get("reads_generator"), config):
        cmd.extend(["--read-type", eid])
    meta_tmp = None
    try:
        if metadata is not None:
            import tempfile

            handle = tempfile.NamedTemporaryFile(
                mode="w", suffix=".csv", delete=False, encoding="utf-8"
            )
            metadata.to_csv(handle.name, index=False)
            handle.close()
            meta_tmp = handle.name
            cmd.extend(["-m", meta_tmp])
        extra = list(config.get("extra_argv") or extra_flags_argv(config.get("extra_flags")))
        cmd.extend(extra)
        if path.suffix.lower() == ".py":
            cmd = [sys.executable, *cmd]
        subprocess.run(cmd, check=True)
    finally:
        if meta_tmp:
            Path(meta_tmp).unlink(missing_ok=True)
    return _collect_fastq(out)


def _collect_fastq(output_dir: Path) -> List[str]:
    found: List[str] = []
    if not output_dir.is_dir():
        return found
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and is_fastq_name(path.name):
            found.append(str(path))
    if not found:
        raise ValueError(f"Custom reads_generator wrote no FASTQ under {output_dir}")
    return found


def get_reads_generator(name: Optional[str]) -> ReadsGenerator:
    kind, canon = resolve_reads_generator(name)
    if kind == "builtin":
        if canon == "iss":
            return IssReadsGenerator()
        return CamisimReadsGenerator(canon)
    lookup_reads_generator(canon)
    return CustomReadsGenerator(canon)


def alias_fastq_with_extra_ids(
    primary_r1: Union[str, Path],
    primary_r2: Union[str, Path],
    extra_ids: Sequence[str],
    *,
    gzip_reads: bool = False,
    read_type: str = "",
) -> List[str]:
    """Copy ``{stem}_R*`` to ``{stem}_{id}_R*`` and optionally tag ``read_type:``."""
    from samovar.camisim import tag_fastq_file

    written = [str(primary_r1), str(primary_r2)]
    r1 = Path(primary_r1)
    r2 = Path(primary_r2)
    if read_type:
        if r1.is_file():
            tag_fastq_file(r1, r1, read_type, "")
        if r2.is_file():
            tag_fastq_file(r2, r2, read_type, "")
    stem = str(r1)
    for suffix in ("_R1.fastq.gz", "_R1.fastq"):
        if stem.endswith(suffix):
            prefix = stem[: -len(suffix)]
            break
    else:
        prefix = str(r1.with_name(r1.name.split("_R1")[0]))
    for eid in extra_ids:
        dest_r1, dest_r2 = fastq_pair_paths(f"{prefix}_{eid}", gzip_reads=gzip_reads)
        if Path(primary_r1).is_file() and dest_r1 != str(primary_r1):
            shutil.copy2(primary_r1, dest_r1)
            written.append(dest_r1)
        if Path(primary_r2).is_file() and dest_r2 != str(primary_r2):
            shutil.copy2(primary_r2, dest_r2)
            written.append(dest_r2)
        tag = eid if not read_type else read_type
        if tag and Path(dest_r1).is_file():
            tag_fastq_file(dest_r1, dest_r1, tag, "")
        if tag and Path(dest_r2).is_file():
            tag_fastq_file(dest_r2, dest_r2, tag, "")
    return written


def simulate_from_sample_tables(
    sample_tables: Dict[str, pd.DataFrame],
    available_genomes: Dict[str, str],
    output_dir: str,
    *,
    annotator_cols: Dict[str, str],
    total_amount: Optional[int] = None,
    model: str = "hiseq",
    read_length: int = 150,
    cpus: int = 2,
    seed: int = 42,
    gzip_reads: bool = False,
    reads_generator: Optional[str] = None,
    extra_flags: Optional[Sequence[str]] = None,
    extra_ids: Optional[Sequence[str]] = None,
    metadata: Optional[pd.DataFrame] = None,
    config: Optional[Dict[str, Any]] = None,
    genome_dir: str = "",
) -> None:
    """ISS (or custom) mix from per-sample abundance tables (prepare/regenerate)."""
    from samovar.table2iss import (
        _annotator_from_n_col,
        _emit_empty_for_annotators,
        _n_columns,
        generate_reads_metagenome,
        iss_cli_extra_flags,
        split_metagenome_to_samples,
    )

    cfg = attach_reads_flags(reads_generator, dict(config or {}))
    sample_names = list(sample_tables)
    annotators = list(annotator_cols.keys()) or ["any"]
    meta_name = cfg.get("metagenome_generator")
    if meta_name:
        from samovar.metagenome_generators import (
            attach_metagenome_flags,
            get_metagenome_generator,
            resolve_metagenome_generator,
        )

        cfg = attach_metagenome_flags(meta_name, cfg)
        mkind, mcanon = resolve_metagenome_generator(meta_name)
        flags = list(
            extra_flags or cfg.get("extra_argv") or extra_flags_argv(cfg.get("extra_flags"))
        )
        if mkind == "custom":
            import tempfile

            gen = get_metagenome_generator(mcanon)
            for annotator_name in annotators:
                with tempfile.TemporaryDirectory() as tmp:
                    csv_path = Path(tmp) / f"{annotator_name}.csv"
                    _write_annotator_abundance_csv(
                        sample_tables, annotator_name, csv_path
                    )
                    job = dict(cfg)
                    job.update(
                        {
                            "abundance_table": str(csv_path),
                            "output_dir": output_dir,
                            "genome_dir": genome_dir or cfg.get("genome_dir"),
                            "stage": "regenerate",
                            "annotator_name": annotator_name,
                            "annotator": annotator_name,
                            "gzip_reads": gzip_reads,
                            "extra_argv": flags,
                            "seed": seed,
                            "model": model,
                        }
                    )
                    gen.generate(metadata, job)
            return
        reads_generator = mcanon

    kind, canon = resolve_reads_generator(reads_generator or cfg.get("reads_generator"))
    flags = list(extra_flags or cfg.get("extra_argv") or extra_flags_argv(cfg.get("extra_flags")))
    ids = list(extra_ids if extra_ids is not None else extra_ids_for_generator(canon, cfg))
    if kind == "builtin" and canon == "hybrid":
        ids = []
    os_makedirs = __import__("os").makedirs
    os_path = __import__("os").path
    rmtree = shutil.rmtree

    if kind == "custom":
        import tempfile

        gen = CustomReadsGenerator(canon)
        for annotator_name in annotators:
            with tempfile.TemporaryDirectory() as tmp:
                csv_path = Path(tmp) / f"{annotator_name}.csv"
                _write_annotator_abundance_csv(
                    sample_tables, annotator_name, csv_path
                )
                job = dict(cfg)
                job.update(
                    {
                        "abundance_table": str(csv_path),
                        "output_dir": output_dir,
                        "genome_dir": genome_dir or cfg.get("genome_dir"),
                        "stage": "regenerate",
                        "annotator_name": annotator_name,
                        "gzip_reads": gzip_reads,
                        "extra_argv": flags,
                        "seed": seed,
                        "model": model,
                    }
                )
                gen.generate(metadata, job)
        return

    os_makedirs(output_dir, exist_ok=True)
    pool_dir = os_path.join(output_dir, ".iss_full")
    os_makedirs(pool_dir, exist_ok=True)
    read_type = _READ_TYPE_FOR.get(canon, "") if kind == "builtin" else ""
    # ISS keeps taxid: in combined-genome headers; ART/wgsim/NanoSim prepare
    # uses the same mixer and then stamps read_type: + extra filename ids.
    with iss_cli_extra_flags(flags):
        for annotator_name in annotators:
            sample_source_counts: Dict[str, Dict[str, int]] = {}
            total_by_taxid: Dict[str, int] = {}
            for sample_name, table in sample_tables.items():
                counts: Dict[str, int] = {}
                n_col = next(
                    (
                        col
                        for col in _n_columns(table)
                        if _annotator_from_n_col(col) == annotator_name
                    ),
                    None,
                )
                if n_col is not None and not table.empty:
                    from samovar.table2iss import _scale_amounts

                    amounts = _scale_amounts(table[n_col].tolist(), total_amount)
                    taxids = table["taxid"].tolist()
                    for taxid, n_reads in zip(taxids, amounts):
                        taxid = str(taxid)
                        if taxid not in available_genomes or n_reads <= 0:
                            continue
                        counts[taxid] = counts.get(taxid, 0) + int(n_reads)
                        total_by_taxid[taxid] = total_by_taxid.get(taxid, 0) + int(
                            n_reads
                        )
                sample_source_counts[sample_name] = counts

            genome_ids = [taxid for taxid, n_reads in total_by_taxid.items() if n_reads > 0]
            if not genome_ids:
                _emit_empty_for_annotators(
                    output_dir, sample_names, [annotator_name], gzip_reads=gzip_reads
                )
                continue

            from samovar.add_annotator import annotator_fastqs_complete

            if annotator_fastqs_complete(
                output_dir, sample_names, annotator_name, gzip_reads=gzip_reads
            ):
                logger.info(
                    "Skipping ISS for %s: regenerated FASTQs already present",
                    annotator_name,
                )
                continue

            r1_full, r2_full = generate_reads_metagenome(
                genome_files=[available_genomes[taxid] for taxid in genome_ids],
                output_dir=pool_dir,
                amount=[total_by_taxid[taxid] for taxid in genome_ids],
                read_length=read_length,
                sample_name="full",
                annotator_name=annotator_name,
                model=model,
                genome_ids=genome_ids,
                cpus=cpus,
                seed=seed + annotators.index(annotator_name),
                gzip_reads=gzip_reads,
            )
            split_metagenome_to_samples(
                r1_full,
                r2_full,
                sample_source_counts,
                output_dir,
                annotator_name,
                gzip_reads=gzip_reads,
            )
            if ids or read_type:
                for sample_name in sample_names:
                    p1, p2 = fastq_pair_paths(
                        os_path.join(output_dir, f"{sample_name}_{annotator_name}"),
                        gzip_reads=gzip_reads,
                    )
                    alias_fastq_with_extra_ids(
                        p1,
                        p2,
                        ids,
                        gzip_reads=gzip_reads,
                        read_type=read_type if read_type != "hybrid" else "",
                    )
    rmtree(pool_dir, ignore_errors=True)


def _write_annotator_abundance_csv(
    sample_tables: Dict[str, pd.DataFrame],
    annotator_name: str,
    dest: Path,
) -> None:
    from samovar.table2iss import _annotator_from_n_col, _n_columns

    frames = []
    for sample, table in sample_tables.items():
        if table is None or table.empty or "taxid" not in table.columns:
            continue
        n_col = next(
            (
                col
                for col in _n_columns(table)
                if _annotator_from_n_col(col) == annotator_name
            ),
            None,
        )
        if n_col is None:
            continue
        piece = table[["taxid", n_col]].copy()
        piece = piece.rename(columns={n_col: f"N_{sample}"})
        frames.append(piece)
    if not frames:
        pd.DataFrame(columns=["taxid"]).to_csv(dest, index=False)
        return
    out = frames[0]
    for piece in frames[1:]:
        out = out.merge(piece, on="taxid", how="outer")
    out = out.fillna(0)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, index=False)


def run_generate_from_yaml(config_path: str) -> List[str]:
    import yaml

    cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8")) or {}
    meta_name = cfg.get("metagenome_generator")
    if meta_name:
        from samovar.metagenome_generators import (
            attach_metagenome_flags,
            get_metagenome_generator,
        )

        cfg = attach_metagenome_flags(meta_name, cfg)
        gen = get_metagenome_generator(meta_name)
    else:
        name = cfg.get("reads_generator") or "iss"
        gen = get_reads_generator(name)
    meta = None
    raw = cfg.get("samples_metadata") or cfg.get("metadata")
    if raw:
        meta = pd.read_csv(raw)
    cfg.setdefault("stage", "generate")
    return gen.generate(meta, cfg)


def write_custom_generate_script(base_dir: str, config_path: str, cores: int = 1) -> str:
    from samovar.paths import (
        python_path,
        repo_root,
        runtime_path_prefix,
        shell_source_install_env_snippet,
        shell_outdir_override_snippet,
    )

    base = Path(base_dir)
    generate_dir = base / ".generate"
    generate_dir.mkdir(parents=True, exist_ok=True)
    pipeline_path = generate_dir / "generate.sh"
    root = repo_root()
    py = python_path()
    tool_path_s = runtime_path_prefix()
    env_snippet = shell_source_install_env_snippet()
    pipeline_content = f"""# Setup
set -e
export SAMOVAR_ROOT="{root}"
{env_snippet}export PATH="${{SAMOVAR_PATH:+$SAMOVAR_PATH:}}{tool_path_s}:{root}/bin:$PATH"
export PYTHONPATH="{root / 'src'}${{PYTHONPATH:+:$PYTHONPATH}}"
PYTHON_PATH="${{PYTHON_PATH:-{py}}}"
if [ -z "$PYTHON_PATH" ] || [ ! -x "$PYTHON_PATH" ]; then
  PYTHON_PATH="$(command -v python3 || command -v python || true)"
fi
PYTHON_PATH="${{PYTHON_PATH:-python3}}"

out_dir="{base}"
{shell_outdir_override_snippet()}
mkdir -p "$out_dir"

"$PYTHON_PATH" -c "from samovar.reads_generators import run_generate_from_yaml; run_generate_from_yaml(r'{config_path}')"
"""
    pipeline_path.write_text(pipeline_content, encoding="utf-8")
    _ = cores
    return str(pipeline_path)
