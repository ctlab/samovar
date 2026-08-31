"""Add a classifier to an existing prepare/exec run without rewriting samovar.sh.

``samovar prepare --add-annotator --kraken2-test 'kraken2 toy'`` keeps the original
Hydra prepare parameters, writes ``.log/samovar_v{N}.sh``, and drops checkpoints
from ``annotate_initial`` onward so exec supplements reports/tables. Annotator
Snakemake already skips existing ``*.out`` files. ISS skips per-annotator FASTQs
that already have reads.
"""

from __future__ import annotations

import re
from argparse import Namespace
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

from samovar.exec_control import CHECKPOINT_STEPS, as_path, clear_checkpoints
from samovar.paths import absolute_path

ACTIVE_PIPELINE = "active_pipeline"
PIPELINE_RE = re.compile(r"^samovar_v(\d+)\.sh$")
KEEP_CHECKPOINTS = {"setup_reads", "qc_initial", "seed_genomes"}
STALE_ABUNDANCE = (
    Path("initial_abundance"),
    Path("regenerated") / ".regenerated_abundance",
)
ISS_DONE = Path("regenerated") / ".process_annotations.done"


def log_dir(output_dir) -> Path:
    return as_path(output_dir) / ".log"


def listed_pipeline_scripts(output_dir) -> List[Path]:
    folder = log_dir(output_dir)
    if not folder.is_dir():
        return []
    found = []
    main = folder / "samovar.sh"
    if main.is_file():
        found.append(main)
    versions = []
    for path in folder.glob("samovar_v*.sh"):
        match = PIPELINE_RE.match(path.name)
        if match:
            versions.append((int(match.group(1)), path))
    versions.sort()
    found.extend(path for _, path in versions)
    return found


def next_pipeline_version(output_dir) -> int:
    highest = 0
    for path in listed_pipeline_scripts(output_dir):
        match = PIPELINE_RE.match(path.name)
        if match:
            highest = max(highest, int(match.group(1)))
    return highest + 1


def write_active_pipeline(output_dir, script_name: str) -> Path:
    dest = log_dir(output_dir) / ACTIVE_PIPELINE
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(Path(script_name).name + "\n", encoding="utf-8")
    return dest


def active_pipeline_script(output_dir) -> Path:
    """Script ``samovar exec`` should run: pointer, else latest ``samovar_vN.sh``, else ``samovar.sh``."""
    folder = log_dir(output_dir)
    pointer = folder / ACTIVE_PIPELINE
    if pointer.is_file():
        name = pointer.read_text(encoding="utf-8").strip()
        if name:
            candidate = folder / Path(name).name
            if candidate.is_file():
                return candidate
    scripts = listed_pipeline_scripts(output_dir)
    if not scripts:
        return folder / "samovar.sh"
    versions = [p for p in scripts if PIPELINE_RE.match(p.name)]
    if versions:
        return versions[-1]
    return scripts[0]


def _plain(value: Any) -> Any:
    if value is None:
        return None
    if hasattr(value, "items"):
        try:
            from omegaconf import OmegaConf

            return OmegaConf.to_container(value, resolve=True)
        except Exception:
            return dict(value)
    return value


def load_prepare_hydra_args(output_dir) -> Dict[str, Any]:
    from samovar.repro import hydra_dir_from_user

    try:
        hydra_root = hydra_dir_from_user(as_path(output_dir))
    except FileNotFoundError:
        return {}
    cfg_path = hydra_root / "config.yaml"
    if not cfg_path.is_file():
        return {}
    from omegaconf import OmegaConf

    snap = OmegaConf.load(cfg_path)
    block = _plain(snap.get("prepare")) or {}
    if not isinstance(block, dict):
        return {}
    args = block.get("args") or {}
    return dict(args) if isinstance(args, dict) else {}


def existing_run_names(output_dir) -> List[str]:
    names: List[str] = []
    init = log_dir(output_dir) / "configs" / "config_init.yaml"
    if init.is_file():
        import yaml

        data = yaml.safe_load(init.read_text(encoding="utf-8")) or {}
        for row in data.get("run_config") or []:
            name = str((row or {}).get("run_name") or "").strip()
            if name:
                names.append(name)
    hydra_args = load_prepare_hydra_args(output_dir)
    for key, value in hydra_args.items():
        if not str(key).startswith("cmd_") or not value:
            continue
        names.append(str(key)[4:])
    seen = set()
    out = []
    for name in names:
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def _as_cmd_cells(value: Any) -> Any:
    if value in (None, "", [], {}):
        return value
    if isinstance(value, str):
        return [[value]]
    if isinstance(value, list) and value and isinstance(value[0], str):
        return [[value[0]]]
    return value


def cli_annotator_keys(args: Any) -> Dict[str, Any]:
    data = vars(args) if not isinstance(args, dict) else args
    out: Dict[str, Any] = {}
    for key, value in data.items():
        if not str(key).startswith("cmd_"):
            continue
        if not value:
            continue
        out[str(key)] = _as_cmd_cells(value)
    for key in ("kraken2", "kaiju", "dummy"):
        value = data.get(key)
        if value:
            out[key] = _as_cmd_cells(value)
    return out


def merge_add_annotator_namespace(args: Any) -> Tuple[Any, List[str]]:
    """Hydra prepare args + new CLI annotators. Returns (namespace, new run names)."""
    out_dir = absolute_path(getattr(args, "output_dir", None) or ".")
    hydra_args = load_prepare_hydra_args(out_dir)
    if not hydra_args and not (log_dir(out_dir) / "samovar.sh").is_file():
        raise ValueError(
            f"--add-annotator needs an existing run with .hydra and .log/samovar.sh under {out_dir}"
        )
    merged = dict(hydra_args)
    merged["output_dir"] = out_dir
    merged.pop("add_annotator", None)
    for key, value in list(merged.items()):
        if str(key).startswith("cmd_") or key in {"kraken2", "kaiju", "dummy"}:
            merged[key] = _as_cmd_cells(value)
    new_keys = cli_annotator_keys(args)
    if not new_keys:
        raise ValueError(
            "--add-annotator requires a new classifier flag, "
            "e.g. --kraken2-test 'kraken2 toy'"
        )
    already = set(existing_run_names(out_dir))
    added: List[str] = []
    for key, value in new_keys.items():
        run_name = key[4:] if key.startswith("cmd_") else key
        if run_name in already or key in hydra_args:
            raise ValueError(
                f"--add-annotator: {run_name} is already in this run "
                f"({', '.join(sorted(already)) or 'see config_init.yaml'})"
            )
        merged[key] = value
        added.append(run_name)
        already.add(run_name)
    ns = Namespace(**merged)
    ns.add_annotator = True
    ns.output_dir = out_dir
    if not getattr(ns, "input_dir", None) and not getattr(ns, "input_config", None):
        initial = as_path(out_dir) / "initial"
        if initial.is_dir():
            ns.input_dir = str(initial)
    return ns, added


def invalidate_add_annotator(output_dir) -> List[str]:
    """Drop checkpoints/tables that must be rebuilt after a new classifier."""
    root = as_path(output_dir)
    removed: List[str] = []
    for name in CHECKPOINT_STEPS:
        if name in KEEP_CHECKPOINTS:
            continue
        for path in clear_checkpoints(root, name):
            removed.append(str(path))
    done = root / ISS_DONE
    if done.is_file():
        done.unlink()
        removed.append(str(done))
    for rel in STALE_ABUNDANCE:
        folder = root / rel
        if not folder.is_dir():
            continue
        for path in folder.glob("*.csv"):
            path.unlink()
            removed.append(str(path))
    return removed


def annotator_fastqs_complete(
    output_dir: str,
    sample_names: Sequence[str],
    annotator_name: str,
    gzip_reads: bool = False,
) -> bool:
    """True when every sample already has a non-empty R1 for this annotator."""
    from samovar.seqio import fastq_has_reads, fastq_pair_paths

    if not sample_names:
        return False
    for sample in sample_names:
        r1, _r2 = fastq_pair_paths(
            Path(output_dir) / f"{sample}_{annotator_name}",
            gzip_reads=gzip_reads,
        )
        alt, _ = fastq_pair_paths(
            Path(output_dir) / f"{sample}_{annotator_name}",
            gzip_reads=not gzip_reads,
        )
        if not (fastq_has_reads(r1) or fastq_has_reads(alt)):
            return False
    return True
