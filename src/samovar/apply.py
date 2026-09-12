"""Apply a previously configured SamovaR pipeline to a new sample.

``samovar apply`` reuses the source run's annotators, QC, export, scoring, and
trained ML reprofiler without fitting. ``--full`` additionally regenerates
abundance tables (the in-run analogue of ModDirect) and refits the configured
reprofiler with ``samovar.reprofilers.run_reprofiler``.
"""

from __future__ import annotations

import argparse
import copy
import os
import random
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import yaml

from samovar.abundance import materialize_observed_abundance, regenerated_abundance_dir
from samovar.abundance_correctors import export_run_dir, is_skipped_export
from samovar.combine_tables import combine_with_cpp
from samovar.exec_control import CHECKPOINT_STEPS, mark_done, resolve_window
from samovar.paths import add_output_dir_argument, repo_root, runtime_path_prefix
from samovar.qc import trim_stage
from samovar.regenerate import stage_regenerate_tables
from samovar.reprofilers import apply_saved_reprofiler, run_reprofiler
from samovar.scorers import run_custom_scorers
from samovar.seqio import has_r1_reads, link_or_copy_reads, list_r1_files
from samovar.table_scorers import stage_score_regenerated_tables
from samovar.viz_annotation import compare_annotations

PathLike = Union[str, os.PathLike]

ISS_STEPS = frozenset(
    {
        "seed_genomes",
        "regenerate_reads",
        "sort_reads",
        "qc_generated",
        "annotate_regenerated",
        "combine_regenerated",
        "viz_regenerated",
    }
)

TABLE_ADAPT_STEPS = (
    "abundance_tables",
    "regenerate_tables",
    "score_regenerated_tables",
)


class ApplyError(ValueError):
    """Invalid apply inputs or pipeline state."""


class MissingPipelineError(ApplyError):
    """Source exec directory is missing required configs."""


class MissingModelError(ApplyError):
    """Normal apply needs a trained ``trained_model.joblib``."""


class IncompatibleInputError(ApplyError):
    """FASTQ input is missing or unusable."""


@dataclass
class PipelineState:
    """Persisted state from a completed (or prepared+executed) SamovaR run."""

    root: Path
    start: str
    end: str
    configs: Dict[str, Path] = field(default_factory=dict)
    init_config: Dict[str, Any] = field(default_factory=dict)
    qc_config: Dict[str, Any] = field(default_factory=dict)
    scoring_config: Dict[str, Any] = field(default_factory=dict)
    export_config: Dict[str, Any] = field(default_factory=dict)
    reprofiling_config: Dict[str, Any] = field(default_factory=dict)
    annotation2iss_config: Dict[str, Any] = field(default_factory=dict)
    model_path: Optional[Path] = None
    regenerated_annotations: Optional[Path] = None
    seed: int = 42


def _load_yaml(path: Optional[PathLike]) -> Dict[str, Any]:
    if not path:
        return {}
    file = Path(path)
    if not file.is_file():
        return {}
    data = yaml.safe_load(file.read_text(encoding="utf-8")) or {}
    return data if isinstance(data, dict) else {}


def _dump_yaml(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def _parse_window_env(path: Path) -> Tuple[str, str]:
    start = end = ""
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines():
            text = line.strip().replace("export ", "", 1)
            if text.startswith("SAMOVAR_START="):
                start = text.split("=", 1)[1].strip().strip('"').strip("'")
            elif text.startswith("SAMOVAR_END="):
                end = text.split("=", 1)[1].strip().strip('"').strip("'")
    return resolve_window(start or None, end or None)


def load_pipeline_state(pipeline_dir: PathLike) -> PipelineState:
    root = Path(pipeline_dir).expanduser().resolve()
    if not root.is_dir():
        raise MissingPipelineError(f"pipeline directory not found: {root}")
    cfg_dir = root / ".log" / "configs"
    configs: Dict[str, Path] = {}
    if cfg_dir.is_dir():
        for path in sorted(cfg_dir.glob("*.yaml")):
            configs[path.stem] = path
    init_path = configs.get("config_init")
    if init_path is None or not init_path.is_file():
        raise MissingPipelineError(
            f"{root} has no .log/configs/config_init.yaml "
            "(run samovar prepare/exec on the source pipeline first)"
        )
    start, end = _parse_window_env(root / ".log" / "window.env")
    repro = _load_yaml(configs.get("config_reprofiling"))
    model = root / "reprofiled_annotations" / "trained_model.joblib"
    if not model.is_file():
        alt = Path(str(repro.get("output_dir") or "")) / "trained_model.joblib"
        model = alt if alt.is_file() else model
    regen = root / "regenerated_annotations"
    if not regen.is_dir():
        regen_path = Path(str(repro.get("regenerated_path") or ""))
        if regen_path.is_file():
            regen = regen_path.parent
        elif regen_path.is_dir():
            regen = regen_path
        else:
            regen = None
    seed = int(repro.get("seed") or 42)
    return PipelineState(
        root=root,
        start=start,
        end=end,
        configs=configs,
        init_config=_load_yaml(init_path),
        qc_config=_load_yaml(configs.get("config_qc")),
        scoring_config=_load_yaml(configs.get("config_scoring")),
        export_config=_load_yaml(configs.get("config_export")),
        reprofiling_config=repro,
        annotation2iss_config=_load_yaml(configs.get("config_annotation2iss")),
        model_path=model if model.is_file() else None,
        regenerated_annotations=regen if regen is not None and Path(regen).exists() else None,
        seed=seed,
    )


def configured_steps(state: PipelineState, *, full: bool) -> List[str]:
    """Steps from setup through the source endpoint, minus apply-mode skips.

    A new sample always starts at ``setup_reads`` even if the source window
    began later. In-silico read regeneration is skipped unless
    ``SAMOVAR_APPLY_ISS=1``. Normal apply skips table regeneration; ``--full``
    keeps those checkpoints so the reprofiler can refit.
    """
    skip_iss = os.environ.get("SAMOVAR_APPLY_ISS", "").strip().lower() not in {
        "1",
        "true",
        "yes",
    }
    _, end = resolve_window(CHECKPOINT_STEPS[0], state.end)
    end_idx = CHECKPOINT_STEPS.index(end)
    out: List[str] = []
    for name in CHECKPOINT_STEPS:
        if CHECKPOINT_STEPS.index(name) > end_idx:
            break
        if skip_iss and name in ISS_STEPS:
            continue
        if not full and name in TABLE_ADAPT_STEPS:
            continue
        out.append(name)
    return out


def _set_seed(seed: int) -> None:
    random.seed(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except Exception:
        pass


def _apply_env(dest: Path) -> Dict[str, str]:
    env = os.environ.copy()
    root = repo_root()
    prefix = runtime_path_prefix()
    env["PATH"] = f"{prefix}:{root / 'bin'}:{env.get('PATH', '')}"
    env["PYTHONPATH"] = f"{root / 'src'}{os.pathsep}{env['PYTHONPATH']}" if env.get("PYTHONPATH") else str(root / "src")
    env["SAMOVAR_ROOT"] = str(root)
    env["SAMOVAR_RUN_DIR"] = str(dest)
    return env


def _retarget_init_config(cfg: Dict[str, Any], dest: Path, cores: int) -> Dict[str, Any]:
    out = copy.deepcopy(cfg)
    trimmed = dest / "initial_trimmed"
    out["r1_dir"] = str(trimmed)
    out["r2_dir"] = str(trimmed)
    out["output_dir"] = str(dest / "initial_reports")
    for run in out.get("run_config") or []:
        if isinstance(run, dict):
            run["threads"] = cores
    return out


def _retarget_annotation2iss(cfg: Dict[str, Any], dest: Path, cores: int) -> Dict[str, Any]:
    out = copy.deepcopy(cfg) if cfg else {}
    out["annotation_dir"] = str(dest / "initial_annotations")
    out["observed_abundance_dir"] = str(dest / "initial_abundance")
    out["abundance_dir"] = str(dest / "regenerated" / ".regenerated_abundance")
    out["output_dir"] = str(dest / "regenerated")
    out["cores"] = cores
    return out


def _run_annotators(config: Dict[str, Any], dest: Path, cores: int, env: Dict[str, str]) -> None:
    cfg_path = dest / ".log" / "configs" / "config_init.yaml"
    _dump_yaml(cfg_path, config)
    snakefile = repo_root() / "workflow" / "annotators" / "Snakefile"
    snakemake = shutil.which("snakemake")
    if snakemake:
        cmd = [
            snakemake,
            "-s",
            str(snakefile),
            "--configfile",
            str(cfg_path),
            "--directory",
            str(dest),
            "--cores",
            str(cores),
        ]
        subprocess.check_call(cmd, env=env)
        return
    from samovar.annotators_wrapper import get_annotator_instance
    from samovar.seqio import find_fastq_mate, list_fastq_samples

    reports = Path(config["output_dir"])
    reports.mkdir(parents=True, exist_ok=True)
    samples = list_fastq_samples(config["r1_dir"])
    if not samples:
        raise IncompatibleInputError(f"no FASTQ samples under {config['r1_dir']}")
    for run_conf in config.get("run_config") or []:
        tool = str(run_conf.get("type") or "")
        annotator = get_annotator_instance(tool, run_conf, config)
        for sample in samples:
            r1 = find_fastq_mate(config["r1_dir"], sample, "R1")
            r2 = find_fastq_mate(config["r2_dir"], sample, "R2")
            if r1 is None:
                raise IncompatibleInputError(f"missing R1 for sample {sample}")
            outs = annotator.get_expected_outputs(sample, str(reports))
            Path(outs[0]).parent.mkdir(parents=True, exist_ok=True)
            cmd = annotator.get_snakemake_shell_cmd(str(r1), str(r2) if r2 else str(r1), outs)
            subprocess.check_call(cmd, shell=True, env=env)


def _export_stage(
    dest: Path,
    stage: str,
    source_rel: str,
    export_cfg: Dict[str, Any],
    *,
    reference: Optional[Path],
    model: Optional[Path],
) -> None:
    if is_skipped_export(export_cfg.get("export") or export_cfg.get("corrector")):
        return
    formats = export_cfg.get("export_formats") or []
    if not formats and not export_cfg.get("export"):
        return
    cfg = dict(export_cfg)
    if reference is not None:
        cfg["reference"] = str(reference)
    if model is not None:
        cfg["model"] = str(model)
        cfg["model_path"] = str(model)
    export_run_dir(dest, stage=stage, source=source_rel, config=cfg)


def _viz_and_score(
    dest: Path,
    annotation_dir: Path,
    plots_dir: Path,
    scoring_cfg: Dict[str, Any],
    stage: str,
    *,
    combined: Optional[Path] = None,
) -> None:
    try:
        compare_annotations(
            annotation_dir=str(annotation_dir),
            output_dir=str(plots_dir),
            csv_file=str(combined) if combined else None,
            show_top=0,
            rank="none",
        )
    except Exception as exc:
        print(f"[apply] visualization skipped ({stage}): {exc}", file=sys.stderr)
    cfg = dict(scoring_cfg)
    cfg["output_dir"] = str(dest)
    run_custom_scorers(dest, config=cfg, stage=stage)


def _copy_regenerated_annotations(source: Path, dest: Path) -> Path:
    target = dest / "regenerated_annotations"
    if target.exists():
        return target
    if source.is_file():
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target / source.name)
        return target
    shutil.copytree(source, target)
    return target


def apply_pipeline(
    input_dir: PathLike,
    pipeline_dir: PathLike,
    output_dir: PathLike,
    *,
    full: bool = False,
    cores: int = 1,
) -> Dict[str, Any]:
    """Run apply. Returns provenance dict written to ``.log/apply.yaml``."""
    src_reads = Path(input_dir).expanduser().resolve()
    dest = Path(output_dir).expanduser().resolve()
    if not src_reads.is_dir():
        raise IncompatibleInputError(f"input directory not found: {src_reads}")
    if not has_r1_reads(src_reads) and not list_r1_files(src_reads):
        raise IncompatibleInputError(
            f"no paired FASTQ (*_R1.fastq) under {src_reads}"
        )
    state = load_pipeline_state(pipeline_dir)
    if dest == state.root:
        raise ApplyError(
            "apply --output_dir must not be the source --pipeline directory "
            "(that would overwrite the trained run)"
        )
    if not full and state.model_path is None:
        raise MissingModelError(
            f"normal apply requires {state.root}/reprofiled_annotations/trained_model.joblib"
        )
    if full and state.regenerated_annotations is None:
        raise ApplyError(
            f"--full requires regenerated annotation tables with known true taxids "
            f"under {state.root}/regenerated_annotations"
        )
    dest.mkdir(parents=True, exist_ok=True)
    (dest / ".log" / "configs").mkdir(parents=True, exist_ok=True)
    _set_seed(state.seed)
    env = _apply_env(dest)
    steps = configured_steps(state, full=full)
    cores = max(1, int(cores or 1))

    scoring_cfg = dict(state.scoring_config)
    scoring_cfg["output_dir"] = str(dest)
    _dump_yaml(dest / ".log" / "configs" / "config_scoring.yaml", scoring_cfg)
    qc_cfg = dict(state.qc_config)
    qc_cfg["output_dir"] = str(dest)
    _dump_yaml(dest / ".log" / "configs" / "config_qc.yaml", qc_cfg)
    export_cfg = dict(state.export_config)
    export_cfg["output_dir"] = str(dest)
    _dump_yaml(dest / ".log" / "configs" / "config_export.yaml", export_cfg)
    a2iss = _retarget_annotation2iss(state.annotation2iss_config, dest, cores)
    _dump_yaml(dest / ".log" / "configs" / "config_annotation2iss.yaml", a2iss)
    repro_cfg = dict(state.reprofiling_config)
    repro_cfg["output_dir"] = str(dest / "reprofiled_annotations")
    repro_cfg["initial_dir"] = str(dest / "initial_annotations")
    repro_cfg["seed"] = state.seed
    _dump_yaml(dest / ".log" / "configs" / "config_reprofiling.yaml", repro_cfg)

    init_cfg = _retarget_init_config(state.init_config, dest, cores)
    reference = (
        Path(state.regenerated_annotations)
        if state.regenerated_annotations is not None
        else None
    )
    applied_model = state.model_path
    trained = False

    for step in steps:
        if step == "setup_reads":
            (dest / "initial").mkdir(parents=True, exist_ok=True)
            link_or_copy_reads(src_reads, dest / "initial")
        elif step == "qc_initial":
            trim_stage(dest, "initial", qc_cfg)
        elif step == "annotate_initial":
            _run_annotators(init_cfg, dest, cores, env)
        elif step == "combine_initial":
            combine_with_cpp(
                str(dest / "initial_reports"),
                str(dest / "initial_annotations"),
                1,
            )
            _export_stage(
                dest,
                "initial",
                "initial_annotations",
                export_cfg,
                reference=reference,
                model=applied_model,
            )
        elif step == "viz_initial":
            _viz_and_score(
                dest,
                dest / "initial_annotations",
                dest / "initial_annotations_plots",
                scoring_cfg,
                "viz_initial",
            )
        elif step == "abundance_tables":
            materialize_observed_abundance(dest)
        elif step == "regenerate_tables":
            stage_regenerate_tables(dest, a2iss)
        elif step == "score_regenerated_tables":
            try:
                stage_score_regenerated_tables(dest, a2iss)
            except FileNotFoundError as exc:
                print(f"[apply] table scoring skipped: {exc}", file=sys.stderr)
        elif step == "reprofile":
            if full:
                copied = _copy_regenerated_annotations(
                    Path(state.regenerated_annotations), dest
                )
                gt_dir = regenerated_abundance_dir(dest)
                name = (
                    repro_cfg.get("reprofiler")
                    or state.reprofiling_config.get("reprofiler")
                    or "ensemble"
                )
                run_reprofiler(
                    name,
                    regenerated_path=copied,
                    ground_truth_dir=gt_dir if gt_dir.is_dir() else None,
                    initial_dir=dest / "initial_annotations",
                    output_dir=dest / "reprofiled_annotations",
                    config=repro_cfg,
                )
                trained = True
                applied_model = dest / "reprofiled_annotations" / "trained_model.joblib"
                reference = copied
            else:
                apply_saved_reprofiler(
                    model_path=state.model_path,
                    initial_dir=dest / "initial_annotations",
                    output_dir=dest / "reprofiled_annotations",
                    config=repro_cfg,
                )
                applied_model = dest / "reprofiled_annotations" / "trained_model.joblib"
            _export_stage(
                dest,
                "reprofiled",
                "reprofiled_annotations",
                export_cfg,
                reference=reference,
                model=applied_model,
            )
        elif step == "viz_reprofiled":
            combined = dest / "reprofiled_annotations" / "combined_annotation_table.csv"
            _viz_and_score(
                dest,
                dest / "reprofiled_annotations",
                dest / "reprofiled_annotations_plots",
                scoring_cfg,
                "viz_reprofiled",
                combined=combined if combined.parent.is_dir() else None,
            )
        else:
            print(f"[apply] skip unimplemented source step {step}", file=sys.stderr)
            continue
        mark_done(dest, step)

    provenance = {
        "mode": "full" if full else "normal",
        "pipeline": str(state.root),
        "input_dir": str(src_reads),
        "output_dir": str(dest),
        "model": str(applied_model) if applied_model else "",
        "source_model": str(state.model_path) if state.model_path else "",
        "seed": state.seed,
        "cores": cores,
        "steps": steps,
        "retrained": trained,
        "window": {"start": state.start, "end": state.end},
    }
    _dump_yaml(dest / ".log" / "apply.yaml", provenance)
    return provenance


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="samovar apply",
        description=(
            "Apply a configured SamovaR pipeline to one new sample. "
            "Normal mode uses the saved ML reprofiler; --full regenerates "
            "abundance tables and refits the configured reprofiler."
        ),
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory with one sample of paired FASTQ (*_R1 / *_R2)",
    )
    parser.add_argument(
        "--pipeline",
        "--pipeline_dir",
        dest="pipeline",
        required=True,
        help="Completed source run directory (prepare/exec outdir)",
    )
    add_output_dir_argument(parser, required=True)
    parser.add_argument(
        "--full",
        action="store_true",
        help="Regenerate tables from the new sample and refit the ML reprofiler",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Threads for annotators (default: 1)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    try:
        apply_pipeline(
            args.input_dir,
            args.pipeline,
            args.output_dir,
            full=bool(args.full),
            cores=args.cores,
        )
    except (
        ApplyError,
        FileNotFoundError,
        subprocess.CalledProcessError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
