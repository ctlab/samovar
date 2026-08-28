"""SparseDOSSA2 abundance-table generator and CV table scorer.

Optional R package: https://github.com/biobakery/SparseDOSSA2
Wiki: https://github.com/biobakery/biobakery/wiki/SparseDOSSA2

Install with ``./install.sh SparseDOSSA2``. Parallelism uses the R ``future``
plan from the wiki (``multisession``; CV: sequential/sequential/multisession).
"""

from __future__ import annotations

import itertools
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import pandas as pd

from samovar.parse_annotators import Annotation
from samovar.paths import load_config, repo_root
from samovar.regenerate import (
    SPARSEDOSSA2_MODE_ALIASES,
    _abundance_table_from_matrix,
    _annotator_name,
    _apply_abundance_scale,
    _count_matrix,
    _filter_taxa,
    _n_samples_or_observed,
    _taxid_columns,
    coerce_seed,
)
from samovar.table_regenerators import apply_extra_flags

PathLike = Union[str, Path]

TEMPLATES = ("Stool", "Vaginal", "IBD")
MODE_ALIASES = SPARSEDOSSA2_MODE_ALIASES
_JOB_SEQ = itertools.count(1)


def _job_tag(prefix: str) -> str:
    return f"{prefix}_{os.getpid()}_{next(_JOB_SEQ)}"


class MissingSparseDOSSA2Error(RuntimeError):
    """R package SparseDOSSA2 is not installed."""


def canonicalize_sparsedossa2_mode(name: Optional[str]) -> Optional[str]:
    key = str(name or "").strip().lower().replace("_", "-")
    return MODE_ALIASES.get(key)


def template_for_mode(mode: str) -> Optional[str]:
    canon = canonicalize_sparsedossa2_mode(mode) or str(mode)
    if canon.endswith("-stool"):
        return "Stool"
    if canon.endswith("-vaginal"):
        return "Vaginal"
    if canon.endswith("-ibd"):
        return "IBD"
    if canon.endswith("-fit") or canon == "sparsedossa2":
        return None
    return None


def rscript_executable() -> str:
    cfg = load_config()
    try:
        from samovar.main_config import iter_tools, tool_path

        spec = iter_tools(cfg).get("Rscript")
        stored = tool_path(spec, "Rscript") if spec else ""
        if stored and Path(stored).expanduser().is_file():
            return str(Path(stored).expanduser())
    except Exception:
        pass
    for key in ("rscript_path", "Rscript"):
        raw = str(cfg.get(key) or "").strip()
        if raw and Path(raw).expanduser().is_file():
            return str(Path(raw).expanduser())
    r_path = str(cfg.get("r_path") or "").strip()
    if r_path:
        candidate = Path(r_path).expanduser()
        if candidate.name in {"R", "R.exe"}:
            sibling = candidate.with_name("Rscript")
            if sibling.is_file():
                return str(sibling)
        rs = candidate.parent / "Rscript"
        if rs.is_file():
            return str(rs)
    found = shutil.which("Rscript")
    if found:
        return found
    raise MissingSparseDOSSA2Error(
        "Rscript is not on PATH. Install R, then ./install.sh SparseDOSSA2."
    )


def driver_script() -> Path:
    env = os.environ.get("SAMOVAR_SPARSEDOSSA2_R", "").strip()
    if env:
        path = Path(env).expanduser()
        if path.is_file():
            return path
    cfg = load_config()
    raw = str(cfg.get("sparsedossa2_r") or "").strip()
    if raw:
        path = Path(raw).expanduser()
        if path.is_file():
            return path
    try:
        from samovar.main_config import iter_tools, tool_path

        spec = iter_tools(cfg).get("sparsedossa2.R")
        stored = tool_path(spec, "sparsedossa2.R") if spec else ""
        if stored and Path(stored).expanduser().is_file():
            return Path(stored).expanduser()
    except Exception:
        pass
    packaged = Path(__file__).with_name("sparsedossa2.R")
    if packaged.is_file():
        return packaged
    fallback = repo_root() / "src" / "samovar" / "sparsedossa2.R"
    if fallback.is_file():
        return fallback
    raise MissingSparseDOSSA2Error("SparseDOSSA2 R driver sparsedossa2.R is missing.")


def sparsedossa2_available() -> bool:
    try:
        rscript = rscript_executable()
    except MissingSparseDOSSA2Error:
        return False
    try:
        proc = subprocess.run(
            [
                rscript,
                "--vanilla",
                "-e",
                'cat(if (requireNamespace("SparseDOSSA2", quietly=TRUE)) "TRUE" else "FALSE")',
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return "TRUE" in (proc.stdout or "")


def require_sparsedossa2() -> None:
    if not sparsedossa2_available():
        raise MissingSparseDOSSA2Error(
            "SparseDOSSA2 is not installed. Run ./install.sh SparseDOSSA2 "
            "(https://github.com/biobakery/SparseDOSSA2)."
        )


def _workers(config: Dict[str, Any]) -> int:
    for key in ("sparsedossa2_workers", "workers", "parallel", "cores"):
        raw = config.get(key)
        if raw in (None, "", False):
            continue
        try:
            return max(1, int(raw))
        except (TypeError, ValueError):
            continue
    env = os.environ.get("SAMOVAR_SPARSEDOSSA2_WORKERS") or os.environ.get(
        "SAMOVAR_WORKERS"
    )
    if env:
        try:
            return max(1, int(env))
        except ValueError:
            pass
    return 1


def _r_env() -> Dict[str, str]:
    env = os.environ.copy()
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        env.setdefault(key, "1")
    return env


def _append_em_controls(argv: List[str], cfg: Dict[str, Any]) -> None:
    if cfg.get("maxit") not in (None, ""):
        argv.extend(["--maxit", str(int(cfg["maxit"]))])
    if cfg.get("max_eval") not in (None, ""):
        argv.extend(["--max-eval", str(int(cfg["max_eval"]))])
    if cfg.get("prec_bits") not in (None, ""):
        argv.extend(["--prec-bits", str(int(cfg["prec_bits"]))])


def _job_workdir(cfg: Optional[Dict[str, Any]] = None) -> Path:
    cfg = cfg or {}
    raw = (
        cfg.get("sd2_workdir")
        or cfg.get("output_dir")
        or os.environ.get("SAMOVAR_SD2_WORKDIR")
        or ""
    )
    if raw:
        path = Path(str(raw)).expanduser() / ".sd2_jobs"
        path.mkdir(parents=True, exist_ok=True)
        return path
    return Path(tempfile.mkdtemp(prefix="samovar_sd2_"))


def run_sparsedossa2_r(
    argv: Sequence[str],
    *,
    timeout: Optional[int] = None,
    config: Optional[Dict[str, Any]] = None,
) -> str:
    require_sparsedossa2()
    _ = config
    cmd = [rscript_executable(), "--vanilla", str(driver_script()), *[str(x) for x in argv]]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=_r_env(),
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"SparseDOSSA2 R driver timed out after {timeout}s"
        ) from exc
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(
            "SparseDOSSA2 R driver failed:\n" + err[-4000:]
        )
    return proc.stdout or ""


def _write_feature_matrix(matrix: pd.DataFrame, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    work = matrix.copy()
    work.index = work.index.astype(str)
    work.to_csv(dest)


def _read_simulated_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def remap_simulated_features(
    simulated: pd.DataFrame,
    taxids: Sequence[str],
) -> pd.DataFrame:
    """Keep taxid rownames when present; otherwise map Feature* onto ranked taxids."""
    wanted = [str(t) for t in taxids]
    sim = simulated.copy()
    sim.index = sim.index.astype(str)
    overlap = [t for t in wanted if t in sim.index]
    if overlap and len(overlap) >= min(3, len(wanted)):
        out = sim.reindex(wanted).fillna(0.0)
        return out
    n = min(len(wanted), sim.shape[0])
    out = pd.DataFrame(0.0, index=wanted, columns=list(sim.columns))
    if n:
        out.iloc[:n, :] = sim.iloc[:n, :].to_numpy()
    return out


def simulate_count_matrix(
    observed: pd.DataFrame,
    *,
    mode: str = "sparsedossa2-fit",
    n_sample: int,
    config: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    cfg = apply_extra_flags(dict(config or {}))
    ranked = list(observed.sum(axis=1).sort_values(ascending=False).index.astype(str))
    template = cfg.get("sparsedossa2_template") or template_for_mode(mode)
    if cfg.get("sparsedossa2_fit"):
        template = None
    n_feature = int(cfg.get("n_feature") or cfg.get("sparsedossa2_n_feature") or max(len(ranked), 1))
    new_features = cfg.get("new_features")
    workers = _workers(cfg)
    seed = coerce_seed(cfg.get("seed", 42))
    verbose = bool(cfg.get("verbose") or cfg.get("sparsedossa2_verbose"))
    work = _job_workdir(cfg)
    tag = _job_tag("sim")
    inp = work / f"{tag}_observed.csv"
    outp = work / f"{tag}_simulated.csv"
    _write_feature_matrix(observed, inp)
    argv = [
        "simulate",
        "--output",
        str(outp),
        "--n-sample",
        str(int(n_sample)),
        "--n-feature",
        str(n_feature),
        "--seed",
        str(seed),
        "--workers",
        str(workers),
    ]
    if verbose:
        argv.append("--verbose")
    if template:
        argv.extend(["--mode", "template", "--template", str(template)])
        if new_features is None:
            argv.extend(["--new-features", "TRUE"])
    else:
        argv.extend(["--mode", "fit", "--fit", "--input", str(inp)])
        if new_features is None:
            argv.extend(["--new-features", "FALSE"])
    if new_features is not None:
        argv.extend(["--new-features", "TRUE" if new_features else "FALSE"])
    if cfg.get("fit_lambda") not in (None, ""):
        argv.extend(["--lambda", str(cfg["fit_lambda"])])
    _append_em_controls(argv, cfg)
    run_sparsedossa2_r(
        argv,
        timeout=int(cfg.get("timeout") or 7200),
        config=cfg,
    )
    simulated = _read_simulated_csv(outp)
    return remap_simulated_features(simulated, ranked)


def _fitcv_argv(
    inp: Path,
    outp: Path,
    cfg: Dict[str, Any],
) -> List[str]:
    folds = int(cfg.get("cv_folds") or cfg.get("K") or 5)
    workers = _workers(cfg)
    seed = coerce_seed(cfg.get("seed", 42))
    verbose = bool(cfg.get("verbose") or cfg.get("sparsedossa2_verbose"))
    lambdas = cfg.get("lambdas") or cfg.get("sparsedossa2_lambdas")
    argv = [
        "fitcv",
        "--input",
        str(inp),
        "--output",
        str(outp),
        "--cv-folds",
        str(max(2, min(folds, 100))),
        "--seed",
        str(seed),
        "--workers",
        str(workers),
    ]
    if verbose:
        argv.append("--verbose")
    if lambdas not in (None, ""):
        argv.extend(["--lambdas", str(lambdas)])
    _append_em_controls(argv, cfg)
    return argv


def _parse_cv_payload(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    score = payload.get("cv_goodness_of_fit")
    try:
        score_f = float(score)
    except (TypeError, ValueError):
        score_f = float("nan")
    payload["cv_goodness_of_fit"] = score_f
    return payload


def concat_aligned_matrices(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    """Feature-align two feature×sample tables and concatenate samples."""
    idx = sorted(set(left.index.astype(str)) | set(right.index.astype(str)))
    a = left.copy()
    b = right.copy()
    a.index = a.index.astype(str)
    b.index = b.index.astype(str)
    a = a.reindex(idx).fillna(0.0)
    b = b.reindex(idx).fillna(0.0)
    b.columns = [f"{c}__b" for c in b.columns]
    return pd.concat([a, b], axis=1)


def fitcv_score_matrix(
    matrix: pd.DataFrame,
    *,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    rows = fitcv_score_jobs([("one", matrix)], config=config)
    return rows["one"]


def fitcv_score_jobs(
    items: Iterable[Tuple[str, pd.DataFrame]],
    *,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Dict[str, Any]]:
    require_sparsedossa2()
    cfg = apply_extra_flags(dict(config or {}))
    work = _job_workdir(cfg)
    timeout = int(cfg.get("timeout") or 7200)
    prepared: List[Tuple[str, Path, List[str]]] = []
    for key, matrix in items:
        tag = _job_tag("cv")
        inp = work / f"{tag}_{key}.csv"
        outp = work / f"{tag}_{key}.json"
        _write_feature_matrix(matrix, inp)
        argv = _fitcv_argv(inp, outp, cfg)
        prepared.append((str(key), outp, argv))

    out: Dict[str, Dict[str, Any]] = {}
    for key, outp, argv in prepared:
        try:
            run_sparsedossa2_r(argv, timeout=timeout, config=cfg)
            out[key] = _parse_cv_payload(outp)
        except Exception as exc:
            out[key] = {
                "cv_goodness_of_fit": float("nan"),
                "error": str(exc),
            }
    return out


def regenerate(
    annotation: Annotation,
    metadata: Optional[pd.DataFrame],
    config: Dict[str, Any],
) -> Dict[str, pd.DataFrame]:
    """``table_reads_generator`` entry point (fit or community template)."""
    _ = metadata
    cfg = apply_extra_flags(dict(config or {}))
    mode = str(
        cfg.get("table_reads_generator")
        or cfg.get("regeneration_mode")
        or "sparsedossa2-fit"
    )
    data = annotation.DataFrame
    if "sample" not in data.columns:
        data = data.copy()
        data["sample"] = "1"
    n_samples = _n_samples_or_observed(data, cfg.get("N"))
    n_reads = cfg.get("N_reads")
    if n_reads is not None:
        n_reads = int(n_reads)
    rescale = bool(cfg.get("rescale_abundance", False))
    threshold = float(cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5)))
    result: Dict[str, pd.DataFrame] = {}
    for col in _taxid_columns(data):
        mat = _count_matrix(data, col, "sample")
        mat = _filter_taxa(mat, threshold, n_reads, rescale=False)
        if mat.empty:
            continue
        simulated = simulate_count_matrix(mat, mode=mode, n_sample=n_samples, config=cfg)
        simulated.columns = [str(c) for c in simulated.columns]
        simulated = _apply_abundance_scale(simulated, n_reads, rescale)
        name = _annotator_name(col)
        result[name] = _abundance_table_from_matrix(simulated, name)
    if not result:
        raise RuntimeError("SparseDOSSA2 produced no per-annotator abundance tables")
    return result


def register_sparsedossa2_tools(cfg: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Write table generators + CV scorer into the install config."""
    from samovar.main_config import set_tool
    from samovar.paths import write_config

    config = dict(cfg or load_config())
    wrapper = str(Path(__file__).resolve())
    try:
        driver = str(driver_script().resolve())
    except MissingSparseDOSSA2Error:
        driver = str(Path(__file__).with_name("sparsedossa2.R"))
    rscript = shutil.which("Rscript") or ""
    try:
        rscript = rscript_executable()
    except MissingSparseDOSSA2Error:
        pass
    set_tool(
        config,
        "sparsedossa2.R",
        path=driver,
        env="",
        workflow="bash",
        group="workflow",
    )
    if rscript:
        set_tool(
            config,
            "Rscript",
            path=rscript,
            env="",
            workflow="bash",
            group="runtime",
        )
    table_flags = {
        "sparsedossa2-fit": "--fit",
        "sparsedossa2-stool": "--template Stool",
        "sparsedossa2-vaginal": "--template Vaginal",
        "sparsedossa2-ibd": "--template IBD",
    }
    for name, flags in table_flags.items():
        set_tool(
            config,
            name,
            path=wrapper,
            env="",
            workflow="bash",
            group="table_reads_generator",
            flags=flags,
        )
    set_tool(
        config,
        "sparsedossa2-cv",
        path=wrapper,
        env="",
        workflow="bash",
        group="scoring",
        flags="--cv-folds 5",
        inputs="regenerated/.regenerated_abundance",
    )
    write_config(config)
    return config
