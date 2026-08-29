"""Abundance-table regenerators: builtins plus imported ``table_reads_generator`` tools.

Each regenerator takes an abundance table (``taxid`` × ``N_<sample>``), an
``Annotation``, or a long per-read table, and returns
``{name: DataFrame(taxid, N_<sample>…)}``.
"""

from __future__ import annotations

import importlib.util
import shlex
import subprocess
import sys
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union

import pandas as pd

from samovar.main_config import (
    flags_target_matches,
    iter_tools,
    merge_flag_strings,
    parse_tool_entry,
    tool_flags,
    tool_path,
)
from samovar.parse_annotators import Annotation
from samovar.paths import load_config
from samovar.regenerate import (
    _n_samples_or_observed,
    coerce_seed,
    max_genomes_from_config,
    regenerate_bootstrap,
    regenerate_camisim,
    regenerate_glm_python,
    regenerate_preserve,
    regenerate_vae,
    resolve_regeneration_mode,
)

TABLE_READS_GENERATOR_GROUP = "table_reads_generator"


class MissingTableRegeneratorError(ValueError):
    """``tools.<name>`` is missing or is not a ``table_reads_generator``."""


def extra_flags_argv(text: Optional[str]) -> List[str]:
    if not text or not str(text).strip():
        return []
    return shlex.split(str(text))


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "t", "on"}


def parse_regeneration_modes(raw: Any) -> List[str]:
    """Split CLI/YAML ``table_reads_generator`` into ordered unique names."""
    if raw in (None, False, ""):
        return []
    if isinstance(raw, (list, tuple)):
        names: List[str] = []
        for item in raw:
            names.extend(parse_regeneration_modes(item))
        return names
    return [piece.strip() for piece in str(raw).replace(",", " ").split() if piece.strip()]


def canonical_regeneration_modes(raw: Any) -> List[str]:
    names: List[str] = []
    seen = set()
    for name in parse_regeneration_modes(raw):
        canon = require_known_regeneration_mode(name)
        if canon in seen:
            continue
        seen.add(canon)
        names.append(canon)
    return names


def flags_apply_to_regenerator(target: str, mode: Optional[str]) -> bool:
    """True if ``--flags TARGET …`` should attach to the current table regenerator."""
    kind, name = resolve_regeneration_mode(mode)
    names = [name, mode]
    if kind == "builtin" and name == "camisim-table":
        names.extend(["camisim", "camisim_table", "cami"])
    if kind == "builtin" and str(name).startswith("sparsedossa2"):
        names.extend(["sparsedossa2", "sd2", "SparseDOSSA2"])
    return flags_target_matches(
        target,
        *names,
        groups=("table_reads_generator", "table_reads", "table"),
    )


def apply_extra_flags(config: Dict[str, Any]) -> Dict[str, Any]:
    """Parse ``extra_flags`` / ``table_reads_generator_flags`` into config keys."""
    cfg = dict(config or {})
    raw = merge_flag_strings(
        cfg.get("extra_flags"),
        cfg.get("table_reads_generator_flags"),
    )
    argv = extra_flags_argv(raw)
    mapping = {
        "--log-mu": ("log_mu", float),
        "--log_mu": ("log_mu", float),
        "--log-sigma": ("log_sigma", float),
        "--log_sigma": ("log_sigma", float),
        "--gauss-mu": ("gauss_mu", float),
        "--gauss_mu": ("gauss_mu", float),
        "--gauss-sigma": ("gauss_sigma", float),
        "--gauss_sigma": ("gauss_sigma", float),
        "--distribution": ("camisim_distribution", str),
        "--mode": ("camisim_distribution", str),
        "--N": ("N", int),
        "--N_reads": ("N_reads", int),
        "--seed": ("seed", int),
        "--template": ("sparsedossa2_template", str),
        "--n-feature": ("n_feature", int),
        "--n_feature": ("n_feature", int),
        "--n-features": ("n_feature", int),
        "--workers": ("sparsedossa2_workers", int),
        "--parallel": ("sparsedossa2_workers", int),
        "--cores": ("cores", int),
        "--cv-folds": ("cv_folds", int),
        "--K": ("cv_folds", int),
        "--lambdas": ("lambdas", str),
        "--lambda": ("fit_lambda", float),
        "--maxit": ("maxit", int),
        "--max-eval": ("max_eval", int),
        "--max_eval": ("max_eval", int),
        "--prec-bits": ("prec_bits", int),
        "--precBits": ("prec_bits", int),
        "--timeout": ("timeout", int),
        "--max-genomes": ("max_genomes", str),
        "--max_genomes": ("max_genomes", str),
    }
    bool_flags = {
        "--fit": ("sparsedossa2_fit", True),
        "--verbose": ("verbose", True),
        "-v": ("verbose", True),
    }
    optional_bool = {
        "--new-features": "new_features",
        "--new_features": "new_features",
    }
    i = 0
    while i < len(argv):
        tok = argv[i]
        key = None
        caster = str
        value = None
        consumed = 1
        if tok.startswith("--") and "=" in tok:
            flag, value = tok.split("=", 1)
            mapped = mapping.get(flag)
            if mapped:
                key, caster = mapped
            elif flag in optional_bool:
                key, caster = optional_bool[flag], _as_bool
            elif flag in bool_flags:
                key, value = bool_flags[flag][0], bool_flags[flag][1]
                caster = lambda x: x
        elif tok in bool_flags:
            key, value = bool_flags[tok][0], bool_flags[tok][1]
            caster = lambda x: x
        elif tok in optional_bool:
            if i + 1 < len(argv) and not str(argv[i + 1]).startswith("-"):
                key, caster, value, consumed = optional_bool[tok], _as_bool, argv[i + 1], 2
            else:
                key, value, caster = optional_bool[tok], True, lambda x: x
        else:
            mapped = mapping.get(tok)
            if mapped and i + 1 < len(argv):
                key, caster = mapped
                value = argv[i + 1]
                consumed = 2
        if key is not None and value is not None:
            try:
                cfg[key] = caster(value)
            except (TypeError, ValueError):
                cfg[key] = value
        i += consumed
    cfg["extra_argv"] = argv
    cfg["extra_flags"] = raw
    return cfg


def attach_regenerator_flags(mode: Optional[str], config: Dict[str, Any]) -> Dict[str, Any]:
    """Merge import-time flags (tools.*[4]) with prepare ``--flags`` / YAML."""
    cfg = dict(config or {})
    kind, name = resolve_regeneration_mode(mode)
    imported = ""
    if kind == "custom":
        spec = lookup_table_regenerator(name)
        imported = tool_flags(spec, name)
    cfg["extra_flags"] = merge_flag_strings(
        imported,
        cfg.get("extra_flags"),
        cfg.get("table_reads_generator_flags"),
    )
    return apply_extra_flags(cfg)


def as_annotation(
    data: Union[Annotation, pd.DataFrame, None],
    annotation_dir: Union[str, Path, None] = None,
) -> Annotation:
    from samovar.abundance import is_abundance_table, load_table_input

    if isinstance(data, Annotation):
        return data
    if annotation_dir is not None and data is None:
        return load_table_input(annotation_dir)
    if data is not None:
        if is_abundance_table(data):
            return Annotation.from_abundance_tables({"table": data})
        return Annotation.from_long_table(data)
    if annotation_dir is not None:
        return load_table_input(annotation_dir)
    return Annotation.from_long_table(pd.DataFrame())


def load_samples_metadata(config: Optional[Dict[str, Any]]) -> Optional[pd.DataFrame]:
    """Optional samples table. Built-ins ignore it; ``None`` is valid."""
    cfg = dict(config or {})
    raw = cfg.get("samples_metadata") or cfg.get("metadata")
    if raw is None or raw is False or str(raw).strip().lower() in {"", "null", "none"}:
        return None
    if isinstance(raw, pd.DataFrame):
        return raw
    path = Path(str(raw)).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"samples_metadata not found: {path}")
    return pd.read_csv(path)


def lookup_table_regenerator(name: str) -> list:
    """Return ``[env, workflow, path, group]`` for an imported table regenerator.

    Raises ``MissingTableRegeneratorError`` if the name is absent from the
    install config or registered under another ``--type``.
    """
    key = str(name or "").strip()
    if not key:
        raise MissingTableRegeneratorError(
            "Empty table_reads_generator name. Import a tool with "
            "`samovar tools import -n NAME --type table`."
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
        raise MissingTableRegeneratorError(
            f"table_reads_generator {key!r} is not in the main install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type table` "
            "before prepare / regenerate."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != TABLE_READS_GENERATOR_GROUP:
        raise MissingTableRegeneratorError(
            f"tools.{matched} has group {group!r}, expected "
            f"{TABLE_READS_GENERATOR_GROUP!r}. Re-import with --type table."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingTableRegeneratorError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return parsed


def require_known_regeneration_mode(mode: Optional[str]) -> str:
    """Canonical builtin name, or custom name after checking the install config."""
    kind, name = resolve_regeneration_mode(mode)
    if kind == "custom":
        lookup_table_regenerator(name)
    return name


class TableRegenerator(ABC):
    """One abundance table in, one (or more) abundance tables out."""

    @abstractmethod
    def run(
        self,
        annotation: Annotation,
        metadata: Optional[pd.DataFrame],
        config: Dict[str, Any],
    ) -> Dict[str, pd.DataFrame]:
        raise NotImplementedError


class DirectTableRegenerator(TableRegenerator):
    def run(self, annotation, metadata, config):
        _ = metadata
        cfg = dict(config or {})
        n_reads = cfg.get("N_reads")
        if n_reads is not None:
            n_reads = int(n_reads)
        return regenerate_preserve(
            annotation,
            n_reads=n_reads,
            rescale=bool(cfg.get("rescale_abundance", False)),
            threshold_amount=float(
                cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5))
            ),
            max_genomes=max_genomes_from_config(cfg),
        )


class BootstrapTableRegenerator(TableRegenerator):
    def run(self, annotation, metadata, config):
        _ = metadata
        cfg = dict(config or {})
        n_reads = cfg.get("N_reads")
        if n_reads is not None:
            n_reads = int(n_reads)
        return regenerate_bootstrap(
            annotation,
            n_samples=_n_samples_or_observed(annotation, cfg.get("N")),
            n_reads=n_reads,
            threshold_amount=float(
                cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5))
            ),
            seed=coerce_seed(cfg.get("seed", 42)),
            rescale=bool(cfg.get("rescale_abundance", False)),
            error_scale=float(cfg.get("bootstrap_error_scale", 0.15)),
            max_genomes=max_genomes_from_config(cfg),
        )


class VaeTableRegenerator(TableRegenerator):
    def run(self, annotation, metadata, config):
        _ = metadata
        cfg = dict(config or {})
        n_reads = cfg.get("N_reads")
        if n_reads is not None:
            n_reads = int(n_reads)
        return regenerate_vae(
            annotation,
            n_samples=_n_samples_or_observed(annotation, cfg.get("N")),
            n_reads=n_reads,
            threshold_amount=float(
                cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5))
            ),
            latent_dim=int(cfg.get("vae_latent_dim", 4)),
            seed=coerce_seed(cfg.get("seed", 42)),
            rescale=bool(cfg.get("rescale_abundance", False)),
            max_genomes=max_genomes_from_config(cfg),
        )


class GlmTableRegenerator(TableRegenerator):
    def run(self, annotation, metadata, config):
        _ = metadata
        cfg = dict(config or {})
        n_reads = cfg.get("N_reads")
        if n_reads is not None:
            n_reads = int(n_reads)
        return regenerate_glm_python(
            annotation,
            n_samples=_n_samples_or_observed(annotation, cfg.get("N")),
            n_reads=n_reads,
            threshold_amount=float(
                cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5))
            ),
            seed=coerce_seed(cfg.get("seed", 42)),
            rescale=bool(cfg.get("rescale_abundance", False)),
            min_cluster_size=int(cfg.get("min_cluster_size", 2)),
            max_cluster_size=int(cfg.get("max_cluster_size", 100)),
            max_genomes=max_genomes_from_config(cfg),
        )


class SamovarRTableRegenerator(TableRegenerator):
    def run(self, annotation, metadata, config):
        from samovar.abundance import input_to_abundance_tables
        from samovar.regenerate import _write_abundance_tables, write_samovar_config_defaults
        from samovar.table2iss import (
            R_REGENERATE_DRIVER,
            _resolve_r_executable,
        )
        import os
        import tempfile
        import yaml

        _ = metadata
        cfg = write_samovar_config_defaults(dict(config or {}))
        out = Path(str(cfg.get("output_dir") or "regenerated_abundance"))
        out.mkdir(parents=True, exist_ok=True)
        tables = input_to_abundance_tables(annotation)
        staging = out / ".samovar_r_input"
        if tables:
            _write_abundance_tables(staging, tables)
            input_dir = str(staging)
        else:
            input_dir = str(cfg.get("annotation_dir") or "")
            if not input_dir:
                raise ValueError(
                    "regeneration_mode='samovar' needs an abundance table or annotation_dir."
                )
        r_cfg = {
            "N": int(cfg.get("N") or 1),
            "N_reads": int(cfg.get("N_reads") or 100),
            "regeneration_mode": "samovar",
            "seed": cfg.get("seed", 42),
            "threshold_amount": cfg.get("threshold_amount", 1e-5),
            "plot_log": False,
        }
        r_path, r_lib_path = _resolve_r_executable()
        env = os.environ.copy()
        if r_lib_path:
            env["R_LIBS"] = r_lib_path
            env["R_LIBS_USER"] = r_lib_path
        cfg_tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, encoding="utf-8"
        )
        yaml.safe_dump(r_cfg, cfg_tmp)
        cfg_tmp.close()
        drv_tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".R", delete=False, encoding="utf-8"
        )
        drv_tmp.write(R_REGENERATE_DRIVER)
        drv_tmp.close()
        cmd = [
            r_path,
            "--vanilla",
            "-s",
            "-f",
            drv_tmp.name,
            "--args",
            "--config",
            cfg_tmp.name,
            "--annotation_dir",
            input_dir,
            "--output_dir",
            str(out),
        ]
        try:
            subprocess.run(cmd, check=True, env=env)
        except FileNotFoundError as exc:
            raise ValueError(
                "regeneration_mode='samovar' uses the optional R regenerator. "
                "Install samovaR (./install.sh R-package) and set r_path."
            ) from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"samovaR regenerator failed (exit {exc.returncode})"
            ) from exc
        finally:
            Path(cfg_tmp.name).unlink(missing_ok=True)
            Path(drv_tmp.name).unlink(missing_ok=True)
        written = {
            p.stem: pd.read_csv(p)
            for p in sorted(out.glob("*.csv"))
            if "taxid" in pd.read_csv(p, nrows=0).columns
        }
        if written:
            return written
        raise ValueError(
            "regeneration_mode='samovar' uses the optional R regenerator. "
            "Install samovaR (./install.sh R-package) and set SAMOVAR_R_REGENERATE."
        )


class CamisimTableRegenerator(TableRegenerator):
    def run(self, annotation, metadata, config):
        _ = metadata
        cfg = apply_extra_flags(dict(config or {}))
        n_reads = cfg.get("N_reads")
        if n_reads is not None:
            n_reads = int(n_reads)
        return regenerate_camisim(
            annotation,
            n_samples=_n_samples_or_observed(annotation, cfg.get("N")),
            n_reads=n_reads,
            threshold_amount=float(
                cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5))
            ),
            seed=coerce_seed(cfg.get("seed", 42)),
            rescale=bool(cfg.get("rescale_abundance", False)),
            distribution=str(
                cfg.get("camisim_distribution")
                or cfg.get("distribution")
                or "differential"
            ),
            log_mu=float(cfg.get("log_mu", 1.0)),
            log_sigma=float(cfg.get("log_sigma", 2.0)),
            gauss_mu=float(cfg.get("gauss_mu", 1.0)),
            gauss_sigma=float(cfg.get("gauss_sigma", 1.0)),
            max_genomes=max_genomes_from_config(cfg),
        )


class CustomTableRegenerator(TableRegenerator):
    """User tool registered as ``tools.<name>`` with group ``table_reads_generator``."""

    def __init__(self, name: str):
        self.name = str(name).strip()

    def run(self, annotation, metadata, config):
        spec = lookup_table_regenerator(self.name)
        path = Path(tool_path(spec, self.name)).expanduser()
        cfg = dict(config or {})
        py_fn = _try_python_callable(path, self.name)
        if py_fn is not None:
            tables = py_fn(annotation, metadata, cfg)
            if not isinstance(tables, dict):
                raise TypeError(
                    f"{path} regenerate() must return dict[str, DataFrame], "
                    f"got {type(tables).__name__}"
                )
            return tables
        return _run_cli_regenerator(path, annotation, metadata, cfg)


def _try_python_callable(path: Path, name: str) -> Optional[Callable]:
    if path.suffix.lower() != ".py" or not path.is_file():
        return None
    try:
        spec = importlib.util.spec_from_file_location(
            f"samovar_custom_regen_{name}", path
        )
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception:
        return None
    fn = getattr(module, "regenerate", None)
    if callable(fn):
        return fn
    cls = getattr(module, "TableRegenerator", None)
    if cls is None:
        return None
    try:
        inst = cls()
    except Exception:
        return None
    run = getattr(inst, "run", None)
    return run if callable(run) else None


def _run_cli_regenerator(
    path: Path,
    annotation: Annotation,
    metadata: Optional[pd.DataFrame],
    config: Dict[str, Any],
) -> Dict[str, pd.DataFrame]:
    annotation_dir = config.get("annotation_dir")
    output_dir = config.get("output_dir")
    if not output_dir:
        raise ValueError(
            f"CLI table_reads_generator {path} needs output_dir "
            "in config (or define regenerate() in a Python module)."
        )
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    from samovar.abundance import dir_looks_like_annotation, input_to_abundance_tables
    from samovar.regenerate import _write_abundance_tables

    if annotation_dir and dir_looks_like_annotation(annotation_dir):
        input_path = str(annotation_dir)
    else:
        staging = out / ".input_abundance"
        tables = input_to_abundance_tables(annotation)
        if not tables:
            raise ValueError(
                f"CLI table_reads_generator {path} needs annotation_dir or an abundance table."
            )
        _write_abundance_tables(staging, tables)
        input_path = str(staging)
    cmd = [str(path), "-i", input_path, "-o", str(out)]
    meta_tmp = None
    try:
        if metadata is not None:
            handle = tempfile.NamedTemporaryFile(
                mode="w", suffix=".csv", delete=False, encoding="utf-8"
            )
            metadata.to_csv(handle.name, index=False)
            handle.close()
            meta_tmp = handle.name
            cmd.extend(["-m", meta_tmp])
        if config.get("seed") is not None:
            cmd.extend(["--seed", str(coerce_seed(config.get("seed")))])
        if config.get("N") not in (None, "", 0):
            cmd.extend(["--N", str(int(config["N"]))])
        if config.get("N_reads") is not None:
            cmd.extend(["--N_reads", str(int(config["N_reads"]))])
        from samovar.regenerate import finite_max_genomes

        limit = finite_max_genomes(config.get("max_genomes"), default_from_env=False)
        if limit is not None:
            cmd.extend(["--max-genomes", str(limit)])
        extra = list(config.get("extra_argv") or extra_flags_argv(config.get("extra_flags")))
        cmd.extend(extra)
        if path.suffix.lower() == ".py":
            cmd = [sys.executable, *cmd]
        subprocess.run(cmd, check=True)
    finally:
        if meta_tmp:
            Path(meta_tmp).unlink(missing_ok=True)
    _ = annotation
    return _read_abundance_dir(out)


def _read_abundance_dir(output_dir: Path) -> Dict[str, pd.DataFrame]:
    tables: Dict[str, pd.DataFrame] = {}
    for csv in sorted(output_dir.glob("*.csv")):
        df = pd.read_csv(csv)
        if "taxid" not in df.columns:
            continue
        tables[csv.stem] = df
    if not tables:
        raise ValueError(
            f"Custom table regenerator wrote no taxid CSVs under {output_dir}"
        )
    return tables


class SparseDOSSA2TableRegenerator(TableRegenerator):
    def __init__(self, mode: str = "sparsedossa2-fit"):
        self.mode = str(mode)

    def run(self, annotation, metadata, config):
        cfg = dict(config or {})
        cfg["table_reads_generator"] = self.mode
        cfg["regeneration_mode"] = self.mode
        from samovar.sparsedossa2 import regenerate as sd2_regenerate

        return sd2_regenerate(annotation, metadata, cfg)


_BUILTIN: Dict[str, type] = {
    "direct": DirectTableRegenerator,
    "bootstrap": BootstrapTableRegenerator,
    "vae": VaeTableRegenerator,
    "glm": GlmTableRegenerator,
    "samovar": SamovarRTableRegenerator,
    "camisim-table": CamisimTableRegenerator,
    "sparsedossa2-fit": SparseDOSSA2TableRegenerator,
    "sparsedossa2-stool": SparseDOSSA2TableRegenerator,
    "sparsedossa2-vaginal": SparseDOSSA2TableRegenerator,
    "sparsedossa2-ibd": SparseDOSSA2TableRegenerator,
}


def get_table_regenerator(mode: Optional[str]) -> TableRegenerator:
    """Factory: builtin mode or imported ``tools.<name>`` custom regenerator."""
    kind, name = resolve_regeneration_mode(mode)
    if kind == "builtin":
        cls = _BUILTIN[name]
        if cls is SparseDOSSA2TableRegenerator:
            return SparseDOSSA2TableRegenerator(name)
        return cls()
    lookup_table_regenerator(name)
    return CustomTableRegenerator(name)
