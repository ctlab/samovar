"""Annotation → abundance-like export, with optional per-taxon correction.

Contract group ``export`` (``samovar tools import --type export``). Built-ins:

* ``logistic`` (default on ``samovar prepare``): on a reference Annotation with
  ``true`` / ``taxID_true``, estimate the per-taxon factor
  ``n_true / n_pred`` (logistic-smoothed precision / recall), then
  **multiply** predicted counts by that factor. Taxa that were never seen as
  predictions, or whose recall is below ``--min-efficiency``, are left unchanged.
  The correction is dropped for an annotator if it would lower abundance R² on
  the same reference.
* ``none`` / ``identity`` / ``raw``: count reads as-is.

Custom tools implement ``export(annotation, dest, config)`` and may read a
second Annotation from ``config["reference"]``.
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from samovar.abundance import (
    annotator_name,
    input_to_abundance_tables,
    load_table_input,
    n_sample_columns,
    normalize_abundance_table,
    taxid_annotation_columns,
    write_abundance_dir,
)
from samovar.annotation_convert import dump_builtin, parse_export_formats
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
from samovar.table_regenerators import extra_flags_argv

PathLike = Union[str, Path]
EXPORT_GROUP = "export"

_BUILTIN_ALIASES = {
    "logistic": "logistic",
    "logit": "logistic",
    "logreg": "logistic",
    "default": "logistic",
    "none": "identity",
    "identity": "identity",
    "raw": "identity",
    "off": "off",
    "skip": "off",
    "false": "off",
    "0": "off",
}

EXPORT_FLAG_GROUPS = (
    "export",
    "abundance_corrector",
    "abundance_correction",
    "corrector",
    "correction",
    "annotation_export",
)

TRUE_ALIASES = (
    "true",
    "taxID_true",
    "taxid_true",
    "taxID_TRUE",
    "truth",
)
UNCLASSIFIED = {"0", "", "nan", "none", "null", "<na>"}

DEFAULT_MIN_EFFICIENCY = 0.05
DEFAULT_MAX_FACTOR = 4.0
DEFAULT_C = 1.0
DEFAULT_MAX_ITER = 500


class MissingExportError(ValueError):
    """``tools.<name>`` is missing or is not ``--type export``."""


def is_skipped_export(name: Optional[str]) -> bool:
    key = str(name or "").strip().lower().replace("-", "_")
    return _BUILTIN_ALIASES.get(key, key) == "off"


def is_identity_export(name: Optional[str]) -> bool:
    key = str(name or "").strip().lower().replace("-", "_")
    return _BUILTIN_ALIASES.get(key, key) == "identity"


def resolve_export(name: Optional[str]) -> Tuple[str, str]:
    key = str(name or "logistic").strip() or "logistic"
    low = key.lower().replace("-", "_")
    if low in _BUILTIN_ALIASES:
        return "builtin", _BUILTIN_ALIASES[low]
    return "custom", key


def require_known_export(name: Optional[str]) -> str:
    kind, canon = resolve_export(name)
    if kind == "custom":
        lookup_export(canon)
        return canon
    return canon


def flags_apply_to_export(target: str, name: Optional[str] = None) -> bool:
    try:
        kind, canon = resolve_export(name)
    except MissingExportError:
        kind, canon = "custom", name
    names = [canon, name]
    if kind == "builtin" and canon == "logistic":
        names.extend(["logit", "logreg", "default"])
    if kind == "builtin" and canon == "identity":
        names.extend(["none", "raw"])
    return flags_target_matches(target, *names, groups=EXPORT_FLAG_GROUPS)


def is_export_flag_target(target: str) -> bool:
    extra = list(_BUILTIN_ALIASES.keys())
    extra.extend(iter_custom_export_names())
    return flags_target_matches(target, *extra, groups=EXPORT_FLAG_GROUPS)


def iter_custom_export_names(tools: Optional[Dict[str, Any]] = None) -> List[str]:
    mapping = tools if tools is not None else iter_tools(load_config())
    names: List[str] = []
    for name, spec in mapping.items():
        parsed = parse_tool_entry(spec, name)
        if str(parsed[3] or "").strip() != EXPORT_GROUP:
            continue
        if name.lower().replace("-", "_") in _BUILTIN_ALIASES:
            continue
        if tool_path(parsed, name):
            names.append(name)
    return names


def lookup_export(name: str) -> list:
    key = str(name or "").strip()
    if not key:
        raise MissingExportError(
            "Empty export name. Import with "
            "`samovar tools import -n NAME --type export` or use logistic|none."
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
        raise MissingExportError(
            f"export {key!r} is not in the install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type export`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != EXPORT_GROUP:
        raise MissingExportError(
            f"tools.{matched} has group {group!r}, expected {EXPORT_GROUP!r}. "
            "Re-import with --type export."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingExportError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return parsed


def attach_export_flags(name: str, config: Dict[str, Any]) -> Dict[str, Any]:
    cfg = dict(config or {})
    kind, canon = resolve_export(name)
    imported = ""
    if kind == "custom":
        spec = lookup_export(canon)
        imported = tool_flags(spec, canon)
    named_map = cfg.get("export_tool_flags") or {}
    named = named_map.get(canon) or named_map.get(name)
    if named is None:
        low = str(canon).lower()
        for key, value in named_map.items():
            if str(key).lower() == low:
                named = value
                break
    cfg["extra_flags"] = merge_flag_strings(
        imported,
        cfg.get("export_flags"),
        named,
        cfg.get("extra_flags"),
    )
    cfg["extra_argv"] = extra_flags_argv(cfg["extra_flags"])
    _apply_export_flag_tokens(cfg)
    cfg["export"] = canon
    return cfg


def _apply_export_flag_tokens(cfg: Dict[str, Any]) -> None:
    argv = list(cfg.get("extra_argv") or [])
    i = 0
    while i < len(argv):
        tok = str(argv[i])
        if tok in {"--min-efficiency", "--min_efficiency", "--floor"} and i + 1 < len(argv):
            cfg["min_efficiency"] = float(argv[i + 1])
            i += 2
            continue
        if tok in {"--max-factor", "--max_factor"} and i + 1 < len(argv):
            cfg["max_factor"] = float(argv[i + 1])
            i += 2
            continue
        if tok in {"--C", "-C"} and i + 1 < len(argv):
            cfg["C"] = float(argv[i + 1])
            i += 2
            continue
        if tok in {"--max-iter", "--max_iter"} and i + 1 < len(argv):
            cfg["max_iter"] = int(argv[i + 1])
            i += 2
            continue
        if tok in {"--n-reads", "--N_reads", "--n_reads"} and i + 1 < len(argv):
            cfg["n_reads"] = int(argv[i + 1])
            i += 2
            continue
        i += 1


def _taxid_str(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "0"
    text = str(value).strip()
    if text.lower() in UNCLASSIFIED:
        return "0"
    try:
        return str(int(float(text)))
    except (TypeError, ValueError):
        return text or "0"


def true_column(frame: pd.DataFrame) -> Optional[str]:
    if frame is None or frame.empty:
        return None
    cols = list(frame.columns)
    for name in TRUE_ALIASES:
        if name in cols:
            return name
    lower = {str(c).lower(): c for c in cols}
    for alias in ("true", "taxid_true", "truth"):
        if alias in lower:
            return lower[alias]
    return None


def normalize_export_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if frame is None or frame.empty:
        return pd.DataFrame()
    work = frame.copy()
    if "read" in work.columns and "seq" not in work.columns:
        work = work.rename(columns={"read": "seq"})
    col = true_column(work)
    if col and col != "true":
        work = work.rename(columns={col: "true"})
    return work


def annotation_frame(data: Any) -> pd.DataFrame:
    if data is None:
        return pd.DataFrame()
    if isinstance(data, Annotation):
        frame = data.DataFrame
        if frame is None or frame.empty:
            return pd.DataFrame()
        work = normalize_export_frame(frame)
        if "true" not in work.columns and getattr(data, "true_annotation", None):
            n = min(len(data.true_annotation), len(work))
            if n:
                work = work.copy()
                work["true"] = list(data.true_annotation)[:n] + [""] * (len(work) - n)
        return work
    if isinstance(data, pd.DataFrame):
        return normalize_export_frame(data)
    return pd.DataFrame()


def as_annotation(data: Any) -> Annotation:
    if isinstance(data, Annotation):
        frame = annotation_frame(data)
        if not frame.empty:
            obj = Annotation.from_long_table(frame)
            obj.abundance_tables = getattr(data, "abundance_tables", None)
            return obj
        return data
    if isinstance(data, pd.DataFrame):
        return Annotation.from_long_table(normalize_export_frame(data))
    if data in (None, "", False):
        return Annotation.from_long_table(pd.DataFrame())
    loaded = load_table_input(data)
    frame = annotation_frame(loaded)
    if not frame.empty:
        obj = Annotation.from_long_table(frame)
        obj.abundance_tables = getattr(loaded, "abundance_tables", None)
        return obj
    return loaded


def _logistic_efficiency(
    true: pd.Series,
    pred: pd.Series,
    *,
    C: float,
    max_iter: int,
    min_efficiency: float,
) -> Dict[str, float]:
    """P(pred == true | true = taxon) via one-hot logistic, else empirical recall."""
    true_s = true.map(_taxid_str)
    pred_s = pred.map(_taxid_str)
    mask = true_s != "0"
    if int(mask.sum()) < 2:
        return {}
    y = (pred_s[mask] == true_s[mask]).astype(int)
    taxa = true_s[mask]
    empirical = y.groupby(taxa).mean().to_dict()
    counts = taxa.value_counts().to_dict()
    rates: Dict[str, float] = {}
    for taxon, rec in empirical.items():
        if taxon == "0":
            continue
        rates[str(taxon)] = float(np.clip(rec, min_efficiency, 1.0))
    if int(y.nunique()) < 2:
        if rates:
            rates["__default__"] = float(np.clip(float(y.mean()), min_efficiency, 1.0))
        return rates
    labels = taxa.astype(str).values.reshape(-1, 1)
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import OneHotEncoder

        try:
            enc = OneHotEncoder(handle_unknown="ignore", sparse_output=True)
        except TypeError:
            enc = OneHotEncoder(handle_unknown="ignore", sparse=True)
        X = enc.fit_transform(labels)
        model = LogisticRegression(
            C=float(C),
            max_iter=int(max_iter),
            solver="lbfgs",
        )
        model.fit(X, y.to_numpy())
        unique = [t for t in enc.categories_[0] if str(t) != "0"]
        if unique:
            probs = model.predict_proba(enc.transform(np.array(unique, dtype=object).reshape(-1, 1)))
            pos = 1 if model.classes_[0] == 0 and len(model.classes_) > 1 else int(np.where(model.classes_ == 1)[0][0]) if 1 in model.classes_ else 0
            for taxon, row in zip(unique, probs):
                rates[str(taxon)] = float(np.clip(row[pos], min_efficiency, 1.0))
    except Exception:
        pass
    if rates:
        weights = [counts.get(t, 1) for t in rates]
        vals = list(rates.values())
        rates["__default__"] = float(
            np.clip(np.average(vals, weights=weights) if weights else np.mean(vals), min_efficiency, 1.0)
        )
    return rates


def fit_taxon_efficiency(
    reference: Any,
    *,
    C: float = DEFAULT_C,
    max_iter: int = DEFAULT_MAX_ITER,
    min_efficiency: float = DEFAULT_MIN_EFFICIENCY,
) -> Dict[str, Dict[str, float]]:
    """``{annotator: {taxid: recall}}`` from a reference Annotation with ``true``."""
    frame = annotation_frame(reference)
    if frame.empty or "true" not in frame.columns:
        return {}
    true = frame["true"]
    out: Dict[str, Dict[str, float]] = {}
    for col in taxid_annotation_columns(frame):
        name = annotator_name(col)
        if name.lower() in {"true", "samovar_confidence"}:
            continue
        rates = _logistic_efficiency(
            true,
            frame[col],
            C=C,
            max_iter=max_iter,
            min_efficiency=min_efficiency,
        )
        if rates:
            out[name] = rates
    return out


def apply_taxon_efficiency(
    tables: Dict[str, pd.DataFrame],
    efficiencies: Dict[str, Dict[str, float]],
    *,
    min_efficiency: float = DEFAULT_MIN_EFFICIENCY,
) -> Dict[str, pd.DataFrame]:
    """Divide predicted counts by per-taxon recall (``count / efficiency``)."""
    if not efficiencies:
        return tables
    corrected: Dict[str, pd.DataFrame] = {}
    for name, table in tables.items():
        work = normalize_abundance_table(table)
        rates = efficiencies.get(name) or efficiencies.get(name.lower()) or {}
        if not rates:
            low = str(name).lower()
            for key, value in efficiencies.items():
                if str(key).lower() == low:
                    rates = value
                    break
        if not rates:
            corrected[name] = work
            continue
        cols = n_sample_columns(work)
        default = float(rates.get("__default__", 1.0))
        scale = []
        for tax in work["taxid"].map(_taxid_str):
            if tax == "0":
                scale.append(1.0)
                continue
            eff = float(rates.get(tax, default))
            scale.append(max(eff, min_efficiency))
        denom = np.array(scale, dtype=float)
        denom[denom <= 0] = min_efficiency
        for col in cols:
            work[col] = pd.to_numeric(work[col], errors="coerce").fillna(0.0) / denom
        corrected[name] = work.reset_index(drop=True)
    return corrected


def maybe_attach_samovar_reference(reference: Annotation, model_path: Optional[PathLike]) -> Annotation:
    """Score regenerated reads with a trained reprofiler so SAMOVAR can be corrected."""
    path = Path(str(model_path or "")).expanduser() if model_path else None
    if path is None or not path.is_file():
        return reference
    frame = annotation_frame(reference)
    if frame.empty or "true" not in frame.columns:
        return reference
    if any(annotator_name(c).lower() == "samovar" for c in taxid_annotation_columns(frame)):
        return reference
    try:
        from samovar.reprofiling import predict_taxid
    except Exception:
        return reference
    try:
        scored = predict_taxid(frame, model_path=str(path))
    except Exception:
        return reference
    if "taxid_SAMOVAR" not in scored.columns:
        return reference
    merged = frame.copy()
    merged["taxid_SAMOVAR"] = scored["taxid_SAMOVAR"].values[: len(merged)]
    obj = Annotation.from_long_table(merged)
    obj.abundance_tables = getattr(reference, "abundance_tables", None)
    return obj


def write_export_formats(
    annotation: Annotation,
    dest: PathLike,
    formats: Sequence[str],
    **kwargs: Any,
) -> Path:
    names = parse_export_formats(*list(formats)) or ["abundance"]
    dest_path = Path(dest)
    if len(names) == 1:
        return dump_builtin(annotation, dest_path, names[0], **kwargs)
    dest_path.mkdir(parents=True, exist_ok=True)
    last = dest_path
    for fmt in names:
        last = dump_builtin(annotation, dest_path / fmt, fmt, **kwargs)
    return last


def export_identity(annotation: Annotation, dest: PathLike, config: Optional[Dict[str, Any]] = None) -> Path:
    cfg = dict(config or {})
    tables = input_to_abundance_tables(annotation)
    tables = {k: normalize_abundance_table(v) for k, v in tables.items()}
    annotation.abundance_tables = tables
    formats = parse_export_formats(cfg.get("to") or cfg.get("formats") or "abundance")
    return write_export_formats(
        annotation,
        dest,
        formats,
        n_reads=cfg.get("n_reads"),
        taxdump=cfg.get("taxdump"),
        taxonomy=cfg.get("taxonomy") or "ncbi",
        extra_flags=cfg.get("extra_flags") or "",
        extra_argv=cfg.get("extra_argv") or [],
    )


def export_logistic(annotation: Annotation, dest: PathLike, config: Optional[Dict[str, Any]] = None) -> Path:
    cfg = attach_export_flags("logistic", dict(config or {}))
    min_eff = float(cfg.get("min_efficiency", DEFAULT_MIN_EFFICIENCY))
    reference = cfg.get("reference")
    model_path = cfg.get("model") or cfg.get("model_path")
    if reference is not None:
        ref = as_annotation(reference)
        ref = maybe_attach_samovar_reference(ref, model_path)
        rates = fit_taxon_efficiency(
            ref,
            C=float(cfg.get("C", DEFAULT_C)),
            max_iter=int(cfg.get("max_iter", DEFAULT_MAX_ITER)),
            min_efficiency=min_eff,
        )
    else:
        rates = {}
    tables = input_to_abundance_tables(annotation)
    tables = apply_taxon_efficiency(tables, rates, min_efficiency=min_eff)
    annotation.abundance_tables = tables
    formats = parse_export_formats(cfg.get("to") or cfg.get("formats") or "abundance")
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if rates:
        meta = dest_path / "taxon_efficiency.json" if dest_path.suffix else dest_path / "taxon_efficiency.json"
        if dest_path.suffix.lower() in {".csv", ".tsv", ".txt"}:
            meta = dest_path.with_name("taxon_efficiency.json")
        try:
            import json

            meta.parent.mkdir(parents=True, exist_ok=True)
            meta.write_text(json.dumps(rates, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        except OSError:
            pass
    return write_export_formats(
        annotation,
        dest,
        formats,
        n_reads=cfg.get("n_reads"),
        taxdump=cfg.get("taxdump"),
        taxonomy=cfg.get("taxonomy") or "ncbi",
        extra_flags=cfg.get("extra_flags") or "",
        extra_argv=cfg.get("extra_argv") or [],
    )


def _load_python_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(f"samovar_export_{name}", path)
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _custom_export(annotation: Annotation, dest: Path, config: Dict[str, Any], spec: list, name: str) -> Path:
    path = Path(tool_path(spec, name))
    if path.suffix.lower() == ".py" and path.is_file():
        module = _load_python_module(path, name)
        fn = getattr(module, "export", None) if module is not None else None
        cls = getattr(module, "Exporter", None) or getattr(module, "AbundanceCorrector", None) if module is not None else None
        if not callable(fn) and cls is not None:
            try:
                inst = cls()
            except Exception:
                inst = None
            if inst is not None:
                fn = getattr(inst, "export", None) or getattr(inst, "run", None)
        if callable(fn):
            written = fn(annotation, dest, config)
            if isinstance(written, dict) and written:
                first = next(iter(written.values()))
                if isinstance(first, pd.DataFrame):
                    annotation.abundance_tables = {
                        str(k): normalize_abundance_table(v) for k, v in written.items()
                    }
                    return write_export_formats(
                        annotation,
                        dest,
                        parse_export_formats(config.get("to") or "abundance"),
                        n_reads=config.get("n_reads"),
                        taxdump=config.get("taxdump"),
                        taxonomy=config.get("taxonomy") or "ncbi",
                    )
            return Path(written) if written else dest
    dest.mkdir(parents=True, exist_ok=True)
    tmp_in = dest / ".annotation.csv"
    frame = annotation_frame(annotation)
    if not frame.empty:
        frame.to_csv(tmp_in, index=False)
    cmd = [str(path), "-i", str(tmp_in if tmp_in.is_file() else dest), "-o", str(dest)]
    reference = config.get("reference")
    if reference is not None:
        ref_path = Path(str(getattr(reference, "path", "") or ""))
        if not ref_path.is_file() and not ref_path.is_dir():
            ref_dir = dest / ".reference"
            ref_ann = as_annotation(reference)
            write_abundance_dir(ref_dir, input_to_abundance_tables(ref_ann) or {})
            frame = annotation_frame(ref_ann)
            if not frame.empty:
                csv = ref_dir / "combined_annotation_table.csv"
                csv.parent.mkdir(parents=True, exist_ok=True)
                frame.to_csv(csv, index=False)
                cmd.extend(["-r", str(csv)])
        else:
            cmd.extend(["-r", str(ref_path)])
    extra = list(config.get("extra_argv") or extra_flags_argv(config.get("extra_flags")))
    cmd.extend(extra)
    subprocess.check_call(cmd)
    return dest


def run_export(
    name: Optional[str],
    annotation: Any,
    dest: PathLike,
    config: Optional[Dict[str, Any]] = None,
) -> Optional[Path]:
    kind, canon = resolve_export(name)
    if canon == "off":
        return None
    cfg = attach_export_flags(canon if kind == "custom" else name or canon, dict(config or {}))
    ann = as_annotation(annotation)
    dest_path = Path(dest)
    if kind == "builtin" and canon == "identity":
        return export_identity(ann, dest_path, cfg)
    if kind == "builtin" and canon == "logistic":
        return export_logistic(ann, dest_path, cfg)
    spec = lookup_export(canon)
    return _custom_export(ann, dest_path, cfg, spec, canon)


def export_run_dir(
    output_dir: PathLike,
    *,
    stage: str,
    source: PathLike,
    dest: Optional[PathLike] = None,
    config: Optional[Dict[str, Any]] = None,
) -> Optional[Path]:
    """Export one annotation directory from a SAMOVAR run."""
    cfg = dict(config or {})
    root = Path(output_dir)
    src = Path(source)
    if not src.is_absolute():
        src = root / src
    out = Path(dest) if dest else root / "exports" / stage
    name = cfg.get("export") or cfg.get("corrector") or "logistic"
    reference = cfg.get("reference")
    if reference in (None, "", False):
        regen = root / "regenerated_annotations"
        if regen.is_dir():
            reference = regen
    model = cfg.get("model") or cfg.get("model_path")
    if model in (None, "", False):
        joblib = root / "reprofiled_annotations" / "trained_model.joblib"
        if joblib.is_file():
            model = joblib
    cfg["reference"] = reference
    cfg["model"] = model
    cfg["stage"] = stage
    return run_export(name, src, out, cfg)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="samovar export-abundance",
        description=(
            "Annotation → abundance-like tables (abundance, kraken2, cami). "
            "Default method: logistic per-taxon recall correction. "
            "Custom tools: samovar tools import --type export."
        ),
    )
    parser.add_argument("-i", "--input", dest="source", required=True, help="Annotation dir or table")
    parser.add_argument("-o", "--output", dest="dest", required=True, help="Destination dir")
    parser.add_argument(
        "-r",
        "--reference",
        default="",
        help="Reference Annotation with taxID_true / true (usually regenerated_annotations)",
    )
    parser.add_argument(
        "--export",
        "--corrector",
        dest="export",
        default="logistic",
        help="logistic (default), none/identity, or an imported --type export name",
    )
    parser.add_argument(
        "--to",
        dest="formats",
        action="append",
        default=None,
        help="Output format (repeatable): abundance, kraken2, cami, …",
    )
    parser.add_argument("--model", dest="model", default="", help="trained_model.joblib for SAMOVAR correction")
    parser.add_argument("--taxdump", default="", help="Taxdump for kraken2 / cami")
    parser.add_argument("--taxonomy", default="ncbi", help="ncbi or gtdb")
    parser.add_argument("--config", dest="config_yaml", default="", help="YAML from prepare (export flags)")
    parser.add_argument("--output_dir", dest="output_dir", default="", help="Run directory (fills reference/model)")
    parser.add_argument("--stage", default="", help="Export stage label (initial/regenerated/reprofiled)")
    parser.add_argument("--n-reads", dest="n_reads", type=int, default=None)
    parser.add_argument(
        "--flags",
        default="",
        help='Extra flags, e.g. "--min-efficiency 0.05 --C 1"',
    )
    return parser


def _load_yaml(path: str) -> Dict[str, Any]:
    if not path:
        return {}
    target = Path(path)
    if not target.is_file():
        return {}
    import yaml

    data = yaml.safe_load(target.read_text(encoding="utf-8")) or {}
    return data if isinstance(data, dict) else {}


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(None if argv is None else list(argv))
    yaml_cfg = _load_yaml(str(args.config_yaml or ""))
    name = args.export or yaml_cfg.get("export") or yaml_cfg.get("corrector") or "logistic"
    formats = parse_export_formats(args.formats, yaml_cfg.get("export_formats"), yaml_cfg.get("to"))
    cfg: Dict[str, Any] = dict(yaml_cfg)
    cfg["export"] = name
    cfg["to"] = formats or ["abundance"]
    cfg["export_flags"] = merge_flag_strings(cfg.get("export_flags"), args.flags)
    cfg["export_tool_flags"] = cfg.get("export_tool_flags") or {}
    if args.n_reads is not None:
        cfg["n_reads"] = args.n_reads
    if args.taxdump:
        cfg["taxdump"] = args.taxdump
    elif yaml_cfg.get("export_taxdump"):
        cfg["taxdump"] = yaml_cfg.get("export_taxdump")
    if args.taxonomy:
        cfg["taxonomy"] = args.taxonomy
    elif yaml_cfg.get("export_taxonomy"):
        cfg["taxonomy"] = yaml_cfg.get("export_taxonomy")
    if args.reference:
        cfg["reference"] = args.reference
    if args.model:
        cfg["model"] = args.model
    root = args.output_dir or yaml_cfg.get("output_dir") or ""
    if root:
        written = export_run_dir(
            root,
            stage=args.stage or "export",
            source=args.source,
            dest=args.dest,
            config=cfg,
        )
    else:
        written = run_export(name, args.source, args.dest, cfg)
    if written is None:
        return 0
    print(written)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
