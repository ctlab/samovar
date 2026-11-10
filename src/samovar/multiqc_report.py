"""Which MultiQC sections a finished run should rebuild.

Section choice is separate from ``viz`` (per-read metrics and sampling) and from
visual style. ``samovar multiqc`` reads saved artifacts only.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Union

SECTIONS = (
    "qc",
    "benchmark",
    "correction",
    "mag",
    "profiles",
    "composition",
    "alpha",
    "beta",
    "runtime",
    "feature_importance",
    "table_scores",
)
STAGES = ("initial", "regenerated", "reprofiled")
FIELDS = ("sections", "include", "exclude", "stages", "annotators")
DEFAULT_SETTINGS = dict(sections="all", include=[], exclude=[], stages="all", annotators="all")
_ALWAYS = {"run_options", "conclusion", "intro", "status"}
_CATALOG = {
    "sections": SECTIONS,
    "include": SECTIONS,
    "exclude": SECTIONS,
    "stages": STAGES,
}

_MISSING = {
    "qc": "no fastp or FastQC report under initial_trimmed or regenerated_trimmed",
    "benchmark": "no score table or plot for the selected stages",
    "correction": "no profile_ratio.json, profile_joint.json, or trained_model.joblib",
    "mag": "no MAG mapping sidecar or MAG_ID column",
    "profiles": "no saved profile table",
    "composition": "no composition plot in this run",
    "alpha": "no alpha-diversity table in this run",
    "beta": "no beta-diversity table in this run",
    "runtime": "no recorded stage or per-sample runtime in this run",
    "feature_importance": "no feature-importance plot",
    "table_scores": "no table_selection.json",
}


def _split_tokens(value: Any) -> List[str]:
    if isinstance(value, str):
        return [part for part in re.split(r"[,\s]+", value.strip()) if part]
    if isinstance(value, (list, tuple)):
        tokens: List[str] = []
        for item in value:
            tokens.extend(_split_tokens(item) if not isinstance(item, str) else _split_tokens(str(item)))
        return tokens
    raise ValueError("multiqc lists must be a string or a list of names")


def _field_value(value: Any, key: str) -> Union[str, List[str]]:
    if value is None or value == "all":
        return "all"
    tokens = _split_tokens(value)
    if tokens == ["all"]:
        return "all"
    catalog = _CATALOG.get(key)
    if catalog is not None:
        unknown = [name for name in tokens if name not in catalog]
        if unknown:
            raise ValueError(
                f"Unknown multiqc {key}: {', '.join(unknown)}. "
                f"Known names: {', '.join(catalog)}"
            )
    return tokens


def parse_multiqc_block(block: Any) -> Dict[str, Any]:
    """Validate one config layer. Omitted fields are left unset."""
    if block is None:
        return {}
    if not isinstance(block, Mapping):
        raise ValueError("multiqc settings must be a mapping")
    unknown = sorted(set(block) - set(FIELDS))
    if unknown:
        raise ValueError(
            f"Unknown multiqc setting(s): {', '.join(unknown)}. "
            f"Known settings: {', '.join(FIELDS)}"
        )
    return {key: _field_value(block[key], key) for key in FIELDS if key in block}


def effective_settings(layers: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """Later layers win. ``include`` adds and ``exclude`` removes sections."""
    merged: Dict[str, Any] = dict(DEFAULT_SETTINGS)
    for layer in layers:
        merged.update(parse_multiqc_block(layer))
    sections = set(SECTIONS) if merged["sections"] == "all" else set(merged["sections"])
    include = [] if merged["include"] == "all" else list(merged["include"])
    if merged["include"] == "all":
        sections = set(SECTIONS)
    sections.update(include)
    exclude = [] if merged["exclude"] == "all" else list(merged["exclude"])
    if merged["exclude"] == "all":
        sections.clear()
    for name in exclude:
        sections.discard(name)
    stages = set(STAGES) if merged["stages"] == "all" else set(merged["stages"])
    annotators = "all" if merged["annotators"] == "all" else [str(name) for name in merged["annotators"]]
    return {
        "sections": sections,
        "include": include,
        "exclude": exclude,
        "stages": stages,
        "annotators": annotators,
    }


def _yaml_mapping(path: Path) -> Mapping[str, Any]:
    if not path.is_file():
        return {}
    import yaml

    loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(loaded, Mapping):
        raise ValueError(f"{path} must contain a mapping")
    return loaded


def load_run_block(path: Path) -> Dict[str, Any]:
    """Read only the ``multiqc`` key from a run snapshot."""
    loaded = _yaml_mapping(path)
    if "multiqc" not in loaded:
        return {}
    return parse_multiqc_block(loaded["multiqc"])


def load_config_file(path: Path) -> Dict[str, Any]:
    """A ``--config`` file is either a snapshot with a ``multiqc`` key or the block itself."""
    loaded = _yaml_mapping(path)
    if "multiqc" in loaded:
        return parse_multiqc_block(loaded["multiqc"])
    if not (set(loaded) & set(FIELDS)):
        return {}
    return parse_multiqc_block(loaded)


def resolve_settings(
    root: Path,
    cli: Optional[Mapping[str, Any]] = None,
    config_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Defaults, then the run's ``.hydra/config.yaml``, then ``--config``, then CLI."""
    layers: List[Mapping[str, Any]] = []
    hydra = Path(root) / ".hydra" / "config.yaml"
    if hydra.is_file():
        layers.append(load_run_block(hydra))
    if config_path is not None:
        layers.append(load_config_file(Path(config_path)))
    if cli:
        layers.append({key: cli[key] for key in FIELDS if key in cli and cli[key] is not None})
    return effective_settings(layers)


def parse_cli(tokens: Sequence[str]) -> tuple:
    """Split report flags from the MultiQC passthrough that follows ``--``."""
    flags = {
        "--sections": "sections",
        "--include": "include",
        "--exclude": "exclude",
        "--stages": "stages",
        "--annotators": "annotators",
        "--config": "config",
    }
    cli: Dict[str, str] = {}
    extra: List[str] = []
    index = 0
    items = list(tokens)
    while index < len(items):
        token = items[index]
        if token == "--":
            extra = items[index + 1 :]
            break
        name, value = token, None
        if token.startswith("--") and "=" in token:
            name, value = token.split("=", 1)
        if name in flags:
            if value is None:
                if index + 1 >= len(items):
                    raise ValueError(f"{name} needs a value")
                value = items[index + 1]
                index += 2
            else:
                index += 1
            cli[flags[name]] = value
            continue
        if token.startswith("-"):
            raise ValueError(f"Unknown multiqc option {token}. Pass MultiQC flags after --")
        raise ValueError(f"Unexpected multiqc argument {token}")
    return cli, extra


def section_of(filename: str) -> str:
    """Which report section a staged file belongs to."""
    name = filename.lower()
    if "run_options" in name:
        return "run_options"
    if "conclusion" in name:
        return "conclusion"
    if "report_status" in name or name.startswith("status_"):
        return "status"
    if name.endswith("_intro_mqc.html"):
        return "intro"
    if "mag_" in name:
        return "mag"
    if "tablescore" in name or "table_score" in name:
        return "table_scores"
    if "composition" in name:
        return "composition"
    if "profile_inventory" in name or "profile_quality" in name or name.startswith("profile"):
        return "profiles"
    if "feature_importance" in name or "importance" in name:
        return "feature_importance"
    if "correct" in name:
        return "correction"
    if "runtime" in name:
        return "runtime"
    if re.search(r"(^|_)alpha(_|$)|shannon|richness", name):
        return "alpha"
    if re.search(r"(^|_)beta(_|$)|pcoa|mantel", name):
        return "beta"
    if "fastp" in name or "fastqc" in name:
        return "qc"
    return "benchmark"


def _matches(label: str, annotators: Sequence[str]) -> bool:
    text = str(label).lower()
    return any(token.lower() in text for token in annotators)


def profile_selected(name, info, annotators):
    """Select complete profiles before numerical calculations, not plotted points."""
    if annotators == "all":
        return True
    labels = [name, *(info.get(field) or "" for field in ("tool", "run_name", "database", "source_identity"))]
    return _matches(" ".join(map(str, labels)), annotators)


def _filter_payload(payload: Dict[str, Any], annotators: Sequence[str]) -> Optional[Dict[str, Any]]:
    if not annotators:
        return None
    title = str(payload.get("section_name") or "")
    categories = payload.get("categories")
    if isinstance(categories, dict):
        kept = {
            key: value
            for key, value in categories.items()
            if _matches(key, annotators)
            or _matches(str((value or {}).get("name", "")), annotators)
        }
        if not kept and not _matches(title, annotators):
            return None
        if kept:
            payload["categories"] = kept
            data = payload.get("data")
            if isinstance(data, dict):
                payload["data"] = {
                    metric: (
                        {key: value for key, value in series.items() if key in kept or _matches(key, annotators)}
                        if isinstance(series, dict)
                        else series
                    )
                    for metric, series in data.items()
                }
        return payload
    data = payload.get("data")
    labels = (payload.get("pconfig") or {}).get("data_labels")
    if isinstance(data, list) and isinstance(labels, list) and len(labels) == len(data):
        kept_data = []
        kept_labels = []
        for series, label in zip(data, labels):
            if _matches(str(label), annotators):
                kept_data.append(series)
                kept_labels.append(label)
        if not kept_data:
            if title and not _matches(title, annotators):
                return None
            return payload
        payload["data"] = kept_data
        payload["pconfig"]["data_labels"] = kept_labels
        return payload
    if isinstance(data, dict) and data and all(isinstance(value, dict) for value in data.values()):
        if all(str(key).isdigit() for key in data):
            rows = {}
            for row in data.values():
                blob = " ".join(
                    str(row.get(field, ""))
                    for field in ("annotator", "tool", "run_name", "profile", "database")
                )
                if not blob.strip() or _matches(blob, annotators):
                    rows[str(len(rows))] = row
            if not rows:
                return None
            payload["data"] = rows
            return payload
        kept = {key: value for key, value in data.items() if _matches(str(key), annotators)}
        if not kept:
            if title and not _matches(title, annotators):
                return None
            return payload
        payload["data"] = kept
        return payload
    if title and not _matches(title, annotators) and payload.get("plot_type") == "heatmap":
        return None
    return payload


def apply_selection(dest: Path, settings: Mapping[str, Any]) -> set:
    """Drop staged files for sections that are off, and annotators that were not asked for."""
    enabled = set(settings["sections"])
    annotators = settings["annotators"]
    produced = set()
    for path in list(dest.iterdir()):
        if not path.is_file():
            continue
        section = section_of(path.name)
        if section in _ALWAYS:
            continue
        if section not in enabled:
            path.unlink()
            continue
        if annotators != "all" and path.suffix == ".json":
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                produced.add(section)
                continue
            if isinstance(payload, dict):
                filtered = _filter_payload(payload, annotators)
                if filtered is None:
                    path.unlink()
                    continue
                path.write_text(json.dumps(filtered, indent=2) + "\n", encoding="utf-8")
        produced.add(section)
    qc = dest / "qc"
    if qc.is_dir() and any(qc.iterdir()):
        if "qc" in enabled:
            produced.add("qc")
        else:
            for path in qc.iterdir():
                path.unlink()
    return produced


def missing_sections(settings: Mapping[str, Any], produced: Iterable[str]) -> List[tuple]:
    have = set(produced)
    return [
        (name, _MISSING[name])
        for name in SECTIONS
        if name in settings["sections"] and name not in have
    ]


def qc_sources(root: Path) -> List[Path]:
    found = []
    for folder in ("initial_trimmed", "regenerated_trimmed"):
        base = Path(root) / folder
        if not base.is_dir():
            continue
        for pattern in ("*.fastp.json", "*.fastp.html", "*_fastqc.html", "*_fastqc.zip"):
            found.extend(path for path in base.glob(pattern) if path.is_file())
    return found


def correction_sources(root: Path) -> List[str]:
    relative = (
        "reprofiled_profiles/profile_ratio.json",
        "reprofiled_profiles/profile_joint.json",
        "reprofiled_annotations/trained_model.joblib",
    )
    return [name for name in relative if (Path(root) / name).is_file()]


def settings_document(settings: Mapping[str, Any]) -> Dict[str, Any]:
    annotators = settings["annotators"]
    return {
        "sections": sorted(settings["sections"]),
        "include": list(settings["include"]),
        "exclude": list(settings["exclude"]),
        "stages": sorted(settings["stages"]),
        "annotators": "all" if annotators == "all" else list(annotators),
        "priority": "CLI overrides --config, which overrides .hydra/config.yaml, which overrides the built-in maximum",
    }
