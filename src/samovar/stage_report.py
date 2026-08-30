"""FastQC/fastp-style JSON summaries plus MultiQC custom-content staging.

Each exec checkpoint writes ``<output>/.log/multiqc/<stage>.samovar.json``.
``bundle_multiqc`` copies native Plotly heatmaps/scatters/bars and score tables
into ``<output>/multiqc_samovar/`` as ``*_mqc.*`` files that stock MultiQC can render.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

from samovar.exec_control import CHECKPOINT_STEPS
from samovar.paths import PACKAGE_VERSION
from samovar.scores import (
    ENSEMBLE_NAME,
    OPAL_DISPLAY,
    SCORE_DISPLAY,
    canonical_annotator_name,
    is_platform_read_type,
)

PathLike = Union[str, os.PathLike]

STAGE_INFO: Dict[str, Dict[str, str]] = {
    "setup_reads": {
        "title": "Read setup",
        "description": (
            "This is the first checkpoint: SAMOVAR puts the metagenome it will "
            "annotate into ``initial/`` as paired FASTQ. That may be real sequencing "
            "data you already had, or InSilicoSeq output from ``samovar generate``. "
            "Later stages never go back to the original genome folder; they only "
            "see these reads. Docs: https://github.com/ctlab/samovar/wiki"
        ),
    },
    "qc_initial": {
        "title": "QC initial reads",
        "description": (
            "Optional imported ``--type QC`` tool trims ``initial/`` into "
            "``initial_trimmed/``. When no QC is configured this step is identity "
            "(reads are already treated as trimmed). Hybrid files may use a "
            "different tool per postfix (``illumina``, ``ont``, …)."
        ),
    },
    "annotate_initial": {
        "title": "Initial annotation",
        "description": (
            "Each configured classifier (Kraken2, Kaiju, and any tools you imported) "
            "is run on the same ``initial_trimmed/`` FASTQ (identity of ``initial/`` "
            "when no QC tool is set). The raw per-tool reports "
            "(``.kreport``, Kaiju ``.out``, and so on) land in ``initial_reports/``. "
            "This is the “what did each database/tool call?” step before anything is merged."
        ),
    },
    "combine_initial": {
        "title": "Combine initial calls",
        "description": (
            "Tool outputs are joined into one long table per sample in "
            "``initial_annotations/``. Columns are ``taxID_<tool>``; when ISS or CAMISIM "
            "headers contain ``taxid:<digits>`` (or a known assembly accession), a "
            "``true`` column is filled so F1 and scores can be computed."
        ),
    },
    "viz_initial": {
        "title": "Raw",
        "description": (
            "Quality of the original community: F1, completeness, and related scores "
            "comparing each annotator (and later SAMOVAR) to the ``true`` taxIDs when "
            "those exist. This is the baseline before regeneration. If ``true`` is "
            "missing, heatmaps of calls remain but F1 is undefined."
        ),
    },
    "abundance_tables": {
        "title": "Observed abundance",
        "description": (
            "Per-read annotations are counted into OTU-style tables "
            "(``taxid`` + ``N_<sample>``) in ``initial_abundance/``. Table regenerators "
            "and ISS community rebuild start from these counts, not from the FASTQ."
        ),
    },
    "regenerate_tables": {
        "title": "Regenerated abundance",
        "description": (
            "One or more ``table_reads_generator`` methods (direct copy, bootstrap, "
            "GLM, SparseDOSSA2, …) write candidate synthetic communities under "
            "``regenerated/.regenerated_abundance/``. The goal is a table that is "
            "statistically close to the observed mix so re-annotation has known labels."
        ),
    },
    "score_regenerated_tables": {
        "title": "Table scoring",
        "description": (
            "When several candidate tables exist, ``--table-score`` (for example "
            "SparseDOSSA2 CV) ranks them against the observed community. The winner is "
            "stored in ``table_selection.json`` and is the only table used for ISS "
            "regeneration. This matters for custom annotators that emit tables rather than reads."
        ),
    },
    "seed_genomes": {
        "title": "Genome library",
        "description": (
            "Assemblies needed to simulate the regenerated community are resolved "
            "(NCBI, the SamovaR catalog, or ``--test-genomes``) into ``genomes/``. "
            "Without this step ISS has nothing to draw reads from."
        ),
    },
    "regenerate_reads": {
        "title": "Community regeneration",
        "description": (
            "ISS (or another configured generator) draws FASTQ from the chosen "
            "abundance table and the seeded genomes. Those reads have known composition, "
            "so the ensemble can train on labels that are not circular with the original sample."
        ),
    },
    "sort_reads": {
        "title": "Sort regenerated reads",
        "description": (
            "Regenerated FASTQ is ordered so re-annotation and the ML step see a "
            "stable read order across tools."
        ),
    },
    "qc_generated": {
        "title": "QC generated reads",
        "description": (
            "Same QC contract as ``qc_initial``, applied to regenerated FASTQ "
            "(``regenerated_trimmed/``). ``--qc-generated`` may differ from "
            "``--qc-initial``; hybrid postfix maps still apply."
        ),
    },
    "annotate_regenerated": {
        "title": "Re-annotation",
        "description": (
            "The same classifiers run on the regenerated FASTQ (``regenerated_reports/``). "
            "These calls, together with the known simulated mix, are the training set "
            "for reprofiling."
        ),
    },
    "combine_regenerated": {
        "title": "Combine re-annotation",
        "description": (
            "Regenerated tool calls are merged like the initial tables "
            "(``regenerated_annotations/``) so each read has a vote from every annotator "
            "plus a true taxID from the simulator."
        ),
    },
    "viz_regenerated": {
        "title": "Regenerated",
        "description": (
            "The same metrics as Raw, but on the in-silico community. Composition is "
            "known, so F1 here measures how the tools behave when the mix is controlled. "
            "Compare this section to Raw to see domain shift between real/ISS input and "
            "the rebuilt community."
        ),
    },
    "reprofile": {
        "title": "SAMOVAR re-profiling",
        "description": (
            "A supervised model (default ensemble, or an imported ``--type ml`` tool) "
            "is trained on regenerated labels and applied to the original samples. "
            "The corrected column is ``taxid_SAMOVAR`` in ``reprofiled_annotations/``."
        ),
    },
    "viz_reprofiled": {
        "title": "Reprofiled",
        "description": (
            "F1 and scores after the ensemble (or custom reprofiler) has updated the "
            "original samples. This is the main “did SAMOVAR help?” comparison against "
            "the single tools in Raw."
        ),
    },
}

STAGE_DIRS = {
    "setup_reads": ("initial",),
    "qc_initial": ("initial_trimmed",),
    "annotate_initial": ("initial_reports",),
    "combine_initial": ("initial_annotations",),
    "viz_initial": ("initial_annotations_plots",),
    "abundance_tables": ("initial_abundance",),
    "regenerate_tables": ("regenerated/.regenerated_abundance",),
    "score_regenerated_tables": ("regenerated/.regenerated_abundance",),
    "seed_genomes": ("genomes",),
    "regenerate_reads": ("regenerated",),
    "sort_reads": ("regenerated",),
    "qc_generated": ("regenerated_trimmed",),
    "annotate_regenerated": ("regenerated_reports",),
    "combine_regenerated": ("regenerated_annotations",),
    "viz_regenerated": ("regenerated_annotations_plots",),
    "reprofile": ("reprofiled_annotations",),
    "viz_reprofiled": ("reprofiled_annotations_plots",),
}

PLOT_SUFFIXES = {".png", ".html", ".svg", ".pdf"}
JSON_NAME = "{stage}.samovar.json"

REPORT_STAGES = ("viz_initial", "viz_regenerated", "viz_reprofiled")
REPORT_SECTION = {
    "viz_initial": {
        "id": "samovar_raw",
        "title": "Raw",
        "description": STAGE_INFO["viz_initial"]["description"],
        "folder": "initial_annotations_plots",
    },
    "viz_regenerated": {
        "id": "samovar_regenerated",
        "title": "Regenerated",
        "description": STAGE_INFO["viz_regenerated"]["description"],
        "folder": "regenerated_annotations_plots",
    },
    "viz_reprofiled": {
        "id": "samovar_reprofiled",
        "title": "Reprofiled",
        "description": STAGE_INFO["viz_reprofiled"]["description"],
        "folder": "reprofiled_annotations_plots",
    },
}

HEATMAP_COLSTOPS = [
    [0.0, "#FFFFFF"],
    [0.1, "#FFFFE5"],
    [0.2, "#F7FCB9"],
    [0.3, "#E6E487"],
    [0.4, "#D9F0A3"],
    [0.5, "#ADDD8E"],
    [0.6, "#78C679"],
    [0.7, "#41AB5D"],
    [0.8, "#238443"],
    [0.9, "#006837"],
    [1.0, "#004529"],
]

ANNOTATOR_COLORS = {
    "kaiju": "#1B9E77",
    "kraken2": "#D95F02",
    "kraken": "#E6AB02",
    "krakenuniq": "#A6761D",
    "centrifuge": "#66A61E",
    "metaphlan": "#1F78B4",
    "metaphlan4": "#1F78B4",
    ENSEMBLE_NAME: "#7570B3",
    "samovar": "#7570B3",
}
ANNOTATOR_COLOR_FALLBACK = (
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#1F78B4",
    "#A6761D",
    "#666666",
)

BAR_METRIC_ORDER = (
    "accuracy",
    "f1",
    "f1_macro",
    "accuracy_purity",
    "f1_purity",
    "completeness",
    "opal_purity",
    "opal_f1",
    "jaccard",
    "bray_curtis",
)

CONCLUSION_METRICS = (
    "f1",
    "accuracy",
    "f1_macro",
    "r2",
    "accuracy_purity",
    "f1_purity",
    "completeness",
    "opal_f1",
)


def as_path(path: PathLike) -> Path:
    return Path(path)


def report_dir(output_dir: PathLike) -> Path:
    return as_path(output_dir) / ".log" / "multiqc"


def multiqc_stage_dir(output_dir: PathLike) -> Path:
    return as_path(output_dir) / "multiqc_samovar"


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _count_files(root: Path, patterns: Sequence[str]) -> int:
    if not root.is_dir():
        return 0
    n = 0
    for pattern in patterns:
        n += sum(1 for p in root.rglob(pattern) if p.is_file())
    return n


def _list_rel(root: Path, suffixes: Iterable[str], limit: int = 80) -> List[str]:
    if not root.is_dir():
        return []
    want = {s.lower() for s in suffixes}
    hits = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.suffix.lower() in want:
            hits.append(str(path.relative_to(root)))
        if len(hits) >= limit:
            break
    return hits


def _safe_json(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_json(v) for v in value]
    try:
        import numpy as np

        if isinstance(value, np.generic):
            return value.item()
    except Exception:
        pass
    return str(value)


def _read_scores_table(plots: Path) -> List[Dict[str, Any]]:
    csv = plots / "quality_scores.csv"
    if not csv.is_file():
        return []
    try:
        import pandas as pd

        table = pd.read_csv(csv)
    except Exception:
        return []
    rows = []
    for _, row in table.iterrows():
        rec = {}
        for key, val in row.items():
            if hasattr(val, "item"):
                try:
                    val = val.item()
                except Exception:
                    pass
            if isinstance(val, float) and val != val:
                continue
            rec[str(key)] = val
        rows.append(rec)
    return rows


def collect_stage_payload(output_dir: PathLike, stage: str) -> Dict[str, Any]:
    root = as_path(output_dir)
    info = STAGE_INFO.get(stage, {"title": stage, "description": ""})
    summary: Dict[str, Any] = {}
    files: Dict[str, Any] = {}
    for name in STAGE_DIRS.get(stage, ()):
        folder = root / name
        files[name] = {
            "exists": folder.is_dir(),
            "n_files": _count_files(folder, ("*",)) if folder.is_dir() else 0,
        }
        if name.endswith("_plots"):
            files[name]["plots"] = _list_rel(folder, PLOT_SUFFIXES)
            scores = _read_scores_table(folder)
            if scores:
                summary["quality_scores"] = scores
                summary["n_annotators"] = len(scores)
            opal = folder / "opal"
            if opal.is_dir():
                files[name]["opal"] = _list_rel(opal, PLOT_SUFFIXES | {".html", ".csv"})
        if name in {"initial", "regenerated"}:
            summary[f"n_{name}_r1"] = _count_files(folder, ("*_R1.fastq", "*_R1.fastq.gz"))
        if name.endswith("_reports"):
            summary["n_reports"] = _count_files(folder, ("*.out", "*.csv", "*.tsv", "*.kreport"))
        if name.endswith("_annotations") and not name.endswith("_plots"):
            summary["n_tables"] = _count_files(folder, ("*.csv",))
        if name == "genomes":
            summary["n_genomes"] = _count_files(folder, ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz"))
    return {
        "id": "samovar",
        "schema_version": 1,
        "samovar_version": PACKAGE_VERSION,
        "stage": stage,
        "section_name": info["title"],
        "description": info["description"],
        "created": _now(),
        "output_dir": str(root.resolve()) if root.exists() else str(root),
        "summary": summary,
        "files": files,
    }


def write_stage_report(output_dir: PathLike, stage: str, extra: Optional[Dict[str, Any]] = None) -> Path:
    dest = report_dir(output_dir)
    dest.mkdir(parents=True, exist_ok=True)
    payload = collect_stage_payload(output_dir, stage)
    if extra:
        payload["extra"] = _safe_json(extra)
    path = dest / JSON_NAME.format(stage=stage)
    path.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
    return path


def write_overview_report(output_dir: PathLike) -> Path:
    root = as_path(output_dir)
    dest = report_dir(root)
    dest.mkdir(parents=True, exist_ok=True)
    stages = []
    for name in CHECKPOINT_STEPS:
        path = dest / JSON_NAME.format(stage=name)
        if path.is_file():
            try:
                stages.append(json.loads(path.read_text(encoding="utf-8")))
            except json.JSONDecodeError:
                stages.append({"stage": name, "error": "invalid json"})
        else:
            stages.append({"stage": name, "missing": True})
    payload = {
        "id": "samovar",
        "schema_version": 1,
        "samovar_version": PACKAGE_VERSION,
        "stage": "overview",
        "section_name": "SAMOVAR report",
        "description": "Run options and quality plots for this SAMOVAR execution.",
        "created": _now(),
        "output_dir": str(root.resolve()) if root.exists() else str(root),
        "stages": stages,
    }
    path = dest / "overview.samovar.json"
    path.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
    return path


def _html_comment(meta: Dict[str, str]) -> str:
    lines = ["<!--"]
    for key, value in meta.items():
        if value is None:
            continue
        lines.append(f"{key}: {json.dumps(value)}")
    lines.append("-->")
    return "\n".join(lines) + "\n"


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_")


def plot_parent(out_dir: PathLike) -> tuple:
    name = as_path(out_dir).name
    for stage, spec in REPORT_SECTION.items():
        if name == spec["folder"]:
            return spec["id"], spec["title"], spec["description"]
    for stage, dirs in STAGE_DIRS.items():
        if name in dirs:
            info = STAGE_INFO[stage]
            return f"samovar_{stage}", info["title"], info["description"]
    return "samovar_plots", "Annotation plots", ""


def report_parent(stage: str, read_type: Optional[str] = None) -> tuple:
    spec = REPORT_SECTION.get(stage, {})
    base_id = spec.get("id") or f"samovar_{stage}"
    title = spec.get("title") or STAGE_INFO.get(stage, {}).get("title") or stage
    desc = spec.get("description") or STAGE_INFO.get(stage, {}).get("description") or ""
    if read_type and is_platform_read_type(read_type):
        label = str(read_type).strip()
        pretty = label.upper() if label in {"ont"} else label.capitalize()
        return (
            f"{base_id}_{_slug(label).lower()}",
            f"{title} — {pretty}",
            f"{desc} Restricted to {pretty} reads.",
            True,
        )
    return base_id, title, desc, False


def _metric_label(key: str) -> str:
    return SCORE_DISPLAY.get(key) or OPAL_DISPLAY.get(key) or key.replace("_", " ")


def _annotator_color(name: str, index: int = 0) -> str:
    key = canonical_annotator_name(name)
    if key in ANNOTATOR_COLORS:
        return ANNOTATOR_COLORS[key]
    low = str(name).strip().lower()
    if low in ANNOTATOR_COLORS:
        return ANNOTATOR_COLORS[low]
    return ANNOTATOR_COLOR_FALLBACK[index % len(ANNOTATOR_COLOR_FALLBACK)]


def _split_plot_stem(stem: str) -> tuple:
    text = stem.replace("_mqc", "")
    if "." in text:
        base, maybe_rt = text.rsplit(".", 1)
        if is_platform_read_type(maybe_rt):
            return base, maybe_rt.lower()
    return text, None


def _is_finite_number(val: Any) -> bool:
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        return False
    return val == val and abs(val) != float("inf")


def _heatmap_pconfig(
    plot_id: str,
    title: str,
    xlab: str,
    ylab: str,
    min_value: Optional[float] = 0,
) -> Dict[str, Any]:
    # Keep FPC / true-vs-pred order; do not treat taxIDs as sample names.
    cfg: Dict[str, Any] = {
        "id": plot_id,
        "title": title,
        "xlab": xlab,
        "ylab": ylab,
        "square": False,
        "xcats_samples": False,
        "ycats_samples": False,
        "cluster_rows": False,
        "cluster_cols": False,
        "colstops": HEATMAP_COLSTOPS,
    }
    if min_value is not None:
        cfg["min"] = min_value
    return cfg


def write_table_mqc(
    rows: Sequence[Dict[str, Any]],
    path: PathLike,
    *,
    section_name: str,
    description: str = "",
    col1_header: str = "Method",
    plot_id: Optional[str] = None,
    parent_id: Optional[str] = None,
    parent_name: Optional[str] = None,
    id_field: str = "mode",
    numeric_fields: Optional[Sequence[str]] = None,
) -> Path:
    """Native MultiQC table (same shape as quality_scores)."""
    dest = as_path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pid, pname, pdesc = plot_parent(dest.parent)
    stem = dest.stem.replace("_mqc", "")
    cid = plot_id or f"{pid}_{stem}"
    data: Dict[str, Dict[str, Any]] = {}
    headers: Dict[str, Dict[str, str]] = {}
    skip = {id_field, "annotator", "scorer", "ok", "error", "details", "mode"}
    for row in rows:
        name = str(row.get(id_field) or row.get("mode") or "method")
        entry: Dict[str, Any] = {}
        keys = list(numeric_fields) if numeric_fields else [
            k for k, v in row.items() if k not in skip and _is_finite_number(v)
        ]
        for key in keys:
            val = row.get(key)
            if not _is_finite_number(val):
                continue
            entry[key] = round(float(val), 4)
            headers.setdefault(key, {"title": _metric_label(key), "format": "{:.3f}"})
        if entry:
            data[name] = entry
    payload = {
        "id": cid,
        "section_id": cid,
        "parent_id": parent_id or pid,
        "parent_name": parent_name or pname,
        "parent_description": pdesc,
        "section_name": section_name,
        "description": description,
        "plot_type": "table",
        "pconfig": {
            "id": f"{cid}_table",
            "title": section_name,
            "col1_header": col1_header,
        },
        "headers": headers,
        "data": data,
    }
    dest.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
    return dest


def write_bargraph_mqc(
    series: Dict[str, Dict[str, float]],
    path: PathLike,
    *,
    section_name: str,
    description: str = "",
    xlab: str = "Score",
    plot_id: Optional[str] = None,
    parent_id: Optional[str] = None,
    parent_name: Optional[str] = None,
    ymin: Optional[float] = None,
    ymax: Optional[float] = None,
) -> Path:
    """Native MultiQC bar graph (grouped, same defaults as score bars)."""
    dest = as_path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pid, pname, pdesc = plot_parent(dest.parent)
    stem = dest.stem.replace("_mqc", "")
    cid = plot_id or f"{pid}_{stem}"
    names: List[str] = []
    seen = set()
    for metric_series in series.values():
        for name in metric_series:
            if name not in seen:
                seen.add(name)
                names.append(name)
    cats = {
        name: {"name": name, "color": _annotator_color(name, i)}
        for i, name in enumerate(names)
    }
    pconfig: Dict[str, Any] = {
        "id": f"{cid}_plot",
        "title": section_name,
        "xlab": xlab,
        "stacking": "group",
        "cpswitch": False,
        "sort_samples": False,
    }
    if ymin is not None:
        pconfig["ymin"] = ymin
    if ymax is not None:
        pconfig["ymax"] = ymax
    payload = {
        "id": cid,
        "section_id": cid,
        "parent_id": parent_id or pid,
        "parent_name": parent_name or pname,
        "parent_description": pdesc,
        "section_name": section_name,
        "description": description,
        "plot_type": "bargraph",
        "categories": cats,
        "pconfig": pconfig,
        "data": series,
    }
    dest.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
    return dest


def write_heatmap_mqc(
    matrix,
    path: PathLike,
    *,
    section_name: str,
    description: str,
    xlab: str,
    ylab: str,
    plot_id: Optional[str] = None,
    parent_id: Optional[str] = None,
    parent_name: Optional[str] = None,
    min_value: Optional[float] = 0,
) -> Path:
    """Native MultiQC heatmap (selectable + ``--export``)."""
    dest = as_path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pid, pname, pdesc = plot_parent(dest.parent)
    stem = dest.stem.replace("_mqc", "")
    cid = plot_id or f"{pid}_{stem}"
    values = []
    for row in matrix.to_numpy(dtype=float).tolist():
        values.append([None if not _is_finite_number(v) else v for v in row])
    payload = {
        "id": cid,
        "section_id": cid,
        "parent_id": parent_id or pid,
        "parent_name": parent_name or pname,
        "parent_description": pdesc,
        "section_name": section_name,
        "description": description,
        "plot_type": "heatmap",
        "xcats": [str(c) for c in matrix.columns],
        "ycats": [str(c) for c in matrix.index],
        "pconfig": _heatmap_pconfig(
            f"{cid}_plot", section_name, xlab, ylab, min_value=min_value
        ),
        "data": values,
    }
    dest.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
    return dest


def write_scatter_mqc(
    table,
    path: PathLike,
    *,
    section_name: str,
    description: str,
    xlab: str,
    ylab: str,
    x: str = "pred_n",
    y: str = "true_n",
    label: str = "true",
    plot_id: Optional[str] = None,
    parent_id: Optional[str] = None,
    parent_name: Optional[str] = None,
) -> Path:
    """Native MultiQC scatter (selectable + ``--export``)."""
    dest = as_path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pid, pname, pdesc = plot_parent(dest.parent)
    stem = dest.stem.replace("_mqc", "")
    cid = plot_id or f"{pid}_{stem}"
    data = {}
    for _, row in table.iterrows():
        key = str(row.get(label, len(data)))
        data[key] = {"x": float(row[x]), "y": float(row[y])}
    payload = {
        "id": cid,
        "section_id": cid,
        "parent_id": parent_id or pid,
        "parent_name": parent_name or pname,
        "parent_description": pdesc,
        "section_name": section_name,
        "description": description,
        "plot_type": "scatter",
        "pconfig": {
            "id": f"{cid}_plot",
            "title": section_name,
            "xlab": xlab,
            "ylab": ylab,
        },
        "data": data,
    }
    dest.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
    return dest


def _normalize_mqc_json(
    payload: Dict[str, Any],
    *,
    parent_id: str,
    parent_name: str,
    parent_description: str,
    stem: str,
) -> Dict[str, Any]:
    """Make IDs unique per stage and fix Plotly defaults that break taxID heatmaps."""
    cid = f"{parent_id}_{stem}"
    payload["id"] = cid
    payload["section_id"] = cid
    payload["parent_id"] = parent_id
    payload["parent_name"] = parent_name
    payload["parent_description"] = parent_description
    section = str(payload.get("section_name") or stem).strip()
    if parent_name and not section.startswith(parent_name):
        payload["section_name"] = f"{parent_name} — {section}"
    pconfig = dict(payload.get("pconfig") or {})
    pconfig["id"] = f"{cid}_plot"
    if payload.get("plot_type") == "heatmap":
        pconfig.setdefault("xcats_samples", False)
        pconfig.setdefault("ycats_samples", False)
        pconfig.setdefault("cluster_rows", False)
        pconfig.setdefault("cluster_cols", False)
        pconfig.setdefault("square", False)
        pconfig.setdefault("colstops", HEATMAP_COLSTOPS)
    if payload.get("plot_type") == "bargraph":
        pconfig.setdefault("stacking", "group")
        pconfig.setdefault("cpswitch", False)
        pconfig.setdefault("sort_samples", False)
    payload["pconfig"] = pconfig
    return payload


def _load_yaml(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        import yaml

        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _short_path(value: Any) -> str:
    text = str(value).strip()
    if not text:
        return ""
    if "/" in text or text.startswith("."):
        return Path(text).name or text
    return text


def collect_run_options(output_dir: PathLike) -> List[tuple]:
    """Flatten this run's generate / annotate / regenerate YAML into display rows."""
    root = as_path(output_dir)
    rows: List[tuple] = []
    seen = set()

    def _add(label: str, value: Any) -> None:
        if value is None or value == "" or value == []:
            return
        key = label.lower()
        if key in seen:
            return
        seen.add(key)
        if isinstance(value, (list, tuple)):
            text = ", ".join(_short_path(v) for v in value if str(v).strip())
        elif isinstance(value, dict):
            text = ", ".join(f"{k}={_short_path(v)}" for k, v in value.items())
        elif isinstance(value, bool):
            text = "yes" if value else "no"
        else:
            text = _short_path(value)
        if text:
            rows.append((label, text))

    generate = _load_yaml(root / ".generate" / "configs" / "camisim.yaml")
    if not generate:
        generate = _load_yaml(root / ".generate" / "configs" / "iss_config.yaml")
    init_cfg = _load_yaml(root / ".log" / "configs" / "config_init.yaml")
    regen_cfg = _load_yaml(root / ".log" / "configs" / "config_annotation2iss.yaml")

    _add("SAMOVAR", PACKAGE_VERSION)
    _add("Simulator", generate.get("simulator"))
    _add("Mode", generate.get("mode"))
    _add("Read types", generate.get("camisim_types"))
    _add("Samples", generate.get("n_samples"))
    _add("Total reads", generate.get("total_reads") or generate.get("n_reads"))
    _add("Size (Gbp)", generate.get("size_gbp"))
    _add("Host fraction", generate.get("host_fraction"))
    _add("ISS model", generate.get("iss_model") or generate.get("model"))
    _add("Read length", generate.get("read_length") or regen_cfg.get("read_length"))
    _add("Host genome", generate.get("host_genome"))
    _add("Generate seed", generate.get("seed"))
    _add("Cores", generate.get("cores") or regen_cfg.get("cores"))
    annotators = init_cfg.get("run_config") or []
    if isinstance(annotators, list) and annotators:
        names = []
        for item in annotators:
            if not isinstance(item, dict):
                continue
            names.append(str(item.get("type") or item.get("run_name") or "annotator"))
        _add("Annotators", names)
    _add("Regeneration mode", regen_cfg.get("regeneration_mode"))
    modes = regen_cfg.get("table_reads_generators")
    if isinstance(modes, list) and len(modes) > 1:
        _add("Table generators", modes)
        _add("Table score", regen_cfg.get("table_score") or "shannon_ks")
    selection = root / "regenerated" / ".regenerated_abundance" / "table_selection.json"
    if selection.is_file():
        try:
            payload = json.loads(selection.read_text(encoding="utf-8"))
        except Exception:
            payload = {}
        if payload.get("winner"):
            _add("Selected table generator", payload.get("winner"))
    _add("Regenerated reads", regen_cfg.get("N_reads"))
    _add("Coverage", regen_cfg.get("coverage"))
    _add("Max genomes", regen_cfg.get("max_genomes"))
    _add("Max genome MB", regen_cfg.get("max_genome_mb"))
    skip = regen_cfg.get("genome_skip_list")
    if skip:
        _add("Genome skip list", skip)
    _add("Regeneration seed", regen_cfg.get("seed"))
    _add("Rescale abundance", regen_cfg.get("rescale_abundance"))
    return rows


def _score_read_type(row: Dict[str, Any]) -> str:
    if "read_type" not in row:
        return "all"
    token = str(row.get("read_type", "all")).strip().lower()
    if token in {"", "nan", "none", "na"}:
        return "all"
    if is_platform_read_type(token):
        return token
    if token == "all":
        return "all"
    return ""


def _dedupe_score_rows(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        raw = str(row.get("annotator", "annotator"))
        name = canonical_annotator_name(raw)
        key = name.lower()
        incoming_vote = raw.strip().lower() == "consensus"
        if key in by:
            existing_vote = str(by[key].get("annotator", "")).strip().lower() == "consensus"
            if incoming_vote and not existing_vote:
                continue
        rec = dict(row)
        rec["annotator"] = name
        by[key] = rec
    return list(by.values())


def _rows_for_read_type(scores: Sequence[Dict[str, Any]], read_type: str) -> List[Dict[str, Any]]:
    want = (read_type or "all").strip().lower()
    picked = [row for row in scores if _score_read_type(row) == want]
    return _dedupe_score_rows(picked)


def _available_read_types(scores: Sequence[Dict[str, Any]]) -> List[str]:
    types = []
    seen = set()
    for row in scores:
        rt = _score_read_type(row)
        if rt and rt != "all" and rt not in seen:
            seen.add(rt)
            types.append(rt)
    return types


def _table_payload(rows: Sequence[Dict[str, Any]]) -> tuple:
    data: Dict[str, Dict[str, Any]] = {}
    headers: Dict[str, Dict[str, str]] = {}
    skip = {"annotator", "read_type"}
    for row in rows:
        name = str(row.get("annotator", "annotator"))
        entry: Dict[str, Any] = {}
        for key, val in row.items():
            if key in skip or not _is_finite_number(val):
                continue
            entry[key] = round(float(val), 4)
            headers.setdefault(key, {"title": _metric_label(key), "format": "{:.3f}"})
        if entry:
            data[name] = entry
    used = {k for rec in data.values() for k in rec}
    headers = {k: v for k, v in headers.items() if k in used}
    return data, headers


def _bars_by_metric(rows: Sequence[Dict[str, Any]]) -> tuple:
    metrics = [
        key
        for key in BAR_METRIC_ORDER
        if any(_is_finite_number(row.get(key)) for row in rows)
    ]
    data: Dict[str, Dict[str, float]] = {}
    annotators: List[str] = []
    seen = set()
    for row in rows:
        name = str(row.get("annotator", "annotator"))
        if name.lower() not in seen:
            seen.add(name.lower())
            annotators.append(name)
    for metric in metrics:
        series = {}
        for row in rows:
            val = row.get(metric)
            if _is_finite_number(val):
                series[str(row["annotator"])] = round(float(val), 4)
        if series:
            data[_metric_label(metric)] = series
    cats = {
        name: {"name": name, "color": _annotator_color(name, i)}
        for i, name in enumerate(annotators)
    }
    return data, cats


def _write_score_plots(
    dest: Path,
    *,
    prefix: str,
    parent_id: str,
    parent_name: str,
    parent_description: str,
    rows: Sequence[Dict[str, Any]],
) -> List[Path]:
    written: List[Path] = []
    if not rows:
        return written
    data, headers = _table_payload(rows)
    if data:
        payload = {
            "id": f"{parent_id}_quality_scores",
            "section_id": f"{parent_id}_quality_scores",
            "parent_id": parent_id,
            "parent_name": parent_name,
            "parent_description": parent_description,
            "section_name": f"{parent_name} — quality scores",
            "description": "",
            "plot_type": "table",
            "pconfig": {
                "id": f"{parent_id}_quality_scores_table",
                "title": f"{parent_name} quality scores",
                "col1_header": "Annotator",
            },
            "headers": headers,
            "data": data,
        }
        out = dest / f"{prefix}_quality_scores_mqc.json"
        out.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
        written.append(out)
    bars, cats = _bars_by_metric(rows)
    if bars:
        bar_payload = {
            "id": f"{parent_id}_score_bars",
            "section_id": f"{parent_id}_score_bars",
            "parent_id": parent_id,
            "parent_name": parent_name,
            "parent_description": parent_description,
            "section_name": f"{parent_name} — score bars",
            "description": "",
            "plot_type": "bargraph",
            "categories": cats,
            "pconfig": {
                "id": f"{parent_id}_score_bars_plot",
                "title": f"{parent_name} scores",
                "xlab": "Score",
                "ymin": 0,
                "ymax": 1,
                "stacking": "group",
                "cpswitch": False,
                "sort_samples": False,
            },
            "data": bars,
        }
        out = dest / f"{prefix}_score_bars_mqc.json"
        out.write_text(json.dumps(_safe_json(bar_payload), indent=2) + "\n", encoding="utf-8")
        written.append(out)
    return written


def _skip_duplicate_ensemble_plot(path: Path, names: Sequence[str]) -> bool:
    stem, _rt = _split_plot_stem(path.stem.replace("_mqc", ""))
    token = stem.split("_", 1)[-1] if "_" in stem else stem
    if token.lower() != "consensus":
        return False
    want = re.sub(r"consensus", ENSEMBLE_NAME, path.name, flags=re.I)
    return want.lower() in {n.lower() for n in names}


def _copy_plot_assets(plots: Path, dest: Path, stage: str, order: int) -> List[Path]:
    written: List[Path] = []
    if not plots.is_dir():
        return written
    sources = sorted(plots.glob("*_mqc.json"))
    names = [p.name for p in sources]
    for src in sources:
        if _skip_duplicate_ensemble_plot(src, names):
            continue
        try:
            payload = json.loads(src.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if not isinstance(payload, dict):
            continue
        stem, read_type = _split_plot_stem(src.stem.replace("_mqc", ""))
        parent_id, parent_name, parent_description, _hidden = report_parent(stage, read_type)
        if "consensus" in stem.lower():
            stem = re.sub(r"consensus", ENSEMBLE_NAME, stem, flags=re.I)
            section = str(payload.get("section_name") or "")
            payload["section_name"] = re.sub(r"consensus", ENSEMBLE_NAME, section, flags=re.I)
        payload = _normalize_mqc_json(
            payload,
            parent_id=parent_id,
            parent_name=parent_name,
            parent_description=parent_description,
            stem=_slug(stem),
        )
        prefix = f"{order:02d}_{_slug(parent_name)}"
        out = dest / f"{prefix}_{src.name}"
        out.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
        written.append(out)
    scores = _read_scores_table(plots)
    types = ["all"] + _available_read_types(scores)
    for read_type in types:
        rows = _rows_for_read_type(scores, read_type)
        if not rows:
            continue
        parent_id, parent_name, parent_description, _hidden = report_parent(
            stage, None if read_type == "all" else read_type
        )
        prefix = f"{order:02d}_{_slug(parent_name)}"
        written.extend(
            _write_score_plots(
                dest,
                prefix=prefix,
                parent_id=parent_id,
                parent_name=parent_name,
                parent_description=parent_description,
                rows=rows,
            )
        )
    return written


def _run_options_html(output_dir: PathLike, hidden_ids: Sequence[str]) -> str:
    wiki = "https://github.com/ctlab/samovar/wiki"
    rows = collect_run_options(output_dir)
    parts = [
        "<p>SAMOVAR takes a metagenome (simulated or real), annotates it with "
        "the tools you configured, rebuilds a community from those calls, and "
        "trains a reprofiler on the known mix. The sections below follow that "
        "order: <strong>Raw</strong> (input), <strong>Regenerated</strong> "
        "(simulated community), <strong>Reprofiled</strong> (ensemble on the original samples).</p>",
        f'<p>Documentation: <a href="{wiki}">{wiki}</a> '
        f'(pipeline overview, <a href="{wiki}/Home">Home</a>).</p>',
        "<p>Options used for this SAMOVAR run:</p>",
    ]
    if rows:
        parts.append('<dl class="dl-horizontal" style="margin-top:0.75em">')
        for label, value in rows:
            parts.append(
                f"<dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd>"
            )
        parts.append("</dl>")
    else:
        parts.append("<p><em>No run-config YAML was found under <code>.log/configs</code> or <code>.generate/configs</code>.</em></p>")
    if hidden_ids:
        parts.append(
            '<p><button type="button" class="btn btn-default btn-sm" id="samovar-toggle-platforms">'
            "Show per-platform plots</button> "
            '<span class="text-muted">Illumina / ONT / other read-type sections are hidden by default.</span></p>'
            "<script>\n"
            "(function(){var b=document.getElementById('samovar-toggle-platforms');"
            "if(!b)return;b.addEventListener('click',function(){"
            "document.body.classList.toggle('samovar-show-platforms');"
            "b.textContent=document.body.classList.contains('samovar-show-platforms')"
            "?'Hide per-platform plots':'Show per-platform plots';});})();"
            "</script>\n"
        )
    return "\n".join(parts)


def _improvement_pct(ensemble: float, baseline: float) -> Optional[float]:
    if not _is_finite_number(ensemble) or not _is_finite_number(baseline):
        return None
    if abs(baseline) < 1e-12:
        return None if abs(ensemble) < 1e-12 else None
    return 100.0 * (float(ensemble) - float(baseline)) / abs(float(baseline))


def _best_tool_row(rows: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    tools = [r for r in rows if not str(r.get("annotator", "")).lower() == ENSEMBLE_NAME.lower()]
    if not tools:
        return None

    def _rank(row: Dict[str, Any]) -> tuple:
        vals = []
        for key in ("f1", "accuracy", "f1_macro", "completeness"):
            val = row.get(key)
            vals.append(float(val) if _is_finite_number(val) else -1e9)
        return tuple(vals)

    return max(tools, key=_rank)


def _conclusion_html(output_dir: PathLike) -> str:
    root = as_path(output_dir)
    chosen_rows: List[Dict[str, Any]] = []
    stage_title = ""
    for stage in ("viz_reprofiled", "viz_regenerated", "viz_initial"):
        folder = REPORT_SECTION[stage]["folder"]
        scores = _read_scores_table(root / folder)
        rows = _rows_for_read_type(scores, "all")
        if rows:
            chosen_rows = rows
            stage_title = REPORT_SECTION[stage]["title"]
            break
    if not chosen_rows:
        return "<p>No quality scores were available to summarise.</p>"
    by_name = {str(r["annotator"]): r for r in chosen_rows}
    ensemble = by_name.get(ENSEMBLE_NAME)
    best = _best_tool_row(chosen_rows)
    parts = [
        f"<p>Summary from the <strong>{html.escape(stage_title)}</strong> stage"
        " (all reads). Metrics that need ground truth are omitted when the "
        "starting metagenome has no usable <code>true</code> taxIDs.</p>"
    ]
    if best:
        parts.append(f"<h4>Best single tool: {html.escape(str(best['annotator']))}</h4><ul>")
        for key in CONCLUSION_METRICS:
            val = best.get(key)
            if _is_finite_number(val):
                parts.append(f"<li>{html.escape(_metric_label(key))}: {float(val):.3f}</li>")
        parts.append("</ul>")
        if not any(_is_finite_number(best.get(k)) for k in CONCLUSION_METRICS):
            n_taxa = best.get("n_taxa")
            n_reads = best.get("n_reads")
            extra = []
            if _is_finite_number(n_reads):
                extra.append(f"{int(n_reads)} reads")
            if _is_finite_number(n_taxa):
                extra.append(f"{int(n_taxa)} taxa")
            if extra:
                parts.append("<p>" + ", ".join(extra) + ".</p>")
    if ensemble and best:
        parts.append("<h4>SAMOVAR ensemble improvement</h4><ul>")
        any_metric = False
        for key in CONCLUSION_METRICS:
            ev = ensemble.get(key)
            bv = best.get(key)
            if not _is_finite_number(ev):
                continue
            any_metric = True
            pct = _improvement_pct(float(ev), float(bv)) if _is_finite_number(bv) else None
            if pct is None:
                parts.append(
                    f"<li>{html.escape(_metric_label(key))}: {float(ev):.3f}</li>"
                )
            else:
                sign = "+" if pct >= 0 else ""
                parts.append(
                    f"<li>{html.escape(_metric_label(key))}: {float(ev):.3f} "
                    f"({sign}{pct:.1f}% vs {html.escape(str(best['annotator']))})</li>"
                )
        if not any_metric:
            parts.append(
                "<li>Ground-truth F1 / R² were not calculated for this run.</li>"
            )
        parts.append("</ul>")
    elif ensemble:
        parts.append("<h4>SAMOVAR ensemble</h4><ul>")
        for key in CONCLUSION_METRICS:
            val = ensemble.get(key)
            if _is_finite_number(val):
                parts.append(f"<li>{html.escape(_metric_label(key))}: {float(val):.3f}</li>")
        parts.append("</ul>")
    else:
        parts.append("<p>SAMOVAR ensemble scores were not present in this run.</p>")
    return "\n".join(parts)


def _hidden_parent_ids(dest: Path) -> List[str]:
    ids = []
    for path in dest.glob("*_mqc.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        parent = str(payload.get("parent_id") or "")
        if parent.count("_") >= 2 and parent not in ids:
            # samovar_raw_illumina has an extra token after the stage id
            base = parent.rsplit("_", 1)[0]
            if base in {spec["id"] for spec in REPORT_SECTION.values()}:
                ids.append(parent)
    return ids


def _write_multiqc_config(dest: Path, hidden_ids: Sequence[str]) -> Path:
    css_path = dest / "samovar_multiqc.css"
    css_rules = [
        "/* SAMOVAR: Dark2-like tool colours; YlGn heatmaps are set in plot JSON. */",
        ".mqc-module-section-first { border-top: 3px solid #1B9E77; }",
    ]
    for hid in hidden_ids:
        css_rules.append(
            f"body:not(.samovar-show-platforms) #mqc-module-section-{hid},"
            f"body:not(.samovar-show-platforms) a[href='#{hid}'] {{ display: none !important; }}"
        )
    css_path.write_text("\n".join(css_rules) + "\n", encoding="utf-8")
    order = ["samovar_run_options"]
    for stage in REPORT_STAGES:
        spec = REPORT_SECTION[stage]
        order.append(spec["id"])
        for hid in hidden_ids:
            if hid.startswith(spec["id"] + "_"):
                order.append(hid)
    order.append("samovar_conclusion")
    lines = [
        f'title: "SAMOVAR"',
        f'subtitle: "Ensemble annotation report (v{PACKAGE_VERSION}). Wiki: https://github.com/ctlab/samovar/wiki"',
        "intro_text: SAMOVAR report",
        f'custom_css_files:',
        f'  - {json.dumps(str(css_path))}',
        "custom_content:",
        "  order:",
    ]
    for item in order:
        lines.append(f"    - {item}")
    cfg = dest / "multiqc_config.yaml"
    cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cfg


def bundle_multiqc(output_dir: PathLike) -> Path:
    """Materialize ``*_mqc.*`` custom content from plot folders and run options."""
    root = as_path(output_dir)
    write_overview_report(root)
    dest = multiqc_stage_dir(root)
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)

    plots = root / REPORT_SECTION["viz_regenerated"]["folder"]
    if not plots.is_dir() or not any(plots.glob("TableScore_*_mqc.json")):
        for candidate in (
            root / "regenerated" / ".regenerated_abundance" / "table_selection.json",
            root / "regenerated" / "table_selection.json",
        ):
            if candidate.is_file():
                try:
                    from samovar.table_scorers import write_table_score_plots

                    write_table_score_plots(
                        json.loads(candidate.read_text(encoding="utf-8")),
                        plots,
                    )
                except Exception:
                    pass
                break

    for i, stage in enumerate(REPORT_STAGES, start=1):
        spec = REPORT_SECTION[stage]
        intro = dest / f"{i:02d}_{_slug(spec['title'])}_intro_mqc.html"
        intro.write_text(
            _html_comment(
                {
                    "id": f"{spec['id']}_intro",
                    "section_name": spec["title"],
                    "description": spec["description"],
                    "parent_id": spec["id"],
                    "parent_name": spec["title"],
                    "parent_description": spec["description"],
                    "plot_type": "html",
                }
            )
            + f"<p>{html.escape(spec['description'])}</p>\n",
            encoding="utf-8",
        )
        folder = root / spec["folder"]
        _copy_plot_assets(folder, dest, stage, i)

    hidden_ids = _hidden_parent_ids(dest)
    options_html = dest / "00_run_options_mqc.html"
    options_html.write_text(
        _html_comment(
            {
                "id": "samovar_run_options",
                "section_name": "Run options",
                "description": "Settings used for this SAMOVAR run.",
                "plot_type": "html",
            }
        )
        + _run_options_html(root, hidden_ids)
        + "\n",
        encoding="utf-8",
    )
    conclusion_html = dest / "99_conclusion_mqc.html"
    conclusion_html.write_text(
        _html_comment(
            {
                "id": "samovar_conclusion",
                "section_name": "Conclusion",
                "description": "Best single-tool stats and SAMOVAR ensemble improvement.",
                "plot_type": "html",
            }
        )
        + _conclusion_html(root)
        + "\n",
        encoding="utf-8",
    )
    _write_multiqc_config(dest, hidden_ids)
    return dest


def run_multiqc(output_dir: PathLike, extra_args: Optional[Sequence[str]] = None) -> int:
    """Invoke the real MultiQC CLI on staged ``*_mqc`` files.

    Extra args are passed through (after an optional ``--``), so you can do
    ``samovar multiqc --output_dir RUN -- --export --module custom_content``.
    Native Plotly plots (heatmap / scatter / bar / table) show MultiQC's
    in-report export checkboxes; ``--export`` / ``-p`` writes png/svg/pdf.
    """
    import subprocess

    from samovar.paths import discover_multiqc

    staged = bundle_multiqc(output_dir)
    out = as_path(output_dir) / "multiqc"
    out.mkdir(parents=True, exist_ok=True)
    cfg = staged / "multiqc_config.yaml"
    extra = list(extra_args or [])
    if extra and extra[0] == "--":
        extra = extra[1:]
    exe = discover_multiqc()
    if exe:
        cmd = [exe, str(staged), "-o", str(out), "-f", "-c", str(cfg)]
    else:
        cmd = [sys.executable, "-m", "multiqc", str(staged), "-o", str(out), "-f", "-c", str(cfg)]
    if extra:
        cmd.extend(extra)
    else:
        cmd.append("--interactive")
    print(" ".join(cmd), flush=True)
    try:
        rc = subprocess.call(cmd)
    except FileNotFoundError:
        rc = 127
    if rc != 0:
        print(
            f"MultiQC was not run (exit {rc}). Stage JSON and *_mqc files are in {staged}. "
            "Optional install: ./install.sh MultiQC",
            file=sys.stderr,
        )
        return 0
    return rc


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m samovar.stage_report",
        description="Write SamovaR JSON summaries and MultiQC custom content",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    emit = sub.add_parser("emit", help="Write JSON for one completed stage")
    emit.add_argument("output_dir")
    emit.add_argument("stage")
    overview = sub.add_parser("overview", help="Write overview JSON from existing stage files")
    overview.add_argument("output_dir")
    bundle = sub.add_parser("bundle", help="Stage *_mqc files for MultiQC")
    bundle.add_argument("output_dir")
    render = sub.add_parser("multiqc", help="Bundle and run the MultiQC CLI")
    render.add_argument("output_dir")
    render.add_argument("multiqc_args", nargs=argparse.REMAINDER, help="Passed to multiqc (use -- first)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(list(argv) if argv is not None else None)
    if args.command == "emit":
        write_stage_report(args.output_dir, args.stage)
        return 0
    if args.command == "overview":
        write_overview_report(args.output_dir)
        return 0
    if args.command == "bundle":
        bundle_multiqc(args.output_dir)
        return 0
    if args.command == "multiqc":
        return run_multiqc(args.output_dir, extra_args=getattr(args, "multiqc_args", None))
    return 2


if __name__ == "__main__":
    sys.exit(main())
