"""FastQC/fastp-style JSON summaries plus MultiQC custom-content staging.

Each exec checkpoint writes ``<output>/.log/multiqc/<stage>.samovar.json``.
``bundle_multiqc`` copies Altair HTML, cnsplots/OPAL PNGs, and score tables into
``<output>/multiqc_samovar/`` as ``*_mqc.*`` files that stock MultiQC can render.
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

PathLike = Union[str, os.PathLike]

STAGE_INFO: Dict[str, Dict[str, str]] = {
    "setup_reads": {
        "title": "Read setup",
        "description": (
            "Source FASTQ is linked or copied into ``initial/``. "
            "ISS generate writes paired reads here when the run starts from genomes."
        ),
    },
    "annotate_initial": {
        "title": "Initial annotation",
        "description": (
            "Each configured classifier (Kraken2, Kaiju, …) labels the real or "
            "ISS reads. Per-tool reports land in ``initial_reports/``."
        ),
    },
    "combine_initial": {
        "title": "Combine initial calls",
        "description": (
            "Tool outputs are merged into one table (``initial_annotations/``) "
            "with a column per annotator and optional ``true`` taxIDs from ISS headers."
        ),
    },
    "viz_initial": {
        "title": "Initial quality plots",
        "description": (
            "Altair HTML for F1 heatmaps, R² scatters, cross-validation, and score bars "
            "on the first annotation pass (``initial_annotations_plots/``)."
        ),
    },
    "seed_genomes": {
        "title": "Genome library",
        "description": (
            "NCBI / cache assemblies needed to re-simulate the community are "
            "seeded into ``genomes/``. Bundled ``data/test_genomes`` stubs are not used "
            "unless this is an ISS toy run."
        ),
    },
    "regenerate_reads": {
        "title": "Community regeneration",
        "description": (
            "annotation2iss rebuilds an in-silico community from the initial "
            "calls (``regenerated/``) so the ensemble can be trained against a known mix."
        ),
    },
    "sort_reads": {
        "title": "Sort regenerated reads",
        "description": "Regenerated FASTQ is sorted so re-annotation and ML see a stable read order.",
    },
    "annotate_regenerated": {
        "title": "Re-annotation",
        "description": "The same classifiers are run on regenerated reads (``regenerated_reports/``).",
    },
    "combine_regenerated": {
        "title": "Combine re-annotation",
        "description": "Regenerated tool calls are merged (``regenerated_annotations/``) for training labels.",
    },
    "viz_regenerated": {
        "title": "Regenerated quality plots",
        "description": (
            "The same Altair suite as the initial pass, now on the simulated "
            "community (``regenerated_annotations_plots/``)."
        ),
    },
    "reprofile": {
        "title": "SAMOVAR re-profiling",
        "description": (
            "Supervised ensemble (``workflow/ML.py``) maps tool votes to a corrected "
            "taxID column ``taxid_SAMOVAR`` in ``reprofiled_annotations/``."
        ),
    },
    "viz_reprofiled": {
        "title": "Re-profiled quality plots",
        "description": (
            "Altair F1 / R² / CV / score charts after the ensemble "
            "(``reprofiled_annotations_plots/``)."
        ),
    },
}

STAGE_DIRS = {
    "setup_reads": ("initial",),
    "annotate_initial": ("initial_reports",),
    "combine_initial": ("initial_annotations",),
    "viz_initial": ("initial_annotations_plots",),
    "seed_genomes": ("genomes",),
    "regenerate_reads": ("regenerated",),
    "sort_reads": ("regenerated",),
    "annotate_regenerated": ("regenerated_reports",),
    "combine_regenerated": ("regenerated_annotations",),
    "viz_regenerated": ("regenerated_annotations_plots",),
    "reprofile": ("reprofiled_annotations",),
    "viz_reprofiled": ("reprofiled_annotations_plots",),
}

PLOT_SUFFIXES = {".png", ".html", ".svg", ".pdf"}
JSON_NAME = "{stage}.samovar.json"


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
        "section_name": "SamovaR pipeline",
        "description": (
            "Ensemble taxonomic annotation, community regeneration, and ML re-profiling. "
            "Each stage below has a Fastp-style JSON summary in ``.log/multiqc/``; "
            "Altair, cnsplots, and OPAL figures are staged for MultiQC."
        ),
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
    for stage, dirs in STAGE_DIRS.items():
        if name in dirs:
            info = STAGE_INFO[stage]
            return f"samovar_{stage}", info["title"], info["description"]
    return "samovar_plots", "Annotation plots", ""


def _heatmap_pconfig(plot_id: str, title: str, xlab: str, ylab: str) -> Dict[str, Any]:
    # Keep FPC / true-vs-pred order; do not treat taxIDs as sample names.
    return {
        "id": plot_id,
        "title": title,
        "xlab": xlab,
        "ylab": ylab,
        "min": 0,
        "square": False,
        "xcats_samples": False,
        "ycats_samples": False,
        "cluster_rows": False,
        "cluster_cols": False,
    }


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
) -> Path:
    """Native MultiQC heatmap (selectable + ``--export``)."""
    dest = as_path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pid, pname, pdesc = plot_parent(dest.parent)
    stem = dest.stem.replace("_mqc", "")
    cid = plot_id or f"{pid}_{stem}"
    values = matrix.to_numpy(dtype=float).tolist()
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
        "pconfig": _heatmap_pconfig(f"{cid}_plot", section_name, xlab, ylab),
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
    if payload.get("plot_type") == "bargraph":
        pconfig.setdefault("stacking", "group")
        pconfig.setdefault("cpswitch", False)
    payload["pconfig"] = pconfig
    return payload


def _altair_caption(stem: str) -> tuple:
    name = stem.lower()
    if name == "scores" or name.endswith("_scores") and "opal" not in name:
        return (
            "Quality scores (Altair)",
            "Interactive grouped bars of per-annotator classification metrics "
            "(accuracy, F1, purity, TN/FN rates).",
        )
    if "opal_scores" in name:
        return (
            "OPAL-style scores (Altair)",
            "Interactive bars for completeness, purity, OPAL F1, Jaccard, L1, and Bray–Curtis.",
        )
    if name.startswith("f1_") or "/f1_" in name:
        tool = stem.split("_", 1)[-1]
        return (
            f"F1 heatmap — {tool} (Altair)",
            "Interactive confusion heatmap of predicted vs true taxa for this annotator.",
        )
    if name.startswith("r2_") or "/r2_" in name:
        tool = stem.split("_", 1)[-1]
        return (
            f"Abundance R² — {tool} (Altair)",
            "Interactive scatter of predicted vs true taxon counts (R² on abundances).",
        )
    if name.startswith("cv_"):
        pair = stem[3:].replace("_", " ")
        return (
            f"Cross-validation — {pair} (Altair)",
            "Interactive heatmap of taxID agreement between two annotators.",
        )
    return (
        f"{stem} (Altair)",
        "Interactive Altair chart written by the SamovaR visualization step.",
    )


def _altair_mqc_html(src: Path, meta: Dict[str, str], vis_id: str) -> str:
    """Embed an Altair HTML chart so Vega still runs inside MultiQC."""
    raw = src.read_text(encoding="utf-8", errors="replace")
    spec_match = re.search(
        r"var spec = (\{.*?\});\s*var embedOpt",
        raw,
        flags=re.S,
    )
    parts = [_html_comment(meta)]
    if spec_match:
        spec = spec_match.group(1)
        parts.append(
            '<script src="https://cdn.jsdelivr.net/npm/vega@6"></script>\n'
            '<script src="https://cdn.jsdelivr.net/npm/vega-lite@6.4.1"></script>\n'
            '<script src="https://cdn.jsdelivr.net/npm/vega-embed@7"></script>\n'
            f'<div id="{html.escape(vis_id)}" style="width:100%;min-height:360px"></div>\n'
            "<script>\n"
            f"vegaEmbed(document.getElementById({json.dumps(vis_id)}), {spec}, "
            '{"mode": "vega-lite"});\n'
            "</script>\n"
        )
        return "".join(parts)
    body = re.search(r"<body[^>]*>(.*)</body>", raw, flags=re.I | re.S)
    blob = body.group(1) if body else raw
    blob = blob.replace('id="vis"', f'id="{vis_id}"').replace("#vis", f"#{vis_id}")
    parts.append(blob)
    parts.append("\n")
    return "".join(parts)


def _copy_plot_assets(plots: Path, dest: Path, stage: str, order: int) -> List[Path]:
    written: List[Path] = []
    if not plots.is_dir():
        return written
    info = STAGE_INFO.get(stage, {"title": stage, "description": ""})
    parent_id = f"samovar_{stage}"
    prefix = f"{order:02d}_{_slug(info['title'])}"
    for src in sorted(plots.glob("*_mqc.json")):
        try:
            payload = json.loads(src.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if not isinstance(payload, dict):
            continue
        payload = _normalize_mqc_json(
            payload,
            parent_id=parent_id,
            parent_name=info["title"],
            parent_description=info["description"],
            stem=src.stem.replace("_mqc", ""),
        )
        out = dest / f"{prefix}_{src.name}"
        out.write_text(json.dumps(_safe_json(payload), indent=2) + "\n", encoding="utf-8")
        written.append(out)
    for src in sorted(plots.glob("*.html")):
        if not src.is_file() or src.name.endswith("_mqc.html"):
            continue
        stem = _slug(src.stem)
        title, caption = _altair_caption(src.stem)
        html_id = f"{parent_id}_{stem}_altair"
        meta = {
            "id": html_id,
            "section_id": html_id,
            "parent_id": parent_id,
            "parent_name": info["title"],
            "parent_description": info["description"],
            "section_name": f"{info['title']} — {title}",
            "description": caption,
            "plot_type": "html",
        }
        out = dest / f"{prefix}_{stem}_altair_mqc.html"
        out.write_text(_altair_mqc_html(src, meta, f"vis_{html_id}"), encoding="utf-8")
        written.append(out)
    scores = _read_scores_table(plots)
    if scores:
        data = {}
        headers: Dict[str, Dict[str, str]] = {}
        skip = {"annotator", "n_reads", "n_taxa"}
        for row in scores:
            name = str(row.get("annotator", "annotator"))
            entry = {}
            for key, val in row.items():
                if key in skip:
                    continue
                if isinstance(val, (int, float)) and not isinstance(val, bool):
                    entry[key] = round(float(val), 4)
                    headers.setdefault(key, {"title": key.replace("_", " "), "format": "{:.3f}"})
            data[name] = entry
        payload = {
            "id": f"{parent_id}_quality_scores",
            "section_id": f"{parent_id}_quality_scores",
            "parent_id": parent_id,
            "parent_name": info["title"],
            "parent_description": info["description"],
            "section_name": f"{info['title']} — quality scores",
            "description": (
                "Per-annotator classification and OPAL-style profile metrics "
                "(same numbers as ``quality_scores.csv`` / ``scores.png``)."
            ),
            "plot_type": "table",
            "pconfig": {
                "id": f"{parent_id}_quality_scores_table",
                "title": f"{info['title']} quality scores",
                "col1_header": "Annotator",
            },
            "headers": headers,
            "data": data,
        }
        out = dest / f"{prefix}_quality_scores_mqc.json"
        out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        written.append(out)
        bar_metrics = [
            key
            for key in (
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
            if any(key in row for row in scores)
        ]
        if bar_metrics:
            bars = {}
            for row in scores:
                name = str(row.get("annotator", "annotator"))
                bars[name] = {
                    key: round(float(row[key]), 4)
                    for key in bar_metrics
                    if isinstance(row.get(key), (int, float))
                }
            bar_payload = {
                "id": f"{parent_id}_score_bars",
                "section_id": f"{parent_id}_score_bars",
                "parent_id": parent_id,
                "parent_name": info["title"],
                "parent_description": info["description"],
                "section_name": f"{info['title']} — score bars",
                "description": "Interactive MultiQC bar graph of the main 0–1 quality metrics.",
                "plot_type": "bargraph",
                "pconfig": {
                    "id": f"{parent_id}_score_bars_plot",
                    "title": f"{info['title']} scores",
                    "ylab": "Score",
                    "ymin": 0,
                    "ymax": 1,
                    "stacking": "group",
                    "cpswitch": False,
                },
                "data": bars,
            }
            out = dest / f"{prefix}_score_bars_mqc.json"
            out.write_text(json.dumps(bar_payload, indent=2) + "\n", encoding="utf-8")
            written.append(out)
    return written


def bundle_multiqc(output_dir: PathLike) -> Path:
    """Materialize ``*_mqc.*`` custom content from stage JSON + plot folders."""
    root = as_path(output_dir)
    write_overview_report(root)
    dest = multiqc_stage_dir(root)
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)

    overview_html = dest / "00_SamovaR_pipeline_mqc.html"
    parts = [
        "<p>SamovaR writes one JSON summary per pipeline stage "
        "(<code>.log/multiqc/*.samovar.json</code>), in the same spirit as FastQC/fastp. "
        "This MultiQC report includes native Plotly heatmaps, scatters, and bars "
        "(tick a plot to export PNG/SVG/PDF) plus Altair HTML for the same charts, "
        "and a short description of every pipeline step.</p><ol>"
    ]
    for i, stage in enumerate(CHECKPOINT_STEPS, start=1):
        info = STAGE_INFO[stage]
        json_path = report_dir(root) / JSON_NAME.format(stage=stage)
        status = "done" if json_path.is_file() else "not yet written"
        parts.append(
            f"<li><strong>{html.escape(info['title'])}</strong> "
            f"(<code>{html.escape(stage)}</code>, {status}) — "
            f"{html.escape(info['description'])}</li>"
        )
        plots_key = next((n for n in STAGE_DIRS.get(stage, ()) if n.endswith("_plots")), None)
        if plots_key:
            _copy_plot_assets(root / plots_key, dest, stage, i)
    parts.append("</ol>")
    overview_html.write_text(
        _html_comment(
            {
                "id": "samovar_overview",
                "section_name": "SamovaR pipeline",
                "description": "Stage map for this run, with Fastp-style JSON under .log/multiqc/.",
                "plot_type": "html",
            }
        )
        + "\n".join(parts)
        + "\n",
        encoding="utf-8",
    )
    config = {
        "title": "SamovaR",
        "subtitle": f"Ensemble annotation report (v{PACKAGE_VERSION})",
        "intro_text": (
            "Interactive MultiQC view of SamovaR: native Plotly plots (selectable "
            "export) plus Altair HTML, and a short description of each pipeline stage."
        ),
        "run_names": True,
    }
    order = ["samovar_overview"] + [f"samovar_{stage}" for stage in CHECKPOINT_STEPS]
    config_text = (
        "title: {title}\n"
        "subtitle: {subtitle}\n"
        "intro_text: |\n"
        "  {intro}\n"
        "custom_content:\n"
        "  order:\n{order}\n"
    ).format(
        title=json.dumps(config["title"]),
        subtitle=json.dumps(config["subtitle"]),
        intro=config["intro_text"],
        order="\n".join(f"    - {item}" for item in order),
    )
    (dest / "multiqc_config.yaml").write_text(config_text, encoding="utf-8")
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
