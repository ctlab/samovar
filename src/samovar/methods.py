"""Fill methods.md / samovar-method from prepare configs (no LLM prose).

Source of truth: checkpoint window, ``.log/configs/*.yaml``, Hydra prepare
args, ``used_citations.bib``, and ``cite/methods/*.txt`` sentence templates.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

from samovar.citations import (
    USED_CITATIONS_NAME,
    _atomic_write_text,
    bibtex_citekey,
    citation_refs_for_tool,
    parse_bibtex_entries,
    resolve_citation_file,
)
from samovar.exec_control import CHECKPOINT_STEPS, resolve_window, step_in_window
from samovar.paths import _code_repo_root, add_output_dir_argument

STAGE_CONTRACT: Dict[str, str] = {
    "setup_reads": "reads_generator",
    "qc_initial": "qc",
    "annotate_initial": "annotator",
    "viz_initial": "scoring",
    "regenerate_tables": "table_reads_generator",
    "score_sample_qc_full": "sample_scoring",
    "filter_sample_qc_full": "sample_filtering",
    "score_regenerated_tables": "table_scoring",
    "score_sample_qc_final": "sample_scoring",
    "filter_sample_qc_final": "sample_filtering",
    "regenerate_reads": "reads_generator",
    "qc_generated": "qc",
    "annotate_regenerated": "annotator",
    "viz_regenerated": "scoring",
    "reprofile": "reprofiler",
    "viz_reprofiled": "scoring",
}

INTERNAL_STAGES = (
    "combine_initial",
    "abundance_tables",
    "seed_genomes",
    "sort_reads",
    "combine_regenerated",
)

_SKIP = frozenset({"", "none", "off", "false", "0", "skip", "no"})
_ENV_ASSIGN = re.compile(r"^export\s+([A-Z0-9_]+)=(.*)$")
_SECTION = re.compile(r"^%\s*---\s*(.+?)\s*---\s*$")


def methods_template_dir() -> Path:
    return _code_repo_root() / "cite" / "methods"


def _load_yaml(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return {}
    return data if isinstance(data, dict) else {}


def _parse_window_env(path: Path) -> Tuple[str, str]:
    start, end = CHECKPOINT_STEPS[0], CHECKPOINT_STEPS[-1]
    if not path.is_file():
        return start, end
    for line in path.read_text(encoding="utf-8").splitlines():
        match = _ENV_ASSIGN.match(line.strip())
        if not match:
            continue
        key, raw = match.group(1), match.group(2).strip().strip('"').strip("'")
        if key == "SAMOVAR_START" and raw:
            start = raw
        elif key == "SAMOVAR_END" and raw:
            end = raw
    return resolve_window(start, end)


def _is_on(value: Any) -> bool:
    return str(value or "").strip().lower() not in _SKIP


def _as_list(value: Any) -> List[str]:
    if value in (None, False, ""):
        return []
    if isinstance(value, (list, tuple)):
        return [str(x).strip() for x in value if str(x).strip()]
    return [p for p in str(value).replace(",", " ").split() if p.strip()]


def join_en(parts: Sequence[str]) -> str:
    items = [str(p).strip() for p in parts if str(p).strip()]
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return f"{items[0]} and {items[1]}"
    return ", ".join(items[:-1]) + ", and " + items[-1]


def latex_escape(text: str) -> str:
    out = str(text or "")
    for src, dest in (
        ("\\", r"\textbackslash{}"),
        ("{", r"\{"),
        ("}", r"\}"),
        ("$", r"\$"),
        ("&", r"\&"),
        ("%", r"\%"),
        ("#", r"\#"),
        ("_", r"\_"),
        ("~", r"\textasciitilde{}"),
        ("^", r"\textasciicircum{}"),
    ):
        out = out.replace(src, dest)
    return out


def db_label(path: Any, name: Any = None) -> str:
    if name not in (None, "", False):
        return str(name).strip()
    text = str(path or "").strip()
    if not text or text == ".":
        return ""
    return Path(text.rstrip("/")).name


def parse_used_citation_keys(bib_text: str) -> Dict[str, List[str]]:
    """Map ``% --- tool ---`` sections in used_citations.bib to citekeys."""
    mapping: Dict[str, List[str]] = {}
    current = ""
    buf: List[str] = []

    def flush() -> None:
        if not current:
            return
        keys = mapping.setdefault(current, [])
        for entry in parse_bibtex_entries("\n".join(buf)):
            key = bibtex_citekey(entry)
            if key and key not in keys:
                keys.append(key)

    for line in str(bib_text or "").splitlines():
        head = _SECTION.match(line.strip())
        if head:
            flush()
            current = head.group(1).strip()
            buf = []
            continue
        buf.append(line)
    flush()
    return mapping


def citation_keys_for_tool(name: str, used_map: Mapping[str, Sequence[str]]) -> List[str]:
    keys: List[str] = []
    for key in used_map.get(name) or []:
        if key not in keys:
            keys.append(str(key))
    if keys:
        return keys
    for ref in citation_refs_for_tool(name):
        path = resolve_citation_file(ref)
        if path is None:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        for entry in parse_bibtex_entries(text):
            cite = bibtex_citekey(entry)
            if cite and cite not in keys:
                keys.append(cite)
    return keys


def load_template(contract: str, root: Optional[Path] = None) -> Dict[str, str]:
    dest = (root or methods_template_dir()) / f"{contract}.txt"
    out = {"contract": contract, "md": "", "tex": ""}
    if not dest.is_file():
        return out
    for line in dest.read_text(encoding="utf-8").splitlines():
        if ":" not in line or line.lstrip().startswith("#"):
            continue
        key, _, rest = line.partition(":")
        out[key.strip()] = rest.strip()
    return out


def _fill(template: str, fields: Mapping[str, str]) -> str:
    class _Map(dict):
        def __missing__(self, key: str) -> str:
            return ""

    try:
        return str(template).format_map(_Map(fields)).strip()
    except (ValueError, IndexError):
        return str(template).strip()


def _cite_md(keys: Sequence[str]) -> str:
    if not keys:
        return ""
    return " [" + "; ".join("@" + k for k in keys) + "]"


def _cite_tex(keys: Sequence[str]) -> str:
    if not keys:
        return ""
    return r" \cite{" + ",".join(keys) + "}"


def _tool_fields(
    tools: Sequence[Mapping[str, str]],
    keys: Sequence[str],
    *,
    mode: str = "",
    contract: str = "",
    stage: str = "",
) -> Dict[str, str]:
    names = [str(t.get("tool") or "") for t in tools if t.get("tool")]
    items: List[str] = []
    for row in tools:
        tool = str(row.get("tool") or "").strip()
        if not tool:
            continue
        db = str(row.get("database") or "").strip()
        if db:
            items.append(f"{tool} on the {db} database")
        else:
            items.append(tool)
    n = len(names) or len(items)
    return {
        "tool": names[0] if names else "",
        "tools": join_en(names),
        "items": join_en(items),
        "n": str(n),
        "verb": "was" if n == 1 else "were",
        "database": str((tools[0] if tools else {}).get("database") or ""),
        "mode": mode,
        "contract": contract,
        "stage": stage,
        "keys": "; ".join(keys),
        "cite_md": _cite_md(keys),
        "cite_tex": _cite_tex(keys),
        "cite": _cite_md(keys),
    }


def _annotators(init: Mapping[str, Any]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for item in init.get("run_config") or []:
        if not isinstance(item, dict):
            continue
        tool = str(item.get("type") or item.get("run_name") or item.get("cmd") or "").strip()
        if tool.endswith("-test"):
            tool = tool[: -len("-test")]
        rows.append(
            {
                "tool": tool,
                "database": db_label(item.get("db_path"), item.get("db_name")),
                "cmd": str(item.get("cmd") or ""),
                "extra": str(item.get("extra") or ""),
            }
        )
    return rows


def _assembly_rows(annotators: Sequence[Mapping[str, str]]) -> Dict[str, List[Dict[str, str]]]:
    from samovar.assembly_profiling import parse_slot_extra

    grouped: Dict[str, List[Dict[str, str]]] = {
        "assembler": [],
        "gene_caller": [],
        "binner": [],
        "binner_qc": [],
        "binner_combine": [],
        "mag_taxonomy": [],
        "aligner": [],
        "mag_quantifier": [],
        "taxon_quantifier": [],
        "read_assigner": [],
    }
    for row in annotators:
        if str(row.get("tool") or "").lower() not in {"assembly", "assembly_profiling"}:
            continue
        slots = parse_slot_extra(str(row.get("extra") or ""))
        if slots.get("assembler"):
            grouped["assembler"].append({"tool": str(slots["assembler"])})
        if slots.get("gene_caller"):
            grouped["gene_caller"].append({"tool": str(slots["gene_caller"])})
        for name in slots.get("binners") or []:
            grouped["binner"].append({"tool": str(name)})
        if slots.get("binner_qc"):
            grouped["binner_qc"].append({"tool": str(slots["binner_qc"])})
        if slots.get("binner_combine"):
            grouped["binner_combine"].append({"tool": str(slots["binner_combine"])})
        if slots.get("mag_taxonomy"):
            grouped["mag_taxonomy"].append({"tool": str(slots["mag_taxonomy"])})
        if slots.get("aligner"):
            grouped["aligner"].append({"tool": str(slots["aligner"])})
        if slots.get("mag_quantifier"):
            grouped["mag_quantifier"].append({"tool": str(slots["mag_quantifier"])})
        if slots.get("taxon_quantifier"):
            grouped["taxon_quantifier"].append({"tool": str(slots["taxon_quantifier"])})
        if slots.get("read_assigner"):
            grouped["read_assigner"].append({"tool": str(slots["read_assigner"])})
    return grouped


def collect_stage_rows(outdir: Path) -> List[Dict[str, Any]]:
    """Enabled stages with tools, contract, citation keys, filled sentences."""
    configs = outdir / ".log" / "configs"
    init = _load_yaml(configs / "config_init.yaml")
    a2iss = _load_yaml(configs / "config_annotation2iss.yaml")
    qc = _load_yaml(configs / "config_qc.yaml")
    scoring = _load_yaml(configs / "config_scoring.yaml")
    repro = _load_yaml(configs / "config_reprofiling.yaml")
    export = _load_yaml(configs / "config_export.yaml")
    start, end = _parse_window_env(outdir / ".log" / "window.env")
    bib_path = outdir / USED_CITATIONS_NAME
    if not bib_path.is_file():
        bib_path = configs / USED_CITATIONS_NAME
    bib_text = bib_path.read_text(encoding="utf-8") if bib_path.is_file() else ""
    used_map = parse_used_citation_keys(bib_text)

    annotators = _annotators(init)
    assembly = _assembly_rows(annotators)
    reads = str(a2iss.get("reads_generator") or "iss").strip() or "iss"
    meta = str(a2iss.get("metagenome_generator") or "").strip()
    table_modes = _as_list(a2iss.get("table_reads_generators") or a2iss.get("table_reads_generator") or a2iss.get("regeneration_mode"))
    table_score = str(a2iss.get("table_score") or "").strip()
    sample_score = str(a2iss.get("sample_score") or "").strip()
    sample_filter = str(a2iss.get("sample_filter") or "").strip()
    qc_initial = str(qc.get("qc_initial") or qc.get("qc") or "").strip()
    qc_generated = str(qc.get("qc_generated") or qc.get("qc") or "").strip()
    scoring_tools = _as_list(scoring.get("scoring_tools"))
    reprofiler = str(repro.get("reprofiler") or "ensemble").strip()
    fi = str(repro.get("feature_importance") or "").strip()
    export_name = str(export.get("export") or export.get("export_corrector") or "").strip()
    export_formats = _as_list(export.get("export_formats"))

    def keys_for(tools: Sequence[Mapping[str, str]]) -> List[str]:
        out: List[str] = []
        for row in tools:
            for key in citation_keys_for_tool(str(row.get("tool") or ""), used_map):
                if key not in out:
                    out.append(key)
        return out

    def pack(
        stage: str,
        contract: str,
        tools: Sequence[Mapping[str, str]],
        *,
        mode: str = "",
    ) -> Optional[Dict[str, Any]]:
        usable = [t for t in tools if _is_on(t.get("tool"))]
        if not usable:
            return None
        tmpl = load_template(contract)
        keys = keys_for(usable)
        fields = _tool_fields(usable, keys, mode=mode, contract=contract, stage=stage)
        return {
            "stage": stage,
            "contract": contract,
            "tools": usable,
            "keys": keys,
            "mode": mode,
            "md": _fill(tmpl.get("md") or "", fields),
            "tex": _fill(tmpl.get("tex") or "", fields),
            "methods": fields.get("items") or fields.get("tools") or "",
        }

    rows: List[Dict[str, Any]] = []
    for stage in CHECKPOINT_STEPS:
        if not step_in_window(stage, start, end):
            continue
        if stage in INTERNAL_STAGES:
            rows.append(
                {
                    "stage": stage,
                    "contract": "pipeline",
                    "tools": [],
                    "keys": [],
                    "mode": "",
                    "md": "",
                    "tex": "",
                    "methods": "",
                }
            )
            continue
        contract = STAGE_CONTRACT.get(stage, "")
        packed: Optional[Dict[str, Any]] = None
        if contract == "annotator":
            packed = pack(stage, "annotator", annotators)
        elif contract == "reads_generator":
            if meta and stage == "setup_reads":
                packed = pack(stage, "metagenome_generator", [{"tool": meta}])
            else:
                packed = pack(stage, "reads_generator", [{"tool": reads}])
        elif contract == "qc":
            name = qc_initial if stage == "qc_initial" else qc_generated
            packed = pack(stage, "qc", [{"tool": name}])
        elif contract == "table_reads_generator":
            packed = pack(
                stage,
                "table_reads_generator",
                [{"tool": m} for m in table_modes],
                mode=", ".join(table_modes),
            )
        elif contract == "table_scoring":
            packed = pack(stage, "table_scoring", [{"tool": table_score}])
        elif contract == "sample_scoring":
            packed = pack(stage, "sample_scoring", [{"tool": sample_score}])
        elif contract == "sample_filtering":
            packed = pack(stage, "sample_filtering", [{"tool": sample_filter}])
        elif contract == "scoring":
            packed = pack(stage, "scoring", [{"tool": n} for n in scoring_tools])
        elif contract == "reprofiler":
            packed = pack(stage, "reprofiler", [{"tool": reprofiler}])
        if packed:
            rows.append(packed)
        elif stage not in INTERNAL_STAGES:
            rows.append(
                {
                    "stage": stage,
                    "contract": contract or "pipeline",
                    "tools": [],
                    "keys": [],
                    "mode": "",
                    "md": "",
                    "tex": "",
                    "methods": "",
                }
            )
        if stage == "reprofile" and _is_on(fi):
            extra = pack("reprofile", "feature_importance", [{"tool": fi}])
            if extra:
                rows.append(extra)
        if stage == "annotate_initial":
            for sub, tools in assembly.items():
                extra = pack("annotate_initial", sub, tools)
                if extra:
                    rows.append(extra)

    if step_in_window("viz_reprofiled", start, end) and _is_on(export_name):
        extra = pack(
            "viz_reprofiled",
            "export",
            [{"tool": export_name}],
            mode=", ".join(export_formats) or "abundance",
        )
        if extra:
            rows.append(extra)

    if step_in_window("setup_reads", start, end):
        extra = pack("setup_reads", "workflow", [{"tool": "snakemake"}])
        if extra:
            rows.append(extra)

    return rows


def render_methods_md(rows: Sequence[Mapping[str, Any]]) -> str:
    lines = [
        "# Methods",
        "",
        "Generated from prepare contracts and `used_citations.bib` (no free-form text).",
        "",
        "## Stages",
        "",
    ]
    for row in rows:
        stage = row.get("stage") or ""
        contract = row.get("contract") or ""
        lines.append(f"### `{stage}` ({contract})")
        lines.append("")
        tools = row.get("tools") or []
        if tools:
            bits = []
            for tool in tools:
                name = str(tool.get("tool") or "")
                db = str(tool.get("database") or "")
                bits.append(f"{name} (database: {db})" if db else name)
            lines.append("- **Tool:** " + "; ".join(bits))
        else:
            lines.append("- **Tool:** —")
        keys = row.get("keys") or []
        lines.append("- **Citation keys:** " + (", ".join(f"`{k}`" for k in keys) if keys else "—"))
        lines.append("")
        sentence = str(row.get("md") or "").strip()
        if sentence:
            lines.append(sentence)
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_methods_tex(rows: Sequence[Mapping[str, Any]]) -> str:
    lines = [
        r"% Generated from prepare contracts and used_citations.bib",
        r"\section{Methods}",
        "",
    ]
    for row in rows:
        stage = latex_escape(str(row.get("stage") or ""))
        contract = latex_escape(str(row.get("contract") or ""))
        lines.append(rf"\subsection{{{stage} ({contract})}}")
        tools = row.get("tools") or []
        if tools:
            bits = []
            for tool in tools:
                name = latex_escape(str(tool.get("tool") or ""))
                db = str(tool.get("database") or "")
                if db:
                    bits.append(f"{name} (database: {latex_escape(db)})")
                else:
                    bits.append(name)
            lines.append("Tool: " + "; ".join(bits) + r".")
        keys = row.get("keys") or []
        if keys:
            lines.append("Citation keys: " + ", ".join(latex_escape(k) for k in keys) + ".")
            lines.append(r"\cite{" + ",".join(keys) + "}")
        sentence = str(row.get("tex") or "").strip()
        if sentence:
            lines.append(sentence)
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_methods_tsv(rows: Sequence[Mapping[str, Any]]) -> str:
    lines = ["contract\tmethods\tcitation"]
    seen = set()
    for row in rows:
        contract = str(row.get("contract") or "")
        methods = str(row.get("methods") or row.get("md") or "").strip()
        citation = "; ".join(row.get("keys") or [])
        if not methods or contract in {"pipeline", ""}:
            continue
        key = (contract, methods, citation)
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"{contract}\t{methods}\t{citation}")
    return "\n".join(lines) + "\n"


def write_methods(outdir: Path) -> Dict[str, str]:
    root = Path(outdir).expanduser()
    rows = collect_stage_rows(root)
    md = render_methods_md(rows)
    tex = render_methods_tex(rows)
    tsv = render_methods_tsv(rows)
    payload = {
        "stages": [
            {
                "stage": r.get("stage"),
                "contract": r.get("contract"),
                "tools": r.get("tools"),
                "keys": r.get("keys"),
            }
            for r in rows
        ]
    }
    written = {
        "methods.md": md,
        "samovar-method.md": md,
        "methods.tex": tex,
        "methods.tsv": tsv,
        "methods.json": json.dumps(payload, indent=2) + "\n",
    }
    configs = root / ".log" / "configs"
    configs.mkdir(parents=True, exist_ok=True)
    for name, text in written.items():
        _atomic_write_text(root / name, text)
        _atomic_write_text(configs / name, text)
    return {name: str(root / name) for name in written}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="samovar methods",
        description="Write methods.md from prepare configs and used_citations.bib",
    )
    add_output_dir_argument(parser, required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)
    out = Path(args.output_dir).expanduser()
    if not out.exists():
        print(f"Error: output dir not found: {out}", file=__import__("sys").stderr)
        return 1
    paths = write_methods(out)
    print(f"Wrote {paths['methods.md']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
