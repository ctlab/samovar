"""Repeat-package contract: manifest, sequential exec, and long summary tables.

Seeds come from SHA-256 of the base seed, the 1-based repeat index and a
namespace, so they do not depend on Python's hash or on the order of execution.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

SCHEMA = 1
SEED_RULE = "sha256-v1"
MAX_REPEAT = 9999
_CHILD = re.compile(r"out_(\d{4})$")


def parse_repeat(value: Any) -> Optional[int]:
    """Missing means a single run. An explicit count, including 1, is a package."""
    if value is None or value is False or value == "":
        return None
    if isinstance(value, bool) or isinstance(value, float) or not str(value).strip().isdigit():
        raise ValueError(f"repeat must be a positive integer, got {value!r}")
    count = int(str(value).strip())
    if count < 1 or count > MAX_REPEAT:
        raise ValueError(f"repeat must be from 1 to {MAX_REPEAT}, got {count}")
    return count


def repeat_seed(base: int, index: int, namespace: str) -> int:
    """Stable seed in ``1 .. 2**31-2``. ``index`` is 1-based."""
    payload = f"{int(base)}:{int(index)}:{namespace}".encode()
    value = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") % (2**31 - 1)
    return value or 1


def manifest_path(root: Path) -> Path:
    return Path(root) / ".log" / "repeats.json"


def load_manifest(root: Path) -> Optional[Dict[str, Any]]:
    path = manifest_path(root)
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema") != SCHEMA:
        raise ValueError(f"unsupported repeat manifest: {path}")
    return data


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def _parent_package(root: Path) -> Optional[Path]:
    if _CHILD.match(root.name) and manifest_path(root.parent).is_file():
        return root.parent
    return None


def _occupied(root: Path, *, allow_template: bool = False) -> bool:
    if not root.exists():
        return False
    for path in root.iterdir():
        if allow_template and path.name in {".generate", ".log"}:
            continue
        return True
    return False


def _rows(count: int, base_seed: int) -> list:
    rows = []
    for index in range(1, count + 1):
        rows.append(
            {
                "repeat_id": f"{index:04d}",
                "out_dir": f"out_{index:04d}",
                "seed": repeat_seed(base_seed, index, "repeat"),
                "state": "planned",
                "stages": {},
            }
        )
    return rows


def _same(existing: Mapping[str, Any], request: Mapping[str, Any], *, identity_required: bool) -> None:
    if int(existing["repeat"]) != int(request["repeat"]) or int(existing["base_seed"]) != int(request["base_seed"]):
        raise ValueError(
            "repeat package does not match this request: "
            f"manifest repeat={existing['repeat']} seed={existing['base_seed']}, "
            f"request repeat={request['repeat']} seed={request['base_seed']}"
        )
    old = existing.get("input") or {}
    new = request["input"]
    if identity_required and (old.get("kind") != new.get("kind") or old.get("identity") != new.get("identity")):
        raise ValueError(
            "repeat package input does not match this request: "
            f"{old.get('kind')} {old.get('identity')} vs {new.get('kind')} {new.get('identity')}"
        )
    if not identity_required and new.get("kind") and (
        old.get("kind") != new.get("kind") or old.get("identity") != new.get("identity")
    ):
        raise ValueError(
            "repeat package input does not match this request: "
            f"{old.get('kind')} vs {new.get('kind')}"
        )


def _ensure_children(root: Path, rows: list) -> None:
    for row in rows:
        (root / row["out_dir"]).mkdir(parents=True, exist_ok=True)


def publish_package(
    root: Path,
    *,
    count: int,
    base_seed: int,
    kind: str,
    identity: str,
    settings: Optional[Mapping[str, Any]] = None,
    adopt_without_identity: bool = False,
    allow_template: bool = False,
) -> Path:
    """Create or resume a package. Raises before replacing an incompatible manifest."""
    root = Path(root)
    if _parent_package(root):
        raise ValueError(f"{root} is already a repeat directory; refusing to nest another package")
    request = {
        "repeat": int(count),
        "base_seed": int(base_seed),
        "input": {"kind": kind, "identity": identity},
    }
    existing = load_manifest(root)
    if existing:
        _same(existing, request, identity_required=not adopt_without_identity)
        _ensure_children(root, existing["repeats"])
        return manifest_path(root)
    if _occupied(root, allow_template=allow_template):
        raise ValueError(f"{root} is not empty and has no repeat manifest; refusing to overwrite")
    payload = {
        "schema": SCHEMA,
        "seed_rule": SEED_RULE,
        "repeat": int(count),
        "base_seed": int(base_seed),
        "input": {"kind": kind, "identity": identity},
        "settings": dict(settings or {}),
        "repeats": _rows(count, int(base_seed)),
    }
    _atomic_json(manifest_path(root), payload)
    _ensure_children(root, payload["repeats"])
    return manifest_path(root)


def generated_identity(args: Any) -> str:
    genome = getattr(args, "genome_dir", None)
    genome_text = str(Path(genome).resolve()) if genome else ""
    abundance = getattr(args, "abundance_table", None)
    abundance_text = str(Path(abundance).resolve()) if abundance else ""
    simulator = str(getattr(args, "reads_generator", None) or getattr(args, "simulator", None) or "iss")
    meta = str(getattr(args, "metagenome_generator", None) or "")
    return "|".join(
        [
            genome_text,
            simulator,
            meta,
            str(getattr(args, "n_samples", None)),
            str(getattr(args, "total_reads", None)),
            abundance_text,
        ]
    )


def prepare_identity(config: Any) -> Optional[tuple]:
    """Reads or profiles can start a package. Anything else must already be generated."""
    parts = []
    kind = None
    if getattr(config, "input_dir", None):
        kind = "reads"
        parts.append("reads:" + str(Path(config.input_dir).resolve()))
    profiles = getattr(config, "profile_inputs", None) or {}
    if profiles:
        kind = kind or "profiles"
        named = []
        for name in sorted(profiles):
            value = profiles[name]
            text = str(Path(value).resolve()) if isinstance(value, (str, Path)) else json.dumps(value, sort_keys=True)
            named.append(f"{name}={text}")
        parts.append("profiles:" + ",".join(named))
    if kind is None:
        return None
    return kind, "|".join(parts)


def record_generator_plan(root: Path, planned: list) -> None:
    """Remember the seed each existing generator will receive. Does not run it."""
    data = load_manifest(root)
    if data is None:
        raise ValueError(f"repeat manifest is missing: {root}")
    by_id = {row["repeat_id"]: row for row in data["repeats"]}
    for item in planned:
        row = by_id[item["repeat_id"]]
        row["generator_seed"] = int(item["generator_seed"])
        row["generator_namespace"] = item["namespace"]
        row["abundance"] = item.get("abundance") or ""
        row["state"] = "planned"
    _atomic_json(manifest_path(root), data)


def set_repeat_state(
    root: Path,
    repeat_id: str,
    state: str,
    error: str = "",
    stage: str = "generate",
) -> None:
    """Record one stage. Does not create reads or truth."""
    data = load_manifest(root)
    if data is None:
        raise ValueError(f"repeat manifest is missing: {root}")
    found = False
    for row in data["repeats"]:
        if row["repeat_id"] != str(repeat_id):
            continue
        stages = dict(row.get("stages") or {})
        stages[stage] = state
        row["stages"] = stages
        if stage == "generate":
            row["state"] = state
        if error:
            row["error"] = error
        elif stage == "pipeline" and state == "succeeded":
            row.pop("error", None)
        found = True
    if not found:
        raise ValueError(f"repeat {repeat_id} is not in {manifest_path(root)}")
    _atomic_json(manifest_path(root), data)


def write_repeat_orchestrator(root: Path, rows: list) -> None:
    """Root generate.sh runs each child's existing generate.sh and keeps going."""
    from samovar.paths import python_path, repo_root, shell_source_install_env_snippet

    root = Path(root)
    lines = [
        "#!/bin/bash",
        "set -u",
        f'export SAMOVAR_ROOT="{repo_root()}"',
        shell_source_install_env_snippet().rstrip(),
        'export PYTHONPATH="' + str(repo_root() / "src") + '${PYTHONPATH:+:$PYTHONPATH}"',
        f'PYTHON_PATH="${{PYTHON_PATH:-{python_path()}}}"',
        f'root="{root.resolve()}"',
        "fail=0",
        "mark()",
        "{",
        '  "$PYTHON_PATH" -m samovar.repeats mark --root "$root" --repeat-id "$1" --state "$2" --error "${3:-}"',
        "}",
    ]
    for row in rows:
        child = root / row["out_dir"]
        lines.extend([
            f'if bash "{child.resolve() / ".generate" / "generate.sh"}" --directory "{child.resolve()}"; then',
            f'  mark {row["repeat_id"]} generated',
            "else",
            f'  mark {row["repeat_id"]} error "generate failed"',
            "  fail=1",
            "fi",
        ])
    lines.append('exit "$fail"')
    script = root / ".generate" / "generate.sh"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.chmod(script, 0o755)


def materialize_repeat_generations(args: Any) -> None:
    """Write one existing generator workflow per repeat. Does not execute it."""
    import argparse

    from samovar.iss_config import setup_iss_test

    root = Path(args.output_dir)
    manifest = load_manifest(root)
    if manifest is None:
        raise ValueError(f"repeat manifest is missing: {root}")
    abundance = str(getattr(args, "abundance_table", None) or "")
    namespace = "reads" if abundance else "community"
    base = int(getattr(args, "seed", None) if getattr(args, "seed", None) is not None else 42)
    planned = []
    for row in manifest["repeats"]:
        index = int(row["repeat_id"])
        seed = repeat_seed(base, index, namespace)
        child = root / row["out_dir"]
        child_args = argparse.Namespace(**vars(args))
        child_args.output_dir = str(child)
        child_args.seed = seed
        child_args.repeat = None
        child_args._repeat_count = None
        child_args._writing_repeat_child = True
        setup_iss_test(child_args)
        planned.append({
            "repeat_id": row["repeat_id"],
            "generator_seed": seed,
            "namespace": namespace,
            "abundance": str(Path(abundance).resolve()) if abundance else "",
        })
    record_generator_plan(root, planned)
    write_repeat_orchestrator(root, manifest["repeats"])


def guard_generate(args: Any) -> Optional[int]:
    """Reject a bad generate target before the generator writes anything."""
    if getattr(args, "_writing_repeat_child", False):
        return None
    count = parse_repeat(getattr(args, "repeat", None))
    root = Path(getattr(args, "output_dir", ""))
    if _parent_package(root):
        raise ValueError(f"{root} is already a repeat directory; refusing to nest another package")
    existing = load_manifest(root) if root.exists() else None
    if count is None:
        if existing:
            raise ValueError(
                f"{root} is a repeat package ({existing['repeat']} repeats); "
                "pass the same --repeat, or choose another output directory"
            )
        return None
    if existing:
        publish_package(
            root,
            count=count,
            base_seed=int(getattr(args, "seed", None) or 42),
            kind="generated",
            identity=generated_identity(args),
        )
    elif _occupied(root):
        raise ValueError(f"{root} is not empty and has no repeat manifest; refusing to overwrite")
    return count


def main(argv: Optional[list] = None) -> None:
    import argparse
    import sys

    parser = argparse.ArgumentParser(prog="samovar.repeats")
    sub = parser.add_subparsers(dest="cmd", required=True)
    mark = sub.add_parser("mark")
    mark.add_argument("--root", required=True)
    mark.add_argument("--repeat-id", required=True)
    mark.add_argument("--state", required=True)
    mark.add_argument("--error", default="")
    mark.add_argument("--stage", default="generate")
    run = sub.add_parser("exec")
    run.add_argument("--root", required=True)
    run.add_argument("--redo", action="store_true")
    run.add_argument("--startpoint", default=None)
    run.add_argument("--endpoint", default=None)
    tables = sub.add_parser("tables")
    tables.add_argument("--root", required=True)
    args = parser.parse_args(argv)
    if args.cmd == "mark":
        try:
            set_repeat_state(Path(args.root), args.repeat_id, args.state, args.error, args.stage)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            raise SystemExit(1) from exc
    elif args.cmd == "exec":
        raise SystemExit(exec_package(
            Path(args.root),
            redo=args.redo,
            startpoint=args.startpoint,
            endpoint=args.endpoint,
        ))
    elif args.cmd == "tables":
        try:
            write_long_tables(Path(args.root))
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            raise SystemExit(1) from exc


def _supplied_repeat_blocked(kind: str, count: int) -> None:
    if count > 1 and kind in {"reads", "profiles"}:
        raise ValueError(
            "resampling existing FASTQ or a saved profile is not supported yet; "
            "--repeat 1 prepares one run and does not create independent replicates"
        )


def write_status_summary(root: Path) -> Path:
    """Counts only. A partial package is incomplete; none succeeded is no data."""
    data = load_manifest(root)
    if data is None:
        raise ValueError(f"repeat manifest is missing: {root}")
    errors = []
    success = 0
    for row in data["repeats"]:
        stages = row.get("stages") or {}
        if stages.get("pipeline") == "succeeded":
            success += 1
            continue
        if stages.get("pipeline") == "failed" or stages.get("generate") in {"error", "failed"} or row.get("error"):
            errors.append({
                "repeat_id": row["repeat_id"],
                "error": row.get("error") or stages.get("pipeline") or stages.get("generate") or "failed",
            })
    expected = len(data["repeats"])
    if success == 0:
        status = "no data"
    elif success < expected:
        status = "incomplete"
    else:
        status = "complete"
    path = Path(root) / "summary" / "status.json"
    _atomic_json(path, {
        "status": status,
        "n_expected": expected,
        "n_success": success,
        "errors": errors,
    })
    return path


_PHASES = ("initial", "regenerated", "reprofiled")
_DESC = ("kind", "taxonomy", "rank", "unit", "basis", "normalization")
_PROFILE_FIELDS = (
    "repeat_id", "sample_id", "seed", "phase", "run_name", *_DESC, "taxon", "value",
)
_RESULT_FIELDS = (
    "repeat_id", "sample_id", "seed", "phase", "run_name", *_DESC,
    "metric", "value", "truth_source", "reason",
)
_SCORE_IDS = {
    "annotator", "ground_truth", "truth_source", "sample", "sample_id", "reason",
    "phase", "run_name", *_DESC,
}


def _atomic_csv(path: Path, rows: List[Mapping[str, Any]], fields: tuple) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value != value:
        return ""
    return str(value)


def _blank_desc() -> Dict[str, str]:
    return {key: "" for key in _DESC}


def _desc(info: Mapping[str, Any]) -> Dict[str, str]:
    return {key: _cell(info.get(key, "")) for key in _DESC}


def _score_rows(rows: List[dict], child: Path, phase: str, repeat_id: str, seed: Any) -> None:
    plots = child / f"{phase}_annotations_plots"
    for name in ("profile_quality_scores.csv", "quality_scores.csv"):
        path = plots / name
        if not path.is_file():
            continue
        with path.open(encoding="utf-8", newline="") as handle:
            for rec in csv.DictReader(handle):
                reason = rec.get("reason") or ""
                base = {
                    "repeat_id": repeat_id,
                    "sample_id": rec.get("sample") or rec.get("sample_id") or "",
                    "seed": seed,
                    "phase": phase,
                    "run_name": rec.get("annotator") or rec.get("run_name") or "",
                    "truth_source": rec.get("ground_truth") or rec.get("truth_source") or "",
                    "reason": reason,
                    **_desc(rec),
                }
                metrics = []
                for key, raw in rec.items():
                    if not key or key in _SCORE_IDS or raw is None or str(raw).strip() == "":
                        continue
                    try:
                        value = float(raw)
                    except ValueError:
                        continue
                    if value != value:
                        continue
                    metrics.append((key, value))
                if not metrics:
                    if reason:
                        rows.append({**base, "metric": "", "value": ""})
                    continue
                for metric, value in metrics:
                    rows.append({**base, "metric": metric, "value": value})


def _runtime_rows(rows: List[dict], child: Path, repeat_id: str, seed: Any) -> None:
    from samovar.runtime import load_sample_runtimes, load_stage_runtimes

    for rec in load_sample_runtimes(child):
        for metric in ("elapsed_s", "max_rss_mb"):
            if rec.get(metric) is None:
                continue
            rows.append({
                "repeat_id": repeat_id,
                "sample_id": rec.get("sample") or "",
                "seed": seed,
                "phase": rec.get("phase") or "",
                "run_name": rec.get("run_name") or "",
                "metric": metric,
                "value": rec[metric],
                "truth_source": "",
                "reason": "",
                **_blank_desc(),
            })
    for rec in load_stage_runtimes(child):
        rows.append({
            "repeat_id": repeat_id,
            "sample_id": "",
            "seed": seed,
            "phase": "",
            "run_name": "",
            "metric": "stage_elapsed_s",
            "value": rec["elapsed_s"],
            "truth_source": "",
            "reason": rec.get("stage") or "",
            **_blank_desc(),
        })


def write_long_tables(root: Path) -> Dict[str, Path]:
    """Rewrite summary long tables from successful child artifacts.

    A missing taxon is 0 only inside one complete abundance namespace.
    Failed repeats, missing samples and missing tools stay absent.
    """
    root = Path(root)
    data = load_manifest(root)
    if data is None:
        raise ValueError(f"repeat manifest is missing: {root}")
    from samovar.abundance import n_sample_columns
    from samovar.profiles import read_profiles

    groups: Dict[tuple, List[dict]] = {}
    profile_rows: List[dict] = []
    result_rows: List[dict] = []
    for row in data["repeats"]:
        if (row.get("stages") or {}).get("pipeline") != "succeeded":
            continue
        child = root / row["out_dir"]
        seed = row.get("seed", "")
        repeat_id = row["repeat_id"]
        for phase in _PHASES:
            tables, metas = read_profiles(child / f"{phase}_profiles")
            for name, table in tables.items():
                info = metas[name]
                key = (phase, name) + tuple(_desc(info)[item] for item in _DESC)
                samples = {}
                for col in n_sample_columns(table):
                    samples[col[2:]] = {
                        str(taxid): value
                        for taxid, value in zip(table["taxid"].astype(str), table[col])
                    }
                groups.setdefault(key, []).append({
                    "repeat_id": repeat_id,
                    "seed": seed,
                    "info": info,
                    "samples": samples,
                })
            _score_rows(result_rows, child, phase, repeat_id, seed)
        _runtime_rows(result_rows, child, repeat_id, seed)

    for key, observations in groups.items():
        phase, run_name = key[0], key[1]
        kind = key[2]
        if kind == "abundance":
            taxa = sorted({
                taxon
                for item in observations
                for values in item["samples"].values()
                for taxon in values
                if taxon != "0" and not str(taxon).startswith("unmapped:")
            })
        else:
            taxa = None
        for item in observations:
            info = _desc(item["info"])
            for sample, values in sorted(item["samples"].items()):
                listed = taxa if taxa is not None else sorted(values)
                for taxon in listed:
                    if taxon not in values and taxa is None:
                        continue
                    profile_rows.append({
                        "repeat_id": item["repeat_id"],
                        "sample_id": sample,
                        "seed": item["seed"],
                        "phase": phase,
                        "run_name": run_name,
                        "taxon": taxon,
                        "value": values.get(taxon, 0),
                        **info,
                    })
                if taxa is not None:
                    for taxon, value in values.items():
                        if taxon == "0" or str(taxon).startswith("unmapped:"):
                            profile_rows.append({
                                "repeat_id": item["repeat_id"],
                                "sample_id": sample,
                                "seed": item["seed"],
                                "phase": phase,
                                "run_name": run_name,
                                "taxon": taxon,
                                "value": value,
                                **info,
                            })
    descriptors: Dict[tuple, tuple] = {}
    for key in groups:
        descriptors.setdefault((key[0], key[1]), key[2:])
    for item in result_rows:
        if item.get("kind") or item["metric"] in {"elapsed_s", "max_rss_mb", "stage_elapsed_s"}:
            continue
        desc = descriptors.get((item["phase"], item["run_name"]))
        if desc is None:
            continue
        item.update(dict(zip(_DESC, desc)))
    profile_rows.sort(key=lambda item: (
        item["repeat_id"], item["phase"], item["run_name"], item["sample_id"], item["taxon"],
    ))
    result_rows.sort(key=lambda item: (
        item["repeat_id"], item["phase"], item["run_name"], item["sample_id"],
        item["metric"], item["truth_source"], item["reason"],
    ))
    summary = root / "summary"
    profiles = summary / "profiles_long.csv"
    results = summary / "results_long.csv"
    _atomic_csv(profiles, profile_rows, _PROFILE_FIELDS)
    _atomic_csv(results, result_rows, _RESULT_FIELDS)
    _atomic_json(summary / "tables.json", {
        "schema": 1,
        "built_from": "saved child artifacts",
        "profiles": "profiles_long.csv",
        "results": "results_long.csv",
        "profile_columns": list(_PROFILE_FIELDS),
        "result_columns": list(_RESULT_FIELDS),
        "zero_fill": "missing taxon is 0 only inside one complete abundance namespace",
    })
    return {"profiles": profiles, "results": results}


def exec_package(
    root: Path,
    *,
    redo: bool = False,
    startpoint: Optional[str] = None,
    endpoint: Optional[str] = None,
) -> int:
    """Run each child through the existing exec, one at a time."""
    import subprocess

    from samovar.paths import repo_root

    root = Path(root)
    manifest = load_manifest(root)
    if manifest is None:
        raise ValueError(f"repeat manifest is missing: {root}")
    from samovar.seqio import has_r1_reads

    failed = 0
    env_base = dict(os.environ)
    exec_script = repo_root() / "bin" / "samovar_exec"
    for row in manifest["repeats"]:
        stages = row.get("stages") or {}
        generate_state = stages.get("generate") or row.get("state")
        child = root / row["out_dir"]
        child_generate = child / ".generate" / "generate.sh"
        if (
            child_generate.is_file()
            and generate_state not in {"generated", "succeeded"}
            and not has_r1_reads(child / "initial")
        ):
            set_repeat_state(root, row["repeat_id"], "running", stage="generate")
            generated = subprocess.run(
                ["bash", str(child_generate), "--directory", str(child)],
                check=False,
            )
            if generated.returncode != 0:
                set_repeat_state(root, row["repeat_id"], "failed", "generate failed", "generate")
                set_repeat_state(root, row["repeat_id"], "failed", "generate failed", "pipeline")
                failed += 1
                continue
            set_repeat_state(root, row["repeat_id"], "succeeded", stage="generate")
            generate_state = "succeeded"
        if generate_state in {"error", "failed"}:
            set_repeat_state(root, row["repeat_id"], "failed", row.get("error") or "generate failed", "pipeline")
            failed += 1
            continue
        if stages.get("pipeline") == "succeeded" and not redo:
            continue
        set_repeat_state(root, row["repeat_id"], "running", stage="pipeline")
        cmd = ["bash", str(exec_script), "--directory", str(child)]
        if redo:
            cmd.append("--redo")
        if startpoint:
            cmd.extend(["--startpoint", startpoint])
        if endpoint:
            cmd.extend(["--endpoint", endpoint])
        result = subprocess.run(cmd, check=False, env=env_base)
        if result.returncode == 0:
            set_repeat_state(root, row["repeat_id"], "succeeded", stage="pipeline")
        else:
            set_repeat_state(
                root, row["repeat_id"], "failed",
                f"exec exited {result.returncode}", "pipeline",
            )
            failed += 1
    write_status_summary(root)
    write_long_tables(root)
    return 1 if failed else 0


def prepare_package(config: Any) -> Optional[Path]:
    """Return the manifest when this prepare is a package. None keeps the single-run layout."""
    if getattr(config, "_repeat_child", False):
        return None
    root = Path(config.output_dir)
    repeat = parse_repeat(getattr(config, "repeat", None))
    existing = load_manifest(root) if root.exists() else None
    if repeat is None:
        if existing is None:
            return None
        repeat = int(existing["repeat"])
    elif existing is not None and int(existing["repeat"]) != int(repeat):
        raise ValueError(
            "repeat package does not match this request: "
            f"manifest repeat={existing['repeat']}, request repeat={repeat}"
        )
    identity = prepare_identity(config)
    own_initial = (root / "initial").resolve()
    generated_here = (
        existing is not None
        and (existing.get("input") or {}).get("kind") == "generated"
        and (
            identity is None
            or (
                getattr(config, "input_dir", None)
                and Path(config.input_dir).resolve() == own_initial
                and not getattr(config, "profile_inputs", None)
            )
        )
    )
    if generated_here or identity is None:
        if existing is None:
            raise ValueError(
                "prepare creates a repeat package only for reads (--input_dir) or "
                "profiles (profile_inputs); a generated package must come from samovar generate"
            )
        kind, ident = existing["input"]["kind"], existing["input"]["identity"]
        adopt = True
    else:
        kind, ident = identity
        adopt = False
    _supplied_repeat_blocked(kind, int(repeat))
    return publish_package(
        root,
        count=repeat,
        base_seed=int(getattr(config, "regeneration_seed", None) or 42),
        kind=kind,
        identity=ident,
        settings={"mode": getattr(config, "mode", "full")},
        adopt_without_identity=adopt,
    )


if __name__ == "__main__":
    main()
