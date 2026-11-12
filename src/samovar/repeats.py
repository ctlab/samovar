"""Repeat-package contract: manifest and out_XXXX directories.

Per-repeat generation and execution are later stages. This module only
records the package. Seeds come from SHA-256 of the base seed, the
1-based repeat index and a namespace, so they do not depend on Python's
hash or on the order of execution.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

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


def guard_generate(args: Any) -> Optional[int]:
    """Reject a bad generate target before the generator writes anything."""
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


def prepare_package(config: Any) -> Optional[Path]:
    """Return the manifest when this prepare is a package. None keeps the single-run layout."""
    root = Path(config.output_dir)
    repeat = parse_repeat(getattr(config, "repeat", None))
    existing = load_manifest(root) if root.exists() else None
    if repeat is None:
        if existing:
            raise ValueError(
                f"{root} is a repeat package ({existing['repeat']} repeats); "
                "pass the same --repeat, or choose another output directory"
            )
        return None
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
    return publish_package(
        root,
        count=repeat,
        base_seed=int(getattr(config, "regeneration_seed", None) or 42),
        kind=kind,
        identity=ident,
        settings={"mode": getattr(config, "mode", "full")},
        adopt_without_identity=adopt,
    )
