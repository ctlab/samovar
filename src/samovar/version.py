"""Package version from ``pyproject.toml`` (single source of truth)."""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Optional

_VERSION_RE = re.compile(
    r'(?m)^version\s*=\s*["\']([^"\']+)["\']\s*$'
)


def pyproject_path() -> Path:
    """Repo ``pyproject.toml`` next to the package (``src/samovar`` → repo root)."""
    return Path(__file__).resolve().parent.parent.parent / "pyproject.toml"


def _version_from_pyproject(path: Path) -> Optional[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    # Prefer stdlib tomllib (3.11+); fall back to a line scan of [project].
    try:
        import tomllib
    except ImportError:  # pragma: no cover - Python 3.10
        tomllib = None  # type: ignore[assignment]
    if tomllib is not None:
        try:
            data = tomllib.loads(text)
            ver = (data.get("project") or {}).get("version")
            if ver:
                return str(ver).strip()
        except Exception:
            pass
    in_project = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_project = stripped == "[project]"
            continue
        if not in_project:
            continue
        match = _VERSION_RE.match(stripped)
        if match:
            return match.group(1).strip()
    return None


@lru_cache(maxsize=1)
def get_version() -> str:
    """Return the package version declared in ``pyproject.toml``.

    Falls back to ``importlib.metadata`` when the tree is installed without a
    nearby ``pyproject.toml``, then to ``"0.0.0"`` only as a last resort.
    """
    from_file = _version_from_pyproject(pyproject_path())
    if from_file:
        return from_file
    try:
        from importlib.metadata import PackageNotFoundError, version

        return version("samovar")
    except Exception:
        return "0.0.0"
