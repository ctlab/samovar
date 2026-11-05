#!/usr/bin/env python3
"""Install the explicitly listed, self-contained SAMOVAR skills (stdlib only)."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile


BUNDLE = Path(__file__).resolve().parent
DESTINATIONS = {
    "codex": ".agents/skills",
    "cursor": ".cursor/skills",
    "claude": ".claude/skills",
    "generic": "samovar-harness/skills",
}
MARKER = ".samovar-harness.json"


def detect_ide():
    explicit = os.environ.get("SAMOVAR_IDE", "")
    if explicit:
        if explicit not in DESTINATIONS:
            raise ValueError("SAMOVAR_IDE must be codex, cursor, claude, or generic")
        return explicit
    # Prefer the active agent over its hosting terminal (e.g. Codex in Cursor).
    if os.environ.get("CLAUDECODE"):
        return "claude"
    if os.environ.get("CODEX_THREAD_ID"):
        return "codex"
    if os.environ.get("CURSOR_AGENT") or os.environ.get("TERM_PROGRAM") == "cursor":
        return "cursor"
    return None


def digest(data):
    return hashlib.sha256(data).hexdigest()


def bundle_files():
    manifest = json.loads((BUNDLE / "bundle.json").read_text())
    if manifest.get("format_version") != 1 or not manifest.get("skills"):
        raise ValueError("Unsupported or empty SAMOVAR bundle manifest")
    result = {}
    for name in manifest["skills"]:
        if not isinstance(name, str) or not re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", name):
            raise ValueError("Invalid skill name in bundle manifest")
        if name in result:
            raise ValueError(f"Duplicate skill in bundle: {name}")
        source = BUNDLE / "skills" / name
        if not (source / "SKILL.md").is_file():
            raise ValueError(f"Missing SKILL.md for {name}")
        files = {}
        # Follow the source reference to shared rules; installed files are regular
        # copies so each skill can be moved without the source checkout.
        for directory, _, names in os.walk(source, followlinks=True):
            for filename in names:
                path = Path(directory) / filename
                if not path.resolve().is_relative_to(BUNDLE):
                    raise ValueError(f"Bundle reference escapes Agents/: {path}")
                relative = path.relative_to(source).as_posix()
                files[relative] = path.read_bytes()
        result[name] = files
    return result


def check_destination(destination):
    if destination.is_symlink():
        raise ValueError(f"Refusing to replace symlink: {destination}")
    if not destination.exists():
        return
    marker = destination / MARKER
    if not marker.is_file():
        raise ValueError(f"Existing skill is not managed by SAMOVAR: {destination}")
    previous = json.loads(marker.read_text())
    if previous.get("format_version") != 1:
        raise ValueError(f"Unsupported ownership marker: {marker}")
    expected = previous["files"]
    actual = {}
    for path in destination.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Locally changed skill (symlink): {path}")
        if path.is_file() and path != marker:
            actual[path.relative_to(destination).as_posix()] = digest(path.read_bytes())
    if actual != expected:
        raise ValueError(f"Locally changed skill; preserve your edits before updating: {destination}")


def install_skill(destination, files):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".samovar-harness-", dir=destination.parent) as temporary:
        staging = Path(temporary) / "new"
        staging.mkdir()
        for relative, data in files.items():
            path = staging / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        marker = {"format_version": 1, "files": {key: digest(value) for key, value in files.items()}}
        (staging / MARKER).write_text(json.dumps(marker, indent=2) + "\n")
        backup = Path(temporary) / "old"
        if destination.exists():
            destination.rename(backup)
        try:
            staging.rename(destination)
        except OSError:
            if backup.exists():
                backup.rename(destination)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ide", choices=["auto", *DESTINATIONS], default="auto")
    parser.add_argument("--project", type=Path, default=Path.cwd(), help="destination project (default: current directory)")
    parser.add_argument("--dry-run", action="store_true", help="validate and show destinations without writing")
    parser.add_argument("--mode", choices=["explicit", "full", "ask"], default="explicit", help="install.sh integration; full skips CI/unknown IDE, ask requires a terminal")
    args = parser.parse_args()
    ide = detect_ide() if args.ide == "auto" else args.ide
    if args.mode == "full" and (os.environ.get("CI") or ide is None):
        print("SAMOVAR harness skipped (CI or no detected IDE). Use skills_install.sh --ide NAME.")
        return
    if args.mode == "ask":
        if not sys.stdin.isatty():
            print("SAMOVAR harness skipped (non-interactive ask). Use skills_install.sh --ide NAME.")
            return
        answer = input(f"Install SAMOVAR harness? [codex/cursor/claude/generic/no] ({ide or 'no'}): ").strip().lower()
        if answer in {"n", "no"} or (not answer and ide is None):
            print("SAMOVAR harness skipped.")
            return
        ide = answer or ide
    if ide not in DESTINATIONS:
        raise ValueError("IDE not detected. Choose --ide codex, cursor, claude, or generic.")
    project = args.project.resolve()
    if not project.is_dir():
        raise ValueError(f"Project directory does not exist: {project}")
    skills = bundle_files()
    # Validate all destinations before changing any skill.
    for name, files in skills.items():
        check_destination(project / DESTINATIONS[ide] / name)
    for name, files in skills.items():
        destination = project / DESTINATIONS[ide] / name
        if not args.dry_run:
            install_skill(destination, files)
        print(f"{'Would install' if args.dry_run else 'Installed'} {name}: {destination}")
    print("Rules are bundled with each skill and apply only to SAMOVAR tasks. No lifecycle hooks enabled.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print(f"SAMOVAR harness error: {error}", file=sys.stderr)
        sys.exit(1)
