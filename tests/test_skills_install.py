"""Portable harness installation must not alter the analysis or user settings."""

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def harness(tmp_path):
    # Exactly the documented publication boundary, with no installed package.
    source = tmp_path / "standalone skills repo"
    source.mkdir()
    shutil.copy2(ROOT / "skills_install.sh", source)
    shutil.copytree(ROOT / "Agents", source / "Agents", symlinks=True)
    project = tmp_path / "analysis with spaces"
    project.mkdir()
    env = os.environ.copy()
    for key in ("SAMOVAR_IDE", "CLAUDECODE", "CODEX_THREAD_ID", "CURSOR_AGENT", "TERM_PROGRAM", "CI"):
        env.pop(key, None)
    env["SAMOVAR_HARNESS_PYTHON"] = sys.executable
    return source, project, env


def run_install(harness, *args):
    source, project, env = harness
    return subprocess.run(
        ["bash", str(source / "skills_install.sh"), *args],
        cwd=project, env=env, input="", text=True, capture_output=True,
    )


@pytest.mark.parametrize("ide,directory", [
    ("codex", ".agents/skills"), ("cursor", ".cursor/skills"),
    ("claude", ".claude/skills"), ("generic", "samovar-harness/skills"),
])
def test_self_contained_bundle_and_repeat_install(harness, ide, directory):
    source, project, _ = harness
    unrelated = project / directory / "unrelated" / "SKILL.md"
    unrelated.parent.mkdir(parents=True)
    unrelated.write_text("user skill")
    result = run_install(harness, "--ide", ide)
    assert result.returncode == 0, result.stderr
    skill = project / directory / "samovar"
    references = skill / "references/rules"
    assert references.is_dir() and not references.is_symlink()
    for rule in (source / "Agents/rules").glob("*.md"):
        assert (references / rule.name).read_bytes() == rule.read_bytes()
    manifest = json.loads((source / "Agents/bundle.json").read_text())
    for name in manifest["skills"]:
        installed = project / directory / name
        assert (installed / "SKILL.md").is_file()
        for target in re.findall(r"\]\(([^)]+)\)", (installed / "SKILL.md").read_text()):
            if "://" not in target:
                assert (installed / target).is_file(), (name, target)
    assert (skill / "SKILL.md").is_file()
    before = {str(path.relative_to(project)): path.read_bytes() for path in project.rglob("*") if path.is_file()}
    assert run_install(harness, "--ide", ide).returncode == 0
    after = {str(path.relative_to(project)): path.read_bytes() for path in project.rglob("*") if path.is_file()}
    assert before == after
    assert unrelated.read_text() == "user skill"
    assert not (project / "AGENTS.md").exists()
    assert not (project / "CLAUDE.md").exists()
    assert not (project / ".claude/settings.json").exists()
    assert not (project / ".cursor/hooks.json").exists()
    # The copied references still work after removing the source checkout.
    shutil.rmtree(source)
    assert "Conditional" in (references / "scientific-integrity-samovar.md").read_text()


def test_update_unmodified_bundle_preserves_foreign_skills(harness):
    source, project, _ = harness
    assert run_install(harness, "--ide", "generic").returncode == 0
    source_skill = source / "Agents/skills/samovar/SKILL.md"
    source_skill.write_text(source_skill.read_text() + "\nUpdated guidance.\n")
    assert run_install(harness, "--ide", "generic").returncode == 0
    assert (project / "samovar-harness/skills/samovar/SKILL.md").read_bytes() == source_skill.read_bytes()


@pytest.mark.parametrize("owned", [False, True])
def test_refuses_foreign_or_locally_changed_skill(harness, owned):
    _, project, _ = harness
    skill = project / ".agents/skills/samovar"
    if owned:
        assert run_install(harness, "--ide", "codex").returncode == 0
    else:
        skill.mkdir(parents=True)
    (skill / "SKILL.md").write_text("my local instructions")
    result = run_install(harness, "--ide", "codex")
    assert result.returncode == 1
    assert "error:" in result.stderr
    assert (skill / "SKILL.md").read_text() == "my local instructions"


def test_dry_run_and_explicit_project(harness, tmp_path):
    _, project, _ = harness
    target = tmp_path / "other project"
    target.mkdir()
    result = run_install(harness, "--ide", "claude", "--project", str(target), "--dry-run")
    assert result.returncode == 0, result.stderr
    assert "Would install" in result.stdout
    assert not list(project.iterdir()) and not list(target.iterdir())
    assert run_install(harness, "--ide", "claude", "--project", str(target)).returncode == 0
    assert (target / ".claude/skills/samovar/SKILL.md").is_file()
    assert not list(project.iterdir())


@pytest.mark.parametrize("variable,value,ide,directory", [
    ("CODEX_THREAD_ID", "test-session", "codex", ".agents/skills"),
    ("CLAUDECODE", "1", "claude", ".claude/skills"),
    ("CURSOR_AGENT", "1", "cursor", ".cursor/skills"),
    ("TERM_PROGRAM", "cursor", "cursor", ".cursor/skills"),
    ("SAMOVAR_IDE", "generic", "generic", "samovar-harness/skills"),
])
def test_full_detects_active_agent(harness, variable, value, ide, directory):
    _, project, env = harness
    env[variable] = value
    result = run_install(harness, "--mode", "full")
    assert result.returncode == 0, result.stderr
    assert (project / directory / "samovar/SKILL.md").is_file()


def test_no_ide_explicit_fails_full_skips(harness):
    _, project, _ = harness
    assert run_install(harness).returncode == 1
    result = run_install(harness, "--mode", "full")
    assert result.returncode == 0 and "skipped" in result.stdout
    assert not list(project.iterdir())


def test_ci_and_noninteractive_ask_skip(harness):
    _, project, env = harness
    env.update(CODEX_THREAD_ID="test", CI="true")
    assert "skipped" in run_install(harness, "--mode", "full").stdout
    assert "non-interactive ask" in run_install(harness, "--mode", "ask").stdout
    assert not list(project.iterdir())


@pytest.mark.parametrize("answer,target", [("claude", ".claude"), ("no", None), ("", None)])
def test_interactive_ask(harness, monkeypatch, answer, target):
    # Run the real main with a terminal and deterministic user input.
    import importlib.util
    source, project, env = harness
    spec = importlib.util.spec_from_file_location("samovar_harness_install", source / "Agents/install.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for key in ("CODEX_THREAD_ID", "CLAUDECODE", "CURSOR_AGENT", "TERM_PROGRAM", "SAMOVAR_IDE"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(sys, "argv", ["install.py", "--mode", "ask", "--project", str(project)])
    monkeypatch.setattr(sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: answer)
    module.main()
    if target:
        assert (project / target / "skills/samovar/SKILL.md").is_file()
    else:
        assert not list(project.iterdir())


@pytest.mark.parametrize("mode,skip", [("full", False), ("ask", False), ("plain", True), ("disabled", True)])
def test_install_sh_harness_wiring(harness, mode, skip):
    source, project, env = harness
    # Exercise the install.sh function without reinstalling Python/bioinformatics
    # tools. This is the actual function invoked at the end of full/ask.
    script = (ROOT / "install.sh").read_text()
    function = script.split("install_agent_harness() {", 1)[1].split("\n}\n", 1)[0]
    env.update(
        ROOT=str(source), INSTALL_CALLER_DIR=str(project), PYTHON_PATH=sys.executable,
        INSTALL_FULL="1" if mode in {"full", "disabled"} else "0",
        INSTALL_ASK="1" if mode == "ask" else "0",
        SAMOVAR_INSTALL_HARNESS="0" if mode == "disabled" else "auto",
        SAMOVAR_IDE="codex",
    )
    result = subprocess.run(
        ["bash", "-c", "install_agent_harness() {" + function + "\n}\ninstall_agent_harness"],
        cwd=project, env=env, input="", text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    installed = (project / ".agents/skills/samovar/SKILL.md").is_file()
    assert installed == (not skip and mode != "ask")  # non-TTY ask must not install


def test_unlisted_skills_are_not_distributed(harness):
    source, project, _ = harness
    private = source / "Agents/skills/private"
    private.mkdir()
    (private / "SKILL.md").write_text("not part of the distribution")
    result = run_install(harness, "--ide", "generic")
    assert result.returncode == 0, result.stderr
    assert not (project / "samovar-harness/skills/private").exists()
    installed = project / "samovar-harness/skills/samovar/.samovar-harness.json"
    assert json.loads(installed.read_text())["format_version"] == 1


def test_distributed_skills_have_validation_records():
    # Release gate from SF-81/#76: a new shipped skill needs recorded scenarios
    # and execution evidence. This checks packaging, not scientific correctness.
    manifest = json.loads((ROOT / "Agents/bundle.json").read_text())
    for name in manifest["skills"]:
        records = ROOT / "Agents/skill-and-hook-testing" / name
        assert (records / "cases.md").is_file(), f"Missing scenarios for {name}"
        assert (records / "PASS.md").is_file(), f"Missing execution evidence for {name}"
