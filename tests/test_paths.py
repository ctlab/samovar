from pathlib import Path

from samovar.paths import (
    PACKAGE_VERSION,
    ncbi_email,
    repo_root,
    resolve_executable,
    user_config_dir,
    workflow_dir,
)


def test_package_version():
    assert PACKAGE_VERSION == "0.10.1"


def test_repo_root_contains_workflow():
    root = repo_root()
    assert (root / "workflow" / "annotators" / "Snakefile").is_file()
    assert workflow_dir() == root / "workflow"


def test_user_config_dir_respects_xdg(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    assert user_config_dir() == tmp_path / "xdg" / "samovar"


def test_ncbi_email_from_env(monkeypatch):
    monkeypatch.setenv("NCBI_EMAIL", "test@samovar.com")
    assert ncbi_email() == "test@samovar.com"


def test_resolve_executable_absolute(tmp_path):
    exe = tmp_path / "kraken2"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    assert resolve_executable(str(exe), tool_key="kraken2") == str(exe.resolve())


def test_resolve_executable_from_path(tmp_path, monkeypatch):
    exe = tmp_path / "kaiju"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    resolved = resolve_executable("kaiju --verbose", tool_key="kaiju")
    assert resolved.startswith(str(exe.resolve()))
    assert resolved.endswith("--verbose")


def test_cli_help_from_other_directory(tmp_path):
    import subprocess

    cli = repo_root() / "bin" / "samovar"
    proc = subprocess.run(
        [str(cli), "help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert proc.returncode == 0, proc.stderr
    assert "prepare" in proc.stdout.lower() or "Usage:" in proc.stdout
