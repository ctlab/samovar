from pathlib import Path

from samovar.paths import (
    PACKAGE_VERSION,
    annotation_regenerate_r,
    iss_executable,
    ncbi_email,
    repo_root,
    resolve_executable,
    user_config_dir,
    workflow_dir,
)


def test_package_version():
    assert PACKAGE_VERSION == "0.10.7"


def test_absolute_path_resolves_relative(tmp_path, monkeypatch):
    from samovar.paths import absolute_path

    monkeypatch.chdir(tmp_path)
    (tmp_path / "reads").mkdir()
    assert absolute_path("reads") == str((tmp_path / "reads").resolve())
    assert absolute_path("/already/abs") == str(Path("/already/abs").resolve())


def test_repo_root_contains_workflow():
    root = repo_root()
    assert (root / "workflow" / "annotators" / "Snakefile").is_file()
    assert workflow_dir() == root / "workflow"
    assert not (root / "DESCRIPTION").exists()
    assert not (root / "R").exists()
    assert not (root / "man").exists()
    assert not (root / "workflow" / "annotation_regenerate.R").exists()


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


def test_iss_executable_uses_config(tmp_path, monkeypatch):
    exe = tmp_path / "iss"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setattr(
        "samovar.paths.load_config",
        lambda: {"iss_path": str(exe)},
    )
    assert iss_executable() == str(exe.resolve())


def test_annotation_regenerate_r_from_env(monkeypatch, tmp_path):
    script = tmp_path / "custom_annotation_regenerate.R"
    script.write_text("# optional R regenerator\n")
    monkeypatch.setenv("SAMOVAR_R_REGENERATE", str(script))
    monkeypatch.delenv("SAMOVAR_ANNOTATION_REGENERATE_R", raising=False)
    assert annotation_regenerate_r() == script


def test_annotation_regenerate_r_missing_when_unbundled(monkeypatch, tmp_path):
    monkeypatch.delenv("SAMOVAR_R_REGENERATE", raising=False)
    monkeypatch.delenv("SAMOVAR_ANNOTATION_REGENERATE_R", raising=False)
    monkeypatch.setattr("samovar.paths.load_config", lambda: {})
    monkeypatch.setattr("samovar.paths.user_config_dir", lambda: tmp_path / "empty_xdg")
    assert annotation_regenerate_r() is None
