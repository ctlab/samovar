from pathlib import Path

from samovar.paths import (
    PACKAGE_VERSION,
    annotation_regenerate_r,
    collect_runtime_path_dirs,
    discover_multiqc,
    discover_opal,
    iss_executable,
    ncbi_email,
    repo_root,
    resolve_executable,
    user_config_dir,
    workflow_dir,
)


def test_package_version():
    assert PACKAGE_VERSION == "0.10.13"


def test_discover_opal_from_config(tmp_path, monkeypatch):
    script = tmp_path / "opal.py"
    script.write_text("#!/usr/bin/env python3\nprint('opal')\n")
    monkeypatch.setattr(
        "samovar.paths.load_config",
        lambda: {"opal_path": str(script)},
    )
    monkeypatch.delenv("SAMOVAR_OPAL_PATH", raising=False)
    monkeypatch.delenv("SAMOVAR_OPAL_BIN", raising=False)
    assert discover_opal() == str(script.resolve())


def test_discover_multiqc_from_config(tmp_path, monkeypatch):
    exe = tmp_path / "multiqc"
    exe.write_text("#!/usr/bin/env python3\nprint('multiqc')\n")
    exe.chmod(0o755)
    monkeypatch.setattr(
        "samovar.paths.load_config",
        lambda: {"multiqc_path": str(exe)},
    )
    monkeypatch.delenv("SAMOVAR_MULTIQC_PATH", raising=False)
    monkeypatch.delenv("SAMOVAR_MULTIQC_BIN", raising=False)
    assert discover_multiqc() == str(exe.resolve())


def test_collect_runtime_path_dirs_includes_opal(tmp_path):
    opal = tmp_path / "bin" / "opal.py"
    opal.parent.mkdir()
    opal.write_text("#!/usr/bin/env python3\n")
    dirs = collect_runtime_path_dirs({"opal_path": str(opal)})
    assert str(opal.parent.resolve()) in dirs


def test_portable_annotator_cmd_warns_when_missing(capsys, monkeypatch):
    from samovar.config import _portable_annotator_cmd

    monkeypatch.setattr("samovar.paths.load_config", lambda: {})
    monkeypatch.setattr("samovar.config.shutil.which", lambda name: None)
    monkeypatch.setattr("samovar.paths.shutil.which", lambda name: None)
    out = _portable_annotator_cmd("definitely_not_a_real_annotator_xyz /db")
    assert out == "definitely_not_a_real_annotator_xyz /db"
    err = capsys.readouterr().err
    assert "Warning:" in err
    assert "definitely_not_a_real_annotator_xyz" in err
    assert "not on PATH" in err


def test_portable_annotator_cmd_silent_when_absolute_exists(tmp_path, capsys):
    from samovar.config import _portable_annotator_cmd

    exe = tmp_path / "kaiju"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    out = _portable_annotator_cmd(f"{exe} /db")
    assert str(exe) in out
    assert "Warning:" not in capsys.readouterr().err


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


def test_collect_runtime_path_dirs_tools_path_and_envs(tmp_path):
    py = tmp_path / "conda" / "bin" / "python"
    py.parent.mkdir(parents=True)
    py.write_text("#!/bin/sh\n")
    py.chmod(0o755)
    kaiju_env = tmp_path / "kaiju_env"
    kaiju_bin = kaiju_env / "bin"
    kaiju_bin.mkdir(parents=True)
    kaiju = kaiju_bin / "kaiju"
    kaiju.write_text("#!/bin/sh\n")
    kaiju.chmod(0o755)
    extra_env = tmp_path / "other_env"
    (extra_env / "bin").mkdir(parents=True)
    k2_prefix = tmp_path / "kraken2_env"
    dirs = collect_runtime_path_dirs(
        {
            "python_path": str(py),
            "path": [str(extra_env)],
            "tools": {"kaiju": str(kaiju)},
            "tool_envs": {"kraken2": str(k2_prefix)},
        }
    )
    assert str(py.parent) in dirs
    assert str(kaiju_bin) in dirs
    assert str(extra_env / "bin") in dirs
    assert str(k2_prefix / "bin") in dirs


def test_collect_runtime_path_dirs_missing_bin_not_doubled():
    dirs = collect_runtime_path_dirs(
        {
            "path": ["/opt/other-env/bin"],
            "tools": {"kraken2": "/opt/kraken/bin/kraken2"},
            "tool_envs": {"kaiju": "/opt/conda/envs/kaiju"},
        }
    )
    assert "/opt/other-env/bin" in dirs
    assert "/opt/other-env/bin/bin" not in dirs
    assert "/opt/conda/envs/kaiju/bin" in dirs
    assert "/opt/kraken/bin" in dirs


def test_resolve_executable_from_path(tmp_path, monkeypatch):
    exe = tmp_path / "kaiju"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setattr("samovar.paths.load_config", lambda: {})
    resolved = resolve_executable("kaiju --verbose", tool_key="kaiju")
    assert resolved.startswith(str(exe.resolve()))
    assert resolved.endswith("--verbose")


def test_resolve_executable_uses_tool_envs(tmp_path, monkeypatch):
    env = tmp_path / "kaiju_env"
    (env / "bin").mkdir(parents=True)
    exe = env / "bin" / "kaiju"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setattr(
        "samovar.paths.load_config",
        lambda: {"tool_envs": {"kaiju": str(env)}},
    )
    assert resolve_executable("kaiju", tool_key="kaiju") == str(exe.resolve())


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
