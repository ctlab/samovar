from pathlib import Path

from samovar.paths import (
    PACKAGE_VERSION,
    annotation_regenerate_r,
    collect_runtime_path_dirs,
    discover_multiqc,
    discover_opal,
    iss_executable,
    load_config,
    ncbi_email,
    recorded_config_path,
    repo_root,
    resolve_executable,
    user_config_dir,
    user_config_path,
    workflow_dir,
    write_config,
    write_install_config_pointer,
)
from samovar.version import get_version, pyproject_path


def test_package_version():
    assert PACKAGE_VERSION == get_version()
    text = pyproject_path().read_text(encoding="utf-8")
    assert f'version = "{PACKAGE_VERSION}"' in text or f"version = '{PACKAGE_VERSION}'" in text


def test_build_config_path_pointer(tmp_path, monkeypatch):
    monkeypatch.delenv("SAMOVAR_CONFIG", raising=False)
    monkeypatch.setenv("SAMOVAR_ROOT", str(tmp_path))
    cfg = tmp_path / "custom" / "config.json"
    cfg.parent.mkdir(parents=True)
    cfg.write_text(
        '{"version": "%s", "root": "%s", "compilers": {"python": "/bin/python3"}, '
        '"API": {}, "genomes": {"test": [], "raw": {}, "processed": {}, "data": {}}, '
        '"databases": {}, "workflows": {}, "tools": {}}\n' % (get_version(), tmp_path)
    )
    pointer = write_install_config_pointer(cfg, root=tmp_path)
    assert pointer is not None
    assert pointer.read_text().strip() == str(cfg.resolve())
    assert recorded_config_path() == cfg.resolve()
    assert user_config_path() == cfg.resolve()
    assert user_config_dir() == cfg.parent
    loaded = load_config()
    assert loaded.get("python_path") == "/bin/python3"


def test_samovar_config_env_overrides_pointer(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_ROOT", str(tmp_path))
    pointed = tmp_path / "pointed.json"
    pointed.write_text('{"compilers": {"python": "/from/pointer"}}\n')
    write_install_config_pointer(pointed, root=tmp_path)
    override = tmp_path / "override.json"
    override.write_text('{"compilers": {"python": "/from/env"}}\n')
    monkeypatch.setenv("SAMOVAR_CONFIG", str(override))
    assert user_config_path() == override
    assert load_config().get("python_path") == "/from/env"


def test_write_config_records_pointer(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_ROOT", str(tmp_path))
    dest = tmp_path / "etc" / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(dest))
    write_config(
        {
            "version": get_version(),
            "root": str(tmp_path),
            "compilers": {"python": "/usr/bin/python3"},
        }
    )
    assert dest.is_file()
    assert (tmp_path / "build" / "config_path").read_text().strip() == str(dest.resolve())
    assert (tmp_path / "build" / "config.json").is_file()


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
    monkeypatch.delenv("SAMOVAR_CONFIG", raising=False)
    monkeypatch.setenv("SAMOVAR_ROOT", str(tmp_path))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    assert user_config_dir() == tmp_path / "xdg" / "samovar"
    assert user_config_path() == tmp_path / "xdg" / "samovar" / "config.json"


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


def test_discover_nanosim_from_tool_envs(tmp_path, monkeypatch):
    from samovar.paths import discover_nanosim

    env = tmp_path / "nanosim"
    bindir = env / "bin"
    bindir.mkdir(parents=True)
    exe = bindir / "simulator.py"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setattr(
        "samovar.paths.load_config",
        lambda: {"tool_envs": {"nanosim": str(env)}},
    )
    monkeypatch.delenv("SAMOVAR_NANOSIM", raising=False)
    monkeypatch.delenv("SAMOVAR_NANOSIM_BIN", raising=False)
    assert discover_nanosim() == str(exe.resolve())


def test_format_install_status_lists_required():
    from samovar.paths import format_install_status, install_status_rows

    text = format_install_status()
    assert "SamovaR tool status" in text
    assert "Required:" in text
    assert "Optional:" in text
    assert "NanoSim" in text
    names = {row["name"] for row in install_status_rows()}
    assert "python" in names
    assert "iss" in names
    assert "CAMISIM" in names
    assert "ART" in names
    assert "seqtk" in names


def test_format_install_status_without_nextflow(monkeypatch):
    """CI images often have no nextflow and an empty nextflow_path in config."""
    import shutil as _shutil

    from samovar.paths import format_install_status, install_status_rows

    monkeypatch.setattr(
        "samovar.paths.load_config",
        lambda: {"nextflow_path": "", "tools": {"nextflow": ""}},
    )
    real_which = _shutil.which

    def fake_which(cmd, path=None):
        if cmd == "nextflow":
            return None
        return real_which(cmd, path=path)

    monkeypatch.setattr("samovar.paths.shutil.which", fake_which)
    text = format_install_status()
    assert "SamovaR tool status" in text
    nxt = next(row for row in install_status_rows() if row["name"] == "Nextflow")
    assert nxt["found"] is False
    assert nxt["path"] == ""


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


def test_output_dir_cli_aliases_and_shell_override(tmp_path):
    import argparse
    import subprocess

    from samovar.paths import add_output_dir_argument, shell_outdir_override_snippet

    parser = argparse.ArgumentParser()
    add_output_dir_argument(parser, required=True)
    assert parser.parse_args(["--directory", "/run/a"]).output_dir == "/run/a"
    assert parser.parse_args(["--outdir", "/run/b"]).output_dir == "/run/b"
    assert parser.parse_args(["--output-dir", "/run/c"]).output_dir == "/run/c"

    script = tmp_path / "samovar.sh"
    script.write_text(
        "#!/bin/bash\nset -e\nout_dir=/baked\n"
        + shell_outdir_override_snippet()
        + 'printf "%s\\n" "$out_dir"\n',
        encoding="utf-8",
    )
    script.chmod(0o755)
    dest = tmp_path / "relocated"
    dest.mkdir()
    result = subprocess.run(
        ["bash", str(script), "--directory", str(dest)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.strip() == str(dest)
