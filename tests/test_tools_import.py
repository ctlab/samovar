import json
import stat
from pathlib import Path

from samovar.annotators_wrapper import CustomAnnotator, get_annotator_instance
from samovar.main_config import parse_tool_entry
from samovar.paths import write_config, update_config
from samovar.tools_import import import_tool, main as import_main


def test_import_tool_writes_config(tmp_path, monkeypatch):
    binary = tmp_path / "mykaiju"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="kaiju",
        tool_type="a",
        env="",
        exec_name="kaiju",
        exec_path=str(binary),
        also_repo_build=False,
    )
    assert spec[2] == str(binary.resolve())
    assert spec[3] == "annotator"
    assert spec[1] == "bash"
    loaded = parse_tool_entry(json.loads(cfg.read_text())["tools"]["kaiju"], "kaiju")
    assert loaded[2] == spec[2]
    assert loaded[3] == "annotator"


def test_import_cli_conda_prefix(tmp_path, monkeypatch):
    prefix = tmp_path / "env"
    exe = prefix / "bin" / "simulator.py"
    exe.parent.mkdir(parents=True)
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    rc = import_main(
        [
            "-n",
            "nanosim",
            "--env",
            "conda",
            "--exec",
            "simulator.py",
            "--exec-path",
            str(prefix),
            "--type",
            "meta",
        ]
    )
    assert rc == 0
    row = json.loads(cfg.read_text())["tools"]["nanosim"]
    spec = parse_tool_entry(row, "nanosim")
    assert spec[0] == "conda"
    assert spec[3] == "metagenome_generator"
    assert Path(spec[2]).resolve() == prefix.resolve()


def test_import_flags_fifth_slot(tmp_path, monkeypatch):
    binary = tmp_path / "myboil.py"
    binary.write_text("def regenerate(annotation, metadata, config):\n    return {}\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="myboil",
        tool_type="table",
        exec_path=str(binary),
        flags="--log-mu 1 --foo bar",
        also_repo_build=False,
    )
    assert spec[3] == "table_reads_generator"
    assert spec[4] == "--log-mu 1 --foo bar"
    raw = json.loads(cfg.read_text())["tools"]["myboil"]
    assert len(raw) == 5
    assert raw[4] == "--log-mu 1 --foo bar"
    spec4 = import_tool(
        name="plain",
        tool_type="table",
        exec_path=str(binary),
        also_repo_build=False,
    )
    assert len(spec4) == 4
    raw4 = json.loads(cfg.read_text())["tools"]["plain"]
    assert len(raw4) == 4


def test_imported_annotator_invokes_binary_not_custom_sh(tmp_path):
    script = tmp_path / "clf"
    script.write_text("#!/bin/sh\nexit 0\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    inst = get_annotator_instance(
        "clf",
        {
            "run_name": "clf",
            "type": "clf",
            "cmd": str(script),
            "db_path": "/db",
            "threads": 4,
        },
        {},
    )
    assert isinstance(inst, CustomAnnotator)
    cmd = inst.get_snakemake_shell_cmd("a_R1.fastq", "a_R2.fastq", ["/tmp/out.out"])
    assert "custom.sh" not in cmd
    assert str(script) in cmd
    assert "-i a_R1.fastq" in cmd
    assert "-d /db" in cmd
    assert "-p clf" not in cmd
