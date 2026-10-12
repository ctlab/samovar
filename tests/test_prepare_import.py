"""samovar prepare --autoimport."""

import argparse
import subprocess
import sys

from samovar.prepare_import import autoimport_mode, offer_prepare_imports


def _offer(args, **kwargs):
    recorded = []
    asked = []
    databases = []

    def import_tool(**kw):
        recorded.append(kw)

    def import_db(**kw):
        databases.append(kw)

    def ask_tool(spec, **_kw):
        asked.append(spec.name)
        return "skipped"

    def no_conda(*_args):
        raise FileNotFoundError("no conda env")

    defaults = dict(
        ask=lambda _prompt: "no",
        interactive=True,
        import_tool_fn=import_tool,
        import_db_fn=import_db,
        ask_tool_fn=ask_tool,
        configured_fn=lambda *_a, **_k: None,
        db_known_fn=lambda *_a: False,
        which_fn=lambda _names: [],
        resolve_conda_fn=no_conda,
        sibling_fn=lambda *_a: "",
    )
    defaults.update(kwargs)
    offer_prepare_imports(args, **defaults)
    return recorded, asked, databases


def test_autoimport_mode():
    assert autoimport_mode(None) == "ask"
    assert autoimport_mode(1) == "auto"
    assert autoimport_mode(0) == "off"


def test_autoimport_0_imports_nothing():
    args = argparse.Namespace(autoimport=0, cmd_kraken2_test=[["kraken2 toy"]])
    recorded, asked, databases = _offer(args, which_fn=lambda _names: ["/usr/bin/kraken2"])
    assert recorded == []
    assert asked == []
    assert databases == []


def test_autoimport_1_records_single_path_hit():
    args = argparse.Namespace(autoimport=1, cmd_kraken2_test=[["kraken2 toy"]])
    recorded, asked, _databases = _offer(
        args, which_fn=lambda _names: ["/usr/bin/kraken2"]
    )
    assert asked == []
    assert recorded[0]["name"] == "kraken2"
    assert recorded[0]["exec_path"] == "/usr/bin/kraken2"


def test_autoimport_1_records_conda_env_when_path_is_empty():
    args = argparse.Namespace(autoimport=1, cmd_kraken2_test=[["kraken2 toy"]])

    def conda(_env, _exe):
        return "/envs/kraken2/bin/kraken2", "/envs/kraken2"

    recorded, asked, _databases = _offer(args, resolve_conda_fn=conda)
    assert asked == []
    assert recorded[0]["env"] == "conda"
    assert recorded[0]["exec_path"] == "/envs/kraken2/bin/kraken2"


def test_autoimport_1_asks_when_nothing_resolves():
    args = argparse.Namespace(autoimport=1, cmd_kraken2_test=[["kraken2 toy"]])
    prompts = []

    def ask(prompt):
        prompts.append(prompt)
        return "no"

    _recorded, asked, databases = _offer(args, ask=ask)
    assert asked == ["kraken2"]
    assert databases == []
    assert any("database toy" in prompt for prompt in prompts)


def test_autoimport_1_records_database_path_without_asking(tmp_path):
    db = tmp_path / "toy"
    db.mkdir()
    binary = tmp_path / "kraken2"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)
    args = argparse.Namespace(autoimport=1, cmd_kraken2_test=[[f"{binary} {db}"]])
    recorded, asked, databases = _offer(args)
    assert asked == []
    assert recorded[0]["exec_path"] == str(binary)
    assert databases[0]["name"] == "toy"
    assert databases[0]["tool"] == "kraken2"


def test_omitted_flag_asks_even_when_binary_is_on_path():
    args = argparse.Namespace(cmd_kraken2_test=[["kraken2 toy"]])
    recorded, asked, _databases = _offer(
        args, which_fn=lambda _names: ["/usr/bin/kraken2"]
    )
    assert recorded == []
    assert asked == ["kraken2"]


def test_builtin_table_generator_is_not_imported():
    args = argparse.Namespace(autoimport=1, table_reads_generator=["direct"])
    recorded, asked, databases = _offer(
        args, which_fn=lambda _names: ["/usr/bin/direct"]
    )
    assert recorded == []
    assert asked == []
    assert databases == []


def test_autoimport_1_without_a_terminal_does_not_ask():
    args = argparse.Namespace(autoimport=1, cmd_kraken2_test=[["kraken2 toy"]])
    recorded, asked, databases = _offer(args, interactive=False)
    assert recorded == []
    assert asked == []
    assert databases == []


def test_prepare_help_documents_autoimport():
    proc = subprocess.run(
        [sys.executable, "bin/samovar_preprocess", "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert "--autoimport" in proc.stdout
    assert "import nothing" in proc.stdout
    rejected = subprocess.run(
        [sys.executable, "bin/samovar_preprocess", "--autoimport", "2", "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert rejected.returncode != 0
