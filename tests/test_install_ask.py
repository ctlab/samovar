import json
from pathlib import Path

from samovar.install_ask import (
    ASK_TOOLS,
    AskTool,
    ask_tool,
    parse_env_answer,
    pick_which,
    record_tool,
)
from samovar.main_config import lookup_tool_record, tool_path
from samovar.paths import load_config


KAIJU = next(tool for tool in ASK_TOOLS if tool.name == "kaiju")


def test_parse_env_answer_formats():
    assert parse_env_answer("n") == "no"
    assert parse_env_answer("N") == "no"
    assert parse_env_answer("no") == "no"
    assert parse_env_answer("base kaiju") == ("base", "", "kaiju")
    assert parse_env_answer("conda kaiju kaiju") == ("conda", "kaiju", "kaiju")
    assert parse_env_answer("conda") is None
    assert parse_env_answer("maybe") is None
    assert parse_env_answer("") is None


def test_pick_which_number():
    hits = ["/opt/a/kaiju", "/opt/b/kaiju"]
    assert pick_which("", hits) == hits[0]
    assert pick_which("y", hits) == hits[0]
    assert pick_which("2", hits) == hits[1]
    assert pick_which("n", hits) is None
    assert pick_which("9", hits) == ""


def test_which_yes_records_that_binary():
    recorded = []

    def ask(prompt):
        assert "1) /usr/bin/kaiju" in prompt
        return "y"

    result = ask_tool(
        KAIJU,
        ask=ask,
        which_fn=lambda names: ["/usr/bin/kaiju"],
        record=lambda spec, path, **kwargs: recorded.append((spec.name, path, kwargs)),
    )
    assert result == "imported-which"
    assert recorded == [("kaiju", "/usr/bin/kaiju", {"env": "", "prefix": "", "exec_name": "kaiju"})]


def test_decline_which_then_conda_env(tmp_path):
    answers = iter(["n", "conda kaiju kaiju"])
    recorded = []

    result = ask_tool(
        KAIJU,
        ask=lambda prompt: next(answers),
        which_fn=lambda names: ["/usr/bin/kaiju"],
        resolve_conda_fn=lambda env, exe: (str(tmp_path / "bin" / exe), str(tmp_path)),
        record=lambda spec, path, **kwargs: recorded.append((path, kwargs)),
    )
    assert result == "imported-env"
    path, kwargs = recorded[0]
    assert path.endswith("/bin/kaiju")
    assert kwargs["env"] == "conda"
    assert kwargs["prefix"] == str(tmp_path)
    assert kwargs["exec_name"] == "kaiju"


def test_no_then_skip_install():
    answers = iter(["n", "n"])
    called = []
    result = ask_tool(
        KAIJU,
        ask=lambda prompt: next(answers),
        which_fn=lambda names: [],
        install=lambda spec: called.append(spec.name),
        record=lambda *args, **kwargs: called.append("record"),
    )
    assert result == "skipped"
    assert called == []


def test_no_then_install_into_samovar():
    answers = iter(["no", "y"])
    seen = []

    def which_fn(names):
        return ["/opt/samovar/bin/kaiju"] if seen else []

    def install(spec):
        seen.append(spec.name)

    recorded = []
    result = ask_tool(
        KAIJU,
        ask=lambda prompt: next(answers),
        which_fn=which_fn,
        install=install,
        locate=lambda spec: ("/opt/samovar/bin/kaiju", "/opt/samovar"),
        record=lambda spec, path, **kwargs: recorded.append((path, kwargs["env"], kwargs["prefix"])),
    )
    assert result == "installed"
    assert seen == ["kaiju"]
    assert recorded == [("/opt/samovar/bin/kaiju", "conda", "/opt/samovar")]


def test_base_exec_updates_config():
    recorded = []
    result = ask_tool(
        KAIJU,
        ask=lambda prompt: "base kaiju",
        which_fn=lambda names: [],
        resolve_base_fn=lambda exe: ("/envs/samovar/bin/kaiju", "/envs/samovar"),
        record=lambda spec, path, **kwargs: recorded.append(kwargs),
    )
    assert result == "imported-env"
    assert recorded[0]["prefix"] == "/envs/samovar"
    assert recorded[0]["exec_name"] == "kaiju"


def test_invalid_then_conda():
    answers = iter(["later", "conda kaiju kaiju"])
    result = ask_tool(
        KAIJU,
        ask=lambda prompt: next(answers),
        which_fn=lambda names: [],
        resolve_conda_fn=lambda env, exe: ("/opt/envs/kaiju/bin/kaiju", "/opt/envs/kaiju"),
        record=lambda *args, **kwargs: None,
    )
    assert result == "imported-env"


def test_record_tool_writes_main_config(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg_path))
    monkeypatch.setenv("SAMOVAR_ROOT", str(tmp_path))
    exe = tmp_path / "envs" / "kaiju" / "bin" / "kaiju"
    exe.parent.mkdir(parents=True)
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    record_tool(
        KAIJU,
        str(exe),
        env="conda",
        prefix=str(tmp_path / "envs" / "kaiju"),
        exec_name="kaiju",
    )
    cfg = load_config()
    rec = lookup_tool_record(cfg, "kaiju")
    assert rec is not None
    assert tool_path(rec, "kaiju") == str(exe.resolve())
    envs = cfg.get("tool_envs") or {}
    assert str(tmp_path / "envs" / "kaiju") in {str(value) for value in envs.values()}
    assert json.loads(cfg_path.read_text(encoding="utf-8"))["tools"]


def test_ask_table_is_unique_and_typed():
    names = [tool.name for tool in ASK_TOOLS]
    assert names.count("kaiju") == 1
    assert len(names) == len(set(names))
    for tool in ASK_TOOLS:
        assert isinstance(tool, AskTool)
        assert tool.install.startswith(("conda:", "optional:"))
        assert tool.which
