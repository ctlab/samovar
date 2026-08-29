import json
from pathlib import Path
from types import SimpleNamespace

from omegaconf import OmegaConf

from samovar.paths import write_config
from samovar.repro import (
    args_to_overrides,
    compare_tools,
    compose_snapshot,
    empty_snapshot,
    export_run,
    load_lazy_install_text,
    record_stage,
    reproduce,
    rewrite_output_dir_argv,
)
from samovar.tools_import import import_tool


def test_load_lazy_install_file_and_command(tmp_path):
    script = tmp_path / "install.sh"
    script.write_text("conda install -y bioconda::kraken2=2.1.3\n")
    assert "kraken2=2.1.3" in load_lazy_install_text(f"@{script}")
    assert load_lazy_install_text("conda install bioconda::kaiju").startswith("conda")
    assert "kraken2" in load_lazy_install_text(str(script))


def test_import_stores_multiline_lazy(tmp_path, monkeypatch):
    binary = tmp_path / "myclf"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    recipe = tmp_path / "build.sh"
    recipe.write_text("#!/bin/bash\nconda create -y -p \"$PREFIX\"\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="myclf",
        tool_type="annotator",
        exec_path=str(binary),
        lazy_install=str(recipe),
        version="0.0.1",
        also_repo_build=False,
    )
    tools = json.loads(cfg.read_text())["tools"]
    rec = tools["myclf:0.0.1"]
    assert "conda create" in rec["lazy-install"]
    assert "PREFIX" in rec["lazy-install"]


def test_record_and_export_slim(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": ["--threads"]}, also_repo_build=False)
    out = tmp_path / "run"
    args = SimpleNamespace(output_dir=str(out), n_samples=2, threads=4)
    dest = record_stage(
        "generate", str(out), args=args, argv=["--n_samples", "2", "--threads", "4"]
    )
    assert dest is not None
    assert (dest / "config.yaml").is_file()
    snap = OmegaConf.load(dest / "config.yaml")
    assert "generate" in list(snap.stages)
    assert snap.generate.args.n_samples == 2
    archive = tmp_path / "run.tgz"
    export_run(out, archive, mode="slim")
    assert archive.is_file()
    code = reproduce(archive, str(tmp_path / "rerun"), install=False, dry_run=True)
    assert code == 0
    snap2 = compose_snapshot(dest, ["generate.n_samples=9"])
    assert int(snap2.generate.n_samples) == 9 or int(snap2.generate.args.n_samples) == 9
    from samovar.repro import patch_argv_overrides

    patched = patch_argv_overrides(
        ["--n_samples", "2", "--threads", "4"], "generate", ["generate.n_samples=9"]
    )
    assert patched[patched.index("--n_samples") + 1] == "9"
    assert (dest / "hydra.yaml").is_file()
    assert (dest / "overrides.yaml").is_file()
    code = reproduce(
        dest,
        str(tmp_path / "rerun2"),
        install=False,
        dry_run=True,
        extra_overrides=["generate.n_samples=9"],
    )
    assert code == 0


def test_export_full_packs_db_and_lazy_install(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    binary = tmp_path / "myclf"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    db = tmp_path / "db"
    db.mkdir()
    (db / "hash.k2d").write_text("x")
    write_config(
        {
            "root": str(tmp_path),
            "tools": {
                "myclf:0.0.1": {
                    "type": "annotator",
                    "exec": {"path": str(binary), "env": "bash"},
                    "lazy-install": "#!/bin/bash\necho rebuild\n",
                    "flags-translate": {"--threads": "--threads"},
                }
            },
            "databases": {"myclf": [["testdb", str(db), "--threads 1"]]},
            "custom-flags": ["--threads"],
        },
        also_repo_build=False,
    )
    out = tmp_path / "run"
    record_stage(
        "generate",
        str(out),
        args=SimpleNamespace(output_dir=str(out), n_samples=1),
        argv=["--n_samples", "1"],
    )
    archive = tmp_path / "full.tgz"
    export_run(out, archive, mode="full")
    import tarfile

    names = tarfile.open(archive, "r:gz").getnames()
    assert any(n.endswith("payload/tools/myclf/lazy-install.sh") for n in names)
    assert any(n.endswith("payload/tools/myclf/meta.json") for n in names)
    assert any("payload/databases/myclf/testdb" in n for n in names)
    assert any(n.endswith("install.sh") for n in names)
    extract = tmp_path / "unfold"
    extract.mkdir()
    with tarfile.open(archive, "r:gz") as tar:
        kwargs = {"filter": "data"} if hasattr(tarfile, "data_filter") else {}
        tar.extractall(extract, **kwargs)
    recipe = extract / "samovar-run" / "payload" / "tools" / "myclf" / "lazy-install.sh"
    assert "rebuild" in recipe.read_text()
    script = (extract / "samovar-run" / "install.sh").read_text()
    assert "--type" in script
    assert "lazy-install-file" in script


def test_rewrite_output_dir_and_overrides():
    argv = rewrite_output_dir_argv(["--genome_dir", "g", "--output_dir", "old"], "/new")
    assert argv[-1] == "/new"
    ov = args_to_overrides("generate", {"n_samples": 1, "threads": 4, "extra": None})
    assert "generate.threads=4" in ov


def test_stages_and_overrides_split():
    from samovar.repro import _stages_and_overrides

    wanted, extra = _stages_and_overrides("generate", ["generate.n_samples=3"])
    assert wanted == ["generate"]
    assert extra == ["generate.n_samples=3"]
    wanted, extra = _stages_and_overrides("generate,prepare", [])
    assert wanted == ["generate", "prepare"]
    wanted, extra = _stages_and_overrides("generate.n_samples=3", [])
    assert wanted is None
    assert extra == ["generate.n_samples=3"]


def test_compare_tools_ok_when_versions_match():
    snap = empty_snapshot()
    snap.tools = {
        "python": {
            "name": "python",
            "version": "",
            "builtin": True,
            "exec": {"path": ""},
            "lazy_install": "",
        }
    }
    rows = compare_tools(snap)
    assert rows
    assert rows[0]["name"] == "python"
