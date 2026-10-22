"""OOP annotators: CustomAnnotator routes, native parsers, taxonomy_engine, features."""

import os
import json
import stat
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

from samovar.annotators_wrapper import (
    CUSTOM_WRAPPER_SCRIPTS,
    CustomAnnotator,
    KaijuAnnotator,
    Kraken2Annotator,
    get_annotator_instance,
)
from samovar.config import setup_pipeline
from samovar.parse_annotators import match_annotation
from samovar.reprofiling import load_read_features, merge_read_features
from samovar.taxonomy_engine import NCBITaxonomyParser

REPO = Path(__file__).resolve().parents[1]
ANNOTATORS = REPO / "src" / "annotators"


def _tiny_fastq(path: Path, n: int = 3, prefix: str = "read") -> None:
    records = []
    for i in range(n):
        records.append(f"@{prefix}{i}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")
    path.write_text("".join(records))


def _custom(tool: str, run_name=None, cmd=None, db_path="/db") -> CustomAnnotator:
    inst = get_annotator_instance(
        tool,
        {
            "run_name": run_name or tool,
            "type": tool,
            "cmd": cmd or tool,
            "db_path": db_path,
            "threads": 2,
        },
        {},
    )
    assert isinstance(inst, CustomAnnotator)
    return inst


@pytest.mark.parametrize(
    "tool",
    ["centrifuge", "metauto", "assembly_hybrid"],
)
def test_factory_routes_custom_tools(tool):
    inst = _custom(tool)
    assert inst.tool_name == tool
    outs = inst.get_expected_outputs("s1", "/tmp/out")
    assert outs[0].endswith(f"s1_{tool}.custom_{tool}.out")
    cmd = inst.get_snakemake_shell_cmd("a_R1.fastq", "a_R2.fastq", outs)
    assert "custom.sh" in cmd
    assert f"-p {tool}" in cmd
    assert "-d /db" in cmd
    # Binary name must not be used as the executable (wrong CLI flags).
    assert not cmd.strip().startswith(tool + " -i")


def test_custom_single_end_omits_mate_flag():
    inst = _custom("centrifuge")
    cmd = inst.get_snakemake_shell_cmd("a_R1.fastq", "", inst.get_expected_outputs("s1", "/tmp/out"))
    assert "-I " not in cmd
    assert "-i a_R1.fastq" in cmd


def test_wrapper_scripts_forward_to_custom_sh():
    for name in CUSTOM_WRAPPER_SCRIPTS:
        if name == "custom.sh":
            continue
        path = ANNOTATORS / name
        assert path.is_file()
        text = path.read_text()
        assert "custom.sh" in text
        assert "-p " in text


def test_custom_annotator_uses_wrapper_path_when_given():
    wrapper = str(ANNOTATORS / "centrifuge.sh")
    inst = _custom("centrifuge", cmd=f"bash {wrapper}")
    cmd = inst.get_snakemake_shell_cmd("r1", "r2", ["/tmp/out.out"])
    assert "centrifuge.sh" in cmd
    assert "-p centrifuge" in cmd


def test_custom_annotator_extra_on_router():
    router = _custom("centrifuge")
    router.extra = "--foo 1"
    cmd = router.get_snakemake_shell_cmd("r1", "r2", ["/tmp/x.out"])
    assert cmd.rstrip().endswith("--foo 1")


def test_custom_sh_accepts_threads_long_option(tmp_path):
    """flags-translate emits --threads; getopts treated that as illegal --."""
    r1 = tmp_path / "a_1.fastq"
    r2 = tmp_path / "a_2.fastq"
    _tiny_fastq(r1)
    _tiny_fastq(r2)
    out = tmp_path / "dummy.out"
    script = ANNOTATORS / "custom.sh"
    proc = subprocess.run(
        [
            "bash",
            str(script),
            "-i",
            str(r1),
            "-I",
            str(r2),
            "-d",
            "/db",
            "-o",
            str(out),
            "-p",
            "dummy",
            "--threads",
            "2",
            "--threads",
            "2",
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    assert "illegal option" not in (proc.stderr + proc.stdout).lower()
    assert out.is_file()


def test_match_annotation_custom_tools():
    assert match_annotation("s.custom_centrifuge.out") == "centrifuge"
    assert match_annotation("s.custom_metauto.out") == "metauto"
    assert match_annotation("s.custom_assembly_hybrid.out") == "assembly_hybrid"
    assert match_annotation("s.kaiju.out") == "kaiju"


def test_native_kraken2_and_kaiju_parse(tmp_path):
    kraken = tmp_path / "k.kraken2.out"
    kraken.write_text(
        "C\tread0\tEscherichia coli (taxid 562)\t100\t562:10\n"
        "U\tread1\tunclassified (taxid 0)\t100\t0:10\n"
    )
    parsed = Kraken2Annotator({"run_name": "k2"}, {}).parse_output(str(kraken))
    assert list(parsed["taxID"].astype(str)) == ["562", "0"]

    kaiju = tmp_path / "k.kaiju.out"
    kaiju.write_text("C\tread0\t562\nU\tread1\t0\n")
    parsed = KaijuAnnotator({"run_name": "kj"}, {}).parse_output(str(kaiju))
    assert list(parsed["taxID"].astype(str)) == ["562", "0"]


def test_custom_parse_output_two_column_and_empty(tmp_path):
    inst = _custom("centrifuge")
    filled = tmp_path / "ok.out"
    filled.write_text("read0\t562\nread1\t1280\n")
    df = inst.parse_output(str(filled))
    assert list(df.columns) == ["seq", "taxID"]
    assert list(df["taxID"].astype(str)) == ["562", "1280"]

    empty = tmp_path / "empty.out"
    empty.write_text("")
    df = inst.parse_output(str(empty))
    assert df.empty
    assert list(df.columns) == ["seq", "taxID"]


def test_custom_sh_empty_fastq_and_unknown_tool(tmp_path):
    r1 = tmp_path / "e_R1.fastq"
    r2 = tmp_path / "e_R2.fastq"
    r1.write_text("")
    r2.write_text("")
    out = tmp_path / "empty.out"
    proc = subprocess.run(
        [
            "bash",
            str(ANNOTATORS / "custom.sh"),
            "-i",
            str(r1),
            "-I",
            str(r2),
            "-d",
            ".",
            "-o",
            str(out),
            "-p",
            "centrifuge",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert out.exists()
    assert out.stat().st_size == 0

    bad = subprocess.run(
        [
            "bash",
            str(ANNOTATORS / "custom.sh"),
            "-i",
            str(r1),
            "-I",
            str(r2),
            "-d",
            ".",
            "-o",
            str(tmp_path / "no.out"),
            "-p",
            "not_a_real_tool",
        ],
        capture_output=True,
        text=True,
    )
    # empty R1 exits 0 before the router; give it one real record
    _tiny_fastq(r1, n=1)
    _tiny_fastq(r2, n=1)
    bad = subprocess.run(
        [
            "bash",
            str(ANNOTATORS / "custom.sh"),
            "-i",
            str(r1),
            "-I",
            str(r2),
            "-d",
            ".",
            "-o",
            str(tmp_path / "no.out"),
            "-p",
            "not_a_real_tool",
        ],
        capture_output=True,
        text=True,
    )
    assert bad.returncode != 0
    assert "Unknown tool" in bad.stdout + bad.stderr


def test_custom_sh_centrifuge_with_mock_binary(tmp_path):
    mock_bin = tmp_path / "bin"
    mock_bin.mkdir()
    centrifuge = mock_bin / "centrifuge"
    centrifuge.write_text(
        "#!/bin/bash\n"
        f'printf "%s\\n" "$@" > "{tmp_path / "cf.argv"}"\n'
        "OUT=''\n"
        "R1=''\n"
        "while [[ $# -gt 0 ]]; do\n"
        "  case $1 in\n"
        "    -S) OUT=$2; shift 2 ;;\n"
        "    -1) R1=$2; shift 2 ;;\n"
        "    *) shift ;;\n"
        "  esac\n"
        "done\n"
        "rid=$(head -1 \"$R1\" | sed 's/^@//;s/ .*//;s|/1$||')\n"
        "printf 'readID\\tseqID\\ttaxID\\tscore\\n' > \"$OUT\"\n"
        "printf '%s\\tref\\t562\\t100\\n' \"$rid\" >> \"$OUT\"\n"
        "printf 'other\\tref\\t0\\t1\\n' >> \"$OUT\"\n"
    )
    centrifuge.chmod(centrifuge.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    from samovar.tools_import import import_tool

    old_cfg = os.environ.get("SAMOVAR_CONFIG")
    old_root = os.environ.get("SAMOVAR_ROOT")
    os.environ["SAMOVAR_CONFIG"] = str(cfg)
    os.environ["SAMOVAR_ROOT"] = str(tmp_path)
    try:
        import_tool(
            name="centrifuge",
            tool_type="annotator",
            exec_path=str(centrifuge),
            also_repo_build=False,
        )
    finally:
        if old_cfg is None:
            os.environ.pop("SAMOVAR_CONFIG", None)
        else:
            os.environ["SAMOVAR_CONFIG"] = old_cfg
        if old_root is None:
            os.environ.pop("SAMOVAR_ROOT", None)
        else:
            os.environ["SAMOVAR_ROOT"] = old_root

    db = tmp_path / "db"
    db.mkdir()
    (db / "idx.1.cf").write_text("")
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    _tiny_fastq(r1, n=2)
    _tiny_fastq(r2, n=2)
    out = tmp_path / "cf.out"
    env = os.environ.copy()
    payload = json.loads(cfg.read_text())
    payload["compilers"] = {"python": sys.executable}
    cfg.write_text(json.dumps(payload))
    env["SAMOVAR_CONFIG"] = str(cfg)
    env["SAMOVAR_ROOT"] = str(tmp_path)
    env["PYTHON_PATH"] = sys.executable
    env["PATH"] = "/usr/bin:/bin"
    proc = subprocess.run(
        [
            "bash",
            str(ANNOTATORS / "custom.sh"),
            "-i",
            str(r1),
            "-I",
            str(r2),
            "-d",
            str(db),
            "-o",
            str(out),
            "-p",
            "centrifuge",
            "-t",
            "1", "--score-min", "17", "-x", "/wrong", "-p99",
        ],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(REPO),
    )
    if proc.returncode != 0:
        pytest.fail(proc.stdout + "\n" + proc.stderr)
    assert str(db / "idx") in proc.stdout
    df = pd.read_table(out, header=None, names=["seq", "taxID"])
    assert (df["taxID"].astype(str) == "562").all()
    assert "read0" in set(df["seq"].astype(str))


    argv = (tmp_path / 'cf.argv').read_text().splitlines()
    assert argv[-2:] == ['--score-min', '17']
    assert argv.count('-x') == 1 and '/wrong' not in argv and '-p99' not in argv

def test_prepare_custom_test_alias(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    _tiny_fastq(reads / "1_full_R1.fastq")
    _tiny_fastq(reads / "1_full_R2.fastq")
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(reads),
            "output_dir": str(tmp_path / "out"),
            "kraken2": None,
            "kaiju": None,
            "dummy": None,
            "cmd_custom-test": [["metauto /path/to/db"]],
        },
    )()
    result = setup_pipeline(args)
    init_cfg = yaml.safe_load(Path(result["configs"]["init_annotator"]).read_text())
    runs = init_cfg["run_config"]
    assert any(r["type"] == "metauto" and r["run_name"] == "metauto" for r in runs)


def test_taxonomy_engine_lineage(tmp_path):
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(
        "562\t|\t561\t|\tspecies\t|\t\n"
        "561\t|\t543\t|\tgenus\t|\t\n"
        "543\t|\t91347\t|\tfamily\t|\t\n"
        "91347\t|\t1\t|\torder\t|\t\n"
        "1\t|\t1\t|\tno rank\t|\t\n"
        "9606\t|\t9605\t|\tspecies\t|\t\n"
        "9605\t|\t1\t|\tgenus\t|\t\n"
    )
    tax = NCBITaxonomyParser(str(nodes))
    assert tax.get_ancestor_by_rank(562, "genus") == 561
    assert tax.get_ancestor_by_rank(562, "species") == 562
    assert tax.get_ancestor_by_rank(562, "family") == 543
    assert tax.get_lineage_ranks(9606)["genus"] == 9605


def test_fastq_annotator_and_feature_merge(tmp_path):
    fq = tmp_path / "reads.fastq"
    _tiny_fastq(fq, n=8, prefix="r")
    out = tmp_path / "features.tsv"
    proc = subprocess.run(
        [
            "python3",
            str(ANNOTATORS / "fastq_annotator.py"),
            str(fq),
            "-o",
            str(out),
            "--clusters",
            "2",
            "--chunk_size",
            "4",
        ],
        capture_output=True,
        text=True,
        cwd=str(REPO),
    )
    if proc.returncode != 0:
        pytest.fail(proc.stdout + "\n" + proc.stderr)
    feats = load_read_features(str(out))
    assert feats is not None
    assert "gc" in feats.columns
    assert len(feats) == 8

    ann = pd.DataFrame({"seq": ["r0", "r1"], "taxid_kaiju": [562, 0], "true": [562, 0]})
    merged = merge_read_features(ann, feats)
    assert "gc" in merged.columns
    assert merged.loc[merged["seq"] == "r0", "gc"].notna().all()


@pytest.mark.parametrize("separator", [[], ["--"]])
def test_custom_sh_metax_uses_configured_taxdump(tmp_path, separator):
    dump = tmp_path / "taxdump"
    dump.mkdir()
    (dump / "nodes.dmp").write_text("1\t|\t1\t|\n")
    (dump / "taxdump.version").write_text("2024-08-14\n")
    db = tmp_path / "metax_db"
    db.mkdir()
    (db / "metax_db.json").write_text("{}\n")
    mock_bin = tmp_path / "bin"
    mock_bin.mkdir()
    metax = mock_bin / "metax"
    log = tmp_path / "metax.argv"
    metax.write_text(
        "#!/bin/bash\n"
        f'printf "%s\\n" "$@" > "{log}"\n'
        "prefix=''\n"
        "prev=''\n"
        "for arg in \"$@\"; do\n"
        "  if [ \"$prev\" = \"--outprefix\" ]; then prefix=\"$arg\"; fi\n"
        "  prev=\"$arg\"\n"
        "done\n"
        "printf 'read1\\t562\\n' > \"${prefix}.classify.txt\"\n"
    )
    metax.chmod(metax.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    cfg.write_text(
        '{"genomes": {"taxdump": "%s", "taxdump_version": "2024-08-14"}, "tools": {}, "databases": {}}\n'
        % dump
    )
    r1 = tmp_path / "a_1.fastq"
    _tiny_fastq(r1)
    out = tmp_path / "sample.metax.out"
    env = os.environ.copy()
    env["PATH"] = f"{mock_bin}{os.pathsep}{env.get('PATH', '')}"
    payload = json.loads(cfg.read_text())
    payload["compilers"] = {"python": sys.executable}
    cfg.write_text(json.dumps(payload))
    env["SAMOVAR_CONFIG"] = str(cfg)
    env["PYTHONPATH"] = f"{REPO / 'src'}{os.pathsep}{env.get('PYTHONPATH', '')}"
    env.pop("METAX_DMP_DIR", None)
    proc = subprocess.run(
        [
            "bash",
            str(ANNOTATORS / "custom.sh"),
            "-i",
            str(r1),
            "-d",
            str(db),
            "-o",
            str(out),
            "-p",
            "metax",
            "-t",
            "1",
            *separator,
            "--min-read-len", "80", "--sequencer", "Nanopore",
            "--threads=99", "-t", "88", "--db", "/wrong", "--is-paired",
            "--outprefix=/wrong", "--native-label", "literal $(no eval)",
        ],
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    argv = log.read_text()
    assert f"--dmp-dir\n{dump}" in argv or f"--dmp-dir {dump}" in argv
    assert (db / "taxonomy").is_symlink()
    assert (db / "taxonomy").resolve() == dump.resolve()
    assert "read1\t562" in out.read_text()

    tokens = log.read_text().splitlines()
    assert tokens.count('--threads') == 1 and tokens[tokens.index('--threads') + 1] == '1'
    assert tokens.count('--db') == 1 and '/wrong' not in tokens
    assert '--is-paired' not in tokens
    assert '--min-read-len' in tokens and '80' in tokens
    assert tokens.count('--sequencer') == 1 and 'Nanopore' in tokens and 'Illumina' not in tokens
    assert tokens[-2:] == ['--native-label', 'literal $(no eval)']
