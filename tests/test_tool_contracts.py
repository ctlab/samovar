"""In→out contracts for each ``samovar tools import`` type.

Default run uses bundled example tools. Pass ``--tool PATH --tool-type GROUP``
(or env SAMOVAR_CONTRACT_TOOL / SAMOVAR_CONTRACT_TYPE) to check a custom dest.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pandas as pd
import pytest

from samovar.abundance import n_sample_columns, normalize_abundance_table
from samovar.main_config import normalize_tool_group
from samovar.paths import repo_root
from samovar.tool_contracts import DEFAULT_TOOLS, load_python_module


def _wanted_type(request) -> str:
    raw = request.config.getoption("--tool-type") or ""
    if not raw:
        import os

        raw = os.environ.get("SAMOVAR_CONTRACT_TYPE") or ""
    return normalize_tool_group(raw) if str(raw).strip() else ""


def _tool_path(request, group: str) -> Path:
    raw = request.config.getoption("--tool")
    if not raw:
        import os

        raw = os.environ.get("SAMOVAR_CONTRACT_TOOL")
    if raw:
        return Path(raw).expanduser().resolve()
    rel = DEFAULT_TOOLS[group]
    return (repo_root() / rel).resolve()


def _skip_if_other_type(request, group: str) -> None:
    wanted = _wanted_type(request)
    if wanted and wanted != group:
        pytest.skip(f"contract type {wanted} != {group}")


def _tiny_abundance():
    return pd.DataFrame(
        {
            "taxid": ["562", "2886930", "9606"],
            "N_1": [10, 4, 2],
            "N_2": [8, 3, 1],
        }
    )


@pytest.fixture
def tiny_fastq(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    rec = "@r0\nACGT\n+\nIIII\n"
    r1.write_text(rec)
    r2.write_text(rec)
    return r1, r2


def test_annotator_contract(request, tmp_path, tiny_fastq):
    _skip_if_other_type(request, "annotator")
    path = _tool_path(request, "annotator")
    assert path.is_file(), f"annotator dest not found: {path}"
    r1, r2 = tiny_fastq
    out = tmp_path / "clf.out"
    if path.suffix.lower() == ".py":
        module = load_python_module(path, "contract_annotator")
        if hasattr(module, "main"):
            rc = module.main(["-i", str(r1), "-I", str(r2), "-d", str(tmp_path), "-o", str(out), "-t", "1"])
            assert rc in (None, 0)
        else:
            out.write_text("r0\t562\n")
        parse = getattr(module, "parse_output", None)
        if callable(parse):
            frame = parse(str(out))
        else:
            frame = pd.read_table(out, header=None)
            frame.columns = ["seq", "taxID"]
    else:
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        import subprocess

        subprocess.check_call(
            [str(path), "-i", str(r1), "-I", str(r2), "-d", str(tmp_path), "-o", str(out), "-t", "1"]
        )
        from samovar.annotators_wrapper import CustomAnnotator

        inst = CustomAnnotator(
            {"run_name": "clf", "type": "clf", "cmd": str(path), "db_path": str(tmp_path), "threads": 1},
            {},
            "clf",
        )
        frame = inst.parse_output(str(out))
    assert "seq" in frame.columns
    assert "taxID" in frame.columns or "taxid" in {c.lower() for c in frame.columns}
    assert not frame.empty


def test_table_regenerator_contract(request, tmp_path):
    _skip_if_other_type(request, "table_reads_generator")
    path = _tool_path(request, "table_reads_generator")
    assert path.suffix.lower() == ".py", "table regenerator contract expects a Python regenerate()"
    module = load_python_module(path, "contract_table")
    fn = getattr(module, "regenerate", None)
    assert callable(fn), f"{path} must define regenerate(data, metadata, config)"
    tables = fn(
        _tiny_abundance(),
        None,
        {"output_dir": str(tmp_path), "seed": 1, "max_genomes": float("inf")},
    )
    assert isinstance(tables, dict) and tables
    table = next(iter(tables.values()))
    assert "taxid" in table.columns
    assert n_sample_columns(normalize_abundance_table(table))
    capped = fn(
        _tiny_abundance(),
        None,
        {"output_dir": str(tmp_path), "seed": 1, "max_genomes": 1},
    )
    capped_table = next(iter(capped.values()))
    assert len(capped_table) <= 1


def test_table_scoring_contract(request):
    _skip_if_other_type(request, "table_scoring")
    path = _tool_path(request, "table_scoring")
    module = load_python_module(path, "contract_table_score")
    obs = _tiny_abundance()
    gen = _tiny_abundance().copy()
    gen["N_1"] = gen["N_1"] + 1
    score_ann = getattr(module, "score_annotator", None)
    score_one = getattr(module, "score_table", None)
    assert callable(score_ann) or callable(score_one), (
        f"{path} needs score_annotator(...) or score_table(observed, generated, config)"
    )
    if callable(score_ann):
        block = score_ann(
            obs,
            {"bootstrap": {"table": gen}},
            "table",
            ["bootstrap"],
            {},
        )
    else:
        block = score_one(obs, gen, {})
    assert isinstance(block, dict)
    if "rank_value" not in block and "candidates" in block:
        cands = block["candidates"]
        assert cands
        assert "rank_value" in cands[0] or "ks_statistic" in cands[0]
    else:
        assert "rank_value" in block or "ks_statistic" in block


def test_scoring_contract(request, tmp_path):
    _skip_if_other_type(request, "scoring")
    path = _tool_path(request, "scoring")
    module = load_python_module(path, "contract_scoring")
    fn = getattr(module, "score", None)
    assert callable(fn), f"{path} must define score(inputs, dest, config)"
    ann = tmp_path / "initial_annotations"
    ann.mkdir()
    (ann / "s.annotation.csv").write_text("seq,taxID_k_0\na,562\n")
    dest = tmp_path / "scores"
    fn([ann], dest, {"stage": "viz_initial"})
    written = [p for p in dest.rglob("*") if p.is_file()] if dest.exists() else []
    assert dest.exists() and written, f"{path} score() must write files under dest"


def test_reads_generator_contract(request, tmp_path):
    _skip_if_other_type(request, "reads_generator")
    path = _tool_path(request, "reads_generator")
    module = load_python_module(path, "contract_reads")
    fn = getattr(module, "generate", None)
    assert callable(fn), f"{path} must define generate(spec, metadata, config)"
    out = tmp_path / "reads"
    paths = fn(
        {
            "output_dir": str(out),
            "abundance_table": str(tmp_path / "t.csv"),
            "genome_dir": str(tmp_path),
            "n_samples": 1,
            "total_reads": 10,
        },
        None,
        {"max_genomes": float("inf")},
    )
    assert isinstance(paths, (list, tuple)) and paths
    assert any(str(p).endswith(".fastq") or str(p).endswith(".fastq.gz") for p in paths)
    assert Path(paths[0]).is_file()


def test_metagenome_generator_contract(request, tmp_path):
    _skip_if_other_type(request, "metagenome_generator")
    path = _tool_path(request, "metagenome_generator")
    module = load_python_module(path, "contract_meta")
    fn = getattr(module, "generate", None)
    assert callable(fn), f"{path} must define generate(spec, metadata, config) like --type reads"
    out = tmp_path / "meta"
    paths = fn(
        {"output_dir": str(out), "n_samples": 1, "total_reads": 10},
        None,
        {"max_genomes": float("inf")},
    )
    assert isinstance(paths, (list, tuple)) and paths
    assert Path(paths[0]).is_file()


def test_reprofiler_contract(request, tmp_path):
    _skip_if_other_type(request, "reprofiler")
    path = _tool_path(request, "reprofiler")
    module = load_python_module(path, "contract_ml")
    fn = getattr(module, "reprofile", None)
    assert callable(fn), f"{path} must define reprofile(regenerated, ground_truth, initial, config)"
    regenerated = pd.DataFrame(
        {
            "seq": [f"r{i}" for i in range(12)],
            "taxid_dummy": [9606, 9606, 562, 562, 9606, 562, 9606, 562, 9606, 562, 9606, 562],
            "length": [50] * 12,
            "true": [9606, 9606, 562, 562, 9606, 562, 9606, 562, 9606, 562, 9606, 562],
        }
    )
    initial = {
        "sample.annotation": pd.DataFrame(
            {"seq": ["a", "b", "c"], "taxid_dummy": [9606, 562, 0], "length": [50, 50, 50]}
        )
    }
    ground = {"dummy": pd.DataFrame({"taxid": [9606, 562], "N_1": [10, 8]})}
    raw = fn(regenerated, ground, initial, {"output_dir": str(tmp_path / "ml"), "seed": 0})
    tables = getattr(raw, "tables", None)
    if tables is None and isinstance(raw, dict):
        tables = raw.get("tables") or raw
    assert tables
    first = next(iter(tables.values())) if isinstance(tables, dict) else tables
    assert isinstance(first, pd.DataFrame)
    assert not first.empty


def test_annotation_converter_contract(request, tmp_path):
    _skip_if_other_type(request, "annotation_converter")
    path = _tool_path(request, "annotation_converter")
    module = load_python_module(path, "contract_annconv")
    dump = getattr(module, "dump", None)
    convert = getattr(module, "convert", None)
    assert callable(dump) or callable(convert), (
        f"{path} must define dump(annotation, dest, config) or convert(src, dest, config)"
    )
    from samovar.parse_annotators import Annotation

    ann = Annotation.from_long_table(
        pd.DataFrame(
            {
                "seq": ["a", "b", "c"],
                "taxID_kaiju_0": ["562", "562", "9606"],
                "sample": ["1", "1", "1"],
            }
        )
    )
    dest = tmp_path / "converted"
    if callable(dump):
        written = dump(ann, dest, {"from": "annotation", "to": "abundance"})
    else:
        written = convert(ann, dest, {"from": "annotation", "to": "abundance"})
    written = Path(written) if written else dest
    assert Path(written).exists() or dest.exists()
    load = getattr(module, "load", None)
    if callable(load):
        roundtrip = load(written if Path(written).exists() else dest, {})
        assert roundtrip is not None


def test_qc_contract(request, tmp_path, tiny_fastq):
    _skip_if_other_type(request, "qc")
    path = _tool_path(request, "qc")
    assert path.is_file(), f"QC dest not found: {path}"
    seq = "ACGT" * 20
    rec = f"@r0\n{seq}\n+\n{'I' * len(seq)}\n"
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    r1.write_text(rec)
    r2.write_text(rec)
    dest_r1 = tmp_path / "out_R1.fastq"
    dest_r2 = tmp_path / "out_R2.fastq"
    if path.suffix.lower() == ".py":
        module = load_python_module(path, "contract_qc")
        fn = getattr(module, "trim", None)
        assert callable(fn), f"{path} must define trim(r1, r2, dest_r1, dest_r2, config)"
        paths = fn(
            str(r1),
            str(r2),
            str(dest_r1),
            str(dest_r2),
            {"min_gc": 0.0, "max_gc": 1.0, "extra_argv": ["--length_required", "1"]},
        )
    else:
        from samovar.qc import apply_qc_executable

        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        paths = apply_qc_executable(
            path,
            r1,
            r2,
            dest_r1,
            dest_r2,
            {
                "extra_argv": ["--length_required", "1"],
                "threads": 1,
            },
            name=path.name,
        )
    assert dest_r1.is_file()
    assert dest_r2.is_file()
    assert Path(paths[0]).is_file()
    text = dest_r1.read_text()
    assert "@" in text
    assert "ACGT" in text

