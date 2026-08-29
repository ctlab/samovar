"""Prepare/exec runs for the CAMI table contract and regenerated from-genomes."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from samovar.combine_tables import combine_with_cpp
from samovar.ground_truth import iter_truth_pairs
from samovar.paths import repo_root
from samovar.paths import test_genomes_dir as bundled_genomes_dir
from samovar.seqio import link_or_copy_reads

REPO = repo_root()
CAMI_TOY = REPO / "tests" / "data" / "cami_toy"
MAPPING = CAMI_TOY / "reads_mapping.tsv.gz"
MAPPING_PLAIN = CAMI_TOY / "reads_mapping.tsv"


def _run(cmd, cwd, env, timeout=900):
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


def _env():
    env = os.environ.copy()
    env["PATH"] = str(REPO / "bin") + os.pathsep + env.get("PATH", "")
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("NCBI_EMAIL", "test@samovar.com")
    return env


def _require_cami_toy():
    mapping = MAPPING if MAPPING.is_file() else MAPPING_PLAIN
    fastqs = list(CAMI_TOY.glob("sample_*_1.fastq.gz")) or list(
        CAMI_TOY.glob("sample_*_1.fastq")
    )
    if len(fastqs) < 3 or not mapping.is_file():
        pytest.skip(f"CAMI toy fixture missing under {CAMI_TOY}")
    return fastqs, mapping


def test_cami_toy_initial_truth_table_parse_genome(tmp_path):
    """CAMI toy FASTQ + mapping table; regenerated combine stays header-parse."""
    _fastqs, mapping = _require_cami_toy()
    out = tmp_path / "cami_gt"
    env = _env()
    samovar = str(REPO / "bin" / "samovar")
    _run(
        [
            samovar,
            "prepare",
            "--input_dir",
            str(CAMI_TOY),
            "--output_dir",
            str(out),
            "--dummy",
            "dummy",
            "--initial-ground-truth-table",
            str(mapping),
            "--regenerated-metagenomes",
            "parse-genome",
            "--endpoint",
            "combine_initial",
            "--no-multiqc",
        ],
        cwd=REPO,
        env=env,
    )
    pipeline = (out / ".log" / "samovar.sh").read_text()
    assert f'--truth-table "{mapping.resolve()}"' in pipeline
    assert "from-genomes" not in pipeline
    assert 'SAMOVAR_REGENERATED_METAGENOMES="parse-genome"' in pipeline

    _run([samovar, "exec", "--output_dir", str(out)], cwd=REPO, env=env, timeout=900)

    tables = list((out / "initial_annotations").glob("*.annotation.csv"))
    assert len(tables) >= 3
    truth = {}
    for seq, tax in iter_truth_pairs(mapping):
        truth[seq] = tax
        if "/" in seq:
            truth[seq.rsplit("/", 1)[0]] = tax
    matched = 0
    for path in tables:
        df = pd.read_csv(path)
        assert "true" in df.columns
        for seq, val in zip(df["seq"].astype(str), df["true"].astype(str)):
            expect = truth.get(seq) or truth.get(f"{seq}/1")
            if expect:
                assert str(int(float(val))) == str(int(float(expect)))
                matched += 1
    assert matched >= 1000


def _norm_true(df: pd.DataFrame) -> pd.Series:
    true = df["true"]
    out = []
    for val in true:
        text = str(val).strip()
        if text in {"", "nan", "None", "0"}:
            out.append("")
            continue
        try:
            out.append(str(int(float(text))))
        except (TypeError, ValueError):
            out.append(text)
    return pd.Series(out, index=df["seq"].astype(str), name="true")


def test_toy_regenerated_ground_truth_table(tmp_path):
    """examples/toy analog: both regenerated-metagenomes modes assign the same true taxids."""
    genomes = bundled_genomes_dir()
    meta = genomes / "meta"
    host = genomes / "host" / "9606.fna"
    gen = tmp_path / "gen"
    env = _env()
    samovar = str(REPO / "bin" / "samovar")
    _run(
        [
            samovar,
            "generate",
            "--genome_dir",
            str(meta),
            "--host_genome",
            str(host),
            "--output_dir",
            str(gen),
            "--n_samples",
            "1",
            "--total_reads",
            "40",
            "--cores",
            "1",
        ],
        cwd=REPO,
        env=env,
    )
    generate_sh = gen / ".generate" / "generate.sh"
    assert generate_sh.is_file()
    _run(["bash", str(generate_sh)], cwd=REPO, env=env, timeout=300)
    assert any((gen / "initial").glob("*_R1.fastq*"))
    trues = {}
    reports_by_mode = {}
    for mode in ("parse-genome", "ground-truth-table"):
        dest = tmp_path / mode.replace("-", "_")
        _run(
            [
                samovar,
                "prepare",
                "--input_dir",
                str(gen / "initial"),
                "--output_dir",
                str(dest),
                "--test-genomes",
                "--dummy",
                "dummy",
                "--regenerated-metagenomes",
                mode,
                "--N_reads",
                "40",
                "--endpoint",
                "combine_regenerated",
                "--scoring",
                "none",
                "--no-multiqc",
            ],
            cwd=REPO,
            env=env,
        )
        pipeline = (dest / ".log" / "samovar.sh").read_text()
        if mode == "ground-truth-table":
            assert "from-genomes" in pipeline
            assert "${SAMOVAR_INJECT_TAXID:-0}" in pipeline
            initial_block = pipeline.split("ckpt_skip combine_initial")[1].split(
                "ckpt_skip viz_initial"
            )[0]
            assert "--truth-table" not in initial_block
        else:
            assert "from-genomes" not in pipeline
            assert "${SAMOVAR_INJECT_TAXID:-1}" in pipeline
        link_or_copy_reads(gen / "initial", dest / "initial")
        _run(
            [samovar, "exec", "--output_dir", str(dest)],
            cwd=REPO,
            env=env,
            timeout=1200,
        )
        regen = list((dest / "regenerated_annotations").glob("*.annotation.csv"))
        assert regen
        trues[mode] = _norm_true(pd.read_csv(regen[0]))
        reports_by_mode[mode] = dest / "regenerated_reports"
        assert (trues[mode] != "").any()

    parse_t = trues["parse-genome"]
    table_t = trues["ground-truth-table"]
    assert sorted(parse_t.tolist()) == sorted(table_t.tolist())
    common = parse_t.index.intersection(table_t.index)
    if len(common) == len(parse_t) == len(table_t) and len(common) > 0:
        pd.testing.assert_series_equal(
            parse_t.loc[common].sort_index(),
            table_t.loc[common].sort_index(),
            check_names=False,
        )

    table_reports = reports_by_mode["ground-truth-table"]
    parse_recombine = tmp_path / "recombine_parse"
    combine_with_cpp(str(table_reports), str(parse_recombine), 2)
    recom = _norm_true(pd.read_csv(next(parse_recombine.glob("*.annotation.csv"))))
    pd.testing.assert_series_equal(
        table_t.sort_index(),
        recom.loc[table_t.index].sort_index(),
        check_names=False,
    )

    gt_path = tmp_path / "ground_truth_table" / ".log" / "regenerated_ground_truth.tsv"
    assert gt_path.is_file()
    assert gt_path.read_text().splitlines()[0] == "seq\ttaxid"
