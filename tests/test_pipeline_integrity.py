"""End-to-end integrity: generate → prepare → exec like samovar_fix2, on bundled genomes.

Uses the real ``iss`` CLI and the dummy annotator (no Kraken/Kaiju index). NCBI is
not required: genomes come from ``data/test_genomes``.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from samovar.paths import repo_root
from samovar.paths import test_genomes_dir as bundled_genomes_dir


def _run(cmd, cwd, env, timeout=600):
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


def test_fix2_style_generate_prepare_exec(tmp_path):
    root = repo_root()
    genomes = bundled_genomes_dir()
    meta = genomes / "meta"
    host = genomes / "host" / "9606.fna"
    assert meta.is_dir(), f"missing {meta}"
    assert host.is_file(), f"missing {host}"

    out = tmp_path / "samovar_fix2_ci"
    env = os.environ.copy()
    env["PATH"] = str(root / "bin") + os.pathsep + env.get("PATH", "")
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("NCBI_EMAIL", "test@samovar.com")
    samovar = str(root / "bin" / "samovar")

    _run(
        [
            samovar,
            "generate",
            "--genome_dir",
            str(meta),
            "--host_genome",
            str(host),
            "--output_dir",
            str(out),
            "--n_samples",
            "1",
            "--total_reads",
            "40",
            "--cores",
            "1",
        ],
        cwd=root,
        env=env,
    )
    assert (out / ".generate" / "generate.sh").is_file()

    _run(
        [
            samovar,
            "prepare",
            "--output_dir",
            str(out),
            "--dummy",
            "dummy",
        ],
        cwd=root,
        env=env,
    )
    assert (out / ".log" / "samovar.sh").is_file()

    _run(
        [samovar, "exec", "--output_dir", str(out)],
        cwd=root,
        env=env,
        timeout=900,
    )

    initial = out / "initial"
    r1 = list(initial.glob("*_R1.fastq")) + list(initial.glob("*_R1.fastq.gz"))
    assert r1, f"ISS did not write reads under {initial}"
    reports = out / "initial_reports"
    assert reports.is_dir()
    assert list(reports.glob("*.csv")) or list(reports.glob("*.out"))
    json_dir = out / ".log" / "multiqc"
    assert (json_dir / "setup_reads.samovar.json").is_file()
    staged = out / "multiqc_samovar"
    assert staged.is_dir()
    assert (staged / "00_SamovaR_pipeline_mqc.html").is_file()
