"""Generate → prepare (dummy annotator) → exec integrity.

B21: both generate contracts must run end-to-end:

- known ground truth: bundled ``data/test_genomes`` keep ``taxid:`` in headers
- unknown ground truth: a copy of those genomes with true taxid stripped
  (headers and numeric / ``Ecoli.fna``-style names)
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd

from samovar.parse_annotators import extract_true_taxid
from samovar.paths import repo_root
from samovar.paths import test_genomes_dir as bundled_genomes_dir
from samovar.seqio import is_fasta_name, open_text

KNOWN_TRUE_TAXIDS = {"562", "4932", "9606", "2886930"}


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


def _env(root: Path):
    env = os.environ.copy()
    env["PATH"] = str(root / "bin") + os.pathsep + env.get("PATH", "")
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("NCBI_EMAIL", "test@samovar.com")
    return env


def _true_tokens(series: pd.Series) -> pd.Series:
    out = []
    for val in series:
        if val is None or (isinstance(val, float) and pd.isna(val)):
            out.append("")
            continue
        text = str(val).strip()
        if text.lower() in {"", "nan", "none", "<na>"}:
            out.append("")
            continue
        try:
            num = str(int(float(text)))
        except (TypeError, ValueError):
            out.append("")
            continue
        out.append("" if num == "0" else num)
    return pd.Series(out, index=series.index)


def _rewrite_fasta_without_true_taxid(src: Path, dest: Path, contig_prefix: str) -> None:
    """Copy nucleotide FASTA, dropping taxid: tokens and recognisable prefixes."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(src) as inn, open_text(dest, "wt") as out:
        for line in inn:
            if line.startswith(">"):
                n += 1
                out.write(f">{contig_prefix}_c{n}\n")
            else:
                out.write(line if line.endswith("\n") else line + "\n")


def copy_test_genomes_without_true_taxid(dest: Path) -> tuple[Path, Path]:
    """Return ``(meta_dir, host_fasta)`` with no recoverable true taxid."""
    genomes = bundled_genomes_dir()
    meta_src = genomes / "meta"
    host_src = genomes / "host" / "9606.fna"
    meta_out = dest / "meta"
    host_out = dest / "host"
    meta_out.mkdir(parents=True)
    host_out.mkdir(parents=True)
    for i, path in enumerate(sorted(meta_src.iterdir())):
        if not path.is_file() or not is_fasta_name(path.name, protein=False):
            continue
        _rewrite_fasta_without_true_taxid(path, meta_out / f"org{i}.fna", f"org{i}")
    host_fa = host_out / "host.fna"
    _rewrite_fasta_without_true_taxid(host_src, host_fa, "host")
    return meta_out, host_fa


def _prepare_exec(
    root: Path,
    env: dict,
    out: Path,
    meta: Path,
    host: Path,
) -> None:
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
    _run(
        [
            samovar,
            "prepare",
            "--output_dir",
            str(out),
            "--test-genomes",
            "--dummy",
            "dummy",
            "--N_reads",
            "40",
            "--scoring",
            "none",
            "--no-multiqc",
        ],
        cwd=root,
        env=env,
    )
    _run([samovar, "exec", "--output_dir", str(out)], cwd=root, env=env, timeout=900)


def _assert_core_outputs(out: Path) -> Path:
    initial = out / "initial"
    r1 = list(initial.glob("*_R1.fastq")) + list(initial.glob("*_R1.fastq.gz"))
    assert r1, f"ISS did not write reads under {initial}"
    reports = out / "initial_reports"
    assert reports.is_dir()
    assert list(reports.glob("*.csv")) or list(reports.glob("*.out"))
    tables = list((out / "initial_annotations").glob("*.annotation.csv"))
    assert tables, f"no combined tables under {out / 'initial_annotations'}"
    return tables[0]


def test_copy_test_genomes_strips_taxid_tokens(tmp_path):
    meta, host = copy_test_genomes_without_true_taxid(tmp_path / "anon")
    files = list(meta.glob("*.fna"))
    assert files
    assert host.is_file()
    for path in files + [host]:
        assert not path.stem.isdigit()
        first = path.open().readline()
        assert first.startswith(">")
        assert "taxid:" not in first.lower()
        token = first[1:].split()[0]
        assert extract_true_taxid(token) == ""
        assert extract_true_taxid(f"{token}_0_0") == ""


def test_integrity_known_ground_truth_after_generate(tmp_path):
    """Bundled test genomes keep taxid: in headers; combine fills ``true``."""
    root = repo_root()
    genomes = bundled_genomes_dir()
    meta = genomes / "meta"
    host = genomes / "host" / "9606.fna"
    assert meta.is_dir(), f"missing {meta}"
    assert host.is_file(), f"missing {host}"
    out = tmp_path / "samovar_known_gt"
    _prepare_exec(root, _env(root), out, meta, host)
    table = _assert_core_outputs(out)
    df = pd.read_csv(table)
    filled = _true_tokens(df["true"])
    assert (filled != "").mean() > 0.9
    assert set(filled[filled != ""]) <= KNOWN_TRUE_TAXIDS
    json_dir = out / ".log" / "multiqc"
    assert (json_dir / "setup_reads.samovar.json").is_file()


def test_integrity_unknown_ground_truth_stripped_genomes(tmp_path):
    """Same generate/prepare/exec path after stripping true taxid from genomes."""
    root = repo_root()
    meta, host = copy_test_genomes_without_true_taxid(tmp_path / "anon_genomes")
    sample = next(meta.glob("*.fna"))
    header = sample.read_text().splitlines()[0]
    assert "taxid:" not in header.lower()
    assert extract_true_taxid(header.lstrip(">")) == ""
    assert extract_true_taxid(f"{header.lstrip('>')}_0_0") == ""

    out = tmp_path / "samovar_unknown_gt"
    _prepare_exec(root, _env(root), out, meta, host)
    table = _assert_core_outputs(out)
    df = pd.read_csv(table)
    filled = _true_tokens(df["true"])
    assert (filled == "").all(), f"expected unknown true, got {sorted(set(filled))}"
