from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from samovar.table2iss import (
    generate_iss_test_samples,
    generate_reads_metagenome,
    process_annotation_tables,
)


def _write_fasta(path, seq_id, seq="ATCG" * 100):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(SeqRecord(Seq(seq), id=seq_id, description=""), path, "fasta")
    return str(path)


def _count_fastq_records(path):
    with open(path) as handle:
        return sum(1 for line in handle if line.startswith("@"))


def _is_iss_cmd(cmd):
    return isinstance(cmd, (list, tuple)) and cmd and Path(str(cmd[0])).name == "iss"


def test_generate_reads_metagenome_single_iss_call(tmp_path):
    genomes = [
        _write_fasta(tmp_path / "genomes" / "562.fa", "orig562"),
        _write_fasta(tmp_path / "genomes" / "9606.fa", "orig9606"),
    ]
    out = tmp_path / "reads"
    calls = []
    real_run = __import__("subprocess").run

    def spy(cmd, **kwargs):
        calls.append(cmd)
        return real_run(cmd, **kwargs)

    with patch("samovar.table2iss.subprocess.run", side_effect=spy):
        generate_reads_metagenome(
            genome_files=genomes,
            output_dir=str(out),
            amount=[8, 4],
            sample_name="merged",
            annotator_name="k1",
            genome_ids=["562", "9606"],
            cpus=1,
            seed=1,
        )

    iss_calls = [c for c in calls if _is_iss_cmd(c)]
    assert len(iss_calls) == 1
    assert "--readcount_file" in iss_calls[0]
    assert iss_calls[0].count("--genomes") == 1
    r1 = out / "merged_k1_R1.fastq"
    assert r1.exists()
    text = r1.read_text()
    assert "taxid:562" in text
    assert "taxid:9606" in text


def test_process_annotation_tables_generates_once_then_splits(tmp_path):
    genome_dir = tmp_path / "genomes"
    _write_fasta(genome_dir / "562.fa", "ecoli")
    _write_fasta(genome_dir / "9606.fa", "human")

    tables = tmp_path / "ann"
    tables.mkdir()
    pd.DataFrame({
        "seqID": ["a", "b", "c", "d"],
        "taxID_k1_0": ["562", "562", "9606", "0"],
        "taxID_k2_0": ["562", "9606", "9606", "562"],
    }).to_csv(tables / "s1.annotation.csv", index=False)
    pd.DataFrame({
        "seqID": ["a", "b"],
        "taxID_k1_0": ["562", "0"],
        "taxID_k2_0": ["9606", "9606"],
    }).to_csv(tables / "s2.annotation.csv", index=False)

    reads = tmp_path / "reads"
    calls = []
    real_run = __import__("subprocess").run

    def spy(cmd, **kwargs):
        calls.append(cmd)
        return real_run(cmd, **kwargs)

    with patch("samovar.table2iss.subprocess.run", side_effect=spy):
        process_annotation_tables(
            table_paths=[
                str(tables / "s1.annotation.csv"),
                str(tables / "s2.annotation.csv"),
            ],
            genome_dir=str(genome_dir),
            output_dir=str(reads),
        )

    iss_calls = [c for c in calls if _is_iss_cmd(c)]
    assert len(iss_calls) == 2
    for annotator in ("k1", "k2"):
        for sample in ("s1", "s2"):
            r1 = reads / f"{sample}_{annotator}_R1.fastq"
            r2 = reads / f"{sample}_{annotator}_R2.fastq"
            assert r1.exists(), r1
            assert r2.exists(), r2
            assert _count_fastq_records(r1) == _count_fastq_records(r2)

    s1_k1 = (reads / "s1_k1_R1.fastq").read_text()
    s2_k1 = (reads / "s2_k1_R1.fastq").read_text()
    assert "taxid:562" in (s1_k1 + s2_k1)
    assert not (reads / ".iss_full").exists()


def test_generate_iss_test_samples_two_iss_calls(tmp_path):
    genome_dir = tmp_path / "meta"
    host = _write_fasta(tmp_path / "host" / "9606.fna", "human")
    _write_fasta(genome_dir / "562.fna", "ecoli")
    _write_fasta(genome_dir / "4932.fna", "yeast")
    out = tmp_path / "initial"
    calls = []
    real_run = __import__("subprocess").run

    def spy(cmd, **kwargs):
        calls.append(cmd)
        return real_run(cmd, **kwargs)

    with patch("samovar.table2iss.subprocess.run", side_effect=spy):
        outputs = generate_iss_test_samples(
            genome_dir=str(genome_dir),
            host_genome=host,
            output_dir=str(out),
            n_samples=2,
            total_reads=16,
            host_fraction=0.25,
            seed=1,
            cpus=1,
        )

    iss_calls = [c for c in calls if _is_iss_cmd(c)]
    assert len(iss_calls) == 2  # host pool + one metagenome pool
    assert len(outputs) == 4
    for sample in ("1", "2"):
        r1 = out / f"{sample}_full_R1.fastq"
        r2 = out / f"{sample}_full_R2.fastq"
        assert r1.exists()
        assert r2.exists()
        assert _count_fastq_records(r1) == _count_fastq_records(r2)
    assert not (out / ".iss_full").exists()


def test_generate_reads_metagenome_raises_on_iss_failure(tmp_path):
    genomes = [_write_fasta(tmp_path / "genomes" / "562.fa", "orig562")]

    def fail(cmd, **kwargs):
        class Result:
            returncode = 1

        return Result()

    with patch("samovar.table2iss.subprocess.run", side_effect=fail):
        with pytest.raises(RuntimeError, match="ISS command failed"):
            generate_reads_metagenome(
                genome_files=genomes,
                output_dir=str(tmp_path / "reads"),
                amount=[10],
                sample_name="merged",
                annotator_name="k1",
                genome_ids=["562"],
            )
