import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from samovar.table2iss import (
    generate_reads_genome,
    generate_reads_metagenome,
    regenerate_metagenome,
    split_metagenome_to_samples,
    split_paired_fastq_by_counts,
    write_combined_genomes,
    write_readcount_file,
)

def create_test_fasta(output_dir, name="test"):
    """Create a test FASTA file."""
    fasta_path = os.path.join(output_dir, f"{name}.fa")
    seq = Seq("ATCG" * 100)  # 400 bp sequence
    record = SeqRecord(seq, id=f"{name}_seq", description="")
    SeqIO.write(record, fasta_path, "fasta")
    return fasta_path

def _fastq_record(read_id, seq="ACGT"):
    return f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n"

def _count_fastq_records(path):
    with open(path) as handle:
        return sum(1 for line in handle if line.startswith("@"))

def test_generate_reads_genome():
    """Test generating reads from a single genome."""
    output_dir = "tests_outs/test_generate_reads_genome"
    os.makedirs(output_dir, exist_ok=True)

    fasta_path = create_test_fasta(output_dir)
    output_path = os.path.join(output_dir, "test_reads")

    generate_reads_genome(
        genome_file=fasta_path,
        output_file=output_path,
        amount=20,
        read_length=150,
        cpus=1,
        seed=1,
    )

    assert os.path.exists(f"{output_path}_R1.fastq")
    assert os.path.exists(f"{output_path}_R2.fastq")

def test_generate_reads_metagenome():
    """Test generating a mixed metagenome from multiple genomes in one ISS call."""
    output_dir = "tests_outs/test_generate_reads_metagenome"
    os.makedirs(output_dir, exist_ok=True)

    fasta_paths = [create_test_fasta(output_dir, f"g{i}") for i in range(3)]

    generate_reads_metagenome(
        genome_files=fasta_paths,
        output_dir=output_dir,
        amount=[20, 20, 20],
        read_length=150,
        sample_name="test_metagenome",
        genome_ids=["g0", "g1", "g2"],
        cpus=1,
        seed=1,
    )

    r1 = os.path.join(output_dir, "test_metagenome_any_R1.fastq")
    r2 = os.path.join(output_dir, "test_metagenome_any_R2.fastq")
    assert os.path.exists(r1)
    assert os.path.exists(r2)
    assert _count_fastq_records(r1) > 0
    assert _count_fastq_records(r1) == _count_fastq_records(r2)

def test_regenerate_metagenome():
    """Test regenerating metagenome reads."""
    output_dir = "tests_outs/test_regenerate_metagenome"
    os.makedirs(output_dir, exist_ok=True)

    fasta_paths = [create_test_fasta(output_dir, f"g{i}") for i in range(3)]

    regenerate_metagenome(
        genome_files=fasta_paths,
        output_dir=output_dir,
        amount=[20, 20, 20],
        read_length=150,
        sample_name="test_regenerated",
        genome_ids=["g0", "g1", "g2"],
        cpus=1,
        seed=1,
    )

    assert os.path.exists(os.path.join(output_dir, "test_regenerated_any_R1.fastq"))
    assert os.path.exists(os.path.join(output_dir, "test_regenerated_any_R2.fastq"))

def test_write_combined_genomes_rewrites_headers(tmp_path):
    fasta_a = create_test_fasta(tmp_path, "562")
    fasta_b = create_test_fasta(tmp_path, "9606")
    dest = tmp_path / "combined.fasta"
    ids = write_combined_genomes([str(fasta_a), str(fasta_b)], str(dest), ["562", "9606"])
    assert ids == ["taxid:562", "taxid:9606"]
    text = dest.read_text()
    assert text.startswith(">taxid:562\n")
    assert ">taxid:9606\n" in text
    assert "test_seq" not in text

def test_write_readcount_file(tmp_path):
    dest = tmp_path / "counts.txt"
    total = write_readcount_file(str(dest), ["562", "9606"], [100, 0])
    assert total == 100
    assert dest.read_text() == "562\t100\n"

def test_split_paired_fastq_by_counts(tmp_path):
    r1 = tmp_path / "full_R1.fastq"
    r2 = tmp_path / "full_R2.fastq"
    r1.write_text("".join(_fastq_record(f"r{i}/1") for i in range(10)))
    r2.write_text("".join(_fastq_record(f"r{i}/2") for i in range(10)))

    out = {
        "1": (str(tmp_path / "1_R1.fastq"), str(tmp_path / "1_R2.fastq")),
        "2": (str(tmp_path / "2_R1.fastq"), str(tmp_path / "2_R2.fastq")),
    }
    written = split_paired_fastq_by_counts(str(r1), str(r2), {"1": 6, "2": 4}, out)
    assert written["1"] + written["2"] == 10
    assert _count_fastq_records(out["1"][0]) == written["1"]
    assert _count_fastq_records(out["1"][1]) == written["1"]
    assert _count_fastq_records(out["2"][0]) == written["2"]

def test_split_metagenome_to_samples_preserves_sources(tmp_path):
    r1 = tmp_path / "full_k1_R1.fastq"
    r2 = tmp_path / "full_k1_R2.fastq"
    records_r1 = []
    records_r2 = []
    for source, n in (("562", 6), ("9606", 4)):
        for i in range(n):
            records_r1.append(_fastq_record(f"{source}_{i}_0/1"))
            records_r2.append(_fastq_record(f"{source}_{i}_0/2"))
    r1.write_text("".join(records_r1))
    r2.write_text("".join(records_r2))

    split_metagenome_to_samples(
        str(r1),
        str(r2),
        {
            "s1": {"562": 4, "9606": 1},
            "s2": {"562": 2, "9606": 3},
        },
        str(tmp_path),
        "k1",
    )

    s1_r1 = (tmp_path / "s1_k1_R1.fastq").read_text()
    s2_r1 = (tmp_path / "s2_k1_R1.fastq").read_text()
    assert s1_r1.count("@562_") == 4
    assert s1_r1.count("@9606_") == 1
    assert s2_r1.count("@562_") == 2
    assert s2_r1.count("@9606_") == 3
    assert _count_fastq_records(tmp_path / "s1_k1_R2.fastq") == 5
    assert _count_fastq_records(tmp_path / "s2_k1_R2.fastq") == 5
