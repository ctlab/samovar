import gzip

from samovar.seqio import (
    concat_fastq_files,
    fastq_pair_paths,
    find_fastq_mate,
    find_genome_file,
    gzip_file,
    gunzip_file,
    has_r1_reads,
    is_fasta_name,
    is_fastq_name,
    list_fasta_files,
    list_fastq_samples,
    open_text,
    processed_genome_path,
    sequence_stem,
    taxid_from_fasta_name,
)


def test_stems_and_taxids():
    assert sequence_stem("562.fna.gz") == "562"
    assert sequence_stem("562-processed.fasta.gz") == "562-processed"
    assert taxid_from_fasta_name("562.fa") == "562"
    assert taxid_from_fasta_name("562.fa.gz") == "562"
    assert taxid_from_fasta_name("562-processed.fasta.gz") == "562"
    assert taxid_from_fasta_name("abc.fa") is None


def test_name_detectors():
    assert is_fasta_name("1.fna.gz", nucleotide=True, protein=False)
    assert is_fasta_name("1.faa.gz", nucleotide=False, protein=True)
    assert not is_fasta_name("1.faa.gz", nucleotide=True, protein=False)
    assert is_fastq_name("s_R1.fastq.gz")
    assert is_fastq_name("s_R1.fq")
    assert not is_fastq_name("s.fasta.gz")


def test_open_text_roundtrip(tmp_path):
    plain = tmp_path / "a.fa"
    plain.write_text(">h\nACGT\n")
    gz = gzip_file(plain)
    assert gz.name.endswith(".gz")
    assert not plain.exists()
    with open_text(gz) as handle:
        assert "ACGT" in handle.read()
    restored = gunzip_file(gz, remove_source=True)
    assert restored.read_text().startswith(">h")


def test_find_genome_prefers_processed_gzip(tmp_path):
    (tmp_path / "562.fa").write_text(">raw\nA\n")
    gz = tmp_path / "562-processed.fasta.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(">p\nACGT\n")
    found = find_genome_file(tmp_path, "562")
    assert found.endswith("562-processed.fasta.gz")


def test_fastq_discovery(tmp_path):
    r1 = tmp_path / "sample_k1_R1.fastq.gz"
    r2 = tmp_path / "sample_k1_R2.fastq.gz"
    with gzip.open(r1, "wt") as handle:
        handle.write("@a\nA\n+\nI\n")
    with gzip.open(r2, "wt") as handle:
        handle.write("@a\nT\n+\nI\n")
    (tmp_path / "other_R1.fq").write_text("@b\nA\n+\nI\n")
    (tmp_path / "other_R2.fq").write_text("@b\nT\n+\nI\n")
    assert list_fastq_samples(tmp_path) == ["other", "sample_k1"]
    assert find_fastq_mate(tmp_path, "sample_k1", "R1") == r1
    assert has_r1_reads(tmp_path)
    dest = tmp_path / "cat.fastq"
    concat_fastq_files([r1, tmp_path / "other_R1.fq"], dest)
    assert dest.read_text().count("@") == 2


def test_fastq_illumina_numeric_mates(tmp_path):
    r1 = tmp_path / "sample_0_1.fastq"
    r2 = tmp_path / "sample_0_2.fastq"
    r1.write_text("@a\nA\n+\nI\n")
    r2.write_text("@a\nT\n+\nI\n")
    gz1 = tmp_path / "run_1.fq.gz"
    gz2 = tmp_path / "run_2.fq.gz"
    with gzip.open(gz1, "wt") as handle:
        handle.write("@b\nA\n+\nI\n")
    with gzip.open(gz2, "wt") as handle:
        handle.write("@b\nT\n+\nI\n")
    samples = list_fastq_samples(tmp_path)
    assert "sample_0" in samples
    assert "run" in samples
    assert find_fastq_mate(tmp_path, "sample_0", "R1") == r1
    assert find_fastq_mate(tmp_path, "sample_0", "R2") == r2
    # ``_R1`` still wins over the ``_1`` suffix of the same filename.
    (tmp_path / "keep_R1.fastq").write_text("@c\nA\n+\nI\n")
    (tmp_path / "keep_R2.fastq").write_text("@c\nT\n+\nI\n")
    assert "keep" in list_fastq_samples(tmp_path)
    assert "keep_R" not in list_fastq_samples(tmp_path)


def test_fastq_ignores_camisim_tech_duplicates(tmp_path):
    (tmp_path / "1_full_R1.fastq").write_text("@a\nA\n+\nI\n")
    (tmp_path / "1_full_R2.fastq").write_text("@a\nT\n+\nI\n")
    (tmp_path / "1_wgsim_R1.fastq").write_text("@b\nA\n+\nI\n")
    (tmp_path / "1_wgsim_R2.fastq").write_text("@b\nT\n+\nI\n")
    (tmp_path / "2_full_R1.fastq").write_text("@c\nA\n+\nI\n")
    assert list_fastq_samples(tmp_path) == ["1_full", "2_full"]


def test_fastq_pair_paths():
    assert fastq_pair_paths("/tmp/s_k1") == ("/tmp/s_k1_R1.fastq", "/tmp/s_k1_R2.fastq")
    r1, r2 = fastq_pair_paths("/tmp/s_k1", gzip_reads=True)
    assert r1.endswith("_R1.fastq.gz")
    assert r2.endswith("_R2.fastq.gz")


def test_processed_genome_path():
    assert processed_genome_path("/g", "562").name == "562-processed.fasta.gz"
    assert processed_genome_path("/g", "562", gzip_genomes=False).name == "562-processed.fasta"


def test_list_fasta_includes_gzip(tmp_path):
    (tmp_path / "1.fna.gz").write_bytes(b"\x1f\x8b")
    (tmp_path / "notes.txt").write_text("nope")
    files = list_fasta_files(tmp_path, nucleotide=True, protein=False)
    assert [p.name for p in files] == ["1.fna.gz"]
