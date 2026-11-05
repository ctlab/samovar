import gzip
import pytest

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
    drop_r2_if_source_single_end,
    list_fastq_samples,
    validate_fastq_directory,
    open_text,
    processed_genome_path,
    sequence_stem,
    taxid_from_fasta_name,
)


def test_stems_and_taxids():
    assert sequence_stem("562.fna.gz") == "562"
    assert sequence_stem("562-processed.fasta.gz") == "562-processed"
    assert sequence_stem("GCF_000819615.1.fa.gz") == "GCF_000819615.1"
    assert taxid_from_fasta_name("562.fa") == "562"
    assert taxid_from_fasta_name("562.fa.gz") == "562"
    assert taxid_from_fasta_name("562-processed.fasta.gz") == "562"
    assert taxid_from_fasta_name("GCF_000819615.1.fa.gz") == "GCF_000819615.1"
    assert taxid_from_fasta_name("abc.fa") == "abc"


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


def test_find_genome_prefers_accession_fa_gz(tmp_path):
    (tmp_path / "GCF_000819615.1.fa").write_text(">raw\nA\n")
    gz = tmp_path / "GCF_000819615.1.fa.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(">p\nACGT\n")
    found = find_genome_file(tmp_path, "GCF_000819615.1")
    assert found.endswith("GCF_000819615.1.fa.gz")


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
    (tmp_path / "bare.fastq").write_text("@c\nA\n+\nI\n")
    with gzip.open(tmp_path / "compressed.fq.gz", "wt") as handle:
        handle.write("@d\nT\n+\nI\n")
    (tmp_path / "orphan_R2.fastq").write_text("@e\nT\n+\nI\n")
    assert list_fastq_samples(tmp_path) == ["bare", "compressed", "other", "sample_k1"]
    assert find_fastq_mate(tmp_path, "bare", "R1") == tmp_path / "bare.fastq"
    assert find_fastq_mate(tmp_path, "bare", "R2") is None
    assert find_fastq_mate(tmp_path, "sample_k1", "R1") == r1
    assert has_r1_reads(tmp_path)
    dest = tmp_path / "cat.fastq"
    concat_fastq_files([r1, tmp_path / "other_R1.fq"], dest)
    assert dest.read_text().count("@") == 2


def test_has_r1_reads_rejects_empty_and_one_byte_stubs(tmp_path):
    from samovar.seqio import fastq_has_reads

    empty = tmp_path / "empty_R1.fastq"
    empty.write_text("")
    stub = tmp_path / "stub_R1.fastq"
    stub.write_bytes(b"x")
    gz = tmp_path / "compressed_empty_R1.fastq.gz"
    with gzip.open(gz, "wb") as handle:
        handle.write(b"")
    assert not fastq_has_reads(empty)
    assert not fastq_has_reads(stub)
    assert not fastq_has_reads(gz)
    assert not has_r1_reads(tmp_path)
    good = tmp_path / "ok_R1.fastq"
    good.write_text("@a\nA\n+\nI\n")
    assert fastq_has_reads(good)
    assert has_r1_reads(tmp_path)


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


def test_validate_fastq_directory_rejects_blank_and_truncated(tmp_path):
    path = tmp_path / "a_R1.fastq"
    path.write_text("@a\nAC\n+\nII\n\n")
    with pytest.raises(ValueError, match="Blank line in FASTQ"):
        validate_fastq_directory(tmp_path)
    path.write_text("@a\nAC\n")
    with pytest.raises(ValueError, match="Truncated FASTQ record"):
        validate_fastq_directory(tmp_path)
    path.write_text("@a\nAC\n+\nII\n")
    validate_fastq_directory(tmp_path)


@pytest.mark.parametrize("other", ["sample_R1.fq", "sample_R1.fastq.gz", "sample_1.fastq", "sample.fastq", "sample.fq.gz"])
def test_fastq_discovery_rejects_colliding_names(tmp_path, other):
    first = tmp_path / "sample_R1.fastq"
    second = tmp_path / other
    first.write_text("@a\nA\n+\nI\n")
    second.write_text("@b\nT\n+\nI\n")
    with pytest.raises(ValueError, match="Multiple R1 FASTQ files for sample 'sample'") as error:
        list_fastq_samples(tmp_path)
    assert str(first) in str(error.value) and str(second) in str(error.value)


def test_fastq_discovery_rejects_colliding_r2(tmp_path):
    (tmp_path / "s_R1.fastq").write_text("@a\nA\n+\nI\n")
    (tmp_path / "s_R2.fastq").write_text("@a\nA\n+\nI\n")
    (tmp_path / "s_2.fq").write_text("@b\nT\n+\nI\n")
    with pytest.raises(ValueError, match="Multiple R2 FASTQ files for sample 's'"):
        list_fastq_samples(tmp_path)


def test_link_reads_rejects_other_source_before_staging(tmp_path):
    from samovar.seqio import link_or_copy_reads

    first, second, out = (tmp_path / name for name in ("first", "second", "initial"))
    first.mkdir()
    second.mkdir()
    (first / "sample_R1.fastq").write_text("@a\nA\n+\nI\n")
    (second / "sample_R1.fastq").write_text("@b\nT\n+\nI\n")
    (second / "another_R1.fastq").write_text("@c\nT\n+\nI\n")
    assert link_or_copy_reads(first, out) == 1
    assert link_or_copy_reads(first, out) == 0
    with pytest.raises(ValueError, match="Sample name collision"):
        link_or_copy_reads(second, out)
    assert not (out / "another_R1.fastq").exists()
    assert (out / "sample_R1.fastq").resolve() == (first / "sample_R1.fastq").resolve()


def test_link_reads_resumes_copied_files(tmp_path, monkeypatch):
    from pathlib import Path
    from samovar.seqio import link_or_copy_reads

    source, out = tmp_path / "source", tmp_path / "out"
    source.mkdir()
    (source / "sample_R1.fastq").write_text("@a\nA\n+\nI\n")

    def no_symlinks(*args, **kwargs):
        raise OSError("Symlinks unavailable")

    monkeypatch.setattr(Path, "symlink_to", no_symlinks)
    assert link_or_copy_reads(source, out) == 1
    assert link_or_copy_reads(source, out) == 0
    (source / "sample_R1.fastq").write_text("@b\nT\n+\nI\n")
    with pytest.raises(ValueError, match="Sample name collision"):
        link_or_copy_reads(source, out)


def test_fastq_pair_paths():
    assert fastq_pair_paths("/tmp/s_k1") == ("/tmp/s_k1_R1.fastq", "/tmp/s_k1_R2.fastq")
    r1, r2 = fastq_pair_paths("/tmp/s_k1", gzip_reads=True)
    assert r1.endswith("_R1.fastq.gz")
    assert r2.endswith("_R2.fastq.gz")


def test_processed_genome_path():
    assert processed_genome_path("/g", "562").name == "562.fa.gz"
    assert processed_genome_path("/g", "GCF_000819615.1").name == "GCF_000819615.1.fa.gz"
    assert processed_genome_path("/g", "562", gzip_genomes=False).name == "562.fa"


def test_drop_r2_when_source_is_single_end(tmp_path):
    initial = tmp_path / "initial"
    regen = tmp_path / "regenerated"
    initial.mkdir()
    regen.mkdir()
    (initial / "sample_0.fq.gz").write_bytes(b"\x1f\x8b")
    (regen / "sample_0_concat_R1.fastq").write_text("@a\nA\n+\nI\n")
    (regen / "sample_0_concat_R2.fastq").write_text("@a\nT\n+\nI\n")
    (regen / "sample_0_metax_R2.fastq.gz").write_bytes(b"\x1f\x8b")
    removed = drop_r2_if_source_single_end(regen)
    assert (regen / "sample_0_concat.fastq").is_file()
    assert not (regen / "sample_0_concat_R1.fastq").exists()
    assert not (regen / "sample_0_concat_R2.fastq").exists()
    assert not (regen / "sample_0_metax_R2.fastq.gz").exists()
    assert {p.name for p in removed} == {
        "sample_0_concat_R2.fastq",
        "sample_0_metax_R2.fastq.gz",
    }
    paired = tmp_path / "paired"
    (paired / "initial").mkdir(parents=True)
    (paired / "regenerated").mkdir()
    (paired / "initial" / "s_R1.fastq").write_text("@a\nA\n+\nI\n")
    (paired / "initial" / "s_R2.fastq").write_text("@a\nT\n+\nI\n")
    kept = paired / "regenerated" / "s_iss_R2.fastq"
    kept.write_text("@a\nT\n+\nI\n")
    assert drop_r2_if_source_single_end(paired / "regenerated") == []
    assert kept.is_file()


def test_list_fasta_includes_gzip(tmp_path):
    (tmp_path / "1.fna.gz").write_bytes(b"\x1f\x8b")
    (tmp_path / "notes.txt").write_text("nope")
    files = list_fasta_files(tmp_path, nucleotide=True, protein=False)
    assert [p.name for p in files] == ["1.fna.gz"]
