"""Minimal reads/metagenome generator for contract tests."""

from pathlib import Path


def generate(spec, metadata, config):
    _ = metadata, config
    out = Path(spec["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    r1 = out / "1_full_R1.fastq"
    r2 = out / "1_full_R2.fastq"
    rec = "@read0 taxid:562\nACGT\n+\nIIII\n"
    r1.write_text(rec)
    r2.write_text(rec)
    return [str(r1), str(r2)]
