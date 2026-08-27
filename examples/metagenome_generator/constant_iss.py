#!/usr/bin/env python3
"""Custom combined metagenome generator: constant (equal) abundance + ISS.

Import::

    samovar tools import -n constant_iss --exec-path examples/metagenome_generator/constant_iss.py \\
        --type meta --flags "--model hiseq"

Then ``samovar generate --metagenome_generator constant_iss`` /
``samovar prepare --metagenome_generator constant_iss``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from samovar.metagenome_generators import constant_abundance_frame, parse_constant_flags
from samovar.table2iss import (
    generate_iss_from_abundance_table,
    generate_iss_test_samples,
    iss_cli_extra_flags,
)


def generate(spec: dict, metadata: Optional[pd.DataFrame], config: dict) -> List[str]:
    _ = metadata
    cfg = dict(config or {})
    parsed = parse_constant_flags(list(cfg.get("extra_argv") or []))
    output_dir = str(spec.get("output_dir") or cfg.get("output_dir"))
    genome_dir = str(spec.get("genome_dir") or cfg.get("genome_dir") or "")
    host_genome = str(spec.get("host_genome") or cfg.get("host_genome") or "")
    host_fraction = spec.get("host_fraction") or cfg.get("host_fraction") or "RANDOM"
    n_samples = int(spec.get("n_samples") or cfg.get("n_samples") or 10)
    total_reads = int(
        parsed.get("n_reads")
        or spec.get("total_reads")
        or cfg.get("total_reads")
        or cfg.get("N_reads")
        or 2000
    )
    seed = int(cfg.get("seed") or 42)
    model = str(parsed.get("model") or cfg.get("model") or "hiseq")
    cpus = int(cfg.get("cores") or cfg.get("cpus") or 2)
    extra = list(parsed.get("iss_extra") or [])
    abundance = spec.get("abundance_table") or cfg.get("abundance_table") or cfg.get("abundance")
    stage = str(spec.get("stage") or cfg.get("stage") or "generate")

    with iss_cli_extra_flags(extra):
        if abundance and Path(str(abundance)).is_file():
            table = pd.read_csv(abundance)
            flat = constant_abundance_frame(table, n_reads=total_reads)
            tmp = Path(output_dir) / ".constant_abundance.csv"
            tmp.parent.mkdir(parents=True, exist_ok=True)
            flat.to_csv(tmp, index=False)
            hf = host_fraction if stage == "generate" else 0
            return generate_iss_from_abundance_table(
                str(tmp),
                genome_dir=genome_dir,
                output_dir=output_dir,
                host_genome=host_genome,
                host_fraction=hf,
                total_reads=total_reads,
                seed=seed,
                model=model,
                cpus=cpus,
                annotator_name=str(spec.get("annotator") or cfg.get("annotator_name") or "full"),
            )
        return generate_iss_test_samples(
            genome_dir=genome_dir,
            host_genome=host_genome,
            output_dir=output_dir,
            n_samples=n_samples,
            total_reads=total_reads,
            host_fraction=host_fraction,
            seed=seed,
            model=model,
            genomes=cfg.get("genomes"),
            cpus=cpus,
            extra_flags=extra,
        )


def _cli(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="constant abundance + ISS metagenome generator")
    p.add_argument("-i", dest="abundance", default=None)
    p.add_argument("-o", dest="output_dir", required=True)
    p.add_argument("-m", dest="metadata", default=None)
    p.add_argument("--genome-dir", dest="genome_dir", default="")
    p.add_argument("--host-genome", dest="host_genome", default="")
    p.add_argument("--host-fraction", dest="host_fraction", default="RANDOM")
    p.add_argument("--n-samples", dest="n_samples", type=int, default=10)
    p.add_argument("--total-reads", dest="total_reads", type=int, default=2000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--model", default="hiseq")
    p.add_argument("--gzip-reads", action="store_true")
    p.add_argument("--stage", default="generate")
    p.add_argument("--annotator", default=None)
    p.add_argument("--read-type", action="append", default=None)
    args, rest = p.parse_known_args(argv)
    meta = pd.read_csv(args.metadata) if args.metadata else None
    spec = {
        "abundance_table": args.abundance,
        "output_dir": args.output_dir,
        "genome_dir": args.genome_dir,
        "host_genome": args.host_genome,
        "host_fraction": args.host_fraction,
        "n_samples": args.n_samples,
        "total_reads": args.total_reads,
        "stage": args.stage,
        "annotator": args.annotator,
        "gzip_reads": args.gzip_reads,
        "extra_ids": args.read_type or [],
    }
    config: Dict[str, Any] = {
        "seed": args.seed,
        "model": args.model,
        "extra_argv": rest,
        "N_reads": args.total_reads,
    }
    generate(spec, meta, config)
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
