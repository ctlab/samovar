"""External profiler example: unchanged Kraken2 -> direct-count taxid profile.

Import with samovar tools import --type profiler --exec-path THIS_FILE --pytest.
Kraken2 itself must be registered in the active SAMOVAR install config.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

import pandas as pd

from samovar.paths import command_argv, resolve_executable
from samovar.profiles import metadata, validate
from samovar.seqio import fastq_has_reads


def parse_report(path):
    counts = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 6:
            raise ValueError("Malformed Kraken2 report: expected six columns")
        percent, clade, direct, rank, taxid, name = fields
        percentage = float(percent)
        if (not math.isfinite(percentage) or not 0 <= percentage <= 100
                or not clade.isdigit() or not direct.isdigit() or not taxid.isdigit()
                or int(direct) > int(clade) or not rank or not name.strip()):
            raise ValueError("Malformed Kraken2 report values")
        if taxid in counts:
            raise ValueError(f"Duplicate Kraken2 taxid: {taxid}")
        counts[taxid] = int(direct)
    if not counts:
        raise ValueError("Empty Kraken2 report for nonempty reads")
    return pd.DataFrame([(tax, count) for tax, count in counts.items() if count > 0],
                        columns=["taxid", "N_sample"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", required=True)
    parser.add_argument("-I", default="")
    parser.add_argument("-d", required=True)
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", type=int, default=1)
    args = parser.parse_args(argv)
    if args.t < 1:
        parser.error("-t must be positive")
    out = Path(args.o)
    scratch = out.parent / ".tmp"
    scratch.mkdir(parents=True, exist_ok=True)
    command = command_argv(resolve_executable("kraken2", tool_key="kraken2"))
    with tempfile.TemporaryDirectory(prefix="kraken2_profile_", dir=scratch) as tmp:
        if fastq_has_reads(args.i):
            report = Path(tmp) / "kraken2.report"
            reads = ["--paired", args.i, args.I] if args.I and fastq_has_reads(args.I) else [args.i]
            subprocess.run([*command, "--db", args.d, "--threads", str(args.t),
                            "--report", str(report), "--output", "/dev/null", *reads], check=True)
            table = parse_report(report)
        else:
            table = pd.DataFrame(columns=["taxid", "N_sample"])
        table, info = validate(table, metadata(dict(rank="taxid", normalization="all")))
        profile = Path(tmp) / "profile.csv"
        sidecar = Path(str(profile) + ".profile.json")
        table.to_csv(profile, index=False)
        sidecar.write_text(json.dumps(info, indent=2) + "\n")
        profile.replace(out)
        sidecar.replace(str(out) + ".profile.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
