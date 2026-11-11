"""Deterministic sample-profile baseline for the profiler contract."""
import argparse
import json
from pathlib import Path

import pandas as pd

from samovar.baselines.constant_taxid import iter_fastq_ids
from samovar.profiles import metadata, validate


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", required=True)
    parser.add_argument("-I", default="")
    parser.add_argument("-o", required=True)
    parser.add_argument("-d", default="")
    parser.add_argument("-t", default=1)
    parser.add_argument("--taxid", default="562")
    args = parser.parse_args(argv)
    count = sum(1 for _ in iter_fastq_ids(args.i))
    frame = pd.DataFrame({"taxid": [args.taxid] if count else [], "N_sample": [count] if count else []})
    frame, meta = validate(frame, metadata())
    Path(args.o).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.o, index=False)
    Path(args.o + ".profile.json").write_text(json.dumps(meta))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
