#!/usr/bin/env python3
"""Custom scoring/viz wrapper: count annotation files under matched inputs.

Import::

    samovar tools import -n counts \\
        --exec-path examples/scoring/count_annotations.py --type scoring \\
        --inputs '*annotations' --flags "--min-files 0"
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Sequence


def score(inputs: Sequence[Path], output_dir: Path, config: dict) -> None:
    extra = list(config.get("extra_argv") or [])
    min_files = 0
    if "--min-files" in extra:
        idx = extra.index("--min-files")
        if idx + 1 < len(extra):
            min_files = int(extra[idx + 1])
    dest = Path(output_dir)
    dest.mkdir(parents=True, exist_ok=True)
    lines = ["path\tn_files"]
    for raw in inputs:
        path = Path(raw)
        if path.is_dir():
            n = sum(1 for p in path.iterdir() if p.is_file())
        elif path.is_file():
            n = 1
        else:
            n = 0
        if n < min_files:
            continue
        lines.append(f"{path}\t{n}")
    (dest / "file_counts.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    stage = config.get("stage") or ""
    if stage:
        (dest / "stage.txt").write_text(str(stage) + "\n", encoding="utf-8")


def _parse_cli(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Count files in scoring inputs")
    parser.add_argument("-i", action="append", dest="inputs", default=[])
    parser.add_argument("-o", dest="output_dir", required=True)
    parser.add_argument("--min-files", type=int, default=0)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = _parse_cli()
    score(
        [Path(p) for p in args.inputs],
        Path(args.output_dir),
        {"extra_argv": ["--min-files", str(args.min_files)]},
    )
