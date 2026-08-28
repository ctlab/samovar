"""Test-only scoring wrapper (same contract as ``samovar tools import --type scoring``)."""

from pathlib import Path


def score(inputs, output_dir, config):
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
