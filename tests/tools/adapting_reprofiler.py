"""Test-only ML wrapper that records when fitting runs (``--type ml``)."""

from pathlib import Path

from samovar.reprofilers import reprofile_linear


def reprofile(regenerated, ground_truth, initial, config):
    dest = Path((config or {}).get("output_dir") or ".")
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "adapted.txt").write_text("retrained\n", encoding="utf-8")
    return reprofile_linear(regenerated, ground_truth, initial, config or {})
