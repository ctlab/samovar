"""Built-in sample filter: keep every sample that has a quality row."""

from __future__ import annotations

import pandas as pd

from samovar.abundance import n_sample_columns


SUPPORTED_PARAMS = ("n", "frac", "sd")


def filter_samples(table, scores, config=None):
    _ = config
    frame = table if isinstance(table, pd.DataFrame) else pd.read_csv(table)
    cols = n_sample_columns(frame)
    if not cols:
        return frame
    keep = cols
    if scores is not None and hasattr(scores, "columns") and "sample" in getattr(scores, "columns", []):
        named = {str(s) for s in scores["sample"].astype(str)}
        matched = [c for c in cols if c[2:] in named or c in named]
        if matched:
            keep = matched
    out_cols = [c for c in frame.columns if c not in cols or c in keep]
    return frame.loc[:, out_cols]
