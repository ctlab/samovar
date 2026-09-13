"""Example ``samovar tools import --type sample-score`` plugin.

Returns 1.0 for every generated sample that has a finite count vector.
"""

from __future__ import annotations

import pandas as pd

from samovar.abundance import abundance_to_matrix, n_sample_columns, normalize_abundance_table


def score_samples(generated, reference, config=None):
    _ = reference, config
    frame = generated
    if not isinstance(frame, pd.DataFrame):
        frame = pd.read_csv(frame)
    if "taxid" in frame.columns or n_sample_columns(frame):
        mat = abundance_to_matrix(normalize_abundance_table(frame))
    else:
        mat = pd.DataFrame()
    rows = []
    for col in mat.columns:
        rows.append({"sample": str(col), "quality": 1.0, "ok": True})
    return pd.DataFrame(rows)
