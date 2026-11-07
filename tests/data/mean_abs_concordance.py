"""Concordance scorer that is not Mantel: mean absolute upper-triangle difference."""

from __future__ import annotations

import numpy as np
import pandas as pd


def score(left, right, config=None):
    _ = config
    left_values = pd.DataFrame(left).to_numpy(dtype=float)
    right_values = pd.DataFrame(right).to_numpy(dtype=float)
    index = np.triu_indices(left_values.shape[0], k=1)
    return {
        "method": "mean_abs",
        "version": "1",
        "statistic": float(np.mean(np.abs(left_values[index] - right_values[index]))),
        "pvalue": None,
        "n_samples": int(left_values.shape[0]),
        "seed": None,
        "permutations": None,
        "status": None,
        "reason": "",
    }
