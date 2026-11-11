"""Joint sample-profile correction through the existing export contract.

Each taxon learns nonnegative weights across profilers, without an intercept:
absent taxa cannot be invented. Model inputs remain typed sample profiles.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from samovar.abundance import n_sample_columns
from samovar.profiles import at_rank, compatible, fractions, metadata, write_profiles

POLICIES = dict(missing_taxa="zero", unseen_taxa="mean", profilers="exact",
                assigned_mass="mean observed", intercept=False)


def _matrix(tables, taxa, samples):
    # taxon × profiler × sample, with identities ordered before NumPy loses labels.
    return np.stack([table.reindex(taxa, fill_value=0)[samples].to_numpy(float)
                     for table in tables.values()], axis=1)


def _inputs(annotation, config, expected=None):
    """Align sample/feature order explicitly; read-derived tables stay separate."""
    infos = getattr(annotation, "profile_metadata", {})
    names = sorted(infos)
    if expected is not None and set(names) != set(expected):
        raise ValueError("Joint profile identities differ from training; retrain for missing/extra profilers")
    if len(names) < 2:
        raise ValueError("Joint correction requires at least two independent profiles")
    tables, metas, samples = {}, {}, None
    for name in names:
        info = metadata(infos[name])
        if expected and info["rank"] not in {expected[name]["rank"], "taxid"}:
            raise ValueError(f"Joint profile rank differs from training: {name}")
        rank = expected[name]["rank"] if expected else (info["rank"] if info["rank"] != "taxid" else "species")
        table, info = at_rank(annotation.abundance_tables[name], info, rank, config.get("taxdump"))
        if info["kind"] != "abundance":
            raise ValueError("Joint correction requires quantitative abundance profiles")
        if expected is not None:
            saved_info = metadata(expected[name])
            if saved_info["kind"] != "abundance":
                raise ValueError("Incompatible joint profile model metadata; retrain")
            compatible(info, saved_info)
            if info.get("database") != expected[name].get("database"):
                raise ValueError(f"Joint profile database differs from training: {name}")
        if metas:
            compatible(info, next(iter(metas.values())))
        cols = sorted(n_sample_columns(table))
        if samples is not None and cols != samples:
            raise ValueError("Joint profile sample identities differ")
        samples = cols
        values = fractions(table, info).set_index("taxid").drop(index="0", errors="ignore")
        if any(str(t).startswith("unmapped:") for t in values.index):
            raise ValueError("Joint correction needs resolved profile taxa")
        tables[name], metas[name] = values[cols], info
    return tables, metas, samples


def export(annotation, dest, config):
    """Fit on regenerated reference/GT, or apply a saved versioned JSON model."""
    from samovar.profile_pipeline import align_truth
    from scipy.optimize import nnls

    cfg = dict(config or {})
    dest = Path(dest)
    path = Path(cfg.get("profile_model") or dest / "profile_joint.json")
    reference = cfg.get("reference")
    if reference is not None:
        from samovar.abundance_correctors import as_annotation
        reference = as_annotation(reference)
        if cfg.get("reference_truth") is None:
            raise ValueError("Joint correction needs regenerated reference_truth")
        training, metas, samples = _inputs(reference, cfg)
        if len(samples) < 2:
            raise ValueError("Joint correction needs at least two regenerated samples")
        first = next(iter(training))
        truth, truth_meta = cfg["reference_truth"]
        # Reuse taxonomy, reference-length and normalization conversions from
        # profile_ratio/scoring. No initial truth participates in fitting.
        _, pm, truth, tm = align_truth(
            training[first].reset_index(), {**metas[first], "unit": "fraction"},
            truth, truth_meta, cfg.get("run_root", "."), cfg)
        compatible(pm, tm)
        if tm["kind"] != "abundance":
            raise ValueError("Presence truth cannot train joint quantitative correction")
        if set(n_sample_columns(truth)) != set(samples):
            raise ValueError("Joint training profile/truth sample identities differ")
        gold = fractions(truth, tm).set_index("taxid").drop(index="0", errors="ignore")
        if any(str(t).startswith("unmapped:") for t in gold.index):
            raise ValueError("Joint correction needs resolved truth taxa")
        names = list(training)
        taxa = sorted(set().union(*(table.index for table in training.values())))
        values = _matrix(training, taxa, samples)
        targets = gold.reindex(taxa, fill_value=0)[samples].to_numpy(float)
        weights = {}
        for i, tax in enumerate(taxa):
            x = values[i].T
            if not np.any(x):
                continue  # No training signal; use the documented unseen-taxon policy.
            weights[tax] = nnls(x, targets[i])[0].tolist()
        model = dict(version=1, method="profile_joint", estimator="nonnegative_least_squares",
                     profiles=metas, weights=weights, training_samples=samples,
                     unrecoverable_taxa=sorted(set(gold.index[gold[samples].sum(axis=1) > 0]) - set(weights)),
                     features=[[name, tax] for tax in weights for name in names],
                     policies=POLICIES)
    elif path.is_file():
        model = json.loads(path.read_text())
    else:
        raise ValueError("Joint correction needs a saved model or regenerated reference")
    if not isinstance(model, dict) or model.get("version") != 1 or model.get("method") != "profile_joint":
        raise ValueError("Incompatible joint profile model; retrain")
    if not isinstance(model.get("profiles"), dict) or not isinstance(model.get("weights"), dict):
        raise ValueError("Incompatible joint profile model schema; retrain")
    names = sorted(model["profiles"])
    expected_features = [[name, tax] for tax in sorted(model["weights"]) for name in names]
    if (model.get("features") != expected_features or model.get("estimator") != "nonnegative_least_squares"
            or model.get("policies") != POLICIES):
        raise ValueError("Incompatible joint profile feature schema; retrain")
    for weight in model["weights"].values():
        if len(weight) != len(names) or not np.isfinite(weight).all() or (np.asarray(weight) < 0).any():
            raise ValueError("Invalid joint profile coefficients")
    tables, metas, samples = _inputs(annotation, cfg, model["profiles"])
    taxa = sorted(set().union(*(table.index for table in tables.values())))
    values = _matrix(tables, taxa, samples)
    means = pd.DataFrame(values.mean(axis=1), index=taxa, columns=samples)
    weights = np.asarray([model["weights"].get(tax, [1 / len(names)] * len(names)) for tax in taxa]).reshape(len(taxa), len(names))
    result = pd.DataFrame(np.einsum("tps,tp->ts", values, weights), index=taxa, columns=samples)
    if not np.isfinite(result.to_numpy()).all():
        raise ValueError("Joint profile predictions are not finite")
    # All-normalized inputs retain their unknown mass. Empty predictions stay
    # empty; neither an intercept nor an imputed taxon can manufacture a hit.
    result = result.div(result.sum().replace(0, np.nan)).mul(means.sum()).fillna(0)
    info = {**metas[names[0]], "unit": "fraction", "identity": "profile_joint",
            "method": "profile_joint", "source_profiles": metas}
    info.pop("database", None)
    info.pop("sample", None)
    info.pop("source", None)
    info["unassigned_by_sample"] = {col[2:]: float(max(0, 1 - result[col].sum())) for col in samples}
    write_profiles(dest, {"profile_joint": result.rename_axis("taxid").reset_index()}, {"profile_joint": info})
    if reference is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(model, indent=2, allow_nan=False) + "\n")
    return dest
