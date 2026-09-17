"""Feature-importance scoring at reprofiler training time."""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

from samovar.feature_importance import (
    extract_native_importance,
    importance_table,
    maybe_score_reprofiler,
    score_feature_importance,
    selected_models,
    unwrap_estimator,
)
from samovar.reprofilers import ReprofileResult, run_reprofiler
from samovar.stage_report import REPORT_STAGES, bundle_multiqc


def _tiny():
    annotation = pd.DataFrame(
        {
            "seq": [f"r{i}" for i in range(16)],
            "taxid_dummy": [9606, 9606, 562, 562] * 4,
            "length": [50] * 16,
            "true": [9606, 9606, 562, 562] * 4,
        }
    )
    initial = {"obs": pd.DataFrame({"taxid": [9606, 562], "N_1": [4, 12]})}
    regenerated = {"obs": pd.DataFrame({"taxid": [9606, 562], "N_1": [8, 8]})}
    return annotation, initial, regenerated


class _FakeHMM:
    emissionprob_ = np.array([[0.1, 0.8, 0.1], [0.7, 0.05, 0.25]])


def test_extract_tree_and_linear_and_hmm():
    annotation, _, _ = _tiny()
    X = annotation[["taxid_dummy", "length"]]
    y = annotation["true"]
    rf = RandomForestClassifier(n_estimators=8, random_state=0).fit(X, y)
    vec, src = extract_native_importance(rf, list(X.columns))
    assert vec is not None and vec.size == 2
    assert "feature_importances_" in src
    lin = LogisticRegression(max_iter=200, random_state=0).fit(X, y)
    coef, src = extract_native_importance(lin, list(X.columns))
    assert coef is not None and coef.size == 2
    assert src == "coef_"
    hmm, src = extract_native_importance(_FakeHMM(), ["a", "b", "c"])
    assert hmm is not None and hmm.size == 3
    assert src == "emissionprob_"
    assert unwrap_estimator(rf) is rf


def test_score_writes_multiqc_json(tmp_path):
    annotation, initial, regenerated = _tiny()
    X = annotation[["taxid_dummy", "length"]]
    y = annotation["true"]
    model = RandomForestClassifier(n_estimators=8, random_state=0).fit(X, y)
    dest = tmp_path / "feature_importance_plots"
    payload = score_feature_importance(
        model,
        annotation,
        initial,
        regenerated,
        {"plot_dir": str(dest), "output_dir": str(tmp_path / "reprofiled_annotations"), "seed": 0},
    )
    assert dest.glob("*_mqc.json")
    bars = next(dest.glob("*_bars_mqc.json"))
    text = bars.read_text()
    assert '"plot_type": "bargraph"' in text
    assert "samovar_feature_importance" in text
    assert payload["n_features"] >= 1
    table = importance_table(
        model=model,
        annotation=annotation,
        initial_abundance=initial,
        regenerated_abundance=regenerated,
        feature_names=["taxid_dummy", "length"],
    )
    assert set(table["feature"]) == {"taxid_dummy", "length"}
    assert table["importance"].sum() > 0


def test_default_best_model_only():
    a = object()
    b = object()
    result = ReprofileResult(tables={}, model=a, models={"RandomForest": a, "AdaBoost": b})
    chosen = selected_models(result, all_models=False)
    assert chosen == [("RandomForest", a)]
    all_models = selected_models(result, all_models=True)
    assert {name for name, _ in all_models} == {"RandomForest", "AdaBoost"}


def test_run_reprofiler_writes_importance(tmp_path):
    annotation, initial_ab, regenerated_ab = _tiny()
    initial = {
        "sample.annotation": pd.DataFrame(
            {
                "seq": ["a", "b", "c"],
                "taxid_dummy": [9606, 562, 0],
                "length": [50, 50, 50],
            }
        )
    }
    dest = tmp_path / "reprofiled_annotations"
    run_reprofiler(
        "rf",
        regenerated=annotation,
        ground_truth=regenerated_ab,
        initial=initial,
        output_dir=dest,
        config={
            "seed": 0,
            "initial_abundance": initial_ab,
            "feature_importance": "builtin",
        },
    )
    plots = tmp_path / "feature_importance_plots"
    assert list(plots.glob("*_mqc.json"))
    staged = bundle_multiqc(tmp_path)
    names = [p.name for p in staged.glob("*_mqc.json")]
    assert any("feature_importance" in n for n in names)
    cfg = (staged / "multiqc_config.yaml").read_text()
    fi = cfg.find("samovar_feature_importance")
    repro = cfg.find("samovar_reprofiled")
    assert fi != -1 and repro != -1 and fi < repro
    assert REPORT_STAGES.index("feature_importance") < REPORT_STAGES.index("viz_reprofiled")


def test_skip_when_none(tmp_path):
    annotation, initial_ab, regenerated_ab = _tiny()
    result = ReprofileResult(tables={}, model=object())
    payloads = maybe_score_reprofiler(
        result,
        annotation=annotation,
        initial_abundance=initial_ab,
        regenerated_abundance=regenerated_ab,
        config={"feature_importance": "none", "plot_dir": str(tmp_path)},
    )
    assert payloads == []
