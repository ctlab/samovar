import unittest
import pandas as pd
import numpy as np
import os
from src.samovar.reprofiling import (
    preprocess_data,
    train_models,
    plot_roc_curves,
    predict_taxid,
    save_model,
    annotator_taxid_columns,
)

import pytest

class TestReprofiling(unittest.TestCase):
    def setUp(self):
        """Set up test data"""
        # Create more diverse sample test data with more samples for better ROC curves
        self.test_data = pd.DataFrame({
            'seq': ['ATCG', 'GCTA', 'TAGC', 'CGAT', 'TACG', 'GATC', 'CTAG', 'AGCT',
                   'TTTT', 'AAAA', 'CCCC', 'GGGG', 'ATAT', 'CGCG', 'TATA', 'GCGC'],
            'taxid_annotator1': ['1', '2', '1', '2', '1', '2', '1', '2',
                                '1', '2', '1', '2', '1', '2', '1', '2'],
            'taxid_annotator2': ['2', '1', '2', '1', '2', '1', '2', '1',
                                '2', '1', '2', '1', '2', '1', '2', '1'],
            'length': [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            'true': [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
        })
        
        # Create test output directory
        import tempfile
        self.output_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.output_dir.cleanup)

    def test_preprocess_data(self):
        """Test data preprocessing"""
        processed_df = preprocess_data(self.test_data.copy())
        
        # Check if seq column is removed
        self.assertNotIn('seq', processed_df.columns)
        
        # Check if taxID columns are converted to numeric
        taxid_cols = [col for col in processed_df.columns if col.startswith('taxid_')]
        for col in taxid_cols:
            self.assertTrue(pd.api.types.is_numeric_dtype(processed_df[col]))
        
        # Check if length is converted to int
        self.assertTrue(pd.api.types.is_integer_dtype(processed_df['length']))
        
        # Check if true is converted to int
        self.assertTrue(pd.api.types.is_integer_dtype(processed_df['true']))

    def test_preprocess_factorizes_mag_id(self):
        df = pd.DataFrame({
            'seq': ['a', 'b', 'c', 'd'],
            'taxid_assembly_0': ['562', '562', '9606', '9606'],
            'feat_assembly_0_MAG_ID': ['magA', 'magB', 'magA', 'magB'],
            'length': [4, 4, 4, 4],
            'true': [562, 562, 9606, 9606],
        })
        processed = preprocess_data(df)
        self.assertTrue(pd.api.types.is_integer_dtype(processed['feat_assembly_0_MAG_ID']))
        self.assertEqual(processed['feat_assembly_0_MAG_ID'].nunique(), 2)

    def test_train_models(self):
        """Test model training functionality"""
        processed_df = preprocess_data(self.test_data.copy())
        best_model, models, metrics, feature_cols = train_models(processed_df, test_size=0.25)  # Adjusted test size
        
        self.assertIn('RandomForest', models)
        self.assertNotIn('AdaBoost', models)
        self.assertIn('RandomForest', metrics)
        
        # Check if best model is selected
        self.assertIsNotNone(best_model)
        
        # Check if metrics are reasonable
        for score in metrics.values():
            self.assertGreaterEqual(score, 0)
            self.assertLessEqual(score, 1)
        
        # Check if feature columns are taxid votes, length, and optional feat_*
        self.assertIsInstance(feature_cols, list)
        self.assertIn('length', feature_cols)
        extra = [col for col in feature_cols if col not in {'length', 'read_type'}]
        self.assertTrue(
            all(col.startswith('taxid_') or col.startswith('feat_') for col in extra)
        )

    def test_predict_taxid(self):
        """Test prediction functionality"""
        # Train model and make predictions
        processed_df = preprocess_data(self.test_data.copy())
        best_model, _, _, feature_cols = train_models(processed_df, test_size=0.25)  # Adjusted test size
        result_df = predict_taxid(self.test_data.copy(), model_path=best_model, feature_cols=feature_cols)
        
        # Check if prediction columns are added
        self.assertIn('taxid_SAMOVAR', result_df.columns)
        self.assertIn('taxid_SAMOVAR_confidence', result_df.columns)
        
        # Check if predictions are within expected range
        self.assertTrue(all(result_df['taxid_SAMOVAR_confidence'] >= 0))
        self.assertTrue(all(result_df['taxid_SAMOVAR_confidence'] <= 1))
        
        # Check if predictions are in the expected classes
        self.assertTrue(all(pred in [1, 2] for pred in result_df['taxid_SAMOVAR']))

    def test_save_and_load_model(self):
        """Test model saving and loading"""
        processed_df = preprocess_data(self.test_data.copy())
        best_model, _, _, feature_cols = train_models(processed_df, test_size=0.25)  # Adjusted test size
        
        # Save model
        model_path = os.path.join(self.output_dir.name, 'test_model.joblib')
        save_model(best_model, model_path)
        
        # Load model and make predictions
        result_df = predict_taxid(self.test_data.copy(), model_path=model_path, feature_cols=feature_cols)
        
        # Check if predictions are made
        self.assertIn('taxid_SAMOVAR', result_df.columns)
        self.assertIn('taxid_SAMOVAR_confidence', result_df.columns)

    def test_predict_production_colnames_not_all_zero(self):
        """taxID_* production columns must yield non-zero SAMOVAR calls."""
        data = pd.DataFrame({
            'seq': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'],
            'taxID_kaiju_0': [9606, 9606, 562, 562, 9606, 562, 0, 562],
            'taxID_kraken2_1': [9606, 562, 562, 9606, 9606, 562, 0, 562],
            'length': [100] * 8,
            'true': [9606, 9606, 562, 562, 9606, 562, 9606, 562],
        })
        result = predict_taxid(data.copy(), model_path=None)
        self.assertFalse((result['taxid_SAMOVAR'] == 0).all())
        self.assertGreater((result['taxid_SAMOVAR'] != 0).sum(), 4)
        self.assertEqual(int(result.loc[0, 'taxid_SAMOVAR']), 9606)
        self.assertEqual(int(result.loc[2, 'taxid_SAMOVAR']), 562)

    def test_annotator_taxid_columns_ignore_samovar(self):
        df = pd.DataFrame(columns=[
            'taxID_kaiju_0', 'taxID_kraken2_1', 'taxid_SAMOVAR',
            'taxid_SAMOVAR_confidence', 'length', 'true'
        ])
        self.assertEqual(
            annotator_taxid_columns(df),
            ['taxID_kaiju_0', 'taxID_kraken2_1'],
        )

    def test_standardize_coalesces_mangled_columns(self):
        data = pd.DataFrame({
            'taxID_kaiju...2': [9606, None, 562],
            'taxID_kaiju...4': [None, 9606, None],
            'taxID_kraken2...3': [9606, 562, 562],
            'length': [10, 10, 10],
            'true': [9606, 9606, 562],
        })
        processed = preprocess_data(data)
        self.assertIn('taxid_kaiju', processed.columns)
        self.assertIn('taxid_kraken2', processed.columns)
        self.assertNotIn('taxid_kaiju...2', processed.columns)
        self.assertEqual(list(processed['taxid_kaiju']), [9606, 9606, 562])

    def test_standardize_keeps_kraken2_databases_apart(self):
        data = pd.DataFrame({
            'taxID_kraken2_0': [562, 9606, 0],
            'taxID_kraken2_1': [9606, 562, 562],
            'feat_assembly_0_MAG_ID': ['m1', 'm2', 'm1'],
            'feat_assembly_1_MAG_ID': ['a', 'b', 'a'],
            'length': [10, 10, 10],
            'true': [562, 9606, 562],
        })
        processed = preprocess_data(data)
        self.assertIn('taxid_kraken2_0', processed.columns)
        self.assertIn('taxid_kraken2_1', processed.columns)
        self.assertNotIn('taxid_kraken2', processed.columns)
        self.assertIn('feat_assembly_0_MAG_ID', processed.columns)
        self.assertIn('feat_assembly_1_MAG_ID', processed.columns)
        self.assertEqual(list(processed['taxid_kraken2_0']), [562, 9606, 0])
        self.assertEqual(list(processed['taxid_kraken2_1']), [9606, 562, 562])

    @pytest.mark.optional
    def test_roc_curves(self):
        """Test ROC curve generation"""
        processed_df = preprocess_data(self.test_data.copy())
        best_model, models, _, feature_cols = train_models(processed_df, test_size=0.25)

        X_test = processed_df[feature_cols]
        y_test = processed_df['true']

        plot_roc_curves(models, X_test, y_test, output_dir=self.output_dir.name)
        self.assertTrue(os.path.exists(os.path.join(self.output_dir.name, 'roc_comparison.png')))

    def test_train_models_with_singleton_class(self):
        df = self.test_data.copy()
        df.loc[0, "true"] = 99
        processed = preprocess_data(df)
        best, models, metrics, cols = train_models(processed, test_size=0.25)
        self.assertIsNotNone(best)
        self.assertTrue(metrics)

    def test_passthrough_reprofile_tables(self):
        from samovar.reprofiling import passthrough_reprofile_tables
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as raw:
            src = Path(raw) / "in"
            dst = Path(raw) / "out"
            src.mkdir()
            pd.DataFrame({"taxid_a": [1, 2], "true": [562, 562]}).to_csv(
                src / "sample.annotation.csv", index=False
            )
            n = passthrough_reprofile_tables(str(src), str(dst))
            self.assertEqual(n, 1)
            out = pd.read_csv(dst / "sample_reprofiled.csv")
            self.assertEqual(list(out["taxid_SAMOVAR"]), [562, 562])

    @pytest.mark.optional
    def test_fit_uses_prepare_cores(self):
        """RandomForest and AdaBoost training call Parallel with prepare's core count."""
        import sklearn.ensemble._forest as forest
        import sklearn.multiclass as multiclass

        from samovar.reprofiling import Parallel as stump_parallel
        from samovar.reprofiling import build_adaboost, fit_logistic

        calls = []

        def _spy(real):
            class Spy:
                def __init__(self, *args, **kwargs):
                    calls.append(kwargs.get("n_jobs"))
                    self._inner = real(*args, **kwargs)

                def __call__(self, iterable):
                    return self._inner(iterable)

            return Spy

        forest_parallel = forest.Parallel
        ovr_parallel = multiclass.Parallel
        forest.Parallel = _spy(forest_parallel)
        multiclass.Parallel = _spy(ovr_parallel)
        import samovar.reprofiling as reprofiling

        saved = reprofiling.Parallel
        reprofiling.Parallel = _spy(stump_parallel)
        try:
            processed = preprocess_data(self.test_data.copy())
            _best, models, _metrics, _cols = train_models(
                processed, test_size=0.25, methods="ensemble", n_jobs=3
            )
            self.assertEqual(models["RandomForest"].n_jobs, 3)
            self.assertEqual(models["AdaBoost"].estimator.n_jobs, 3)
            rng = np.random.default_rng(0)
            frame = pd.DataFrame(
                {
                    "taxid_a": rng.integers(1, 4, 40),
                    "taxid_b": rng.integers(1, 4, 40),
                    "length": 100,
                    "true": rng.integers(0, 3, 40),
                }
            )
            logistic = fit_logistic(
                frame[["taxid_a", "taxid_b", "length"]],
                frame["true"],
                n_jobs=3,
                max_iter=50,
                solver="lbfgs",
            )
            self.assertEqual(logistic.n_jobs, 3)
            plain = build_adaboost(n_jobs=1, random_state=0, n_estimators=8)
            parallel = build_adaboost(n_jobs=3, random_state=0, n_estimators=8)
            X = processed[["taxid_annotator1", "taxid_annotator2", "length"]]
            y = processed["true"]
            plain.fit(X, y)
            parallel.fit(X, y)
            self.assertTrue(np.array_equal(plain.predict(X), parallel.predict(X)))
        finally:
            forest.Parallel = forest_parallel
            multiclass.Parallel = ovr_parallel
            reprofiling.Parallel = saved
        self.assertIn(3, calls)
        self.assertGreaterEqual(calls.count(3), 2)

    @pytest.mark.optional
    def test_gpu_mode_checks_each_model(self):
        import sys
        import types

        from sklearn.ensemble import RandomForestClassifier as SkForest

        import samovar.reprofiling as reprofiling
        from samovar.reprofiling import (
            GpuNotSupported,
            build_adaboost,
            build_random_forest,
            gpu_training_support,
            trains_on_gpu,
        )

        cpu = build_random_forest(n_estimators=2, n_jobs=1)
        self.assertIsInstance(cpu, SkForest)
        self.assertFalse(trains_on_gpu(cpu))
        with self.assertRaises(GpuNotSupported) as ada:
            build_adaboost(gpu=True, n_estimators=2)
        self.assertIn("CPU", str(ada.exception))

        class CudaForest:
            def __init__(self, **kwargs):
                self.kwargs = kwargs

        CudaForest.__module__ = "cuml.ensemble"

        class CpuLookalike:
            def __init__(self, **kwargs):
                self.kwargs = kwargs

        CpuLookalike.__module__ = "sklearn.ensemble"

        saved_load = reprofiling._load_gpu_class
        saved_cuda = reprofiling.cuda_device_available
        added = []
        for mod_name in ("cuml", "cuml.ensemble", "cuml.linear_model"):
            if mod_name not in sys.modules:
                sys.modules[mod_name] = types.ModuleType(mod_name)
                added.append(mod_name)
        reprofiling._load_gpu_class = lambda module, name: CudaForest
        reprofiling.cuda_device_available = lambda: True
        try:
            self.assertIsNone(gpu_training_support("RandomForest"))
            self.assertIsNotNone(gpu_training_support("AdaBoost"))
            model = build_random_forest(gpu=True, n_estimators=4, n_jobs=2)
            self.assertTrue(trains_on_gpu(model))
            self.assertEqual(model.kwargs["n_estimators"], 4)
            reprofiling._load_gpu_class = lambda module, name: CpuLookalike
            with self.assertRaises(GpuNotSupported) as rejected:
                build_random_forest(gpu=True, n_estimators=2)
            self.assertIn("does not train on GPU", str(rejected.exception))
        finally:
            reprofiling._load_gpu_class = saved_load
            reprofiling.cuda_device_available = saved_cuda
            for mod_name in added:
                sys.modules.pop(mod_name, None)


if __name__ == '__main__':
    unittest.main() 