import unittest
import pandas as pd
import numpy as np
import os
from src.samovar.reprofiling import (
    preprocess_data,
    train_models,
    plot_roc_curves,
    predict_taxid,
    save_model
)

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
        os.makedirs('tests_outs', exist_ok=True)

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

    def test_train_models(self):
        """Test model training functionality"""
        processed_df = preprocess_data(self.test_data.copy())
        best_model, models, metrics, feature_cols = train_models(processed_df, test_size=0.25)  # Adjusted test size
        
        # Check if we got both models
        self.assertIn('RandomForest', models)
        self.assertIn('AdaBoost', models)
        
        # Check if metrics are calculated
        self.assertIn('RandomForest', metrics)
        self.assertIn('AdaBoost', metrics)
        
        # Check if best model is selected
        self.assertIsNotNone(best_model)
        
        # Check if metrics are reasonable
        for score in metrics.values():
            self.assertGreaterEqual(score, 0)
            self.assertLessEqual(score, 1)
        
        # Check if feature columns are correct
        self.assertIsInstance(feature_cols, list)
        self.assertIn('length', feature_cols)
        self.assertTrue(all(col.startswith('taxid_') for col in feature_cols if col != 'length'))

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
        model_path = 'tests_outs/test_model.joblib'
        save_model(best_model, model_path)
        
        # Load model and make predictions
        result_df = predict_taxid(self.test_data.copy(), model_path=model_path, feature_cols=feature_cols)
        
        # Check if predictions are made
        self.assertIn('taxid_SAMOVAR', result_df.columns)
        self.assertIn('taxid_SAMOVAR_confidence', result_df.columns)

    def test_roc_curves(self):
        """Test ROC curve generation"""
        processed_df = preprocess_data(self.test_data.copy())
        best_model, models, _, feature_cols = train_models(processed_df, test_size=0.25)  # Adjusted test size
        
        # Get test data
        X_test = processed_df[feature_cols]
        y_test = processed_df['true']
        
        # Plot ROC curves
        plot_roc_curves(models, X_test, y_test)
        
        # Check if plot file exists
        self.assertTrue(os.path.exists('tests_outs/roc_comparison.png'))

if __name__ == '__main__':
    unittest.main() 