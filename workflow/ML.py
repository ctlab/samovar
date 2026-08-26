import argparse
import os
import pandas as pd
from samovar.parse_annotators import Annotation, match_annotation
from samovar.reprofiling import (
    train_models,
    predict_taxid,
    save_model,
    plot_roc_curves,
    annotator_taxid_columns,
    preprocess_data,
    merge_read_features,
    load_read_features,
    passthrough_reprofile_tables,
)

def process_sample(sample_file, output_dir, model=None, label_encoder=None, classify_unclassified=False, features_df=None):
    """Process a single sample with the trained model.
    
    Args:
        sample_file (str): Path to the sample file
        output_dir (str): Directory to save output
        model: Trained model to use for prediction
        label_encoder: Label encoder for the model
        classify_unclassified (bool): If False, sequences with all taxid fields = 0 will not be classified
        features_df: Optional per-read biological features to merge on seq
    """
    print(f"Processing {sample_file}...")
    df = pd.read_csv(sample_file)
    df = merge_read_features(df, features_df)

    # Fill missing predictor fields with 0, but never coerce missing true
    # taxIDs to 0 (0 means unclassified, not "unknown truth").
    fill_cols = [c for c in df.columns if c != "true"]
    df[fill_cols] = df[fill_cols].fillna(0)
    
    # Make predictions
    result_df = predict_taxid(df, model_path=None if model is None else model)
    
    # If classify_unclassified is False, set prediction to 0 for sequences with all taxid fields = 0
    if not classify_unclassified:
        taxid_cols = annotator_taxid_columns(df)
        if taxid_cols:
            numeric = df[taxid_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
            unclassified_mask = (numeric == 0).all(axis=1)
            result_df.loc[unclassified_mask, "taxid_SAMOVAR"] = 0
            result_df.loc[unclassified_mask, "taxid_SAMOVAR_confidence"] = 0
    
    # Save results
    output_file = os.path.join(output_dir, f"{os.path.basename(sample_file).split('.')[0]}_reprofiled.csv")
    result_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--reprofiling_dir", "-r", type=str, required=True,
                    help="Directory containing files to be reprofiled")
parser.add_argument("--validation_file", "-v", type=str, required=True,
                    help="File containing validation data")
parser.add_argument("--output_dir", "-o", type=str, required=True,
                    help="Directory to save output files")
parser.add_argument("--classify-unclassified", action="store_true",
                    help="If set, will attempt to classify sequences that have all taxid fields = 0")
parser.add_argument("--features", "-f", type=str, required=False,
                    help="TSV/CSV of per-read features from src/annotators/fastq_annotator.py")
parser.add_argument("--seed", type=int, default=42,
                    help="Random seed for train/test split and classifiers")
args = parser.parse_args()

# Create output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Load validation data
if not os.path.exists(args.validation_file):
    raise FileNotFoundError(f"Validation file not found: {args.validation_file}")

validation_df = pd.read_csv(args.validation_file)
features_df = load_read_features(getattr(args, "features", None))
if features_df is not None:
    print(f"\n[INFO] Loading features from {args.features}...")
    before = validation_df.shape[1]
    validation_df = merge_read_features(validation_df, features_df)
    print(f"[INFO] Merged features: columns increased from {before} to {validation_df.shape[1]}")
if "seq" in validation_df.columns:
    validation_df.drop("seq", axis=1, inplace=True)

# Drop rows with NaN in true column
initial_rows = len(validation_df)
validation_df = validation_df.dropna(subset=['true'])
dropped_rows = initial_rows - len(validation_df)
print(f"\nDropped {dropped_rows} rows with NaN in 'true' column")
print(f"Remaining rows for training: {len(validation_df)}")

# Fill NaN values with 0 in predictors only; do not coerce missing true taxIDs to 0.
fill_cols = [c for c in validation_df.columns if c != "true"]
validation_df[fill_cols] = validation_df[fill_cols].fillna(0)

# Train model
print("\nTraining model on validation data...")
if len(validation_df) == 0:
    print("No validation rows available; skipping ML reprofiling.")
    n = passthrough_reprofile_tables(
        args.reprofiling_dir, args.output_dir, features_df=features_df
    )
    print(f"Wrote {n} passthrough reprofile table(s) (no model).")
    raise SystemExit(0)

try:
    best_model, models, metrics, feature_cols = train_models(
        validation_df, random_state=int(args.seed)
    )
except ValueError as exc:
    print(f"Skipping ML reprofiling: {exc}")
    n = passthrough_reprofile_tables(
        args.reprofiling_dir, args.output_dir, features_df=features_df
    )
    print(f"Wrote {n} passthrough reprofile table(s) (taxid_SAMOVAR copies true).")
    raise SystemExit(0)

# Print model performance
print("\nModel performance:")
for name, score in metrics.items():
    print(f"{name}: {score:.3f}")

# Save best model
model_path = os.path.join(args.output_dir, "trained_model.joblib")
save_model(best_model, model_path)
print(f"\nBest model saved to {model_path}")

# Plot ROC curves
print("\nGenerating ROC curves...")
processed_validation = preprocess_data(validation_df.copy())
X_test = processed_validation.reindex(columns=feature_cols, fill_value=0)
y_test = processed_validation['true']
plot_roc_curves(models, X_test, y_test, output_dir=args.output_dir)
print(f"ROC curves saved to {args.output_dir}/roc_comparison.png")

# Process each file in the reprofiling directory
print("\nProcessing files for reprofiling...")
for filename in os.listdir(args.reprofiling_dir):
    if filename.endswith('.csv'):
        sample_file = os.path.join(args.reprofiling_dir, filename)
        process_sample(
            sample_file, 
            args.output_dir, 
            model=best_model,
            classify_unclassified=args.classify_unclassified,
            features_df=features_df)