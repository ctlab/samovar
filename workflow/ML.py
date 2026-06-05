import argparse
import os

import pandas as pd

from samovar.parse_annotators import Annotation, match_annotation
from samovar.reprofiling import plot_roc_curves, predict_taxid, save_model, train_models


def process_sample(
    sample_file, output_dir, model=None, classify_unclassified=False, features_df=None
):
    print(f"Processing {sample_file}...")
    df = pd.read_csv(sample_file, sep="\t")
    df.columns = [col.lower() for col in df.columns]

    # [FIX]: Merge biological features into the sample dataframe before prediction
    if features_df is not None:
        df = pd.merge(df, features_df, on="seq", how="left")

    df = df.fillna(0)

    # Remove old samovar columns if they exist
    cols_to_drop = [c for c in df.columns if "samovar" in c.lower()]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)

    result_df = predict_taxid(df, model_path=None if model is None else model)

    if not classify_unclassified:
        taxid_cols = [col for col in df.columns if col.startswith("taxid_")]
        unclassified_mask = (df[taxid_cols] == 0).all(axis=1)
        result_df.loc[unclassified_mask, "taxid_SAMOVAR"] = 0
        result_df.loc[unclassified_mask, "taxid_SAMOVAR_confidence"] = 0

    output_file = os.path.join(
        output_dir, f"{os.path.basename(sample_file).split('.')[0]}_reprofiled.csv"
    )
    result_df.to_csv(output_file, index=False, sep="\t")
    print(f"Results saved to {output_file}")


parser = argparse.ArgumentParser()
parser.add_argument(
    "--reprofiling_dir",
    "-r",
    type=str,
    required=True,
    help="Directory containing files to be reprofiled",
)
parser.add_argument(
    "--validation_file",
    "-v",
    type=str,
    required=True,
    help="File containing validation data",
)
parser.add_argument(
    "--output_dir", "-o", type=str, required=True, help="Directory to save output files"
)
parser.add_argument(
    "--classify-unclassified",
    action="store_true",
    help="If set, will attempt to classify sequences that have all taxid fields = 0",
)
# [FIX]: Added features argument
parser.add_argument(
    "--features",
    "-f",
    type=str,
    required=False,
    help="TSV file containing extracted biological features",
)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

if not os.path.exists(args.validation_file):
    raise FileNotFoundError(f"Validation file not found: {args.validation_file}")

# Load validation data
validation_df = pd.read_csv(args.validation_file, sep="\t")
validation_df.columns = [col.lower() for col in validation_df.columns]

# Load and merge features
features_df = None
if args.features and os.path.exists(args.features):
    print(f"\n[INFO] Loading features from {args.features}...")
    features_df = pd.read_csv(args.features, sep="\t")

    # Rename 'read_id' from fastq_annotator to 'seq' to match validation_df
    if "read_id" in features_df.columns:
        features_df = features_df.rename(columns={"read_id": "seq"})

    initial_shape = validation_df.shape
    validation_df = pd.merge(validation_df, features_df, on="seq", how="left")
    print(
        f"[INFO] Merged features: columns increased from {initial_shape[1]} to {validation_df.shape[1]}"
    )

# Clean up training data
initial_rows = len(validation_df)
validation_df = validation_df.dropna(subset=["true"])
dropped_rows = initial_rows - len(validation_df)
print(f"\nDropped {dropped_rows} rows with NaN in 'true' column")
print(f"Remaining rows for training: {len(validation_df)}")

# [FIX]: Drop 'seq' AFTER merging features, because train_models doesn't need string IDs
if "seq" in validation_df.columns:
    validation_df = validation_df.drop("seq", axis=1)

validation_df = validation_df.fillna(0)

# Train model
print("\nTraining model on validation data...")
# Assuming train_models unpacks to 5 variables in your updated samovar.reprofiling
best_model, models, metrics, X_test_processed, y_test_processed = train_models(
    validation_df
)

print("\nModel performance:")
for name, score in metrics.items():
    print(f"{name}: {score:.3f}")

model_path = os.path.join(args.output_dir, "trained_model.joblib")
save_model(best_model, model_path)
print(f"\nBest model saved to {model_path}")

print("\nGenerating ROC curves...")
plot_roc_curves(models, X_test_processed, y_test_processed, output_dir=args.output_dir)
print(f"ROC curves saved to {args.output_dir}/roc_comparison.png")

# Process files for reprofiling (inference)
print("\nProcessing files for reprofiling...")
for filename in os.listdir(args.reprofiling_dir):
    if filename.endswith(".csv"):
        sample_file = os.path.join(args.reprofiling_dir, filename)
        process_sample(
            sample_file,
            args.output_dir,
            model=best_model,
            classify_unclassified=args.classify_unclassified,
            features_df=features_df,  # Pass features to be merged inside!
        )
