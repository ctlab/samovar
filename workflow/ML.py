import argparse
import os
import pandas as pd
from samovar.parse_annotators import Annotation, match_annotation
from samovar.reprofiling import train_models, predict_taxid, save_model, plot_roc_curves

def process_sample(sample_file, output_dir, model=None, label_encoder=None):
    """Process a single sample with the trained model."""
    print(f"Processing {sample_file}...")
    df = pd.read_csv(sample_file)
    
    # Fill NaN values with 0 in all columns
    df = df.fillna(0)
    
    # Make predictions
    result_df = predict_taxid(df, model_path=None if model is None else model)
    
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
args = parser.parse_args()

# Create output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Load validation data
if not os.path.exists(args.validation_file):
    raise FileNotFoundError(f"Validation file not found: {args.validation_file}")

validation_df = pd.read_csv(args.validation_file)
validation_df.drop("seq", axis=1, inplace=True)

# Drop rows with NaN in true column
initial_rows = len(validation_df)
validation_df = validation_df.dropna(subset=['true'])
dropped_rows = initial_rows - len(validation_df)
print(f"\nDropped {dropped_rows} rows with NaN in 'true' column")
print(f"Remaining rows for training: {len(validation_df)}")

# Fill NaN values with 0 in all other columns
validation_df = validation_df.fillna(0)

# Train model
print("\nTraining model on validation data...")
best_model, models, metrics, feature_cols = train_models(validation_df)

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
validation_df.columns = [col.lower() for col in validation_df.columns]
X_test = validation_df[feature_cols]
y_test = validation_df['true']
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
            model=best_model)