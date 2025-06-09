import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, auc
from sklearn.multiclass import OneVsRestClassifier
import matplotlib.pyplot as plt
import joblib
import os
import re
from collections import Counter

def standardize_taxid_columns(df):
    """
    Standardize taxID column names to format taxid_[toolname].
    Example: taxID_kaiju_1 -> taxid_kaiju
    
    Args:
        df (pd.DataFrame): Input dataframe
        
    Returns:
        pd.DataFrame: Dataframe with standardized column names
    """
    df = df.copy()
    # Find all taxID columns
    taxid_cols = [col for col in df.columns if col.lower().startswith('taxid_')]
    
    # Create mapping for new column names
    new_names = {}
    for col in taxid_cols:
        # Extract tool name using regex
        match = re.search(r'taxid_([^_]+)', col.lower())
        if match:
            tool_name = match.group(1)
            new_name = f'taxid_{tool_name}'
            new_names[col] = new_name
    
    # Rename columns
    df = df.rename(columns=new_names)
    return df

def preprocess_data(df):
    """
    Preprocess the input dataframe for ML training.
    
    Args:
        df (pd.DataFrame): Input dataframe with columns seq, taxID_*, length, true
        
    Returns:
        pd.DataFrame: Preprocessed dataframe ready for ML
    """
    # Create a copy to avoid modifying the original
    df = df.copy()
    
    # Standardize taxID column names
    df = standardize_taxid_columns(df)
    
    # Remove seq column as it's not needed for prediction
    if 'seq' in df.columns:
        df.drop('seq', axis=1, inplace=True)
    
    # Remove sample column if present
    if 'sample' in df.columns:
        df.drop('sample', axis=1, inplace=True)
    
    # Convert taxID columns to categorical and handle NA values
    taxid_cols = [col for col in df.columns if col.startswith('taxid_')]
    for col in taxid_cols:
        # First convert to string and handle NA values
        df[col] = df[col].fillna('0').astype(str)
        # Remove any non-numeric characters and convert to numeric
        df[col] = df[col].str.extract(r'(\d+)', expand=False).fillna('0')
        # Convert to numeric
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    # Convert length to int if present
    if 'length' in df.columns:
        df['length'] = pd.to_numeric(df['length'], errors='coerce').fillna(0).astype(int)
    
    # Convert true to numeric and handle NaN values
    if 'true' in df.columns:
        # First convert to numeric
        df['true'] = pd.to_numeric(df['true'], errors='coerce')
        # Remove rows with NaN in true column
        df = df.dropna(subset=['true'])
        # Convert to int after removing NaN
        df['true'] = df['true'].astype(int)
    
    return df

def train_models(df, test_size=0.2, random_state=42):
    """
    Train RandomForest and AdaBoost models for taxID prediction.
    
    Args:
        df (pd.DataFrame): Input dataframe
        test_size (float): Proportion of data to use for testing
        random_state (int): Random seed for reproducibility
        
    Returns:
        tuple: (best_model, models_dict, metrics_dict, feature_cols)
    """
    # Preprocess the data first
    df_processed = preprocess_data(df)
    
    # Select features for training (exclude seq and true)
    feature_cols = sorted([col for col in df_processed.columns 
                         if col not in ['seq', 'true']])
    
    # Ensure all required features are present
    required_features = ['length'] + [col for col in feature_cols if col.startswith('taxid_')]
    missing_features = [f for f in required_features if f not in feature_cols]
    if missing_features:
        raise ValueError(f"Missing required features: {missing_features}")
    
    # Prepare features and target
    X = df_processed[feature_cols]
    y = df_processed['true']
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )
    
    # Initialize models
    models = {
        'RandomForest': RandomForestClassifier(random_state=random_state, n_estimators=100),
        'AdaBoost': AdaBoostClassifier(random_state=random_state, n_estimators=100)
    }
    
    # Train and evaluate models
    metrics = {}
    best_score = -np.inf
    best_model = None
    
    for name, model in models.items():
        model.fit(X_train, y_train)
        score = model.score(X_test, y_test)
        metrics[name] = score
        
        if score > best_score:
            best_score = score
            best_model = model
    
    return best_model, models, metrics, feature_cols

def plot_roc_curves(models, X_test, y_test, output_dir='tests_outs'):
    """
    Plot ROC curves for all models using micro-averaging.
    
    Args:
        models (dict): Dictionary of trained models
        X_test (pd.DataFrame): Test features
        y_test (np.array): Test target
        output_dir (str): Directory to save the plot
    """
    plt.figure(figsize=(10, 6))
    
    for name, model in models.items():
        # Get prediction probabilities
        y_score = model.predict_proba(X_test)
        
        # Convert to binary classification using micro-averaging
        y_binary = pd.get_dummies(y_test).values
        fpr, tpr, _ = roc_curve(y_binary.ravel(), y_score.ravel())
        roc_auc = auc(fpr, tpr)
        
        # Plot ROC curve for this model
        plt.plot(fpr, tpr, 
                label=f'{name} (AUC = {roc_auc:.2f})')
    
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for Different Models')
    plt.legend(loc="lower right")
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'roc_comparison.png'))
    plt.close()

def predict_taxid(df, model_path=None, feature_cols=None):
    """
    Predict taxID and confidence for input dataframe.
    
    Args:
        df (pd.DataFrame): Input dataframe
        model_path (str or model): Path to saved model or model object
        feature_cols (list): List of feature columns in the correct order
        
    Returns:
        pd.DataFrame: Input dataframe with added prediction columns
    """
    # Create a copy to avoid modifying the original
    df = df.copy()
    
    # Preprocess data
    df_processed = preprocess_data(df)
    
    # Use provided feature columns or get them from the dataframe
    if feature_cols is None:
        feature_cols = sorted([col for col in df_processed.columns 
                             if col not in ['seq', 'true']])
    
    # Ensure all required features are present
    required_features = ['length'] + [col for col in feature_cols if col.startswith('taxid_')]
    missing_features = [f for f in required_features if f not in df_processed.columns]
    if missing_features:
        raise ValueError(f"Missing required features: {missing_features}")
    
    # Load or use provided model
    if isinstance(model_path, str) and os.path.exists(model_path):
        model = joblib.load(model_path)
    elif model_path is not None:
        model = model_path  # Use provided model object
    else:
        # Train new model if no model provided
        model, _, _, feature_cols = train_models(df_processed)
    
    # Make predictions using the same feature order as training
    X = df_processed[feature_cols]
    predictions = model.predict(X)
    probabilities = model.predict_proba(X)
    
    # Get the most common non-zero prediction for each row
    final_predictions = []
    final_confidences = []
    
    for i, pred in enumerate(predictions):
        # Get all non-zero taxID values for this row
        taxid_values = [val for val in X.iloc[i] if val > 0 and str(val).startswith('taxid_')]
        
        if not taxid_values:  # If no non-zero taxIDs
            final_predictions.append(pred)  # Use model prediction
            final_confidences.append(np.max(probabilities[i]))  # Use model confidence
        else:
            # Count occurrences of each non-zero taxID
            counts = Counter(taxid_values)
            most_common = counts.most_common(1)[0][0]
            
            final_predictions.append(most_common)
            final_confidences.append(0.0)  # Set confidence to 0 for these cases
    
    # Create a mapping from processed index to original index
    processed_to_original = dict(zip(range(len(df_processed)), df_processed.index))
    
    # Initialize arrays for the full length of original dataframe
    full_predictions = np.zeros(len(df), dtype=object)
    full_confidences = np.zeros(len(df))
    
    # Fill in predictions and confidences for processed rows
    for i, (pred, conf) in enumerate(zip(final_predictions, final_confidences)):
        orig_idx = processed_to_original[i]
        full_predictions[orig_idx] = pred
        full_confidences[orig_idx] = conf
    
    # Add predictions to dataframe
    df['taxid_SAMOVAR'] = full_predictions
    df['taxid_SAMOVAR_confidence'] = full_confidences
    
    return df

def save_model(model, path):
    """
    Save trained model to disk.
    
    Args:
        model: Trained model object
        path (str): Path to save model
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    joblib.dump(model, path)
