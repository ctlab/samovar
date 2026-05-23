import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import joblib
import os
import re

def standardize_taxid_columns(df):
    df = df.copy()
    taxid_cols = [col for col in df.columns if col.lower().startswith('taxid_')]
    new_names = {}
    for col in taxid_cols:
        match = re.search(r'taxid_([^_]+)', col.lower())
        if match:
            tool_name = match.group(1)
            new_names[col] = f'taxid_{tool_name}'
            
    df = df.rename(columns=new_names)
    df = df.loc[:, ~df.columns.duplicated()]
    return df

def preprocess_data(df):
    df = df.copy()
    df = standardize_taxid_columns(df)
    
    # Safely drop non-feature columns if they exist
    cols_to_drop = ['seq', 'sample']
    for col in cols_to_drop:
        if col in df.columns:
            df.drop(col, axis=1, inplace=True)
            
    # Safely convert taxid columns to numeric
    for col_name in df.columns:
        if col_name.lower().startswith('taxid_'):
            series = df[col_name].fillna('0').astype(str)
            extracted = series.str.extract(r'(\d+)', expand=False).fillna('0')
            df[col_name] = pd.to_numeric(extracted, errors='coerce').fillna(0)
            
    if 'length' in df.columns:
        df['length'] = pd.to_numeric(df['length'], errors='coerce').fillna(0).astype(int)
        
    if 'true' in df.columns:
        df['true'] = pd.to_numeric(df['true'], errors='coerce')
        df = df.dropna(subset=['true'])
        df['true'] = df['true'].astype(int)
        
    return df

def train_models(df, test_size=0.2, random_state=42):
    # df already contains all features (merged in ML.py)
    df_processed = preprocess_data(df)
    
    # All columns except 'true' are our n-dimensional features (TaxIDs, GC, entropy, etc.)
    feature_cols = sorted([col for col in df_processed.columns if col != 'true'])
    print(f"[INFO] Initializing training with feature set: {feature_cols}")
    
    X = df_processed[feature_cols]
    y = df_processed['true']
    
    if len(y.unique()) > 1:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state, stratify=y)
    else:
        print("[WARNING] Only 1 class found in validation dataset. Stratification disabled.")
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    
    models = {
        'RandomForest': RandomForestClassifier(random_state=random_state, n_estimators=100),
        'AdaBoost': AdaBoostClassifier(random_state=random_state, n_estimators=100)
    }
    
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
            
    best_model.feature_cols_ = feature_cols
    return best_model, models, metrics, X_test, y_test

def plot_roc_curves(models, X_test, y_test, output_dir='tests_outs'):
    plt.figure(figsize=(10, 6))
    unique_classes = np.unique(y_test)
    
    if len(unique_classes) <= 1:
        print(f"[WARNING] Cannot plot ROC curves: Only one class ({unique_classes}) is present in the test set.")
        plt.text(0.5, 0.5, 'ROC AUC requires >= 2 classes\n(Only 1 class found in validation set)',
                 horizontalalignment='center', verticalalignment='center', fontsize=12)
    else:
        for name, model in models.items():
            try:
                y_score = model.predict_proba(X_test)
                classes = model.classes_
                y_binary = label_binarize(y_test, classes=classes)
                
                if y_binary.shape[1] == 1:
                    y_binary = np.hstack((1 - y_binary, y_binary))
                    
                fpr, tpr, _ = roc_curve(y_binary.ravel(), y_score.ravel())
                roc_auc = auc(fpr, tpr)
                plt.plot(fpr, tpr, label=f'{name} (AUC = {roc_auc:.2f})')
            except Exception as e:
                print(f"[ERROR] Failed to plot ROC for {name}: {e}")
        
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for Different Models')
    plt.legend(loc="lower right")
    
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'roc_comparison.png'))
    plt.close()

def predict_taxid(df, model_path=None, feature_cols=None):
    # df already contains all features (merged in ML.py)
    df_processed = preprocess_data(df)
    
    if isinstance(model_path, str) and os.path.exists(model_path):
        model = joblib.load(model_path)
    elif model_path is not None:
        model = model_path
    else:
        model, _, _, X_test, y_test = train_models(df_processed)
        
    if hasattr(model, 'feature_cols_'):
        feature_cols = model.feature_cols_
    elif feature_cols is None:
        feature_cols = sorted([col for col in df_processed.columns if col != 'true'])
        
    missing_features = [f for f in feature_cols if f not in df_processed.columns]
    if missing_features:
        raise ValueError(f"Missing required features: {missing_features}")
        
    X = df_processed[feature_cols]
    predictions = model.predict(X)
    probabilities = model.predict_proba(X)
    
    full_predictions = pd.Series(0, index=df.index, dtype=object)
    full_confidences = pd.Series(0.0, index=df.index, dtype=float)
    
    full_predictions.loc[df_processed.index] = predictions
    full_confidences.loc[df_processed.index] = np.max(probabilities, axis=1)
    
    df['taxid_SAMOVAR'] = full_predictions
    df['taxid_SAMOVAR_confidence'] = full_confidences
    return df

def save_model(model, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    joblib.dump(model, path)