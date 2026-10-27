import pandas as pd
import numpy as np
import importlib
import inspect
from contextlib import contextmanager
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, auc
from sklearn.multiclass import OneVsRestClassifier
from sklearn.tree import DecisionTreeClassifier
import matplotlib.pyplot as plt
import joblib
import os
import re
import sys
from typing import Any, Optional
from collections import Counter

from samovar.memory import parallel_jobs, predict_threads, proba_chunk_rows, read_csv_memmap

def _tool_name_from_column(col: str) -> Optional[str]:
    """Extract annotator name from taxid/taxID column labels, never from feat_."""
    name = str(col)
    low = name.lower()
    if low.startswith("feat_") or "confidence" in low:
        return None
    match = re.match(r"taxid[_./]*([a-z][a-z0-9]*)", low)
    if not match:
        return None
    tool = match.group(1)
    if tool in {"samovar", "true", "length", "sample", "read", "feat"}:
        return None
    return tool


def merge_read_features(df: pd.DataFrame, features_df: pd.DataFrame) -> pd.DataFrame:
    """Join per-read features (fastq_annotator) onto an annotation table."""
    if features_df is None or features_df.empty:
        return df
    feats = features_df.copy()
    if "read_id" in feats.columns and "seq" not in feats.columns:
        feats = feats.rename(columns={"read_id": "seq"})
    if "seq" not in df.columns or "seq" not in feats.columns:
        return df
    return pd.merge(df, feats, on="seq", how="left")


def load_read_features(path: Optional[str]) -> Optional[pd.DataFrame]:
    """Load a TSV/CSV of per-read features, or None if missing."""
    if not path or not os.path.exists(path):
        return None
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    if sep != ",":
        return pd.read_csv(path, sep=sep)
    from samovar.memory import memmap_min_bytes, read_csv_memmap

    if os.path.getsize(path) < memmap_min_bytes():
        return pd.read_csv(path, sep=sep)
    return read_csv_memmap(path)


def annotator_taxid_columns(df: pd.DataFrame) -> list:
    """Return annotator taxid columns, ignoring SAMOVAR outputs and confidence."""
    return [col for col in df.columns if _tool_name_from_column(col) is not None]


def standardize_taxid_columns(df):
    """
    Standardize taxID column names to format taxid_[toolname].
    Example: taxID_kaiju_1 -> taxid_kaiju
    Duplicate columns for the same tool (taxID_kaiju...2) are coalesced.
    """
    df = df.copy()
    grouped = {}
    for col in list(df.columns):
        tool = _tool_name_from_column(col)
        if tool is None:
            continue
        grouped.setdefault(tool, []).append(col)

    for tool, cols in grouped.items():
        new_name = f"taxid_{tool}"
        numeric = df[cols].apply(pd.to_numeric, errors="coerce")
        coalesced = numeric.bfill(axis=1).iloc[:, 0]
        drop_cols = [c for c in cols if c != new_name]
        if drop_cols:
            df = df.drop(columns=drop_cols)
        df[new_name] = coalesced
    return df

class IncompatiblePreprocessingError(ValueError):
    """A saved model cannot safely encode the supplied features."""


def _category_values(values, column):
    values = values.astype("string").fillna("").astype(str)
    if column == "MAG_ID" or column.endswith("_MAG_ID"):
        # Contract alignment can fill an absent MAG feature with numeric zero.
        values = values.replace({"0": "", "0.0": ""})
    return values


def preprocess_data(df, drop_missing_true: bool = True, *, preprocessing=None):
    """
    Preprocess the input dataframe for ML training.

    Args:
        df (pd.DataFrame): Input dataframe with columns seq, taxID_*, length, true
        drop_missing_true: Drop rows with NaN in ``true`` (training only; inference
            on initial annotations keeps rows when ``true`` is unknown).
        preprocessing: Fitted schema from ``model.samovar_preprocessing_``.
            Without it, fit dictionaries and retain them in dataframe attrs for
            ``train_models(preprocess_data(raw))`` callers. Missing and unseen
            categories use -1; dictionaries are never extended during prediction.

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
    
    # Convert length to int if present; dummy/custom tables often omit it
    if 'length' not in df.columns:
        df['length'] = 0
    df['length'] = pd.to_numeric(df['length'], errors='coerce').fillna(0).astype(int)

    encoded = df.attrs.get("samovar_preprocessing")
    if encoded is not None and preprocessing is not None and encoded != preprocessing:
        raise IncompatiblePreprocessingError("Features were encoded with a different model's dictionaries")
    schema = preprocessing if preprocessing is not None else encoded
    if schema is None:
        features = {}
        for col in df.columns:
            if col == "true":
                continue
            numeric = pd.to_numeric(df[col], errors="coerce")
            categorical = col == "read_type" or col == "MAG_ID" or col.endswith("_MAG_ID")
            if not categorical and (
                numeric.notna().all()
                or numeric.notna().sum() >= max(1, int(0.9 * len(df)))
            ):
                features[col] = None
            else:
                values = _category_values(df[col], col)
                features[col] = {value: i for i, value in enumerate(sorted(set(values) - {""}))}
        schema = {"version": 1, "features": features}
    if schema.get("version") != 1:
        raise IncompatiblePreprocessingError("Unsupported model preprocessing version; retrain the reprofiler")
    for col, categories in schema["features"].items():
        if col not in df.columns:
            df[col] = 0 if categories is None else -1
        elif encoded is None:
            if categories is None:
                df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
            else:
                df[col] = _category_values(df[col], col).map(categories).fillna(-1).astype(int)
    df = df[[*schema["features"], *(["true"] if "true" in df.columns else [])]]
    df.attrs["samovar_preprocessing"] = schema
    
    # Convert true to numeric and handle NaN values
    if 'true' in df.columns:
        df['true'] = pd.to_numeric(df['true'], errors='coerce')
        if drop_missing_true:
            df = df.dropna(subset=['true'])
            df['true'] = df['true'].astype(int)
    
    return df


def preprocess_for_model(df, model, *, drop_missing_true=True):
    """Transform with saved dictionaries; legacy taxid/length models remain usable."""
    schema = getattr(model, "samovar_preprocessing_", None)
    names = getattr(model, "feature_names_in_", None)
    if schema is None:
        df = standardize_taxid_columns(df.copy())
        if names is not None:
            df = df.reindex(columns=[*names, *(["true"] if "true" in df.columns else [])])
    processed = preprocess_data(df, drop_missing_true, preprocessing=schema)
    if schema is None and any(
        col != "length" and not col.startswith("taxid_")
        for col in processed.attrs["samovar_preprocessing"]["features"]
    ):
        # A numeric-only CSV chunk cannot prove a legacy feat_* was numeric at
        # training time: it might be a subset of string categories like "001".
        raise IncompatiblePreprocessingError(
            "Legacy model has no preprocessing schema for auxiliary features. "
            "Retrain the reprofiler from regenerated annotations before prediction."
        )
    if names is not None:
        processed = processed.reindex(columns=[*names, *(["true"] if "true" in processed.columns else [])])
    return processed


def training_n_jobs(n_jobs=None, config=None) -> int:
    """Thread count from ``samovar prepare --cores`` (never ``-1`` / all CPUs)."""
    raw = n_jobs
    if raw is None and isinstance(config, dict):
        raw = config.get("cores", config.get("threads", config.get("n_jobs", 1)))
    try:
        jobs = int(raw if raw is not None else 1)
    except (TypeError, ValueError):
        jobs = 1
    return max(1, jobs)


class GpuNotSupported(RuntimeError):
    """GPU mode was requested, but this model's ``fit`` cannot run on CUDA."""


# cuML estimators train on the GPU. sklearn's fit for these models does not.
_GPU_CLASSES = {
    "RandomForest": ("cuml.ensemble", "RandomForestClassifier"),
    "LogisticRegression": ("cuml.linear_model", "LogisticRegression"),
}


def use_gpu(gpu=None, config=None) -> bool:
    if gpu is not None:
        return bool(gpu)
    if isinstance(config, dict):
        return bool(config.get("gpu") or config.get("device") == "cuda")
    return False


def cuda_device_available() -> bool:
    ensure_cuda_toolkit_path()
    try:
        cupy = importlib.import_module("cupy")
        return int(cupy.cuda.runtime.getDeviceCount()) > 0
    except Exception:
        pass
    try:
        torch = importlib.import_module("torch")
        return bool(torch.cuda.is_available())
    except Exception:
        return False


def trains_on_gpu(model) -> bool:
    """True only when this estimator's implementation runs ``fit`` on CUDA.

    A ``device`` attribute or a ``cuml`` class counts. sklearn ensembles and
    logistic regression stay on CPU even if a GPU is present.
    """
    module = str(getattr(type(model), "__module__", "") or "")
    if module.startswith("cuml"):
        return True
    device = getattr(model, "device", None)
    return device is not None and "cuda" in str(device).lower()


def gpu_training_support(kind: str) -> Optional[str]:
    """``None`` when ``kind`` can train on a GPU, otherwise the reason it cannot."""
    if kind == "AdaBoost":
        return (
            "AdaBoostClassifier.fit is a CPU loop: each round refits a sklearn "
            "decision stump on reweighted rows, and there is no CUDA implementation."
        )
    spec = _GPU_CLASSES.get(kind)
    if spec is None:
        return f"{kind} has no GPU training backend"
    module, _cls_name = spec
    try:
        importlib.import_module(module)
    except ImportError:
        return (
            f"{kind} GPU training needs {module} (cuML). "
            "The sklearn estimator's fit does not use the GPU."
        )
    if not cuda_device_available():
        return "CUDA device is not available"
    return None


def _load_gpu_class(module: str, name: str):
    return getattr(importlib.import_module(module), name)


def _supported_kwargs(cls, kwargs: dict) -> dict:
    params = inspect.signature(cls.__init__).parameters
    if any(param.kind == inspect.Parameter.VAR_KEYWORD for param in params.values()):
        return dict(kwargs)
    allowed = {key for key in params if key != "self"}
    picked = {key: value for key, value in kwargs.items() if key in allowed}
    if "output_type" in allowed and "output_type" not in picked:
        picked["output_type"] = "numpy"
    return picked


def ensure_cuda_toolkit_path() -> None:
    """Point CuPy at the conda CUDA toolkit before it is imported.

    CuPy's conda lookup uses ``CONDA_PREFIX/targets/<arch>``, and it caches
    that path on first import. ``samovar exec`` often runs the env interpreter
    without ``conda activate``, so ``CONDA_PREFIX`` is empty and NVRTC cannot
    find ``cuda.h``.
    """
    prefix = getattr(sys, "prefix", "") or ""
    if prefix and not os.environ.get("CONDA_PREFIX"):
        os.environ["CONDA_PREFIX"] = prefix
    current = os.environ.get("CUDA_PATH") or ""
    if current and os.path.isfile(os.path.join(current, "include", "cuda.h")):
        return
    roots = [
        os.path.join(prefix, "targets", "x86_64-linux"),
        os.path.join(prefix, "targets", "x86_64-linux-gnu"),
        prefix,
    ]
    for root in roots:
        if os.path.isfile(os.path.join(root, "include", "cuda.h")):
            os.environ["CUDA_PATH"] = root
            return


def build_gpu_estimator(kind: str, **kwargs):
    """Construct a CUDA estimator and reject anything that would still fit on CPU."""
    ensure_cuda_toolkit_path()
    reason = gpu_training_support(kind)
    if reason:
        raise GpuNotSupported(reason)
    module, cls_name = _GPU_CLASSES[kind]
    cls = _load_gpu_class(module, cls_name)
    model = cls(**_supported_kwargs(cls, kwargs))
    if not trains_on_gpu(model):
        raise GpuNotSupported(
            f"{kind} was built as {type(model).__module__}.{type(model).__name__}, "
            "which does not train on GPU."
        )
    return model


@contextmanager
def blas_threads(n_jobs: int):
    """Cap BLAS/OpenMP threads for solvers that do not take ``n_jobs`` (lbfgs)."""
    limit = training_n_jobs(n_jobs)
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        yield
        return
    with threadpool_limits(limits=limit):
        yield


class _ParallelDepth1Tree(BaseEstimator, ClassifierMixin):
    """Depth-1 tree whose split search runs in ``n_jobs`` threads.

    A depth-1 tree uses one feature. Searching each feature independently and
    keeping the lowest weighted impurity matches ``DecisionTreeClassifier(max_depth=1)``,
    which is the AdaBoost base learner. Boosting rounds stay sequential; this
    parallelizes the work inside each round. sklearn 1.7 AdaBoost no longer
    accepts ``n_jobs`` because that parameter did not affect ``fit``.
    """

    def __init__(self, n_jobs: int = 1, random_state=None):
        self.n_jobs = n_jobs
        self.random_state = random_state

    def fit(self, X, y, sample_weight=None):
        matrix = np.asarray(X)
        n_features = int(matrix.shape[1]) if matrix.ndim == 2 else 0
        if n_features < 1:
            raise ValueError("AdaBoost stump needs at least one feature")
        jobs = training_n_jobs(self.n_jobs)
        random_state = self.random_state

        def _one(index: int):
            tree = DecisionTreeClassifier(max_depth=1, random_state=random_state)
            column = matrix[:, [index]]
            tree.fit(column, y, sample_weight=sample_weight)
            layout = tree.tree_
            left = int(layout.children_left[0])
            right = int(layout.children_right[0])
            if left < 0:
                score = float(layout.impurity[0])
            else:
                weights = layout.weighted_n_node_samples
                score = float(
                    (weights[left] * layout.impurity[left] + weights[right] * layout.impurity[right])
                    / max(float(weights[0]), 1e-12)
                )
            return score, index, tree

        # Cython split search releases the GIL, so threads overlap the training work.
        found = Parallel(n_jobs=jobs, prefer="threads")(
            delayed(_one)(index) for index in range(n_features)
        )
        _score, index, tree = min(found, key=lambda row: (row[0], row[1]))
        importances = np.zeros(n_features, dtype=float)
        if int(tree.tree_.feature[0]) >= 0:
            importances[index] = 1.0
        self.tree_ = tree
        self.feature_index_ = int(index)
        self.feature_importances_ = importances
        self.classes_ = tree.classes_
        self.n_features_in_ = n_features
        return self

    def predict(self, X):
        column = np.asarray(X)[:, [self.feature_index_]]
        return self.tree_.predict(column)

    def predict_proba(self, X):
        column = np.asarray(X)[:, [self.feature_index_]]
        return self.tree_.predict_proba(column)


def build_random_forest(n_jobs=1, random_state=42, n_estimators=100, gpu=False):
    """Forest whose ``fit`` builds trees via ``Parallel(n_jobs=...)`` or cuML on GPU."""
    jobs = training_n_jobs(n_jobs)
    if use_gpu(gpu):
        return build_gpu_estimator(
            "RandomForest",
            n_estimators=n_estimators,
            random_state=random_state,
            n_streams=jobs,
            n_jobs=jobs,
        )
    return RandomForestClassifier(
        random_state=random_state,
        n_estimators=n_estimators,
        n_jobs=jobs,
    )


def build_adaboost(n_jobs=1, random_state=42, n_estimators=100, gpu=False) -> AdaBoostClassifier:
    """AdaBoost. ``n_jobs>1`` parallelizes each round's stump search. No GPU fit."""
    if use_gpu(gpu):
        raise GpuNotSupported(gpu_training_support("AdaBoost"))
    jobs = training_n_jobs(n_jobs)
    if jobs <= 1:
        return AdaBoostClassifier(random_state=random_state, n_estimators=n_estimators)
    return AdaBoostClassifier(
        estimator=_ParallelDepth1Tree(n_jobs=jobs, random_state=random_state),
        random_state=random_state,
        n_estimators=n_estimators,
    )


def fit_logistic(X, y, *, n_jobs=1, gpu=False, **kwargs) -> Any:
    """Fit logistic regression using ``n_jobs`` on the training itself.

    Multiclass with more than one thread uses one-vs-rest so each class is a
    separate ``fit`` (sklearn's multinomial lbfgs ignores ``n_jobs``). Binary
    lbfgs is one problem; BLAS/OpenMP threads are capped at ``n_jobs``.
    ``gpu=True`` fits a cuML model on CUDA instead of sklearn.
    """
    from sklearn.linear_model import LogisticRegression

    if use_gpu(gpu):
        model = build_gpu_estimator("LogisticRegression", **kwargs)
        model.fit(X, y)
        return model
    nbytes = int(getattr(X, "shape", (0, 0))[0] or 0) * int(getattr(X, "shape", (0, 1))[1] or 0) * 8
    jobs = parallel_jobs(training_n_jobs(n_jobs), nbytes=nbytes)
    labels = pd.Series(np.asarray(y)).dropna()
    n_classes = int(labels.nunique())
    if jobs > 1 and n_classes > 2:
        model = OneVsRestClassifier(LogisticRegression(**kwargs), n_jobs=jobs)
        model.fit(X, y)
        return model
    model = LogisticRegression(**kwargs)
    with blas_threads(jobs):
        model.fit(X, y)
    return model


def train_models(df, test_size=0.2, random_state=42, methods=None, n_jobs=1, gpu=False):
    """
    Train RandomForest and/or AdaBoost models for taxID prediction.

    Args:
        df (pd.DataFrame): Input dataframe
        test_size (float): Proportion of data to use for testing
        random_state (int): Random seed for reproducibility
        methods: ``random_forest`` (default), ``adaboost``, ``ensemble`` (both, pick best),
            or an iterable of those names.
        n_jobs: Threads for training. ``samovar prepare --cores`` supplies this.
        gpu: Train on CUDA. Each model is checked; sklearn estimators are rejected.

    Returns:
        tuple: (best_model, models_dict, metrics_dict, feature_cols)
    """
    # Preprocess the data first
    df_processed = preprocess_data(df)
    if df_processed.empty:
        raise ValueError("No training rows available after preprocessing")

    # Select features for training (exclude seq and true)
    feature_cols = sorted([col for col in df_processed.columns 
                         if col not in ['seq', 'true']])

    if 'length' not in feature_cols:
        df_processed['length'] = 0
        feature_cols.append('length')
        feature_cols = sorted(feature_cols)
    
    # Ensure all required features are present
    required_features = ['length'] + [col for col in feature_cols if col.startswith('taxid_')]
    missing_features = [f for f in required_features if f not in feature_cols]
    if missing_features:
        raise ValueError(f"Missing required features: {missing_features}")
    
    # Prepare features and target
    X = df_processed[feature_cols]
    y = df_processed['true']
    min_class = int(os.environ.get("SAMOVAR_MIN_CLASS_COUNT", "0") or 0)
    if min_class > 1:
        counts = y.value_counts()
        keep_labels = counts[counts >= min_class].index
        if len(keep_labels) >= 2:
            mask = y.isin(keep_labels)
            dropped = int((~mask).sum())
            print(
                f"[rfc] skip rare profiles: drop {dropped} rows in classes with n<{min_class}; "
                f"keep {len(keep_labels)} classes",
                flush=True,
            )
            X = X.loc[mask]
            y = y.loc[mask]
        else:
            print(
                f"[rfc] skip rare profiles: not applied, fewer than 2 classes have n>={min_class}",
                flush=True,
            )
    
    # Split data. Stratify only when every class has ≥2 members (sklearn).
    split_kwargs = {"test_size": test_size, "random_state": random_state}
    class_counts = y.value_counts()
    if len(class_counts) > 1 and int(class_counts.min()) >= 2:
        split_kwargs["stratify"] = y
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, **split_kwargs
    )

    wanted = _normalize_sklearn_methods(methods)
    if use_gpu(gpu):
        ensure_cuda_toolkit_path()
    matrix_bytes = int(getattr(X_train, "shape", (0, 0))[0] or 0) * int(
        getattr(X_train, "shape", (0, 1))[1] or 0
    ) * 8
    jobs = parallel_jobs(training_n_jobs(n_jobs), nbytes=matrix_bytes)
    on_gpu = use_gpu(gpu)
    if on_gpu:
        for name in wanted:
            reason = gpu_training_support(name)
            if reason:
                raise GpuNotSupported(reason)
    builders = {
        "RandomForest": lambda: build_random_forest(
            n_jobs=jobs, random_state=random_state, gpu=on_gpu
        ),
        "AdaBoost": lambda: build_adaboost(
            n_jobs=jobs, random_state=random_state, gpu=on_gpu
        ),
    }
    # Build only the requested estimators. AdaBoost has no CUDA fit, so
    # constructing it whenever gpu=True rejects a RandomForest-only run.
    models = {name: builders[name]() for name in wanted}
    
    # Train and evaluate models
    metrics = {}
    best_score = -np.inf
    best_model = None
    
    for name, model in models.items():
        model.fit(X_train, y_train)
        model.samovar_preprocessing_ = df_processed.attrs["samovar_preprocessing"]
        score = model.score(X_test, y_test)
        metrics[name] = score
        
        if score > best_score:
            best_score = score
            best_model = model
    
    return best_model, models, metrics, feature_cols


def _normalize_sklearn_methods(methods) -> list:
    alias = {
        "ensemble": ("RandomForest", "AdaBoost"),
        "samovar": ("RandomForest", "AdaBoost"),
        "default": ("RandomForest", "AdaBoost"),
        "both": ("RandomForest", "AdaBoost"),
        "random_forest": ("RandomForest",),
        "rf": ("RandomForest",),
        "randomforest": ("RandomForest",),
        "adaboost": ("AdaBoost",),
        "ada": ("AdaBoost",),
        "randomforestclassifier": ("RandomForest",),
        "adaboostclassifier": ("AdaBoost",),
    }
    if methods is None or methods is False:
        names = ["random_forest"]
    elif isinstance(methods, str):
        names = [methods]
    else:
        names = list(methods)
    out = []
    seen = set()
    for raw in names:
        key = str(raw or "ensemble").strip().lower().replace("-", "_")
        mapped = alias.get(key)
        if mapped is None:
            raise ValueError(
                f"Unknown built-in reprofiler {raw!r}. Use ensemble, random_forest, or adaboost."
            )
        for name in mapped:
            if name not in seen:
                seen.add(name)
                out.append(name)
    return out or ["RandomForest", "AdaBoost"]


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

def _predict_in_chunks(model, X: pd.DataFrame):
    """Predict one block at a time so the reads × classes matrix is not retained."""
    n = len(X)
    predictions = np.empty(n, dtype=object)
    confidences = np.zeros(n, dtype=float)
    if n == 0:
        return predictions, confidences
    classes = getattr(model, "classes_", None)
    n_classes = len(classes) if classes is not None else 2
    chunk = proba_chunk_rows(n_classes, int(X.shape[1]))
    annotator_cols = annotator_taxid_columns(X)
    requested = getattr(model, "n_jobs", 1) or 1
    with predict_threads(model, requested):
        for start in range(0, n, chunk):
            stop = min(n, start + chunk)
            block = X.iloc[start:stop]
            pred = np.asarray(model.predict(block))
            proba = np.asarray(model.predict_proba(block))
            conf = proba.max(axis=1) if proba.size else np.zeros(stop - start)
            del proba
            values = block[annotator_cols] if annotator_cols else None
            for i, label in enumerate(pred):
                votes = []
                if values is not None:
                    for col in annotator_cols:
                        try:
                            value = int(values.iloc[i][col])
                        except (TypeError, ValueError):
                            continue
                        if value > 0:
                            votes.append(value)
                confidence = float(conf[i])
                if os.environ.get("SAMOVAR_RFC_CORRECTION") == "1":
                    min_conf = float(os.environ.get("SAMOVAR_MIN_CONFIDENCE", "0.5"))
                    if not votes or confidence < min_conf:
                        predictions[start + i] = 0
                        confidences[start + i] = confidence
                    else:
                        predictions[start + i] = int(label)
                        confidences[start + i] = confidence
                    continue
                if votes:
                    top_value, top_count = Counter(votes).most_common(1)[0]
                    if top_count == len(votes):
                        predictions[start + i] = top_value
                        confidences[start + i] = 1.0 if len(votes) > 1 else confidence
                        continue
                predictions[start + i] = int(label)
                confidences[start + i] = confidence
    return predictions, confidences


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
    
    # Load or use provided model
    if isinstance(model_path, (str, os.PathLike)) and os.path.exists(model_path):
        model = joblib.load(model_path)
    elif model_path is not None:
        model = model_path
    else:
        model, _, _, feature_cols = train_models(df)

    df_processed = preprocess_for_model(df, model, drop_missing_true=False)
    if feature_cols is None:
        feature_cols = sorted(col for col in df_processed.columns if col != 'true')

    if hasattr(model, "feature_names_in_"):
        feature_cols = list(model.feature_names_in_)
    if "length" not in df_processed.columns:
        df_processed["length"] = 0
    X = df_processed.reindex(columns=feature_cols, fill_value=0)
    final_predictions, final_confidences = _predict_in_chunks(model, X)
    
    # Chunked read_csv keeps the file index (20000, 20001, ...), so labels are
    # not positions in this dataframe.
    processed_to_original = dict(zip(range(len(df_processed)), df_processed.index))
    label_to_pos = {label: pos for pos, label in enumerate(df.index)}

    full_predictions = np.zeros(len(df), dtype=object)
    full_confidences = np.zeros(len(df))

    for i, (pred, conf) in enumerate(zip(final_predictions, final_confidences)):
        pos = label_to_pos[processed_to_original[i]]
        full_predictions[pos] = pred
        full_confidences[pos] = conf
    
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


def passthrough_reprofile_tables(
    reprofiling_dir: str,
    output_dir: str,
    features_df=None,
) -> int:
    """Write ``*_reprofiled.csv`` without a model so viz_reprofiled can run.

    ``taxid_SAMOVAR`` copies ``true`` when present, otherwise 0.
    """
    os.makedirs(output_dir, exist_ok=True)
    written = 0
    names = sorted(os.listdir(reprofiling_dir))
    per_sample = [
        name
        for name in names
        if name.endswith(".csv") and not name.startswith("combined_annotation_table")
    ]
    if not per_sample:
        for name in names:
            if name.startswith("combined_annotation_table") and name.endswith(".csv"):
                written += _passthrough_combined(
                    os.path.join(reprofiling_dir, name), output_dir, features_df
                )
        return written
    for filename in per_sample:
        sample_file = os.path.join(reprofiling_dir, filename)
        df = pd.read_csv(sample_file)
        _write_passthrough_frame(df, output_dir, os.path.basename(sample_file).split(".")[0], features_df)
        written += 1
    return written


def _write_passthrough_frame(df, output_dir, stem, features_df):
    df = merge_read_features(df, features_df)
    if "sample" in df.columns:
        df = df.drop(columns=["sample"])
    if "taxid_SAMOVAR" not in df.columns:
        if "true" in df.columns:
            df["taxid_SAMOVAR"] = df["true"]
        else:
            df["taxid_SAMOVAR"] = 0
        df["taxid_SAMOVAR_confidence"] = 0
    output_file = os.path.join(output_dir, f"{stem}_reprofiled.csv")
    df.to_csv(output_file, index=False)


def _passthrough_combined(path, output_dir, features_df) -> int:
    seen = set()
    n = 0
    for chunk in pd.read_csv(path, chunksize=20_000):
        if "sample" not in chunk.columns:
            groups = [(os.path.basename(path).split(".")[0], chunk)]
        else:
            groups = list(chunk.groupby("sample", sort=False))
        for sample, group in groups:
            key = str(sample)
            stem = key[:-11] if key.endswith("_reprofiled") else key
            out = os.path.join(output_dir, f"{stem}_reprofiled.csv")
            frame = merge_read_features(group.drop(columns=["sample"], errors="ignore"), features_df)
            if "taxid_SAMOVAR" not in frame.columns:
                if "true" in frame.columns:
                    frame["taxid_SAMOVAR"] = frame["true"]
                else:
                    frame["taxid_SAMOVAR"] = 0
                frame["taxid_SAMOVAR_confidence"] = 0
            first = key not in seen
            frame.to_csv(out, index=False, mode="w" if first else "a", header=first)
            seen.add(key)
            if first:
                n += 1
    return n
