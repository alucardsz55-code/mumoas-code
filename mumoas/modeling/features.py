from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, StandardScaler


@dataclass(frozen=True)
class DesignMatrix:
    features: pd.DataFrame
    conditions: pd.DataFrame
    labels: np.ndarray
    patient_ids: pd.Series
    metadata: dict[str, Any]
    encoders: dict[str, Any]


def _numeric_frame(features: pd.DataFrame) -> pd.DataFrame:
    numeric = features.select_dtypes(include=[np.number]).copy()
    if numeric.empty:
        raise ValueError("Feature table must contain at least one numeric column")
    return numeric.replace([np.inf, -np.inf], np.nan).fillna(0.0)


def prune_correlated_features(
    features: pd.DataFrame, threshold: float = 0.9
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Drop later columns from highly correlated Spearman feature pairs."""
    numeric = _numeric_frame(features)
    if not 0 <= threshold <= 1:
        raise ValueError("Correlation threshold must be between 0 and 1")

    ranked = numeric.rank(method="average")
    corr = ranked.corr(method="pearson").abs().fillna(0.0)
    dropped: list[str] = []
    dropped_set: set[str] = set()

    columns = list(numeric.columns)
    for i, left in enumerate(columns):
        if left in dropped_set:
            continue
        for right in columns[i + 1 :]:
            if right in dropped_set:
                continue
            if corr.loc[left, right] > threshold:
                dropped.append(right)
                dropped_set.add(right)

    kept = [column for column in columns if column not in dropped_set]
    return numeric.loc[:, kept].copy(), kept, dropped


def _vif_for_column(values: np.ndarray, column_index: int) -> float:
    y = values[:, column_index]
    if np.nanstd(y) <= 1e-12:
        return np.inf

    x = np.delete(values, column_index, axis=1)
    if x.shape[1] == 0:
        return 1.0

    x = np.column_stack([np.ones(x.shape[0]), x])
    if x.shape[0] <= x.shape[1]:
        return 1.0
    try:
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        fitted = x @ beta
    except np.linalg.LinAlgError:
        return np.inf

    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot <= 1e-12:
        return np.inf
    r_squared = 1.0 - ss_res / ss_tot
    if not np.isfinite(r_squared):
        return np.inf
    r_squared = min(max(r_squared, 0.0), 1.0)
    if r_squared >= 1.0 - 1e-12:
        return np.inf
    return 1.0 / (1.0 - r_squared)


def prune_vif(
    features: pd.DataFrame, threshold: float = 10.0
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Iteratively remove the feature with the largest variance inflation factor."""
    if threshold <= 0:
        raise ValueError("VIF threshold must be positive")

    numeric = _numeric_frame(features)
    kept = list(numeric.columns)
    dropped: list[str] = []

    while len(kept) > 1:
        values = numeric.loc[:, kept].to_numpy(dtype=float)
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
        vifs = np.array([_vif_for_column(values, idx) for idx in range(len(kept))])
        worst_index = int(np.argmax(vifs))
        worst_vif = float(vifs[worst_index])
        if worst_vif <= threshold:
            break
        dropped.append(kept.pop(worst_index))

    return numeric.loc[:, kept].copy(), kept, dropped


def _make_one_hot_encoder() -> OneHotEncoder:
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def _encoded_column_names(encoder: OneHotEncoder, columns: list[str]) -> list[str]:
    if not columns:
        return []
    try:
        return list(encoder.get_feature_names_out(columns))
    except AttributeError:
        return list(encoder.get_feature_names(columns))


def build_design_matrix(
    clinical: pd.DataFrame,
    imaging: pd.DataFrame,
    patient_col: str,
    label_col: str,
    categorical_cols: list[str] | tuple[str, ...],
    continuous_cols: list[str] | tuple[str, ...],
    label_order: list[str] | tuple[str, ...] | None = None,
) -> DesignMatrix:
    """Join clinical and imaging tables into model-ready feature and condition matrices."""
    if patient_col not in clinical.columns or patient_col not in imaging.columns:
        raise ValueError(f"Both tables must include patient column '{patient_col}'")
    if label_col not in clinical.columns:
        raise ValueError(f"Clinical table must include label column '{label_col}'")

    categorical = [column for column in categorical_cols if column in clinical.columns]
    continuous = [column for column in continuous_cols if column in clinical.columns]
    missing_categorical = sorted(set(categorical_cols) - set(categorical))
    missing_continuous = sorted(set(continuous_cols) - set(continuous))

    merged = clinical.merge(imaging, on=patient_col, how="inner", validate="one_to_one")
    if merged.empty:
        raise ValueError("Clinical and imaging tables have no overlapping patients")

    imaging_columns = [
        column
        for column in imaging.columns
        if column != patient_col and pd.api.types.is_numeric_dtype(imaging[column])
    ]
    if not imaging_columns and not continuous and not categorical:
        raise ValueError("No usable numeric or categorical predictors were found")

    numeric_columns = continuous + imaging_columns
    numeric_scaler = StandardScaler()
    numeric_values = pd.DataFrame(index=merged.index)
    if numeric_columns:
        numeric_raw = merged.loc[:, numeric_columns].apply(pd.to_numeric, errors="coerce")
        numeric_raw = numeric_raw.replace([np.inf, -np.inf], np.nan).fillna(numeric_raw.median())
        numeric_raw = numeric_raw.fillna(0.0)
        numeric_values = pd.DataFrame(
            numeric_scaler.fit_transform(numeric_raw),
            columns=numeric_columns,
            index=merged.index,
        )

    categorical_encoder = _make_one_hot_encoder()
    categorical_values = pd.DataFrame(index=merged.index)
    if categorical:
        categorical_raw = merged.loc[:, categorical].fillna("missing").astype(str)
        encoded = categorical_encoder.fit_transform(categorical_raw)
        categorical_values = pd.DataFrame(
            encoded,
            columns=_encoded_column_names(categorical_encoder, categorical),
            index=merged.index,
        )

    feature_table = pd.concat([numeric_values, categorical_values], axis=1)
    condition_columns = continuous + list(categorical_values.columns)
    condition_table = pd.concat(
        [numeric_values.loc[:, continuous] if continuous else pd.DataFrame(index=merged.index), categorical_values],
        axis=1,
    )

    observed_labels = merged[label_col].astype(str)
    if label_order:
        label_classes = list(label_order)
        missing_from_order = sorted(set(observed_labels) - set(label_classes))
        if missing_from_order:
            raise ValueError(
                f"Observed labels are not present in configured label_order: {missing_from_order}"
            )
    else:
        label_classes = sorted(observed_labels.dropna().unique().tolist())
    if not label_classes:
        raise ValueError(f"Label column '{label_col}' does not contain any usable labels")
    label_to_index = {label: index for index, label in enumerate(label_classes)}
    labels = observed_labels.map(label_to_index)
    if labels.isna().any():
        missing_labels = sorted(observed_labels[labels.isna()].unique().tolist())
        raise ValueError(f"Could not encode labels: {missing_labels}")
    labels_array = labels.astype(int).to_numpy()

    metadata = {
        "feature_columns": list(feature_table.columns),
        "condition_columns": condition_columns,
        "continuous_columns": continuous,
        "categorical_columns": categorical,
        "imaging_columns": imaging_columns,
        "label_classes": label_classes,
        "missing_categorical_columns": missing_categorical,
        "missing_continuous_columns": missing_continuous,
    }
    encoders = {
        "numeric_scaler": numeric_scaler,
        "categorical_encoder": categorical_encoder,
        "label_to_index": label_to_index,
    }

    return DesignMatrix(
        features=feature_table.astype(float),
        conditions=condition_table.astype(float),
        labels=labels_array,
        patient_ids=merged[patient_col].reset_index(drop=True),
        metadata=metadata,
        encoders=encoders,
    )
