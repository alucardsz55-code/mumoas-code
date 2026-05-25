from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve
from sklearn.preprocessing import label_binarize


def probability_columns(table: pd.DataFrame, classes: list[str] | None = None) -> list[str]:
    if classes is not None:
        expected = [f"prob_{label}" for label in classes]
        missing = [column for column in expected if column not in table.columns]
        if missing:
            raise ValueError(f"Missing probability columns: {missing}")
        return expected
    return [column for column in table.columns if column.startswith("prob_")]


def decode_prediction_labels(
    probabilities: pd.DataFrame,
    classes: list[str] | None = None,
    prefix: str = "prob_",
) -> pd.Series:
    prob_cols = probability_columns(probabilities, classes)
    if not prob_cols:
        raise ValueError("No probability columns found; expected columns named like 'prob_<class>'.")
    labels = [column.removeprefix(prefix) for column in prob_cols]
    predicted_index = probabilities[prob_cols].to_numpy(dtype=float).argmax(axis=1)
    return pd.Series([labels[idx] for idx in predicted_index], index=probabilities.index, name="predicted_class")


def confusion_matrix_table(
    y_true: pd.Series | list[str],
    y_pred: pd.Series | list[str],
    labels: list[str] | None = None,
) -> pd.DataFrame:
    true = pd.Series(y_true)
    predicted = pd.Series(y_pred)
    if labels is None:
        labels = sorted(set(true.dropna().tolist()) | set(predicted.dropna().tolist()), key=str)
    matrix = confusion_matrix(true, predicted, labels=labels)
    rows = []
    for row_label, counts in zip(labels, matrix):
        for column_label, count in zip(labels, counts):
            rows.append({"true_label": row_label, "predicted_label": column_label, "n": int(count)})
    return pd.DataFrame(rows)


def multiclass_auc_table(
    y_true: pd.Series | list[str],
    probabilities: pd.DataFrame,
    classes: list[str] | None = None,
) -> pd.DataFrame:
    true = pd.Series(y_true)
    if classes is None:
        classes = [column.removeprefix("prob_") for column in probability_columns(probabilities)]
    rows = []
    for label in classes:
        column = f"prob_{label}"
        if column not in probabilities.columns:
            raise ValueError(f"Missing probability column '{column}'.")
        binary_true = (true == label).astype(int)
        n_positive = int(binary_true.sum())
        n_negative = int(len(binary_true) - n_positive)
        auc = np.nan
        if n_positive > 0 and n_negative > 0:
            auc = float(roc_auc_score(binary_true, probabilities[column].to_numpy(dtype=float)))
        rows.append(
            {
                "class": label,
                "auc": auc,
                "n_positive": n_positive,
                "n_negative": n_negative,
            }
        )
    return pd.DataFrame(rows)


def one_vs_rest_roc_table(
    y_true: pd.Series | list[str],
    probabilities: pd.DataFrame,
    classes: list[str] | None = None,
) -> pd.DataFrame:
    true = pd.Series(y_true)
    if classes is None:
        classes = [column.removeprefix("prob_") for column in probability_columns(probabilities)]
    rows = []
    for label in classes:
        column = f"prob_{label}"
        if column not in probabilities.columns:
            raise ValueError(f"Missing probability column '{column}'.")
        binary_true = (true == label).astype(int)
        if binary_true.nunique() < 2:
            continue
        fpr, tpr, thresholds = roc_curve(binary_true, probabilities[column].to_numpy(dtype=float))
        for fp_rate, tp_rate, threshold in zip(fpr, tpr, thresholds):
            rows.append(
                {
                    "class": label,
                    "fpr": float(fp_rate),
                    "tpr": float(tp_rate),
                    "threshold": float(threshold),
                }
            )
    return pd.DataFrame(rows)


def macro_micro_auc(
    y_true: pd.Series | list[str],
    probabilities: pd.DataFrame,
    classes: list[str] | None = None,
) -> pd.DataFrame:
    true = pd.Series(y_true)
    if classes is None:
        classes = [column.removeprefix("prob_") for column in probability_columns(probabilities)]
    prob_cols = [f"prob_{label}" for label in classes]
    missing = [column for column in prob_cols if column not in probabilities.columns]
    if missing:
        raise ValueError(f"Missing probability columns: {missing}")

    y_binary = label_binarize(true, classes=classes)
    if y_binary.shape[1] == 1:
        y_binary = np.column_stack([1 - y_binary[:, 0], y_binary[:, 0]])
    scores = probabilities[prob_cols].to_numpy(dtype=float)

    rows = []
    for average in ("macro", "micro"):
        auc = np.nan
        try:
            auc = float(roc_auc_score(y_binary, scores, average=average))
        except ValueError:
            pass
        rows.append({"average": average, "auc": auc})
    return pd.DataFrame(rows)
