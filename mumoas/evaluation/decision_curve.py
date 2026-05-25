from __future__ import annotations

import numpy as np
import pandas as pd


def _validate_thresholds(thresholds: list[float] | np.ndarray) -> np.ndarray:
    values = np.asarray(thresholds, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("thresholds must be a non-empty one-dimensional sequence.")
    if np.any(values <= 0) or np.any(values >= 1):
        raise ValueError("Decision-curve thresholds must be between 0 and 1.")
    return values


def net_benefit_curve(
    y_true: pd.Series | list[int] | list[str],
    y_pred: pd.Series | list[float],
    thresholds: list[float] | np.ndarray,
    positive_label: str | int | None = None,
) -> pd.DataFrame:
    true = pd.Series(y_true)
    if positive_label is None:
        values = set(true.dropna().unique().tolist())
        if not values.issubset({0, 1, False, True}):
            raise ValueError("positive_label is required when y_true is not binary 0/1.")
        binary = true.astype(int).to_numpy()
    else:
        binary = (true == positive_label).astype(int).to_numpy()
    predicted = np.asarray(y_pred, dtype=float)
    if len(binary) != len(predicted):
        raise ValueError("y_true and y_pred must have the same length.")
    if len(binary) == 0:
        raise ValueError("y_true and y_pred must not be empty.")
    if np.any(predicted < 0) or np.any(predicted > 1):
        raise ValueError("Predicted probabilities must be in [0, 1].")

    thresholds_array = _validate_thresholds(thresholds)
    n = len(binary)
    prevalence = float(binary.mean()) if n else 0.0
    rows = []
    for threshold in thresholds_array:
        weight = threshold / (1.0 - threshold)
        predicted_positive = predicted >= threshold
        tp = float(((binary == 1) & predicted_positive).sum())
        fp = float(((binary == 0) & predicted_positive).sum())
        model_nb = tp / n - fp / n * weight
        treat_all_nb = prevalence - (1.0 - prevalence) * weight
        rows.extend(
            [
                {"threshold": float(threshold), "strategy": "model", "net_benefit": float(model_nb)},
                {"threshold": float(threshold), "strategy": "treat_all", "net_benefit": float(treat_all_nb)},
                {"threshold": float(threshold), "strategy": "treat_none", "net_benefit": 0.0},
            ]
        )
    return pd.DataFrame(rows)
