from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2


def _binary_series(y_true: pd.Series | list[int] | list[str], positive_label: str | int | None = None) -> pd.Series:
    true = pd.Series(y_true)
    if positive_label is None:
        values = set(true.dropna().unique().tolist())
        if not values.issubset({0, 1, False, True}):
            raise ValueError("positive_label is required when y_true is not binary 0/1.")
        return true.astype(int)
    return (true == positive_label).astype(int)


def calibration_table(
    y_true: pd.Series | list[int] | list[str],
    y_pred: pd.Series | list[float],
    n_bins: int = 10,
    positive_label: str | int | None = None,
) -> pd.DataFrame:
    if n_bins < 1:
        raise ValueError("n_bins must be at least 1.")
    observed = _binary_series(y_true, positive_label)
    predicted = pd.Series(y_pred, dtype=float)
    if len(observed) != len(predicted):
        raise ValueError("y_true and y_pred must have the same length.")
    if len(observed) == 0:
        raise ValueError("y_true and y_pred must not be empty.")
    if predicted.lt(0).any() or predicted.gt(1).any():
        raise ValueError("Predicted probabilities must be in [0, 1].")

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins = pd.cut(predicted, bins=edges, include_lowest=True, labels=False)
    grouped = (
        pd.DataFrame({"bin": bins, "predicted": predicted, "observed": observed})
        .dropna(subset=["bin"])
        .groupby("bin", observed=True)
        .agg(mean_predicted=("predicted", "mean"), observed_rate=("observed", "mean"), n=("observed", "size"))
        .reset_index()
    )
    grouped["bin"] = grouped["bin"].astype(int) + 1
    return grouped[["bin", "mean_predicted", "observed_rate", "n"]]


def hosmer_lemeshow_summary(
    y_true: pd.Series | list[int] | list[str],
    y_pred: pd.Series | list[float],
    n_bins: int = 10,
    positive_label: str | int | None = None,
) -> pd.DataFrame:
    if n_bins < 2:
        raise ValueError("n_bins must be at least 2 for a Hosmer-Lemeshow summary.")
    observed = _binary_series(y_true, positive_label)
    predicted = pd.Series(y_pred, dtype=float)
    if len(observed) != len(predicted):
        raise ValueError("y_true and y_pred must have the same length.")
    if len(observed) == 0:
        raise ValueError("y_true and y_pred must not be empty.")
    if predicted.lt(0).any() or predicted.gt(1).any():
        raise ValueError("Predicted probabilities must be in [0, 1].")
    table = pd.DataFrame({"observed": observed, "predicted": predicted})
    table["bin"] = pd.qcut(table["predicted"].rank(method="first"), q=min(n_bins, len(table)), labels=False) + 1
    grouped = table.groupby("bin", observed=True).agg(
        observed_events=("observed", "sum"),
        expected_events=("predicted", "sum"),
        n=("observed", "size"),
    )
    grouped["observed_nonevents"] = grouped["n"] - grouped["observed_events"]
    grouped["expected_nonevents"] = grouped["n"] - grouped["expected_events"]

    eps = np.finfo(float).eps
    statistic = (
        ((grouped["observed_events"] - grouped["expected_events"]) ** 2 / grouped["expected_events"].clip(lower=eps))
        + (
            (grouped["observed_nonevents"] - grouped["expected_nonevents"]) ** 2
            / grouped["expected_nonevents"].clip(lower=eps)
        )
    ).sum()
    degrees = max(int(grouped.shape[0]) - 2, 1)
    return pd.DataFrame(
        [
            {
                "n_groups": int(grouped.shape[0]),
                "statistic": float(statistic),
                "df": degrees,
                "p_value": float(chi2.sf(statistic, degrees)),
            }
        ]
    )
