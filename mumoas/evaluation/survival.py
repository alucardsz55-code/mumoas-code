from __future__ import annotations

import pandas as pd


def prepare_survival_table(
    clinical: pd.DataFrame,
    predictions: pd.DataFrame | None = None,
    patient_col: str = "patient_id",
    time_col: str = "efs_time",
    event_col: str = "efs_event",
    prediction_group_col: str = "predicted_class",
    output_group_col: str = "predicted_group",
) -> pd.DataFrame:
    required = [patient_col, time_col, event_col]
    missing = [column for column in required if column not in clinical.columns]
    if missing:
        raise ValueError(f"Missing survival columns in clinical table: {missing}")

    table = clinical[required].copy()
    table[time_col] = pd.to_numeric(table[time_col], errors="coerce")
    table[event_col] = pd.to_numeric(table[event_col], errors="coerce")
    if table[[time_col, event_col]].isna().any().any():
        raise ValueError(f"Survival columns '{time_col}' and '{event_col}' must be numeric and non-missing.")
    if (table[time_col] <= 0).any():
        raise ValueError(f"Survival time column '{time_col}' must be positive.")
    if not set(table[event_col].unique().tolist()).issubset({0, 1, 0.0, 1.0}):
        raise ValueError(f"Survival event column '{event_col}' must contain 0/1 values.")

    if predictions is not None:
        required_prediction = [patient_col, prediction_group_col]
        missing_prediction = [column for column in required_prediction if column not in predictions.columns]
        if missing_prediction:
            raise ValueError(f"Missing prediction columns for survival grouping: {missing_prediction}")
        groups = predictions[required_prediction].rename(columns={prediction_group_col: output_group_col})
        table = table.merge(groups, on=patient_col, how="left")
    elif output_group_col not in table.columns:
        table[output_group_col] = "all"
    return table


def fit_kaplan_meier_by_group(
    survival_table: pd.DataFrame,
    time_col: str,
    event_col: str,
    group_col: str = "predicted_group",
) -> dict[str, object]:
    try:
        from lifelines import KaplanMeierFitter
    except ImportError as exc:
        raise ImportError("Install lifelines to fit Kaplan-Meier curves.") from exc

    fits: dict[str, object] = {}
    for group, group_table in survival_table.groupby(group_col, dropna=False):
        fitter = KaplanMeierFitter()
        fitter.fit(group_table[time_col], event_observed=group_table[event_col], label=str(group))
        fits[str(group)] = fitter
    return fits
