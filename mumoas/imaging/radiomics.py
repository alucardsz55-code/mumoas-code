from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd
from pandas.api.types import is_numeric_dtype


def select_feature_columns(
    table: pd.DataFrame,
    patient_col: str = "patient_id",
    preserve_columns: Sequence[str] | None = None,
) -> list[str]:
    """Return patient id, explicitly preserved columns, and numeric feature columns."""
    if patient_col not in table.columns:
        raise ValueError(f"Patient identifier column '{patient_col}' is missing from feature table")

    preserved = [col for col in (preserve_columns or []) if col != patient_col]
    missing_preserved = [col for col in preserved if col not in table.columns]
    if missing_preserved:
        raise ValueError(f"Preserved columns are missing from feature table: {missing_preserved}")

    numeric_features = [
        col
        for col in table.columns
        if col != patient_col and col not in preserved and is_numeric_dtype(table[col])
    ]
    if not numeric_features:
        raise ValueError("Feature table does not contain any numeric feature columns")

    return [patient_col, *preserved, *numeric_features]


def validate_feature_table(
    table: pd.DataFrame,
    patient_col: str = "patient_id",
    preserve_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Keep patient identifiers and usable numeric features; drop text metadata by default."""
    columns = select_feature_columns(table, patient_col=patient_col, preserve_columns=preserve_columns)
    features = table.loc[:, columns].copy()
    if features[patient_col].isna().any():
        raise ValueError(f"Patient identifier column '{patient_col}' contains missing values")
    if features[patient_col].duplicated().any():
        raise ValueError(f"Patient identifier column '{patient_col}' contains duplicate values")
    return features


def merge_sequence_feature_tables(
    tables: Mapping[str, pd.DataFrame] | Sequence[pd.DataFrame],
    patient_col: str = "patient_id",
    preserve_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Validate and outer-merge sequence-specific radiomics tables by patient id."""
    if isinstance(tables, Mapping):
        named_tables = list(tables.items())
    else:
        named_tables = [(f"seq{idx + 1}", table) for idx, table in enumerate(tables)]
    if not named_tables:
        raise ValueError("At least one feature table is required")

    merged: pd.DataFrame | None = None
    for name, table in named_tables:
        validated = validate_feature_table(
            table, patient_col=patient_col, preserve_columns=preserve_columns
        )
        rename_map = {
            col: f"{name}_{col}"
            for col in validated.columns
            if col != patient_col and col not in (preserve_columns or [])
        }
        validated = validated.rename(columns=rename_map)
        merged = validated if merged is None else merged.merge(validated, on=patient_col, how="outer")

    if merged is None:
        raise ValueError("At least one feature table is required")
    return merged
