from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from scipy.stats import spearmanr


def fdr_bh(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values aligned to the input order."""
    p = pd.to_numeric(pd.Series(p_values), errors="coerce").to_numpy(dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return adjusted

    valid_p = p[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n_tests = ranked.size
    raw = ranked * n_tests / np.arange(1, n_tests + 1)
    monotonic = np.minimum.accumulate(raw[::-1])[::-1]
    monotonic = np.clip(monotonic, 0.0, 1.0)

    valid_adjusted = np.empty_like(valid_p)
    valid_adjusted[order] = monotonic
    adjusted[np.where(valid)[0]] = valid_adjusted
    return adjusted


def _numeric_feature_columns(table: pd.DataFrame, patient_col: str, exclude: set[str]) -> list[str]:
    return [
        column
        for column in table.columns
        if column != patient_col and column not in exclude and is_numeric_dtype(table[column])
    ]


def _require_columns(table: pd.DataFrame, columns: Sequence[str], table_name: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def spearman_program_associations(
    imaging: pd.DataFrame,
    programs: pd.DataFrame,
    patient_col: str,
    program_cols: Sequence[str],
) -> pd.DataFrame:
    """Correlate numeric imaging features with molecular program scores."""
    _require_columns(imaging, [patient_col], "imaging")
    _require_columns(programs, [patient_col, *program_cols], "programs")

    program_set = set(program_cols)
    feature_cols = _numeric_feature_columns(imaging, patient_col=patient_col, exclude=program_set)
    if not feature_cols:
        raise ValueError("imaging does not contain numeric feature columns")

    merged = imaging[[patient_col, *feature_cols]].merge(
        programs[[patient_col, *program_cols]],
        on=patient_col,
        how="inner",
    )
    rows: list[dict[str, float | str]] = []
    for feature in feature_cols:
        for program in program_cols:
            pair = merged[[feature, program]].dropna()
            rho = np.nan
            p_value = np.nan
            if len(pair) >= 2 and pair[feature].nunique() > 1 and pair[program].nunique() > 1:
                result = spearmanr(pair[feature], pair[program])
                rho = float(getattr(result, "statistic", result[0]))
                p_value = float(getattr(result, "pvalue", result[1]))
            rows.append(
                {
                    "feature": feature,
                    "program": program,
                    "rho": rho,
                    "p_value": p_value,
                }
            )

    associations = pd.DataFrame(rows, columns=["feature", "program", "rho", "p_value"])
    associations["fdr"] = fdr_bh(associations["p_value"].to_numpy())
    return associations[["feature", "program", "rho", "p_value", "fdr"]]


def correlation_matrix(associations: pd.DataFrame, value_col: str = "rho") -> pd.DataFrame:
    _require_columns(associations, ["feature", "program", value_col], "associations")
    return associations.pivot(index="feature", columns="program", values=value_col)


def top_associations(associations: pd.DataFrame, n: int = 10, by_abs: bool = True) -> pd.DataFrame:
    _require_columns(associations, ["rho"], "associations")
    if n < 1:
        raise ValueError("n must be at least 1")

    ranked = associations.copy()
    sort_col = "_abs_rho" if by_abs else "rho"
    if by_abs:
        ranked[sort_col] = ranked["rho"].abs()
    return ranked.sort_values(sort_col, ascending=False).drop(
        columns=[sort_col], errors="ignore"
    ).head(n)
