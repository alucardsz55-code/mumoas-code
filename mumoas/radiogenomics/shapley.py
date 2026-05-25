from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd


def compute_shap_if_available(
    model: Any,
    x: Any,
    feature_names: Sequence[str],
) -> pd.DataFrame | dict[str, str]:
    """Compute SHAP values when the optional dependency is installed."""
    try:
        import shap
    except ImportError:
        return {
            "status": "unavailable",
            "reason": "Optional dependency 'shap' is not installed.",
        }

    explainer = shap.Explainer(model, x)
    values = explainer(x)
    shap_values = getattr(values, "values", values)
    array = np.asarray(shap_values)
    if array.ndim > 2:
        array = array.reshape(array.shape[0], -1)

    columns = list(feature_names)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.shape[1] != len(columns):
        columns = [f"shap_{idx}" for idx in range(array.shape[1])]
    return pd.DataFrame(array, columns=columns)
