import importlib.util
from types import SimpleNamespace

import pandas as pd
import pytest

from mumoas.evaluation.calibration import calibration_table
from mumoas.evaluation.classification import multiclass_auc_table
from mumoas.evaluation.decision_curve import net_benefit_curve


def _load_evaluation_script():
    spec = importlib.util.spec_from_file_location(
        "evaluate_predictions",
        "scripts/06_evaluate_predictions.py",
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_multiclass_auc_table_returns_one_row_per_class_for_three_classes():
    y_true = pd.Series(["A", "B", "C", "A", "B", "C"])
    probabilities = pd.DataFrame(
        {
            "prob_A": [0.8, 0.1, 0.1, 0.7, 0.2, 0.1],
            "prob_B": [0.1, 0.7, 0.2, 0.2, 0.6, 0.2],
            "prob_C": [0.1, 0.2, 0.7, 0.1, 0.2, 0.7],
        }
    )

    auc = multiclass_auc_table(y_true, probabilities, classes=["A", "B", "C"])

    assert list(auc["class"]) == ["A", "B", "C"]
    assert len(auc) == 3


def test_configured_class_order_requires_all_probability_columns():
    module = _load_evaluation_script()
    cfg = SimpleNamespace(
        modeling={"class_label_order": ["A", "B", "C"]},
        columns=SimpleNamespace(class_label_order=[]),
    )
    predictions = pd.DataFrame({"prob_A": [0.7], "prob_C": [0.3]})

    with pytest.raises(ValueError, match="prob_B"):
        module._class_order(cfg, predictions)


def test_calibration_table_returns_expected_columns():
    table = calibration_table(
        y_true=[0, 1, 1, 0],
        y_pred=[0.1, 0.4, 0.8, 0.9],
        n_bins=2,
    )

    assert list(table.columns) == ["bin", "mean_predicted", "observed_rate", "n"]


def test_net_benefit_curve_returns_requested_thresholds():
    curve = net_benefit_curve(
        y_true=[0, 1, 1, 0],
        y_pred=[0.1, 0.4, 0.8, 0.9],
        thresholds=[0.2, 0.5],
    )

    assert sorted(curve["threshold"].unique().tolist()) == [0.2, 0.5]
