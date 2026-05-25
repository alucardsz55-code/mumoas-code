from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mumoas.config import class_label_order, load_config
from mumoas.evaluation import (
    calibration_table,
    confusion_matrix_table,
    multiclass_auc_table,
    net_benefit_curve,
    prepare_survival_table,
)
from mumoas.evaluation.classification import decode_prediction_labels, probability_columns
from mumoas.io import read_table, write_table


TRAIN_FIRST_MESSAGE = (
    "Missing prediction input. Run "
    "`python scripts/02_train_cvae_rcb.py --config configs/public_example.yaml` first, "
    "or set a prediction path in the config."
)


def _configured_prediction_path(cfg) -> Path | None:
    for key in ("predictions", "prediction", "cvae_predictions"):
        path = cfg.paths.get(key)
        if path is not None:
            return Path(path)
    modeling_path = cfg.modeling.get("predictions_path") or cfg.modeling.get("prediction_path")
    return Path(modeling_path) if modeling_path else None


def _resolve_prediction_path(cfg) -> Path:
    default_path = cfg.project.output_dir / "modeling" / "cvae_predictions.csv"
    if default_path.exists():
        return default_path
    configured = _configured_prediction_path(cfg)
    if configured is not None and configured.exists():
        return configured
    checked = [str(default_path)]
    if configured is not None:
        checked.append(str(configured))
    raise FileNotFoundError(f"{TRAIN_FIRST_MESSAGE} Checked: {', '.join(checked)}")


def _class_order(cfg, predictions: pd.DataFrame) -> list[str]:
    configured_order = cfg.modeling.get("class_label_order") or cfg.columns.class_label_order
    if configured_order:
        classes = list(configured_order)
        missing = [f"prob_{label}" for label in classes if f"prob_{label}" not in predictions.columns]
        if missing:
            raise ValueError(
                "Configured class_label_order requires matching prediction probability columns; "
                f"missing: {missing}"
            )
        return classes
    return [column.removeprefix("prob_") for column in probability_columns(predictions)]


def _merge_truth(cfg, predictions: pd.DataFrame) -> pd.DataFrame:
    patient_col = cfg.columns.patient_id
    label_col = cfg.modeling.get("target", cfg.columns.rcb_class)
    if label_col in predictions.columns:
        return predictions.copy()
    if patient_col not in predictions.columns:
        raise ValueError(f"Predictions must contain '{patient_col}' when true labels are read from clinical data.")
    clinical = read_table(cfg.paths["clinical"])
    required = [patient_col, label_col]
    missing = [column for column in required if column not in clinical.columns]
    if missing:
        raise ValueError(f"Clinical table is missing required columns for evaluation: {missing}")
    return predictions.merge(clinical[required], on=patient_col, how="left", validate="one_to_one")


def _binary_target_for_curves(table: pd.DataFrame, label_col: str, classes: list[str]) -> tuple[pd.Series, pd.Series]:
    predicted = table["predicted_class"] if "predicted_class" in table.columns else decode_prediction_labels(table, classes)
    confidence = table[[f"prob_{label}" for label in classes]].max(axis=1)
    correctness = (table[label_col] == predicted).astype(int)
    return correctness, confidence


def _decision_thresholds() -> np.ndarray:
    return np.round(np.linspace(0.05, 0.95, 19), 2)


def run(config_path: str | Path) -> dict[str, Path]:
    cfg = load_config(config_path)
    prediction_path = _resolve_prediction_path(cfg)
    predictions = read_table(prediction_path)
    classes = _class_order(cfg, predictions)
    if not classes:
        raise ValueError("Predictions must contain probability columns named like 'prob_<class>'.")

    table = _merge_truth(cfg, predictions)
    label_col = cfg.modeling.get("target", cfg.columns.rcb_class)
    if table[label_col].isna().any():
        raise ValueError("Some predictions could not be matched to true labels for evaluation.")

    if "predicted_class" not in table.columns:
        table["predicted_class"] = decode_prediction_labels(table, classes)

    output_dir = cfg.project.output_dir / "evaluation"
    outputs = {
        "multiclass_auc": output_dir / "multiclass_auc.csv",
        "confusion_matrix": output_dir / "confusion_matrix.csv",
        "calibration": output_dir / "calibration.csv",
        "decision_curve": output_dir / "decision_curve.csv",
    }

    write_table(multiclass_auc_table(table[label_col], table, classes), outputs["multiclass_auc"])
    write_table(confusion_matrix_table(table[label_col], table["predicted_class"], classes), outputs["confusion_matrix"])

    curve_true, curve_score = _binary_target_for_curves(table, label_col, classes)
    write_table(calibration_table(curve_true, curve_score, n_bins=10), outputs["calibration"])
    write_table(net_benefit_curve(curve_true, curve_score, thresholds=_decision_thresholds()), outputs["decision_curve"])

    if cfg.columns.efs_time and cfg.columns.efs_event:
        clinical = read_table(cfg.paths["clinical"])
        survival_output = output_dir / "survival_input.csv"
        try:
            survival_table = prepare_survival_table(
                clinical=clinical,
                predictions=table[[cfg.columns.patient_id, "predicted_class"]],
                patient_col=cfg.columns.patient_id,
                time_col=cfg.columns.efs_time,
                event_col=cfg.columns.efs_event,
            )
        except ValueError:
            survival_table = None
        if survival_table is not None:
            write_table(survival_table, survival_output)
            outputs["survival_input"] = survival_output

    for name, path in outputs.items():
        print(f"Wrote {name}: {path}")
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate public MuMoAS cVAE prediction exports.")
    parser.add_argument("--config", required=True, help="Path to a MuMoAS YAML configuration file.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        run(args.config)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
