from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mumoas.config import class_label_order, load_config
from mumoas.io import read_table, write_table
from mumoas.modeling.features import build_design_matrix
from mumoas.modeling.train import export_metrics, predict_cvae, save_cvae_state, train_cvae


MAKE_SYNTHETIC_COMMAND = "python scripts/00_make_synthetic_data.py --config configs/public_example.yaml"
EXTRACT_IMAGING_COMMAND = "python scripts/01_extract_imaging_features.py --config configs/public_example.yaml"


def _validated_imaging_path(cfg) -> Path:
    return cfg.project.output_dir / "imaging" / "imaging_features_validated.csv"


def _missing_training_inputs_message(missing: list[Path], configured_imaging: Path | None = None) -> str:
    checked = [str(path) for path in missing]
    if configured_imaging is not None:
        checked.append(str(configured_imaging))
    return (
        "Missing cVAE training input(s): "
        f"{checked}. Generate the public synthetic inputs with `{MAKE_SYNTHETIC_COMMAND}`, "
        f"then validate imaging features with `{EXTRACT_IMAGING_COMMAND}` before training."
    )


def _imaging_features_path(cfg) -> Path:
    validated = _validated_imaging_path(cfg)
    if validated.exists():
        return validated
    configured = cfg.paths["imaging_features"]
    if configured.exists():
        return configured
    raise FileNotFoundError(_missing_training_inputs_message([validated], configured_imaging=configured))


def _read_required_table(path: Path):
    if not path.exists():
        raise FileNotFoundError(_missing_training_inputs_message([path]))
    return read_table(path)


def run(config_path: str | Path) -> dict[str, Path]:
    cfg = load_config(config_path)
    clinical = _read_required_table(cfg.paths["clinical"])
    imaging = read_table(_imaging_features_path(cfg))

    modeling_cfg = cfg.modeling
    cvae_cfg = modeling_cfg.get("cvae", {})
    design = build_design_matrix(
        clinical=clinical,
        imaging=imaging,
        patient_col=cfg.columns.patient_id,
        label_col=modeling_cfg.get("target", cfg.columns.rcb_class),
        categorical_cols=modeling_cfg.get("categorical_columns", cfg.columns.categorical),
        continuous_cols=modeling_cfg.get("continuous_columns", cfg.columns.continuous),
        label_order=class_label_order(cfg),
    )

    model, metrics = train_cvae(
        features=design.features,
        conditions=design.conditions,
        labels=design.labels,
        latent_dim=int(cvae_cfg.get("latent_dim", 8)),
        hidden_dims=tuple(cvae_cfg.get("hidden_dims", (32, 16))),
        learning_rate=float(cvae_cfg.get("learning_rate", 1e-3)),
        batch_size=int(cvae_cfg.get("batch_size", 16)),
        epochs=int(cvae_cfg.get("epochs", 5)),
        validation_fraction=float(cvae_cfg.get("validation_fraction", 0.25)),
        seed=int(cfg.project.seed),
        num_classes=len(design.metadata["label_classes"]),
    )
    metrics["design"] = design.metadata

    output_dir = cfg.project.output_dir / "modeling"
    predictions = predict_cvae(
        model,
        design.features,
        design.conditions,
        patient_ids=design.patient_ids,
        patient_col=cfg.columns.patient_id,
        label_names=design.metadata["label_classes"],
    )

    predictions_path = output_dir / "cvae_predictions.csv"
    metrics_path = output_dir / "cvae_metrics.json"
    state_path = output_dir / "cvae_state.pt"
    write_table(predictions, predictions_path)
    export_metrics(metrics, metrics_path)
    save_cvae_state(
        model,
        state_path,
        metadata={
            "feature_columns": design.metadata["feature_columns"],
            "condition_columns": design.metadata["condition_columns"],
            "label_classes": design.metadata["label_classes"],
            "cvae": cvae_cfg,
        },
    )

    print(f"Wrote cVAE predictions: {predictions_path}")
    print(f"Wrote cVAE metrics: {metrics_path}")
    print(f"Wrote cVAE state: {state_path}")
    return {
        "predictions": predictions_path,
        "metrics": metrics_path,
        "state": state_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train the public MuMoAS cVAE RCB model.")
    parser.add_argument("--config", required=True, help="Path to a MuMoAS YAML configuration file.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        run(args.config)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(
            f"{exc}\nFor the public example, run `{MAKE_SYNTHETIC_COMMAND}` and "
            f"`{EXTRACT_IMAGING_COMMAND}` first.",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
