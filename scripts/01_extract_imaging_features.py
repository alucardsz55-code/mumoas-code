from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mumoas.config import load_config
from mumoas.io import read_table, write_table
from mumoas.imaging.radiomics import validate_feature_table


MAKE_SYNTHETIC_COMMAND = "python scripts/00_make_synthetic_data.py --config configs/public_example.yaml"
MANIFEST_PATH_KEYS = ("image_manifest", "mask_manifest", "imaging_manifest")
DEFAULT_PATH_COLUMNS_BY_KEY = {
    "image_manifest": ("image_path",),
    "mask_manifest": ("mask_path",),
    "imaging_manifest": ("image_path", "mask_path"),
}


def _missing_imaging_features_message(path: Path) -> str:
    return (
        f"Missing configured imaging feature table: {path}. "
        f"Generate the public synthetic inputs first with `{MAKE_SYNTHETIC_COMMAND}`."
    )


def _read_required_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(_missing_imaging_features_message(path))
    return read_table(path)


def _validate_manifest(path: Path, required_columns: tuple[str, ...]) -> pd.DataFrame:
    manifest = read_table(path)
    missing = [column for column in required_columns if column not in manifest.columns]
    if missing:
        raise ValueError(f"Manifest {path} is missing required path columns: {missing}")
    print(
        f"Validated manifest {path} columns {list(required_columns)}. "
        "Private image processing requires controlled local data and is not run by this public script."
    )
    return manifest


def _manifest_path_columns(cfg) -> dict[str, tuple[str, ...]]:
    configured = cfg.imaging.get("manifest_path_columns", {})
    columns_by_key = {}
    for key in MANIFEST_PATH_KEYS:
        columns_by_key[key] = tuple(configured.get(key, DEFAULT_PATH_COLUMNS_BY_KEY[key]))
    return columns_by_key


def run(config_path: str | Path) -> Path:
    cfg = load_config(config_path)
    imaging_features_path = cfg.paths["imaging_features"]
    features = _read_required_table(imaging_features_path)
    validated = validate_feature_table(features, patient_col=cfg.columns.patient_id)

    output_path = cfg.project.output_dir / "imaging" / "imaging_features_validated.csv"
    write_table(validated, output_path)

    for key, required_columns in _manifest_path_columns(cfg).items():
        manifest_path = cfg.paths.get(key)
        if manifest_path is not None and manifest_path.exists():
            _validate_manifest(manifest_path, required_columns)

    print(f"Wrote validated imaging features: {output_path}")
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate public MuMoAS imaging feature tables."
    )
    parser.add_argument("--config", required=True, help="Path to a MuMoAS YAML configuration file.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        run(args.config)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(f"{exc}\nFor the public example, run `{MAKE_SYNTHETIC_COMMAND}` first.", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
