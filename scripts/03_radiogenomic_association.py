from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mumoas.config import load_config
from mumoas.io import read_table, write_table
from mumoas.radiogenomics.association import correlation_matrix, spearman_program_associations


MAKE_SYNTHETIC_COMMAND = "python scripts/00_make_synthetic_data.py --config configs/public_example.yaml"
EXTRACT_IMAGING_COMMAND = "python scripts/01_extract_imaging_features.py --config configs/public_example.yaml"


def _missing_inputs_message(missing: list[Path]) -> str:
    return (
        "Missing radiogenomic input(s): "
        f"{[str(path) for path in missing]}. Generate public synthetic inputs with "
        f"`{MAKE_SYNTHETIC_COMMAND}` and validate imaging features with `{EXTRACT_IMAGING_COMMAND}` "
        "before running radiogenomic associations."
    )


def _imaging_features_path(cfg) -> Path:
    validated_path = cfg.project.output_dir / "imaging" / "imaging_features_validated.csv"
    if validated_path.exists():
        return validated_path
    configured_path = cfg.paths["imaging_features"]
    if configured_path.exists():
        return configured_path
    raise FileNotFoundError(_missing_inputs_message([validated_path, configured_path]))


def _read_required_table(path: Path):
    if not path.exists():
        raise FileNotFoundError(_missing_inputs_message([path]))
    return read_table(path)


def _program_columns(cfg) -> list[str]:
    configured = cfg.radiogenomics.get("program_columns")
    if configured:
        return list(configured)
    return list(cfg.signatures)


def run(config_path: str | Path) -> dict[str, Path]:
    cfg = load_config(config_path)
    imaging = read_table(_imaging_features_path(cfg))
    programs = _read_required_table(cfg.paths["mp_scores"])
    program_columns = _program_columns(cfg)

    associations = spearman_program_associations(
        imaging=imaging,
        programs=programs,
        patient_col=cfg.columns.patient_id,
        program_cols=program_columns,
    )
    matrix = correlation_matrix(associations)

    output_dir = cfg.project.output_dir / "radiogenomics"
    associations_path = output_dir / "mp_feature_associations.csv"
    matrix_path = output_dir / "mp_feature_correlation_matrix.csv"
    write_table(associations, associations_path)
    write_table(matrix.reset_index(), matrix_path)

    print(f"Wrote radiogenomic associations: {associations_path}")
    print(f"Wrote radiogenomic correlation matrix: {matrix_path}")
    return {
        "associations": associations_path,
        "correlation_matrix": matrix_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute public MuMoAS radiogenomic feature-program associations."
    )
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
