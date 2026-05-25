from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mumoas.config import MuMoASConfig, class_label_order, load_config
from mumoas.io import set_reproducible_seed, write_table


CONFIG_HELP = "Check --config and run: python scripts/00_make_synthetic_data.py --config configs/public_example.yaml"


def _patient_frame(cfg: MuMoASConfig, n_per_class: int, rng: np.random.Generator) -> pd.DataFrame:
    classes = class_label_order(cfg)
    records = []
    for class_index, rcb_class in enumerate(classes):
        for offset in range(n_per_class):
            patient_number = class_index * n_per_class + offset + 1
            records.append(
                {
                    cfg.columns.patient_id: f"P{patient_number:03d}",
                    cfg.columns.cohort: "training" if offset < int(n_per_class * 0.7) else "validation",
                    cfg.columns.rcb_class: rcb_class,
                    "age": int(rng.integers(31, 76)),
                    "tumor_size_cm": round(float(rng.normal(2.6 + 0.35 * class_index, 0.55)), 2),
                    "ki67_percent": round(float(rng.normal(28 + 7 * class_index, 8)), 1),
                    "menopausal_status": rng.choice(["pre", "post"]),
                    "clinical_stage": rng.choice(["II", "III"], p=[0.65, 0.35]),
                    "hr_status": rng.choice(["positive", "negative"], p=[0.7, 0.3]),
                    "her2_status": rng.choice(["positive", "negative"], p=[0.25, 0.75]),
                    cfg.columns.efs_time or "efs_time": round(float(rng.gamma(5.0 - 0.4 * class_index, 6.0)), 1),
                    cfg.columns.efs_event or "efs_event": int(rng.random() < 0.08 + 0.08 * class_index),
                }
            )
    return pd.DataFrame.from_records(records)


def _imaging_features(clinical: pd.DataFrame, cfg: MuMoASConfig, rng: np.random.Generator) -> pd.DataFrame:
    class_rank = clinical[cfg.columns.rcb_class].map(
        {label: index for index, label in enumerate(class_label_order(cfg))}
    )
    table = clinical[[cfg.columns.patient_id]].copy()
    for feature_index in range(1, 13):
        noise = rng.normal(0, 1, len(clinical))
        table[f"rad_feature_{feature_index:02d}"] = np.round(
            noise + class_rank.to_numpy() * (0.12 * feature_index), 4
        )
    return table


def _mp_scores(clinical: pd.DataFrame, cfg: MuMoASConfig, rng: np.random.Generator) -> pd.DataFrame:
    class_rank = clinical[cfg.columns.rcb_class].map(
        {label: index for index, label in enumerate(class_label_order(cfg))}
    ).to_numpy()
    table = clinical[[cfg.columns.patient_id]].copy()
    for index, signature in enumerate(cfg.signatures, start=1):
        direction = -1 if signature in {"MP1", "MP4"} else 1
        table[signature] = np.round(rng.normal(0, 0.35, len(clinical)) + direction * class_rank * 0.35, 4)
    return table


def _signature_genes(cfg: MuMoASConfig) -> list[str]:
    genes: list[str] = []
    for members in cfg.signatures.values():
        for gene in members:
            if gene not in genes:
                genes.append(gene)
    return genes


def _omics_table(
    clinical: pd.DataFrame,
    cfg: MuMoASConfig,
    rng: np.random.Generator,
    scale: float,
) -> pd.DataFrame:
    genes = _signature_genes(cfg)
    class_rank = clinical[cfg.columns.rcb_class].map(
        {label: index for index, label in enumerate(class_label_order(cfg))}
    ).to_numpy()
    table = clinical[[cfg.columns.patient_id]].copy()
    for gene_index, gene in enumerate(genes):
        trend = (gene_index % 4 - 1.5) * 0.08 * class_rank
        table[gene] = np.round(rng.lognormal(mean=scale + trend, sigma=0.25), 4)
    return table


def _single_cell_tables(
    clinical: pd.DataFrame,
    cfg: MuMoASConfig,
    rng: np.random.Generator,
    cells_per_patient: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = _signature_genes(cfg)
    cell_id_col = cfg.single_cell.get("cell_id_column", "cell_id")
    sample_id_col = cfg.single_cell.get("sample_id_column", "sample_id")
    cell_type_col = cfg.single_cell.get("cell_type_column", "cell_type")
    cell_types = ["T cell", "B cell", "myeloid", "epithelial", "stromal"]
    expression_records = []
    obs_records = []
    for _, row in clinical.iterrows():
        for cell_offset in range(cells_per_patient):
            cell_id = f"{row[cfg.columns.patient_id]}_C{cell_offset + 1:03d}"
            cell_type = rng.choice(cell_types)
            expression = {
                gene: int(rng.poisson(2 + (gene in cfg.signatures.get("MP4", [])) * (cell_type == "T cell") * 3))
                for gene in genes
            }
            expression_records.append({cell_id_col: cell_id, **expression})
            obs_records.append(
                {
                    cell_id_col: cell_id,
                    cfg.columns.patient_id: row[cfg.columns.patient_id],
                    sample_id_col: f"S_{row[cfg.columns.patient_id]}",
                    cfg.columns.cohort: row[cfg.columns.cohort],
                    cfg.columns.rcb_class: row[cfg.columns.rcb_class],
                    cell_type_col: cell_type,
                }
            )
    return pd.DataFrame(expression_records), pd.DataFrame(obs_records)


def make_synthetic_data(cfg: MuMoASConfig, n_per_class: int = 20, cells_per_patient: int = 8) -> list[Path]:
    set_reproducible_seed(cfg.project.seed)
    rng = np.random.default_rng(cfg.project.seed)
    clinical = _patient_frame(cfg, n_per_class=n_per_class, rng=rng)

    outputs = {
        "clinical": clinical,
        "imaging_features": _imaging_features(clinical, cfg, rng),
        "mp_scores": _mp_scores(clinical, cfg, rng),
        "bulk_expression": _omics_table(clinical, cfg, rng, scale=2.0),
        "proteomics": _omics_table(clinical, cfg, rng, scale=1.2),
    }
    single_cell_expression, single_cell_obs = _single_cell_tables(clinical, cfg, rng, cells_per_patient)
    outputs["single_cell_expression"] = single_cell_expression
    outputs["single_cell_obs"] = single_cell_obs

    written = []
    for key, table in outputs.items():
        written.append(write_table(table, cfg.paths[key]))
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create deterministic MuMoAS public synthetic CSV inputs.")
    parser.add_argument("--config", default="configs/public_example.yaml", help="Path to the public YAML config.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        cfg = load_config(args.config)
        written = make_synthetic_data(cfg)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(f"{exc}\n{CONFIG_HELP}", file=sys.stderr)
        raise SystemExit(1) from exc
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
