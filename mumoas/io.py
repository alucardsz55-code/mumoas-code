from __future__ import annotations

import random
from pathlib import Path

import numpy as np
import pandas as pd


def ensure_directory(path: str | Path) -> Path:
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def read_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    if table_path.suffix == ".csv":
        return pd.read_csv(table_path)
    if table_path.suffix == ".tsv":
        return pd.read_csv(table_path, sep="\t")
    raise ValueError(f"Unsupported table suffix '{table_path.suffix}' for {table_path}")


def write_table(table: pd.DataFrame, path: str | Path) -> Path:
    table_path = Path(path)
    ensure_directory(table_path.parent)
    if table_path.suffix == ".csv":
        table.to_csv(table_path, index=False)
        return table_path
    if table_path.suffix == ".tsv":
        table.to_csv(table_path, sep="\t", index=False)
        return table_path
    raise ValueError(f"Unsupported table suffix '{table_path.suffix}' for {table_path}")


def set_reproducible_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
