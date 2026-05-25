from pathlib import Path

import pandas as pd

from mumoas.config import load_config
from mumoas.io import ensure_directory, read_table, write_table


def test_load_public_config():
    cfg = load_config("configs/public_example.yaml")

    assert cfg.project.name == "MuMoAS public Python model workflow"
    assert cfg.columns.patient_id == "patient_id"
    assert cfg.columns.rcb_class == "rcb_class"
    assert cfg.imaging["manifest_path_columns"]["imaging_manifest"] == ["image_path", "mask_path"]
    assert cfg.radiogenomics["program_columns"] == ["MP1", "MP2", "MP3", "MP4", "MP5", "MP6"]
    assert cfg.single_cell["cell_id_column"] == "cell_id"
    assert cfg.multiomics["ranking_positive_group"] == "RCB-III"


def test_table_roundtrip(tmp_path: Path):
    output_dir = ensure_directory(tmp_path / "tables")
    table_path = output_dir / "one_row.csv"
    expected = pd.DataFrame({"patient_id": ["P001"], "value": [1.25]})

    write_table(expected, table_path)
    observed = read_table(table_path)

    pd.testing.assert_frame_equal(observed, expected)
