import numpy as np
import pandas as pd
import torch

from mumoas.imaging.foundation import CompactMambaTransformerEncoder
from mumoas.imaging.patches import extract_center_patches
from mumoas.imaging.radiomics import validate_feature_table


def test_extract_center_patches_shape():
    volume = np.ones((8, 32, 32), dtype=np.float32)
    mask = np.zeros_like(volume)
    mask[:, 8:24, 8:24] = 1

    patches = extract_center_patches(volume, mask, patch_size=16, max_patches=4)

    assert patches.shape == (4, 1, 16, 16)


def test_validate_feature_table_keeps_patient_and_numeric_features():
    table = pd.DataFrame(
        {"patient_id": ["P001"], "feat_a": [0.5], "text_metadata": ["drop"]}
    )

    features = validate_feature_table(table, patient_col="patient_id")

    assert list(features.columns) == ["patient_id", "feat_a"]


def test_compact_encoder_forward_shape():
    model = CompactMambaTransformerEncoder(in_channels=3, embed_dim=16, output_dim=8)
    x = torch.randn(2, 3, 32, 32)

    y = model(x)

    assert y.shape == (2, 8)
