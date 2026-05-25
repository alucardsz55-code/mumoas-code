import numpy as np
import pandas as pd
import torch

from mumoas.modeling.cvae import ConditionalVAE, cvae_loss
from mumoas.modeling.features import build_design_matrix, prune_correlated_features, prune_vif
from mumoas.modeling.train import train_cvae


def test_prune_correlated_features_drops_duplicate_signal():
    features = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5],
            "b": [1, 2, 3, 4, 5],
            "c": [5, 1, 4, 2, 3],
        }
    )

    pruned, kept, dropped = prune_correlated_features(features, threshold=0.9)

    assert list(pruned.columns) == ["a", "c"]
    assert kept == ["a", "c"]
    assert dropped == ["b"]


def test_prune_vif_returns_subset_for_identity_matrix():
    features = pd.DataFrame(np.eye(4), columns=["a", "b", "c", "d"])

    pruned, kept, dropped = prune_vif(features, threshold=10.0)

    assert set(kept).issubset(features.columns)
    assert list(pruned.columns) == kept
    assert dropped == []


def test_conditional_vae_forward_and_loss_are_well_shaped():
    torch.manual_seed(20260525)
    model = ConditionalVAE(input_dim=6, condition_dim=3, latent_dim=4, num_classes=3)
    x = torch.randn(5, 6)
    condition = torch.eye(3)[torch.tensor([0, 1, 2, 0, 1])].float()
    labels = torch.tensor([0, 1, 2, 0, 1])

    output = model(x, condition)
    loss = cvae_loss(output, x, labels, beta=0.1, class_weight=1.0)

    assert output.logits.shape == (5, 3)
    assert loss.total.item() > 0
    assert loss.reconstruction.item() >= 0
    assert loss.kl.item() >= 0
    assert loss.classification.item() >= 0


def test_build_design_matrix_uses_configured_label_order():
    clinical = pd.DataFrame(
        {
            "patient_id": ["P1", "P2", "P3"],
            "rcb_class": ["RCB-III", "RCB-0", "RCB-II"],
            "age": [50, 45, 60],
            "cohort": ["A", "A", "B"],
        }
    )
    imaging = pd.DataFrame(
        {
            "patient_id": ["P1", "P2", "P3"],
            "feat_a": [0.1, 0.2, 0.3],
        }
    )

    design = build_design_matrix(
        clinical=clinical,
        imaging=imaging,
        patient_col="patient_id",
        label_col="rcb_class",
        categorical_cols=["cohort"],
        continuous_cols=["age"],
        label_order=["RCB-0", "RCB-I", "RCB-II", "RCB-III"],
    )

    assert design.metadata["label_classes"] == ["RCB-0", "RCB-I", "RCB-II", "RCB-III"]
    assert design.labels.tolist() == [3, 0, 2]


def test_train_cvae_keeps_configured_class_count_when_class_absent():
    features = pd.DataFrame(np.random.default_rng(1).normal(size=(6, 3)))
    conditions = pd.DataFrame(np.random.default_rng(2).normal(size=(6, 2)))
    labels = np.array([0, 1, 1, 0, 1, 0])

    model, _ = train_cvae(
        features=features,
        conditions=conditions,
        labels=labels,
        latent_dim=2,
        hidden_dims=(4,),
        epochs=1,
        batch_size=3,
        validation_fraction=0.0,
        seed=1,
        num_classes=4,
    )

    assert model.num_classes == 4
