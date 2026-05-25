from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset

from mumoas.io import ensure_directory
from mumoas.modeling.cvae import ConditionalVAE, cvae_loss


def split_train_validation(
    labels: np.ndarray,
    validation_fraction: float = 0.25,
    seed: int = 20260525,
) -> tuple[np.ndarray, np.ndarray]:
    indices = np.arange(len(labels))
    if len(indices) < 2 or validation_fraction <= 0:
        return indices, np.array([], dtype=int)

    labels = np.asarray(labels)
    class_counts = np.bincount(labels) if np.issubdtype(labels.dtype, np.integer) else np.array([])
    stratify = labels if class_counts.size and np.all(class_counts[class_counts > 0] >= 2) else None

    try:
        train_idx, val_idx = train_test_split(
            indices,
            test_size=validation_fraction,
            random_state=seed,
            stratify=stratify,
        )
    except ValueError:
        train_idx, val_idx = train_test_split(
            indices,
            test_size=validation_fraction,
            random_state=seed,
            stratify=None,
        )
    return np.sort(train_idx), np.sort(val_idx)


def _to_tensor(values: pd.DataFrame | np.ndarray, dtype: torch.dtype = torch.float32) -> torch.Tensor:
    if isinstance(values, pd.DataFrame):
        array = values.to_numpy(dtype=np.float32)
    else:
        array = np.asarray(values, dtype=np.float32)
    return torch.as_tensor(array, dtype=dtype)


def _make_loader(
    x: torch.Tensor,
    condition: torch.Tensor,
    labels: torch.Tensor,
    indices: np.ndarray,
    batch_size: int,
    seed: int,
    shuffle: bool,
) -> DataLoader:
    dataset = TensorDataset(x[indices], condition[indices], labels[indices])
    generator = torch.Generator()
    generator.manual_seed(seed)
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, generator=generator)


def train_cvae(
    features: pd.DataFrame | np.ndarray,
    conditions: pd.DataFrame | np.ndarray,
    labels: np.ndarray,
    latent_dim: int = 8,
    hidden_dims: tuple[int, ...] | list[int] = (32, 16),
    learning_rate: float = 1e-3,
    batch_size: int = 16,
    epochs: int = 5,
    validation_fraction: float = 0.25,
    seed: int = 20260525,
    beta: float = 1.0,
    class_weight: float = 1.0,
    classification_weight: float | None = None,
    num_classes: int | None = None,
) -> tuple[ConditionalVAE, dict[str, Any]]:
    if classification_weight is not None:
        class_weight = classification_weight
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(1)

    x = _to_tensor(features)
    condition = _to_tensor(conditions)
    y = torch.as_tensor(np.asarray(labels), dtype=torch.long)
    inferred_classes = int(y.max().item()) + 1
    if num_classes is None:
        num_classes = inferred_classes
    if num_classes < inferred_classes:
        raise ValueError(
            f"num_classes={num_classes} is smaller than the encoded labels require ({inferred_classes})"
        )

    model = ConditionalVAE(
        input_dim=x.shape[1],
        condition_dim=condition.shape[1],
        latent_dim=latent_dim,
        num_classes=num_classes,
        hidden_dims=hidden_dims,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    train_idx, val_idx = split_train_validation(labels, validation_fraction, seed)
    train_loader = _make_loader(x, condition, y, train_idx, batch_size, seed, shuffle=True)

    history: list[dict[str, float]] = []
    for epoch in range(int(epochs)):
        model.train()
        epoch_losses: list[float] = []
        for batch_x, batch_condition, batch_y in train_loader:
            optimizer.zero_grad()
            output = model(batch_x, batch_condition)
            losses = cvae_loss(output, batch_x, batch_y, beta=beta, class_weight=class_weight)
            losses.total.backward()
            optimizer.step()
            epoch_losses.append(float(losses.total.detach().item()))

        row = {"epoch": float(epoch + 1), "train_total_loss": float(np.mean(epoch_losses))}
        if val_idx.size:
            val_metrics = evaluate_cvae(model, x, condition, y, val_idx, beta, class_weight)
            row.update({f"validation_{key}": value for key, value in val_metrics.items()})
        history.append(row)

    metrics = {
        "history": history,
        "train_indices": train_idx.tolist(),
        "validation_indices": val_idx.tolist(),
    }
    if val_idx.size:
        metrics["validation"] = evaluate_cvae(model, x, condition, y, val_idx, beta, class_weight)
    metrics["training"] = evaluate_cvae(model, x, condition, y, train_idx, beta, class_weight)
    return model, metrics


def evaluate_cvae(
    model: ConditionalVAE,
    x: torch.Tensor,
    condition: torch.Tensor,
    labels: torch.Tensor,
    indices: np.ndarray,
    beta: float = 1.0,
    class_weight: float = 1.0,
    classification_weight: float | None = None,
) -> dict[str, float]:
    if classification_weight is not None:
        class_weight = classification_weight
    if len(indices) == 0:
        return {}
    model.eval()
    with torch.no_grad():
        output = model(x[indices], condition[indices])
        losses = cvae_loss(output, x[indices], labels[indices], beta=beta, class_weight=class_weight)
        predicted = torch.argmax(output.logits, dim=1).cpu().numpy()
        observed = labels[indices].cpu().numpy()
    return {
        "total_loss": float(losses.total.item()),
        "reconstruction_loss": float(losses.reconstruction.item()),
        "kl_loss": float(losses.kl.item()),
        "classification_loss": float(losses.classification.item()),
        "accuracy": float(accuracy_score(observed, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(observed, predicted)),
    }


def predict_cvae(
    model: ConditionalVAE,
    features: pd.DataFrame | np.ndarray,
    conditions: pd.DataFrame | np.ndarray,
    patient_ids: pd.Series | list[str] | np.ndarray | None = None,
    patient_col: str = "patient_id",
    label_names: list[str] | None = None,
) -> pd.DataFrame:
    model.eval()
    x = _to_tensor(features)
    condition = _to_tensor(conditions)
    with torch.no_grad():
        output = model(x, condition)
        probabilities = torch.softmax(output.logits, dim=1).cpu().numpy()
        predicted_index = probabilities.argmax(axis=1)

    if label_names is None:
        label_names = [str(idx) for idx in range(probabilities.shape[1])]
    table = pd.DataFrame(
        probabilities,
        columns=[f"prob_{label}" for label in label_names],
    )
    table.insert(0, "predicted_class", [label_names[idx] for idx in predicted_index])
    if patient_ids is not None:
        table.insert(0, patient_col, list(patient_ids))
    return table


def export_metrics(metrics: dict[str, Any], path: str | Path) -> Path:
    output_path = Path(path)
    ensure_directory(output_path.parent)
    output_path.write_text(json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8")
    return output_path


def save_cvae_state(
    model: ConditionalVAE,
    path: str | Path,
    metadata: dict[str, Any] | None = None,
) -> Path:
    output_path = Path(path)
    ensure_directory(output_path.parent)
    torch.save({"state_dict": model.state_dict(), "metadata": metadata or {}}, output_path)
    return output_path
