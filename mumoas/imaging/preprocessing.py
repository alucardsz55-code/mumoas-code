from __future__ import annotations

from typing import Sequence

import numpy as np


def zscore_inside_mask(volume: np.ndarray, mask: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Z-score a volume using only mask-positive voxels for the summary statistics."""
    volume_array = np.asarray(volume, dtype=np.float32)
    mask_array = np.asarray(mask)
    if volume_array.shape != mask_array.shape:
        raise ValueError(
            f"volume and mask must have identical shapes, got {volume_array.shape} and {mask_array.shape}"
        )

    values = volume_array[mask_array > 0]
    if values.size == 0:
        raise ValueError("mask does not contain any positive voxels for z-score normalization")
    mean = float(values.mean())
    std = float(values.std())
    if std < eps:
        std = 1.0
    return ((volume_array - mean) / std).astype(np.float32, copy=False)


def crop_to_mask(volume: np.ndarray, mask: np.ndarray, margin: int | Sequence[int] = 0) -> tuple[np.ndarray, np.ndarray]:
    """Crop a volume and mask to the mask-positive bounding box."""
    volume_array = np.asarray(volume)
    mask_array = np.asarray(mask)
    if volume_array.shape != mask_array.shape:
        raise ValueError(
            f"volume and mask must have identical shapes, got {volume_array.shape} and {mask_array.shape}"
        )
    if volume_array.ndim != 3:
        raise ValueError(f"volume and mask must be 3D arrays shaped (z, y, x), got ndim={volume_array.ndim}")

    coords = np.argwhere(mask_array > 0)
    if coords.size == 0:
        raise ValueError("mask does not contain any positive voxels for cropping")

    if isinstance(margin, int):
        margins = (margin, margin, margin)
    else:
        margins = tuple(int(value) for value in margin)
    if len(margins) != 3 or any(value < 0 for value in margins):
        raise ValueError("margin must be a non-negative integer or a sequence of three non-negative integers")

    starts = np.maximum(coords.min(axis=0) - np.asarray(margins), 0)
    stops = np.minimum(coords.max(axis=0) + np.asarray(margins) + 1, volume_array.shape)
    slices = tuple(slice(int(start), int(stop)) for start, stop in zip(starts, stops))
    return volume_array[slices], mask_array[slices]


def resample_volume_interface(
    volume: np.ndarray,
    spacing: Sequence[float] | None = None,
    target_spacing: Sequence[float] | None = None,
) -> np.ndarray:
    """Public interface for controlled local image resampling backends."""
    volume_array = np.asarray(volume)
    if target_spacing is None or spacing is None or tuple(spacing) == tuple(target_spacing):
        return volume_array.copy()

    raise RuntimeError(
        "Real image resampling is not bundled with the public MuMoAS repository. "
        "Install and configure a controlled local imaging backend such as SimpleITK, "
        "MONAI, or scipy.ndimage, then adapt resample_volume_interface for your data governance environment."
    )
