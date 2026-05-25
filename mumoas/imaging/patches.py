from __future__ import annotations

import numpy as np


def _center_crop_with_padding(image: np.ndarray, center_yx: tuple[float, float], size: int) -> np.ndarray:
    half = size // 2
    center_y = int(round(center_yx[0]))
    center_x = int(round(center_yx[1]))
    y0 = center_y - half
    x0 = center_x - half
    y1 = y0 + size
    x1 = x0 + size

    pad_top = max(0, -y0)
    pad_left = max(0, -x0)
    pad_bottom = max(0, y1 - image.shape[0])
    pad_right = max(0, x1 - image.shape[1])

    padded = np.pad(
        image,
        ((pad_top, pad_bottom), (pad_left, pad_right)),
        mode="constant",
        constant_values=0,
    )
    y0 += pad_top
    y1 += pad_top
    x0 += pad_left
    x1 += pad_left
    return padded[y0:y1, x0:x1]


def _deterministic_indices(n_items: int, max_items: int) -> list[int]:
    if n_items >= max_items:
        return np.linspace(0, n_items - 1, num=max_items, dtype=int).tolist()
    return [idx % n_items for idx in range(max_items)]


def extract_center_patches(
    volume: np.ndarray, mask: np.ndarray, patch_size: int, max_patches: int
) -> np.ndarray:
    """Extract deterministic 2D tumor-centered patches from a 3D image volume."""
    volume_array = np.asarray(volume)
    mask_array = np.asarray(mask)
    if volume_array.shape != mask_array.shape:
        raise ValueError(
            f"volume and mask must have identical shapes, got {volume_array.shape} and {mask_array.shape}"
        )
    if volume_array.ndim != 3:
        raise ValueError(f"volume and mask must be 3D arrays shaped (z, y, x), got ndim={volume_array.ndim}")
    if patch_size <= 0 or max_patches <= 0:
        raise ValueError("patch_size and max_patches must be positive integers")

    positive_slices = np.flatnonzero(mask_array.reshape(mask_array.shape[0], -1).any(axis=1))
    if positive_slices.size == 0:
        raise ValueError("mask does not contain any positive voxels")

    patches: list[np.ndarray] = []
    for z_index in positive_slices:
        slice_mask = mask_array[z_index] > 0
        y_coords, x_coords = np.nonzero(slice_mask)
        center = (float(y_coords.mean()), float(x_coords.mean()))
        patch = _center_crop_with_padding(volume_array[z_index], center, patch_size)
        patches.append(patch.astype(np.float32, copy=False))

    selected = [patches[idx] for idx in _deterministic_indices(len(patches), max_patches)]
    return np.stack(selected, axis=0)[:, np.newaxis, :, :].astype(np.float32, copy=False)
