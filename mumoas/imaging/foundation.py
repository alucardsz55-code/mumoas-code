from __future__ import annotations

import torch
from torch import nn


class GatedStateMixingBlock(nn.Module):
    """Small gated token-mixing block used as a public Mamba-like scaffold."""

    def __init__(self, embed_dim: int) -> None:
        super().__init__()
        self.norm = nn.LayerNorm(embed_dim)
        self.mix = nn.Conv1d(embed_dim, embed_dim, kernel_size=3, padding=1, groups=embed_dim)
        self.gate = nn.Linear(embed_dim, embed_dim)
        self.value = nn.Linear(embed_dim, embed_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x
        x_norm = self.norm(x)
        mixed = self.mix(x_norm.transpose(1, 2)).transpose(1, 2)
        gated = torch.sigmoid(self.gate(x_norm)) * self.value(mixed)
        return residual + gated


class CompactMambaTransformerEncoder(nn.Module):
    """Lightweight public imaging encoder scaffold.

    Manuscript training data and weights are not distributed in this public repository.
    This module provides a callable architecture for reviewer inspection and smoke tests.
    """

    def __init__(
        self,
        in_channels: int = 1,
        embed_dim: int = 32,
        output_dim: int = 16,
        patch_size: int = 4,
        depth: int = 2,
        num_heads: int = 4,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        if embed_dim % num_heads != 0:
            raise ValueError("embed_dim must be divisible by num_heads")

        self.patch_embed = nn.Conv2d(
            in_channels, embed_dim, kernel_size=patch_size, stride=patch_size
        )
        self.state_blocks = nn.ModuleList(
            GatedStateMixingBlock(embed_dim) for _ in range(depth)
        )
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embed_dim,
            nhead=num_heads,
            dim_feedforward=embed_dim * 4,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
            norm_first=True,
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=depth)
        self.norm = nn.LayerNorm(embed_dim)
        self.projection = nn.Linear(embed_dim, output_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        tokens = self.patch_embed(x).flatten(2).transpose(1, 2)
        for block in self.state_blocks:
            tokens = block(tokens)
        encoded = self.transformer(tokens)
        pooled = self.norm(encoded).mean(dim=1)
        return self.projection(pooled)
