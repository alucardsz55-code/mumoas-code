from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn
from torch.nn import functional as F


@dataclass
class ConditionalVAEOutput:
    reconstruction: torch.Tensor
    mu: torch.Tensor
    logvar: torch.Tensor
    z: torch.Tensor
    logits: torch.Tensor


@dataclass
class ConditionalVAELoss:
    total: torch.Tensor
    reconstruction: torch.Tensor
    kl: torch.Tensor
    classification: torch.Tensor


class ConditionalVAE(nn.Module):
    """Small conditional variational autoencoder with an auxiliary classifier."""

    def __init__(
        self,
        input_dim: int,
        condition_dim: int,
        latent_dim: int,
        num_classes: int,
        hidden_dims: tuple[int, ...] | list[int] = (64, 32),
    ) -> None:
        super().__init__()
        if input_dim <= 0 or latent_dim <= 0 or num_classes <= 1:
            raise ValueError("input_dim, latent_dim, and num_classes must be positive")
        if condition_dim < 0:
            raise ValueError("condition_dim cannot be negative")

        self.input_dim = input_dim
        self.condition_dim = condition_dim
        self.latent_dim = latent_dim
        self.num_classes = num_classes

        encoder_dims = [input_dim + condition_dim, *hidden_dims]
        self.encoder = _make_mlp(encoder_dims)
        encoder_out = encoder_dims[-1]
        self.fc_mu = nn.Linear(encoder_out, latent_dim)
        self.fc_logvar = nn.Linear(encoder_out, latent_dim)

        decoder_dims = [latent_dim + condition_dim, *reversed(hidden_dims), input_dim]
        self.decoder = _make_mlp(decoder_dims, final_activation=False)
        self.classifier = nn.Linear(latent_dim + condition_dim, num_classes)

    def encode(self, x: torch.Tensor, condition: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        encoded = self.encoder(torch.cat([x, condition], dim=1))
        return self.fc_mu(encoded), self.fc_logvar(encoded)

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        if not self.training:
            return mu
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z: torch.Tensor, condition: torch.Tensor) -> torch.Tensor:
        return self.decoder(torch.cat([z, condition], dim=1))

    def forward(self, x: torch.Tensor, condition: torch.Tensor) -> ConditionalVAEOutput:
        mu, logvar = self.encode(x, condition)
        z = self.reparameterize(mu, logvar)
        reconstruction = self.decode(z, condition)
        logits = self.classifier(torch.cat([z, condition], dim=1))
        return ConditionalVAEOutput(
            reconstruction=reconstruction,
            mu=mu,
            logvar=logvar,
            z=z,
            logits=logits,
        )


def _make_mlp(dims: list[int], final_activation: bool = True) -> nn.Sequential:
    layers: list[nn.Module] = []
    for idx in range(len(dims) - 1):
        layers.append(nn.Linear(dims[idx], dims[idx + 1]))
        if idx < len(dims) - 2 or final_activation:
            layers.append(nn.ReLU())
    return nn.Sequential(*layers)


def cvae_loss(
    output: ConditionalVAEOutput,
    target: torch.Tensor,
    labels: torch.Tensor,
    beta: float = 1.0,
    class_weight: float = 1.0,
    classification_weight: float | None = None,
) -> ConditionalVAELoss:
    if classification_weight is not None:
        class_weight = classification_weight
    reconstruction = F.mse_loss(output.reconstruction, target, reduction="mean")
    kl = -0.5 * torch.mean(1 + output.logvar - output.mu.pow(2) - output.logvar.exp())
    classification = F.cross_entropy(output.logits, labels)
    total = reconstruction + beta * kl + class_weight * classification
    return ConditionalVAELoss(
        total=total,
        reconstruction=reconstruction,
        kl=kl,
        classification=classification,
    )
