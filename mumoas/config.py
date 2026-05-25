from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

DEFAULT_CLASS_LABEL_ORDER = ["RCB-0", "RCB-I", "RCB-II", "RCB-III"]


@dataclass(frozen=True)
class ProjectConfig:
    name: str
    seed: int
    output_dir: Path


@dataclass(frozen=True)
class ColumnConfig:
    patient_id: str
    cohort: str
    rcb_class: str
    efs_time: str | None = None
    efs_event: str | None = None
    categorical: list[str] = field(default_factory=list)
    continuous: list[str] = field(default_factory=list)
    class_label_order: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class MuMoASConfig:
    project: ProjectConfig
    columns: ColumnConfig
    paths: dict[str, Path]
    modeling: dict[str, Any]
    signatures: dict[str, list[str]]
    imaging: dict[str, Any] = field(default_factory=dict)
    radiogenomics: dict[str, Any] = field(default_factory=dict)
    single_cell: dict[str, Any] = field(default_factory=dict)
    multiomics: dict[str, Any] = field(default_factory=dict)


def _resolve_path(base_dir: Path, value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path.resolve()
    return (base_dir / path).resolve()


def class_label_order(cfg: MuMoASConfig) -> list[str]:
    modeling_order = cfg.modeling.get("class_label_order")
    if modeling_order:
        return list(modeling_order)
    if cfg.columns.class_label_order:
        return list(cfg.columns.class_label_order)
    return list(DEFAULT_CLASS_LABEL_ORDER)


def load_config(path: str | Path) -> MuMoASConfig:
    config_path = Path(path)
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Configuration file {config_path} must contain a YAML mapping")

    base_dir = config_path.parent.parent
    project_data = data["project"]
    project = ProjectConfig(
        name=project_data["name"],
        seed=int(project_data.get("seed", 20260525)),
        output_dir=_resolve_path(base_dir, project_data["output_dir"]),
    )

    paths = {key: _resolve_path(base_dir, value) for key, value in data["paths"].items()}
    columns = ColumnConfig(**data["columns"])
    signatures = data.get("signatures", {})
    if not isinstance(signatures, dict):
        raise ValueError("'signatures' must be a mapping of signature name to genes")
    modeling = dict(data.get("modeling") or {})
    if "response_composites" not in modeling and "response_composites" in data:
        modeling["response_composites"] = data["response_composites"]

    return MuMoASConfig(
        project=project,
        columns=columns,
        paths=paths,
        modeling=modeling,
        signatures={key: list(value) for key, value in signatures.items()},
        imaging=dict(data.get("imaging") or {}),
        radiogenomics=dict(data.get("radiogenomics") or {}),
        single_cell=dict(data.get("single_cell") or modeling.get("single_cell") or {}),
        multiomics=dict(data.get("multiomics") or {}),
    )
