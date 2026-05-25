from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mumoas.config import load_config
from mumoas.io import ensure_directory


PYTHON_STAGES = [
    ("make synthetic data", Path("scripts/00_make_synthetic_data.py")),
    ("extract imaging features", Path("scripts/01_extract_imaging_features.py")),
    ("train cVAE RCB model", Path("scripts/02_train_cvae_rcb.py")),
    ("radiogenomic association", Path("scripts/03_radiogenomic_association.py")),
    ("evaluate predictions", Path("scripts/06_evaluate_predictions.py")),
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run public MuMoAS smoke checks on synthetic Python inputs and R code parsing.")
    parser.add_argument("--config", default="configs/public_example.yaml", help="Path to the Python YAML configuration file.")
    parser.add_argument("--r-config", default="analysis_r/config/mumoas_sc_multiomics_template.yml", help="Path to the R YAML configuration template.")
    parser.add_argument("--skip-r", action="store_true", help="Skip the R parse/config smoke check.")
    return parser


def run_command(command: list[str], stage_name: str) -> None:
    result = subprocess.run(command, cwd=REPO_ROOT)
    if result.returncode != 0:
        print(f"Smoke test failed at stage: {stage_name}", file=sys.stderr)
        raise SystemExit(result.returncode)


def run_python_stages(config: Path) -> None:
    for stage_name, script_path in PYTHON_STAGES:
        run_command([sys.executable, str(REPO_ROOT / script_path), "--config", str(config)], stage_name)


def run_r_smoke(r_config: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        print("Rscript was not found; skipping R smoke parse check.", file=sys.stderr)
        return
    run_command([rscript, str(REPO_ROOT / "analysis_r/smoke_parse.R"), "--config", str(r_config)], "R analysis parse/config check")


def run(config_path: str | Path, r_config_path: str | Path, skip_r: bool = False) -> Path:
    config = Path(config_path)
    if not config.is_absolute():
        config = (REPO_ROOT / config).resolve()
    r_config = Path(r_config_path)
    if not r_config.is_absolute():
        r_config = (REPO_ROOT / r_config).resolve()
    cfg = load_config(config)
    run_python_stages(config)
    if not skip_r:
        run_r_smoke(r_config)
    completion_path = ensure_directory(cfg.project.output_dir) / "public_smoke_test_complete.txt"
    completion_path.write_text("MuMoAS public smoke test completed.\n", encoding="utf-8")
    print(f"Wrote completion marker: {completion_path}")
    return completion_path


def main() -> None:
    args = build_parser().parse_args()
    try:
        run(args.config, args.r_config, skip_r=args.skip_r)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(f"{exc}\nCheck configuration paths and rerun scripts/run_public_smoke_test.py.", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
