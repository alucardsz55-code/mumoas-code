from pathlib import Path


def test_public_smoke_runner_exists():
    assert Path("scripts/run_public_smoke_test.py").exists()


def test_public_smoke_runner_uses_retained_python_and_r_parse_stages():
    runner_text = Path("scripts/run_public_smoke_test.py").read_text(encoding="utf-8")
    expected = [
        "scripts/00_make_synthetic_data.py",
        "scripts/01_extract_imaging_features.py",
        "scripts/02_train_cvae_rcb.py",
        "scripts/03_radiogenomic_association.py",
        "scripts/06_evaluate_predictions.py",
        "analysis_r/smoke_parse.R",
    ]
    removed = [
        "scripts/" + "04_" + "single_cell_" + "metaprograms.py",
        "scripts/" + "05_" + "multiomics_" + "validation.py",
    ]
    for item in expected:
        assert item in runner_text
    for item in removed:
        assert item not in runner_text
