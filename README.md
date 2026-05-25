# MuMoAS Public Review Code

This repository contains public journal-review code for the MuMoAS study without controlled study data, private paths, trained weights, patient identifiers, or manuscript figure outputs.

## Repository Structure

- `mumoas/`: Python package for MRI feature validation, radiogenomic association, cVAE modeling, plotting, and model evaluation.
- `scripts/`: Python command-line stages for synthetic table generation, imaging feature checks, cVAE training, radiogenomic association, and prediction evaluation.
- `analysis_r/`: R workflow for single-cell RNA-seq, epithelial meta-programs, inferCNV, AUCell/UCell-style scoring, CellChat, bulk RNA-seq, proteomics, enrichment, RCB score figures, and external biological validation.
- `configs/`: Public Python configuration template.
- `analysis_r/config/`: Public R configuration template with required site-name and country-origin manifest fields.
- `tests/`: Python tests for the retained modeling, imaging, radiogenomic, configuration, and evaluation code.

## Python Model Code

Use Python 3.10 or newer in an isolated environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

On Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
pip install -e .
```

The public synthetic table generator can be used to exercise the Python interfaces only:

```bash
python scripts/00_make_synthetic_data.py --config configs/public_example.yaml
python scripts/01_extract_imaging_features.py --config configs/public_example.yaml
python scripts/02_train_cvae_rcb.py --config configs/public_example.yaml
python scripts/03_radiogenomic_association.py --config configs/public_example.yaml
python scripts/06_evaluate_predictions.py --config configs/public_example.yaml
```

The public smoke test runs those retained Python stages on deterministic synthetic inputs and parses the R workflow/configuration:

```bash
python scripts/run_public_smoke_test.py --config configs/public_example.yaml --r-config analysis_r/config/mumoas_sc_multiomics_template.yml
```

The synthetic outputs demonstrate table contracts and software execution. They do not reproduce manuscript results, effect sizes, biological conclusions, or model performance.

## R Single-Cell And Multi-Omics Code

Install the R packages listed in `analysis_r/R_packages.txt`, copy `analysis_r/config/mumoas_sc_multiomics_template.yml` to a private location, and replace placeholder paths with controlled data paths.

Run the full R workflow from the repository root:

```bash
Rscript analysis_r/run_all.R --config analysis_r/config/mumoas_sc_multiomics_template.yml
```

Run an individual step:

```bash
Rscript analysis_r/scripts/06_geneNMF_metaprograms.R --config analysis_r/config/mumoas_sc_multiomics_template.yml
```

The R workflow contains the 18 single-cell and multi-omics analysis steps used for the biological figures: Seurat clustering and proportions, marker plots, diagonal volcano plots, inferCNV, OR heatmaps, geneNMF MP discovery, KEGG, Hallmark enrichment, MP-RCB AUCell scoring and high/low grouping, CellChat source networks and ligand-receptor chord diagrams, GSVA/ssGSEA and fgsea statistics, external validation statistics, validation scRNA-seq scoring, MP violins, RCB score figures, and enrichment heatmaps.

Main-text biological panels are explicitly mapped by `main_figure_panel_map()` in `analysis_r/R/mumoas_sc_multiomics.R`; running `analysis_r/run_all.R` writes `outputs/r_analysis/main_figure_panel_code_map.tsv`. The main-text biological figures covered by the R workflow are Figure 5 and Figure 6.

## Data Boundaries

The repository does not include raw MRI data, segmentation masks, clinical tables, omics matrices, single-cell matrices, Seurat objects, inferCNV objects, CellChat objects, trained weights, checkpoints, private configurations, manuscript files, or generated figures.

Private runs should use local controlled configuration files outside version control. The R configuration template includes source-manifest fields for `site_name` and `country` so each analysis step can be linked to its data source.

## Tests

Run retained Python tests with:

```bash
pytest
```
