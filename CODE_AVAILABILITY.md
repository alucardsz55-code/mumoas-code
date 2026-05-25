# Code Availability

This repository provides public journal-review code for the MuMoAS study. It exposes the computational workflow, data interfaces, and analysis scripts without releasing controlled study data.

## Shared Code

- Python code under `mumoas/` and `scripts/` covers MRI feature validation, radiogenomic association, cVAE-based RCB prediction, and prediction evaluation.
- R code under `analysis_r/` covers the single-cell RNA-seq and multi-omics analyses used for the biological interpretation and validation figures.

## Data Access Boundaries

The public repository does not include raw MRI data, segmentation masks, image manifests, clinical records, patient-level metadata, bulk RNA-seq matrices, proteomic matrices, single-cell matrices, Seurat objects, inferCNV objects, CellChat objects, model checkpoints, trained weights, or generated manuscript outputs.

Reproducing manuscript results requires controlled access to the approved MRI, clinical, single-cell, transcriptomic, and proteomic datasets. Public synthetic Python inputs are provided only to demonstrate table contracts and code execution.

## Expected Data Interfaces

Python configuration uses `configs/public_example.yaml` as a public schema template for model-related tables:

- Clinical table: one row per patient with patient ID, cohort, RCB class, survival variables, categorical covariates, and continuous covariates.
- Imaging feature table: one row per patient with numeric MRI-derived features.
- Metaprogram score table: one row per patient with MP score columns for radiogenomic association.
- Prediction table: patient ID plus probability columns named `prob_<class>`.

R configuration uses `analysis_r/config/mumoas_sc_multiomics_template.yml` as a public schema template:

- scRNA-seq: controlled Seurat RDS files or sample manifests pointing to count matrices, with sample ID, site name, country, RCB group, predicted group, and cell-type metadata.
- inferCNV: gene-order file and non-malignant reference groups, or a controlled inferCNV object from an approved run.
- Bulk RNA-seq and proteomics: gene-by-sample or protein-by-sample matrices plus sample metadata containing sample ID, pCR group, RCB group, response group, batch/source, site name, and country.
- Model validation: prediction table with truth labels and class probability columns.
- Data-source manifests: one manifest per modality or cohort with `site_name` and `country` columns to document the source of each analysis step.

The R workflow explicitly maps main-text Figure 5 and Figure 6 panels to the generating analysis steps through `main_figure_panel_map()`. Supplementary biological analyses such as inferCNV score plots, diagonal volcano plots, and validation scRNA-seq feature maps are retained as numbered workflow outputs.

## Smoke Test

The public smoke test is a lightweight execution check. It runs the retained Python synthetic workflow and parses the R analysis code/configuration. It is not intended to reproduce manuscript results:

```bash
python scripts/run_public_smoke_test.py --config configs/public_example.yaml --r-config analysis_r/config/mumoas_sc_multiomics_template.yml
```

## Review Scope

Reviewers can inspect the complete code paths and configuration schemas in this repository. Running the manuscript-level analyses requires private configuration files pointing to controlled data locations. No private path or identifier should be committed to this repository.
