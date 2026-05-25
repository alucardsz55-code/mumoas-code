# MuMoAS R Analysis Workflow

Run R analyses from the repository root after copying the template configuration to a private path and replacing input paths with controlled data locations:

```bash
Rscript analysis_r/run_all.R --config analysis_r/config/mumoas_sc_multiomics_template.yml
```

Individual scripts are available under `analysis_r/scripts/`:

1. `01_dimension_clustering_lineages.R`: Seurat preprocessing, UMAP, major lineages, epithelial subclustering, cell proportions.
2. `02_marker_gene_dotplot.R`: canonical marker dot plots.
3. `03_diagonal_volcano.R`: epithelial response-group differential-expression scatter plot.
4. `04_infercnv_cnv_score.R`: inferCNV and CNV score plots.
5. `05_odds_ratio_heatmap.R`: cluster or cell-state odds-ratio heatmap across RCB groups.
6. `06_geneNMF_metaprograms.R`: geneNMF discovery of MP1-MP6.
7. `07_kegg_enrichment.R`: KEGG enrichment for MP gene sets.
8. `08_hallmark_gene_sets.R`: Hallmark overlap heatmap for MP gene sets.
9. `09_mp_rcb_signatures.R`: MP-RCB0 and MP-RCB3 scoring and UMAP feature maps.
10. `10_aucell_scores.R`: AUCell calculation for MP1-MP6 and MP-RCB0/MP-RCB3 response gene sets.
11. `11_rcb_group_scores.R`: sample-level MP scores by RCB group.
12. `12_cellchat_interactions.R`: CellChat interaction networks, source-state network plots, outgoing/incoming pathway dots, and ligand-receptor chord diagrams.
13. `13_gsea_visualization.R`: bulk and proteomic GSVA/ssGSEA, limma, fgsea, and GSEA panel plots.
14. `14_external_model_validation_stats.R`: external model AUC, confusion matrix, and calibration tables.
15. `15_validation_scrna_workflow.R`: validation scRNA-seq scoring and plots.
16. `16_mp_score_violin.R`: MP score violin plots.
17. `17_rcb_score_figures.R`: bulk RNA-seq and proteomic RCB score figures.
18. `18_enrichment_heatmap.R`: Hallmark enrichment heatmap.

`main_figure_panel_map()` in `analysis_r/R/mumoas_sc_multiomics.R` marks main-text Figure 5 and Figure 6 panels and the code step that generates each one.

The configuration includes required data-source manifest columns for site name and country of origin. Private manifests, matrices, Seurat objects, inferCNV objects, CellChat objects, and generated figures are excluded from the public repository.
