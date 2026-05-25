`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) y else x
}

mumoas_pkg <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) stop(paste("Missing R packages:", paste(missing, collapse = ", ")), call. = FALSE)
}

cfg_get <- function(x, path, default = NULL) {
  value <- x
  for (key in path) {
    if (is.null(value) || is.null(value[[key]])) return(default)
    value <- value[[key]]
  }
  value
}

mumoas_read_config <- function(path) {
  mumoas_pkg("yaml")
  cfg <- yaml::read_yaml(path)
  cfg$.config_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  cfg$.base_dir <- normalizePath(cfg_get(cfg, c("project", "base_dir"), getwd()), winslash = "/", mustWork = FALSE)
  cfg$.output_dir <- mumoas_path(cfg, cfg_get(cfg, c("project", "output_dir"), "outputs/r_analysis"))
  dir.create(cfg$.output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg
}

mumoas_args_config <- function(default = "analysis_r/config/mumoas_sc_multiomics_template.yml") {
  args <- commandArgs(trailingOnly = TRUE)
  index <- match("--config", args)
  if (!is.na(index) && length(args) >= index + 1) args[[index + 1]] else default
}

mumoas_path <- function(cfg, path) {
  if (is.null(path) || identical(path, "")) return(NULL)
  if (grepl("^[A-Za-z]:[/\\\\]|^/", path)) return(normalizePath(path, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(cfg$.base_dir, path), winslash = "/", mustWork = FALSE)
}

mumoas_output <- function(cfg, ...) {
  path <- file.path(cfg$.output_dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

read_table_auto <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!file.exists(path)) stop(paste("Input file not found:", path), call. = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) return(as.data.frame(data.table::fread(path, data.table = FALSE, check.names = FALSE)))
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  utils::read.table(path, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_table_auto <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  invisible(path)
}

save_plot_pdf <- function(plot, path, width = 8, height = 6) {
  mumoas_pkg("ggplot2")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, device = grDevices::cairo_pdf, limitsize = FALSE)
  invisible(path)
}

save_plot_png <- function(plot, path, width = 8, height = 6, dpi = 300) {
  mumoas_pkg("ggplot2")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  invisible(path)
}

get_assay_matrix <- function(object, assay = NULL, slot = "data") {
  mumoas_pkg("Seurat")
  assay <- assay %||% Seurat::DefaultAssay(object)
  out <- tryCatch(Seurat::GetAssayData(object, assay = assay, layer = slot), error = function(e) NULL)
  if (is.null(out)) out <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  out
}

config_gene_sets <- function(cfg) {
  sets <- cfg_get(cfg, c("signatures"), list())
  sets <- sets[!vapply(sets, is.null, logical(1))]
  lapply(sets, unique)
}

response_programs <- function(cfg) {
  cfg_get(cfg, c("response_programs"), list())
}

marker_sets <- function(cfg, dataset = "discovery") {
  cfg_get(cfg, c("single_cell", dataset, "marker_sets"), cfg_get(cfg, c("single_cell", "marker_sets"), list()))
}

load_seurat_from_manifest <- function(cfg, dataset = "discovery") {
  mumoas_pkg(c("Seurat", "Matrix"))
  sc <- cfg_get(cfg, c("single_cell", dataset), list())
  manifest_path <- mumoas_path(cfg, sc$sample_manifest)
  manifest <- read_table_auto(manifest_path)
  counts_column <- sc$counts_column %||% "counts_path"
  sample_column <- sc$sample_id_column %||% "sample_id"
  min_cells <- cfg_get(cfg, c("single_cell", "preprocess", "min_cells"), 3)
  min_features <- cfg_get(cfg, c("single_cell", "preprocess", "min_features"), 200)
  objects <- lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, , drop = FALSE]
    counts_path <- mumoas_path(cfg, row[[counts_column]])
    counts <- if (dir.exists(counts_path)) Seurat::Read10X(counts_path) else read_table_auto(counts_path)
    if (is.data.frame(counts)) {
      rownames(counts) <- counts[[1]]
      counts <- as.matrix(counts[, -1, drop = FALSE])
    }
    object <- Seurat::CreateSeuratObject(counts = counts, min.cells = min_cells, min.features = min_features, project = as.character(row[[sample_column]]))
    for (column in setdiff(colnames(row), counts_column)) object[[column]] <- row[[column]]
    object
  })
  Reduce(function(x, y) merge(x, y), objects)
}

load_seurat_dataset <- function(cfg, dataset = "discovery") {
  mumoas_pkg("Seurat")
  sc <- cfg_get(cfg, c("single_cell", dataset), list())
  rds <- mumoas_path(cfg, sc$seurat_rds)
  if (!is.null(rds) && file.exists(rds)) return(readRDS(rds))
  load_seurat_from_manifest(cfg, dataset)
}

prepare_seurat_object <- function(object, cfg, dataset = "discovery") {
  mumoas_pkg(c("Seurat", "ggplot2"))
  set.seed(cfg_get(cfg, c("project", "seed"), 20260525))
  pp <- cfg_get(cfg, c("single_cell", "preprocess"), list())
  if (!"percent.mt" %in% colnames(object@meta.data)) object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = pp$mitochondrial_pattern %||% "^MT-")
  if (!"percent.rb" %in% colnames(object@meta.data)) object[["percent.rb"]] <- Seurat::PercentageFeatureSet(object, pattern = pp$ribosomal_pattern %||% "^RP[SL]|^RPLP|^RPSA")
  if (!"percent.hb" %in% colnames(object@meta.data)) object[["percent.hb"]] <- Seurat::PercentageFeatureSet(object, pattern = pp$hemoglobin_pattern %||% "^HB[^(P)]")
  object <- subset(object, subset = nFeature_RNA >= (pp$nfeature_min %||% 200) & nFeature_RNA <= (pp$nfeature_max %||% Inf) & percent.mt <= (pp$percent_mt_max %||% 25))
  method <- pp$normalization %||% "LogNormalize"
  if (identical(method, "SCT")) {
    object <- Seurat::SCTransform(object, vars.to.regress = pp$vars_to_regress %||% c("percent.mt", "percent.rb"), verbose = FALSE)
  } else {
    object <- Seurat::NormalizeData(object, verbose = FALSE)
    object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = pp$nfeatures %||% 2000, verbose = FALSE)
    object <- Seurat::ScaleData(object, vars.to.regress = pp$vars_to_regress %||% NULL, verbose = FALSE)
  }
  npcs <- pp$npcs %||% 40
  object <- Seurat::RunPCA(object, npcs = npcs, verbose = FALSE)
  reduction <- "pca"
  if (identical(pp$integration, "harmony")) {
    mumoas_pkg("harmony")
    batch_column <- pp$batch_column %||% "orig.ident"
    object <- harmony::RunHarmony(object, group.by.vars = batch_column, reduction.use = "pca", dims.use = seq_len(npcs), verbose = FALSE)
    reduction <- "harmony"
  }
  dims <- seq_len(pp$dims %||% min(30, npcs))
  object <- Seurat::RunUMAP(object, reduction = reduction, dims = dims, verbose = FALSE)
  object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims, verbose = FALSE)
  object <- Seurat::FindClusters(object, resolution = pp$resolution %||% 0.8, verbose = FALSE)
  object
}

load_or_prepare_seurat <- function(cfg, dataset = "discovery") {
  object <- load_seurat_dataset(cfg, dataset)
  run_preprocess <- cfg_get(cfg, c("single_cell", dataset, "run_preprocess"), FALSE)
  has_umap <- "umap" %in% names(object@reductions)
  if (isTRUE(run_preprocess) || !has_umap) object <- prepare_seurat_object(object, cfg, dataset)
  out <- mumoas_output(cfg, "single_cell", dataset, paste0(dataset, "_seurat_processed.rds"))
  saveRDS(object, out)
  object
}

metadata_umap <- function(object) {
  mumoas_pkg("Seurat")
  emb <- as.data.frame(Seurat::Embeddings(object, "umap"))
  colnames(emb)[seq_len(2)] <- c("UMAP_1", "UMAP_2")
  cbind(object@meta.data, emb[colnames(object), , drop = FALSE])
}

plot_umap <- function(object, group_col, label_col = group_col, point_size = 0.08, label = TRUE, title = NULL) {
  mumoas_pkg(c("ggplot2", "ggrepel"))
  df <- metadata_umap(object)
  df[[group_col]] <- factor(df[[group_col]])
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$UMAP_1, .data$UMAP_2, color = .data[[group_col]])) +
    ggplot2::geom_point(size = point_size, alpha = 0.9) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2", color = NULL, title = title) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.key.height = grid::unit(0.35, "cm"), legend.key.width = grid::unit(0.35, "cm"))
  if (isTRUE(label)) {
    centers <- stats::aggregate(df[, c("UMAP_1", "UMAP_2")], by = list(label = df[[label_col]]), FUN = median)
    p <- p + ggrepel::geom_text_repel(data = centers, ggplot2::aes(UMAP_1, UMAP_2, label = .data$label), inherit.aes = FALSE, size = 3.5, fontface = "bold", color = "black")
  }
  p
}

plot_celltype_proportion <- function(object, sample_col, group_col, celltype_col) {
  mumoas_pkg(c("ggplot2", "dplyr"))
  md <- object@meta.data
  md[[sample_col]] <- factor(md[[sample_col]], levels = unique(md[[sample_col]][order(md[[group_col]])]))
  df <- md |>
    dplyr::count(.data[[sample_col]], .data[[group_col]], .data[[celltype_col]], name = "n") |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::mutate(prop = .data$n / sum(.data$n)) |>
    dplyr::ungroup()
  ggplot2::ggplot(df, ggplot2::aes(.data[[sample_col]], .data$prop, fill = .data[[celltype_col]])) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::facet_grid(. ~ .data[[group_col]], scales = "free_x", space = "free_x") +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(x = NULL, y = "Proportion", fill = "celltype") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank())
}

plot_marker_dot <- function(object, features, group_col) {
  mumoas_pkg(c("Seurat", "ggplot2"))
  Seurat::Idents(object) <- object@meta.data[[group_col]]
  Seurat::DotPlot(object, features = features) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank()) +
    ggplot2::scale_color_gradientn(colours = grDevices::colorRampPalette(c("#2c7fb8", "white", "#d81b60"))(100))
}

subcluster_cells <- function(object, cfg, labels, output_name) {
  mumoas_pkg("Seurat")
  celltype_col <- cfg_get(cfg, c("single_cell", "celltype_column"), "celltype")
  cells <- rownames(object@meta.data)[object@meta.data[[celltype_col]] %in% labels]
  sub <- subset(object, cells = cells)
  sub <- Seurat::NormalizeData(sub, verbose = FALSE)
  sub <- Seurat::FindVariableFeatures(sub, verbose = FALSE)
  sub <- Seurat::ScaleData(sub, verbose = FALSE)
  sub <- Seurat::RunPCA(sub, npcs = 40, verbose = FALSE)
  sub <- Seurat::RunUMAP(sub, dims = 1:30, verbose = FALSE)
  sub <- Seurat::FindNeighbors(sub, dims = 1:30, verbose = FALSE)
  sub <- Seurat::FindClusters(sub, resolution = cfg_get(cfg, c("single_cell", "epithelial_resolution"), 1.0), verbose = FALSE)
  saveRDS(sub, mumoas_output(cfg, "single_cell", output_name, paste0(output_name, "_subclustered.rds")))
  sub
}

score_seurat_gene_sets <- function(object, gene_sets, assay = NULL, slot = "data") {
  mumoas_pkg("Matrix")
  mat <- get_assay_matrix(object, assay, slot)
  for (name in names(gene_sets)) {
    genes <- intersect(gene_sets[[name]], rownames(mat))
    object[[name]] <- if (length(genes) == 0) NA_real_ else Matrix::colMeans(mat[genes, , drop = FALSE])
  }
  object
}

composite_scores <- function(meta, cfg) {
  programs <- response_programs(cfg)
  for (name in names(programs)) {
    pos <- unlist(programs[[name]]$positive %||% character())
    neg <- unlist(programs[[name]]$negative %||% character())
    score <- rep(0, nrow(meta))
    if (length(pos) > 0) score <- score + rowMeans(meta[, intersect(pos, colnames(meta)), drop = FALSE], na.rm = TRUE)
    if (length(neg) > 0) score <- score - rowMeans(meta[, intersect(neg, colnames(meta)), drop = FALSE], na.rm = TRUE)
    meta[[name]] <- score
  }
  meta
}

plot_score_violin <- function(meta, scores, group_col, palette = NULL) {
  mumoas_pkg(c("ggplot2", "tidyr", "dplyr", "ggpubr"))
  df <- meta |>
    dplyr::select(dplyr::all_of(c(group_col, scores))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(scores), names_to = "program", values_to = "score")
  ggplot2::ggplot(df, ggplot2::aes(.data[[group_col]], .data$score, fill = .data[[group_col]])) +
    ggplot2::geom_violin(trim = FALSE, linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.2) +
    ggplot2::facet_wrap(~program, scales = "free_y", ncol = min(4, length(scores))) +
    ggpubr::stat_compare_means(size = 2.5) +
    ggplot2::labs(x = NULL, y = "Score", fill = NULL) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
}

plot_score_density_umap <- function(object, scores) {
  mumoas_pkg(c("ggplot2", "patchwork"))
  df <- metadata_umap(object)
  plots <- lapply(scores, function(score) {
    ggplot2::ggplot(df, ggplot2::aes(.data$UMAP_1, .data$UMAP_2, color = .data[[score]])) +
      ggplot2::geom_point(size = 0.08) +
      ggplot2::scale_color_viridis_c(option = "magma", na.value = "grey90") +
      ggplot2::labs(title = score, x = "UMAP1", y = "UMAP2", color = NULL) +
      ggplot2::theme_void(base_size = 9) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black", color = NA), plot.title = ggplot2::element_text(color = "black", hjust = 0.5))
  })
  patchwork::wrap_plots(plots, ncol = min(4, length(plots)))
}

find_threshold <- function(score, label) {
  if (length(unique(stats::na.omit(label))) < 2 || !requireNamespace("pROC", quietly = TRUE)) return(stats::median(score, na.rm = TRUE))
  roc <- pROC::roc(label, score, quiet = TRUE, direction = "auto")
  as.numeric(pROC::coords(roc, "best", ret = "threshold", transpose = FALSE)[1, 1])
}

main_figure_panel_map <- function() {
  data.frame(
    figure = c(rep("Figure 5", 9), rep("Figure 6", 19)),
    panel = c("5a", "5b", "5c", "5d", "5e", "5f", "5g", "5h", "5i", "6a", "6b", "6c", "6d", "6e", "6f", "6g", "6h", "6i", "6j", "6k", "6l", "6m", "6n", "6o", "6p", "6q", "6r", "6s"),
    analysis_step = c(
      "01 lineage UMAP", "01 cell-type proportions", "02 marker dot plot", "01 epithelial subclustering", "06 geneNMF Jaccard heatmap", "16 MP1-MP6 score violin", "09 MP feature maps", "06 MP annotation table", "18 Hallmark overlap heatmap",
      "09 AUCell MP-RCB scoring and 16 violin plot", "05 OR heatmap for MP-RCB0", "05 OR heatmap for MP-RCB3", "12 CellChat MP-RCB3 source network", "12 CellChat SPP1-CD44 chord", "12 CellChat VEGFA-VEGFR2 chord", "12 CellChat TIGIT-NECTIN2 chord", "12 CellChat MP-RCB0 source network", "12 CellChat TNC-SDC1 chord", "12 CellChat THBS2-SDC1 chord", "12 CellChat GRN-SORT1 chord", "13 RNA-seq GSEA MP-RCB0", "13 RNA-seq GSEA MP-RCB3", "13 proteomic GSEA MP-RCB0", "13 proteomic GSEA MP-RCB3", "17 in-house RNA pCR score validation", "17 in-house RNA RCB score validation", "17 external RNA pCR score validation", "17 external RNA RCB score validation"
    ),
    script = c(
      "analysis_r/scripts/01_dimension_clustering_lineages.R", "analysis_r/scripts/01_dimension_clustering_lineages.R", "analysis_r/scripts/02_marker_gene_dotplot.R", "analysis_r/scripts/01_dimension_clustering_lineages.R", "analysis_r/scripts/06_geneNMF_metaprograms.R", "analysis_r/scripts/16_mp_score_violin.R", "analysis_r/scripts/09_mp_rcb_signatures.R", "analysis_r/scripts/06_geneNMF_metaprograms.R", "analysis_r/scripts/18_enrichment_heatmap.R",
      "analysis_r/scripts/09_mp_rcb_signatures.R", "analysis_r/scripts/05_odds_ratio_heatmap.R", "analysis_r/scripts/05_odds_ratio_heatmap.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/12_cellchat_interactions.R", "analysis_r/scripts/13_gsea_visualization.R", "analysis_r/scripts/13_gsea_visualization.R", "analysis_r/scripts/13_gsea_visualization.R", "analysis_r/scripts/13_gsea_visualization.R", "analysis_r/scripts/17_rcb_score_figures.R", "analysis_r/scripts/17_rcb_score_figures.R", "analysis_r/scripts/17_rcb_score_figures.R", "analysis_r/scripts/17_rcb_score_figures.R"
    ),
    stringsAsFactors = FALSE
  )
}

write_main_figure_panel_map <- function(cfg) {
  write_table_auto(main_figure_panel_map(), mumoas_output(cfg, "main_figure_panel_code_map.tsv"))
}

load_rdata_object <- function(path, preferred = NULL) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  names <- ls(env)
  if (!is.null(preferred) && preferred %in% names) return(env[[preferred]])
  env[[names[[1]]]]
}

format_hallmark_name <- function(pathway) {
  out <- gsub("^HALLMARK_", "", pathway)
  out <- gsub("_", " ", out)
  tools::toTitleCase(tolower(out))
}

top_p_to_overlap_table <- function(top_p) {
  records <- list()
  k <- 1
  for (program in names(top_p)) {
    df <- as.data.frame(top_p[[program]])
    if (!all(c("pathway", "overlap", "size") %in% colnames(df))) next
    df$program <- program
    df$pathway_label <- vapply(df$pathway, format_hallmark_name, character(1))
    df$overlap <- as.numeric(df$overlap)
    df$size <- as.numeric(df$size)
    df$gene_ratio <- ifelse(df$size == 0, 0, df$overlap / df$size)
    records[[k]] <- df
    k <- k + 1
  }
  if (length(records) == 0) data.frame() else do.call(rbind, records)
}

extract_overlap_genes <- function(top_p, programs) {
  genes <- character()
  for (program in programs) {
    item <- top_p[[program]]
    if (is.null(item) && grepl("^[0-9]+$", as.character(program))) item <- top_p[[as.integer(program)]]
    if (is.null(item) || is.null(item$overlapGenes)) next
    genes <- c(genes, unique(unlist(item$overlapGenes, use.names = FALSE)))
  }
  unique(stats::na.omit(genes))
}

derived_response_gene_sets <- function(cfg, fallback_sets = config_gene_sets(cfg)) {
  programs <- response_programs(cfg)
  top_p_path <- mumoas_path(cfg, cfg_get(cfg, c("single_cell", "geneNMF", "hallmark_gsea_rdata"), NULL))
  top_p <- NULL
  if (!is.null(top_p_path) && file.exists(top_p_path)) top_p <- load_rdata_object(top_p_path, "top_p")
  out <- list()
  for (name in names(programs)) {
    overlap_programs <- programs[[name]]$derive_from_hallmark_overlap
    if (!is.null(top_p) && length(overlap_programs) > 0) out[[name]] <- extract_overlap_genes(top_p, unlist(overlap_programs))
    if (is.null(out[[name]]) || length(out[[name]]) == 0) {
      pos <- unlist(programs[[name]]$positive %||% character())
      out[[name]] <- unique(unlist(fallback_sets[intersect(pos, names(fallback_sets))], use.names = FALSE))
    }
  }
  out
}

compute_aucell_scores <- function(object, gene_sets) {
  mumoas_pkg("AUCell")
  rankings <- AUCell::AUCell_buildRankings(get_assay_matrix(object, slot = "counts"), plotStats = FALSE, verbose = FALSE)
  auc <- AUCell::AUCell_calcAUC(gene_sets, rankings, verbose = FALSE)
  scores <- as.data.frame(t(as.matrix(AUCell::getAUC(auc))), check.names = FALSE)
  scores$cell <- rownames(scores)
  scores
}

add_aucell_scores_to_seurat <- function(object, gene_sets) {
  scores <- compute_aucell_scores(object, gene_sets)
  rownames(scores) <- scores$cell
  score_cols <- setdiff(colnames(scores), "cell")
  object@meta.data[rownames(scores), score_cols] <- scores[rownames(object@meta.data), score_cols, drop = FALSE]
  object
}

plot_auc_histogram <- function(meta, score_col, group_col, positive_group, threshold) {
  mumoas_pkg("ggplot2")
  ggplot2::ggplot(meta, ggplot2::aes(.data[[score_col]])) +
    ggplot2::geom_histogram(bins = 60, fill = "#2b8cbe", color = "white", alpha = 0.85) +
    ggplot2::geom_density(color = "#08589e", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = threshold, color = "#d7191c", linewidth = 0.9) +
    ggplot2::labs(x = paste0(score_col, " AUCell score"), y = "Frequency", title = paste0(score_col, " high/low threshold")) +
    ggplot2::theme_classic(base_size = 9)
}

apply_response_aucell <- function(object, cfg) {
  base_sets <- config_gene_sets(cfg)
  response_sets <- derived_response_gene_sets(cfg, base_sets)
  if (length(response_sets) > 0) object <- add_aucell_scores_to_seurat(object, response_sets)
  group_col <- cfg_get(cfg, c("single_cell", "group_column"), "rcb_group")
  thresholds <- data.frame(program = character(), threshold = numeric(), positive_group = character())
  for (program in names(response_sets)) {
    positive_group <- cfg_get(cfg, c("response_programs", program, "positive_group"), NULL)
    label <- if (!is.null(positive_group) && group_col %in% colnames(object@meta.data)) object@meta.data[[group_col]] == positive_group else NULL
    threshold <- if (is.null(label)) stats::median(object@meta.data[[program]], na.rm = TRUE) else find_threshold(object@meta.data[[program]], label)
    object@meta.data[[paste0(program, "_group")]] <- ifelse(object@meta.data[[program]] >= threshold, paste0(program, "_Hi"), paste0(program, "_Lo"))
    thresholds <- rbind(thresholds, data.frame(program = program, threshold = threshold, positive_group = positive_group %||% "median"))
  }
  attr(object, "mumoas_response_thresholds") <- thresholds
  object
}

plot_hallmark_overlap_heatmap <- function(table) {
  mumoas_pkg(c("ggplot2", "dplyr"))
  if ("pathway_label" %in% colnames(table)) table$pathway <- table$pathway_label else table$pathway <- gsub("^HALLMARK_", "", table$pathway)
  ggplot2::ggplot(table, ggplot2::aes(.data$pathway, .data$program, fill = .data$gene_ratio)) +
    ggplot2::geom_tile(color = "grey85", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#325EA1", mid = "#FBF7EB", high = "#E7648E", midpoint = stats::median(table$gene_ratio, na.rm = TRUE)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Overlap ratio", title = "GSEA overlapped gene ratio of HALLMARK gene set") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
}

read_expression_matrix <- function(path) {
  table <- read_table_auto(path)
  rownames(table) <- make.unique(as.character(table[[1]]))
  mat <- as.matrix(table[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  mat
}

log2_if_needed <- function(mat) {
  qx <- stats::quantile(mat, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  if ((qx[[5]] > 100) || ((qx[[6]] - qx[[1]]) > 50 && qx[[2]] > 0)) mat <- log2(pmax(mat, 0) + 1)
  mat
}

combat_matrix <- function(mat, sample_meta, sample_col, batch_col) {
  if (is.null(batch_col) || !batch_col %in% colnames(sample_meta)) return(mat)
  mumoas_pkg("sva")
  sample_meta <- sample_meta[match(colnames(mat), sample_meta[[sample_col]]), , drop = FALSE]
  sva::ComBat(mat, batch = sample_meta[[batch_col]], mod = NULL, par.prior = TRUE)
}

matrix_gene_set_scores <- function(mat, gene_sets, scale_rows = TRUE) {
  if (isTRUE(scale_rows)) mat <- t(scale(t(mat)))
  scores <- sapply(gene_sets, function(genes) {
    genes <- intersect(genes, rownames(mat))
    if (length(genes) == 0) rep(NA_real_, ncol(mat)) else colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
  })
  scores <- as.data.frame(scores, check.names = FALSE)
  scores$sample_id <- rownames(scores)
  scores
}

rank_features_limma <- function(mat, sample_meta, sample_col, group_col, contrast) {
  mumoas_pkg("limma")
  sample_meta <- sample_meta[match(colnames(mat), sample_meta[[sample_col]]), , drop = FALSE]
  group <- factor(sample_meta[[group_col]])
  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(mat, design)
  cm <- limma::makeContrasts(contrasts = contrast, levels = design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, cm))
  table <- limma::topTable(fit, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  table$feature <- rownames(table)
  table
}

hallmark_overlap <- function(gene_sets, universe) {
  mumoas_pkg(c("msigdbr", "dplyr"))
  msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") |>
    dplyr::select(gs_name, gene_symbol) |>
    dplyr::distinct()
  hallmarks <- split(msig$gene_symbol, msig$gs_name)
  records <- list()
  k <- 1
  for (program in names(gene_sets)) {
    genes <- unique(intersect(gene_sets[[program]], universe))
    for (pathway in names(hallmarks)) {
      pathway_genes <- unique(intersect(hallmarks[[pathway]], universe))
      overlap <- length(intersect(genes, pathway_genes))
      p <- stats::fisher.test(matrix(c(overlap, length(genes) - overlap, length(pathway_genes) - overlap, length(universe) - length(genes) - length(pathway_genes) + overlap), nrow = 2))$p.value
      records[[k]] <- data.frame(program = program, pathway = pathway, overlap = overlap, gene_ratio = overlap / max(length(genes), 1), p_value = p)
      k <- k + 1
    }
  }
  out <- do.call(rbind, records)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

kegg_enrichment <- function(gene_sets) {
  mumoas_pkg(c("clusterProfiler", "org.Hs.eg.db"))
  records <- list()
  for (program in names(gene_sets)) {
    ids <- suppressMessages(clusterProfiler::bitr(gene_sets[[program]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
    if (nrow(ids) == 0) next
    kk <- suppressMessages(clusterProfiler::enrichKEGG(gene = ids$ENTREZID, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1))
    tab <- as.data.frame(kk)
    if (nrow(tab) == 0) next
    tab$program <- program
    records[[program]] <- tab
  }
  if (length(records) == 0) data.frame() else do.call(rbind, records)
}

msigdb_hallmark_sets <- function() {
  mumoas_pkg(c("msigdbr", "dplyr"))
  msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") |>
    dplyr::select(gs_name, gene_symbol) |>
    dplyr::distinct()
  split(msig$gene_symbol, msig$gs_name)
}

plot_enrichment_bar <- function(table, program, top_n = 10) {
  mumoas_pkg(c("ggplot2", "dplyr"))
  df <- table |>
    dplyr::filter(.data$program == !!program) |>
    dplyr::arrange(.data$pvalue) |>
    dplyr::slice_head(n = top_n)
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())
  df$fold_enrichment <- vapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)) / vapply(strsplit(df$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  df$Description <- factor(df$Description, levels = rev(df$Description))
  ggplot2::ggplot(df, ggplot2::aes(.data$fold_enrichment, .data$Description, fill = -log10(.data$pvalue))) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "#c51b7d") +
    ggplot2::labs(x = "Fold Enrichment", y = NULL, fill = "-log10(pvalue)", title = paste("KEGG pathways", program)) +
    ggplot2::theme_classic(base_size = 9)
}

plot_enrichment_heatmap <- function(table, value_col = "gene_ratio") {
  mumoas_pkg(c("ggplot2", "dplyr", "tidyr"))
  df <- table |>
    dplyr::group_by(.data$program) |>
    dplyr::slice_min(.data$fdr, n = 10, with_ties = FALSE) |>
    dplyr::ungroup()
  df$pathway <- gsub("^HALLMARK_", "", df$pathway)
  ggplot2::ggplot(df, ggplot2::aes(.data$pathway, .data$program, fill = .data[[value_col]])) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d81b60", midpoint = stats::median(df[[value_col]], na.rm = TRUE)) +
    ggplot2::labs(x = NULL, y = NULL, fill = value_col) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

gsva_score_matrix <- function(mat, gene_sets, method = "ssgsea", kcdf = "Gaussian") {
  mumoas_pkg("GSVA")
  if (utils::packageVersion("GSVA") >= "1.50.0" && exists("ssgseaParam", asNamespace("GSVA"))) {
    param <- if (identical(method, "ssgsea")) GSVA::ssgseaParam(mat, gene_sets) else GSVA::gsvaParam(mat, gene_sets, kcdf = kcdf)
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(mat, gene_sets, method = method, kcdf = kcdf, abs.ranking = identical(method, "ssgsea"), verbose = FALSE)
}

limma_table_for_scores <- function(scores, sample_meta, sample_col, group_col, contrast) {
  mumoas_pkg("limma")
  sample_meta <- sample_meta[match(colnames(scores), sample_meta[[sample_col]]), , drop = FALSE]
  group <- factor(sample_meta[[group_col]])
  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(scores, design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, limma::makeContrasts(contrasts = contrast, levels = design)))
  table <- limma::topTable(fit, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  table$pathway <- rownames(table)
  table
}

plot_gsva_bar <- function(table) {
  mumoas_pkg(c("ggplot2", "dplyr"))
  df <- table
  df$pathway <- gsub("^HALLMARK_", "", df$pathway)
  df$threshold <- factor(ifelse(df$t >= 2, "Up", ifelse(df$t <= -2, "Down", "Not significant")), levels = c("Up", "Down", "Not significant"))
  df <- df[order(df$t), , drop = FALSE]
  df$pathway <- factor(df$pathway, levels = df$pathway)
  ggplot2::ggplot(df, ggplot2::aes(.data$pathway, .data$t, fill = .data$threshold)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Up" = "#2166AC", "Down" = "#67A65B", "Not significant" = "#cccccc")) +
    ggplot2::geom_hline(yintercept = c(-2, 2), color = "white", linetype = 2, linewidth = 0.3) +
    ggplot2::labs(x = NULL, y = "t value of GSVA score", fill = NULL) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
}

fgsea_results <- function(rank_table, gene_sets, score_col = "t") {
  mumoas_pkg("fgsea")
  ranks <- rank_table[[score_col]]
  names(ranks) <- rank_table$feature
  ranks <- sort(stats::na.omit(ranks), decreasing = TRUE)
  out <- fgsea::fgsea(pathways = gene_sets, stats = ranks, minSize = 5, maxSize = 500, eps = 0)
  as.data.frame(out)
}

plot_gsea_panel <- function(rank_table, gene_sets, pathways, title = NULL, score_col = "t") {
  mumoas_pkg(c("fgsea", "ggplot2", "patchwork"))
  ranks <- rank_table[[score_col]]
  names(ranks) <- rank_table$feature
  ranks <- sort(stats::na.omit(ranks), decreasing = TRUE)
  plots <- lapply(pathways, function(pathway) {
    genes <- gene_sets[[pathway]]
    if (is.null(genes)) return(ggplot2::ggplot() + ggplot2::theme_void())
    fgsea::plotEnrichment(genes, ranks) +
      ggplot2::labs(title = pathway, x = NULL, y = "Running enrichment score") +
      ggplot2::theme_classic(base_size = 8)
  })
  patchwork::wrap_plots(plots, ncol = 1) + patchwork::plot_annotation(title = title)
}

cellchat_source_network_plot <- function(cellchat, source_celltype, count = TRUE, show_interaction = TRUE, color_set = NULL, celltype_order = NULL, celltype_size = FALSE, flipped = TRUE) {
  mumoas_pkg(c("reshape2", "ggplot2", "ggraph", "tidygraph", "dplyr", "igraph", "CellChat"))
  mat <- as.data.frame(if (isTRUE(count)) cellchat@net$count else cellchat@net$weight)
  mat <- mat[order(mat[, source_celltype], decreasing = TRUE), , drop = FALSE]
  group_size <- as.data.frame(table(cellchat@idents))
  df <- data.frame(from = rep(source_celltype, nrow(mat)), to = rownames(mat), inter = mat[, source_celltype])
  nodes <- data.frame(name = df$to, inter = df$inter)
  if (isTRUE(celltype_size)) {
    colnames(group_size) <- c("name", "size")
    nodes <- merge(nodes, group_size, by = "name", all = FALSE)
  } else {
    nodes$size <- 5
  }
  nodes <- nodes[order(nodes$inter, decreasing = TRUE), , drop = FALSE]
  color_use <- color_set %||% c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767", "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416")[seq_len(nrow(mat))]
  names(color_use) <- celltype_order %||% df$to
  net <- tidygraph::tbl_graph(nodes = nodes, edges = df[c("from", "to", "inter")])
  p <- ggraph::ggraph(net, layout = "igraph", algorithm = "circle") +
    ggraph::geom_edge_bend(ggplot2::aes(edge_width = inter), strength = 0.2, alpha = 0.8, flipped = flipped, edge_color = "#A9AAAA", n = 50, show.legend = FALSE) +
    ggraph::geom_edge_loop(ggplot2::aes(edge_width = inter), colour = "#A9AAAA", alpha = 0.5, show.legend = FALSE) +
    ggraph::scale_edge_width_continuous(range = c(0, 5))
  if (isTRUE(show_interaction)) {
    p + ggraph::geom_node_point(ggplot2::aes(size = inter, colour = inter)) +
      ggraph::geom_node_point(ggplot2::aes(size = inter), show.legend = FALSE, shape = 21, colour = "black", stroke = 1.2) +
      ggraph::geom_node_text(ggplot2::aes(x = x * 1.06, y = y * 1.06, label = name), hjust = 0, fontface = "bold", size = 3) +
      ggplot2::scale_color_gradientn(colors = grDevices::colorRampPalette(c("#2166AC", "#90C0DC", "white", "#EF8C65", "#B2182B"))(100)) +
      ggplot2::scale_size_continuous(range = c(1, 10)) +
      ggraph::theme_graph()
  } else {
    p + ggraph::geom_node_point(ggplot2::aes(size = size, colour = name), show.legend = FALSE) +
      ggraph::geom_node_point(ggplot2::aes(size = size), show.legend = FALSE, shape = 21, colour = "black", stroke = 1.2) +
      ggraph::geom_node_text(ggplot2::aes(x = x * 1.06, y = y * 1.06, label = name), hjust = 0, size = 3) +
      ggplot2::scale_color_manual(values = color_use) +
      ggplot2::scale_size_continuous(range = c(1, 15)) +
      ggraph::theme_graph()
  }
}

cellchat_communication_table <- function(cellchat) {
  mumoas_pkg("CellChat")
  as.data.frame(CellChat::subsetCommunication(cellchat))
}

plot_cellchat_lr_chord <- function(cellchat, pairs, sources = NULL, targets = NULL) {
  mumoas_pkg(c("circlize", "dplyr"))
  tab <- cellchat_communication_table(cellchat)
  tab$pair <- paste(tab$ligand, tab$receptor, sep = "-")
  tab <- tab[tab$pair %in% pairs | tab$interaction_name %in% pairs, , drop = FALSE]
  if (!is.null(sources)) tab <- tab[tab$source %in% sources, , drop = FALSE]
  if (!is.null(targets)) tab <- tab[tab$target %in% targets, , drop = FALSE]
  if (nrow(tab) == 0) return(invisible(NULL))
  circ <- tab[, c("source", "target", "prob")]
  colnames(circ) <- c("from", "to", "value")
  circlize::circos.clear()
  circlize::chordDiagram(circ, transparency = 0.25, directional = 1, direction.type = c("arrows", "diffHeight"), annotationTrack = "grid")
  invisible(tab)
}

plot_cellchat_role_dot <- function(cellchat, mode = c("outgoing", "incoming")) {
  mumoas_pkg(c("ggplot2", "reshape2"))
  mode <- match.arg(mode)
  prob <- cellchat@netP$prob
  pathways <- dimnames(prob)[[3]]
  cells <- dimnames(prob)[[1]]
  mat <- sapply(pathways, function(pathway) {
    x <- prob[, , pathway, drop = FALSE]
    if (identical(mode, "outgoing")) rowSums(x[, , 1, drop = FALSE]) else colSums(x[, , 1, drop = FALSE])
  })
  df <- reshape2::melt(mat)
  colnames(df) <- c("celltype", "pathway", "contribution")
  ggplot2::ggplot(df, ggplot2::aes(.data$pathway, .data$celltype, size = .data$contribution, color = .data$celltype)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::labs(x = NULL, y = NULL, size = "Contribution", title = paste(tools::toTitleCase(mode), "communication patterns")) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "right")
}

odds_ratio_table <- function(meta, feature_col, group_col) {
  groups <- unique(meta[[group_col]])
  features <- unique(meta[[feature_col]])
  records <- list()
  k <- 1
  for (feature in features) {
    for (group in groups) {
      a <- sum(meta[[feature_col]] == feature & meta[[group_col]] == group)
      b <- sum(meta[[feature_col]] == feature & meta[[group_col]] != group)
      c <- sum(meta[[feature_col]] != feature & meta[[group_col]] == group)
      d <- sum(meta[[feature_col]] != feature & meta[[group_col]] != group)
      ft <- stats::fisher.test(matrix(c(a, b, c, d), 2, 2))
      records[[k]] <- data.frame(feature = as.character(feature), group = as.character(group), or = unname(ft$estimate), p_value = ft$p.value)
      k <- k + 1
    }
  }
  out <- do.call(rbind, records)
  out$fdr <- stats::p.adjust(out$p_value, "BH")
  out
}

plot_or_heatmap <- function(table) {
  mumoas_pkg("ggplot2")
  table$stars <- cut(table$fdr, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))
  ggplot2::ggplot(table, ggplot2::aes(.data$group, .data$feature, fill = pmin(.data$or, 3))) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$stars), color = "#d7191c", size = 2.5) +
    ggplot2::scale_fill_gradient(low = "#3b4cc0", high = "#fdae61") +
    ggplot2::labs(x = NULL, y = NULL, fill = "OR") +
    ggplot2::theme_classic(base_size = 9)
}

diagonal_volcano <- function(object, subset_col, subset_value, group_col, group_low, group_high, p_cut = 0.05, logfc_cut = 0.5) {
  mumoas_pkg(c("Seurat", "ggplot2", "ggrepel"))
  Seurat::Idents(object) <- object@meta.data[[subset_col]]
  sub <- subset(object, idents = subset_value)
  Seurat::Idents(sub) <- sub@meta.data[[group_col]]
  de <- Seurat::FindMarkers(sub, ident.1 = group_high, ident.2 = group_low, logfc.threshold = 0, min.pct = 0.1, verbose = FALSE)
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(de))[1]
  avg <- as.data.frame(log1p(Seurat::AverageExpression(sub, group.by = group_col, assays = Seurat::DefaultAssay(sub), verbose = FALSE)[[Seurat::DefaultAssay(sub)]]))
  avg$gene <- rownames(avg)
  avg <- avg[rownames(de), , drop = FALSE]
  tab <- cbind(gene = rownames(de), de, avg[, c(group_low, group_high), drop = FALSE])
  tab$direction <- ifelse(tab$p_val_adj < p_cut & tab[[fc_col]] >= logfc_cut, group_high, ifelse(tab$p_val_adj < p_cut & tab[[fc_col]] <= -logfc_cut, group_low, "ns"))
  labels <- head(tab$gene[order(tab$p_val_adj)], 12)
  p <- ggplot2::ggplot(tab, ggplot2::aes(.data[[group_low]], .data[[group_high]], color = .data$direction)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.8) +
    ggrepel::geom_text_repel(data = tab[tab$gene %in% labels, ], ggplot2::aes(label = .data$gene), size = 2.5, color = "black") +
    ggplot2::scale_color_manual(values = c(stats::setNames(c("#66bfbf", "#c06c84", "grey80"), c(group_low, group_high, "ns")))) +
    ggplot2::labs(x = group_low, y = group_high, color = NULL, title = subset_value) +
    ggplot2::theme_classic(base_size = 9)
  list(table = tab, plot = p)
}

cnv_score_table <- function(infercnv_object, seurat_object) {
  expr <- infercnv_object@expr.data
  score <- data.frame(cell = colnames(expr), cnv_score = Matrix::colMeans((expr - 1) ^ 2))
  md <- seurat_object@meta.data
  md$cell <- rownames(md)
  merge(score, md, by = "cell")
}

plot_cnv_scores <- function(table, x_col, fill_col = x_col) {
  mumoas_pkg("ggplot2")
  ggplot2::ggplot(table, ggplot2::aes(.data[[x_col]], .data$cnv_score, fill = .data[[fill_col]])) +
    ggplot2::geom_violin(scale = "width", linewidth = 0.1, color = NA) +
    ggplot2::labs(x = NULL, y = "CNV score", fill = NULL) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
}

run_01_dimension_clustering_lineages <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "discovery")
  celltype_col <- cfg_get(cfg, c("single_cell", "celltype_column"), "celltype")
  sample_col <- cfg_get(cfg, c("single_cell", "sample_id_column"), "sample_id")
  group_col <- cfg_get(cfg, c("single_cell", "group_column"), "rcb_group")
  p1 <- plot_umap(object, celltype_col, title = "Major lineages")
  p2 <- plot_celltype_proportion(object, sample_col, group_col, celltype_col)
  save_plot_pdf(p1, mumoas_output(cfg, "single_cell", "01_lineage_umap.pdf"), 7, 6)
  save_plot_pdf(p2, mumoas_output(cfg, "single_cell", "01_celltype_proportions.pdf"), 10, 4)
  epithelial_labels <- cfg_get(cfg, c("single_cell", "epithelial_labels"), "Epithelial")
  epi <- subcluster_cells(object, cfg, epithelial_labels, "epithelial")
  save_plot_pdf(plot_umap(epi, "seurat_clusters", title = "Epithelial subclusters"), mumoas_output(cfg, "single_cell", "01_epithelial_umap.pdf"), 7, 6)
  invisible(list(discovery = object, epithelial = epi))
}

run_02_marker_genes <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "discovery")
  celltype_col <- cfg_get(cfg, c("single_cell", "celltype_column"), "celltype")
  features <- marker_sets(cfg, "discovery")
  p <- plot_marker_dot(object, features, celltype_col)
  save_plot_pdf(p, mumoas_output(cfg, "single_cell", "02_marker_dotplot.pdf"), 10, 4)
  invisible(p)
}

run_03_diagonal_volcano <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "discovery")
  cmp <- cfg_get(cfg, c("single_cell", "volcano"), list())
  result <- diagonal_volcano(object, cmp$subset_col %||% cfg_get(cfg, c("single_cell", "celltype_column"), "celltype"), cmp$subset_value %||% "Epithelial", cmp$group_col %||% cfg_get(cfg, c("single_cell", "group_column"), "rcb_group"), cmp$group_low %||% "RCB0", cmp$group_high %||% "RCB3", cmp$p_cut %||% 0.05, cmp$logfc_cut %||% 0.5)
  write_table_auto(result$table, mumoas_output(cfg, "single_cell", "03_diagonal_volcano.tsv"))
  save_plot_pdf(result$plot, mumoas_output(cfg, "single_cell", "03_diagonal_volcano.pdf"), 5, 4)
  invisible(result)
}

run_04_infercnv_cnv_score <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "discovery")
  inf <- cfg_get(cfg, c("single_cell", "infercnv"), list())
  infercnv_rds <- mumoas_path(cfg, inf$infercnv_rds %||% "outputs/r_analysis/single_cell/infercnv/run.final.infercnv_obj")
  if (isTRUE(inf$run)) {
    mumoas_pkg(c("infercnv", "rjags"))
    ann <- data.frame(cell = colnames(object), group = object@meta.data[[inf$annotation_column %||% cfg_get(cfg, c("single_cell", "celltype_column"), "celltype")]])
    ann_path <- mumoas_output(cfg, "single_cell", "infercnv", "annotations.tsv")
    utils::write.table(ann, ann_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = get_assay_matrix(object, slot = "counts"), annotations_file = ann_path, gene_order_file = mumoas_path(cfg, inf$gene_order_file), delim = "\t", ref_group_names = inf$ref_group_names)
    obj <- infercnv::run(obj, cutoff = inf$cutoff %||% 0.1, out_dir = mumoas_output(cfg, "single_cell", "infercnv"), cluster_by_groups = TRUE, denoise = TRUE, HMM = isTRUE(inf$HMM), no_prelim_plot = TRUE, write_expr_matrix = TRUE, num_threads = inf$num_threads %||% 4)
    saveRDS(obj, infercnv_rds)
  }
  if (!file.exists(infercnv_rds)) stop("inferCNV object not found; set single_cell.infercnv.run true or provide infercnv_rds", call. = FALSE)
  cnv <- cnv_score_table(readRDS(infercnv_rds), object)
  write_table_auto(cnv, mumoas_output(cfg, "single_cell", "04_cnv_scores.tsv"))
  save_plot_pdf(plot_cnv_scores(cnv, cfg_get(cfg, c("single_cell", "celltype_column"), "celltype")), mumoas_output(cfg, "single_cell", "04_cnv_score_by_celltype.pdf"), 8, 4)
  invisible(cnv)
}

run_05_odds_ratio <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "discovery")
  meta <- object@meta.data
  feature_col <- cfg_get(cfg, c("single_cell", "or_feature_column"), "seurat_clusters")
  group_col <- cfg_get(cfg, c("single_cell", "group_column"), "rcb_group")
  tab <- odds_ratio_table(meta, feature_col, group_col)
  write_table_auto(tab, mumoas_output(cfg, "single_cell", "05_cluster_group_odds_ratio.tsv"))
  save_plot_pdf(plot_or_heatmap(tab), mumoas_output(cfg, "single_cell", "05_cluster_group_odds_ratio_heatmap.pdf"), 6, 10)
  invisible(tab)
}

run_06_geneNMF_metaprograms <- function(cfg) {
  mumoas_pkg(c("GeneNMF", "Seurat", "ComplexHeatmap", "viridis"))
  object <- load_or_prepare_seurat(cfg, "discovery")
  epithelial_labels <- cfg_get(cfg, c("single_cell", "epithelial_labels"), "Epithelial")
  celltype_col <- cfg_get(cfg, c("single_cell", "celltype_column"), "celltype")
  epi <- subset(object, cells = rownames(object@meta.data)[object@meta.data[[celltype_col]] %in% epithelial_labels])
  sample_col <- cfg_get(cfg, c("single_cell", "sample_id_column"), "sample_id")
  nmf <- cfg_get(cfg, c("single_cell", "geneNMF"), list())
  programs <- GeneNMF::multiNMF(Seurat::SplitObject(epi, split.by = sample_col), assay = nmf$assay %||% "RNA", k = unlist(nmf$k_range %||% 4:12), min.exp = nmf$min_exp %||% 0.05, nfeatures = nmf$nfeatures %||% 2000)
  mp <- GeneNMF::getMetaPrograms(programs, metric = nmf$metric %||% "jaccard", weight.explained = nmf$weight_explained %||% 0.7, min.confidence = nmf$min_confidence %||% 0.1, specificity.weight = nmf$specificity_weight %||% 1, nMP = nmf$nMP %||% 6, max.genes = nmf$max_genes %||% 100)
  save(mp, file = mumoas_output(cfg, "single_cell", "06_geneNMF_metaprograms.RData"))
  genes <- data.frame(gene = unlist(mp$metaprograms.genes), program = rep(names(mp$metaprograms.genes), lengths(mp$metaprograms.genes)))
  write_table_auto(genes, mumoas_output(cfg, "single_cell", "06_metaprogram_genes.tsv"))
  annotation <- cfg_get(cfg, c("single_cell", "mp_annotations"), list())
  if (length(annotation) > 0) write_table_auto(do.call(rbind, lapply(names(annotation), function(x) data.frame(program = x, pathway = annotation[[x]]$pathway %||% NA_character_, markers = paste(annotation[[x]]$markers %||% character(), collapse = ";")))), mumoas_output(cfg, "single_cell", "06_metaprogram_annotations.tsv"))
  if (exists("runGSEA", mode = "function")) {
    top_p <- lapply(mp$metaprograms.genes, function(program) runGSEA(program, universe = rownames(epi), category = "H"))
    save(top_p, file = mumoas_output(cfg, "single_cell", "06_metaprogram_hallmark_gsea.RData"))
    top_p_table <- top_p_to_overlap_table(top_p)
    if (nrow(top_p_table) > 0) write_table_auto(top_p_table, mumoas_output(cfg, "single_cell", "06_metaprogram_hallmark_gsea.tsv"))
  }
  if (!is.null(mp$programs.similarity)) {
    grDevices::pdf(mumoas_output(cfg, "single_cell", "06_metaprogram_similarity_heatmap.pdf"), width = 7, height = 7)
    ComplexHeatmap::Heatmap(mp$programs.similarity, name = "Jaccard", col = viridis::viridis(100, option = "A", direction = -1), show_row_names = FALSE, show_column_names = FALSE)
    grDevices::dev.off()
  }
  invisible(mp)
}

run_07_kegg_enrichment <- function(cfg) {
  gene_sets <- config_gene_sets(cfg)
  tab <- kegg_enrichment(gene_sets)
  write_table_auto(tab, mumoas_output(cfg, "multiomics", "07_kegg_enrichment.tsv"))
  for (program in unique(tab$program)) save_plot_pdf(plot_enrichment_bar(tab, program), mumoas_output(cfg, "multiomics", paste0("07_kegg_", program, ".pdf")), 7, 4)
  invisible(tab)
}

run_08_hallmark_gene_sets <- function(cfg) {
  gene_sets <- config_gene_sets(cfg)
  universe <- unique(unlist(gene_sets))
  top_p_path <- mumoas_path(cfg, cfg_get(cfg, c("single_cell", "geneNMF", "hallmark_gsea_rdata"), NULL) %||% "outputs/r_analysis/single_cell/06_metaprogram_hallmark_gsea.RData")
  tab <- if (!is.null(top_p_path) && file.exists(top_p_path)) top_p_to_overlap_table(load_rdata_object(top_p_path, "top_p")) else hallmark_overlap(gene_sets, universe)
  write_table_auto(tab, mumoas_output(cfg, "multiomics", "08_hallmark_overlap.tsv"))
  p <- if ("pathway_label" %in% colnames(tab)) plot_hallmark_overlap_heatmap(tab) else plot_enrichment_heatmap(tab, "gene_ratio")
  save_plot_pdf(p, mumoas_output(cfg, "multiomics", "08_hallmark_overlap_heatmap.pdf"), 12, 4)
  invisible(tab)
}

run_09_mp_rcb_signatures <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "discovery")
  gene_sets <- config_gene_sets(cfg)
  object <- score_seurat_gene_sets(object, gene_sets)
  object <- apply_response_aucell(object, cfg)
  meta <- object@meta.data
  group_col <- cfg_get(cfg, c("single_cell", "group_column"), "rcb_group")
  thresholds <- attr(object, "mumoas_response_thresholds")
  if (!is.null(thresholds) && nrow(thresholds) > 0) {
    write_table_auto(thresholds, mumoas_output(cfg, "single_cell", "09_mp_rcb_aucell_thresholds.tsv"))
    for (i in seq_len(nrow(thresholds))) {
      program <- thresholds$program[[i]]
      p <- plot_auc_histogram(object@meta.data, program, group_col, thresholds$positive_group[[i]], thresholds$threshold[[i]])
      save_plot_pdf(p, mumoas_output(cfg, "single_cell", paste0("09_", program, "_aucell_histogram.pdf")), 5, 4)
    }
  }
  saveRDS(object, mumoas_output(cfg, "single_cell", "09_mp_rcb_scored_seurat.rds"))
  scores <- c(names(gene_sets), names(response_programs(cfg)))
  group_cols <- paste0(names(response_programs(cfg)), "_group")
  write_table_auto(cbind(cell = rownames(object@meta.data), object@meta.data[, intersect(c(scores, group_cols), colnames(object@meta.data)), drop = FALSE]), mumoas_output(cfg, "single_cell", "09_cell_program_scores.tsv"))
  save_plot_pdf(plot_score_density_umap(object, intersect(scores, colnames(object@meta.data))), mumoas_output(cfg, "single_cell", "09_program_feature_maps.pdf"), 11, 7)
  invisible(object)
}

run_10_aucell <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "discovery")
  gene_sets <- c(config_gene_sets(cfg), derived_response_gene_sets(cfg))
  auc <- compute_aucell_scores(object, gene_sets)
  write_table_auto(auc, mumoas_output(cfg, "single_cell", "10_aucell_scores.tsv"))
  invisible(auc)
}

run_11_rcb_group_scores <- function(cfg) {
  object <- run_09_mp_rcb_signatures(cfg)
  sample_col <- cfg_get(cfg, c("single_cell", "sample_id_column"), "sample_id")
  group_col <- cfg_get(cfg, c("single_cell", "group_column"), "rcb_group")
  scores <- c(names(config_gene_sets(cfg)), names(response_programs(cfg)))
  md <- object@meta.data
  sample_scores <- stats::aggregate(md[, intersect(scores, colnames(md)), drop = FALSE], by = list(sample_id = md[[sample_col]], group = md[[group_col]]), FUN = mean, na.rm = TRUE)
  write_table_auto(sample_scores, mumoas_output(cfg, "single_cell", "11_sample_program_scores.tsv"))
  save_plot_pdf(plot_score_violin(md, intersect(scores, colnames(md)), group_col), mumoas_output(cfg, "single_cell", "11_cell_program_violin_by_rcb.pdf"), 12, 7)
  invisible(sample_scores)
}

run_12_cellchat <- function(cfg) {
  mumoas_pkg(c("CellChat", "Seurat"))
  object <- run_09_mp_rcb_signatures(cfg)
  cc <- cfg_get(cfg, c("single_cell", "cellchat"), list())
  group_col <- cc$group_column %||% paste0(names(response_programs(cfg))[1], "_group")
  data_input <- as.matrix(get_assay_matrix(object, slot = "data"))
  meta <- data.frame(labels = object@meta.data[[group_col]], row.names = rownames(object@meta.data))
  chat <- CellChat::createCellChat(object = data_input, meta = meta, group.by = "labels")
  CellChat::CellChatDB.use <- CellChat::CellChatDB.human
  chat@DB <- CellChat::CellChatDB.use
  chat <- CellChat::subsetData(chat)
  chat <- CellChat::identifyOverExpressedGenes(chat)
  chat <- CellChat::identifyOverExpressedInteractions(chat)
  chat <- CellChat::computeCommunProb(chat)
  chat <- CellChat::filterCommunication(chat, min.cells = cc$min_cells %||% 10)
  chat <- CellChat::computeCommunProbPathway(chat)
  chat <- CellChat::aggregateNet(chat)
  saveRDS(chat, mumoas_output(cfg, "single_cell", "12_cellchat.rds"))
  grDevices::pdf(mumoas_output(cfg, "single_cell", "12_cellchat_interaction_count.pdf"), width = 7, height = 7)
  CellChat::netVisual_circle(chat@net$count, vertex.weight = as.numeric(table(chat@idents)), weight.scale = TRUE, label.edge = FALSE)
  grDevices::dev.off()
  sources <- unlist(cc$source_cell_states %||% unique(as.character(chat@idents)))
  for (source in sources) {
    save_plot_pdf(cellchat_source_network_plot(chat, source_celltype = source, count = TRUE, show_interaction = TRUE, flipped = TRUE), mumoas_output(cfg, "single_cell", paste0("12_cellchat_source_", gsub("[^A-Za-z0-9]+", "_", source), ".pdf")), 4.2, 3.6)
  }
  if (!is.null(cc$ligand_receptor_pairs)) {
    for (pair in unlist(cc$ligand_receptor_pairs)) {
      grDevices::pdf(mumoas_output(cfg, "single_cell", paste0("12_cellchat_lr_", gsub("[^A-Za-z0-9]+", "_", pair), ".pdf")), width = 4.5, height = 4.5)
      plot_cellchat_lr_chord(chat, pairs = pair, sources = cc$lr_sources %||% NULL, targets = cc$lr_targets %||% NULL)
      grDevices::dev.off()
    }
  }
  save_plot_pdf(plot_cellchat_role_dot(chat, "outgoing"), mumoas_output(cfg, "single_cell", "12_cellchat_outgoing_patterns.pdf"), 12, 4)
  save_plot_pdf(plot_cellchat_role_dot(chat, "incoming"), mumoas_output(cfg, "single_cell", "12_cellchat_incoming_patterns.pdf"), 12, 4)
  invisible(chat)
}

run_13_gsea_visualization <- function(cfg) {
  mo <- cfg_get(cfg, c("multiomics"), list())
  meta <- read_table_auto(mumoas_path(cfg, mo$sample_metadata))
  assays <- cfg_get(cfg, c("multiomics", "assays"), list())
  pathway_sets <- msigdb_hallmark_sets()
  selected_pathways <- unlist(mo$gsea_pathways %||% names(pathway_sets)[seq_len(min(6, length(pathway_sets)))])
  results <- list()
  for (assay in names(assays)) {
    mat <- log2_if_needed(read_expression_matrix(mumoas_path(cfg, assays[[assay]]$matrix)))
    mat <- combat_matrix(mat, meta, mo$sample_id_column %||% "sample_id", mo$batch_column)
    sample_col <- mo$sample_id_column %||% "sample_id"
    group_col <- assays[[assay]]$group_column %||% mo$group_column %||% "response_group"
    contrast <- assays[[assay]]$contrast %||% "RCB3-RCB0_1"
    rank <- rank_features_limma(mat, meta, sample_col, group_col, contrast)
    write_table_auto(rank, mumoas_output(cfg, "multiomics", paste0("13_", assay, "_limma_rank.tsv")))
    gsva_method <- assays[[assay]]$gsva_method %||% "ssgsea"
    gsva_scores <- gsva_score_matrix(mat, pathway_sets, method = gsva_method)
    gsva_limma <- limma_table_for_scores(gsva_scores, meta, sample_col, group_col, contrast)
    write_table_auto(data.frame(pathway = rownames(gsva_scores), gsva_scores, check.names = FALSE), mumoas_output(cfg, "multiomics", paste0("13_", assay, "_", gsva_method, "_scores.tsv")))
    write_table_auto(gsva_limma, mumoas_output(cfg, "multiomics", paste0("13_", assay, "_", gsva_method, "_limma.tsv")))
    save_plot_pdf(plot_gsva_bar(gsva_limma), mumoas_output(cfg, "multiomics", paste0("13_", assay, "_", gsva_method, "_barplot.pdf")), 8, 7)
    fgsea_tab <- fgsea_results(rank, pathway_sets, score_col = "t")
    write_table_auto(fgsea_tab, mumoas_output(cfg, "multiomics", paste0("13_", assay, "_fgsea.tsv")))
    save_plot_pdf(plot_gsea_panel(rank, pathway_sets, intersect(selected_pathways, names(pathway_sets)), title = paste(assay, "GSEA"), score_col = "t"), mumoas_output(cfg, "multiomics", paste0("13_", assay, "_gsea_panel.pdf")), 7, max(3, length(intersect(selected_pathways, names(pathway_sets))) * 1.4))
    results[[assay]] <- list(rank = rank, gsva_limma = gsva_limma, fgsea = fgsea_tab)
  }
  invisible(results)
}

run_14_model_external_validation_stats <- function(cfg) {
  mumoas_pkg(c("pROC", "ggplot2", "dplyr"))
  mv <- cfg_get(cfg, c("model_validation"), list())
  pred <- read_table_auto(mumoas_path(cfg, mv$predictions))
  truth_col <- mv$truth_column %||% "truth"
  class_cols <- grep("^prob_", colnames(pred), value = TRUE)
  classes <- sub("^prob_", "", class_cols)
  auc <- data.frame(class = classes, auc = vapply(seq_along(classes), function(i) as.numeric(pROC::auc(pROC::roc(pred[[truth_col]] == classes[[i]], pred[[class_cols[[i]]]], quiet = TRUE))), numeric(1)))
  pred$predicted <- classes[max.col(pred[, class_cols, drop = FALSE])]
  confusion <- as.data.frame.matrix(table(pred[[truth_col]], pred$predicted))
  confusion$truth <- rownames(confusion)
  write_table_auto(auc, mumoas_output(cfg, "model_validation", "14_auc.tsv"))
  write_table_auto(confusion, mumoas_output(cfg, "model_validation", "14_confusion_matrix.tsv"))
  calibration <- lapply(class_cols, function(col) {
    df <- data.frame(prob = pred[[col]], obs = as.integer(pred[[truth_col]] == sub("^prob_", "", col)))
    df$bin <- cut(df$prob, breaks = seq(0, 1, length.out = 11), include.lowest = TRUE)
    out <- df |>
      dplyr::group_by(.data$bin) |>
      dplyr::summarise(mean_predicted = mean(.data$prob), observed_rate = mean(.data$obs), n = dplyr::n(), .groups = "drop")
    out$class <- sub("^prob_", "", col)
    out
  })
  calibration <- do.call(rbind, calibration)
  write_table_auto(calibration, mumoas_output(cfg, "model_validation", "14_calibration.tsv"))
  invisible(list(auc = auc, confusion = confusion, calibration = calibration))
}

run_15_validation_scrna <- function(cfg) {
  object <- load_or_prepare_seurat(cfg, "validation")
  gene_sets <- config_gene_sets(cfg)
  object <- score_seurat_gene_sets(object, gene_sets)
  object <- apply_response_aucell(object, cfg)
  saveRDS(object, mumoas_output(cfg, "single_cell_validation", "15_validation_scored_seurat.rds"))
  celltype_col <- cfg_get(cfg, c("single_cell", "celltype_column"), "celltype")
  group_col <- cfg_get(cfg, c("single_cell", "validation_group_column"), cfg_get(cfg, c("single_cell", "group_column"), "rcb_group"))
  save_plot_pdf(plot_umap(object, celltype_col, title = "Validation scRNA-seq"), mumoas_output(cfg, "single_cell_validation", "15_validation_umap.pdf"), 7, 6)
  save_plot_pdf(plot_score_violin(object@meta.data, intersect(c(names(gene_sets), names(response_programs(cfg))), colnames(object@meta.data)), group_col), mumoas_output(cfg, "single_cell_validation", "15_validation_program_violin.pdf"), 12, 7)
  save_plot_pdf(plot_score_density_umap(object, intersect(c(names(gene_sets), names(response_programs(cfg))), colnames(object@meta.data))), mumoas_output(cfg, "single_cell_validation", "15_validation_feature_maps.pdf"), 12, 7)
  invisible(object)
}

run_16_mp_score_violin <- function(cfg) {
  object <- run_09_mp_rcb_signatures(cfg)
  group_col <- cfg_get(cfg, c("single_cell", "group_column"), "rcb_group")
  scores <- intersect(c(names(config_gene_sets(cfg)), names(response_programs(cfg))), colnames(object@meta.data))
  p <- plot_score_violin(object@meta.data, scores, group_col)
  save_plot_pdf(p, mumoas_output(cfg, "single_cell", "16_mp_score_violin.pdf"), 12, 7)
  invisible(p)
}

run_17_rcb_figures <- function(cfg) {
  mo <- cfg_get(cfg, c("multiomics"), list())
  meta <- read_table_auto(mumoas_path(cfg, mo$sample_metadata))
  assays <- cfg_get(cfg, c("multiomics", "assays"), list())
  gene_sets <- config_gene_sets(cfg)
  score_sets <- c(gene_sets, derived_response_gene_sets(cfg, gene_sets))
  outputs <- list()
  for (assay in names(assays)) {
    mat <- log2_if_needed(read_expression_matrix(mumoas_path(cfg, assays[[assay]]$matrix)))
    mat <- combat_matrix(mat, meta, mo$sample_id_column %||% "sample_id", mo$batch_column)
    scores <- matrix_gene_set_scores(mat, score_sets)
    scores <- merge(scores, meta, by.x = "sample_id", by.y = mo$sample_id_column %||% "sample_id")
    write_table_auto(scores, mumoas_output(cfg, "multiomics", paste0("17_", assay, "_program_scores.tsv")))
    group_cols <- unique(c(assays[[assay]]$group_column, assays[[assay]]$validation_group_columns, mo$group_column %||% "response_group"))
    for (group_col in group_cols) {
      if (!group_col %in% colnames(scores)) next
      p <- plot_score_violin(scores, intersect(c(names(score_sets), names(response_programs(cfg))), colnames(scores)), group_col)
      save_plot_pdf(p, mumoas_output(cfg, "multiomics", paste0("17_", assay, "_", group_col, "_score_violin.pdf")), 12, 7)
    }
    outputs[[assay]] <- scores
  }
  invisible(outputs)
}

run_18_enrichment_heatmap <- function(cfg) {
  gene_sets <- config_gene_sets(cfg)
  top_p_path <- mumoas_path(cfg, cfg_get(cfg, c("single_cell", "geneNMF", "hallmark_gsea_rdata"), NULL) %||% "outputs/r_analysis/single_cell/06_metaprogram_hallmark_gsea.RData")
  tab <- if (!is.null(top_p_path) && file.exists(top_p_path)) top_p_to_overlap_table(load_rdata_object(top_p_path, "top_p")) else hallmark_overlap(gene_sets, unique(unlist(gene_sets)))
  write_table_auto(tab, mumoas_output(cfg, "multiomics", "18_hallmark_enrichment_for_heatmap.tsv"))
  p <- if ("pathway_label" %in% colnames(tab)) plot_hallmark_overlap_heatmap(tab) else plot_enrichment_heatmap(tab, "gene_ratio")
  save_plot_pdf(p, mumoas_output(cfg, "multiomics", "18_hallmark_enrichment_heatmap.pdf"), 12, 4)
  invisible(tab)
}

run_all_sc_multiomics <- function(cfg) {
  write_main_figure_panel_map(cfg)
  run_01_dimension_clustering_lineages(cfg)
  run_02_marker_genes(cfg)
  run_03_diagonal_volcano(cfg)
  if (isTRUE(cfg_get(cfg, c("single_cell", "infercnv", "enabled"), FALSE))) run_04_infercnv_cnv_score(cfg)
  run_05_odds_ratio(cfg)
  if (isTRUE(cfg_get(cfg, c("single_cell", "geneNMF", "enabled"), FALSE))) run_06_geneNMF_metaprograms(cfg)
  run_07_kegg_enrichment(cfg)
  run_08_hallmark_gene_sets(cfg)
  run_09_mp_rcb_signatures(cfg)
  run_10_aucell(cfg)
  run_11_rcb_group_scores(cfg)
  if (isTRUE(cfg_get(cfg, c("single_cell", "cellchat", "enabled"), FALSE))) run_12_cellchat(cfg)
  run_13_gsea_visualization(cfg)
  if (!is.null(cfg_get(cfg, c("model_validation", "predictions"), NULL))) run_14_model_external_validation_stats(cfg)
  if (!is.null(cfg_get(cfg, c("single_cell", "validation", "seurat_rds"), NULL)) || !is.null(cfg_get(cfg, c("single_cell", "validation", "sample_manifest"), NULL))) run_15_validation_scrna(cfg)
  run_16_mp_score_violin(cfg)
  run_17_rcb_figures(cfg)
  run_18_enrichment_heatmap(cfg)
  invisible(TRUE)
}
