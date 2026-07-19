suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

setClass("KeyMixin", contains = "VIRTUAL", slots = list(key = "character"))
setClass("Assay", contains = "KeyMixin",
         slots = c(counts = "AnyMatrix", data = "AnyMatrix", scale.data = "matrix",
                   assay.orig = "OptionalCharacter", var.features = "vector",
                   meta.features = "data.frame", misc = "OptionalList"))
setClass("Seurat",
         slots = c(assays = "list", meta.data = "data.frame", active.assay = "character",
                   active.ident = "factor", graphs = "list", neighbors = "list",
                   reductions = "list", images = "list", project.name = "character",
                   misc = "list", version = "package_version", commands = "list",
                   tools = "list"))

workspace_out <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"
plan_dir <- "D:/OC_spatiogenomics/infercnv/integrated_oc_plan_analysis"
tables_dir <- file.path(plan_dir, "tables")
figures_dir <- file.path(plan_dir, "figures")
scripts_dir <- file.path(plan_dir, "scripts")
source_dir <- file.path(plan_dir, "source_outputs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

message("Copying existing source outputs...")
files_to_copy <- c(
  "infercnv_cell_to_subclone_k5.csv",
  "infercnv_subcluster_to_subclone_k5.csv",
  "infercnv_subclone_cnv_regions_k5.csv",
  "infercnv_subclone_cnv_state_summary_k5.csv",
  "infercnv_subclone_heatmap.png",
  "infercnv_subclone_heatmap.pdf",
  "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv",
  "integrated_oc_Bcell_subtypes_mapped_from_integratedocBcells.csv",
  "integrated_oc_interaction_group_cell_counts.csv",
  "integrated_oc_mapping_summary.csv",
  "lr_group_pair_summary_connectome_like.csv",
  "lr_top_interactions_involving_cnv_subclones.csv",
  "lr_top_interactions_cnv_subclone_to_celltype.csv",
  "lr_top_interactions_celltype_to_cnv_subclone.csv",
  "lr_top_distinct_LR_pairs_involving_cnv_subclones.csv",
  "lr_group_pair_total_score_heatmap.png",
  "lr_group_pair_total_score_heatmap.pdf",
  "lr_cnv_subclone_celltype_directional_heatmaps.pdf",
  "lr_top40_subclone_involving_pairs.png",
  "map_subtypes_and_run_lr_interactions.R",
  "plot_infercnv_subclone_heatmap.R"
)
for (f in files_to_copy) {
  src <- file.path(workspace_out, f)
  if (file.exists(src)) file.copy(src, file.path(source_dir, f), overwrite = TRUE)
}

message("Loading integrated_oc...")
obj <- readRDS("D:/OC_spatiogenomics/infercnv/integrated_oc.RData")
expr <- obj@assays$RNA@data
meta0 <- obj@meta.data
meta0$cell_integrated_oc <- rownames(meta0)
meta <- read.csv(file.path(workspace_out, "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv"), stringsAsFactors = FALSE)
meta <- meta[match(rownames(meta0), meta$cell_integrated_oc), ]
if (!all(meta$cell_integrated_oc == rownames(meta0))) stop("Metadata order mismatch.")

message("Exporting metadata, PCA and UMAP coordinates...")
pca <- obj@reductions$pca@cell.embeddings[, 1:30, drop = FALSE]
umap <- obj@reductions$umap@cell.embeddings
write.csv(data.frame(cell = rownames(meta0), pca), file.path(tables_dir, "integrated_oc_pca30_for_knn.csv"), row.names = FALSE)
write.csv(data.frame(cell = rownames(meta0), umap, meta[, c("cell_type", "cnv_subclone", "interaction_group", "sample_type", "seurat_clusters")]),
          file.path(tables_dir, "integrated_oc_umap_with_cnv_and_interaction_groups.csv"), row.names = FALSE)
write.csv(meta, file.path(tables_dir, "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv"), row.names = FALSE)

message("Summarising clone composition...")
clone_cells <- !is.na(meta$cnv_subclone)
clone_counts <- as.data.frame(table(meta$cnv_subclone[clone_cells]), stringsAsFactors = FALSE)
colnames(clone_counts) <- c("cnv_subclone", "n_cells")
clone_counts$fraction_of_cnv_cells <- clone_counts$n_cells / sum(clone_counts$n_cells)
write.csv(clone_counts, file.path(tables_dir, "cnv_subclone_cell_counts.csv"), row.names = FALSE)

sample_clone <- as.data.frame(table(meta$sample_type[clone_cells], meta$cnv_subclone[clone_cells]), stringsAsFactors = FALSE)
colnames(sample_clone) <- c("sample_type", "cnv_subclone", "n_cells")
sample_clone$total_in_sample <- ave(sample_clone$n_cells, sample_clone$sample_type, FUN = sum)
sample_clone$fraction_in_sample <- ifelse(sample_clone$total_in_sample > 0, sample_clone$n_cells / sample_clone$total_in_sample, NA)
write.csv(sample_clone, file.path(tables_dir, "cnv_subclone_composition_by_sample_type.csv"), row.names = FALSE)

subcluster_clone <- read.csv(file.path(workspace_out, "infercnv_subcluster_to_subclone_k5.csv"), stringsAsFactors = FALSE)
write.csv(subcluster_clone, file.path(tables_dir, "infercnv_subcluster_to_subclone_k5.csv"), row.names = FALSE)
cnv_regions <- read.csv(file.path(workspace_out, "infercnv_subclone_cnv_regions_k5.csv"), stringsAsFactors = FALSE)
region_summary <- aggregate(list(n_regions = cnv_regions$cnv_name),
                            by = list(cnv_subclone = cnv_regions$cnv_subclone, chr = cnv_regions$chr, event = cnv_regions$event),
                            FUN = length)
write.csv(region_summary, file.path(tables_dir, "cnv_subclone_cnv_event_counts_by_chr.csv"), row.names = FALSE)

message("Computing functional module scores for CNV subclones...")
gene_sets <- list(
  EMT = c("VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH2", "ITGA5", "COL1A1", "COL1A2", "SPARC", "TAGLN", "ACTA2"),
  Hypoxia = c("VEGFA", "CA9", "EGLN3", "BNIP3", "LDHA", "SLC2A1", "PDK1", "NDRG1", "ADM", "ENO1", "ALDOA"),
  KRAS_activation = c("DUSP6", "FOS", "JUN", "EGR1", "ETV4", "ETV5", "SPRY2", "SPRY4", "CCND1", "MYC"),
  Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "UBE2C", "BIRC5"),
  Stemness_epithelial = c("PROM1", "ALDH1A1", "SOX2", "NANOG", "EPCAM", "KRT8", "KRT18", "KRT19", "MUC16"),
  Immune_modulatory = c("MIF", "CD274", "PDCD1LG2", "LGALS9", "CD47", "TGFB1", "HLA-A", "HLA-B", "HLA-C", "B2M"),
  M2_related = c("CSF1", "CCL2", "CCL5", "SPP1", "TGFB1", "MIF", "APOE", "GAS6"),
  LAPTM5_axis = c("LAPTM5")
)
score_mat <- matrix(NA_real_, nrow = ncol(expr), ncol = length(gene_sets), dimnames = list(colnames(expr), names(gene_sets)))
for (nm in names(gene_sets)) {
  genes <- intersect(gene_sets[[nm]], rownames(expr))
  if (length(genes) > 0) {
    score_mat[, nm] <- Matrix::colMeans(expr[genes, , drop = FALSE])
  }
}
score_df <- data.frame(cell_integrated_oc = rownames(meta0), score_mat[rownames(meta0), , drop = FALSE], stringsAsFactors = FALSE)
score_df$cnv_subclone <- meta$cnv_subclone
score_df$interaction_group <- meta$interaction_group
write.csv(score_df, file.path(tables_dir, "cell_level_function_module_scores.csv"), row.names = FALSE)

clone_score_summary <- do.call(rbind, lapply(names(gene_sets), function(score) {
  x <- score_df[!is.na(score_df$cnv_subclone), c("cnv_subclone", score)]
  colnames(x) <- c("cnv_subclone", "score")
  out <- aggregate(score ~ cnv_subclone, data = x, FUN = function(v) c(mean = mean(v, na.rm = TRUE), median = median(v, na.rm = TRUE)))
  data.frame(cnv_subclone = out$cnv_subclone, module = score, mean_score = out$score[, "mean"], median_score = out$score[, "median"])
}))
clone_score_summary <- clone_score_summary[order(clone_score_summary$module, -clone_score_summary$mean_score), ]
write.csv(clone_score_summary, file.path(tables_dir, "cnv_subclone_function_module_score_summary.csv"), row.names = FALSE)

key_genes <- intersect(c("LAPTM5", "KRAS", "MIF", "SPP1", "CSF1", "CCL2", "TGFB1", "VEGFA", "CD47", "EPCAM", "KRT8", "VIM", "FN1", "MKI67"), rownames(expr))
key_expr <- t(as.matrix(expr[key_genes, rownames(meta0)[clone_cells], drop = FALSE]))
key_expr_df <- data.frame(cnv_subclone = meta$cnv_subclone[clone_cells], key_expr, check.names = FALSE)
key_summary <- do.call(rbind, lapply(key_genes, function(g) {
  out <- aggregate(key_expr_df[[g]], by = list(cnv_subclone = key_expr_df$cnv_subclone), FUN = function(v) c(mean = mean(v), pct = mean(v > 0)))
  data.frame(cnv_subclone = out$cnv_subclone, gene = g, mean_expr = out$x[, "mean"], pct_expr = out$x[, "pct"])
}))
write.csv(key_summary, file.path(tables_dir, "cnv_subclone_key_gene_expression_summary.csv"), row.names = FALSE)

message("Plotting functional score heatmap...")
wide <- reshape(clone_score_summary[, c("cnv_subclone", "module", "mean_score")], idvar = "cnv_subclone", timevar = "module", direction = "wide")
rownames(wide) <- wide$cnv_subclone
wide$cnv_subclone <- NULL
colnames(wide) <- sub("^mean_score\\.", "", colnames(wide))
scaled <- scale(as.matrix(wide))
pdf(file.path(figures_dir, "cnv_subclone_function_module_score_heatmap.pdf"), width = 8, height = 5)
Heatmap(t(scaled), name = "z-score", col = colorRamp2(c(-1.5, 0, 1.5), c("#2166AC", "white", "#B2182B")),
        cluster_rows = TRUE, cluster_columns = TRUE, column_title = "CNV subclones", row_title = "Functional modules")
dev.off()
png(file.path(figures_dir, "cnv_subclone_function_module_score_heatmap.png"), width = 1600, height = 1000, res = 200)
Heatmap(t(scaled), name = "z-score", col = colorRamp2(c(-1.5, 0, 1.5), c("#2166AC", "white", "#B2182B")),
        cluster_rows = TRUE, cluster_columns = TRUE, column_title = "CNV subclones", row_title = "Functional modules")
dev.off()

message("Plotting clone composition...")
p <- ggplot(sample_clone, aes(x = sample_type, y = fraction_in_sample, fill = cnv_subclone)) +
  geom_col(color = "grey30", size = 0.1) +
  theme_bw(base_size = 11) +
  labs(x = "sample_type", y = "Fraction within CNV cells", fill = "CNV subclone",
       title = "CNV subclone composition by sample type")
ggsave(file.path(figures_dir, "cnv_subclone_composition_by_sample_type.png"), p, width = 7, height = 4.5, dpi = 220)
ggsave(file.path(figures_dir, "cnv_subclone_composition_by_sample_type.pdf"), p, width = 7, height = 4.5)

message("Preparing focused LR axis tables...")
lr_all <- read.csv(file.path(workspace_out, "lr_top_interactions_involving_cnv_subclones.csv"), stringsAsFactors = FALSE)
focus_axes <- list(
  MIF_CD74_CXCR4 = list(lig = c("MIF"), rec = c("CD74", "CXCR4")),
  SPP1_CD44_integrin = list(lig = c("SPP1"), rec = c("CD44", "ITGAV", "ITGB1", "ITGB5")),
  CSF1_CSF1R = list(lig = c("CSF1"), rec = c("CSF1R")),
  CCL2_CCR2 = list(lig = c("CCL2"), rec = c("CCR2")),
  CCL5_CCR5 = list(lig = c("CCL5"), rec = c("CCR5")),
  TGFB_TGFBR = list(lig = c("TGFB1", "TGFB2", "TGFB3"), rec = c("TGFBR1", "TGFBR2", "TGFBR3")),
  IL6_IL6R = list(lig = c("IL6"), rec = c("IL6R", "IL6ST")),
  GAS6_AXL = list(lig = c("GAS6"), rec = c("AXL")),
  CD47_SIRPA = list(lig = c("CD47"), rec = c("SIRPA")),
  APOE_LRP1 = list(lig = c("APOE"), rec = c("LRP1")),
  VEGFA_VEGFR = list(lig = c("VEGFA"), rec = c("KDR", "FLT1", "FLT4")),
  CXCL12_CXCR4 = list(lig = c("CXCL12"), rec = c("CXCR4"))
)
focus <- do.call(rbind, lapply(names(focus_axes), function(axis) {
  z <- focus_axes[[axis]]
  y <- lr_all[lr_all$ligand %in% z$lig & lr_all$receptor %in% z$rec, , drop = FALSE]
  if (nrow(y) == 0) return(NULL)
  y$focus_axis <- axis
  y
}))
if (!is.null(focus)) {
  focus <- focus[order(focus$focus_axis, -focus$score), ]
  write.csv(focus, file.path(tables_dir, "focused_LR_axes_involving_cnv_subclones.csv"), row.names = FALSE)
}

message("Done prepare script.")
