suppressPackageStartupMessages({
  library(Matrix)
  library(matrixStats)
  library(limma)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(corrplot)
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

base_dir <- "D:/OC_spatiogenomics/infercnv/CNV_expression_joint_analysis_integrated_oc"
tables_dir <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")
scripts_dir <- file.path(base_dir, "scripts")
source_dir <- file.path(base_dir, "source")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

workspace_out <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"
prior_dir <- "D:/OC_spatiogenomics/infercnv/integrated_oc_plan_analysis"
infer_dir <- "D:/OC_spatiogenomics/infercnv/infercnv_Other_vs_Immune_subcluster"

canonical_cell <- function(x) sub("^([^_]+)_\\1_", "\\1_", x)

message("Copying source files...")
copy_candidates <- c(
  file.path(workspace_out, "infercnv_cell_to_subclone_k5.csv"),
  file.path(workspace_out, "infercnv_subcluster_to_subclone_k5.csv"),
  file.path(workspace_out, "infercnv_subclone_cnv_regions_k5.csv"),
  file.path(workspace_out, "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv"),
  file.path(workspace_out, "lr_top_distinct_LR_pairs_involving_cnv_subclones.csv"),
  file.path(prior_dir, "integrated_oc_infercnv_plan_task_breakdown_and_results.md"),
  "D:/Downloads/CNV_expression_joint_analysis_integrated_oc (1).docx"
)
for (f in copy_candidates[file.exists(copy_candidates)]) file.copy(f, file.path(source_dir, basename(f)), overwrite = TRUE)

message("Loading integrated_oc and metadata...")
obj <- readRDS("D:/OC_spatiogenomics/infercnv/integrated_oc.RData")
expr <- obj@assays$RNA@data
counts <- obj@assays$RNA@counts
meta0 <- obj@meta.data
meta0$cell_integrated_oc <- rownames(meta0)
meta0$cell_canonical <- canonical_cell(rownames(meta0))

meta_file <- file.path(workspace_out, "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[match(rownames(meta0), meta$cell_integrated_oc), ]
if (!all(meta$cell_integrated_oc == rownames(meta0))) stop("metadata alignment failed")

cell_map <- read.csv(file.path(workspace_out, "infercnv_cell_to_subclone_k5.csv"), stringsAsFactors = FALSE)
cell_map$cell_canonical <- canonical_cell(cell_map$cell)
meta$infercnv_cell_group_name <- cell_map$cell_group_name[match(meta$cell_canonical, cell_map$cell_canonical)]
meta$sample_id <- sub("_.*$", "", meta$cell_canonical)
meta$CNV_clone <- meta$cnv_subclone
tumor_idx <- !is.na(meta$CNV_clone)
tumor_cells <- meta$cell_integrated_oc[tumor_idx]
message("Tumor/CNV cells: ", length(tumor_cells))
write.csv(meta, file.path(tables_dir, "integrated_oc_metadata_with_CNV_clone_and_subtypes.csv"), row.names = FALSE)

message("Constructing HMM gene-level CNV matrix from pred_cnv_genes.dat...")
genes_used <- read.delim(file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat"),
                         stringsAsFactors = FALSE, check.names = FALSE)
genes_used$gene <- rownames(genes_used)
colnames(genes_used)[1:3] <- c("chr", "start", "end")

pred_genes <- read.delim(file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"),
                         stringsAsFactors = FALSE)
pred_genes <- pred_genes[grepl("^observation\\.", pred_genes$cell_group_name), ]
obs_groups <- sort(unique(cell_map$cell_group_name))
gene_universe <- intersect(genes_used$gene, rownames(expr))
cnv_gene_state <- matrix(0, nrow = length(gene_universe), ncol = length(obs_groups),
                         dimnames = list(gene_universe, obs_groups))
pg <- pred_genes[pred_genes$gene %in% gene_universe & pred_genes$cell_group_name %in% obs_groups, ]
cnv_gene_state[cbind(pg$gene, pg$cell_group_name)] <- pg$state - 3
saveRDS(cnv_gene_state, file.path(tables_dir, "cnv_gene_state_matrix_gene_by_infercnv_subcluster.rds"))

message("CNV burden and chromosome scores...")
group_burden <- colMeans(abs(cnv_gene_state))
meta$CNV_burden <- group_burden[meta$infercnv_cell_group_name]
clone_burden <- aggregate(CNV_burden ~ CNV_clone, data = meta[tumor_idx, ], FUN = function(v) c(mean = mean(v), median = median(v), sd = sd(v)))
clone_burden <- data.frame(CNV_clone = clone_burden$CNV_clone,
                           mean = clone_burden$CNV_burden[, "mean"],
                           median = clone_burden$CNV_burden[, "median"],
                           sd = clone_burden$CNV_burden[, "sd"])
write.csv(clone_burden, file.path(tables_dir, "CNV_clone_burden_summary.csv"), row.names = FALSE)
write.csv(meta[, c("cell_integrated_oc", "cell_canonical", "CNV_clone", "infercnv_cell_group_name", "CNV_burden", "sample_id", "sample_type", "interaction_group")],
          file.path(tables_dir, "cell_level_CNV_burden_and_clone_metadata.csv"), row.names = FALSE)

go <- genes_used[match(rownames(cnv_gene_state), genes_used$gene), ]
chr_levels <- unique(go$chr)
chr_group_score <- sapply(chr_levels, function(chr) {
  genes <- rownames(cnv_gene_state)[go$chr == chr]
  colMeans(cnv_gene_state[genes, , drop = FALSE])
})
chr_group_score <- t(chr_group_score)
clone_levels <- sort(unique(meta$CNV_clone[tumor_idx]))
chr_clone_score <- sapply(clone_levels, function(cl) {
  grps <- intersect(unique(meta$infercnv_cell_group_name[meta$CNV_clone == cl]), colnames(chr_group_score))
  rowMeans(chr_group_score[, grps, drop = FALSE])
})
write.csv(chr_clone_score, file.path(tables_dir, "chromosome_level_CNV_score_by_clone.csv"))

pdf(file.path(figures_dir, "chromosome_level_CNV_score_by_clone_heatmap.pdf"), width = 7.5, height = 8)
Heatmap(chr_clone_score, name = "CNV score", col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        cluster_rows = FALSE, cluster_columns = TRUE, row_title = "Chromosome", column_title = "CNV clone")
dev.off()
png(file.path(figures_dir, "chromosome_level_CNV_score_by_clone_heatmap.png"), width = 1500, height = 1600, res = 220)
Heatmap(chr_clone_score, name = "CNV score", col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        cluster_rows = FALSE, cluster_columns = TRUE, row_title = "Chromosome", column_title = "CNV clone")
dev.off()

p_burden <- ggplot(meta[tumor_idx, ], aes(x = CNV_clone, y = CNV_burden, fill = CNV_clone)) +
  geom_violin(scale = "width", color = "grey30", size = 0.15) +
  geom_boxplot(width = 0.14, outlier.size = 0.2) +
  theme_bw(base_size = 11) +
  labs(x = "CNV clone", y = "CNV burden: mean abs(HMM state - 3)", title = "CNV burden by clone") +
  theme(legend.position = "none")
ggsave(file.path(figures_dir, "CNV_burden_by_clone_violin.png"), p_burden, width = 7, height = 4.8, dpi = 220)
ggsave(file.path(figures_dir, "CNV_burden_by_clone_violin.pdf"), p_burden, width = 7, height = 4.8)

message("Signature and TF target activity scores...")
gene_sets <- list(
  EMT = c("VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH2", "ITGA5", "COL1A1", "COL1A2", "SPARC", "TAGLN", "ACTA2"),
  Hypoxia = c("VEGFA", "CA9", "EGLN3", "BNIP3", "LDHA", "SLC2A1", "PDK1", "NDRG1", "ADM", "ENO1", "ALDOA"),
  KRAS_UP = c("DUSP6", "FOS", "JUN", "EGR1", "ETV4", "ETV5", "SPRY2", "SPRY4", "CCND1", "MYC", "AREG", "EREG"),
  KRAS_DN = c("DUSP1", "TOB1", "BTG2", "KLF6", "JUNB", "FOSB"),
  LAPTM5_axis = c("LAPTM5"),
  Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "UBE2C", "BIRC5"),
  Immune_modulatory = c("MIF", "CD274", "PDCD1LG2", "LGALS9", "CD47", "TGFB1", "HLA-A", "HLA-B", "HLA-C", "B2M")
)
tf_sets <- list(
  SNAI1 = c("VIM", "FN1", "CDH2", "ZEB1", "ZEB2", "MMP2", "MMP9"),
  SNAI2 = c("VIM", "FN1", "CDH2", "ITGA5", "MMP14", "SERPINE1"),
  ZEB1 = c("VIM", "ITGA5", "COL1A1", "COL1A2", "SPARC"),
  HIF1A = c("VEGFA", "CA9", "LDHA", "SLC2A1", "PDK1", "BNIP3", "NDRG1"),
  EPAS1 = c("VEGFA", "EGLN3", "ADM", "NDRG1", "SLC2A1"),
  JUN = c("JUN", "FOS", "MMP1", "MMP3", "IL6", "DUSP1"),
  FOS = c("FOS", "JUN", "DUSP1", "EGR1", "IER2"),
  ATF3 = c("ATF3", "JUN", "FOS", "DUSP1", "EGR1"),
  EGR1 = c("EGR1", "FOS", "JUN", "DUSP1", "BTG2"),
  STAT3 = c("STAT3", "SOCS3", "JUNB", "BCL3", "IRF1", "IL6ST"),
  NFKB1 = c("NFKB1", "NFKBIA", "TNFAIP3", "BIRC3", "ICAM1", "CXCL8"),
  RELA = c("RELA", "NFKBIA", "TNFAIP3", "BIRC3", "ICAM1", "CXCL8"),
  MYC = c("MYC", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "ODC1"),
  TP53 = c("TP53", "CDKN1A", "MDM2", "BBC3", "BAX", "GADD45A"),
  FOXM1 = c("FOXM1", "MKI67", "TOP2A", "BIRC5", "UBE2C", "CCNB1", "CCNB2"),
  KLF4 = c("KLF4", "CDKN1A", "JUNB", "BTG2", "EGR1")
)
all_sets <- c(gene_sets, setNames(tf_sets, paste0("TF_", names(tf_sets))))
score_mat <- matrix(NA_real_, nrow = length(tumor_cells), ncol = length(all_sets), dimnames = list(tumor_cells, names(all_sets)))
for (nm in names(all_sets)) {
  genes <- intersect(all_sets[[nm]], rownames(expr))
  if (length(genes) > 0) score_mat[, nm] <- Matrix::colMeans(expr[genes, tumor_cells, drop = FALSE])
}
score_df <- data.frame(cell_integrated_oc = tumor_cells,
                       CNV_clone = meta$CNV_clone[tumor_idx],
                       CNV_burden = meta$CNV_burden[tumor_idx],
                       sample_id = meta$sample_id[tumor_idx],
                       score_mat, check.names = FALSE)
write.csv(score_df, file.path(tables_dir, "tumor_cell_CNV_burden_signature_TF_scores.csv"), row.names = FALSE)

clone_function_summary <- do.call(rbind, lapply(colnames(score_mat), function(module) {
  v <- score_df[, c("CNV_clone", module)]
  colnames(v) <- c("CNV_clone", "score")
  out <- aggregate(score ~ CNV_clone, data = v, FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
  data.frame(CNV_clone = out$CNV_clone, feature = module, mean_score = out$score[, "mean"], median_score = out$score[, "median"], sd_score = out$score[, "sd"])
}))
write.csv(clone_function_summary, file.path(tables_dir, "CNV_clone_functional_and_TF_activity_summary.csv"), row.names = FALSE)

feature_clone_mat <- reshape(clone_function_summary[, c("CNV_clone", "feature", "mean_score")], idvar = "CNV_clone", timevar = "feature", direction = "wide")
rownames(feature_clone_mat) <- feature_clone_mat$CNV_clone
feature_clone_mat$CNV_clone <- NULL
colnames(feature_clone_mat) <- sub("^mean_score\\.", "", colnames(feature_clone_mat))
feature_clone_mat <- as.matrix(feature_clone_mat)
pdf(file.path(figures_dir, "CNV_clone_signature_and_TF_activity_heatmap.pdf"), width = 9, height = 9)
Heatmap(t(scale(feature_clone_mat)), name = "z-score", col = colorRamp2(c(-1.5, 0, 1.5), c("#2166AC", "white", "#B2182B")),
        cluster_rows = TRUE, cluster_columns = TRUE, row_title = "Signature / TF target score", column_title = "CNV clone")
dev.off()
png(file.path(figures_dir, "CNV_clone_signature_and_TF_activity_heatmap.png"), width = 1800, height = 1800, res = 220)
Heatmap(t(scale(feature_clone_mat)), name = "z-score", col = colorRamp2(c(-1.5, 0, 1.5), c("#2166AC", "white", "#B2182B")),
        cluster_rows = TRUE, cluster_columns = TRUE, row_title = "Signature / TF target score", column_title = "CNV clone")
dev.off()

message("Clone marker analysis: one-vs-rest expression statistics...")
expr_tumor <- expr[, tumor_cells, drop = FALSE]
groups <- factor(meta$CNV_clone[tumor_idx], levels = clone_levels)
avg_by_clone <- sapply(clone_levels, function(cl) Matrix::rowMeans(expr_tumor[, groups == cl, drop = FALSE]))
pct_by_clone <- sapply(clone_levels, function(cl) Matrix::rowMeans(expr_tumor[, groups == cl, drop = FALSE] > 0))
avg_rest <- sapply(clone_levels, function(cl) Matrix::rowMeans(expr_tumor[, groups != cl, drop = FALSE]))
pct_rest <- sapply(clone_levels, function(cl) Matrix::rowMeans(expr_tumor[, groups != cl, drop = FALSE] > 0))
logfc <- log2(avg_by_clone + 0.1) - log2(avg_rest + 0.1)
marker_list <- list()
for (cl in clone_levels) {
  cand <- rownames(expr_tumor)[logfc[, cl] > 0.25 & pct_by_clone[, cl] >= 0.10]
  cand <- head(cand[order(logfc[cand, cl], decreasing = TRUE)], 2000)
  pvals <- vapply(cand, function(g) {
    suppressWarnings(wilcox.test(as.numeric(expr_tumor[g, groups == cl]), as.numeric(expr_tumor[g, groups != cl]))$p.value)
  }, numeric(1))
  marker_list[[cl]] <- data.frame(
    cluster = cl, gene = cand, avg_log2FC = logfc[cand, cl], pct.1 = pct_by_clone[cand, cl],
    pct.2 = pct_rest[cand, cl], p_val = pvals, p_val_adj = p.adjust(pvals, "BH"),
    stringsAsFactors = FALSE
  )
}
markers <- do.call(rbind, marker_list)
markers <- markers[order(markers$cluster, markers$p_val_adj, -markers$avg_log2FC), ]
write.csv(markers, file.path(tables_dir, "CNV_clone_DEG_single_cell_wilcoxon_one_vs_rest.csv"), row.names = FALSE)

top_markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(x) head(x[order(x$p_val_adj, -x$avg_log2FC), ], 20)))
write.csv(top_markers, file.path(tables_dir, "CNV_clone_top20_markers_per_clone.csv"), row.names = FALSE)
top_genes <- unique(top_markers$gene)
avg_top <- avg_by_clone[top_genes, , drop = FALSE]
pdf(file.path(figures_dir, "CNV_clone_top_marker_average_expression_heatmap.pdf"), width = 8, height = 12)
Heatmap(t(scale(t(avg_top))), name = "z-score", col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
        cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 5),
        column_title = "CNV clone", row_title = "Top marker genes")
dev.off()
png(file.path(figures_dir, "CNV_clone_top_marker_average_expression_heatmap.png"), width = 1600, height = 2400, res = 220)
Heatmap(t(scale(t(avg_top))), name = "z-score", col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
        cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 5),
        column_title = "CNV clone", row_title = "Top marker genes")
dev.off()

message("Pseudo-bulk DEG by sample_id x CNV clone using limma-trend...")
pb_group <- paste(meta$sample_id[tumor_idx], meta$CNV_clone[tumor_idx], sep = "|")
pb_levels <- unique(pb_group)
pb_counts <- matrix(0, nrow = nrow(counts), ncol = length(pb_levels), dimnames = list(rownames(counts), pb_levels))
for (pg_name in pb_levels) {
  cells <- tumor_cells[pb_group == pg_name]
  pb_counts[, pg_name] <- Matrix::rowSums(counts[, cells, drop = FALSE])
}
pb_meta <- data.frame(pb_group = pb_levels, stringsAsFactors = FALSE)
parts <- do.call(rbind, strsplit(pb_meta$pb_group, "\\|"))
pb_meta$sample_id <- parts[, 1]
pb_meta$CNV_clone <- parts[, 2]
lib_size <- colSums(pb_counts)
keep_gene <- rowSums(pb_counts >= 10) >= 2
pb_use <- pb_counts[keep_gene, , drop = FALSE]
logcpm <- log2(t(t(pb_use + 1) / (lib_size + 1) * 1e6))
design <- model.matrix(~ 0 + CNV_clone + sample_id, data = pb_meta)
colnames(design) <- make.names(colnames(design))
if (qr(design)$rank < ncol(design)) {
  design <- model.matrix(~ 0 + CNV_clone, data = pb_meta)
  colnames(design) <- make.names(colnames(design))
}
fit <- lmFit(logcpm, design)
fit <- eBayes(fit, trend = TRUE)
coef_cols <- grep("^CNV_clone", colnames(design), value = TRUE)
contrast_defs <- c()
ref_coef <- paste0("CNV_clone", clone_levels[1])
for (cl in clone_levels[-1]) {
  cc <- paste0("CNV_clone", cl)
  if (cc %in% coef_cols && ref_coef %in% coef_cols) contrast_defs <- c(contrast_defs, paste0(cc, "-", ref_coef))
}
if (length(contrast_defs) > 0) {
  names(contrast_defs) <- paste0(clone_levels[-1], "_vs_", clone_levels[1])
  cont <- makeContrasts(contrasts = contrast_defs, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont), trend = TRUE)
  pb_deg_all <- do.call(rbind, lapply(colnames(cont), function(coef) {
    tt <- topTable(fit2, coef = coef, number = Inf, sort.by = "P")
    tt$gene <- rownames(tt)
    tt$contrast <- coef
    tt
  }))
  write.csv(pb_deg_all, file.path(tables_dir, "pseudo_bulk_DEG_limma_clone_vs_Subclone_01.csv"), row.names = FALSE)
}
write.csv(pb_meta, file.path(tables_dir, "pseudo_bulk_sample_clone_metadata.csv"), row.names = FALSE)

message("CNV-expression dosage correlation at inferCNV subcluster level...")
group_expr_avg <- sapply(obs_groups, function(g) {
  cells <- meta$cell_integrated_oc[meta$infercnv_cell_group_name == g]
  Matrix::rowMeans(expr[gene_universe, cells, drop = FALSE])
})
common_dosage_genes <- intersect(rownames(cnv_gene_state), rownames(group_expr_avg))
dosage_res <- do.call(rbind, lapply(common_dosage_genes, function(g) {
  x <- as.numeric(cnv_gene_state[g, obs_groups])
  y <- as.numeric(group_expr_avg[g, obs_groups])
  if (sd(x) == 0 || sd(y) == 0) return(data.frame(gene = g, rho = NA_real_, p_value = NA_real_, cnv_sd = sd(x), expr_sd = sd(y)))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(gene = g, rho = as.numeric(ct$estimate), p_value = ct$p.value, cnv_sd = sd(x), expr_sd = sd(y))
}))
dosage_res$padj <- p.adjust(dosage_res$p_value, method = "BH")
dosage_res <- dosage_res[order(-dosage_res$rho), ]
write.csv(dosage_res, file.path(tables_dir, "CNV_expression_dosage_correlation_by_infercnv_subcluster.csv"), row.names = FALSE)
dosage_sig <- dosage_res[!is.na(dosage_res$padj) & dosage_res$padj < 0.05 & dosage_res$rho > 0.25, ]
write.csv(dosage_sig, file.path(tables_dir, "CNV_expression_dosage_genes_positive_rho_padj005.csv"), row.names = FALSE)

top_dosage_genes <- head(dosage_sig$gene, 50)
if (length(top_dosage_genes) >= 2) {
  dosage_heat <- group_expr_avg[top_dosage_genes, obs_groups, drop = FALSE]
  col_split <- cell_map$cnv_subclone[match(obs_groups, cell_map$cell_group_name)]
  pdf(file.path(figures_dir, "top_CNV_expression_dosage_gene_expression_by_subcluster.pdf"), width = 11, height = 10)
  Heatmap(t(scale(t(dosage_heat))), name = "expr z", col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
          cluster_rows = TRUE, cluster_columns = TRUE, column_split = col_split,
          row_names_gp = gpar(fontsize = 6), column_title = "inferCNV subclusters split by CNV clone")
  dev.off()
  png(file.path(figures_dir, "top_CNV_expression_dosage_gene_expression_by_subcluster.png"), width = 2200, height = 2000, res = 220)
  Heatmap(t(scale(t(dosage_heat))), name = "expr z", col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
          cluster_rows = TRUE, cluster_columns = TRUE, column_split = col_split,
          row_names_gp = gpar(fontsize = 6), column_title = "inferCNV subclusters split by CNV clone")
  dev.off()
}

message("CNV-transcriptome coupling correlations...")
model_df <- score_df
model_df$LAPTM5_expr <- if ("LAPTM5" %in% rownames(expr)) as.numeric(expr["LAPTM5", model_df$cell_integrated_oc]) else NA_real_
vars <- c("CNV_burden", "EMT", "Hypoxia", "KRAS_UP", "KRAS_DN", "LAPTM5_axis", "LAPTM5_expr", "Proliferation", "Immune_modulatory",
          "TF_SNAI2", "TF_ATF3", "TF_HIF1A", "TF_STAT3", "TF_MYC", "TF_FOXM1", "TF_RELA")
vars <- intersect(vars, colnames(model_df))
cor_mat <- cor(model_df[, vars], method = "spearman", use = "pairwise.complete.obs")
write.csv(cor_mat, file.path(tables_dir, "CNV_transcriptome_coupling_spearman_correlation_matrix.csv"))
cor_pairs <- do.call(rbind, combn(vars, 2, function(v) {
  ct <- suppressWarnings(cor.test(model_df[[v[1]]], model_df[[v[2]]], method = "spearman", exact = FALSE))
  data.frame(var1 = v[1], var2 = v[2], rho = as.numeric(ct$estimate), p_value = ct$p.value)
}, simplify = FALSE))
cor_pairs$padj <- p.adjust(cor_pairs$p_value, "BH")
cor_pairs <- cor_pairs[order(cor_pairs$padj, -abs(cor_pairs$rho)), ]
write.csv(cor_pairs, file.path(tables_dir, "CNV_transcriptome_coupling_spearman_pairs.csv"), row.names = FALSE)
pdf(file.path(figures_dir, "CNV_transcriptome_coupling_correlation_heatmap.pdf"), width = 8, height = 8)
corrplot(cor_mat, method = "color", type = "upper", tl.col = "black", tl.cex = 0.7, mar = c(0,0,2,0), title = "CNV-transcriptome coupling")
dev.off()
png(file.path(figures_dir, "CNV_transcriptome_coupling_correlation_heatmap.png"), width = 1600, height = 1600, res = 220)
corrplot(cor_mat, method = "color", type = "upper", tl.col = "black", tl.cex = 0.7, mar = c(0,0,2,0), title = "CNV-transcriptome coupling")
dev.off()

clone_summary <- merge(clone_burden[, c("CNV_clone", "mean", "median")], 
                       reshape(clone_function_summary[, c("CNV_clone", "feature", "mean_score")], idvar = "CNV_clone", timevar = "feature", direction = "wide"),
                       by = "CNV_clone", all = TRUE)
colnames(clone_summary)[colnames(clone_summary) == "mean"] <- "CNV_burden_mean"
colnames(clone_summary)[colnames(clone_summary) == "median"] <- "CNV_burden_median"
write.csv(clone_summary, file.path(tables_dir, "CNV_clone_functional_summary.csv"), row.names = FALSE)

message("Done CNV-expression joint analysis.")
