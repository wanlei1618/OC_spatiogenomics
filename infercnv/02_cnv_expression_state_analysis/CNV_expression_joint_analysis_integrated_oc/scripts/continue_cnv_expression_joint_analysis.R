suppressPackageStartupMessages({
  library(Matrix)
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
workspace_out <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"

message("Loading integrated_oc...")
obj <- readRDS("D:/OC_spatiogenomics/infercnv/integrated_oc.RData")
expr <- obj@assays$RNA@data
meta <- read.csv(file.path(tables_dir, "integrated_oc_metadata_with_CNV_clone_and_subtypes.csv"), stringsAsFactors = FALSE)
cnv_gene_state <- readRDS(file.path(tables_dir, "cnv_gene_state_matrix_gene_by_infercnv_subcluster.rds"))
cell_map <- read.csv(file.path(workspace_out, "infercnv_cell_to_subclone_k5.csv"), stringsAsFactors = FALSE)

message("Computing subcluster-level expression averages for dosage correlation...")
obs_groups <- colnames(cnv_gene_state)
genes <- intersect(rownames(cnv_gene_state), rownames(expr))
group_expr_avg <- matrix(NA_real_, nrow = length(genes), ncol = length(obs_groups), dimnames = list(genes, obs_groups))
for (g in obs_groups) {
  cells <- intersect(meta$cell_integrated_oc[meta$infercnv_cell_group_name == g], colnames(expr))
  if (length(cells) > 0) {
    group_expr_avg[, g] <- Matrix::rowMeans(expr[genes, cells, drop = FALSE])
  }
}
valid_groups <- colnames(group_expr_avg)[colSums(!is.na(group_expr_avg)) > 0]
group_expr_avg <- group_expr_avg[, valid_groups, drop = FALSE]
cnv_use <- cnv_gene_state[genes, valid_groups, drop = FALSE]

message("Running gene-wise CNV-expression dosage correlations...")
dosage_res <- do.call(rbind, lapply(genes, function(gene) {
  x <- as.numeric(cnv_use[gene, ])
  y <- as.numeric(group_expr_avg[gene, ])
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 8 || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(data.frame(gene = gene, rho = NA_real_, p_value = NA_real_, n_groups = sum(keep), cnv_sd = sd(x[keep]), expr_sd = sd(y[keep])))
  }
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.frame(gene = gene, rho = as.numeric(ct$estimate), p_value = ct$p.value, n_groups = sum(keep), cnv_sd = sd(x[keep]), expr_sd = sd(y[keep]))
}))
dosage_res$padj <- p.adjust(dosage_res$p_value, method = "BH")
dosage_res <- dosage_res[order(-dosage_res$rho), ]
write.csv(dosage_res, file.path(tables_dir, "CNV_expression_dosage_correlation_by_infercnv_subcluster.csv"), row.names = FALSE)
dosage_sig <- dosage_res[!is.na(dosage_res$padj) & dosage_res$padj < 0.05 & dosage_res$rho > 0.25, ]
write.csv(dosage_sig, file.path(tables_dir, "CNV_expression_dosage_genes_positive_rho_padj005.csv"), row.names = FALSE)

top_dosage_genes <- head(dosage_sig$gene, 50)
if (length(top_dosage_genes) >= 2) {
  dosage_heat <- group_expr_avg[top_dosage_genes, , drop = FALSE]
  col_split <- cell_map$cnv_subclone[match(colnames(dosage_heat), cell_map$cell_group_name)]
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

message("CNV-transcriptome coupling model...")
score_df <- read.csv(file.path(tables_dir, "tumor_cell_CNV_burden_signature_TF_scores.csv"), stringsAsFactors = FALSE, check.names = FALSE)
if ("LAPTM5" %in% rownames(expr)) {
  score_df$LAPTM5_expr <- as.numeric(expr["LAPTM5", score_df$cell_integrated_oc])
}
vars <- c("CNV_burden", "EMT", "Hypoxia", "KRAS_UP", "KRAS_DN", "LAPTM5_axis", "LAPTM5_expr", "Proliferation", "Immune_modulatory",
          "TF_SNAI2", "TF_ATF3", "TF_HIF1A", "TF_STAT3", "TF_MYC", "TF_FOXM1", "TF_RELA")
vars <- intersect(vars, colnames(score_df))
cor_mat <- cor(score_df[, vars], method = "spearman", use = "pairwise.complete.obs")
write.csv(cor_mat, file.path(tables_dir, "CNV_transcriptome_coupling_spearman_correlation_matrix.csv"))
cor_pairs <- do.call(rbind, combn(vars, 2, function(v) {
  ct <- suppressWarnings(cor.test(score_df[[v[1]]], score_df[[v[2]]], method = "spearman", exact = FALSE))
  data.frame(var1 = v[1], var2 = v[2], rho = as.numeric(ct$estimate), p_value = ct$p.value)
}, simplify = FALSE))
cor_pairs$padj <- p.adjust(cor_pairs$p_value, "BH")
cor_pairs <- cor_pairs[order(cor_pairs$padj, -abs(cor_pairs$rho)), ]
write.csv(cor_pairs, file.path(tables_dir, "CNV_transcriptome_coupling_spearman_pairs.csv"), row.names = FALSE)

pdf(file.path(figures_dir, "CNV_transcriptome_coupling_correlation_heatmap.pdf"), width = 8, height = 8)
corrplot(cor_mat, method = "color", type = "upper", tl.col = "black", tl.cex = 0.7, mar = c(0, 0, 2, 0), title = "CNV-transcriptome coupling")
dev.off()
png(file.path(figures_dir, "CNV_transcriptome_coupling_correlation_heatmap.png"), width = 1600, height = 1600, res = 220)
corrplot(cor_mat, method = "color", type = "upper", tl.col = "black", tl.cex = 0.7, mar = c(0, 0, 2, 0), title = "CNV-transcriptome coupling")
dev.off()

clone_burden <- read.csv(file.path(tables_dir, "CNV_clone_burden_summary.csv"), stringsAsFactors = FALSE)
clone_func <- read.csv(file.path(tables_dir, "CNV_clone_functional_and_TF_activity_summary.csv"), stringsAsFactors = FALSE)
wide <- reshape(clone_func[, c("CNV_clone", "feature", "mean_score")], idvar = "CNV_clone", timevar = "feature", direction = "wide")
clone_summary <- merge(clone_burden[, c("CNV_clone", "mean", "median")], wide, by = "CNV_clone", all = TRUE)
colnames(clone_summary)[colnames(clone_summary) == "mean"] <- "CNV_burden_mean"
colnames(clone_summary)[colnames(clone_summary) == "median"] <- "CNV_burden_median"
write.csv(clone_summary, file.path(tables_dir, "CNV_clone_functional_summary.csv"), row.names = FALSE)

message("Done continuation.")
