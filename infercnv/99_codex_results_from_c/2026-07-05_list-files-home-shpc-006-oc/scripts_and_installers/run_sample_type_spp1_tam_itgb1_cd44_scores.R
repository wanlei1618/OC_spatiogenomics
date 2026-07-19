suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

root_dir <- "D:/OC_spatiogenomics/infercnv"
out_dir <- file.path(root_dir, "sample_type_SPP1_TAM_ITGB1_CD44_scores")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
script_dir <- file.path(out_dir, "scripts")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(root_dir, "integrated_oc_plan_analysis/tables/integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv")
object_path <- file.path(root_dir, "integrated_oc.RData")

spp1_tam_genes <- c("SPP1", "CD68", "CD163", "CD14", "LST1", "TYROBP",
                    "C1QA", "C1QB", "C1QC", "APOE", "MRC1", "MSR1",
                    "FCGR3A", "ITGAM", "CSF1R")
itgb1_cd44_genes <- c("ITGB1", "CD44")

read_seurat_object <- function(path) {
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(obj)) return(obj)
  env <- new.env()
  load(path, envir = env)
  env[[ls(env)[1]]]
}

score_signature_zmean <- function(expr_mat, genes, cells) {
  genes_present <- intersect(genes, rownames(expr_mat))
  if (length(genes_present) == 0) {
    return(list(score = rep(NA_real_, length(cells)), genes_present = character()))
  }
  sub <- as.matrix(expr_mat[genes_present, cells, drop = FALSE])
  z <- t(scale(t(sub)))
  z[!is.finite(z)] <- 0
  list(score = colMeans(z, na.rm = TRUE), genes_present = genes_present)
}

extract_gene_expr <- function(expr_mat, genes, cells) {
  out <- data.frame(cell_integrated_oc = cells, stringsAsFactors = FALSE)
  for (g in genes) {
    out[[paste0(g, "_expr")]] <- if (g %in% rownames(expr_mat)) as.numeric(expr_mat[g, cells, drop = TRUE]) else NA_real_
  }
  out
}

summary_by_sample_type <- function(d, score_col, analysis_label) {
  d <- d[is.finite(d[[score_col]]) & !is.na(d$sample_type) & d$sample_type != "", , drop = FALSE]
  if (nrow(d) == 0) return(data.frame())
  rows <- lapply(sort(unique(d$sample_type)), function(st) {
    x <- d[d$sample_type == st, score_col]
    data.frame(analysis = analysis_label, sample_type = st, n = length(x),
               mean = mean(x), median = median(x), sd = sd(x),
               q25 = as.numeric(quantile(x, 0.25)),
               q75 = as.numeric(quantile(x, 0.75)))
  })
  rbindlist(rows)
}

kruskal_test <- function(d, score_col, analysis_label) {
  d <- d[is.finite(d[[score_col]]) & !is.na(d$sample_type) & d$sample_type != "", , drop = FALSE]
  if (length(unique(d$sample_type)) < 2) {
    return(data.frame(analysis = analysis_label, statistic = NA_real_, p_value = NA_real_, n_groups = length(unique(d$sample_type)), n_cells = nrow(d)))
  }
  kt <- kruskal.test(as.formula(paste(score_col, "~ sample_type")), data = d)
  data.frame(analysis = analysis_label, statistic = unname(kt$statistic), p_value = kt$p.value,
             n_groups = length(unique(d$sample_type)), n_cells = nrow(d))
}

pairwise_wilcox <- function(d, score_col, analysis_label) {
  d <- d[is.finite(d[[score_col]]) & !is.na(d$sample_type) & d$sample_type != "", , drop = FALSE]
  groups <- sort(unique(d$sample_type))
  if (length(groups) < 2) return(data.frame())
  rows <- list()
  idx <- 1
  for (i in seq_len(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      x <- d[d$sample_type == g1, score_col]
      y <- d[d$sample_type == g2, score_col]
      wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
      rows[[idx]] <- data.frame(analysis = analysis_label, group1 = g1, group2 = g2,
                                n1 = length(x), n2 = length(y),
                                median1 = median(x), median2 = median(y),
                                delta_median_group2_minus_group1 = median(y) - median(x),
                                p_value = wt$p.value)
      idx <- idx + 1
    }
  }
  out <- rbindlist(rows)
  out$p_adj_BH <- p.adjust(out$p_value, method = "BH")
  out
}

plot_score <- function(d, score_col, title, subtitle, file_prefix) {
  d <- d[is.finite(d[[score_col]]) & !is.na(d$sample_type) & d$sample_type != "", , drop = FALSE]
  d$sample_type <- factor(d$sample_type, levels = names(sort(table(d$sample_type), decreasing = TRUE)))
  p <- ggplot(d, aes(x = sample_type, y = .data[[score_col]], fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.82) +
    geom_jitter(width = 0.18, size = 0.35, alpha = 0.22, color = "grey20") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2.4, fill = "white") +
    labs(title = title, subtitle = subtitle, x = "sample_type", y = score_col) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, paste0(file_prefix, ".png")), p, width = 7.2, height = 5.2, dpi = 220)
  ggsave(file.path(fig_dir, paste0(file_prefix, ".pdf")), p, width = 7.2, height = 5.2)
}

obj <- read_seurat_object(object_path)
meta <- fread(metadata_path, data.table = FALSE)
expr <- GetAssayData(obj, assay = "RNA", slot = "data")

common_cells <- intersect(colnames(expr), meta$cell_integrated_oc)
meta <- meta[match(common_cells, meta$cell_integrated_oc), , drop = FALSE]
rownames(meta) <- meta$cell_integrated_oc

myeloid_flag <- grepl("Myeloid|Macro", meta$interaction_group) |
  meta$cell_type %in% c("Macrophages", "Monocytes", "DC")
myeloid_cells <- meta$cell_integrated_oc[myeloid_flag]

tumor_flag <- !is.na(meta$cnv_subclone) & meta$cnv_subclone != ""
tumor_cells <- meta$cell_integrated_oc[tumor_flag]

spp1_score <- score_signature_zmean(expr, spp1_tam_genes, myeloid_cells)
spp1_df <- meta[myeloid_cells, c("cell_integrated_oc", "sample_type", "cell_type", "myeloid_subtype", "interaction_group"), drop = FALSE]
spp1_df$SPP1_TAM_score <- spp1_score$score
spp1_expr <- extract_gene_expr(expr, c("SPP1", "CD68", "CD163", "APOE", "C1QA", "C1QB", "C1QC"), myeloid_cells)
spp1_df <- merge(spp1_df, spp1_expr, by = "cell_integrated_oc", all.x = TRUE, sort = FALSE)

itgb_score <- score_signature_zmean(expr, itgb1_cd44_genes, tumor_cells)
itgb_df <- meta[tumor_cells, c("cell_integrated_oc", "sample_type", "cell_type", "cnv_subclone", "interaction_group"), drop = FALSE]
itgb_df$ITGB1_CD44_tumor_score <- itgb_score$score
itgb_expr <- extract_gene_expr(expr, c("ITGB1", "CD44", "EPCAM", "KRT8", "PAX8", "MUC16"), tumor_cells)
itgb_df <- merge(itgb_df, itgb_expr, by = "cell_integrated_oc", all.x = TRUE, sort = FALSE)

fwrite(spp1_df, file.path(tab_dir, "cell_level_SPP1_TAM_score_myeloid_by_sample_type.csv"))
fwrite(itgb_df, file.path(tab_dir, "cell_level_ITGB1_CD44_tumor_score_CNV_cells_by_sample_type.csv"))

gene_coverage <- rbind(
  data.frame(score = "SPP1_TAM_score", requested_gene = spp1_tam_genes, present = spp1_tam_genes %in% rownames(expr)),
  data.frame(score = "ITGB1_CD44_tumor_score", requested_gene = itgb1_cd44_genes, present = itgb1_cd44_genes %in% rownames(expr))
)
fwrite(gene_coverage, file.path(tab_dir, "score_gene_coverage.csv"))

summary_tab <- rbind(
  summary_by_sample_type(spp1_df, "SPP1_TAM_score", "myeloid_SPP1_TAM_score"),
  summary_by_sample_type(itgb_df, "ITGB1_CD44_tumor_score", "CNV_tumor_ITGB1_CD44_score")
)
fwrite(summary_tab, file.path(tab_dir, "sample_type_score_summary.csv"))

kruskal_tab <- rbind(
  kruskal_test(spp1_df, "SPP1_TAM_score", "myeloid_SPP1_TAM_score"),
  kruskal_test(itgb_df, "ITGB1_CD44_tumor_score", "CNV_tumor_ITGB1_CD44_score")
)
kruskal_tab$p_adj_BH <- p.adjust(kruskal_tab$p_value, method = "BH")
fwrite(kruskal_tab, file.path(tab_dir, "sample_type_score_kruskal_wallis_tests.csv"))

wilcox_tab <- rbind(
  pairwise_wilcox(spp1_df, "SPP1_TAM_score", "myeloid_SPP1_TAM_score"),
  pairwise_wilcox(itgb_df, "ITGB1_CD44_tumor_score", "CNV_tumor_ITGB1_CD44_score")
)
fwrite(wilcox_tab, file.path(tab_dir, "sample_type_score_pairwise_wilcoxon_tests.csv"))

plot_score(spp1_df, "SPP1_TAM_score",
           "SPP1_TAM_score across sample_type",
           paste0("Myeloid/TAM cells; genes present: ", paste(spp1_score$genes_present, collapse = ", ")),
           "SPP1_TAM_score_by_sample_type_myeloid")
plot_score(itgb_df, "ITGB1_CD44_tumor_score",
           "ITGB1_CD44_tumor_score across sample_type",
           paste0("CNV tumor cells; genes present: ", paste(itgb_score$genes_present, collapse = ", ")),
           "ITGB1_CD44_tumor_score_by_sample_type_CNV_tumor")

report_path <- file.path(out_dir, "sample_type_SPP1_TAM_ITGB1_CD44_score_report.md")
sink(report_path)
cat("# sample_type association for SPP1_TAM and ITGB1_CD44 tumor scores\n\n")
cat("Input object: `", object_path, "`\n\n", sep = "")
cat("Metadata: `", metadata_path, "`\n\n", sep = "")
cat("## Score definitions\n\n")
cat("- `SPP1_TAM_score`: mean z-scored RNA log-normalized expression of present TAM genes in myeloid/macrophage cells.\n")
cat("- `ITGB1_CD44_tumor_score`: mean z-scored RNA log-normalized expression of ITGB1 and CD44 in CNV-subclone tumor cells.\n\n")
cat("## Gene coverage\n\n")
print(gene_coverage)
cat("\n## sample_type summary\n\n")
print(summary_tab)
cat("\n## Kruskal-Wallis tests\n\n")
print(kruskal_tab)
cat("\n## Pairwise Wilcoxon tests\n\n")
print(wilcox_tab)
cat("\n## Output files\n\n")
cat("- `tables/cell_level_SPP1_TAM_score_myeloid_by_sample_type.csv`\n")
cat("- `tables/cell_level_ITGB1_CD44_tumor_score_CNV_cells_by_sample_type.csv`\n")
cat("- `tables/sample_type_score_summary.csv`\n")
cat("- `tables/sample_type_score_kruskal_wallis_tests.csv`\n")
cat("- `tables/sample_type_score_pairwise_wilcoxon_tests.csv`\n")
cat("- `figures/SPP1_TAM_score_by_sample_type_myeloid.png`\n")
cat("- `figures/ITGB1_CD44_tumor_score_by_sample_type_CNV_tumor.png`\n")
sink()

file.copy(normalizePath("run_sample_type_spp1_tam_itgb1_cd44_scores.R", winslash = "/", mustWork = FALSE),
          file.path(script_dir, "run_sample_type_spp1_tam_itgb1_cd44_scores.R"), overwrite = TRUE)
