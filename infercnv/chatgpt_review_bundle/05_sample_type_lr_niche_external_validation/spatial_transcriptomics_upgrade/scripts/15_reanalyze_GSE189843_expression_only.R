options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

local_root <- "D:/OC_spatiogenomics/spatial_data/spatial_analysis_correction_v2"
extract_dir <- "D:/OC_spatiogenomics/spatial_data/raw/GSE189843/extracted"
audit <- fread(file.path(local_root, "spatial_sample_audit.csv"))
sample_audit <- audit[dataset_id == "GSE189843", .(
  dataset_id, sample_id, response_group, included_by_author,
  available_in_matrix, final_include, exclusion_reason,
  analysis_level = "EXPRESSION_ONLY",
  statistical_unit = "sample/patient"
)]
fwrite(sample_audit, file.path(local_root, "GSE189843_sample_audit.csv"),
       na = "NA")
signatures <- fread(file.path(local_root, "updated_spatial_signatures.csv"))
get_signature <- function(name) {
  signatures[signature_name == name, unique(gene)]
}
needed <- unique(c(
  get_signature("macrophage_identity"),
  get_signature("SPP1_macrophage_program"),
  get_signature("C1QC_macrophage_control"),
  get_signature("epithelial_identity"),
  get_signature("CNV_supported_malignant_epithelial_signature"),
  get_signature("Subclone02_04_KRAS_hypoxia"), "ITGB1", "CD44"
))
zscore <- function(v) {
  s <- sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(v)))
  (v - mean(v, na.rm = TRUE)) / s
}
score_signature <- function(norm, genes) {
  genes <- intersect(genes, rownames(norm))
  if (length(genes) < 2L) return(rep(NA_real_, ncol(norm)))
  z <- vapply(genes, function(g) zscore(as.numeric(norm[g, ])),
              numeric(ncol(norm)))
  rowMeans(z, na.rm = TRUE)
}
high_top <- function(v, fraction = .25) {
  v >= quantile(v, 1 - fraction, na.rm = TRUE, names = FALSE)
}
read_one <- function(sid) {
  matrix_file <- list.files(
    extract_dir, pattern = paste0("^", sid, ".*matrix.*\\.mtx$"),
    full.names = TRUE
  )[[1L]]
  feature_file <- list.files(
    extract_dir, pattern = paste0("^", sid, ".*features.*\\.tsv$"),
    full.names = TRUE
  )[[1L]]
  barcode_file <- list.files(
    extract_dir, pattern = paste0("^", sid, ".*barcodes.*\\.tsv$"),
    full.names = TRUE
  )[[1L]]
  features <- fread(feature_file, header = FALSE)
  barcodes <- fread(barcode_file, header = FALSE)[[1L]]
  gene_names <- make.unique(as.character(features[[2L]]))
  m <- as(readMM(matrix_file), "CsparseMatrix")
  rownames(m) <- gene_names
  colnames(m) <- barcodes
  keep <- intersect(needed, rownames(m))
  lib <- Matrix::colSums(m)
  norm <- t(t(m[keep, , drop = FALSE]) / pmax(lib, 1)) * 1e4
  norm@x <- log1p(norm@x)
  gene_expr <- function(g) {
    if (g %in% rownames(norm)) as.numeric(norm[g, ]) else
      rep(NA_real_, ncol(norm))
  }
  macrophage <- score_signature(norm, get_signature("macrophage_identity"))
  spp1 <- score_signature(norm, get_signature("SPP1_macrophage_program"))
  c1qc <- score_signature(norm, get_signature("C1QC_macrophage_control"))
  epithelial <- score_signature(norm, get_signature("epithelial_identity"))
  malignant <- score_signature(
    norm, get_signature("CNV_supported_malignant_epithelial_signature")
  )
  subclone <- score_signature(
    norm, get_signature("Subclone02_04_KRAS_hypoxia")
  )
  itgb1 <- gene_expr("ITGB1")
  cd44 <- gene_expr("CD44")
  epi_base <- epithelial >= median(epithelial, na.rm = TRUE)
  malignant_base <- malignant >= median(malignant, na.rm = TRUE)
  spp1_high <- macrophage >= median(macrophage, na.rm = TRUE) &
    high_top(spp1)
  itgb1_high <- epi_base & malignant_base & high_top(itgb1)
  cd44_high <- epi_base & malignant_base & high_top(cd44)
  cor_itgb1 <- suppressWarnings(cor(spp1, itgb1, method = "spearman",
                                    use = "pairwise.complete.obs"))
  cor_cd44 <- suppressWarnings(cor(spp1, cd44, method = "spearman",
                                   use = "pairwise.complete.obs"))
  data.table(
    sample_id = sid, n_spots = ncol(m),
    median_SPP1_program = median(spp1, na.rm = TRUE),
    median_C1QC_program = median(c1qc, na.rm = TRUE),
    median_ITGB1 = median(itgb1, na.rm = TRUE),
    median_CD44 = median(cd44, na.rm = TRUE),
    median_malignant_epithelial_score = median(malignant, na.rm = TRUE),
    median_Subclone02_04_KRAS_hypoxia_score =
      median(subclone, na.rm = TRUE),
    SPP1_ITGB1_sample_internal_spearman = cor_itgb1,
    SPP1_CD44_sample_internal_spearman = cor_cd44,
    high_SPP1_spot_fraction = mean(spp1_high),
    high_ITGB1_receiver_spot_fraction = mean(itgb1_high),
    high_CD44_receiver_spot_fraction = mean(cd44_high),
    analysis_level = "EXPRESSION_ONLY",
    statistical_unit = "sample/patient"
  )
}

rows <- lapply(sample_audit[final_include == TRUE, sample_id], function(sid) {
  message("Expression-only scoring: ", sid)
  read_one(sid)
})
summary <- rbindlist(rows, fill = TRUE)
summary <- merge(
  sample_audit[, .(sample_id, response_group)], summary,
  by = "sample_id", all.y = TRUE
)
fwrite(summary,
       file.path(local_root, "GSE189843_expression_summary_by_sample.csv"),
       na = "NA")

metrics <- setdiff(names(summary), c(
  "sample_id", "response_group", "n_spots", "analysis_level",
  "statistical_unit"
))
set.seed(20260730)
comparison <- rbindlist(lapply(metrics, function(metric) {
  excellent <- summary[response_group == "Excellent", get(metric)]
  poor <- summary[response_group == "Poor", get(metric)]
  if (!length(excellent) || !length(poor)) return(NULL)
  wt <- suppressWarnings(wilcox.test(excellent, poor, exact = TRUE))
  pair <- outer(poor, excellent, "-")
  cliffs <- (sum(pair > 0) - sum(pair < 0)) / length(pair)
  boot <- replicate(
    2000L,
    median(sample(poor, replace = TRUE)) -
      median(sample(excellent, replace = TRUE))
  )
  data.table(
    metric,
    n_excellent = length(excellent), n_poor = length(poor),
    excellent_median = median(excellent),
    poor_median = median(poor),
    median_difference_poor_minus_excellent =
      median(poor) - median(excellent),
    exact_wilcoxon_p = wt$p.value,
    cliffs_delta_poor_vs_excellent = cliffs,
    bootstrap_ci_low = quantile(boot, .025, names = FALSE),
    bootstrap_ci_high = quantile(boot, .975, names = FALSE),
    analysis_level = "EXPRESSION_ONLY",
    statistical_unit = "sample/patient"
  )
}))
comparison[, BH_FDR := p.adjust(exact_wilcoxon_p, "BH")]
fwrite(comparison,
       file.path(local_root, "GSE189843_response_group_comparison.csv"),
       na = "NA")

plot_metrics <- c(
  "median_SPP1_program", "median_ITGB1", "median_CD44",
  "median_malignant_epithelial_score",
  "SPP1_ITGB1_sample_internal_spearman",
  "SPP1_CD44_sample_internal_spearman"
)
plot_dt <- melt(
  summary,
  id.vars = c("sample_id", "response_group"),
  measure.vars = plot_metrics,
  variable.name = "metric", value.name = "value"
)
p <- ggplot(plot_dt, aes(response_group, value, color = response_group)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +
  geom_point(position = position_jitter(width = .08), size = 1.8) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "GSE189843 expression-only sample-level comparison",
    x = NULL, y = "sample-level value"
  )
ggsave(file.path(local_root, "GSE189843_expression_response_dotplot.png"),
       p, width = 11, height = 6, dpi = 180)
message("GSE189843 expression-only patient-level analysis complete")
