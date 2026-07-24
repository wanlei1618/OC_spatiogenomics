options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
out_dir <- file.path(data_root, "research_validation_independent_cnv")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(out_dir, "copykat_defined_bias_by_patient.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

extract_assay5 <- function(path) {
  x <- readRDS(path)
  assay <- x@assays[[x@active.assay]]
  m <- attr(assay, "layers")$counts
  dimnames(m) <- list(
    rownames(attr(assay, "features"))[seq_len(nrow(m))],
    rownames(attr(assay, "cells"))[seq_len(ncol(m))]
  )
  list(counts = as(m, "dgCMatrix"), metadata = x@meta.data)
}
full <- extract_assay5(file.path(
  data_root, "GSE154600", "objects", "GSE154600_preannotation.rds"
))
stability <- fread(file.path(
  out_dir, "GSE154600_copykat_stability_by_cell_v2.csv.gz"
))
ids <- intersect(stability$cell_id, colnames(full$counts))
counts <- full$counts[, ids, drop = FALSE]
metadata <- as.data.table(full$metadata, keep.rownames = "cell_id")
cell <- merge(
  stability,
  metadata[, .(cell_id, nCount_RNA, nFeature_RNA, percent.mt)],
  by = "cell_id", all.x = TRUE
)

genes <- c("EPCAM", "KRT8", "KRT18", "KRT19", "CD44", "ITGB1")
raw <- counts[intersect(genes, rownames(counts)), , drop = FALSE]
lib <- pmax(Matrix::colSums(counts), 1)
norm <- raw
if (length(norm@x))
  norm@x <- log1p(norm@x * rep(1e4 / lib, diff(norm@p)))
expr <- function(g) {
  if (g %in% rownames(norm)) as.numeric(norm[g, ]) else
    rep(NA_real_, ncol(norm))
}
expr_dt <- data.table(cell_id = colnames(counts))
for (g in genes) expr_dt[[g]] <- expr(g)
expr_dt[, epithelial_marker_score :=
          rowMeans(.SD, na.rm = TRUE), .SDcols = c("EPCAM", "KRT8", "KRT18", "KRT19")]
cell <- merge(cell, expr_dt, by = "cell_id", all.x = TRUE)

metrics <- c(
  "nCount_RNA", "nFeature_RNA", "percent.mt", "EPCAM", "KRT8",
  "KRT18", "KRT19", "epithelial_marker_score", "CD44", "ITGB1"
)
summarize_metric <- function(v, prefix) {
  ans <- c(
    median = median(v, na.rm = TRUE),
    q25 = as.numeric(quantile(v, .25, na.rm = TRUE, names = FALSE)),
    q75 = as.numeric(quantile(v, .75, na.rm = TRUE, names = FALSE))
  )
  names(ans) <- paste(prefix, names(ans), sep = "_")
  as.list(ans)
}
by_group <- cell[, {
  ans <- list(n_cells = .N)
  for (metric in metrics)
    ans <- c(ans, summarize_metric(get(metric), metric))
  ans
}, by = .(dataset_id, patient_id, stability_class)]
fwrite(by_group, out, na = "NA")

bias_one <- function(d, patient) {
  a <- d[stability_class == "STABLE_ANEUPLOID"]
  n <- d[stability_class == "MOSTLY_NOT_DEFINED"]
  enough <- nrow(a) >= 20L && nrow(n) >= 20L
  count_ratio <- if (enough)
    median(a$nCount_RNA, na.rm = TRUE) / pmax(median(n$nCount_RNA, na.rm = TRUE), 1) else NA_real_
  feature_ratio <- if (enough)
    median(a$nFeature_RNA, na.rm = TRUE) / pmax(median(n$nFeature_RNA, na.rm = TRUE), 1) else NA_real_
  marker_shift <- if (enough)
    median(a$epithelial_marker_score, na.rm = TRUE) -
      median(n$epithelial_marker_score, na.rm = TRUE) else NA_real_
  status <- if (!enough) "INSUFFICIENT_DATA" else
    if (count_ratio > 2 || feature_ratio > 1.5)
      "HIGHER_DEPTH_IN_DEFINED_CALLS" else
        if (count_ratio < .5 || feature_ratio < (1 / 1.5))
          "LOWER_DEPTH_IN_NOT_DEFINED" else
            if (abs(marker_shift) > .5) "MARKER_COMPOSITION_SHIFT" else
              "NO_MAJOR_TECHNICAL_SHIFT"
  data.table(
    dataset_id = "GSE154600",
    patient_id = patient,
    n_stable_aneuploid = nrow(a),
    n_mostly_not_defined = nrow(n),
    stable_to_notdefined_nCount_median_ratio = count_ratio,
    stable_to_notdefined_nFeature_median_ratio = feature_ratio,
    epithelial_marker_score_median_difference = marker_shift,
    copykat_defined_bias_status = status,
    potential_detection_depth_selection_bias =
      isTRUE(count_ratio > 2 || feature_ratio > 1.5),
    inference_scope =
      "descriptive_technical_audit_not_patient_level_biological_evidence"
  )
}
summary <- rbindlist(lapply(
  unique(cell$patient_id),
  function(pt) bias_one(cell[patient_id == pt], pt)
))
overall_status <- if (any(summary$potential_detection_depth_selection_bias))
  "HIGHER_DEPTH_IN_DEFINED_CALLS" else
    if (all(summary$copykat_defined_bias_status == "INSUFFICIENT_DATA"))
      "INSUFFICIENT_DATA" else "NO_MAJOR_TECHNICAL_SHIFT"
summary <- rbind(
  summary,
  data.table(
    dataset_id = "GSE154600", patient_id = "ALL",
    n_stable_aneuploid = sum(summary$n_stable_aneuploid),
    n_mostly_not_defined = sum(summary$n_mostly_not_defined),
    stable_to_notdefined_nCount_median_ratio = NA_real_,
    stable_to_notdefined_nFeature_median_ratio = NA_real_,
    epithelial_marker_score_median_difference = NA_real_,
    copykat_defined_bias_status = overall_status,
    potential_detection_depth_selection_bias =
      any(summary$potential_detection_depth_selection_bias),
    inference_scope =
      "patient_stratified_descriptive_technical_audit"
  )
)
fwrite(summary, file.path(out_dir, "copykat_defined_bias_summary.csv"), na = "NA")

plot_dt <- cell[
  stability_class %in% c("STABLE_ANEUPLOID", "MOSTLY_NOT_DEFINED"),
  .(patient_id, stability_class, nCount_RNA, nFeature_RNA)
]
plot_dt <- melt(
  plot_dt,
  id.vars = c("patient_id", "stability_class"),
  measure.vars = c("nCount_RNA", "nFeature_RNA"),
  variable.name = "metric", value.name = "value"
)
p <- ggplot(plot_dt, aes(stability_class, value, fill = stability_class)) +
  geom_boxplot(outlier.shape = NA, width = .7) +
  scale_y_log10() +
  facet_grid(metric ~ patient_id, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "CopyKAT defined versus mostly-not-defined technical depth",
    x = NULL, y = "value (log10 scale)"
  )
ggsave(file.path(out_dir, "copykat_defined_vs_notdefined_qc.png"),
       p, width = 12, height = 6, dpi = 180)
message("CopyKAT defined-cell technical selection audit complete")
