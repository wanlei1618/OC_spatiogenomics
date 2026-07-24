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
out_dir <- file.path(data_root, "research_spatial_transition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cnv <- fread(file.path(out_dir, "GSE154600_calibrated_cnv_by_cell.csv.gz"))
counts <- readRDS(file.path(
  data_root, "research_validation_independent_cnv",
  "GSE154600_complete_final_epithelial_counts.rds"
))
assignment <- fread(file.path(
  data_root, "diagnostics_v2_marker_ready_cleaned", "GSE154600",
  "cleaned_cell_assignments.csv.gz"
))
object <- readRDS(file.path(
  data_root, "GSE154600", "objects", "GSE154600_preannotation.rds"
))
meta <- as.data.table(object@meta.data, keep.rownames = "cell_id")
rm(object)

cells <- intersect(cnv$cell_id, colnames(counts))
stopifnot(length(cells) == nrow(cnv))
counts <- counts[, cells, drop = FALSE]
lib <- Matrix::colSums(counts)
genes <- c("EPCAM", "KRT8", "KRT18", "KRT19")
present <- intersect(genes, rownames(counts))
log_cpm <- log1p(t(t(counts[present, , drop = FALSE]) / pmax(lib, 1)) * 1e4)
epithelial_score <- Matrix::colMeans(log_cpm)
qc <- data.table(
  cell_id = cells, library_size = as.numeric(lib),
  epithelial_score = as.numeric(epithelial_score)
)
for (g in genes) {
  qc[[paste0(g, "_log1p_cpm")]] <- if (g %in% present) {
    as.numeric(log_cpm[g, ])
  } else {
    NA_real_
  }
}

x <- merge(cnv, assignment[, .(cell_id, final_cluster)],
           by = "cell_id", all.x = TRUE)
qc_cols <- intersect(c("cell_id", "nCount_RNA", "nFeature_RNA", "percent.mt"),
                     names(meta))
x <- merge(x, meta[, ..qc_cols], by = "cell_id", all.x = TRUE)
x <- merge(x, qc, by = "cell_id", all.x = TRUE)
x[, evidence_group := fcase(
  integrated_calibrated_cnv_evidence == "CALIBRATED_DUAL_METHOD_SUPPORT",
  "dual_support",
  integrated_calibrated_cnv_evidence == "COPYKAT_ONLY_SUPPORT",
  "copykat_only",
  integrated_calibrated_cnv_evidence == "INFERCNV_ONLY_SUPPORT",
  "infercnv_only",
  default = "neither_or_not_evaluable"
)]

cluster <- x[, .(
  n_cells = .N,
  n_dual_support = sum(evidence_group == "dual_support"),
  n_copykat_only = sum(evidence_group == "copykat_only"),
  n_infercnv_only = sum(evidence_group == "infercnv_only"),
  n_neither_or_not_evaluable = sum(evidence_group == "neither_or_not_evaluable"),
  fraction_dual_support = mean(evidence_group == "dual_support"),
  fraction_copykat_only = mean(evidence_group == "copykat_only"),
  fraction_infercnv_only = mean(evidence_group == "infercnv_only"),
  fraction_neither_or_not_evaluable =
    mean(evidence_group == "neither_or_not_evaluable"),
  median_library_size = as.numeric(median(library_size, na.rm = TRUE)),
  median_nFeature_RNA = as.numeric(median(nFeature_RNA, na.rm = TRUE)),
  median_percent_mt = as.numeric(median(percent.mt, na.rm = TRUE)),
  median_epithelial_score =
    as.numeric(median(epithelial_score, na.rm = TRUE))
), by = .(dataset_id, patient_id, final_cluster)]
cluster[, cnv_discordance_class := fcase(
  n_cells < 20L, "INSUFFICIENT_CELLS",
  fraction_dual_support >= .5, "CONCORDANT_CNV_SUPPORT",
  fraction_copykat_only >= .5, "COPYKAT_DOMINANT",
  fraction_infercnv_only >= .5, "INFERCNV_DOMINANT",
  fraction_neither_or_not_evaluable >= .5, "LOW_CNV_OR_DIPLOID_LIKE",
  default = "METHOD_DISCORDANT"
)]
cluster[, `:=`(
  dual_method_fraction = fraction_dual_support,
  copykat_only_fraction = fraction_copykat_only,
  infercnv_only_fraction = fraction_infercnv_only,
  neither_fraction = fraction_neither_or_not_evaluable,
  median_nCount = median_library_size,
  median_nFeature = median_nFeature_RNA
)]
setorder(cluster, patient_id, -fraction_copykat_only, final_cluster)
fwrite(cluster, file.path(out_dir, "cnv_evidence_by_epithelial_cluster.csv"),
       na = "NA")

t77 <- copy(cluster[patient_id == "T77"])
t77[, copykat_only_rank := frank(-fraction_copykat_only, ties.method = "min")]
t77[, t77_copykat_only_concentration :=
      if (sum(n_copykat_only) > 0) {
        n_copykat_only / sum(n_copykat_only)
      } else {
        rep(NA_real_, .N)
      }]
fwrite(t77, file.path(out_dir, "T77_cnv_discordance_audit.csv"), na = "NA")

plot_dt <- melt(
  cluster,
  id.vars = c("patient_id", "final_cluster", "cnv_discordance_class"),
  measure.vars = c(
    "fraction_dual_support", "fraction_copykat_only",
    "fraction_infercnv_only", "fraction_neither_or_not_evaluable"
  ),
  variable.name = "evidence", value.name = "fraction"
)
p <- ggplot(
  plot_dt,
  aes(reorder(final_cluster, fraction), fraction, fill = evidence)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~patient_id, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Calibrated CNV evidence by epithelial cluster",
    x = "final cluster", y = "cell fraction", fill = "evidence"
  )
ggsave(file.path(out_dir, "cnv_evidence_by_patient_cluster.png"),
       p, width = 13, height = 10, dpi = 180)
message("CNV method-discordance audit complete")
