options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
out_dir <- file.path(data_root, "research_spatial_transition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(
  data_root, "research_validation_independent_cnv",
  "GSE154600_complete_final_epithelial_counts.rds"
))
cnv <- fread(file.path(out_dir, "GSE154600_calibrated_cnv_by_cell.csv.gz"))
genes <- c("CD44", "ITGB1", "ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9")
present <- intersect(genes, rownames(counts))

tiers <- rbind(
  cnv[integrated_calibrated_cnv_evidence ==
        "CALIBRATED_DUAL_METHOD_SUPPORT",
      .(patient_id, cell_id,
        evidence_tier = "CALIBRATED_DUAL_METHOD_SUPPORT")],
  cnv[stability_class == "STABLE_ANEUPLOID",
      .(patient_id, cell_id,
        evidence_tier = "COPYKAT_STABLE_ANEUPLOID_SENSITIVITY")]
)

pb <- tiers[, {
  ids <- intersect(cell_id, colnames(counts))
  if (!length(ids)) return(data.table())
  m <- counts[present, ids, drop = FALSE]
  summed <- Matrix::rowSums(m)
  lib <- sum(summed)
  data.table(
    gene = genes,
    n_cells = length(ids),
    raw_pseudobulk_count = as.numeric(summed[match(genes, names(summed))]),
    library_size = lib,
    cpm = as.numeric(summed[match(genes, names(summed))]) / pmax(lib, 1) * 1e6,
    log2_cpm_plus1 =
      log2(as.numeric(summed[match(genes, names(summed))]) /
             pmax(lib, 1) * 1e6 + 1),
    detection_fraction = vapply(
      genes,
      function(g) if (g %in% present) mean(m[g, ] > 0) else NA_real_,
      numeric(1)
    )
  )
}, by = .(patient_id, evidence_tier)]
pb[is.na(raw_pseudobulk_count), `:=`(
  raw_pseudobulk_count = 0, cpm = 0, log2_cpm_plus1 = 0
)]
setnames(pb, "detection_fraction", "descriptive_detection_fraction")
fwrite(pb, file.path(out_dir, "receiver_pseudobulk_by_patient.csv"), na = "NA")

metric_fraction <- function(ids) {
  if (!length(ids)) {
    return(data.table(
      receiver_metric = c("CD44", "ITGB1", "ITGB1_any_alpha",
                          "CD44_ITGB1_dual"),
      descriptive_detection_fraction = NA_real_
    ))
  }
  positive <- function(g) {
    if (g %in% rownames(counts)) as.vector(counts[g, ids] > 0) else
      rep(FALSE, length(ids))
  }
  cd44 <- positive("CD44")
  itgb1 <- positive("ITGB1")
  alpha <- Reduce(`|`, lapply(c("ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9"),
                              positive))
  data.table(
    receiver_metric = c("CD44", "ITGB1", "ITGB1_any_alpha",
                        "CD44_ITGB1_dual"),
    descriptive_detection_fraction = c(
      mean(cd44), mean(itgb1), mean(itgb1 & alpha), mean(cd44 & itgb1)
    )
  )
}

primary <- tiers[evidence_tier == "CALIBRATED_DUAL_METHOD_SUPPORT"]
overall <- primary[, {
  ids <- intersect(cell_id, colnames(counts))
  ans <- metric_fraction(ids)
  ans[, primary_n_cells := length(ids)]
  ans
}, by = patient_id]
depth <- primary[, {
  ids <- intersect(cell_id, colnames(counts))
  if (length(ids) < 3L) {
    cbind(
      data.table(depth_tertile = "NOT_EVALUABLE", n_cells = length(ids)),
      metric_fraction(character())
    )
  } else {
    lib <- Matrix::colSums(counts[, ids, drop = FALSE])
    rank <- frank(lib, ties.method = "average")
    tert <- cut(
      rank, breaks = c(-Inf, length(ids) / 3, 2 * length(ids) / 3, Inf),
      labels = c("low", "middle", "high")
    )
    rbindlist(lapply(levels(tert), function(tt) {
      keep <- ids[tert == tt]
      cbind(
        data.table(depth_tertile = tt, n_cells = length(keep)),
        metric_fraction(keep)
      )
    }))
  }
}, by = patient_id]
depth <- merge(
  depth, overall,
  by = c("patient_id", "receiver_metric"), all.x = TRUE,
  suffixes = c("", "_overall")
)
depth[, depth_robustness := {
  vals <- descriptive_detection_fraction[
    depth_tertile %in% c("low", "middle", "high")]
  overall <- unique(descriptive_detection_fraction_overall)[1]
  if (is.na(overall) || length(vals) < 3L || primary_n_cells[1] < 9L) {
    "NOT_EVALUABLE"
  } else if ((overall >= .10 && sum(vals >= .10, na.rm = TRUE) >= 2L) ||
             (overall < .10 && sum(vals > 0, na.rm = TRUE) >= 2L)) {
    "DEPTH_ROBUST"
  } else if (vals[3] > 0 && sum(vals[1:2] > 0, na.rm = TRUE) == 0L) {
    "DEPTH_SENSITIVE"
  } else {
    "DEPTH_INCONCLUSIVE"
  }
}, by = .(patient_id, receiver_metric)]
fwrite(depth, file.path(out_dir, "receiver_depth_sensitivity_by_patient.csv"),
       na = "NA")

repro <- unique(depth[, .(
  patient_id, receiver_metric, primary_n_cells,
  detection_fraction = descriptive_detection_fraction_overall,
  depth_robustness
)])
summary <- repro[, {
  n_support <- sum(detection_fraction >= .10, na.rm = TRUE)
  n_detect <- sum(detection_fraction > 0, na.rm = TRUE)
  n_eval <- sum(depth_robustness != "NOT_EVALUABLE")
  status <- if (n_eval < 2L) {
    "NOT_EVALUABLE"
  } else if (n_support >= 2L &&
             all(depth_robustness[detection_fraction >= .10] ==
                   "DEPTH_ROBUST")) {
    "REPLICATED_ROBUST"
  } else if (n_support >= 2L) {
    "REPLICATED_BUT_DEPTH_SENSITIVE"
  } else if (n_support == 1L) {
    "SINGLE_PATIENT_SUPPORT"
  } else if (n_detect >= 2L) {
    "REPLICATED_LOW"
  } else {
    "NOT_DETECTED"
  }
  .(
    n_patients_evaluable = n_eval,
    n_patients_detection_fraction_ge_0_10 = n_support,
    n_patients_detected = n_detect,
    receiver_reproducibility_status = status,
    interpretation =
      "SPP1-associated malignant epithelial ITGB1/CD44 expression context"
  )
}, by = receiver_metric]
fwrite(summary,
       file.path(out_dir, "receiver_reproducibility_calibrated.csv"), na = "NA")
message("Receiver pseudobulk and depth-robustness analysis complete")
