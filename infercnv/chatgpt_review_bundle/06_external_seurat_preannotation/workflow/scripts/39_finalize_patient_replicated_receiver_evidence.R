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
out_dir <- file.path(data_root, "research_validation_independent_cnv")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(out_dir, "receiver_evidence_by_patient_final.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

counts <- readRDS(file.path(
  out_dir, "GSE154600_complete_final_epithelial_counts.rds"
))
consensus <- fread(file.path(
  out_dir, "GSE154600_copykat_infercnv_consensus_by_cell.csv.gz"
))
bias <- fread(file.path(out_dir, "copykat_defined_bias_summary.csv"))
bias <- bias[patient_id != "ALL"]
genes <- c("CD44", "ITGB1", "ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9")
raw <- counts[intersect(genes, rownames(counts)), , drop = FALSE]
det <- function(g) {
  if (g %in% rownames(raw)) as.numeric(raw[g, ] > 0) else
    rep(NA_real_, ncol(raw))
}
cell <- data.table(cell_id = colnames(counts))
for (g in genes) cell[[paste0(g, "_positive")]] <- det(g)
cell <- merge(
  consensus[, .(
    patient_id, cell_id, stability_class, infercnv_status,
    integrated_cnv_evidence
  )],
  cell, by = "cell_id", all.x = TRUE
)

support <- function(n, fraction) fcase(
  n < 20, "NOT_EVALUABLE",
  is.na(fraction) | fraction == 0, "NOT_DETECTED",
  fraction >= .10, "SUPPORTED",
  default = "DETECTED_LOW"
)
alpha <- c("ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9")
receiver <- rbindlist(lapply(unique(cell$patient_id), function(pt) {
  d <- cell[patient_id == pt]
  dual <- d[integrated_cnv_evidence == "DUAL_METHOD_MALIGNANT_SUPPORT"]
  if (nrow(dual) >= 20L) {
    z <- dual
    tier <- "DUAL_METHOD_MALIGNANT_SUPPORT"
    tier_method <- "CopyKAT_stable_aneuploid_plus_inferCNV_high"
  } else {
    z <- d[stability_class == "STABLE_ANEUPLOID"]
    tier <- "STABLE_COPYKAT_ANEUPLOID_SINGLE_METHOD_FALLBACK"
    tier_method <- "CopyKAT_stable_aneuploid_single_method"
  }
  alpha_cols <- paste0(alpha, "_positive")
  alpha_mat <- as.matrix(z[, ..alpha_cols])
  any_alpha <- if (nrow(z)) rowSums(alpha_mat > 0, na.rm = TRUE) > 0 else logical()
  alpha_frac <- if (nrow(z)) colMeans(alpha_mat > 0, na.rm = TRUE) else
    setNames(rep(NA_real_, length(alpha)), alpha_cols)
  cd44 <- if (nrow(z)) mean(z$CD44_positive > 0, na.rm = TRUE) else NA_real_
  itgb1 <- if (nrow(z)) mean(z$ITGB1_positive > 0, na.rm = TRUE) else NA_real_
  itgb1_alpha <- if (nrow(z))
    mean(z$ITGB1_positive > 0 & any_alpha, na.rm = TRUE) else NA_real_
  dual_frac <- if (nrow(z))
    mean(z$CD44_positive > 0 & z$ITGB1_positive > 0, na.rm = TRUE) else NA_real_
  b <- bias[patient_id == pt][1]
  data.table(
    dataset_id = "GSE154600",
    patient_id = pt,
    receiver_tier = tier,
    receiver_tier_method = tier_method,
    n_dual_method_cells_available = nrow(dual),
    n_receiver_cells = nrow(z),
    CD44_positive_fraction = cd44,
    ITGB1_positive_fraction = itgb1,
    ITGB1_any_alpha_copositive_fraction = itgb1_alpha,
    CD44_ITGB1_copositive_fraction = dual_frac,
    ITGA4_positive_fraction = alpha_frac[["ITGA4_positive"]],
    ITGA5_positive_fraction = alpha_frac[["ITGA5_positive"]],
    ITGAV_positive_fraction = alpha_frac[["ITGAV_positive"]],
    ITGA8_positive_fraction = alpha_frac[["ITGA8_positive"]],
    ITGA9_positive_fraction = alpha_frac[["ITGA9_positive"]],
    dominant_alpha_partner = if (all(is.na(alpha_frac)))
      NA_character_ else alpha[which.max(alpha_frac)],
    CD44_patient_status = support(nrow(z), cd44),
    ITGB1_patient_status = support(nrow(z), itgb1),
    ITGB1_alpha_patient_status = support(nrow(z), itgb1_alpha),
    dual_CD44_ITGB1_patient_status = support(nrow(z), dual_frac),
    potential_detection_depth_selection_bias =
      b$potential_detection_depth_selection_bias,
    receiver_bias_note = if (isTRUE(b$potential_detection_depth_selection_bias))
      "potential detection-depth selection bias" else
        "no task-threshold depth selection bias detected",
    malignant_receiver_expression_context = "DESCRIPTIVE_ONLY"
  )
}), fill = TRUE)
fwrite(receiver, out, na = "NA")

reproducibility <- function(fractions, n_cells) {
  evaluable <- n_cells >= 20 & !is.na(fractions)
  x <- fractions[evaluable]
  if (length(x) < 2L) return(list(
    status = "NOT_EVALUABLE", n_evaluable = length(x),
    n_supported = sum(x >= .10), n_low = sum(x > 0 & x < .10),
    n_not_detected = sum(x == 0), heterogeneity = FALSE
  ))
  category <- ifelse(x >= .10, "SUPPORTED",
                     ifelse(x > 0, "DETECTED_LOW", "NOT_DETECTED"))
  ns <- sum(category == "SUPPORTED")
  nl <- sum(category == "DETECTED_LOW")
  nn <- sum(category == "NOT_DETECTED")
  status <- if (ns >= 2L) "REPLICATED_SUPPORTED" else
    if (ns == 1L) "SINGLE_PATIENT_SUPPORTED" else
      if (nl >= 2L) "REPLICATED_LOW" else
        if (nn == length(x)) "NOT_DETECTED" else "HETEROGENEOUS"
  list(
    status = status, n_evaluable = length(x), n_supported = ns,
    n_low = nl, n_not_detected = nn,
    heterogeneity = length(unique(category)) > 1L
  )
}
metric_map <- c(
  CD44 = "CD44_positive_fraction",
  ITGB1 = "ITGB1_positive_fraction",
  ITGB1_alpha = "ITGB1_any_alpha_copositive_fraction",
  dual_CD44_ITGB1 = "CD44_ITGB1_copositive_fraction"
)
summary <- rbindlist(lapply(names(metric_map), function(metric) {
  ans <- reproducibility(receiver[[metric_map[[metric]]]],
                         receiver$n_receiver_cells)
  data.table(
    dataset_id = "GSE154600",
    receiver_metric = metric,
    patient_reproducibility = ans$status,
    n_evaluable_patients = ans$n_evaluable,
    n_supported_patients = ans$n_supported,
    n_detected_low_patients = ans$n_low,
    n_not_detected_patients = ans$n_not_detected,
    patient_heterogeneity_present = ans$heterogeneity,
    malignant_receiver_expression_context = "DESCRIPTIVE_ONLY",
    evidence_limit = paste(
      "no sample-level association statistics;",
      "no spatial proximity evidence;",
      "no direct ligand-receptor or causal evidence"
    )
  )
}))
fwrite(summary, file.path(out_dir, "receiver_reproducibility_summary.csv"),
       na = "NA")
message("Patient-replicated receiver evidence finalized")
