options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
out_dir <- file.path(data_root, "research_spatial_transition")

fpr <- fread(file.path(out_dir, "infercnvpy_negative_control_fpr_by_patient.csv"))
cnv_patient <- fread(file.path(out_dir, "GSE154600_calibrated_cnv_by_patient.csv"))
t77 <- fread(file.path(out_dir, "T77_cnv_discordance_audit.csv"))
receiver <- fread(file.path(out_dir, "receiver_reproducibility_calibrated.csv"))
manifest <- fread(file.path(out_dir, "spatial_dataset_pilot_manifest.csv"))
spatial <- fread(file.path(out_dir, "spatial_permutation_results.csv"))

first <- spatial[adjacency_order == "first_order" &
                   !is.na(observed_expected_ratio) & !is.na(empirical_p)]
classify_spatial <- function(primary_test, control_test) {
  x <- first[spatial_test == primary_test]
  ctrl <- first[spatial_test == control_test]
  if (nrow(x) < 2L) return("NOT_EVALUABLE")
  n_support <- sum(x$observed_expected_ratio > 1 & x$empirical_p < .05)
  matched <- merge(
    x[, .(dataset_id, sample_id, primary_oe = observed_expected_ratio)],
    ctrl[, .(dataset_id, sample_id, control_oe = observed_expected_ratio)],
    by = c("dataset_id", "sample_id")
  )
  not_stronger <- nrow(matched) >= 2L &&
    median(matched$primary_oe, na.rm = TRUE) <=
      median(matched$control_oe, na.rm = TRUE)
  if (not_stronger) {
    "NO_SPECIFIC_ENRICHMENT"
  } else if (n_support >= 2L) {
    "REPLICATED_SPATIAL_SUPPORT"
  } else if (n_support == 1L) {
    "SINGLE_SAMPLE_SUPPORT"
  } else if (sum(x$observed_expected_ratio <= 1, na.rm = TRUE) >
             nrow(x) / 2) {
    "NOT_SUPPORTED"
  } else {
    "NO_SPECIFIC_ENRICHMENT"
  }
}

itgb1_status <- classify_spatial(
  "SPP1_macrophage_to_ITGB1_receiver",
  "C1QC_macrophage_to_ITGB1_receiver"
)
cd44_status <- classify_spatial(
  "SPP1_macrophage_to_CD44_receiver",
  "all_macrophage_to_epithelial"
)
itgb1_receiver <- receiver[receiver_metric == "ITGB1",
                           receiver_reproducibility_status][1]
cd44_receiver <- receiver[receiver_metric == "CD44",
                          receiver_reproducibility_status][1]
rank_status <- function(x) {
  match(x, c(
    "NOT_EVALUABLE", "NOT_SUPPORTED", "NO_SPECIFIC_ENRICHMENT",
    "SINGLE_SAMPLE_SUPPORT", "REPLICATED_SPATIAL_SUPPORT"
  ))
}
priority <- if (rank_status(itgb1_status) >= rank_status(cd44_status)) {
  "ITGB1"
} else {
  "CD44"
}
worth_wetlab <- any(c(itgb1_status, cd44_status) %in%
                      c("REPLICATED_SPATIAL_SUPPORT", "SINGLE_SAMPLE_SUPPORT")) &&
  any(c(itgb1_receiver, cd44_receiver) %in%
        c("REPLICATED_ROBUST", "REPLICATED_BUT_DEPTH_SENSITIVE"))

matrix <- rbindlist(list(
  data.table(
    evidence_domain = "held-out infercnvpy calibration",
    target = fpr$patient_id,
    evidence_status = fpr$infercnv_threshold_status,
    quantitative_summary = sprintf(
      "median test FPR=%.4f; %d/3 splits <=0.05",
      fpr$median_test_control_fpr, fpr$n_splits_fpr_le_0_05
    ),
    interpretation = fifelse(
      fpr$fpr_calibration_pass,
      "calibrated infercnvpy usable in the primary integrated tier",
      "infercnvpy retained only as sensitivity evidence"
    )
  ),
  receiver[, .(
    evidence_domain = "receiver expression robustness",
    target = receiver_metric,
    evidence_status = receiver_reproducibility_status,
    quantitative_summary = sprintf(
      "%d evaluable patients; %d with detection fraction >=0.10",
      n_patients_evaluable, n_patients_detection_fraction_ge_0_10
    ),
    interpretation
  )],
  data.table(
    evidence_domain = "spatial adjacency pilot",
    target = c("SPP1_to_ITGB1", "SPP1_to_CD44"),
    evidence_status = c(itgb1_status, cd44_status),
    quantitative_summary = c(
      sprintf("%d/%d samples OE>1 and empirical p<0.05",
              first[
                spatial_test == "SPP1_macrophage_to_ITGB1_receiver",
                sum(observed_expected_ratio > 1 & empirical_p < .05)
              ],
              first[spatial_test ==
                      "SPP1_macrophage_to_ITGB1_receiver", .N]),
      sprintf("%d/%d samples OE>1 and empirical p<0.05",
              first[
                spatial_test == "SPP1_macrophage_to_CD44_receiver",
                sum(observed_expected_ratio > 1 & empirical_p < .05)
              ],
              first[spatial_test ==
                      "SPP1_macrophage_to_CD44_receiver", .N])
    ),
    interpretation = paste(
      "Spatial proximity is associative and does not establish direct",
      "SPP1-receptor binding or causal signaling."
    )
  )
), fill = TRUE)
fwrite(matrix, file.path(out_dir, "spatial_validation_evidence_matrix.csv"),
       na = "NA")

selected <- manifest[selected_for_pilot == TRUE,
                     paste0(dataset_id, "/", sample_id)]
top_t77 <- head(t77[order(-n_copykat_only)], 5L)
fmt_spatial <- function(test) {
  x <- first[spatial_test == test]
  if (!nrow(x)) return("not evaluable")
  paste(sprintf(
    "%s/%s OE=%.3f, p=%.4f",
    x$dataset_id, x$sample_id, x$observed_expected_ratio, x$empirical_p
  ), collapse = "; ")
}
fpr_lines <- paste(sprintf(
  "%s: median FPR %.4f (%d/3 splits <=0.05), %s",
  fpr$patient_id, fpr$median_test_control_fpr, fpr$n_splits_fpr_le_0_05,
  fpr$infercnv_threshold_status
), collapse = "\n- ")
t77_lines <- paste(sprintf(
  "%s: %d CopyKAT-only cells (%.1f%% of cluster; %.1f%% of T77 CopyKAT-only)",
  top_t77$final_cluster, top_t77$n_copykat_only,
  100 * top_t77$fraction_copykat_only,
  100 * top_t77$t77_copykat_only_concentration
), collapse = "\n- ")
receiver_lines <- paste(sprintf(
  "%s: %s",
  receiver$receiver_metric, receiver$receiver_reproducibility_status
), collapse = "\n- ")

report <- c(
  "# Final calibrated CNV and spatial pilot decision",
  "",
  "## 1. Held-out immune negative-control FPR",
  "",
  paste0("- ", fpr_lines),
  "",
  "## 2. T77 CNV method discordance",
  "",
  paste0("- ", t77_lines),
  "",
  paste0(
    "The cluster table separates sequencing-depth and epithelial-marker ",
    "summaries from method calls; it does not reinterpret fixed clusters. ",
    "Clusters 4, 5 and 6 account for ",
    sprintf("%.1f%%", 100 * sum(
      t77[final_cluster %chin% c(
        "Epithelial__A_uncorrected__4",
        "Epithelial__A_uncorrected__5",
        "Epithelial__A_uncorrected__6"
      ), n_copykat_only]
    ) / sum(t77$n_copykat_only)),
    " of T77 CopyKAT-only cells. Clusters 4 and 5 retain higher median depth ",
    "and epithelial scores than clusters 6 and 7, so the discordance is ",
    "cluster-associated and cannot be attributed uniformly to low depth."
  ),
  "",
  "## 3. Receiver pseudobulk and depth strata",
  "",
  paste0("- ", receiver_lines),
  "",
  paste0(
    "The primary tier uses only CALIBRATED_DUAL_METHOD_SUPPORT cells. ",
    "CopyKAT-stable cells are reported separately as sensitivity evidence."
  ),
  "",
  "## 4. Spatial datasets entering the pilot",
  "",
  paste0("- ", selected),
  "",
  "## 5. SPP1 macrophage proximity to ITGB1 receiver",
  "",
  paste0("Status: **", itgb1_status, "**. ", fmt_spatial(
    "SPP1_macrophage_to_ITGB1_receiver"
  )),
  "",
  "## 6. SPP1 macrophage proximity to CD44 receiver",
  "",
  paste0("Status: **", cd44_status, "**. ", fmt_spatial(
    "SPP1_macrophage_to_CD44_receiver"
  )),
  "",
  "## 7. Comparison with C1QC macrophage control",
  "",
  fmt_spatial("C1QC_macrophage_to_ITGB1_receiver"),
  "",
  paste0(
    "SPP1-to-ITGB1 is stronger than the matched C1QC control in two of the ",
    "three SPP1-supporting samples, but C1QC is stronger in GSM6177617. ",
    "The decision status explicitly downgrades an SPP1 result if it is not ",
    "stronger than the matched C1QC or general-macrophage control."
  ),
  "",
  "## 8. Wet-lab decision",
  "",
  paste0(
    "Proceed to a focused wet-lab follow-up: **",
    if (worth_wetlab) "yes" else "not yet", "**. Priority receiver: **",
    priority, "**."
  ),
  "",
  "## 9. Interpretation boundary",
  "",
  paste0(
    "The supported wording is an SPP1-macrophage spatial proximity to a ",
    "CNV-supported epithelial ITGB1/CD44 expression context. The analysis ",
    "does not prove direct SPP1 binding, receptor activation, causal tumor ",
    "progression, or tumor specificity of SPP1 macrophages."
  ),
  "",
  "## 10. Remaining limitations",
  "",
  paste0(
    "Visium spots are mixtures, inferred signature scores are not cell ",
    "identity calls, samples rather than spots are biological replicates, ",
    "and the pilot contains only the locally available technically usable ",
    "datasets. Spatial association requires orthogonal validation."
  )
)
writeLines(report, file.path(out_dir, "FINAL_SPATIAL_PILOT_DECISION.md"),
           useBytes = TRUE)
message("Final calibrated CNV and spatial pilot decision generated")
