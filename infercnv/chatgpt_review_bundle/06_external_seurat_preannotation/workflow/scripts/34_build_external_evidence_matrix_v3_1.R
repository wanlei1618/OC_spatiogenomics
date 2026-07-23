options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v6 <- file.path(data_root, "diagnostics_v6_malignant_receiver_validation")
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v61, "external_scrna_evidence_matrix_v3_1.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

m <- fread(file.path(v6, "external_scrna_evidence_matrix_v3.csv"))
drop <- intersect(c(
  "n_malignant_high_confidence", "n_malignant_supportive",
  "CD44_receiver_status", "ITGB1_receiver_status",
  "ITGB1_alpha_partner_status", "dual_receiver_status"
), names(m))
m[, (drop) := NULL]
new_count_cols <- c(
  "n_final_epithelial", "n_copykat_available", "n_copykat_submitted",
  "n_copykat_defined_any", "n_stable_malignant_supportive",
  "n_stable_diploid_like", "n_copykat_unstable",
  "n_copykat_mostly_not_defined"
)
for (col in new_count_cols) m[, (col) := NA_integer_]
m[, `:=`(
  copykat_target_coverage = NA_real_,
  copykat_defined_call_rate = NA_real_,
  copykat_stability_status = "NOT_APPLICABLE",
  CD44_receiver_status_stable_malignant = "NOT_EVALUABLE",
  ITGB1_receiver_status_stable_malignant = "NOT_EVALUABLE",
  ITGB1_alpha_partner_status_stable_malignant = "NOT_EVALUABLE",
  dual_receiver_status_stable_malignant = "NOT_EVALUABLE",
  C1QC_reference_program_support_status = "NOT_APPLICABLE",
  FOLR2_reference_program_support_status = "NOT_APPLICABLE"
)]

st <- fread(file.path(v61, "GSE154600_copykat_stability_summary.csv"))[1]
rec <- fread(file.path(v61, "stable_malignant_receiver_context_summary.csv"))[
  receiver_tier == "stable_malignant_supportive"
][1]
i <- which(m$dataset_id == "GSE154600")
m[i, `:=`(
  n_final_epithelial = st$n_final_epithelial,
  n_copykat_available = st$n_available_for_copykat,
  n_copykat_submitted = st$n_submitted_any_run,
  n_copykat_defined_any = st$n_defined_in_at_least_1_run,
  n_stable_malignant_supportive = st$n_stable_aneuploid,
  n_stable_diploid_like = st$n_stable_diploid,
  n_copykat_unstable = st$n_unstable,
  n_copykat_mostly_not_defined = st$n_mostly_not_defined,
  copykat_target_coverage = st$target_input_coverage,
  copykat_defined_call_rate = st$any_defined_call_rate,
  copykat_stability_status = if (st$n_stable_aneuploid > 0)
    "STABLE_SINGLE_METHOD_EVIDENCE_AVAILABLE" else
      if (st$n_unstable > 0) "MALIGNANCY_SINGLE_METHOD_UNSTABLE" else "NOT_EVALUABLE",
  malignant_epithelial_available = if (st$n_stable_aneuploid > 0)
    "MALIGNANT_SUPPORTIVE_STABLE_SINGLE_METHOD" else
      if (st$n_unstable > 0) "MALIGNANCY_SINGLE_METHOD_UNSTABLE" else "NOT_AVAILABLE",
  malignancy_method = "CopyKAT_three_reference_seeds_single_method",
  CD44_receiver_status_stable_malignant = rec$CD44_receiver_status,
  ITGB1_receiver_status_stable_malignant = rec$ITGB1_receiver_status,
  ITGB1_alpha_partner_status_stable_malignant = rec$ITGB1_alpha_partner_status,
  dual_receiver_status_stable_malignant = rec$dual_receiver_status,
  sender_receiver_context_status = if (
    rec$ITGB1_receiver_status == "SUPPORTED" &
      rec$ITGB1_alpha_partner_status == "SUPPORTED"
  ) "SUPPORTED_STABLE_SINGLE_METHOD_MALIGNANCY" else
    "PARTIAL_RECEIVER_SUPPORT_STABLE_SINGLE_METHOD_MALIGNANCY",
  main_limitation = paste(
    "Stable CopyKAT calls remain single-method evidence without inferCNV;",
    "some final epithelial cells are absent from all four lineage count inputs"
  )
)]

ref <- fread(file.path(v61, "GSE154763_reference_state_summary_v2.csv"))
ev <- ref[evaluable == TRUE]
replicate_status <- function(x, positive) {
  n <- sum(x == positive)
  if (n >= 2) "DETECTED_REPLICATED" else if (n == 1) "DETECTED_SINGLE_PATIENT" else "LOW_OR_NOT_DETECTED"
}
program_status <- function(x) {
  n <- sum(x == "REFERENCE_PROGRAM_SUPPORTED")
  if (n >= 2) "SUPPORTED_REPLICATED_REFERENCE" else
    if (n == 1) "SUPPORTED_SINGLE_REFERENCE_PATIENT" else "REFERENCE_PROGRAM_LOW"
}
j <- which(m$dataset_id == "GSE154763")
m[j, `:=`(
  SPP1_transcript_detection_status =
    replicate_status(ev$SPP1_transcript_detection, "TRANSCRIPT_DETECTED"),
  SPP1_program_support_status = program_status(ev$SPP1_reference_program_support),
  SPP1_threshold_robustness = unique(ev$SPP1_reference_threshold_robustness)[1],
  C1QC_core_detection_status =
    replicate_status(ev$C1QC_core_transcript_detection, "CORE_TRANSCRIPT_DETECTED"),
  C1QC_program_support_status = program_status(ev$C1QC_reference_program_support),
  FOLR2_transcript_detection_status =
    replicate_status(ev$FOLR2_transcript_detection, "TRANSCRIPT_DETECTED"),
  FOLR2_program_support_status = program_status(ev$FOLR2_reference_program_support),
  C1QC_reference_program_support_status = program_status(ev$C1QC_reference_program_support),
  FOLR2_reference_program_support_status = program_status(ev$FOLR2_reference_program_support),
  main_limitation = paste(
    "Author-normalized myeloid reference; transcript detection uses explicit expression > 0,",
    "while reference program support is reported separately"
  )
)]
m[, tumor_specificity_status := "NOT_ESTABLISHED"]
fwrite(m, out, na = "NA")
message("External scRNA evidence matrix v3.1 complete")
