options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v6 <- file.path(data_root, "diagnostics_v6_malignant_receiver_validation")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v6, "external_scrna_evidence_matrix_v3.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

roles <- data.table(
  dataset_id = c("GSE154600", "GSE158722", "GSE147082", "GSE151214", "GSE154763"),
  dataset_role = c(
    "primary_tumor_ecosystem", "malignant_fluid_tumor_ecosystem",
    "tumor_sensitivity_validation", "normal_fallopian_tube_reference",
    "author_annotated_myeloid_reference"
  ),
  evidence_level = c("primary", "secondary", "secondary", "reference", "reference")
)
state <- fread(file.path(v6, "state_detection_by_patient.csv"))
sens <- unique(fread(file.path(v6, "state_threshold_sensitivity.csv"))[
  minimum_macrophages == 20 &
    transcript_positive_fraction == .10 & program_cell_fraction == .10,
  .(dataset_id, SPP1_threshold_robustness)
])
mal <- fread(file.path(v6, "malignancy_summary_by_patient.csv"))
recv <- fread(file.path(v6, "malignant_epithelial_receiver_context.csv"))

rep_status <- function(n_eval, n_positive, positive_label, single_label,
                       absent_label = "LOW_OR_NOT_DETECTED") {
  if (!n_eval) return("NOT_EVALUABLE")
  if (n_positive >= 2) return(positive_label)
  if (n_positive == 1) return(single_label)
  absent_label
}

state_ds <- state[, {
  ev <- evaluable == TRUE
  n_eval <- sum(ev)
  spp1_tx <- sum(ev & SPP1_transcript_detection == "TRANSCRIPT_DETECTED")
  spp1_prog <- sum(ev & SPP1_program_support == "PROGRAM_SUPPORTED")
  spp1_joint <- sum(ev & SPP1_joint_status == "TRANSCRIPT_AND_PROGRAM_PRESENT")
  c1core <- sum(ev & C1QC_core_detection == "CORE_DETECTED")
  c1prog <- sum(ev & C1QC_program_support == "PROGRAM_SUPPORTED")
  ftx <- sum(ev & FOLR2_transcript_detection == "TRANSCRIPT_DETECTED")
  fprog <- sum(ev & FOLR2_program_support == "PROGRAM_SUPPORTED")
  list(
    n_high_confidence_macrophages = sum(n_macrophages),
    n_evaluable_patients = n_eval,
    SPP1_transcript_detection_status = rep_status(
      n_eval, spp1_tx, "DETECTED_REPLICATED", "DETECTED_SINGLE_PATIENT"
    ),
    SPP1_program_support_status = rep_status(
      n_eval, spp1_prog, "SUPPORTED_REPLICATED", "SUPPORTED_SINGLE_PATIENT",
      "PROGRAM_LOW"
    ),
    SPP1_cross_patient_reproducibility = rep_status(
      n_eval, spp1_joint, "REPLICATED", "SUPPORTIVE_SINGLE_PATIENT", "NOT_REPLICATED"
    ),
    SPP1_relative_enrichment_status = if (!n_eval) "NOT_EVALUABLE" else
      if (sum(ev & SPP1_relative_enrichment == "RELATIVELY_ENRICHED") >= 2)
        "RELATIVELY_ENRICHED_MULTIPLE_PATIENTS" else
          if (any(ev & SPP1_relative_enrichment == "RELATIVELY_ENRICHED"))
            "RELATIVELY_ENRICHED_SINGLE_PATIENT" else "NOT_RELATIVELY_ENRICHED",
    C1QC_core_detection_status = rep_status(
      n_eval, c1core, "DETECTED_REPLICATED", "DETECTED_SINGLE_PATIENT"
    ),
    C1QC_program_support_status = rep_status(
      n_eval, c1prog, "SUPPORTED_REPLICATED", "SUPPORTED_SINGLE_PATIENT",
      "PROGRAM_LOW"
    ),
    FOLR2_transcript_detection_status = rep_status(
      n_eval, ftx, "DETECTED_REPLICATED", "DETECTED_SINGLE_PATIENT"
    ),
    FOLR2_program_support_status = rep_status(
      n_eval, fprog, "SUPPORTED_REPLICATED", "SUPPORTED_SINGLE_PATIENT",
      "PROGRAM_LOW"
    )
  )
}, by = dataset_id]
state_ds <- merge(state_ds, sens, by = "dataset_id", all.x = TRUE)

mal_ds <- mal[, .(
  malignant_epithelial_available = if (sum(n_malignant_high_confidence) > 0)
    "MALIGNANT_HIGH_CONFIDENCE_AVAILABLE" else
      if (sum(n_malignant_supportive) > 0) "MALIGNANT_SUPPORTIVE_SINGLE_METHOD" else
        "NOT_AVAILABLE",
  malignancy_method = paste(unique(malignancy_method), collapse = ";"),
  n_malignant_high_confidence = sum(n_malignant_high_confidence),
  n_malignant_supportive = sum(n_malignant_supportive)
), by = dataset_id]

aggregate_receiver <- function(d, field) {
  z <- d[[field]]
  if (any(z == "SUPPORTED", na.rm = TRUE)) "SUPPORTED" else
    if (any(z == "DETECTED_LOW", na.rm = TRUE)) "DETECTED_LOW" else
      if (all(z == "NOT_EVALUABLE" | is.na(z))) "NOT_EVALUABLE" else "NOT_DETECTED"
}
recv_ds <- recv[receiver_tier == "malignant_high_confidence", {
  list(
    CD44_receiver_status = aggregate_receiver(.SD, "CD44_receiver_support"),
    ITGB1_receiver_status = aggregate_receiver(.SD, "ITGB1_receiver_support"),
    ITGB1_alpha_partner_status = aggregate_receiver(.SD, "ITGB1_alpha_partner_support"),
    dual_receiver_status = aggregate_receiver(.SD, "dual_CD44_ITGB1_support")
  )
}, by = dataset_id]
supportive <- recv[receiver_tier == "malignant_supportive", {
  list(
    CD44_supportive = aggregate_receiver(.SD, "CD44_receiver_support"),
    ITGB1_supportive = aggregate_receiver(.SD, "ITGB1_receiver_support"),
    alpha_supportive = aggregate_receiver(.SD, "ITGB1_alpha_partner_support"),
    dual_supportive = aggregate_receiver(.SD, "dual_CD44_ITGB1_support")
  )
}, by = dataset_id]
recv_ds <- merge(recv_ds, supportive, by = "dataset_id", all = TRUE)
recv_ds[
  CD44_receiver_status == "NOT_EVALUABLE" & CD44_supportive != "NOT_EVALUABLE",
  `:=`(
    CD44_receiver_status = paste0(CD44_supportive, "_MALIGNANT_SUPPORTIVE_TIER"),
    ITGB1_receiver_status = paste0(ITGB1_supportive, "_MALIGNANT_SUPPORTIVE_TIER"),
    ITGB1_alpha_partner_status = paste0(alpha_supportive, "_MALIGNANT_SUPPORTIVE_TIER"),
    dual_receiver_status = paste0(dual_supportive, "_MALIGNANT_SUPPORTIVE_TIER")
  )
]
recv_ds[, c("CD44_supportive", "ITGB1_supportive", "alpha_supportive", "dual_supportive") := NULL]

ans <- Reduce(function(x, y) merge(x, y, by = "dataset_id", all.x = TRUE),
              list(roles, state_ds, mal_ds, recv_ds))
ans[is.na(malignant_epithelial_available), `:=`(
  malignant_epithelial_available = "NOT_APPLICABLE_REFERENCE",
  malignancy_method = "NOT_APPLICABLE",
  n_malignant_high_confidence = 0L,
  n_malignant_supportive = 0L
)]
for (col in c(
  "CD44_receiver_status", "ITGB1_receiver_status",
  "ITGB1_alpha_partner_status", "dual_receiver_status"
)) ans[is.na(get(col)), (col) := "NOT_EVALUABLE"]
ans[, sender_receiver_context_status := fcase(
  SPP1_cross_patient_reproducibility == "REPLICATED" &
    grepl("MALIGNANT_SUPPORTIVE", malignant_epithelial_available) &
    grepl("SUPPORTED", ITGB1_receiver_status) &
    grepl("SUPPORTED", ITGB1_alpha_partner_status),
  "SUPPORTED_WITH_SINGLE_METHOD_MALIGNANCY",
  grepl("reference", dataset_role, ignore.case = TRUE),
  "REFERENCE_BACKGROUND_ONLY",
  default = "NOT_EVALUABLE_OR_NOT_ESTABLISHED"
)]
ans[, tumor_specificity_status := "NOT_ESTABLISHED"]
ans[, main_limitation := fcase(
  dataset_id == "GSE154600",
  "CopyKAT supports aneuploid epithelial cells, but inferCNV is unavailable; no two-method high-confidence malignancy consensus",
  dataset_id == "GSE147082",
  "PT-2834 cluster 4/7 are CopyKAT diploid-like; only 15 broad epithelial candidates and no malignant receiver",
  dataset_id == "GSE158722",
  "platform identity is unavailable, so same-platform formal malignancy reference and malignant receiver are not evaluable",
  dataset_id == "GSE151214",
  "normal fallopian-tube reference only; excluded from tumor effect aggregation",
  dataset_id == "GSE154763",
  "author-normalized, author-annotation-driven myeloid reference; no raw-count tumor receiver analysis"
)]
setcolorder(ans, c(
  "dataset_id", "dataset_role", "evidence_level",
  "n_high_confidence_macrophages", "n_evaluable_patients",
  "SPP1_transcript_detection_status", "SPP1_program_support_status",
  "SPP1_cross_patient_reproducibility", "SPP1_relative_enrichment_status",
  "SPP1_threshold_robustness", "C1QC_core_detection_status",
  "C1QC_program_support_status", "FOLR2_transcript_detection_status",
  "FOLR2_program_support_status", "malignant_epithelial_available",
  "malignancy_method", "n_malignant_high_confidence", "n_malignant_supportive",
  "CD44_receiver_status", "ITGB1_receiver_status",
  "ITGB1_alpha_partner_status", "dual_receiver_status",
  "sender_receiver_context_status", "tumor_specificity_status", "main_limitation"
))
fwrite(ans, out, na = "NA")
message("External scRNA evidence matrix v3 complete")
