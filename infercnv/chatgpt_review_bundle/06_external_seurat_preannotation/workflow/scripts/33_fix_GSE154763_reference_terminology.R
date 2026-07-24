options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v4 <- file.path(data_root, "diagnostics_v4_cross_dataset_validation")
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v61, "GSE154763_reference_state_summary_v2.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

x <- fread(file.path(v4, "GSE154763_refined", "author_annotation_with_refined_state.csv.gz"))
x <- x[final_cell_type == "Macrophage"]
required <- c(
  "SPP1_expression", "C1QA_expression", "C1QB_expression",
  "C1QC_expression", "FOLR2_expression"
)
available <- required %in% names(x)
names(available) <- required
rank_pct <- function(v) frank(v, ties.method = "average") / length(v)
x[, `:=`(
  SPP1_program_pct = rank_pct(SPP1_program),
  C1QC_program_pct = rank_pct(C1QC_program),
  FOLR2_program_pct = rank_pct(FOLR2_program)
)]
x[, SPP1_transcript_detected := SPP1_expression > 0]
x[, C1QC_core_transcript_detected :=
    (C1QA_expression > 0) + (C1QB_expression > 0) + (C1QC_expression > 0) >= 2]
x[, FOLR2_transcript_detected := FOLR2_expression > 0]
x[, SPP1_reference_program_cell :=
    grepl("SPP1", cell_type_original, ignore.case = TRUE) | SPP1_program_pct >= .75]
x[, C1QC_reference_program_cell :=
    grepl("C1Q", cell_type_original, ignore.case = TRUE) | C1QC_program_pct >= .75]
x[, FOLR2_reference_program_cell :=
    grepl("FOLR2", cell_type_original, ignore.case = TRUE) | FOLR2_program_pct >= .75]

summary <- x[, .(
  n_macrophages = .N,
  expression_basis = "author_normalized_expression",
  SPP1_transcript_positive_fraction = mean(SPP1_transcript_detected),
  SPP1_reference_program_fraction = mean(SPP1_reference_program_cell),
  C1QC_core_transcript_positive_fraction = mean(C1QC_core_transcript_detected),
  C1QC_reference_program_fraction = mean(C1QC_reference_program_cell),
  FOLR2_transcript_positive_fraction = mean(FOLR2_transcript_detected),
  FOLR2_reference_program_fraction = mean(FOLR2_reference_program_cell)
), by = .(
  patient_id = as.character(patient), sample_id = as.character(library_id)
)]
summary[, dataset_id := "GSE154763"]
setcolorder(summary, c("dataset_id", "patient_id", "sample_id"))
summary[, evaluable := n_macrophages >= 20]
summary[, SPP1_transcript_detection := fcase(
  !evaluable, "NOT_EVALUABLE",
  SPP1_transcript_positive_fraction >= .10, "TRANSCRIPT_DETECTED",
  default = "LOW_OR_NOT_DETECTED"
)]
summary[, SPP1_reference_program_support := fcase(
  !evaluable, "NOT_EVALUABLE",
  SPP1_reference_program_fraction >= .10, "REFERENCE_PROGRAM_SUPPORTED",
  default = "REFERENCE_PROGRAM_LOW"
)]
summary[, C1QC_core_transcript_detection := fcase(
  !evaluable, "NOT_EVALUABLE",
  C1QC_core_transcript_positive_fraction >= .10, "CORE_TRANSCRIPT_DETECTED",
  default = "LOW_OR_NOT_DETECTED"
)]
summary[, C1QC_reference_program_support := fcase(
  !evaluable, "NOT_EVALUABLE",
  C1QC_reference_program_fraction >= .10, "REFERENCE_PROGRAM_SUPPORTED",
  default = "REFERENCE_PROGRAM_LOW"
)]
summary[, FOLR2_transcript_detection := fcase(
  !evaluable, "NOT_EVALUABLE",
  FOLR2_transcript_positive_fraction >= .10, "TRANSCRIPT_DETECTED",
  default = "LOW_OR_NOT_DETECTED"
)]
summary[, FOLR2_reference_program_support := fcase(
  !evaluable, "NOT_EVALUABLE",
  FOLR2_reference_program_fraction >= .10, "REFERENCE_PROGRAM_SUPPORTED",
  default = "REFERENCE_PROGRAM_LOW"
)]
thresholds <- c(.05, .10, .20)
consistent <- all(vapply(thresholds, function(th) {
  sum(summary$evaluable &
        summary$SPP1_transcript_positive_fraction >= th &
        summary$SPP1_reference_program_fraction >= th) >= 2
}, logical(1)))
summary[, SPP1_reference_threshold_robustness :=
           if (sum(evaluable) == 0) "reference_not_evaluable" else
             if (consistent) "reference_consistent" else "reference_threshold_sensitive"]
fwrite(summary, out, na = "NA")

audit <- data.table(
  state = c("SPP1", "C1QC", "FOLR2"),
  transcript_expression_columns = c(
    "SPP1_expression",
    "C1QA_expression;C1QB_expression;C1QC_expression",
    "FOLR2_expression"
  ),
  expression_columns_available = c(
    available[["SPP1_expression"]],
    all(available[c("C1QA_expression", "C1QB_expression", "C1QC_expression")]),
    available[["FOLR2_expression"]]
  ),
  transcript_detection_definition = c(
    "SPP1_expression > 0",
    "at least 2 of C1QA/C1QB/C1QC expression > 0",
    "FOLR2_expression > 0"
  ),
  reference_program_field = c(
    "SPP1_reference_program_support",
    "C1QC_reference_program_support",
    "FOLR2_reference_program_support"
  ),
  terminology_decision = c(
    "true transcript detection plus separate reference program support",
    "true core transcript detection plus separate reference program support",
    "true transcript detection plus separate reference program support"
  )
)
fwrite(audit, file.path(v61, "GSE154763_reference_terminology_audit.csv"), na = "NA")
message("GSE154763 reference terminology and direct-expression detection fixed")
