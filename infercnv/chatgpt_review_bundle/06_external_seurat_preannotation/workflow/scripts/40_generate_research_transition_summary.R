options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
out_dir <- file.path(data_root, "research_validation_independent_cnv")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
report_path <- file.path(out_dir, "FINAL_EXTERNAL_SCRNA_RESEARCH_SUMMARY.md")
if (file.exists(report_path) && !replace_generated)
  stop("Output exists: ", report_path)

trace <- fread(file.path(out_dir, "missing_epithelial_count_trace.csv"))
coverage <- fread(file.path(out_dir, "GSE154600_final_count_coverage_v2.csv"))
stability <- fread(file.path(
  out_dir, "GSE154600_copykat_stability_by_patient_v2.csv"
))
bias <- fread(file.path(out_dir, "copykat_defined_bias_summary.csv"))
validation <- fread(file.path(
  out_dir, "GSE154600_infercnv_validation_by_patient.csv"
))
receiver <- fread(file.path(out_dir, "receiver_evidence_by_patient_final.csv"))
receiver_summary <- fread(file.path(
  out_dir, "receiver_reproducibility_summary.csv"
))
old <- fread(file.path(v61, "external_scrna_evidence_matrix_v3_1.csv"))

get_repro <- function(metric) {
  receiver_summary[receiver_metric == metric, patient_reproducibility][1]
}
gse_total <- validation[, .(
  infercnv_high_n = sum(n_infercnv_high),
  dual_method_support_n = sum(n_dual_method_support)
)]
st_total <- stability[, .(
  n_final = sum(n_final_epithelial),
  n_stable = sum(n_stable_aneuploid)
)]
coverage_total <- coverage[, sum(n_counts_available_final) / sum(n_final_epithelial)]
bias_status <- bias[patient_id == "ALL", copykat_defined_bias_status][1]

matrix <- old[, .(
  dataset_id,
  dataset_role,
  SPP1_cross_patient_reproducibility,
  final_epithelial_count_coverage = NA_real_,
  copykat_stable_aneuploid_n = NA_integer_,
  infercnv_high_n = NA_integer_,
  dual_method_support_n = NA_integer_,
  copykat_defined_bias_status = "NOT_APPLICABLE",
  CD44_patient_reproducibility = "NOT_APPLICABLE",
  ITGB1_patient_reproducibility = "NOT_APPLICABLE",
  ITGB1_alpha_patient_reproducibility = "NOT_APPLICABLE",
  dual_CD44_ITGB1_patient_reproducibility = "NOT_APPLICABLE",
  malignant_receiver_expression_context = "NOT_APPLICABLE",
  spatial_validation_priority = "UNCHANGED_NOT_PRIORITIZED",
  tumor_specificity_status,
  main_limitation
)]
i <- which(matrix$dataset_id == "GSE154600")
matrix[i, `:=`(
  final_epithelial_count_coverage = coverage_total,
  copykat_stable_aneuploid_n = st_total$n_stable,
  infercnv_high_n = gse_total$infercnv_high_n,
  dual_method_support_n = gse_total$dual_method_support_n,
  copykat_defined_bias_status = bias_status,
  CD44_patient_reproducibility = get_repro("CD44"),
  ITGB1_patient_reproducibility = get_repro("ITGB1"),
  ITGB1_alpha_patient_reproducibility = get_repro("ITGB1_alpha"),
  dual_CD44_ITGB1_patient_reproducibility =
    get_repro("dual_CD44_ITGB1"),
  malignant_receiver_expression_context = "DESCRIPTIVE_ONLY",
  spatial_validation_priority = "HIGH",
  main_limitation = paste(
    "Independent CNV used infercnvpy fallback because standard inferCNV R",
    "was unavailable; receiver expression lacks sample-level association,",
    "spatial proximity, direct ligand-receptor, and causal evidence"
  )
)]
fwrite(matrix, file.path(out_dir, "external_scrna_research_evidence_matrix.csv"),
       na = "NA")

fmt <- function(x, d = 3) {
  ifelse(is.na(x), "NA", format(round(x, d), trim = TRUE))
}
trace_by_patient <- trace[, .(
  n_traced = .N,
  n_recovered = sum(recoverable),
  n_found_raw = sum(found_in_raw_matrix),
  n_found_full_object = sum(found_in_full_object)
), by = patient_id]
trace_lines <- trace_by_patient[, paste0(
  "- ", patient_id, ": traced=", n_traced,
  "; recovered=", n_recovered,
  "; found in full object=", n_found_full_object,
  "; found in GEO feature-barcode matrix=", n_found_raw, "."
)]
coverage_lines <- coverage[, paste0(
  "- ", patient_id, ": ", n_counts_available_final, "/",
  n_final_epithelial, " (", fmt(final_count_coverage), "); recovered=",
  n_counts_recovered, "."
)]
infer_lines <- validation[, paste0(
  "- ", patient_id, ": high CNV=", n_infercnv_high,
  "; stable CopyKAT aneuploid=", n_copykat_stable_aneuploid,
  "; dual-method=", n_dual_method_support,
  "; CopyKAT-only=", n_copykat_only,
  "; inferCNV-only=", n_infercnv_only,
  "; concordance among stable CopyKAT=", fmt(copykat_infercnv_concordance), "."
)]
receiver_lines <- receiver[, paste0(
  "- ", patient_id, ": tier=", receiver_tier,
  "; n=", n_receiver_cells,
  "; CD44=", fmt(CD44_positive_fraction), " (", CD44_patient_status, ")",
  "; ITGB1=", fmt(ITGB1_positive_fraction), " (", ITGB1_patient_status, ")",
  "; ITGB1-alpha=", fmt(ITGB1_any_alpha_copositive_fraction),
  " (", ITGB1_alpha_patient_status, ")",
  "; CD44/ITGB1=", fmt(CD44_ITGB1_copositive_fraction),
  " (", dual_CD44_ITGB1_patient_status, ")."
)]
repro_lines <- receiver_summary[, paste0(
  "- ", receiver_metric, ": ", patient_reproducibility,
  "; evaluable patients=", n_evaluable_patients,
  "; supported=", n_supported_patients,
  "; low=", n_detected_low_patients,
  "; not detected=", n_not_detected_patients, "."
)]

lines <- c(
  "# Final external scRNA research transition summary",
  "",
  "## 1. Recovery of 511 final epithelial cells",
  "",
  trace_lines,
  "",
  paste0(
    "All ", nrow(trace), " previously unsubmitted final epithelial cells were ",
    "recovered from the authoritative GSE154600 preannotation Assay5 raw-count ",
    "layer. Their GEO feature-barcode source was also traced. The prior omission ",
    "was caused by incomplete lineage-strategy matrices, not by absent raw counts."
  ),
  "",
  "## 2. Final count coverage",
  "",
  coverage_lines,
  "",
  paste0(
    "Final GSE154600 epithelial count coverage is ",
    sum(coverage$n_counts_available_final), "/", sum(coverage$n_final_epithelial),
    " (", fmt(coverage_total), "). T76 is ",
    coverage[patient_id == "T76", n_counts_available_final], "/",
    coverage[patient_id == "T76", n_final_epithelial], " (",
    fmt(coverage[patient_id == "T76", final_count_coverage]), ")."
  ),
  "",
  "## 3. CopyKAT defined-cell technical bias",
  "",
  paste0(
    "Patient-stratified descriptive audit status: `", bias_status,
    "`. Cell-level tests are treated only as technical diagnostics, not as ",
    "patient-level biological replication."
  ),
  "",
  "## 4. Independent CNV support",
  "",
  infer_lines,
  "",
  paste0(
    "Across patients, inferCNV-high cells=", gse_total$infercnv_high_n,
    "; dual-method supportive cells=", gse_total$dual_method_support_n, "."
  ),
  "",
  "The standard inferCNV R package was unavailable because the local Bioconductor installation could not validate through the Windows security channel. The task-authorized infercnvpy fallback was therefore used for continuous CNV signal validation. No HMM subclone calls and no `confirmed malignant` wording are used.",
  "",
  "## 5. Patient-replicated receiver expression",
  "",
  receiver_lines,
  "",
  repro_lines,
  "",
  "`malignant_receiver_expression_context` is `DESCRIPTIVE_ONLY`: there is no sample-level association statistic, spatial proximity evidence, direct ligand-receptor evidence, or causal evidence.",
  "",
  "## 6. Research transition",
  "",
  "GSE154600 has `spatial_validation_priority = HIGH`. The external scRNA evidence is sufficient to prioritize a spatial validation study of SPP1 macrophages and CNV-supported epithelial ITGB1/CD44 context, but it does not itself establish spatial interaction, direct receptor engagement, tumor specificity, or causality.",
  "",
  "External scRNA preprocessing stops here: no v6.2/v6.3 threshold versions, additional CopyKAT seeds, SPP1 threshold changes, or new cell-type cleanup are introduced."
)
writeLines(lines, report_path)
message("Final external scRNA research transition summary complete")
