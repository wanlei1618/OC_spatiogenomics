options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v6 <- file.path(data_root, "diagnostics_v6_malignant_receiver_validation")
out <- file.path(v6, "FINAL_EXTERNAL_SCRNA_FREEZE_REPORT.md")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

matrix <- fread(file.path(v6, "external_scrna_evidence_matrix_v3.csv"))
sample <- fread(file.path(v6, "state_detection_by_sample.csv"))
sens <- unique(fread(file.path(v6, "state_threshold_sensitivity.csv"))[
  minimum_macrophages == 20 &
    transcript_positive_fraction == .10 & program_cell_fraction == .10,
  .(dataset_id, main_conclusion, SPP1_threshold_robustness)
])
cyt <- fread(file.path(v6, "GSE147082_cluster6_tcr_nk_evidence.csv"))
cnv <- fread(file.path(v6, "malignancy_summary_by_patient.csv"))
pt <- fread(file.path(v6, "GSE147082_PT2834_formal_cnv_cluster_summary.csv"))
recv <- fread(file.path(v6, "malignant_epithelial_receiver_context.csv"))

fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", format(round(x, digits), trim = TRUE))
matrix_lines <- apply(matrix, 1, function(x) paste0(
  "- ", x[["dataset_id"]], ": SPP1 transcript=", x[["SPP1_transcript_detection_status"]],
  "; program=", x[["SPP1_program_support_status"]],
  "; reproducibility=", x[["SPP1_cross_patient_reproducibility"]],
  "; malignancy=", x[["malignant_epithelial_available"]],
  "; receiver context=", x[["sender_receiver_context_status"]], "."
))
p04 <- sample[dataset_id == "GSE158722" & sample_id == "P04_Time3"][1]
sub2 <- cyt[subcluster == "2"][1]
g154 <- cnv[dataset_id == "GSE154600", .(
  supportive = sum(n_malignant_supportive),
  high_confidence = sum(n_malignant_high_confidence),
  evaluated = sum(n_copykat_evaluated)
)]
g147 <- cnv[dataset_id == "GSE147082", .(
  supportive = sum(n_malignant_supportive),
  high_confidence = sum(n_malignant_high_confidence),
  evaluated = sum(n_copykat_evaluated)
)]
r154 <- recv[
  dataset_id == "GSE154600" & receiver_tier == "malignant_supportive" &
    n_receiver_cells >= 20,
  .(
    n_samples = .N,
    median_CD44_positive = median(CD44_positive_fraction, na.rm = TRUE),
    median_ITGB1_positive = median(ITGB1_positive_fraction, na.rm = TRUE),
    median_dual_positive = median(CD44_ITGB1_copositive_fraction, na.rm = TRUE),
    median_ITGB1_alpha_copositive =
      median(ITGB1_alpha_copositive_fraction, na.rm = TRUE)
  )
]

lines <- c(
  "# Final external scRNA evidence freeze report",
  "",
  "## Frozen scope",
  "",
  "This report freezes the external scRNA preprocessing and evidence definitions after targeted count-based state detection, GSE147082 cluster 6 correction, formal CopyKAT assessment, and malignancy-stratified receiver analysis. No dataset-wide QC, clustering, or broad annotation was rerun.",
  "",
  "## Dataset-level evidence",
  "",
  matrix_lines,
  "",
  "SPP1 transcript detection, companion-program support, and within-dataset relative enrichment are separate endpoints. The top-quartile SPP1-high fraction is retained only as a relative indicator and never defines transcript presence.",
  "",
  "## GSE158722 P04_Time3",
  "",
  paste0(
    "P04_Time3 contains ", p04$n_macrophages, " high-confidence macrophages; ",
    fmt(100 * p04$SPP1_positive_fraction, 1), "% detect SPP1 transcript and ",
    fmt(100 * p04$SPP1_program_cell_fraction, 1),
    "% meet the SPP1-plus-companion program definition. Its final status is `",
    p04$SPP1_joint_status, "`, not SPP1-negative."
  ),
  "",
  "## Threshold sensitivity",
  "",
  paste0(
    "- ", sens$dataset_id, ": ", sens$main_conclusion,
    "; robustness=", sens$SPP1_threshold_robustness, "."
  ),
  "",
  "The formal analysis uses at least 20 macrophages, transcript-positive fraction at least 0.10, and program-cell fraction at least 0.10. Sensitivity spans 10/20/30 cells and 0.05/0.10/0.20 transcript and program thresholds.",
  "",
  "## GSE147082 cluster 6",
  "",
  paste0(
    "Cycling is stored only in `cell_state`. Subcluster 2 has CD3/TCR co-positive fraction ",
    fmt(sub2$CD3_TCR_copositive_fraction), " and NK-marker-positive fraction ",
    fmt(sub2$NCR1_NCAM1_FCER1G_KLRD1_positive_fraction),
    "; it is frozen as `", sub2$final_cell_type,
    "` with `cell_state=NK_like_cytotoxic`, `patient_enriched=TRUE`, and confidence `Review`."
  ),
  "",
  "## Formal malignancy assessment",
  "",
  paste0(
    "No reusable target-specific inferCNV or CopyKAT result was found. New patient-internal CopyKAT runs evaluated ",
    g154$evaluated, " GSE154600 epithelial candidates and ", g147$evaluated,
    " GSE147082 PT-2834 candidates. GSE154600 has ", g154$supportive,
    " single-method malignant-supportive cells and ", g154$high_confidence,
    " two-method high-confidence cells. GSE147082 has ", g147$supportive,
    " malignant-supportive and ", g147$high_confidence,
    " high-confidence cells."
  ),
  "",
  paste0(
    "PT-2834 cluster 4: `", pt[target_group == "cluster_4", final_label],
    "` (CopyKAT aneuploid fraction ",
    fmt(pt[target_group == "cluster_4", aneuploid_fraction]), ")."
  ),
  paste0(
    "PT-2834 cluster 7: `", pt[target_group == "cluster_7", final_label],
    "` (CopyKAT aneuploid fraction ",
    fmt(pt[target_group == "cluster_7", aneuploid_fraction]), ")."
  ),
  "",
  "The previous CNV-like intensity ratio is audit-only and has no role in malignancy classification. CopyKAT diploid calls are not promoted to `DIPLOID_SUPPORTIVE` because inferCNV is unavailable. Likewise, CopyKAT aneuploid calls remain `MALIGNANT_SUPPORTIVE`, not two-method `MALIGNANT_HIGH_CONFIDENCE`.",
  "",
  "GSE158722 malignancy is `NOT_EVALUABLE`: platform identity cannot be reliably recovered, so a same-platform reference cannot be selected without fabrication.",
  "",
  "## Malignant receiver context",
  "",
  paste0(
    "Among evaluable GSE154600 malignant-supportive samples (n=", r154$n_samples,
    "), median positive fractions are CD44=", fmt(r154$median_CD44_positive),
    ", ITGB1=", fmt(r154$median_ITGB1_positive),
    ", CD44/ITGB1 dual=", fmt(r154$median_dual_positive),
    ", and ITGB1/alpha-integrin co-positive=", fmt(r154$median_ITGB1_alpha_copositive), "."
  ),
  "",
  "CD44, ITGB1, ITGB1-alpha partner, and dual support are reported separately. ITGB1 expression alone is not evidence of a complete functional integrin receptor, and expression co-occurrence does not establish direct SPP1 binding.",
  "",
  "## Interpretation boundary",
  "",
  "- SPP1 transcript and companion-program support recur across multiple patients, with patient heterogeneity.",
  "- GSE151214 and GSE154763 remain reference-only and are excluded from tumor effect aggregation.",
  "- `tumor_specificity_status` is `NOT_ESTABLISHED` for every dataset.",
  "- The data do not establish tumor-specific SPP1 macrophages, direct SPP1-ITGB1 binding, receptor activation, causality, or spatial contact.",
  "- Spatial assays and/or wet-lab perturbation are required for mechanism validation.",
  "",
  "## Frozen outputs",
  "",
  "- `state_detection_by_sample.csv` and `state_detection_by_patient.csv`",
  "- `state_threshold_sensitivity.csv`",
  "- `GSE147082_cluster6_final_cell_annotation.csv` and `GSE147082_cluster6_tcr_nk_evidence.csv`",
  "- `existing_malignancy_results_audit.csv`, `malignancy_summary_by_patient.csv`, and local per-cell consensus",
  "- `malignant_epithelial_receiver_context.csv`",
  "- `external_scrna_evidence_matrix_v3.csv`",
  "",
  "External scRNA preprocessing and evidence definitions are frozen at v6; subsequent work should proceed to spatial and mechanistic validation rather than another ordinary threshold revision."
)
writeLines(lines, out)
message("Final external scRNA freeze report complete")
