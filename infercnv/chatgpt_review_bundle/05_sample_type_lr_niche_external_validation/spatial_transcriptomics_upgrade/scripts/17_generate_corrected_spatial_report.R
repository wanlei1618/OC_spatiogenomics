options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

local_root <- "D:/OC_spatiogenomics/spatial_data/spatial_analysis_correction_v2"
audit <- fread(file.path(local_root, "spatial_sample_audit.csv"))
neighborhood <- fread(file.path(
  local_root, "GSE203612_corrected_neighborhood_results.csv"
))
expr <- fread(file.path(
  local_root, "GSE189843_expression_summary_by_sample.csv"
))
expr_cmp <- fread(file.path(
  local_root, "GSE189843_response_group_comparison.csv"
))
gse211 <- fread(file.path(local_root, "GSE211956_status.csv"))
signature <- fread(file.path(local_root, "updated_spatial_signatures.csv"))
evidence_tier <- unique(signature[
  signature_name == "CNV_supported_malignant_epithelial_signature",
  evidence_tier
])[[1L]]

primary <- neighborhood[
  threshold == "top25" & neighbor_definition == "hex6"
]
get_result <- function(sample_value, comparison_value) {
  primary[
    primary$sample_id == sample_value &
      primary$comparison == comparison_value
  ][1L]
}
axis_status <- function(receiver) {
  main_name <- paste0("SPP1_sender_to_", receiver, "_receiver")
  c1_name <- paste0("C1QC_sender_to_", receiver, "_receiver")
  rows <- primary[comparison == main_name]
  specific <- vapply(rows$sample_id, function(sid) {
    main <- get_result(sid, main_name)
    c1 <- get_result(sid, c1_name)
    general <- get_result(sid, "all_macrophage_to_all_epithelial")
    isTRUE(main$observed_expected_ratio > 1) &&
      isTRUE(main$empirical_p < .05) &&
      isTRUE(main$observed_expected_ratio > c1$observed_expected_ratio) &&
      isTRUE(main$observed_expected_ratio > general$observed_expected_ratio)
  }, logical(1))
  if (nrow(rows) < 2L) return("NOT_EVALUABLE")
  if (all(specific)) return("REPLICATED_SPATIAL_SUPPORT")
  if (any(rows$observed_expected_ratio > 1) &&
      any(rows$observed_expected_ratio <= 1)) {
    return("SPATIALLY_HETEROGENEOUS")
  }
  if (sum(specific) == 1L) return("LIMITED_SPATIAL_SUPPORT")
  "NO_SPECIFIC_SPATIAL_SUPPORT"
}
itgb1_status <- axis_status("ITGB1")
cd44_status <- axis_status("CD44")

matrix_203 <- rbindlist(lapply(c("GSM6177614", "GSM6177617"), function(sid) {
  it <- get_result(sid, "SPP1_sender_to_ITGB1_receiver")
  cd <- get_result(sid, "SPP1_sender_to_CD44_receiver")
  c1it <- get_result(sid, "C1QC_sender_to_ITGB1_receiver")
  c1cd <- get_result(sid, "C1QC_sender_to_CD44_receiver")
  general <- get_result(sid, "all_macrophage_to_all_epithelial")
  data.table(
    dataset_id = "GSE203612", sample_id = sid,
    tumor_type_verified = "ovarian carcinoma",
    analysis_level = "COORDINATE_AWARE_VISIUM",
    coordinates_verified = TRUE,
    SPP1_ITGB1_effect_direction = fifelse(
      it$observed_expected_ratio > 1, "ENRICHED", "NOT_ENRICHED"
    ),
    SPP1_ITGB1_enrichment_ratio = it$observed_expected_ratio,
    SPP1_ITGB1_empirical_p = it$empirical_p,
    SPP1_CD44_effect_direction = fifelse(
      cd$observed_expected_ratio > 1, "ENRICHED", "NOT_ENRICHED"
    ),
    SPP1_CD44_enrichment_ratio = cd$observed_expected_ratio,
    SPP1_CD44_empirical_p = cd$empirical_p,
    C1QC_ITGB1_control_result = sprintf(
      "OE=%.4f;p=%.4f", c1it$observed_expected_ratio, c1it$empirical_p
    ),
    C1QC_CD44_control_result = sprintf(
      "OE=%.4f;p=%.4f", c1cd$observed_expected_ratio, c1cd$empirical_p
    ),
    general_macrophage_epithelial_control_result = sprintf(
      "OE=%.4f;p=%.4f", general$observed_expected_ratio,
      general$empirical_p
    ),
    ITGB1_dataset_conclusion = itgb1_status,
    CD44_dataset_conclusion = cd44_status,
    malignant_signature_evidence_tier = evidence_tier,
    statistical_unit = "independent ovarian section",
    main_limitation = "Visium spots are mixtures; n=2 sections"
  )
}))
matrix_189 <- expr[, .(
  dataset_id = "GSE189843", sample_id,
  tumor_type_verified = "high-grade serous ovarian carcinoma",
  analysis_level = "EXPRESSION_ONLY",
  coordinates_verified = FALSE,
  SPP1_ITGB1_effect_direction = fifelse(
    SPP1_ITGB1_sample_internal_spearman > 0, "POSITIVE_EXPRESSION_ASSOCIATION",
    "NONPOSITIVE_EXPRESSION_ASSOCIATION"
  ),
  SPP1_ITGB1_enrichment_ratio = NA_real_,
  SPP1_ITGB1_empirical_p = NA_real_,
  SPP1_CD44_effect_direction = fifelse(
    SPP1_CD44_sample_internal_spearman > 0, "POSITIVE_EXPRESSION_ASSOCIATION",
    "NONPOSITIVE_EXPRESSION_ASSOCIATION"
  ),
  SPP1_CD44_enrichment_ratio = NA_real_,
  SPP1_CD44_empirical_p = NA_real_,
  C1QC_ITGB1_control_result = "EXPRESSION_SUMMARY_ONLY",
  C1QC_CD44_control_result = "EXPRESSION_SUMMARY_ONLY",
  general_macrophage_epithelial_control_result = "NOT_EVALUABLE_WITHOUT_COORDINATES",
  ITGB1_dataset_conclusion = "EXPRESSION_LEVEL_SUPPORT_ONLY",
  CD44_dataset_conclusion = "EXPRESSION_LEVEL_SUPPORT_ONLY",
  malignant_signature_evidence_tier = evidence_tier,
  statistical_unit = "sample/patient",
  main_limitation = "coordinates not authoritatively verified"
)]
matrix_211 <- data.table(
  dataset_id = "GSE211956", sample_id = NA_character_,
  tumor_type_verified = "NOT_VERIFIED_IN_FROZEN_REGISTRY",
  analysis_level = gse211$status,
  coordinates_verified = FALSE,
  SPP1_ITGB1_effect_direction = "NOT_EVALUABLE",
  SPP1_ITGB1_enrichment_ratio = NA_real_,
  SPP1_ITGB1_empirical_p = NA_real_,
  SPP1_CD44_effect_direction = "NOT_EVALUABLE",
  SPP1_CD44_enrichment_ratio = NA_real_,
  SPP1_CD44_empirical_p = NA_real_,
  C1QC_ITGB1_control_result = "NOT_EVALUABLE",
  C1QC_CD44_control_result = "NOT_EVALUABLE",
  general_macrophage_epithelial_control_result = "NOT_EVALUABLE",
  ITGB1_dataset_conclusion = "NOT_EVALUABLE",
  CD44_dataset_conclusion = "NOT_EVALUABLE",
  malignant_signature_evidence_tier = evidence_tier,
  statistical_unit = "not analyzed",
  main_limitation = gse211$reason
)
evidence <- rbindlist(list(matrix_203, matrix_189, matrix_211), fill = TRUE)
fwrite(evidence, file.path(local_root, "corrected_spatial_evidence_matrix.csv"),
       na = "NA")

fmt_axis <- function(axis) {
  test <- paste0("SPP1_sender_to_", axis, "_receiver")
  x <- primary[comparison == test]
  paste(sprintf(
    "%s OE=%.3f, empirical p=%.4f", x$sample_id,
    x$observed_expected_ratio, x$empirical_p
  ), collapse = "; ")
}
fmt_control <- function(test) {
  x <- primary[comparison == test]
  paste(sprintf(
    "%s OE=%.3f, p=%.4f", x$sample_id,
    x$observed_expected_ratio, x$empirical_p
  ), collapse = "; ")
}
comparison_focus <- expr_cmp[metric %chin% c(
  "median_SPP1_program", "median_ITGB1", "median_CD44",
  "SPP1_ITGB1_sample_internal_spearman",
  "SPP1_CD44_sample_internal_spearman"
)]
expression_lines <- paste(sprintf(
  "%s: poor-minus-excellent median %.4f, Cliff delta %.3f, exact p %.4f",
  comparison_focus$metric,
  comparison_focus$median_difference_poor_minus_excellent,
  comparison_focus$cliffs_delta_poor_vs_excellent,
  comparison_focus$exact_wilcoxon_p
), collapse = "\n- ")
priority <- if (itgb1_status %chin% c(
  "REPLICATED_SPATIAL_SUPPORT", "LIMITED_SPATIAL_SUPPORT"
)) "ITGB1" else if (cd44_status %chin% c(
  "REPLICATED_SPATIAL_SUPPORT", "LIMITED_SPATIAL_SUPPORT"
)) "CD44" else "neither until further validation"

report <- c(
  "# Corrected spatial validation report",
  "",
  "## 1. Sample audit and GSM6177618 exclusion",
  "",
  paste0(
    "Only GSM6177614 and GSM6177617 are coordinate-aware ovarian sections. ",
    "GSM6177618 is verified PDAC and is excluded from coordinate analysis, ",
    "expression analysis, pooled ovarian conclusions, figures and evidence ",
    "matrices (`PDAC_NOT_OVARIAN`)."
  ),
  "",
  "## 2. GSE203612 SPP1 to ITGB1",
  "",
  paste0(fmt_axis("ITGB1"), ". Dataset conclusion: **", itgb1_status, "**."),
  "",
  "## 3. GSE203612 SPP1 to CD44",
  "",
  paste0(fmt_axis("CD44"), ". Dataset conclusion: **", cd44_status, "**."),
  "",
  "## 4. C1QC and general macrophage controls",
  "",
  paste0("C1QC to ITGB1: ", fmt_control(
    "C1QC_sender_to_ITGB1_receiver"
  )),
  "",
  paste0("C1QC to CD44: ", fmt_control(
    "C1QC_sender_to_CD44_receiver"
  )),
  "",
  paste0("General macrophage to epithelial: ", fmt_control(
    "all_macrophage_to_all_epithelial"
  )),
  "",
  "## 5. GSE189843 expression-only evidence",
  "",
  paste0(
    nrow(expr), " author-included samples with available matrices were ",
    "retained (6 Excellent, 6 Poor). Coordinates remain unverified; these ",
    "results are sample/patient-level expression summaries and must not be ",
    "called spatial colocalization, proximity or neighborhood enrichment."
  ),
  "",
  paste0("- ", expression_lines),
  "",
  "## 6. Optional GSE211956 replication",
  "",
  paste0(
    "**", gse211$status, "**: ", gse211$reason,
    ". No new download or coordinate analysis was performed."
  ),
  "",
  "## 7. Wet-lab priority",
  "",
  paste0(
    "Priority: **", priority, "**, based on the corrected two-section ",
    "specificity comparison rather than expression-product maps."
  ),
  "",
  "## 8. Interpretation boundary",
  "",
  paste0(
    "The analysis can describe enrichment, depletion or heterogeneity of ",
    "SPP1 macrophage programs near ITGB1/CD44-positive CNV-supported ",
    "epithelial expression regions. It does not prove direct SPP1 binding, ",
    "receptor activation, causal tumor progression or tumor-specificity of ",
    "SPP1 macrophages. No output is represented as COMMOT or optimal ",
    "transport."
  )
)
writeLines(report,
           file.path(local_root, "CORRECTED_SPATIAL_VALIDATION_REPORT.md"),
           useBytes = TRUE)
message("Corrected spatial evidence matrix and report generated")
