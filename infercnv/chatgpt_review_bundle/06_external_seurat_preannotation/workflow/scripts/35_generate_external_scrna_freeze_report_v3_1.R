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
out <- file.path(v61, "FINAL_EXTERNAL_SCRNA_FREEZE_REPORT_V3_1.md")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

cov <- fread(file.path(v61, "GSE154600_copykat_target_coverage_audit.csv"))
st <- fread(file.path(v61, "GSE154600_copykat_stability_by_patient.csv"))
runs <- fread(file.path(v61, "copykat_stability_run_status.csv"))
rec <- fread(file.path(v61, "stable_malignant_receiver_context_by_patient.csv"))
rec_sum <- fread(file.path(v61, "stable_malignant_receiver_context_summary.csv"))[
  receiver_tier == "stable_malignant_supportive"
][1]
ref <- fread(file.path(v61, "GSE154763_reference_terminology_audit.csv"))
pt <- fread(file.path(v6, "GSE147082_PT2834_formal_cnv_cluster_summary.csv"))
total <- st[, lapply(.SD, sum, na.rm = TRUE), .SDcols = patterns("^n_")]
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", format(round(x, d), trim = TRUE))

patient_lines <- merge(cov, st, by = c("dataset_id", "patient_id"))
patient_lines <- patient_lines[, paste0(
  "- ", patient_id, ": final=", n_final_epithelial.x,
  "; available/submitted=", n_available_for_copykat, "/", n_submitted_any_run,
  " (coverage ", fmt(target_input_coverage), ")",
  "; any defined=", n_defined_in_at_least_1_run,
  " (rate ", fmt(any_defined_call_rate), ")",
  "; stable aneuploid=", n_stable_aneuploid,
  "; stable diploid=", n_stable_diploid,
  "; unstable=", n_unstable,
  "; mostly not defined=", n_mostly_not_defined,
  "; not submitted=", n_not_submitted, "."
)]
seed_ranges <- runs[, .(
  aneuploid_min = min(n_aneuploid), aneuploid_max = max(n_aneuploid),
  diploid_min = min(n_diploid), diploid_max = max(n_diploid),
  not_defined_min = min(n_not_defined), not_defined_max = max(n_not_defined)
), by = patient_id]
seed_lines <- seed_ranges[, paste0(
  "- ", patient_id, ": aneuploid ", aneuploid_min, "-", aneuploid_max,
  "; diploid ", diploid_min, "-", diploid_max,
  "; not defined ", not_defined_min, "-", not_defined_max, "."
)]
receiver_lines <- rec[receiver_tier == "stable_malignant_supportive", paste0(
  "- ", patient_id, ": n=", n_receiver_cells,
  "; CD44=", fmt(CD44_positive_fraction),
  " (", CD44_receiver_support, ")",
  "; ITGB1=", fmt(ITGB1_positive_fraction),
  " (", ITGB1_receiver_support, ")",
  "; ITGB1/any-alpha=", fmt(ITGB1_any_alpha_copositive_fraction),
  " (", ITGB1_alpha_partner_support, ")",
  "; dual=", fmt(CD44_ITGB1_copositive_fraction),
  " (", dual_CD44_ITGB1_support, ")",
  "; dominant alpha=", dominant_alpha_partner, "."
)]
lines <- c(
  "# Final external scRNA freeze report v3.1",
  "",
  "## Scope",
  "",
  "v6.1 audits GSE154600 CopyKAT target coverage, repeats same-patient immune-reference sampling with seeds 20260718/20260719/20260720, builds a stable single-method malignancy layer, recalculates receiver expression, and corrects GSE154763 reference terminology. v6 outputs remain unchanged; no full QC, clustering, broad annotation, or SPP1 threshold was rerun.",
  "",
  "## CopyKAT coverage and stability",
  "",
  patient_lines,
  "",
  paste0(
    "Across GSE154600, ", total$n_final_epithelial, " final epithelial cells were audited; ",
    total$n_available_for_copykat, " were found in at least one of the four specified lineage count inputs and submitted; ",
    total$n_defined_in_at_least_1_run, " received a defined call at least once; ",
    total$n_stable_aneuploid, " were stable aneuploid; ",
    total$n_stable_diploid, " were stable diploid-like; ",
    total$n_unstable, " were discordant; ",
    total$n_mostly_not_defined, " were mostly not defined; and ",
    total$n_not_submitted, " were not submitted because counts were absent."
  ),
  "",
  "## T76 coverage correction",
  "",
  paste0(
    "v6 submitted 37 of 125 T76 final epithelial cells. Cross-lineage collection still finds 37 of 125 (coverage ",
    fmt(st[patient_id == "T76", target_input_coverage]),
    "); the remaining 88 are absent from all four task-specified lineage count matrices. Coverage therefore did not numerically increase, but the omission is now explicit and classified as `NOT_SUBMITTED`, rather than silently excluded."
  ),
  "",
  "## Three-seed reference stability",
  "",
  seed_lines,
  "",
  paste0(
    "Reference resampling produced ", total$n_stable_aneuploid,
    " stable aneuploid and ", total$n_stable_diploid,
    " stable diploid-like calls, with ", total$n_unstable,
    " discordant cells. Seed-to-seed totals varied modestly but did not overturn the dominant patient-level pattern."
  ),
  "",
  "Stable CopyKAT aneuploid is `MALIGNANT_SUPPORTIVE_STABLE`, which remains single-method supportive evidence. It is not double-method high-confidence malignancy because inferCNV is unavailable.",
  "",
  "## Stable malignant receiver",
  "",
  receiver_lines,
  "",
  paste0(
    "Across evaluable patients, stable malignant-supportive receiver status is CD44=",
    rec_sum$CD44_receiver_status, ", ITGB1=", rec_sum$ITGB1_receiver_status,
    ", ITGB1-alpha partner=", rec_sum$ITGB1_alpha_partner_status,
    ", and CD44/ITGB1 dual=", rec_sum$dual_receiver_status, "."
  ),
  "",
  "ITGB1 expression alone does not establish a complete functional receptor, and co-expression does not establish direct SPP1 binding.",
  "",
  "## GSE147082 retained interpretation",
  "",
  paste0("- Cluster 4 remains `", pt[target_group == "cluster_4", final_label], "`."),
  paste0("- Cluster 7 remains `", pt[target_group == "cluster_7", final_label], "`."),
  "",
  "The prior single CopyKAT diploid result is retained as single-method evidence and does not over-upgrade or force a malignant/normal label.",
  "",
  "## GSE154763 terminology correction",
  "",
  paste0(
    "- ", ref$state, ": ", ref$terminology_decision,
    "; expression columns available=", ref$expression_columns_available, "."
  ),
  "",
  "SPP1, C1QC-core, and FOLR2 transcript detection now use their explicit normalized expression columns (>0). Reference program support is reported separately. Threshold robustness is termed `reference_consistent`, `reference_threshold_sensitive`, or `reference_not_evaluable`, rather than raw-count `stable_replicated`.",
  "",
  "## Frozen evidence",
  "",
  "- SPP1/C1QC/FOLR2 v6 raw-count thresholds and conclusions remain frozen.",
  "- GSE154600 CopyKAT coverage, defined-call rate, three-seed stability, and stable receiver context are frozen at v3.1.",
  "- GSE151214 remains normal-reference only; GSE154763 remains author-normalized myeloid-reference only.",
  "- `tumor_specificity_status` remains `NOT_ESTABLISHED` for every dataset.",
  "",
  "## Remaining validation",
  "",
  "inferCNV or another independent formal CNV method is required for double-method high-confidence malignancy. Spatial assays are required to establish sender-receiver proximity. Wet-lab assays are required to establish direct binding, receptor activation, and causality. The available data do not establish tumor-specific SPP1 macrophages or SPP1-driven progression.",
  "",
  "External scRNA evidence is frozen at v6.1; subsequent work should move to independent CNV, spatial, and mechanistic validation rather than another ordinary threshold revision."
)
writeLines(lines, out)
message("Final external scRNA freeze report v3.1 complete")
