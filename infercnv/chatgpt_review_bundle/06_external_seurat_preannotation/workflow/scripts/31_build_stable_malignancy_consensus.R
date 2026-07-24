options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v61, "GSE154600_copykat_stability_by_cell.csv.gz")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

manifest <- fread(file.path(v61, "GSE154600_copykat_target_manifest.csv.gz"))
pred <- fread(file.path(v61, "copykat_stability_raw_predictions.csv.gz"))
runs <- fread(file.path(v61, "copykat_stability_run_status.csv"))
completed <- runs[run_status == "COMPLETED", .(n_runs_attempted = .N), by = patient_id]

calls <- pred[, .(
  n_runs_defined = sum(copykat_status %in% c("aneuploid", "diploid")),
  n_aneuploid_calls = sum(copykat_status == "aneuploid"),
  n_diploid_calls = sum(copykat_status == "diploid")
), by = .(patient_id, cell_id)]
cell <- merge(
  manifest[, .(
    dataset_id, patient_id, cell_id, available_for_copykat,
    source_lineage_input, n_source_matches
  )],
  completed, by = "patient_id", all.x = TRUE
)
cell <- merge(cell, calls, by = c("patient_id", "cell_id"), all.x = TRUE)
for (col in c("n_runs_attempted", "n_runs_defined", "n_aneuploid_calls", "n_diploid_calls"))
  set(cell, which(is.na(cell[[col]])), col, 0L)
cell[, n_not_defined := pmax(n_runs_attempted - n_runs_defined, 0L)]
cell[, aneuploid_fraction_among_defined := fifelse(
  n_runs_defined > 0, n_aneuploid_calls / n_runs_defined, NA_real_
)]
cell[, diploid_fraction_among_defined := fifelse(
  n_runs_defined > 0, n_diploid_calls / n_runs_defined, NA_real_
)]
cell[, defined_call_rate := fifelse(
  n_runs_attempted > 0, n_runs_defined / n_runs_attempted, NA_real_
)]
cell[, stability_class := fcase(
  !available_for_copykat | n_runs_attempted == 0, "NOT_SUBMITTED",
  n_runs_defined >= 2 & aneuploid_fraction_among_defined >= .80, "STABLE_ANEUPLOID",
  n_runs_defined >= 2 & diploid_fraction_among_defined >= .80, "STABLE_DIPLOID",
  n_runs_defined >= 2, "UNSTABLE_DISCORDANT",
  default = "MOSTLY_NOT_DEFINED"
)]
cell[, malignancy_evidence := fcase(
  stability_class == "STABLE_ANEUPLOID", "MALIGNANT_SUPPORTIVE_STABLE",
  stability_class == "STABLE_DIPLOID", "DIPLOID_LIKE_STABLE",
  stability_class == "UNSTABLE_DISCORDANT", "DISCORDANT_SINGLE_METHOD",
  default = "NOT_EVALUABLE"
)]
fwrite(cell, out, compress = "gzip", na = "NA")

by_patient <- cell[, .(
  n_final_epithelial = .N,
  n_available_for_copykat = sum(available_for_copykat),
  n_submitted_any_run = sum(available_for_copykat & n_runs_attempted > 0),
  n_defined_in_at_least_1_run = sum(n_runs_defined >= 1),
  n_defined_in_at_least_2_runs = sum(n_runs_defined >= 2),
  n_stable_aneuploid = sum(stability_class == "STABLE_ANEUPLOID"),
  n_stable_diploid = sum(stability_class == "STABLE_DIPLOID"),
  n_unstable = sum(stability_class == "UNSTABLE_DISCORDANT"),
  n_mostly_not_defined = sum(stability_class == "MOSTLY_NOT_DEFINED"),
  n_not_submitted = sum(stability_class == "NOT_SUBMITTED")
), by = .(dataset_id, patient_id)]
by_patient[, `:=`(
  target_input_coverage = n_available_for_copykat / n_final_epithelial,
  any_defined_call_rate = n_defined_in_at_least_1_run / n_submitted_any_run,
  stable_defined_call_rate = n_defined_in_at_least_2_runs / n_submitted_any_run,
  stable_aneuploid_fraction_all_final = n_stable_aneuploid / n_final_epithelial,
  stable_aneuploid_fraction_defined = fifelse(
    n_defined_in_at_least_2_runs > 0,
    n_stable_aneuploid / n_defined_in_at_least_2_runs, NA_real_
  )
)]
fwrite(by_patient, file.path(v61, "GSE154600_copykat_stability_by_patient.csv"), na = "NA")

summary <- by_patient[, lapply(.SD, sum, na.rm = TRUE), .SDcols = patterns("^n_")]
summary[, `:=`(
  dataset_id = "GSE154600",
  target_input_coverage = n_available_for_copykat / n_final_epithelial,
  any_defined_call_rate = n_defined_in_at_least_1_run / n_submitted_any_run,
  stable_defined_call_rate = n_defined_in_at_least_2_runs / n_submitted_any_run,
  stable_aneuploid_fraction_all_final = n_stable_aneuploid / n_final_epithelial,
  stable_aneuploid_fraction_defined = n_stable_aneuploid / n_defined_in_at_least_2_runs
)]
setcolorder(summary, "dataset_id")
fwrite(summary, file.path(v61, "GSE154600_copykat_stability_summary.csv"), na = "NA")

p1dt <- melt(
  by_patient,
  id.vars = "patient_id",
  measure.vars = c("target_input_coverage", "any_defined_call_rate", "stable_defined_call_rate"),
  variable.name = "metric", value.name = "fraction"
)
p1 <- ggplot(p1dt, aes(patient_id, fraction, fill = metric)) +
  geom_col(position = "dodge") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  labs(title = "GSE154600 CopyKAT target coverage and call rates", x = NULL, y = "fraction")
ggsave(file.path(v61, "01_copykat_coverage_and_call_rate.png"), p1,
       width = 9, height = 5, dpi = 180)

consistency <- pred[, .(
  aneuploid_fraction = mean(copykat_status == "aneuploid"),
  diploid_fraction = mean(copykat_status == "diploid"),
  not_defined_fraction = mean(!copykat_status %in% c("aneuploid", "diploid"))
), by = .(patient_id, seed)]
p2dt <- melt(consistency, id.vars = c("patient_id", "seed"),
             variable.name = "call", value.name = "fraction")
p2 <- ggplot(p2dt, aes(factor(seed), fraction, fill = call)) +
  geom_col() +
  facet_wrap(~patient_id) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "CopyKAT call composition across three reference seeds",
       x = "seed", y = "submitted target fraction")
ggsave(file.path(v61, "02_copykat_three_seed_consistency.png"), p2,
       width = 11, height = 5.5, dpi = 180)

p3dt <- melt(
  by_patient,
  id.vars = "patient_id",
  measure.vars = c(
    "n_stable_aneuploid", "n_stable_diploid", "n_unstable",
    "n_mostly_not_defined", "n_not_submitted"
  ),
  variable.name = "stability", value.name = "n_cells"
)
p3 <- ggplot(p3dt, aes(patient_id, n_cells, fill = stability)) +
  geom_col(position = "fill") +
  theme_bw() +
  labs(title = "Stable CopyKAT classification among final epithelial cells",
       x = NULL, y = "fraction")
ggsave(file.path(v61, "03_copykat_stability_composition.png"), p3,
       width = 9, height = 5, dpi = 180)
message("Stable single-method malignancy consensus complete")
