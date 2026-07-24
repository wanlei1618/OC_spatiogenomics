#!/usr/bin/env Rscript

required <- c("yaml", "data.table", "jsonlite", "digest")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else "."
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(cfg$project$repo_root, winslash = "/", mustWork = TRUE)
local_root <- file.path(data_root, "diagnostics_v2")
repo_package <- file.path(repo_root, "infercnv", "chatgpt_review_bundle",
                          "06_external_seurat_preannotation", "diagnostics_v2")
dir.create(repo_package, recursive = TRUE, showWarnings = FALSE)

read_if <- function(path) if (file.exists(path)) data.table::fread(path, showProgress = FALSE) else data.table::data.table()
mt147 <- read_if(file.path(local_root, "GSE147082", "01_mt_audit", "mt_feature_detection.csv"))
mt158 <- read_if(file.path(local_root, "GSE158722", "01_mt_audit", "mt_feature_detection.csv"))
ret147 <- read_if(file.path(local_root, "GSE147082", "01_mt_audit", "qc_cell_retention_repaired.csv"))
ret158 <- read_if(file.path(local_root, "GSE158722", "01_mt_audit", "qc_cell_retention_repaired.csv"))
avail158 <- read_if(file.path(local_root, "GSE158722", "01_mt_audit",
                              "raw_qc_patient_mt_availability.csv"))
forensic <- read_if(file.path(local_root, "00_forensic", "current_cluster_sample_metrics.csv"))
status <- read_if(file.path(local_root, "strategy_comparison", "strategy_run_status.csv"))
recs <- read_if(file.path(local_root, "strategy_comparison",
                          "recommended_strategy_by_dataset_and_lineage.csv"))
comparison <- read_if(file.path(local_root, "strategy_comparison",
                                "batch_strategy_comparison.csv"))
diag154 <- read_if(file.path(local_root, "GSE154600", "02_dominance",
                             "cluster_dominance_diagnostic_table.csv"))
diag158 <- read_if(file.path(local_root, "GSE158722", "02_dominance",
                             "cluster_dominance_diagnostic_table.csv"))

strong_lines <- function(dataset_id) {
  dataset_id_value <- dataset_id
  z <- forensic[dataset_id == dataset_id_value & dominance_label == "strong_sample_dominance"]
  if (!nrow(z)) return("- None identified.")
  paste0("- Cluster ", z$cluster, ": ", z$dominant_sample,
         " (", sprintf("%.1f%%", 100 * z$dominant_sample_fraction),
         ", n=", z$n_cells, ")")
}
strategy_lines <- if (nrow(recs)) {
  paste0("- ", recs$dataset_id, " / ", recs$lineage, ": `",
         recs$recommended_strategy, "` - ", recs$rationale)
} else "- Strategy evaluation is unavailable."

interpretation_lines <- function(dt, label) {
  if (!nrow(dt) || !"likely_interpretation" %in% names(dt)) return("- Not available.")
  counts <- dt[, .(clusters = paste(cluster, collapse = ", "), n = .N),
               by = likely_interpretation]
  paste0("- ", label, " / `", counts$likely_interpretation, "` (n=", counts$n,
         "): clusters ", counts$clusters)
}

tumor_candidate_lines <- function(dt, label) {
  if (!nrow(dt) || !"likely_interpretation" %in% names(dt)) return("- Not available.")
  z <- dt[likely_interpretation == "likely_patient_specific_tumor_state"]
  if (!nrow(z)) return(paste0("- ", label, ": no cluster passed the conservative candidate rule."))
  paste0("- ", label, " cluster ", z$cluster, ": dominant sample ", z$dominant_sample,
         " (", sprintf("%.1f%%", 100 * z$dominant_sample_fraction),
         "); retained as uncorrected epithelial/tumor-state candidate.")
}

strategy_effect_lines <- function() {
  if (!nrow(recs) || !nrow(comparison)) return("- Strategy effect metrics unavailable.")
  out <- list()
  for (i in seq_len(nrow(recs))) {
    dataset_value <- recs$dataset_id[[i]]
    lineage_value <- recs$lineage[[i]]
    strategy_value <- recs$recommended_strategy[[i]]
    selected <- comparison[dataset_id == dataset_value & lineage == lineage_value &
                             strategy == strategy_value]
    baseline <- comparison[dataset_id == dataset_value & lineage == lineage_value &
                             strategy == "A_uncorrected"]
    if (!nrow(selected) || !nrow(baseline)) next
    out[[length(out) + 1L]] <- paste0(
      "- ", dataset_value, " / ", lineage_value, " / `", strategy_value, "`: ",
      "median dominant-sample fraction ", signif(baseline$median_cluster_dominant_sample_fraction[[1L]], 3),
      " -> ", signif(selected$median_cluster_dominant_sample_fraction[[1L]], 3),
      "; entropy ", signif(baseline$mean_cluster_normalized_entropy[[1L]], 3),
      " -> ", signif(selected$mean_cluster_normalized_entropy[[1L]], 3),
      "; ARI vs A=", signif(selected$ARI_vs_uncorrected[[1L]], 3),
      "; rare-population retention=", signif(selected$rare_population_retention[[1L]], 3), "."
    )
  }
  if (length(out)) unlist(out, use.names = FALSE) else "- No comparable completed strategy rows."
}

rpca_lines <- function() {
  completed <- comparison[strategy == "C_RPCA"]
  effect <- if (nrow(completed)) paste0(
    "- Completed RPCA / ", completed$dataset_id, " / ", completed$lineage,
    ": median dominant-sample fraction=",
    signif(completed$median_cluster_dominant_sample_fraction, 3),
    "; entropy=", signif(completed$mean_cluster_normalized_entropy, 3),
    "; ARI vs A=", signif(completed$ARI_vs_uncorrected, 3),
    "; rare-population retention=", signif(completed$rare_population_retention, 3), "."
  ) else character()
  failed <- status[strategy == "C_RPCA" & grepl("^FAILED", get("status"))]
  blockers <- if (nrow(failed)) paste0(
    "- RPCA blocker / ", failed$dataset_id, " / ", failed$lineage,
    ": ", failed$message
  ) else character()
  c(effect, blockers)
}

report <- c(
  "# Final diagnostic report", "",
  "## Scope and safeguards", "",
  "Only GSE147082, GSE154600, and GSE158722 were processed. GSE154763 was never run. Existing stage 06 `results/` and local result folders were not overwritten. No final `cell_type` or `cell_subtype` was assigned. All marker tests used the RNA assay.", "",
  "## 1-3. Why percent.mt was zero, whether it is usable, and feature sources", "",
  paste0("- GSE147082: the prepared row names use the R `make.names` dot form (`MT.*`), while the old workflow searched only `^MT-`. The repaired calculation explicitly used ",
         if (nrow(mt147)) mt147$n_mt_features_used[[1L]] else "audited", " mitochondrial features. Availability: ",
         if (nrow(mt147)) mt147$status[[1L]] else "not recorded", ". Retained after repaired mt QC: ",
         if (nrow(ret147)) paste0(ret147$n_after_repaired_mt_qc[[1L]], "/", ret147$n_current_singlets[[1L]]) else "not recorded", "."),
  paste0("- GSE158722: the prepared common-gene matrix contains no mitochondrial features, but original per-patient raw files were audited. Patients with `MT-` features were explicitly recalculated from their own raw counts; patients lacking credible mt features retain `NA` and skip mt filtering. No genes were added to the prepared matrix. Availability: ",
         if (nrow(mt158)) mt158$status[[1L]] else "not recorded", "; cell-level available fraction: ",
         if (nrow(mt158)) signif(mt158$fraction_cells_percent_mt_available[[1L]], 4) else "not recorded",
         if (nrow(avail158) && any(avail158$percent_mt_available == FALSE, na.rm = TRUE))
           paste0("; unavailable raw-source patients: ",
                  paste(avail158[percent_mt_available == FALSE, patient_id], collapse = ", ")) else "",
         "."), "",
  "## 4. Strong sample-dominant clusters in GSE154600", "",
  strong_lines("GSE154600"), "",
  "## 5. Strong patient/timepoint-dominant clusters in GSE158722", "",
  strong_lines("GSE158722"), "",
  "## 6. Biological versus technical interpretation", "",
  "Dominance is an audit flag, not a deletion rule. Strong clusters carrying epithelial/tumor programs (for example EPCAM/KRT/WFDC2/MSLN) are retained as candidate patient-specific malignant states. Clusters whose dominance coincides with abnormal nCount, nFeature, percent.mt, or doublet score are flagged as technical/QC suspects. Shared immune and stromal lineages are evaluated with correction sensitivity analyses.", "",
  interpretation_lines(diag154, "GSE154600"),
  interpretation_lines(diag158, "GSE158722"), "",
  "## 7-8. Harmony/RPCA effects and lineage-preserving choices", "",
  if (nrow(status)) paste0("- Strategy status counts: ", paste(names(table(status$status)), table(status$status), collapse = "; ")) else "- Strategy run status unavailable.",
  strategy_lines, "",
  strategy_effect_lines(), "",
  "### RPCA completed comparisons and blockers", "",
  rpca_lines(), "",
  "## 9. Patient-specific tumor clusters to retain", "",
  tumor_candidate_lines(diag154, "GSE154600"),
  tumor_candidate_lines(diag158, "GSE158722"),
  "All epithelial-like uncorrected clusters remain available as the primary view. Corrected epithelial results, when successful, are parallel shared-state sensitivity outputs and never replace the uncorrected result.", "",
  "## 10. Version for subsequent manual cell-type review", "",
  "Use the strategy listed per dataset and non-epithelial lineage in `strategy_comparison/recommended_strategy_by_dataset_and_lineage.csv`. For epithelial-like cells, use `A_uncorrected`. The blank `manual_annotation_template.csv` files are the handoff point.", "",
  "## 11. Remaining limitations", "",
  "- iLISI and kBET remain skipped when their packages are unavailable; their status is non-blocking and recorded.",
  "- GSE158722 current-object PC/UMAP metadata audit may be memory-guarded; newly generated lineage strategy embeddings provide the sensitivity views.",
  "- Broad lineages are provisional diagnostic strata, not final annotations.",
  "- Marker overlap is a lineage-level gene-set sensitivity metric; cluster identities are not assumed to map one-to-one after correction.", "",
  "- The task description flags a possible T61/T77 identity conflict, but T61 is absent from the supplied GSE154600 result metadata. This package therefore retains source sample IDs and makes no resolved patient-identity claim for T77.", "",
  "## Reproducibility", "",
  "Run `workflow/scripts/run_diagnostics_v2.ps1`. Large RDS/count payloads remain local under `diagnostics_v2/objects` and are excluded from GitHub."
)
writeLines(report, file.path(local_root, "final_diagnostic_report.md"), useBytes = TRUE)

readme <- c(
  "# diagnostics_v2 review package", "",
  "This small package contains forensic audits, mitochondrial-QC decisions, sample-dominance tables, A/B/C strategy summaries, RNA marker tables, plots, logs, and blank manual annotation templates.", "",
  "Start with [`final_diagnostic_report.md`](final_diagnostic_report.md), then inspect:", "",
  "- `00_forensic/current_result_snapshot.md`",
  "- `GSE147082/01_mt_audit/mt_qc_decision.md`",
  "- `GSE158722/01_mt_audit/mt_qc_decision.md`",
  "- `GSE154600/02_dominance/cluster_dominance_diagnostic_table.csv`",
  "- `GSE158722/02_dominance/cluster_dominance_diagnostic_table.csv`",
  "- `GSE158722/02_dominance/patient_timepoint_cluster_correspondence.csv`",
  "- `strategy_comparison/batch_strategy_decision_report.md`",
  "- each dataset's `05_markers/top20_markers_per_cluster.csv` and blank `manual_annotation_template.csv`.", "",
  "Large RDS files, count matrices, per-cell embeddings, and per-cell QC/lineage assignments are intentionally local-only."
)
writeLines(readme, file.path(local_root, "README.md"), useBytes = TRUE)

objects <- read_if(file.path(local_root, "00_forensic", "object_inventory.csv"))
summary <- list(
  task = "patient dominance and mitochondrial QC diagnostics_v2",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  datasets = c("GSE147082", "GSE154600", "GSE158722"),
  excluded_dataset = "GSE154763",
  original_stage06_overwritten = FALSE,
  final_cell_type_assigned = FALSE,
  marker_assay = "RNA",
  object_inventory = if (nrow(objects)) as.data.frame(objects) else NULL,
  strategy_status = if (nrow(status)) as.data.frame(status) else NULL
)
jsonlite::write_json(summary, file.path(local_root, "run_summary.json"),
                     auto_unbox = TRUE, pretty = TRUE, na = "null")
capture.output(sessionInfo(), file = file.path(local_root, "sessionInfo.txt"))

# Copy only review-sized artifacts. Per-cell tables, raw QC parts, large marker
# matrices, and all local RDS inputs remain outside Git.
files <- list.files(local_root, recursive = TRUE, full.names = TRUE)
rel <- substring(normalizePath(files, winslash = "/", mustWork = FALSE),
                 nchar(normalizePath(local_root, winslash = "/")) + 2L)
allowed <- grepl("\\.(csv|csv\\.gz|md|pdf|png|txt|json|log)$", rel, ignore.case = TRUE)
excluded <- grepl(paste(c(
  "(^|/)objects/", "raw_qc_parts/", "qc_metadata_mt_repaired",
  "provisional_broad_lineage_assignments", "cell_embedding_and_clusters",
  "evaluation_reduction", "all_cluster_markers", "cluster_average_expression"
), collapse = "|"), rel, ignore.case = TRUE)
size_ok <- file.info(files)$size <= 15 * 1024^2
selected <- files[allowed & !excluded & size_ok]
selected_rel <- rel[allowed & !excluded & size_ok]
for (i in seq_along(selected)) {
  target <- file.path(repo_package, selected_rel[[i]])
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(selected[[i]], target, overwrite = TRUE)) stop("Failed to copy ", selected[[i]])
  if (grepl("\\.(txt|log|md|json)$", target, ignore.case = TRUE)) {
    cleaned <- tryCatch(suppressWarnings({
      lines <- readLines(target, warn = FALSE, encoding = "UTF-8")
      sub("[[:blank:]]+$", "", lines)
    }), error = function(e) NULL)
    if (!is.null(cleaned)) writeLines(cleaned, target, useBytes = TRUE)
  }
}
target_files <- file.path(repo_package, selected_rel)
manifest <- data.table::data.table(
  path = gsub("\\\\", "/", selected_rel),
  size_bytes = file.info(target_files)$size,
  sha256 = vapply(target_files, digest::digest, character(1),
                  algo = "sha256", file = TRUE),
  checksum_algorithm = "SHA-256"
)
write_csv(manifest, file.path(repo_package, "review_package_manifest.csv"))
message("Review package generated at ", repo_package)
