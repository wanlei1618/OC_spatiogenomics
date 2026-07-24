#!/usr/bin/env Rscript

required <- c("yaml", "data.table", "cluster")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else "."
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
datasets <- split_arg(arg_value("--datasets", "GSE154600,GSE158722"))
assert_datasets(datasets, c("GSE154600", "GSE158722"))
seed <- as.integer(cfg$project$random_seed)

comb2 <- function(x) x * (x - 1) / 2

adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n < 2L) return(NA_real_)
  a <- sum(comb2(tab)); b <- sum(comb2(rowSums(tab))); c <- sum(comb2(colSums(tab)))
  expected <- b * c / comb2(n)
  denom <- (b + c) / 2 - expected
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (a - expected) / denom
}

normalized_mi <- function(x, y) {
  tab <- table(x, y)
  pxy <- tab / sum(tab); px <- rowSums(pxy); py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(i) {
    r <- nz[i, 1L]; c <- nz[i, 2L]
    pxy[r, c] * log(pxy[r, c] / (px[[r]] * py[[c]]))
  }, numeric(1)))
  hx <- -sum(px[px > 0] * log(px[px > 0])); hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (!is.finite(hx + hy) || hx + hy == 0) return(NA_real_)
  2 * mi / (hx + hy)
}

mean_silhouette <- function(mat, labels, max_cells = 3000L) {
  ok <- complete.cases(mat) & !is.na(labels)
  mat <- as.matrix(mat[ok, , drop = FALSE]); labels <- as.character(labels[ok])
  if (nrow(mat) < 3L || data.table::uniqueN(labels) < 2L) return(NA_real_)
  set.seed(seed)
  if (nrow(mat) > max_cells) {
    idx <- sample.int(nrow(mat), max_cells)
    mat <- mat[idx, , drop = FALSE]; labels <- labels[idx]
  }
  groups <- as.integer(factor(labels))
  if (any(table(groups) < 2L)) {
    keep <- groups %in% as.integer(names(table(groups)[table(groups) >= 2L]))
    mat <- mat[keep, , drop = FALSE]; groups <- as.integer(factor(labels[keep]))
  }
  if (nrow(mat) < 3L || length(unique(groups)) < 2L) return(NA_real_)
  mean(cluster::silhouette(groups, stats::dist(mat))[, "sil_width"])
}

numeric_eta_squared <- function(cluster, value) {
  ok <- !is.na(cluster) & is.finite(value)
  cluster <- as.character(cluster[ok]); value <- value[ok]
  if (length(value) < 3L || data.table::uniqueN(cluster) < 2L) return(NA_real_)
  grand <- mean(value)
  ss_between <- sum(vapply(split(value, cluster), function(z) length(z) * (mean(z) - grand)^2,
                           numeric(1)))
  ss_total <- sum((value - grand)^2)
  if (!is.finite(ss_total) || ss_total == 0) return(NA_real_)
  ss_between / ss_total
}

weighted_purity <- function(cluster, label) {
  tab <- table(cluster, label, useNA = "no")
  if (!sum(tab)) return(NA_real_)
  sum(apply(tab, 1L, max)) / sum(tab)
}

calc_rare_population_retention <- function(baseline, corrected,
                                           rare_fraction = 0.05, min_cells = 20L) {
  joined <- merge(baseline, corrected, by = "cell_id",
                  suffixes = c("_base", "_strategy"))
  if (!nrow(joined)) return(NA_real_)
  base_counts <- table(joined$cluster_base)
  rare <- names(base_counts)[base_counts >= min_cells &
                               base_counts / sum(base_counts) <= rare_fraction]
  if (!length(rare)) return(NA_real_)
  # Mean best Jaccard overlap quantifies whether each rare uncorrected cluster
  # remains a discrete population rather than being swallowed after correction.
  mean(vapply(rare, function(base_cluster) {
    base_cells <- joined$cell_id[joined$cluster_base == base_cluster]
    candidates <- split(joined$cell_id, joined$cluster_strategy)
    max(vapply(candidates, function(strategy_cells) {
      length(intersect(base_cells, strategy_cells)) /
        length(union(base_cells, strategy_cells))
    }, numeric(1)))
  }, numeric(1)))
}

scope_files <- list()
for (dataset_id in datasets) {
  root <- file.path(data_root, "diagnostics_v2", dataset_id, "04_strategies")
  files <- list.files(root, pattern = "evaluation_reduction\\.csv\\.gz$",
                      recursive = TRUE, full.names = TRUE)
  scope_files <- c(scope_files, files)
}
if (!length(scope_files)) stop("No completed strategy reduction files from step 05")

records <- list(); qc_records <- list(); cell_tables <- list()
for (path in scope_files) {
  dt <- data.table::fread(path, showProgress = FALSE)
  dataset_id <- unique(dt$dataset_id)[[1L]]
  scope <- unique(dt$analysis_scope_lineage)[[1L]]
  strategy <- unique(dt$strategy)[[1L]]
  key <- paste(dataset_id, scope, strategy, sep = "||")
  cell_tables[[key]] <- dt[, .(cell_id, cluster)]
  reduction_cols <- grep("^(PC_|harmony_|integratedRPCA_|integrated\\.rpca_|rpca_|[A-Za-z]+_[0-9]+$)",
                         names(dt), value = TRUE)
  if (!length(reduction_cols)) reduction_cols <- tail(names(dt), as.integer(cfg$analysis$dims_use))
  reduction_cols <- reduction_cols[vapply(dt[, ..reduction_cols], is.numeric, logical(1))]
  reduction_cols <- head(reduction_cols, as.integer(cfg$analysis$dims_use))
  mat <- as.matrix(dt[, ..reduction_cols])
  dom <- dominance_metrics(dt$cluster, dt$sample_id, dataset_id, strategy, scope)
  cluster_sil <- mean_silhouette(mat, dt$cluster)
  sample_sil <- mean_silhouette(mat, dt$sample_id)
  lineage_sil <- if (data.table::uniqueN(dt$provisional_broad_lineage) > 1L) {
    mean_silhouette(mat, dt$provisional_broad_lineage)
  } else NA_real_
  author_purity <- if (dataset_id == "GSE154600") {
    weighted_purity(dt$cluster, dt$author_label_source)
  } else NA_real_
  broad_purity <- weighted_purity(dt$cluster, dt$provisional_broad_lineage)
  ilisi_status <- if (requireNamespace("lisi", quietly = TRUE)) "AVAILABLE_NOT_RUN_FULL_MEMORY_GUARD" else "SKIPPED_PACKAGE_UNAVAILABLE"
  kbet_status <- if (requireNamespace("kBET", quietly = TRUE)) "AVAILABLE_NOT_RUN_FULL_MEMORY_GUARD" else "SKIPPED_PACKAGE_UNAVAILABLE"
  records[[key]] <- data.frame(
    dataset_id = dataset_id, lineage = scope, strategy = strategy,
    n_cells = nrow(dt), n_clusters = data.table::uniqueN(dt$cluster),
    median_cluster_dominant_sample_fraction = stats::median(dom$dominant_sample_fraction),
    mean_cluster_normalized_entropy = mean(dom$normalized_shannon_entropy, na.rm = TRUE),
    mean_cluster_effective_sample_number = mean(dom$effective_sample_number, na.rm = TRUE),
    sample_silhouette = sample_sil,
    broad_lineage_silhouette = lineage_sil,
    cluster_embedding_coherence = cluster_sil,
    broad_lineage_cluster_purity = broad_purity,
    author_celltype_cluster_purity = author_purity,
    cell_retention_vs_baseline = 1,
    iLISI = NA_real_, iLISI_status = ilisi_status,
    kBET = NA_real_, kBET_status = kbet_status,
    marker_coherence_status = "embedding_and_broad_marker_lineage_proxy; RNA markers rerun_in_step_07",
    stringsAsFactors = FALSE
  )
  expected_numeric <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "doublet_score")
  for (variable in expected_numeric) {
    available <- variable %in% names(dt)
    qc_records[[length(qc_records) + 1L]] <- data.frame(
      dataset_id, lineage = scope, strategy, variable,
      association_type = "eta_squared",
      association_value = if (available) numeric_eta_squared(dt$cluster, as.numeric(dt[[variable]])) else NA_real_,
      status = if (available) "COMPLETE" else "SKIPPED_METADATA_UNAVAILABLE"
    )
  }
  expected_categorical <- c("sample_id", "patient_id", "timepoint", "treatment")
  for (variable in expected_categorical) {
    available <- variable %in% names(dt) && any(!is.na(dt[[variable]]) & nzchar(as.character(dt[[variable]])))
    qc_records[[length(qc_records) + 1L]] <- data.frame(
      dataset_id, lineage = scope, strategy, variable,
      association_type = "cramers_v",
      association_value = if (available) cramers_v(dt$cluster, dt[[variable]]) else NA_real_,
      status = if (available) "COMPLETE" else "SKIPPED_METADATA_UNAVAILABLE"
    )
  }
}

comparison <- data.table::rbindlist(records, fill = TRUE)
comparison[, `:=`(ARI_vs_uncorrected = NA_real_, NMI_vs_uncorrected = NA_real_,
                  rare_population_retention = NA_real_)]
for (i in seq_len(nrow(comparison))) {
  row <- comparison[i]
  base_key <- paste(row$dataset_id, row$lineage, "A_uncorrected", sep = "||")
  key <- paste(row$dataset_id, row$lineage, row$strategy, sep = "||")
  if (!is.null(cell_tables[[base_key]]) && !is.null(cell_tables[[key]])) {
    joined <- merge(cell_tables[[base_key]], cell_tables[[key]], by = "cell_id",
                    suffixes = c("_base", "_strategy"))
    comparison[i, ARI_vs_uncorrected := adjusted_rand(joined$cluster_base, joined$cluster_strategy)]
    comparison[i, NMI_vs_uncorrected := normalized_mi(joined$cluster_base, joined$cluster_strategy)]
    comparison[i, rare_population_retention :=
                 calc_rare_population_retention(cell_tables[[base_key]], cell_tables[[key]])]
  }
}
qc_assoc <- data.table::rbindlist(qc_records, fill = TRUE)
qc_summary <- qc_assoc[association_type == "eta_squared",
                       .(max_cluster_qc_eta_squared = max(association_value, na.rm = TRUE)),
                       by = .(dataset_id, lineage, strategy)]
qc_summary[!is.finite(max_cluster_qc_eta_squared), max_cluster_qc_eta_squared := NA_real_]
comparison <- merge(comparison, qc_summary,
                    by = c("dataset_id", "lineage", "strategy"), all.x = TRUE)

comparison[, balanced_score :=
             mean_cluster_normalized_entropy - median_cluster_dominant_sample_fraction -
             pmax(sample_silhouette, 0, na.rm = TRUE) +
             0.5 * data.table::fcoalesce(ARI_vs_uncorrected, 1) +
             0.5 * data.table::fcoalesce(broad_lineage_cluster_purity, 0) -
             0.25 * (1 - data.table::fcoalesce(rare_population_retention, 1)) -
             0.25 * data.table::fcoalesce(max_cluster_qc_eta_squared, 0)]

recommendations <- comparison[lineage != "Combined_broad_lineages", {
  if (.BY$lineage == "Epithelial_like") {
    chosen <- "A_uncorrected"
    rationale <- paste(
      "Preserve the uncorrected patient-specific epithelial structure as the primary result.",
      "Harmony/RPCA remain parallel shared-state sensitivity outputs when completed."
    )
  } else {
    valid <- which(is.finite(balanced_score))
    chosen <- if (length(valid)) strategy[[valid[[which.max(balanced_score[valid])]]]] else "A_uncorrected"
    rationale <- paste0("Highest balanced score among completed strategies; score combines sample mixing, ",
                        "cluster preservation, broad-lineage purity, and QC association. Score=",
                        signif(max(balanced_score, na.rm = TRUE), 4), ".")
  }
  list(recommended_strategy = chosen, rationale = rationale,
       preserve_uncorrected_parallel = .BY$lineage == "Epithelial_like")
}, by = .(dataset_id, lineage)]

out <- file.path(data_root, "diagnostics_v2", "strategy_comparison")
write_csv(comparison, file.path(out, "batch_strategy_comparison.csv"))
write_csv(qc_assoc, file.path(out, "cluster_qc_and_metadata_association.csv"))
write_csv(recommendations, file.path(out, "recommended_strategy_by_dataset_and_lineage.csv"))

report <- c(
  "# Batch strategy decision report", "",
  "Strategies were compared within each provisional broad lineage. Lower sample dominance and sample silhouette were balanced against cluster/lineage preservation; mixing alone was never the selection rule.", "",
  "Epithelial-like cells always retain the uncorrected result as the primary patient-specific tumor-state view. Corrected epithelial results are sensitivity views only.", "",
  paste0("- iLISI: ", paste(unique(comparison$iLISI_status), collapse = "; ")),
  paste0("- kBET: ", paste(unique(comparison$kBET_status), collapse = "; ")),
  "- RNA marker coherence is finalized in step 07; this step uses embedding coherence and broad-marker lineage purity as an explicitly labeled proxy.", "",
  "## Recommendations", "",
  apply(recommendations, 1L, function(x) paste0("- ", x[["dataset_id"]], " / ", x[["lineage"]],
                                                ": `", x[["recommended_strategy"]], "` - ", x[["rationale"]]))
)
writeLines(report, file.path(out, "batch_strategy_decision_report.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo_06.txt"))
message("Batch strategy evaluation complete")
