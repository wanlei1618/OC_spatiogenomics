#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}

split_arg <- function(x) {
  z <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  z[nzchar(z)]
}

required <- c("data.table", "digest", "jsonlite", "SeuratObject")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                           FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

data_root <- normalizePath(
  arg_value("--data-root",
            "D:/OC_spatiogenomics/infercnv/external_seurat_preannotation"),
  winslash = "/", mustWork = TRUE
)
project_root <- normalizePath(
  arg_value("--project-root", "D:/OC_spatiogenomics"),
  winslash = "/", mustWork = TRUE
)
datasets <- split_arg(arg_value(
  "--datasets", "GSE147082,GSE154600,GSE158722"
))
allowed <- c("GSE147082", "GSE154600", "GSE158722")
if (!length(datasets) || any(!datasets %in% allowed)) {
  stop("Datasets must be a non-empty subset of: ", paste(allowed, collapse = ", "))
}
max_object_load_gb <- as.numeric(arg_value("--max-object-load-gb", "2.5"))
if (!is.finite(max_object_load_gb) || max_object_load_gb <= 0) {
  stop("--max-object-load-gb must be positive")
}

output_root <- file.path(data_root, "diagnostics_v2", "00_forensic")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  data.table::fwrite(x, path, na = "NA")
}

safe_sha256 <- function(path) {
  tryCatch(
    digest::digest(file = path, algo = "sha256", serialize = FALSE),
    error = function(e) NA_character_
  )
}

read_status <- function(dataset_id) {
  path <- file.path(data_root, dataset_id, "logs", "run_status.json")
  if (!file.exists(path)) {
    return(list(path = path, status = NA_character_, n_after_qc = NA_real_,
                finished_at = NA_character_))
  }
  x <- tryCatch(jsonlite::read_json(path, simplifyVector = TRUE),
                error = function(e) list())
  list(
    path = path,
    status = as.character(x$status %||% NA_character_),
    n_after_qc = as.numeric(x$n_after_qc %||% NA_real_),
    finished_at = as.character(x$finished_at %||% NA_character_)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

find_object_candidates <- function(dataset_id) {
  name <- paste0(dataset_id, "_preannotation.rds")
  expected <- file.path(data_root, dataset_id, "objects", name)
  if (file.exists(expected)) return(normalizePath(expected, winslash = "/"))
  hits <- list.files(project_root, pattern = paste0("^", name, "$"),
                     recursive = TRUE, full.names = TRUE,
                     ignore.case = FALSE, include.dirs = FALSE)
  unique(normalizePath(hits, winslash = "/", mustWork = FALSE))
}

object_rows <- list()
selected_objects <- setNames(rep(NA_character_, length(datasets)), datasets)

for (dataset_id in datasets) {
  expected <- normalizePath(
    file.path(data_root, dataset_id, "objects",
              paste0(dataset_id, "_preannotation.rds")),
    winslash = "/", mustWork = FALSE
  )
  candidates <- find_object_candidates(dataset_id)
  status <- read_status(dataset_id)
  if (!length(candidates)) {
    object_rows[[length(object_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, candidate_path = NA_character_,
      size_bytes = NA_real_, modified_time = NA_character_, sha256 = NA_character_,
      exact_expected_path = FALSE, selected = FALSE,
      selection_reason = "no_candidate_found", run_status_file = status$path,
      run_status = status$status, run_status_n_after_qc = status$n_after_qc,
      run_status_finished_at = status$finished_at
    )
    next
  }
  exact <- candidates == expected
  selected <- if (sum(exact) == 1L) which(exact) else integer()
  if (length(selected)) selected_objects[[dataset_id]] <- candidates[[selected]]
  for (i in seq_along(candidates)) {
    info <- file.info(candidates[[i]])
    reason <- if (i %in% selected) {
      "exact_expected_path"
    } else if (length(candidates) > 1L && !any(exact)) {
      "multiple_candidates_no_automatic_selection"
    } else {
      "not_selected"
    }
    message("Hashing object: ", candidates[[i]])
    object_rows[[length(object_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      candidate_path = candidates[[i]],
      size_bytes = unname(info$size),
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
      sha256 = safe_sha256(candidates[[i]]),
      exact_expected_path = exact[[i]],
      selected = i %in% selected,
      selection_reason = reason,
      run_status_file = status$path,
      run_status = status$status,
      run_status_n_after_qc = status$n_after_qc,
      run_status_finished_at = status$finished_at
    )
  }
}

object_inventory <- data.table::rbindlist(object_rows, fill = TRUE)
write_csv(object_inventory, file.path(output_root, "object_inventory.csv"))

result_hash_rows <- list()
for (dataset_id in datasets) {
  dataset_root <- file.path(data_root, dataset_id)
  files <- list.files(dataset_root, recursive = TRUE, full.names = TRUE,
                      include.dirs = FALSE)
  files <- files[!grepl("[/\\\\]objects[/\\\\]", files)]
  files <- files[!grepl("[/\\\\]diagnostics_v2[/\\\\]", files)]
  for (path in sort(files)) {
    info <- file.info(path)
    message("Hashing result: ", path)
    result_hash_rows[[length(result_hash_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      relative_path = substring(normalizePath(path, winslash = "/"),
                                nchar(data_root) + 2L),
      size_bytes = unname(info$size),
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
      sha256 = safe_sha256(path)
    )
  }
}
current_hashes <- data.table::rbindlist(result_hash_rows, fill = TRUE)
write_csv(current_hashes, file.path(output_root, "current_result_sha256.csv"))

cramers_v <- function(tab) {
  tab <- as.matrix(tab)
  n <- sum(tab)
  if (!length(tab) || n <= 0 || nrow(tab) < 2L || ncol(tab) < 2L) return(NA_real_)
  expected <- outer(rowSums(tab), colSums(tab)) / n
  ok <- expected > 0
  chi <- sum(((tab[ok] - expected[ok]) ^ 2) / expected[ok])
  denom <- n * min(nrow(tab) - 1L, ncol(tab) - 1L)
  if (denom <= 0) NA_real_ else sqrt(chi / denom)
}

sample_annotations <- function(dataset_id, sample_ids) {
  patient <- sample_ids
  timepoint <- rep(NA_character_, length(sample_ids))
  treatment <- rep(NA_character_, length(sample_ids))
  if (dataset_id == "GSE158722") {
    patient <- sub("_.*$", "", sample_ids)
    timepoint <- sub("^[^_]+_", "", sample_ids)
    timepoint[timepoint == sample_ids] <- NA_character_
    treatment[grepl("pre", timepoint, ignore.case = TRUE)] <- "pre"
    treatment[grepl("post", timepoint, ignore.case = TRUE)] <- "post"
  }
  data.frame(sample_id = sample_ids, patient_id = patient,
             timepoint = timepoint, treatment = treatment)
}

aggregate_axis <- function(cs, annotation, axis) {
  map <- annotation[, c("sample_id", axis), drop = FALSE]
  names(map)[2L] <- "axis_value"
  z <- merge(cs, map, by = "sample_id", all.x = TRUE)
  z <- z[!is.na(z$axis_value) & nzchar(z$axis_value), , drop = FALSE]
  if (!nrow(z) || length(unique(z$axis_value)) < 2L) return(NA_real_)
  tab <- xtabs(n_cells ~ seurat_cluster + axis_value, z)
  cramers_v(tab)
}

metric_rows <- list()
for (dataset_id in datasets) {
  path <- file.path(data_root, dataset_id, "02_clustering",
                    "cluster_by_sample_counts.csv")
  if (!file.exists(path)) next
  cs <- data.table::fread(path)
  names(cs) <- sub("^analysis_", "", names(cs))
  if (!all(c("sample_id", "seurat_cluster", "n_cells") %in% names(cs))) {
    stop(dataset_id, ": unexpected cluster-by-sample columns")
  }
  cs$n_cells <- as.numeric(cs$n_cells)
  all_samples <- sort(unique(as.character(cs$sample_id)))
  annotation <- sample_annotations(dataset_id, all_samples)
  sample_tab <- xtabs(n_cells ~ seurat_cluster + sample_id, cs)
  sample_v <- cramers_v(sample_tab)
  patient_v <- aggregate_axis(cs, annotation, "patient_id")
  timepoint_v <- aggregate_axis(cs, annotation, "timepoint")
  treatment_v <- aggregate_axis(cs, annotation, "treatment")
  total_sample_count <- length(all_samples)
  for (cluster_id in sort(unique(as.character(cs$seurat_cluster)))) {
    z <- cs[as.character(cs$seurat_cluster) == cluster_id & cs$n_cells > 0, ]
    n <- sum(z$n_cells)
    p <- z$n_cells / n
    dominant_i <- which.max(z$n_cells)
    entropy <- if (total_sample_count > 1L) {
      -sum(p * log(p)) / log(total_sample_count)
    } else {
      NA_real_
    }
    dominant_fraction <- p[[dominant_i]]
    dominance_label <- if (n >= 100 && dominant_fraction >= 0.80) {
      "strong_sample_dominance"
    } else if (dominant_fraction >= 0.60) {
      "moderate_sample_dominance"
    } else {
      "mixed"
    }
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      cluster = cluster_id,
      n_cells = n,
      n_samples = nrow(z),
      total_samples_in_dataset = total_sample_count,
      dominant_sample = as.character(z$sample_id[[dominant_i]]),
      dominant_sample_n = z$n_cells[[dominant_i]],
      dominant_sample_fraction = dominant_fraction,
      normalized_shannon_entropy = entropy,
      simpson_diversity = 1 - sum(p ^ 2),
      effective_sample_number = 1 / sum(p ^ 2),
      dominance_label = dominance_label,
      cluster_sample_cramers_v = sample_v,
      cluster_patient_cramers_v = patient_v,
      cluster_timepoint_cramers_v = timepoint_v,
      cluster_treatment_cramers_v = treatment_v
    )
  }
}
cluster_metrics <- data.table::rbindlist(metric_rows, fill = TRUE)
write_csv(cluster_metrics,
          file.path(output_root, "current_cluster_sample_metrics.csv"))

eta_squared <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- as.factor(group[ok])
  if (length(x) < 3L || nlevels(group) < 2L || stats::var(x) == 0) return(NA_real_)
  mu <- mean(x)
  between <- sum(tapply(x, group, function(v) length(v) * (mean(v) - mu) ^ 2))
  total <- sum((x - mu) ^ 2)
  if (total <= 0) NA_real_ else between / total
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

metadata_column <- function(md, candidates) {
  z <- candidates[candidates %in% names(md)]
  if (length(z)) md[[z[[1L]]]] else rep(NA, nrow(md))
}

pc_rows <- list()
for (dataset_id in datasets) {
  path <- selected_objects[[dataset_id]]
  if (is.na(path) || !file.exists(path)) {
    pc_rows[[length(pc_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, status = "NO_SELECTED_OBJECT", PC = NA_integer_,
      variance_explained = NA_real_, sample_eta_squared = NA_real_,
      patient_eta_squared = NA_real_, timepoint_eta_squared = NA_real_,
      treatment_eta_squared = NA_real_, nCount_spearman = NA_real_,
      nFeature_spearman = NA_real_, percent_mt_spearman = NA_real_,
      doublet_score_spearman = NA_real_, message = "No unambiguous object candidate"
    )
    next
  }
  size_gb <- file.info(path)$size / 1024 ^ 3
  if (size_gb > max_object_load_gb) {
    pc_rows[[length(pc_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, status = "SKIPPED_MEMORY_GUARD", PC = NA_integer_,
      variance_explained = NA_real_, sample_eta_squared = NA_real_,
      patient_eta_squared = NA_real_, timepoint_eta_squared = NA_real_,
      treatment_eta_squared = NA_real_, nCount_spearman = NA_real_,
      nFeature_spearman = NA_real_, percent_mt_spearman = NA_real_,
      doublet_score_spearman = NA_real_,
      message = sprintf("Object %.2f GB exceeds %.2f GB load guard", size_gb,
                        max_object_load_gb)
    )
    next
  }
  message("Loading object for PC audit: ", path)
  obj <- tryCatch(readRDS(path), error = identity)
  if (inherits(obj, "error")) {
    pc_rows[[length(pc_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, status = "OBJECT_LOAD_FAILED", PC = NA_integer_,
      variance_explained = NA_real_, sample_eta_squared = NA_real_,
      patient_eta_squared = NA_real_, timepoint_eta_squared = NA_real_,
      treatment_eta_squared = NA_real_, nCount_spearman = NA_real_,
      nFeature_spearman = NA_real_, percent_mt_spearman = NA_real_,
      doublet_score_spearman = NA_real_, message = conditionMessage(obj)
    )
    next
  }
  if (!"pca" %in% names(obj@reductions)) {
    pc_rows[[length(pc_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, status = "PCA_NOT_FOUND", PC = NA_integer_,
      variance_explained = NA_real_, sample_eta_squared = NA_real_,
      patient_eta_squared = NA_real_, timepoint_eta_squared = NA_real_,
      treatment_eta_squared = NA_real_, nCount_spearman = NA_real_,
      nFeature_spearman = NA_real_, percent_mt_spearman = NA_real_,
      doublet_score_spearman = NA_real_, message = "Object lacks pca reduction"
    )
    rm(obj); gc()
    next
  }
  emb <- SeuratObject::Embeddings(obj[["pca"]])
  md <- obj@meta.data[rownames(emb), , drop = FALSE]
  sdev <- obj[["pca"]]@stdev
  variance <- sdev ^ 2 / sum(sdev ^ 2)
  n_pc <- min(30L, ncol(emb), length(variance))
  sample_id <- metadata_column(md, c("analysis_sample_id", "sample_id"))
  patient_id <- metadata_column(md, c("patient_id"))
  timepoint <- metadata_column(md, c("timepoint"))
  treatment <- metadata_column(md, c("treatment"))
  n_count <- as.numeric(metadata_column(md, c("nCount_RNA")))
  n_feature <- as.numeric(metadata_column(md, c("nFeature_RNA")))
  percent_mt <- as.numeric(metadata_column(md, c("percent.mt")))
  doublet <- as.numeric(metadata_column(md,
                                        c("scDblFinder.score", "doublet_score")))
  for (j in seq_len(n_pc)) {
    x <- emb[, j]
    pc_rows[[length(pc_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, status = "COMPLETE", PC = j,
      variance_explained = variance[[j]],
      sample_eta_squared = eta_squared(x, sample_id),
      patient_eta_squared = eta_squared(x, patient_id),
      timepoint_eta_squared = eta_squared(x, timepoint),
      treatment_eta_squared = eta_squared(x, treatment),
      nCount_spearman = safe_cor(x, n_count),
      nFeature_spearman = safe_cor(x, n_feature),
      percent_mt_spearman = safe_cor(x, percent_mt),
      doublet_score_spearman = safe_cor(x, doublet),
      message = ""
    )
  }
  rm(obj, emb, md); gc()
}
pc_metrics <- data.table::rbindlist(pc_rows, fill = TRUE)
write_csv(pc_metrics,
          file.path(output_root, "current_pc_sample_association.csv"))

strong <- cluster_metrics[dominance_label == "strong_sample_dominance",
                          .N, by = dataset_id]
moderate <- cluster_metrics[dominance_label == "moderate_sample_dominance",
                            .N, by = dataset_id]
summary_lines <- c(
  "# Current stage 06 forensic snapshot",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "This is a read-only snapshot. Existing stage 06 results and RDS objects were not modified.",
  "",
  "## Object selection",
  ""
)
for (ds_id in datasets) {
  z <- object_inventory[dataset_id == ds_id]
  selected <- z[selected == TRUE]
  line <- if (nrow(selected)) {
    sprintf("- %s: `%s` (%0.2f GB), SHA-256 `%s`", ds_id,
            selected$candidate_path[[1L]], selected$size_bytes[[1L]] / 1024 ^ 3,
            selected$sha256[[1L]])
  } else {
    sprintf("- %s: no unambiguous object selected", ds_id)
  }
  summary_lines <- c(summary_lines, line)
}
summary_lines <- c(summary_lines, "", "## Current sample-dominance audit", "")
for (ds_id in datasets) {
  s <- strong[dataset_id == ds_id, N]
  m <- moderate[dataset_id == ds_id, N]
  if (!length(s)) s <- 0L
  if (!length(m)) m <- 0L
  z <- cluster_metrics[dataset_id == ds_id]
  v <- if (nrow(z)) unique(z$cluster_sample_cramers_v)[[1L]] else NA_real_
  summary_lines <- c(summary_lines,
                     sprintf("- %s: %d strong, %d moderate clusters; cluster-sample Cramer's V = %.4f",
                             ds_id, s, m, v))
}
summary_lines <- c(summary_lines, "", "## PC audit status", "")
for (ds_id in datasets) {
  z <- pc_metrics[dataset_id == ds_id]
  status <- paste(unique(z$status), collapse = ", ")
  msg <- paste(unique(z$message[nzchar(z$message)]), collapse = "; ")
  summary_lines <- c(summary_lines,
                     paste0("- ", ds_id, ": ", status,
                            if (nzchar(msg)) paste0(" (", msg, ")") else ""))
}
summary_lines <- c(
  summary_lines, "",
  "Dominance labels are audit flags, not deletion rules. Biological interpretation requires markers, lineage evidence, and QC context."
)
writeLines(summary_lines,
           file.path(output_root, "current_result_snapshot.md"), useBytes = TRUE)

message("Forensic audit complete: ", output_root)
