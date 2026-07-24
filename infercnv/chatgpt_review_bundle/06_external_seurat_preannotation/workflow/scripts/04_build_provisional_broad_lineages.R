#!/usr/bin/env Rscript

required <- c("yaml", "data.table", "Matrix", "SingleCellExperiment",
              "SummarizedExperiment")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else "."
source(file.path(script_dir, "_diagnostics_v2_common.R"))

z <- read_diagnostics_config()
cfg <- z$cfg
set.seed(as.integer(cfg$project$random_seed))
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
prepared_root <- normalizePath(cfg$project$prepared_input_root, winslash = "/", mustWork = TRUE)
datasets <- split_arg(arg_value("--datasets", "GSE154600,GSE158722"))
assert_datasets(datasets, c("GSE154600", "GSE158722"))

author_to_broad <- function(x) {
  y <- rep(NA_character_, length(x))
  x <- toupper(as.character(x))
  y[grepl("EPI|TUM|CANC|OVAR|MAL", x)] <- "Epithelial_like"
  y[grepl("MYE|MAC|MONO|DEND|DC", x)] <- "Myeloid_like"
  y[grepl("(^|[^A-Z])T([^A-Z]|$)|NK", x)] <- "T_NK_like"
  y[grepl("B CELL|B_CELL|PLAS|^B$", x)] <- "B_Plasma_like"
  y[grepl("FIB|STROM|MES", x)] <- "Fibroblast_like"
  y[grepl("ENDO|VASC", x)] <- "Endothelial_like"
  y[grepl("CYCL|PROLIF", x)] <- "Cycling_like"
  y
}

score_modules <- function(counts, marker_sets) {
  genes_upper <- toupper(rownames(counts))
  lib <- Matrix::colSums(counts)
  lib[!is.finite(lib) | lib <= 0] <- 1
  score_list <- lapply(names(marker_sets), function(lineage) {
    idx <- which(genes_upper %in% toupper(marker_sets[[lineage]]))
    if (!length(idx)) return(rep(NA_real_, ncol(counts)))
    scaled <- counts[idx, , drop = FALSE] %*%
      Matrix::Diagonal(x = 10000 / lib)
    Matrix::colMeans(log1p(scaled))
  })
  scores <- do.call(cbind, score_list)
  colnames(scores) <- names(marker_sets)
  scores
}

for (dataset_id in datasets) {
  message("===== ", dataset_id, " =====")
  ds_cfg <- cfg$datasets[[dataset_id]]
  prepared_path <- file.path(prepared_root, ds_cfg$prepared_sce)
  out <- file.path(data_root, "diagnostics_v2", dataset_id, "03_broad_lineages")
  object_out <- file.path(data_root, "diagnostics_v2", "objects", dataset_id,
                          "lineage_inputs")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(object_out, recursive = TRUE, showWarnings = FALSE)

  message("Reading prepared counts: ", prepared_path)
  sce <- readRDS(prepared_path)
  counts <- SummarizedExperiment::assay(sce, "counts")
  sce_md <- as.data.frame(SummarizedExperiment::colData(sce))
  sce_md$cell_id <- colnames(counts)
  qc <- read_repaired_qc(data_root, dataset_id)
  qc <- qc[diagnostics_v2_qc_pass %in% TRUE]
  idx <- match(qc$cell_id, sce_md$cell_id)
  if (anyNA(idx)) stop(dataset_id, ": ", sum(is.na(idx)), " QC cells missing from prepared SCE")
  counts <- counts[, idx, drop = FALSE]
  sce_md <- sce_md[idx, , drop = FALSE]
  rownames(sce_md) <- sce_md$cell_id
  rm(sce); gc()

  qidx <- match(sce_md$cell_id, qc$cell_id)
  for (column in setdiff(names(qc), names(sce_md))) sce_md[[column]] <- qc[[column]][qidx]
  if (!"sample_id" %in% names(sce_md)) sce_md$sample_id <- qc$sample_id[qidx]
  if (!"patient_id" %in% names(sce_md)) sce_md$patient_id <- sce_md$sample_id
  if (!"timepoint" %in% names(sce_md)) sce_md$timepoint <- NA_character_

  scores <- score_modules(counts, broad_marker_sets)
  scores_for_rank <- scores
  scores_for_rank[!is.finite(scores_for_rank)] <- -Inf
  top_index <- max.col(scores_for_rank, ties.method = "first")
  top_score <- scores_for_rank[cbind(seq_len(nrow(scores_for_rank)), top_index)]
  second_score <- apply(scores_for_rank, 1L, function(x) sort(x, decreasing = TRUE)[[2L]])
  module_lineage <- colnames(scores_for_rank)[top_index]
  module_lineage[!is.finite(top_score) | top_score < 0.05 |
                   (top_score - second_score) < 0.03] <- "Uncertain"

  author_source <- rep(NA_character_, nrow(sce_md))
  for (column in intersect(c("celltype", "hpca.celltype", "encode.celltype"), names(sce_md))) {
    empty <- is.na(author_source) | !nzchar(author_source)
    author_source[empty] <- as.character(sce_md[[column]][empty])
  }
  author_lineage <- if (dataset_id == "GSE154600") author_to_broad(author_source) else rep(NA_character_, nrow(sce_md))
  final_lineage <- module_lineage
  use_author <- !is.na(author_lineage)
  final_lineage[use_author] <- author_lineage[use_author]
  margin <- top_score - second_score
  confidence <- data.table::fcase(
    use_author, "high_author_supported",
    final_lineage == "Uncertain", "low_uncertain",
    margin >= 0.25, "high_module_margin",
    margin >= 0.10, "moderate_module_margin",
    default = "low_module_margin"
  )
  prepared_genes_upper <- toupper(rownames(counts))
  evidence_by_lineage <- vapply(names(broad_marker_sets), function(lin) {
    present <- broad_marker_sets[[lin]][
      toupper(broad_marker_sets[[lin]]) %in% prepared_genes_upper
    ]
    paste(present, collapse = ";")
  }, character(1))
  evidence <- unname(evidence_by_lineage[final_lineage])
  evidence[is.na(evidence)] <- ""

  assignments <- data.table::data.table(
    dataset_id = dataset_id,
    cell_id = sce_md$cell_id,
    sample_id = as.character(sce_md$sample_id),
    patient_id = as.character(sce_md$patient_id),
    timepoint = as.character(sce_md$timepoint),
    provisional_broad_lineage = final_lineage,
    broad_lineage_confidence = confidence,
    broad_lineage_marker_evidence = evidence,
    module_top_score = top_score,
    module_score_margin = margin,
    author_broad_lineage = author_lineage,
    author_label_source = author_source,
    nCount_RNA = as.numeric(sce_md$nCount_RNA %||% Matrix::colSums(counts)),
    nFeature_RNA = as.numeric(sce_md$nFeature_RNA %||% Matrix::colSums(counts > 0)),
    percent.mt = as.numeric(sce_md$percent.mt %||% rep(NA_real_, ncol(counts))),
    doublet_score = as.numeric(sce_md$scDblFinder.score %||% rep(NA_real_, ncol(counts)))
  )
  for (j in seq_len(ncol(scores))) assignments[[paste0("score_", colnames(scores)[[j]])]] <- scores[, j]
  write_csv_gz(assignments, file.path(out, "provisional_broad_lineage_assignments.csv.gz"))

  summary <- assignments[, .(
    n_cells = .N,
    n_samples = data.table::uniqueN(sample_id),
    median_module_top_score = stats::median(module_top_score, na.rm = TRUE),
    median_score_margin = stats::median(module_score_margin, na.rm = TRUE),
    author_supported_fraction = mean(!is.na(author_broad_lineage))
  ), by = .(dataset_id, provisional_broad_lineage)]
  write_csv(summary, file.path(out, "provisional_broad_lineage_summary.csv"))

  sampling_rows <- list()
  eligible <- setdiff(unique(assignments$provisional_broad_lineage), "Uncertain")
  for (lineage in eligible) {
    line_md <- assignments[provisional_broad_lineage == lineage]
    if (nrow(line_md) < as.integer(cfg$analysis$broad_lineage_min_cells)) next
    chosen <- balanced_cells(line_md, as.integer(cfg$analysis$max_cells_per_lineage_strategy),
                             as.integer(cfg$project$random_seed))
    cidx <- match(chosen, colnames(counts))
    payload <- list(
      dataset_id = dataset_id,
      lineage = lineage,
      counts = counts[, cidx, drop = FALSE],
      metadata = as.data.frame(line_md[match(chosen, cell_id)]),
      creation_note = "Prepared counts were subset only; no gene or count values were changed."
    )
    path <- file.path(object_out, paste0(lineage, "_strategy_input.rds"))
    saveRDS(payload, path, compress = TRUE)
    sampling_rows[[length(sampling_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, provisional_broad_lineage = lineage,
      n_cells_available = nrow(line_md), n_cells_strategy_input = length(chosen),
      n_samples_available = data.table::uniqueN(line_md$sample_id),
      sampling = if (nrow(line_md) > length(chosen)) "balanced_cap_by_sample" else "all_cells",
      object_path = normalizePath(path, winslash = "/", mustWork = FALSE)
    )
    rm(payload); gc()
  }
  combined_md <- assignments[provisional_broad_lineage != "Uncertain"]
  if (nrow(combined_md) >= as.integer(cfg$analysis$broad_lineage_min_cells)) {
    chosen <- balanced_cells(combined_md,
                             as.integer(cfg$analysis$max_cells_per_lineage_strategy),
                             as.integer(cfg$project$random_seed))
    cidx <- match(chosen, colnames(counts))
    payload <- list(
      dataset_id = dataset_id,
      lineage = "Combined_broad_lineages",
      counts = counts[, cidx, drop = FALSE],
      metadata = as.data.frame(combined_md[match(chosen, cell_id)]),
      creation_note = paste(
        "Balanced cross-lineage sensitivity input.",
        "Prepared counts were subset only; no gene or count values were changed."
      )
    )
    path <- file.path(object_out, "Combined_broad_lineages_strategy_input.rds")
    saveRDS(payload, path, compress = TRUE)
    sampling_rows[[length(sampling_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      provisional_broad_lineage = "Combined_broad_lineages",
      n_cells_available = nrow(combined_md), n_cells_strategy_input = length(chosen),
      n_samples_available = data.table::uniqueN(combined_md$sample_id),
      sampling = if (nrow(combined_md) > length(chosen)) "balanced_cap_by_sample" else "all_cells",
      object_path = normalizePath(path, winslash = "/", mustWork = FALSE)
    )
    rm(payload); gc()
  }
  write_csv(data.table::rbindlist(sampling_rows, fill = TRUE),
            file.path(out, "lineage_strategy_input_audit.csv"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  rm(counts, sce_md, qc, scores, scores_for_rank, assignments); gc()
}

message("Provisional broad-lineage construction complete")
