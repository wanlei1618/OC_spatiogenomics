#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args) >= 1) args[[1]] else "D:/OC_spatiogenomics/infercnv/external_cell_annotations"
data_root <- normalizePath(data_root, winslash = "/", mustWork = FALSE)
raw_dir <- file.path(data_root, "raw", "GSE154600")
processed_dir <- file.path(data_root, "processed")
audit_dir <- file.path(data_root, "audit")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

sce_files <- list.files(raw_dir, pattern = "sample[0-9]+_sce\\.rds$", full.names = TRUE)
rows <- list()
for (path in sce_files) {
  sample_id <- sub("_sce\\.rds$", "", basename(path))
  sample_id <- sub("^sample", "T", sample_id)
  status <- "ok"
  message <- ""
  n_cells <- NA_integer_
  n_genes <- NA_integer_
  out_path <- file.path(processed_dir, paste0(tools::file_path_sans_ext(basename(path)), "_coldata.csv.gz"))
  tryCatch({
    sce <- readRDS(path)
    n_cells <- ncol(sce)
    n_genes <- nrow(sce)
    cd <- as.data.frame(SummarizedExperiment::colData(sce))
    cd$cell_id_original <- colnames(sce)
    cd$sample_id_author <- sample_id
    utils::write.csv(cd, gzfile(out_path), row.names = FALSE)
  }, error = function(e) {
    status <<- "failed"
    message <<- conditionMessage(e)
  })
  rows[[length(rows) + 1]] <- data.frame(
    dataset_id = "GSE154600",
    author_sample_id = sample_id,
    source_file = path,
    coldata_file = out_path,
    n_cells = n_cells,
    n_genes = n_genes,
    status = status,
    message = message,
    stringsAsFactors = FALSE
  )
}
audit <- do.call(rbind, rows)
utils::write.csv(audit, file.path(audit_dir, "GSE154600_sce_extract_audit.csv"), row.names = FALSE)
if (!nrow(audit) || any(audit$status == "failed")) {
  quit(status = 1)
}
