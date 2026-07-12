#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args) >= 1) args[[1]] else "D:/OC_spatiogenomics/infercnv/external_cell_annotations"
log_dir <- file.path(data_root, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  "GSE147082 is marked secondary_reannotation_only.",
  "This lightweight implementation preserves expression barcodes when available and labels confidence as secondary_placeholder.",
  "It is restricted to sensitivity analysis and is not used as author-original primary evidence."
), file.path(log_dir, "GSE147082_secondary_annotation_note.txt"))
