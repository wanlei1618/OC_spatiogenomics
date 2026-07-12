#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args) >= 1) args[[1]] else "D:/OC_spatiogenomics/infercnv/external_cell_annotations"
log_dir <- file.path(data_root, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  "GSE151214 is classified as normal_fallopian_tube_reference.",
  "It is excluded from tumor-ecosystem LR meta-analysis.",
  "Normal epithelial CD44/ITGB1 background should be interpreted as reference context only."
), file.path(log_dir, "GSE151214_normal_reference_note.txt"))
