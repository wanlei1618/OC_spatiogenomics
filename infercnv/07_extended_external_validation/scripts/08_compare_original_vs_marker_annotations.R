#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE)
workflow <- if (grepl("infercnv/07_extended_external_validation", root, fixed = TRUE)) root else file.path(root, "infercnv/07_extended_external_validation")
tables <- file.path(workflow, "tables")
conf <- file.path(tables, "original_vs_marker_confusion.csv")
if (!file.exists(conf)) {
  write.csv(data.frame(note = "No marker-score v1 comparison table available."), conf, row.names = FALSE)
}
