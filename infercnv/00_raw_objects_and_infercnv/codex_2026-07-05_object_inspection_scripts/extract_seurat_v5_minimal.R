suppressPackageStartupMessages(library(Matrix))

setClass("KeyMixin", contains = "VIRTUAL", slots = list(key = "character"))
setClass("LogMap", contains = "matrix")
setClass(
  "StdAssay",
  contains = c("VIRTUAL", "KeyMixin"),
  slots = c(
    layers = "list",
    cells = "LogMap",
    features = "LogMap",
    default = "integer",
    assay.orig = "character",
    meta.data = "data.frame",
    misc = "list"
  )
)
setClass("Assay5", contains = "StdAssay")
setClass(
  "Seurat",
  slots = c(
    assays = "list",
    meta.data = "data.frame",
    active.assay = "character",
    active.ident = "factor",
    graphs = "list",
    neighbors = "list",
    reductions = "list",
    images = "list",
    project.name = "character",
    misc = "list",
    version = "package_version",
    commands = "list",
    tools = "list"
  )
)

obj <- readRDS("D:/OC_spatiogenomics/infercnv/integratedocTcells.RData")
cat("object:", class(obj), "\n")
cat("metadata:", nrow(obj@meta.data), "cells x", ncol(obj@meta.data), "columns\n")
print(colnames(obj@meta.data))
print(table(obj@meta.data$cell_type, useNA = "ifany"))
rna <- obj@assays$RNA
cat("RNA class:", class(rna), "\n")
cat("layers:", names(rna@layers), "\n")
for (nm in names(rna@layers)) {
  cat(nm, paste(dim(rna@layers[[nm]]), collapse = "x"), class(rna@layers[[nm]]), "\n")
}
cat("features map:", dim(rna@features), "\n")
cat("cells map:", dim(rna@cells), "\n")
cat("first cells:", paste(head(rownames(obj@meta.data), 3), collapse = ", "), "\n")
