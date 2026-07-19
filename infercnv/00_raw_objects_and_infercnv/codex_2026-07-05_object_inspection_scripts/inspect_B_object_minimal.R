suppressPackageStartupMessages(library(Matrix))
setClass("KeyMixin", contains = "VIRTUAL", slots = list(key = "character"))
setClass("LogMap", contains = "matrix")
setClass("StdAssay", contains = c("VIRTUAL", "KeyMixin"),
         slots = c(layers = "list", cells = "LogMap", features = "LogMap",
                   default = "integer", assay.orig = "character",
                   meta.data = "data.frame", misc = "list"))
setClass("Assay5", contains = "StdAssay")
setClass("Seurat", slots = c(assays = "list", meta.data = "data.frame",
                             active.assay = "character", active.ident = "factor",
                             graphs = "list", neighbors = "list", reductions = "list",
                             images = "list", project.name = "character", misc = "list",
                             version = "package_version", commands = "list", tools = "list"))
obj <- readRDS("D:/OC_spatiogenomics/infercnv/integratedocBcells.RData")
cat("meta dim:", nrow(obj@meta.data), ncol(obj@meta.data), "\n")
print(colnames(obj@meta.data))
print(table(obj@meta.data$cell_type, useNA = "ifany"))
cat("head cells:\n")
print(head(rownames(obj@meta.data)))
