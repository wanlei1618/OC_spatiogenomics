suppressPackageStartupMessages({
  library(Matrix)
})

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
  "Assay",
  contains = "KeyMixin",
  slots = c(
    counts = "AnyMatrix",
    data = "AnyMatrix",
    scale.data = "matrix",
    assay.orig = "OptionalCharacter",
    var.features = "vector",
    meta.features = "data.frame",
    misc = "OptionalList"
  )
)
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

obj <- readRDS("D:/OC_spatiogenomics/infercnv/integrated_oc.RData")
cat("class:", class(obj), "\n")
cat("assays:", names(obj@assays), "\n")
cat("active.assay:", obj@active.assay, "\n")
cat("meta dim:", dim(obj@meta.data), "\n")
cat("meta cols:", paste(colnames(obj@meta.data), collapse=", "), "\n")
print(table(obj@meta.data$cell_type, useNA="ifany"))
for (assay_name in names(obj@assays)) {
  assay <- obj@assays[[assay_name]]
  cat("\nassay", assay_name, "class", class(assay), "\n")
  if ("Assay" %in% class(assay)) {
    x <- slot(assay, "data")
    cat("data dim", dim(x), "class", class(x), "\n")
    cat("first genes", paste(head(rownames(x)), collapse=", "), "\n")
    x <- slot(assay, "counts")
    cat("counts dim", dim(x), "class", class(x), "\n")
  }
  if ("Assay5" %in% class(assay)) {
    cat("layers", paste(names(assay@layers), collapse=", "), "\n")
    for (nm in names(assay@layers)) {
      cat("layer", nm, "dim", dim(assay@layers[[nm]]), "class", class(assay@layers[[nm]]), "\n")
    }
  }
}
