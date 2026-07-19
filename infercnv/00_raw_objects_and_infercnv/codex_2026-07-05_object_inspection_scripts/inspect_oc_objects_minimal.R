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

describe_obj <- function(path) {
  cat("\n====", path, "====\n")
  obj <- readRDS(path)
  cat("class:", class(obj), "\n")
  cat("meta dim:", nrow(obj@meta.data), ncol(obj@meta.data), "\n")
  print(colnames(obj@meta.data))
  for (nm in intersect(c("cell_type", "sample_type", "orig.ident", "seurat_clusters"), colnames(obj@meta.data))) {
    cat("\n", nm, "\n")
    print(head(sort(table(obj@meta.data[[nm]], useNA = "ifany"), decreasing = TRUE), 30))
  }
  cat("assays:", names(obj@assays), "active:", obj@active.assay, "\n")
  for (assay in names(obj@assays)) {
    a <- obj@assays[[assay]]
    cat("assay", assay, "class", class(a), "\n")
    layers <- tryCatch(a@layers, error = function(e) NULL)
    if (!is.null(layers)) {
      cat("layers:", names(layers), "\n")
      for (ln in names(layers)) cat(" ", ln, paste(dim(layers[[ln]]), collapse = "x"), class(layers[[ln]]), "\n")
    } else {
      counts <- tryCatch(a@counts, error = function(e) NULL)
      data <- tryCatch(a@data, error = function(e) NULL)
      scale.data <- tryCatch(a@scale.data, error = function(e) NULL)
      if (!is.null(counts)) cat("counts:", paste(dim(counts), collapse = "x"), "\n")
      if (!is.null(data)) cat("data:", paste(dim(data), collapse = "x"), "\n")
      if (!is.null(scale.data)) cat("scale.data:", paste(dim(scale.data), collapse = "x"), "\n")
    }
  }
  invisible(NULL)
}

files <- c(
  "D:/OC_spatiogenomics/infercnv/integrated_oc.RData",
  "D:/OC_spatiogenomics/infercnv/integratedocTcells.RData",
  "D:/OC_spatiogenomics/infercnv/integratedocMyecells.RData"
)
for (f in files[file.exists(files)]) describe_obj(f)
