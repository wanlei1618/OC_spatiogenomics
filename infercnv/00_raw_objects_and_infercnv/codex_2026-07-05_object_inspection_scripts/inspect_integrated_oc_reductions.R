suppressPackageStartupMessages(library(Matrix))

setClass("KeyMixin", contains = "VIRTUAL", slots = list(key = "character"))
setClass("Assay", contains = "KeyMixin",
         slots = c(counts = "AnyMatrix", data = "AnyMatrix", scale.data = "matrix",
                   assay.orig = "OptionalCharacter", var.features = "vector",
                   meta.features = "data.frame", misc = "OptionalList"))
setClass("Seurat",
         slots = c(assays = "list", meta.data = "data.frame", active.assay = "character",
                   active.ident = "factor", graphs = "list", neighbors = "list",
                   reductions = "list", images = "list", project.name = "character",
                   misc = "list", version = "package_version", commands = "list",
                   tools = "list"))
obj <- readRDS("D:/OC_spatiogenomics/infercnv/integrated_oc.RData")
cat("reductions:", names(obj@reductions), "\n")
for (nm in names(obj@reductions)) {
  r <- obj@reductions[[nm]]
  cat("reduction", nm, "class", class(r), "\n")
  emb <- tryCatch(r@cell.embeddings, error = function(e) NULL)
  if (!is.null(emb)) {
    cat(" embeddings:", paste(dim(emb), collapse = "x"), "\n")
    print(head(colnames(emb)))
  }
}
