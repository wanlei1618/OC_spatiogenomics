suppressPackageStartupMessages(library(Matrix))

setClass("KeyMixin", contains = "VIRTUAL", slots = list(key = "character"))
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
meta <- obj@meta.data
b <- meta$cell_type == "B cells"
cat("B cells:", sum(b), "\n")
print(sort(table(meta$seurat_clusters[b]), decreasing = TRUE))

markers <- c(
  "TCL1A", "TXNIP", "FCER2", "MS4A1", "IGHG1", "IGHG2", "IGHG3", "IGHG4",
  "MZB1", "JCHAIN", "XBP1", "HSPA1B", "HSP90AA1", "DNAJB1", "CD83",
  "TNFRSF13B", "CD70", "FCRLA", "LAPTM5", "BTG1"
)
expr <- obj@assays$RNA@data
markers <- intersect(markers, rownames(expr))
clusters <- sort(unique(meta$seurat_clusters[b]))
avg <- matrix(NA_real_, nrow = length(markers), ncol = length(clusters), dimnames = list(markers, clusters))
pct <- avg
for (cl in clusters) {
  cells <- rownames(meta)[b & meta$seurat_clusters == cl]
  m <- expr[markers, cells, drop = FALSE]
  avg[, as.character(cl)] <- Matrix::rowMeans(m)
  pct[, as.character(cl)] <- Matrix::rowMeans(m > 0)
}
out_dir <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"
write.csv(avg, file.path(out_dir, "integrated_oc_Bcells_marker_avg_by_full_cluster.csv"))
write.csv(pct, file.path(out_dir, "integrated_oc_Bcells_marker_pct_by_full_cluster.csv"))
print(round(avg, 2))
