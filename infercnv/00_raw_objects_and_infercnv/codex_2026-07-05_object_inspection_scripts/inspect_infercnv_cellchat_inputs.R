cat("R version:", R.version.string, "\n")
cat(".libPaths:\n")
print(.libPaths())

infer_dir <- "D:/OC_spatiogenomics/infercnv/infercnv_Other_vs_Immune_subcluster"
cg_file <- file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings")
regions_file <- file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat")
obj_file <- "D:/OC_spatiogenomics/infercnv/integratedocTcells.RData"

cat("\nInstalled relevant packages:\n")
ip <- installed.packages()
pat <- "^(Seurat|SeuratObject|CellChat|NMF|future|patchwork|ggplot2|igraph|dplyr|Matrix)$"
print(ip[grep(pat, rownames(ip)), c("Package", "Version"), drop = FALSE])

cat("\ninferCNV cell_groupings:\n")
cg <- read.delim(cg_file, stringsAsFactors = FALSE)
print(dim(cg))
print(head(cg))
print(length(unique(cg$cell_group_name)))
print(head(sort(table(cg$cell_group_name), decreasing = TRUE), 25))

cat("\ninferCNV regions:\n")
regions <- read.delim(regions_file, stringsAsFactors = FALSE)
print(dim(regions))
print(head(regions))
print(table(regions$state, useNA = "ifany"))

cat("\nObject file first bytes:\n")
con <- file(obj_file, "rb")
print(readBin(con, "raw", n = 32))
close(con)

cat("\nTrying readRDS:\n")
obj <- tryCatch(readRDS(obj_file), error = function(e) e)
print(class(obj))
if (inherits(obj, "error")) {
  print(obj$message)
} else {
  print(class(obj))
  if (inherits(obj, "Seurat")) {
    print(dim(obj))
    print(head(colnames(obj@meta.data), 100))
    if ("cell_type" %in% colnames(obj@meta.data)) {
      print(table(obj@meta.data$cell_type, useNA = "ifany"))
    }
    print(head(rownames(obj@meta.data)))
  } else {
    print(str(obj, max.level = 2))
  }
}
