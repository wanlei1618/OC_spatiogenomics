library(SeuratObject)
obj <- readRDS("D:/OC_spatiogenomics/infercnv/integratedocTcells.RData")
cat("class:", class(obj), "\n")
cat("nmeta:", nrow(obj@meta.data), "\n")
cat("assays:", names(obj@assays), "\n")
rna <- obj@assays$RNA
cat("RNA class:", class(rna), "\n")
cat("RNA slots:", slotNames(rna), "\n")
cat("layer names:", names(rna@layers), "\n")
for (nm in names(rna@layers)) {
  cat("layer", nm, "dim", paste(dim(rna@layers[[nm]]), collapse = "x"), "\n")
}
cat("features length:", length(rna@features), "\n")
cat("cells length:", length(rna@cells), "\n")
cat("meta columns:\n")
print(colnames(obj@meta.data))
print(table(obj@meta.data$cell_type, useNA = "ifany"))
