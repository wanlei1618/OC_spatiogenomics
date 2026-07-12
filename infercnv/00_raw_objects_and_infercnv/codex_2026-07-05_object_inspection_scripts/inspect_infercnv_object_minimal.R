obj_file <- "D:/OC_spatiogenomics/infercnv/infercnv_Other_vs_Immune_subcluster/01_incoming_data.infercnv_obj"
obj <- readRDS(obj_file)
cat("class:", class(obj), "\n")
cat("type:", typeof(obj), "\n")
print(attributes(obj))
