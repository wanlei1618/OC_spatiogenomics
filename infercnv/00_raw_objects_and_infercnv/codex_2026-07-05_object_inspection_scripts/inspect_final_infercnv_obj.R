obj_file <- "D:/OC_spatiogenomics/infercnv/infercnv_Other_vs_Immune_subcluster/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.infercnv_obj"
setClass(
  "infercnv",
  slots = c(
    expr.data = "matrix",
    count.data = "matrix",
    gene_order = "data.frame",
    reference_grouped_cell_indices = "list",
    observation_grouped_cell_indices = "list",
    tumor_subclusters = "list",
    options = "list",
    .hspike = "ANY"
  )
)
obj <- readRDS(obj_file)
cat("class:", class(obj), "\n")
cat("slots:", slotNames(obj), "\n")
for (s in slotNames(obj)) {
  x <- tryCatch(slot(obj, s), error = function(e) NULL)
  cat("\nSLOT", s, "class", class(x), "\n")
  if (is.matrix(x) || inherits(x, "Matrix")) cat("dim", paste(dim(x), collapse = "x"), "\n")
  else if (is.data.frame(x)) { cat("dim", paste(dim(x), collapse = "x"), "\n"); print(head(x)) }
  else if (is.list(x)) cat("length", length(x), "names", paste(head(names(x)), collapse = ","), "\n")
  else print(utils::head(x))
}
