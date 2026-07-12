suppressPackageStartupMessages(library(Matrix))

out_dir <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
infer_dir <- "D:/OC_spatiogenomics/infercnv/infercnv_Other_vs_Immune_subcluster"

cg <- read.delim(file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"), stringsAsFactors = FALSE)
regions <- read.delim(file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat"), stringsAsFactors = FALSE)

obs_groups <- sort(unique(regions$cell_group_name[grepl("^observation\\.", regions$cell_group_name)]))
regions_obs <- regions[regions$cell_group_name %in% obs_groups, , drop = FALSE]
region_key <- with(regions_obs, paste(chr, start, end, sep = ":"))
all_regions <- sort(unique(region_key))
mat <- matrix(0, nrow = length(obs_groups), ncol = length(all_regions), dimnames = list(obs_groups, all_regions))
signed_state <- function(state) ifelse(state < 3, -1, ifelse(state > 3, 1, 0))
for (i in seq_len(nrow(regions_obs))) {
  mat[regions_obs$cell_group_name[i], region_key[i]] <- signed_state(regions_obs$state[i])
}
mat <- mat[, colSums(abs(mat)) > 1, drop = FALSE]
pca <- prcomp(mat, center = TRUE, scale. = FALSE)
pcs <- pca$x[, seq_len(min(10, ncol(pca$x))), drop = FALSE]
set.seed(1)
k <- 5
km <- kmeans(pcs, centers = k, nstart = 100, iter.max = 100)

group_sizes <- table(cg$cell_group_name)
cell_counts_by_cluster <- tapply(as.integer(group_sizes[rownames(mat)]), km$cluster, sum)
ranked_clusters <- names(sort(cell_counts_by_cluster, decreasing = TRUE))
label_map <- setNames(sprintf("Subclone_%02d", seq_along(ranked_clusters)), ranked_clusters)
subclone <- unname(label_map[as.character(km$cluster)])

group_map <- data.frame(
  infercnv_subcluster = rownames(mat),
  kmeans_cluster = km$cluster,
  cnv_subclone = subclone,
  n_cells = as.integer(group_sizes[rownames(mat)]),
  stringsAsFactors = FALSE
)
group_map <- group_map[order(group_map$cnv_subclone, group_map$infercnv_subcluster), ]

cell_map <- cg[grepl("^observation\\.", cg$cell_group_name), , drop = FALSE]
cell_map$cnv_subclone <- group_map$cnv_subclone[match(cell_map$cell_group_name, group_map$infercnv_subcluster)]
cell_map$cell_for_integratedocTcells_style <- sub("^([^_]+)_(.+)$", "\\1_\\1_\\2", cell_map$cell)
cell_map <- cell_map[, c("cell", "cell_for_integratedocTcells_style", "cell_group_name", "cnv_subclone")]

regions_obs$cnv_subclone <- group_map$cnv_subclone[match(regions_obs$cell_group_name, group_map$infercnv_subcluster)]
regions_obs$event <- ifelse(regions_obs$state < 3, "loss", "gain")

state_summary <- aggregate(
  list(n_regions = regions_obs$cnv_name),
  by = list(cnv_subclone = regions_obs$cnv_subclone, state = regions_obs$state, event = regions_obs$event),
  FUN = length
)
state_summary <- state_summary[order(state_summary$cnv_subclone, state_summary$state), ]

write.csv(group_map, file.path(out_dir, "infercnv_subcluster_to_subclone_k5.csv"), row.names = FALSE)
write.csv(cell_map, file.path(out_dir, "infercnv_cell_to_subclone_k5.csv"), row.names = FALSE)
write.csv(regions_obs, file.path(out_dir, "infercnv_subclone_cnv_regions_k5.csv"), row.names = FALSE)
write.csv(state_summary, file.path(out_dir, "infercnv_subclone_cnv_state_summary_k5.csv"), row.names = FALSE)

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

tobj <- readRDS("D:/OC_spatiogenomics/infercnv/integratedocTcells.RData")
tmeta <- tobj@meta.data
tmeta$cell_integratedocTcells <- rownames(tmeta)
tmeta$cell_infercnv_style <- sub("^[^_]+_", "", rownames(tmeta))
tmeta$infercnv_cell_group_name <- cg$cell_group_name[match(tmeta$cell_infercnv_style, cg$cell)]
tmeta$infercnv_group_type <- ifelse(grepl("^observation\\.", tmeta$infercnv_cell_group_name), "observation",
                                    ifelse(grepl("^reference\\.", tmeta$infercnv_cell_group_name), "reference", NA))
tmeta$cnv_subclone_if_observation <- cell_map$cnv_subclone[match(tmeta$cell_infercnv_style, cell_map$cell)]
write.csv(tmeta, file.path(out_dir, "integratedocTcells_metadata_with_infercnv_match.csv"), row.names = FALSE)

sink(file.path(out_dir, "analysis_notes_subclone_cellchat_status.txt"))
cat("inferCNV subclone reclustering\n")
cat("Method: HMM predicted CNV regions -> signed event matrix (-1 loss, 0 absent/neutral, +1 gain), PCA, k-means k=5.\n")
cat("Observation subclusters:", length(obs_groups), "\n")
cat("Observation cells:", nrow(cell_map), "\n\n")
cat("Subclone cell counts:\n")
print(sort(table(cell_map$cnv_subclone), decreasing = TRUE))
cat("\nT cell object metadata rows:", nrow(tmeta), "\n")
cat("T cell_type counts:\n")
print(table(tmeta$cell_type, useNA = "ifany"))
cat("\nT cells matched to inferCNV groups:\n")
print(table(tmeta$infercnv_group_type, useNA = "ifany"))
cat("\nNote: integratedocTcells contains T/NK cells only. It has no tumor observation cells/expression, so a full CellChat run between tumor subclones and T cell_type requires the full OO/OC Seurat object containing malignant cells, or a separate tumor expression matrix with matching barcodes.\n")
sink()

cat("done\n")
print(sort(table(cell_map$cnv_subclone), decreasing = TRUE))
print(table(tmeta$infercnv_group_type, useNA = "ifany"))
