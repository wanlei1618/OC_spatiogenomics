suppressPackageStartupMessages(library(cluster))

infer_dir <- "D:/OC_spatiogenomics/infercnv/infercnv_Other_vs_Immune_subcluster"
out_dir <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cell_group_file <- file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings")
regions_file <- file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat")

cg <- read.delim(cell_group_file, stringsAsFactors = FALSE)
regions <- read.delim(regions_file, stringsAsFactors = FALSE)

obs_groups <- sort(unique(regions$cell_group_name[grepl("^observation\\.", regions$cell_group_name)]))
regions_obs <- regions[regions$cell_group_name %in% obs_groups, , drop = FALSE]

region_key <- with(regions_obs, paste(chr, start, end, sep = ":"))
all_regions <- sort(unique(region_key))
mat <- matrix(0, nrow = length(obs_groups), ncol = length(all_regions),
              dimnames = list(obs_groups, all_regions))

signed_state <- function(state) {
  ifelse(state < 3, state - 3, ifelse(state > 3, state - 3, 0))
}
for (i in seq_len(nrow(regions_obs))) {
  g <- regions_obs$cell_group_name[i]
  key <- paste(regions_obs$chr[i], regions_obs$start[i], regions_obs$end[i], sep = ":")
  mat[g, key] <- signed_state(regions_obs$state[i])
}

d <- dist(mat, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

k_values <- 4:7
sil <- sapply(k_values, function(k) {
  cl <- cutree(hc, k = k)
  mean(silhouette(cl, d)[, "sil_width"])
})
best_k <- k_values[which.max(sil)]
cl <- cutree(hc, k = best_k)

sizes <- sort(table(cl), decreasing = TRUE)
rank_map <- setNames(sprintf("Subclone_%02d", seq_along(sizes)), names(sizes))
subclone <- unname(rank_map[as.character(cl)])

group_map <- data.frame(
  infercnv_subcluster = names(cl),
  cnv_subclone = subclone,
  n_cells = as.integer(table(factor(cg$cell_group_name, levels = names(cl)))),
  stringsAsFactors = FALSE
)
group_map <- group_map[order(group_map$cnv_subclone, group_map$infercnv_subcluster), ]

cell_map <- cg[grepl("^observation\\.", cg$cell_group_name), , drop = FALSE]
cell_map$cnv_subclone <- group_map$cnv_subclone[match(cell_map$cell_group_name, group_map$infercnv_subcluster)]
cell_map$cell_for_seurat_v5_style <- sub("^([^_]+)_(.+)$", "\\1_\\1_\\2", cell_map$cell)
cell_map <- cell_map[, c("cell", "cell_for_seurat_v5_style", "cell_group_name", "cnv_subclone")]

regions_obs$cnv_subclone <- group_map$cnv_subclone[match(regions_obs$cell_group_name, group_map$infercnv_subcluster)]
subclone_summary <- aggregate(
  list(n_regions = regions_obs$cnv_name),
  by = list(cnv_subclone = regions_obs$cnv_subclone, state = regions_obs$state),
  FUN = length
)
subclone_summary <- subclone_summary[order(subclone_summary$cnv_subclone, subclone_summary$state), ]

top_regions <- regions_obs
top_regions$event <- ifelse(top_regions$state < 3, "loss", "gain")
top_regions$region_label <- paste(top_regions$chr, top_regions$start, top_regions$end, sep = ":")
top_regions <- top_regions[order(top_regions$cnv_subclone, top_regions$state, top_regions$chr, top_regions$start), ]

write.csv(group_map, file.path(out_dir, "infercnv_subcluster_to_4_7_subclone.csv"), row.names = FALSE)
write.csv(cell_map, file.path(out_dir, "infercnv_cell_to_subclone.csv"), row.names = FALSE)
write.csv(subclone_summary, file.path(out_dir, "infercnv_subclone_cnv_state_summary.csv"), row.names = FALSE)
write.csv(top_regions, file.path(out_dir, "infercnv_subclone_cnv_regions_long.csv"), row.names = FALSE)

png(file.path(out_dir, "infercnv_subclone_hclust.png"), width = 1800, height = 1000, res = 150)
plot(hc, labels = FALSE, main = sprintf("inferCNV observation subclusters reclustered, k=%d", best_k),
     xlab = "inferCNV observation subclusters", sub = "")
rect.hclust(hc, k = best_k, border = 2:(best_k + 1))
dev.off()

cat("best_k", best_k, "\n")
print(sil)
cat("observation subclusters", length(obs_groups), "\n")
cat("observation cells", nrow(cell_map), "\n")
print(table(cell_map$cnv_subclone))
