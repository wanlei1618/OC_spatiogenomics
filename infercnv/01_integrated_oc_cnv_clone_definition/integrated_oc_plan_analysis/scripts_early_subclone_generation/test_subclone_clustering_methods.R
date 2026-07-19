infer_dir <- "D:/OC_spatiogenomics/infercnv/infercnv_Other_vs_Immune_subcluster"
cg <- read.delim(file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"), stringsAsFactors = FALSE)
regions <- read.delim(file.path(infer_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat"), stringsAsFactors = FALSE)
obs_groups <- sort(unique(regions$cell_group_name[grepl("^observation\\.", regions$cell_group_name)]))
regions_obs <- regions[regions$cell_group_name %in% obs_groups, ]
key <- with(regions_obs, paste(chr, start, end, sep = ":"))
all_regions <- sort(unique(key))
mat <- matrix(0, nrow = length(obs_groups), ncol = length(all_regions), dimnames = list(obs_groups, all_regions))
signed_state <- function(state) ifelse(state < 3, -1, ifelse(state > 3, 1, 0))
for (i in seq_len(nrow(regions_obs))) {
  mat[regions_obs$cell_group_name[i], key[i]] <- signed_state(regions_obs$state[i])
}
group_sizes <- table(cg$cell_group_name)
obs_sizes <- as.integer(group_sizes[rownames(mat)])
obs_sizes[is.na(obs_sizes)] <- 0

mat2 <- mat[, colSums(abs(mat)) > 1, drop = FALSE]
pca <- prcomp(mat2, center = TRUE, scale. = FALSE)
pcs <- pca$x[, seq_len(min(10, ncol(pca$x))), drop = FALSE]
set.seed(1)
for (k in 4:7) {
  km <- kmeans(pcs, centers = k, nstart = 100, iter.max = 100)
  cat("\nkmeans k=", k, "subclusters\n")
  print(table(km$cluster))
  cat("cells\n")
  print(tapply(obs_sizes, km$cluster, sum))
}
cat("\nhclust average binary k sizes\n")
d <- dist(pcs)
hc <- hclust(d, method = "average")
for (k in 4:7) {
  cl <- cutree(hc, k = k)
  cat("\nhclust k=", k, "subclusters\n")
  print(table(cl))
  cat("cells\n")
  print(tapply(obs_sizes, cl, sum))
}
