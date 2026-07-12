suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

out_dir <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"

cell_map <- read.csv(
  file.path(out_dir, "infercnv_cell_to_subclone_k5.csv"),
  stringsAsFactors = FALSE
)
regions <- read.csv(
  file.path(out_dir, "infercnv_subclone_cnv_regions_k5.csv"),
  stringsAsFactors = FALSE
)

regions$region_id <- paste(regions$chr, regions$start, regions$end, sep = ":")
regions$region_label <- paste(regions$chr, regions$start, regions$end, sep = "_")

obs_groups <- unique(cell_map$cell_group_name)
region_ids <- unique(regions$region_id)

group_region_mat <- matrix(
  0,
  nrow = length(region_ids),
  ncol = length(obs_groups),
  dimnames = list(region_ids, obs_groups)
)

for (i in seq_len(nrow(regions))) {
  group_region_mat[regions$region_id[i], regions$cell_group_name[i]] <- regions$state[i] - 3
}

cnv_matrix_obs <- group_region_mat[, cell_map$cell_group_name, drop = FALSE]
colnames(cnv_matrix_obs) <- cell_map$cell

row_info <- unique(regions[, c("region_id", "chr", "start", "end", "region_label")])
row_info <- row_info[match(rownames(cnv_matrix_obs), row_info$region_id), ]
rownames(cnv_matrix_obs) <- row_info$region_label

subclone <- factor(cell_map$cnv_subclone)
sample_id <- sub("_.*$", "", cell_map$cell)

subclone_cols <- c(
  Subclone_01 = "#4C78A8",
  Subclone_02 = "#F58518",
  Subclone_03 = "#54A24B",
  Subclone_04 = "#E45756",
  Subclone_05 = "#B279A2",
  Subclone_06 = "#72B7B2",
  Subclone_07 = "#FF9DA6"
)
subclone_cols <- subclone_cols[levels(subclone)]

sample_cols <- setNames(
  circlize::rand_color(length(unique(sample_id)), luminosity = "bright"),
  sort(unique(sample_id))
)

ha <- HeatmapAnnotation(
  Subclone = subclone,
  Sample = sample_id,
  col = list(Subclone = subclone_cols, Sample = sample_cols),
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(3, "mm")
)

col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2, 3),
  c("#2166AC", "#67A9CF", "white", "#F4A582", "#B2182B", "#67001F")
)

chr_num <- suppressWarnings(as.integer(sub("^chr", "", row_info$chr)))
chr_num[is.na(chr_num)] <- match(row_info$chr[is.na(chr_num)], unique(row_info$chr[is.na(chr_num)])) + 100
row_order <- order(chr_num, row_info$start)
column_order <- order(subclone, sample_id, cell_map$cell_group_name, cell_map$cell)

cnv_matrix_obs <- cnv_matrix_obs[row_order, column_order, drop = FALSE]
row_info <- row_info[row_order, , drop = FALSE]
subclone <- subclone[column_order]
sample_id <- sample_id[column_order]

ha <- HeatmapAnnotation(
  Subclone = subclone,
  Sample = sample_id,
  col = list(Subclone = subclone_cols, Sample = sample_cols),
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(3, "mm")
)

ht <- Heatmap(
  cnv_matrix_obs,
  name = "CNV",
  col = col_fun,
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = row_info$chr,
  column_split = subclone,
  column_order = seq_len(ncol(cnv_matrix_obs)),
  row_title = "CNV regions",
  column_title = "Cells grouped by inferCNV subclone",
  use_raster = TRUE,
  raster_quality = 2,
  heatmap_legend_param = list(
    title = "HMM state\n(state - 3)",
    at = c(-2, -1, 0, 1, 2, 3),
    labels = c("deep loss", "loss", "neutral/none", "gain", "amp", "high amp")
  )
)

pdf(file.path(out_dir, "infercnv_subclone_heatmap.pdf"), width = 13, height = 8)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(file.path(out_dir, "infercnv_subclone_heatmap.png"), width = 2600, height = 1600, res = 200)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

write.csv(
  data.frame(
    cnv_subclone = names(sort(table(subclone), decreasing = TRUE)),
    n_cells = as.integer(sort(table(subclone), decreasing = TRUE))
  ),
  file.path(out_dir, "infercnv_subclone_heatmap_cell_counts.csv"),
  row.names = FALSE
)

message("Done: infercnv_subclone_heatmap.pdf and infercnv_subclone_heatmap.png")
