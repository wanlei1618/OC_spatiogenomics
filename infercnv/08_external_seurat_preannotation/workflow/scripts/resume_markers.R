#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1]]
}

config_path <- arg_value("--config", "config/five_external_datasets.yaml")
dataset_id <- arg_value("--dataset")
if (is.null(dataset_id) || !nzchar(dataset_id)) stop("--dataset is required")

required <- c("yaml", "Seurat", "dplyr", "ggplot2", "jsonlite")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                           FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

cfg <- yaml::read_yaml(config_path)
set.seed(as.integer(cfg$project$random_seed %||% 20260713))
root <- file.path(cfg$project$output_root, dataset_id)
object_path <- file.path(root, "objects", paste0(dataset_id, "_preannotation.rds"))
status_path <- file.path(root, "logs", "run_status.json")
if (!file.exists(object_path)) stop("Checkpoint not found: ", object_path)

write_csv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE)
}
write_status <- function(x) {
  jsonlite::write_json(x, status_path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null")
}

message(dataset_id, ": loading preannotation checkpoint")
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
Idents(obj) <- "seurat_clusters"

status <- list(dataset_id = dataset_id,
               status = "RESUMED_FROM_PREANNOTATION_CHECKPOINT",
               resumed_at = as.character(Sys.time()))
write_status(status)

message(dataset_id, ": running FindAllMarkers for ", length(levels(Idents(obj))), " clusters")
markers <- FindAllMarkers(
  obj, assay = "RNA",
  only.pos = isTRUE(cfg$markers$only_pos),
  test.use = cfg$markers$test_use,
  min.pct = cfg$markers$min_pct,
  logfc.threshold = cfg$markers$logfc_threshold,
  return.thresh = cfg$markers$adjusted_p_threshold,
  verbose = TRUE
)
if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
write_csv_gz(markers, file.path(root, "03_markers", "all_cluster_markers.csv.gz"))

for (n in unlist(cfg$markers$export_top_n)) {
  top <- markers |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) |>
    slice_head(n = n) |>
    ungroup()
  write.csv(top, file.path(root, "03_markers",
                           paste0("top", n, "_markers_per_cluster.csv")),
            row.names = FALSE)
}

message(dataset_id, ": exporting average expression and marker plots")
avg <- AverageExpression(obj, assays = "RNA", slot = "data",
                         group.by = "seurat_clusters", verbose = FALSE)$RNA
avg <- data.frame(gene = rownames(avg), as.data.frame(avg), check.names = FALSE)
write_csv_gz(avg, file.path(root, "03_markers",
                            "cluster_average_expression.csv.gz"))

panel <- c("EPCAM","KRT8","KRT18","KRT19","MSLN","WFDC2",
           "PTPRC","CD3D","CD3E","TRBC1","NKG7","GNLY",
           "CD79A","MS4A1","MZB1","JCHAIN",
           "LYZ","LST1","TYROBP","FCER1G","C1QA","C1QB","C1QC","SPP1",
           "COL1A1","COL1A2","DCN","COL3A1",
           "PECAM1","VWF","KDR","MKI67","TOP2A")
panel <- panel[panel %in% rownames(obj)]
if (length(panel)) {
  p <- DotPlot(obj, features = panel, group.by = "seurat_clusters") + RotatedAxis()
  ggsave(file.path(root, "03_markers", "broad_marker_dotplot.pdf"),
         plot = p, width = 16, height = 8)
}

top10 <- markers |>
  group_by(cluster) |>
  arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) |>
  slice_head(n = 10) |>
  pull(gene) |>
  unique()
top10 <- top10[top10 %in% rownames(obj)]
if (length(top10)) {
  heatmap_cells <- WhichCells(obj, downsample = 200,
                              seed = cfg$project$random_seed)
  heatmap_obj <- subset(obj, cells = heatmap_cells)
  p <- DoHeatmap(heatmap_obj, features = top10,
                 group.by = "seurat_clusters", raster = TRUE) + NoLegend()
  ggsave(file.path(root, "03_markers", "top_marker_heatmap.pdf"),
         plot = p, width = 14, height = 12)
  rm(heatmap_obj)
}

cc <- read.csv(file.path(root, "02_clustering", "cluster_cell_counts.csv"),
               check.names = FALSE)
cc$seurat_cluster <- as.character(cc$seurat_cluster)
top_text <- markers |>
  mutate(cluster = as.character(cluster)) |>
  group_by(cluster) |>
  arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) |>
  slice_head(n = 20) |>
  summarise(top_markers = paste(gene, collapse = ";"), .groups = "drop")
manual <- cc |>
  left_join(top_text, by = c("seurat_cluster" = "cluster")) |>
  mutate(dataset_id = dataset_id, cell_type_manual = "",
         cell_subtype_manual = "", confidence = "", notes = "") |>
  select(dataset_id, seurat_cluster, n_cells, top_markers,
         cell_type_manual, cell_subtype_manual, confidence, notes)
write.csv(manual,
          file.path(root, "04_manual_annotation", "manual_annotation_template.csv"),
          row.names = FALSE)

capture.output(sessionInfo(), file = file.path(root, "logs", "sessionInfo.txt"))
retention <- read.csv(file.path(root, "01_qc", "qc_cell_retention.csv"))
status$status <- "PREANNOTATION_COMPLETE_WAITING_FOR_MANUAL_CELLTYPE"
status$finished_at <- as.character(Sys.time())
status$n_input <- retention$n_input[[1]]
status$n_after_qc <- retention$n_after_doublet_filter[[1]]
status$primary_resolution <- as.numeric(cfg$analysis$primary_resolution)
status$n_clusters <- length(unique(obj$seurat_clusters))
status$n_marker_rows <- nrow(markers)
status$doublet_status <- retention$doublet_status[[1]]
write_status(status)
message(dataset_id, ": marker resume complete")
