#!/usr/bin/env Rscript

# Transfer single-cell reference labels and CNV-subclone identities to spatial spots.
#
# Usage:
#   Rscript 04_reference_mapping_to_cnv_niches.R \
#     D:/OC_spatiogenomics/spatial_data \
#     D:/OC_spatiogenomics/infercnv/integrated_oc.RData \
#     D:/OC_spatiogenomics/infercnv/integrated_oc_plan_analysis/tables/integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv
#
# The script maps broad reference states rather than claiming one cell per Visium
# spot. Prediction scores should be interpreted as spot composition/state evidence.

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Provide: spatial_root, integrated_oc.RData, mapped_metadata.csv")
}

root <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
reference_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)

processed_dir <- file.path(root, "processed")
result_dir <- file.path(root, "results", "reference_mapping")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

reference <- tryCatch({
  loaded_names <- load(reference_path)
  seurat_names <- loaded_names[vapply(loaded_names, function(x) inherits(get(x), "Seurat"), logical(1))]
  if (length(seurat_names) == 0) stop("No Seurat object found in RData workspace")
  get(seurat_names[[1]])
}, error = function(e) {
  obj <- readRDS(reference_path)
  if (!inherits(obj, "Seurat")) stop("Reference file is neither an RData workspace nor an RDS Seurat object: ", conditionMessage(e))
  obj
})

mapped <- fread(metadata_path)
barcode_candidates <- c("cell_integrated_oc", "cell_canonical", "cell", "cell_id", "barcode", "Cell", "Barcode", colnames(mapped)[1])
barcode_col <- barcode_candidates[barcode_candidates %in% colnames(mapped)][1]
if (is.na(barcode_col)) {
  stop("Cannot identify cell barcode column in mapped metadata")
}
mapped_barcodes <- as.character(mapped[[barcode_col]])
match_idx <- match(colnames(reference), mapped_barcodes)

clone_candidates <- c("cnv_subclone", "CNV_clone", "CNV_subclone")
clone_col <- clone_candidates[clone_candidates %in% colnames(mapped)][1]
group_candidates <- c("interaction_group", "cell_type", "celltype")
group_col <- group_candidates[group_candidates %in% colnames(mapped)][1]
if (is.na(clone_col) || is.na(group_col)) {
  stop("Mapped metadata must contain a CNV-clone column and an interaction/cell-type column")
}

reference$cnv_subclone_curated <- mapped[[clone_col]][match_idx]
reference$interaction_group_curated <- mapped[[group_col]][match_idx]

reference$reference_label <- as.character(reference$interaction_group_curated)
is_tumor <- !is.na(reference$cnv_subclone_curated) & reference$cnv_subclone_curated != ""
reference$reference_label[is_tumor] <- paste0("CNV_", reference$cnv_subclone_curated[is_tumor])
reference$reference_label[is.na(reference$reference_label) | reference$reference_label == ""] <- "Other"

keep_labels <- c(
  "CNV_Subclone_01", "CNV_Subclone_02", "CNV_Subclone_03",
  "CNV_Subclone_04", "CNV_Subclone_05"
)
keep <- reference$reference_label %in% keep_labels |
  grepl("Myeloid|Macro|Monocyte|DC", reference$reference_label, ignore.case = TRUE)

reference <- subset(reference, cells = colnames(reference)[keep])
DefaultAssay(reference) <- "RNA"
reference <- NormalizeData(reference, verbose = FALSE)
reference <- FindVariableFeatures(reference, nfeatures = 3000, verbose = FALSE)
reference <- ScaleData(reference, features = VariableFeatures(reference), verbose = FALSE)
reference <- RunPCA(reference, features = VariableFeatures(reference), npcs = 30, verbose = FALSE)

spatial_objects <- readRDS(file.path(processed_dir, "spatial_objects_curated_scored.rds"))
prediction_rows <- list()

for (sample_id in names(spatial_objects)) {
  query <- spatial_objects[[sample_id]]
  if (ncol(query) < 20) next

  DefaultAssay(query) <- "Spatial"
  query <- NormalizeData(query, verbose = FALSE)
  query <- FindVariableFeatures(query, nfeatures = 3000, verbose = FALSE)

  transfer_status <- "seurat_anchor_transfer"
  predictions <- tryCatch({
    anchors <- FindTransferAnchors(
      reference = reference,
      query = query,
      reference.assay = "RNA",
      query.assay = "Spatial",
      reduction = "pcaproject",
      dims = 1:30,
      features = intersect(VariableFeatures(reference), rownames(query)),
      verbose = FALSE
    )
    TransferData(
      anchorset = anchors,
      refdata = reference$reference_label,
      dims = 1:30,
      verbose = FALSE
    )
  }, error = function(e) {
    transfer_status <<- paste0("fallback_score_based: ", conditionMessage(e))
    target <- query$Target_Subclone02_04_score
    myeloid <- query$SPP1_myeloid_score
    z_target <- as.numeric(scale(target))
    z_myeloid <- as.numeric(scale(myeloid))
    z_target[!is.finite(z_target)] <- 0
    z_myeloid[!is.finite(z_myeloid)] <- 0
    p_target <- plogis(z_target)
    p_myeloid <- plogis(z_myeloid)
    p_other <- pmax(0.05, 1 - pmax(p_target, p_myeloid))
    pred <- data.frame(
      predicted.id = ifelse(p_target >= p_myeloid & p_target >= 0.6, "CNV_Subclone_02_04_like",
                            ifelse(p_myeloid >= 0.6, "SPP1_myeloid_like", "Uncertain")),
      prediction.score.max = pmax(p_target, p_myeloid, p_other),
      prediction.score.CNV_Subclone_02 = p_target / 2,
      prediction.score.CNV_Subclone_04 = p_target / 2,
      prediction.score.SPP1_myeloid_like = p_myeloid,
      row.names = colnames(query)
    )
    pred
  })
  query <- AddMetaData(query, predictions)

  target_cols <- intersect(
    c("prediction.score.CNV_Subclone_02", "prediction.score.CNV_Subclone_04"),
    colnames(query@meta.data)
  )
  if (length(target_cols) > 0) {
    query$prediction.score.CNV_Subclone_02_04 <- rowSums(
      query@meta.data[, target_cols, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    query$prediction.score.CNV_Subclone_02_04 <- NA_real_
  }

  prediction_rows[[length(prediction_rows) + 1]] <- data.table(
    dataset = unique(query$dataset),
    sample_id = sample_id,
    barcode = colnames(query),
    predicted_label = query$predicted.id,
    prediction_max_score = query$prediction.score.max,
    prediction_score_CNV_Subclone_02_04 = query$prediction.score.CNV_Subclone_02_04,
    SPP1_myeloid_score = query$SPP1_myeloid_score,
    Target_Subclone02_04_score = query$Target_Subclone02_04_score,
    transfer_status = transfer_status
  )
  spatial_objects[[sample_id]] <- query
}

saveRDS(
  spatial_objects,
  file.path(processed_dir, "spatial_objects_curated_scored_reference_mapped.rds")
)
fwrite(
  rbindlist(prediction_rows, fill = TRUE),
  file.path(result_dir, "spatial_reference_mapping_predictions.csv.gz")
)

message(sprintf("Reference mapping completed: %s", result_dir))
