#!/usr/bin/env Rscript

# Build curated Seurat objects for GSE203612 and GSE189843.
#
# Usage:
#   Rscript 02_build_spatial_objects.R D:/OC_spatiogenomics/spatial_data
#
# Requirements:
#   Seurat >= 4, Matrix, data.table, jsonlite
#
# Outputs:
#   processed/spatial_objects_curated.rds
#   processed/spatial_qc_raw_summary.csv
#   processed/spatial_object_build_log.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = FALSE) else
  "D:/OC_spatiogenomics/spatial_data"

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                    winslash = "/", mustWork = FALSE))
manifest_path <- file.path(script_dir, "..", "metadata", "spatial_sample_manifest.csv")
manifest <- fread(manifest_path)
manifest <- manifest[toupper(include_in_ovarian_analysis) == "TRUE"]

processed_dir <- file.path(root, "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

log_rows <- list()
objects <- list()

append_log <- function(dataset, sample_id, status, n_spots = NA_integer_, message = "") {
  log_rows[[length(log_rows) + 1]] <<- data.table(
    dataset = dataset,
    sample_id = sample_id,
    status = status,
    n_spots = n_spots,
    message = message
  )
}

find_one <- function(root_dir, pattern) {
  hits <- list.files(root_dir, pattern = pattern, recursive = TRUE, full.names = TRUE,
                     ignore.case = TRUE)
  if (length(hits) == 0) {
    stop(sprintf("No file matching '%s' under %s", pattern, root_dir))
  }
  if (length(hits) > 1) {
    warning(sprintf("Multiple files match '%s'; using %s", pattern, hits[[1]]))
  }
  hits[[1]]
}

read_10x_h5_rhdf5 <- function(path) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package hdf5r is unavailable and rhdf5 is not installed; cannot read 10x HDF5")
  }
  base <- "matrix"
  data <- rhdf5::h5read(path, file.path(base, "data"))
  indices <- rhdf5::h5read(path, file.path(base, "indices"))
  indptr <- rhdf5::h5read(path, file.path(base, "indptr"))
  shape <- as.integer(rhdf5::h5read(path, file.path(base, "shape")))
  barcodes <- as.character(rhdf5::h5read(path, file.path(base, "barcodes")))
  feature_names <- as.character(rhdf5::h5read(path, file.path(base, "features", "name")))
  feature_ids <- as.character(rhdf5::h5read(path, file.path(base, "features", "id")))
  genes <- make.unique(ifelse(nzchar(feature_names), feature_names, feature_ids))
  mat <- Matrix::sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = shape
  )
  rownames(mat) <- genes
  colnames(mat) <- barcodes
  mat
}

read_mtx_triplet <- function(matrix_file, feature_file, barcode_file) {
  counts <- Matrix::readMM(matrix_file)
  features <- fread(feature_file, header = FALSE)
  barcodes <- fread(barcode_file, header = FALSE)
  gene_col <- if (ncol(features) >= 2) 2 else 1
  rownames(counts) <- make.unique(as.character(features[[gene_col]]))
  colnames(counts) <- as.character(barcodes[[1]])
  counts
}

read_visium_positions <- function(path) {
  pos <- fread(path, header = FALSE)
  if (ncol(pos) >= 6) {
    setnames(pos, names(pos)[1:6], c("barcode", "in_tissue", "array_row", "array_col", "coord_y", "coord_x"))
  } else {
    stop("Unexpected tissue position format: ", path)
  }
  pos
}

build_gse203612 <- function(sample_id) {
  sample_dir <- file.path(root, "raw", "GSE203612", sample_id)
  required <- c(
    file.path(sample_dir, "filtered_feature_bc_matrix.h5"),
    file.path(sample_dir, "spatial", "tissue_positions_list.csv"),
    file.path(sample_dir, "spatial", "scalefactors_json.json"),
    file.path(sample_dir, "spatial", "tissue_hires_image.png")
  )
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop(sprintf("Missing GSE203612 files: %s", paste(missing, collapse = "; ")))
  }

  counts <- read_10x_h5_rhdf5(file.path(sample_dir, "filtered_feature_bc_matrix.h5"))
  positions <- read_visium_positions(file.path(sample_dir, "spatial", "tissue_positions_list.csv"))
  counts <- counts[, intersect(colnames(counts), positions$barcode), drop = FALSE]
  object <- CreateSeuratObject(
    counts = counts,
    assay = "Spatial",
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  md_pos <- positions[match(colnames(object), barcode)]
  object$in_tissue <- md_pos$in_tissue
  object$array_row <- md_pos$array_row
  object$array_col <- md_pos$array_col
  object$coord_x <- md_pos$coord_x
  object$coord_y <- md_pos$coord_y
  object$dataset <- "GSE203612"
  object$sample_id <- sample_id
  object$coordinate_status <- "available"
  object$analysis_level <- "coordinate_aware"
  object$clinical_group <- "not_available"
  object
}

build_gse189843 <- function(sample_id, clinical_group) {
  extract_dir <- file.path(root, "raw", "GSE189843", "extracted")
  matrix_file <- find_one(extract_dir, paste0("^", sample_id, ".*matrix.*\\.mtx$"))
  feature_file <- find_one(extract_dir, paste0("^", sample_id, ".*features.*\\.tsv$"))
  barcode_file <- find_one(extract_dir, paste0("^", sample_id, ".*barcodes.*\\.tsv$"))

  counts <- read_mtx_triplet(matrix_file, feature_file, barcode_file)
  object <- CreateSeuratObject(
    counts = counts,
    assay = "Spatial",
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  object$dataset <- "GSE189843"
  object$sample_id <- sample_id
  object$coordinate_status <- "not_released_in_GEO_supplement"
  object$analysis_level <- "expression_only"
  object$clinical_group <- clinical_group
  object
}

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  sample_id <- row$sample_id
  dataset <- row$dataset

  tryCatch({
    object <- if (dataset == "GSE203612") {
      build_gse203612(sample_id)
    } else if (dataset == "GSE189843") {
      build_gse189843(sample_id, row$clinical_group)
    } else {
      stop(sprintf("Unsupported dataset: %s", dataset))
    }

    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    objects[[sample_id]] <- object
    append_log(dataset, sample_id, "built", ncol(object), "")
  }, error = function(e) {
    append_log(dataset, sample_id, "failed", NA_integer_, conditionMessage(e))
  })
}

qc_summary <- rbindlist(lapply(names(objects), function(sample_id) {
  object <- objects[[sample_id]]
  data.table(
    dataset = unique(object$dataset),
    sample_id = sample_id,
    clinical_group = unique(object$clinical_group),
    coordinate_status = unique(object$coordinate_status),
    n_spots = ncol(object),
    median_nCount = median(object$nCount_Spatial, na.rm = TRUE),
    median_nFeature = median(object$nFeature_Spatial, na.rm = TRUE),
    median_percent_mt = median(object$percent.mt, na.rm = TRUE),
    q95_percent_mt = unname(quantile(object$percent.mt, 0.95, na.rm = TRUE))
  )
}), fill = TRUE)

fwrite(rbindlist(log_rows, fill = TRUE),
       file.path(processed_dir, "spatial_object_build_log.csv"))

if (length(objects) == 0) {
  stop("No spatial object was built. Inspect processed/spatial_object_build_log.csv.")
}

saveRDS(objects, file.path(processed_dir, "spatial_objects_curated.rds"))
fwrite(qc_summary, file.path(processed_dir, "spatial_qc_raw_summary.csv"))

message(sprintf("Built %d curated spatial objects under %s", length(objects), processed_dir))
