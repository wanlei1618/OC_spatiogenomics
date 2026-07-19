#!/usr/bin/env Rscript

# Curated scoring and coordinate-aware spatial statistics.
#
# This script deliberately separates:
#   1) expression-level association (all included samples), and
#   2) coordinate-aware neighborhood analysis (GSE203612 ovarian samples only).
#
# GSE189843 is never used for coordinate-based claims because GEO does not
# release tissue_positions/scalefactors files in the supplementary archive.

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
config_path <- if (length(args) >= 1) args[[1]] else file.path(script_dir, "..", "config", "spatial_config.yml")
config <- yaml::read_yaml(config_path)

root <- normalizePath(config$project_root, winslash = "/", mustWork = FALSE)
processed_dir <- file.path(root, "processed")
result_dir <- file.path(root, "results", "spatial_curated")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(as.integer(config$seed))
objects <- readRDS(file.path(processed_dir, "spatial_objects_curated.rds"))

get_assay_data_compat <- function(object, assay = "Spatial", slot = "data") {
  tryCatch(
    GetAssayData(object, assay = assay, layer = slot),
    error = function(e) GetAssayData(object, assay = assay, slot = slot)
  )
}

score_gene_set <- function(data_mat, genes) {
  present <- intersect(genes, rownames(data_mat))
  if (length(present) < 2) {
    warning(sprintf("Only %d genes found for score: %s", length(present), paste(genes, collapse = ",")))
  }
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(data_mat)))
  }
  x <- as.matrix(data_mat[present, , drop = FALSE])
  gene_sd <- apply(x, 1, sd)
  gene_sd[!is.finite(gene_sd) | gene_sd == 0] <- 1
  z <- sweep(sweep(x, 1, rowMeans(x), "-"), 1, gene_sd, "/")
  colMeans(z, na.rm = TRUE)
}

get_gene_expression <- function(data_mat, gene) {
  if (!gene %in% rownames(data_mat)) {
    return(rep(NA_real_, ncol(data_mat)))
  }
  as.numeric(data_mat[gene, ])
}

adaptive_filter <- function(object) {
  md <- object@meta.data
  count_col <- grep("^nCount_", colnames(md), value = TRUE)[1]
  feature_col <- grep("^nFeature_", colnames(md), value = TRUE)[1]
  if (is.na(count_col) || is.na(feature_col)) {
    stop("Cannot identify nCount/nFeature metadata columns")
  }

  keep <- md[[feature_col]] >= as.numeric(config$qc$min_features) &
    md[[count_col]] >= as.numeric(config$qc$min_counts) &
    md$percent.mt <= as.numeric(config$qc$max_percent_mt)

  subset(object, cells = rownames(md)[keep])
}

get_coordinates <- function(object) {
  if (all(c("coord_x", "coord_y") %in% colnames(object@meta.data))) {
    return(data.table(
      barcode = colnames(object),
      coord_x = as.numeric(object$coord_x),
      coord_y = as.numeric(object$coord_y)
    ))
  }
  if (length(Images(object)) == 0) {
    return(NULL)
  }
  image_name <- Images(object)[1]
  coords <- tryCatch(
    GetTissueCoordinates(object, image = image_name),
    error = function(e) NULL
  )
  if (is.null(coords) || nrow(coords) == 0) {
    return(NULL)
  }

  coords <- as.data.frame(coords)
  coords$barcode <- rownames(coords)
  preferred <- list(c("row", "col"), c("imagerow", "imagecol"), c("x", "y"))
  chosen <- NULL
  for (pair in preferred) {
    if (all(pair %in% colnames(coords))) {
      chosen <- pair
      break
    }
  }
  if (is.null(chosen)) {
    numeric_cols <- names(coords)[vapply(coords, is.numeric, logical(1))]
    if (length(numeric_cols) < 2) {
      return(NULL)
    }
    chosen <- numeric_cols[1:2]
  }
  data.table(
    barcode = coords$barcode,
    coord_x = as.numeric(coords[[chosen[1]]]),
    coord_y = as.numeric(coords[[chosen[2]]])
  )
}

knn_indices <- function(coords, k = 6L) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for coordinate-aware kNN. Install with install.packages('RANN').")
  }
  k_use <- min(as.integer(k) + 1L, nrow(coords))
  nn <- RANN::nn2(as.matrix(coords), k = k_use)$nn.idx
  if (ncol(nn) <= 1) {
    stop("Too few spots for kNN analysis")
  }
  nn[, -1, drop = FALSE]
}

permutation_neighbor_test <- function(source_high, target_high, nn_idx, n_perm = 1000L) {
  source_idx <- which(source_high)
  if (length(source_idx) == 0 || sum(target_high) == 0) {
    return(list(
      observed = NA_real_, expected = mean(target_high), enrichment = NA_real_,
      empirical_p = NA_real_, n_edges = 0L
    ))
  }

  target_neighbor <- target_high[nn_idx[source_idx, , drop = FALSE]]
  observed <- mean(target_neighbor)
  expected <- mean(target_high)
  permuted <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    permuted_target <- sample(target_high, replace = FALSE)
    permuted[b] <- mean(permuted_target[nn_idx[source_idx, , drop = FALSE]])
  }
  p <- (sum(permuted >= observed) + 1) / (n_perm + 1)
  list(
    observed = observed,
    expected = expected,
    enrichment = ifelse(expected > 0, observed / expected, NA_real_),
    empirical_p = p,
    n_edges = length(target_neighbor)
  )
}

correlation_rows <- list()
neighborhood_rows <- list()
spot_rows <- list()
qc_rows <- list()

for (sample_id in names(objects)) {
  object <- objects[[sample_id]]
  dataset <- unique(object$dataset)
  clinical_group <- unique(object$clinical_group)
  coordinate_status <- unique(object$coordinate_status)
  analysis_level <- unique(object$analysis_level)

  n_before <- ncol(object)
  object <- adaptive_filter(object)
  n_after <- ncol(object)
  if (n_after < 20) {
    warning(sprintf("%s has fewer than 20 spots after QC; skipping", sample_id))
    next
  }

  DefaultAssay(object) <- "Spatial"
  object <- NormalizeData(
    object,
    normalization.method = config$normalization$method,
    scale.factor = as.numeric(config$normalization$scale_factor),
    verbose = FALSE
  )
  data_mat <- get_assay_data_compat(object, assay = "Spatial", slot = "data")

  spp1_myeloid <- score_gene_set(data_mat, unlist(config$gene_sets$spp1_myeloid))
  target_0204 <- score_gene_set(data_mat, unlist(config$gene_sets$target_subclone_02_04))
  kras_hypoxia <- score_gene_set(data_mat, unlist(config$gene_sets$kras_hypoxia))
  spp1_expr <- get_gene_expression(data_mat, "SPP1")
  cd44_expr <- get_gene_expression(data_mat, "CD44")
  itgb1_expr <- get_gene_expression(data_mat, "ITGB1")

  object$SPP1_myeloid_score <- spp1_myeloid
  object$Target_Subclone02_04_score <- target_0204
  object$KRAS_hypoxia_score <- kras_hypoxia
  object$SPP1_CD44_expr_product <- spp1_expr * cd44_expr
  object$SPP1_ITGB1_expr_product <- spp1_expr * itgb1_expr

  cor_s <- suppressWarnings(cor.test(spp1_myeloid, target_0204, method = "spearman", exact = FALSE))
  cor_p <- suppressWarnings(cor.test(spp1_myeloid, target_0204, method = "pearson"))

  correlation_rows[[length(correlation_rows) + 1]] <- data.table(
    dataset = dataset,
    sample_id = sample_id,
    clinical_group = clinical_group,
    coordinate_status = coordinate_status,
    analysis_level = analysis_level,
    n_spots_raw = n_before,
    n_spots_qc = n_after,
    pearson_r = unname(cor_p$estimate),
    pearson_p = cor_p$p.value,
    spearman_r = unname(cor_s$estimate),
    spearman_p = cor_s$p.value,
    mean_SPP1_CD44_score = mean(object$SPP1_CD44_expr_product, na.rm = TRUE),
    mean_SPP1_ITGB1_score = mean(object$SPP1_ITGB1_expr_product, na.rm = TRUE)
  )

  coord_dt <- get_coordinates(object)
  if (!is.null(coord_dt) && coordinate_status == "available") {
    common <- intersect(colnames(object), coord_dt$barcode)
    coord_dt <- coord_dt[match(common, barcode)]
    md <- object@meta.data[common, , drop = FALSE]
    nn <- knn_indices(coord_dt[, .(coord_x, coord_y)], k = config$spatial_statistics$k_neighbors)

    source_cut <- quantile(md$SPP1_myeloid_score, 1 - config$spatial_statistics$top_fraction,
                           na.rm = TRUE)
    target_cut <- quantile(md$Target_Subclone02_04_score, 1 - config$spatial_statistics$top_fraction,
                           na.rm = TRUE)
    source_high <- md$SPP1_myeloid_score >= source_cut
    target_high <- md$Target_Subclone02_04_score >= target_cut

    perm <- permutation_neighbor_test(
      source_high = source_high,
      target_high = target_high,
      nn_idx = nn,
      n_perm = as.integer(config$spatial_statistics$permutations)
    )

    neighborhood_rows[[length(neighborhood_rows) + 1]] <- data.table(
      dataset = dataset,
      sample_id = sample_id,
      clinical_group = clinical_group,
      coordinate_status = coordinate_status,
      k_neighbors = ncol(nn),
      top_fraction = as.numeric(config$spatial_statistics$top_fraction),
      observed_neighbor_fraction = perm$observed,
      expected_global_fraction = perm$expected,
      enrichment_ratio = perm$enrichment,
      empirical_p = perm$empirical_p,
      n_source_spots = sum(source_high),
      n_target_spots = sum(target_high),
      n_neighbor_edges = perm$n_edges,
      analysis_validity = "coordinate_aware_exploratory"
    )
  }

  coord_for_spots <- get_coordinates(object)
  spot_dt <- data.table(
    dataset = dataset,
    sample_id = sample_id,
    barcode = colnames(object),
    clinical_group = clinical_group,
    coordinate_status = coordinate_status,
    SPP1_myeloid_score = object$SPP1_myeloid_score,
    Target_Subclone02_04_score = object$Target_Subclone02_04_score,
    KRAS_hypoxia_score = object$KRAS_hypoxia_score,
    SPP1_CD44_expr_product = object$SPP1_CD44_expr_product,
    SPP1_ITGB1_expr_product = object$SPP1_ITGB1_expr_product
  )
  if (!is.null(coord_for_spots)) {
    spot_dt <- merge(spot_dt, coord_for_spots, by = "barcode", all.x = TRUE)
  } else {
    spot_dt[, `:=`(coord_x = NA_real_, coord_y = NA_real_)]
  }
  spot_rows[[length(spot_rows) + 1]] <- spot_dt

  qc_rows[[length(qc_rows) + 1]] <- data.table(
    dataset = dataset,
    sample_id = sample_id,
    n_spots_raw = n_before,
    n_spots_qc = n_after,
    retained_fraction = n_after / n_before
  )
  objects[[sample_id]] <- object
}

cor_dt <- rbindlist(correlation_rows, fill = TRUE)
cor_dt[, pearson_q_within_dataset := p.adjust(pearson_p, method = "BH"), by = dataset]
cor_dt[, spearman_q_within_dataset := p.adjust(spearman_p, method = "BH"), by = dataset]
setcolorder(cor_dt, c(
  "dataset", "sample_id", "clinical_group", "coordinate_status", "analysis_level",
  "n_spots_raw", "n_spots_qc", "pearson_r", "pearson_p", "pearson_q_within_dataset",
  "spearman_r", "spearman_p", "spearman_q_within_dataset",
  "mean_SPP1_CD44_score", "mean_SPP1_ITGB1_score"
))

neigh_dt <- rbindlist(neighborhood_rows, fill = TRUE)
if (nrow(neigh_dt) > 0) {
  neigh_dt[, empirical_q := p.adjust(empirical_p, method = "BH")]
}

fwrite(cor_dt, file.path(result_dir, "spatial_correlation_curated.csv"))
fwrite(neigh_dt, file.path(result_dir, "spatial_neighborhood_enrichment_curated.csv"))
fwrite(rbindlist(spot_rows, fill = TRUE),
       file.path(result_dir, "spatial_spot_scores_curated.csv.gz"))
fwrite(rbindlist(qc_rows, fill = TRUE),
       file.path(result_dir, "spatial_qc_filtered_summary.csv"))

saveRDS(objects, file.path(processed_dir, "spatial_objects_curated_scored.rds"))
message(sprintf("Curated spatial outputs written to %s", result_dir))
