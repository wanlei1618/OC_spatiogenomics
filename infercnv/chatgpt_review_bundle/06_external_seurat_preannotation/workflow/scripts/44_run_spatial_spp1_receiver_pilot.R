options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(RANN)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
out_dir <- file.path(data_root, "research_spatial_transition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gse211956_root <- Sys.glob(file.path(
  "D:/OC_spatiogenomics", "*", "ovarian_spatial_geo", "GSE211956", "suppl"
))
stopifnot(length(gse211956_root) == 1L)
gse203612_object <- paste0(
  "D:/OC_spatiogenomics/spatial_data/processed/",
  "spatial_objects_curated_scored_reference_mapped.rds"
)

manifest <- rbindlist(list(
  data.table(
    dataset_id = "GSE203612",
    sample_id = c("GSM6177614", "GSM6177617"),
    platform = "10x_Visium",
    ovarian_tumor_relevant = TRUE,
    counts_available = TRUE,
    coordinates_available = TRUE,
    pathology_or_tissue_mask_available = TRUE,
    cell_resolved_or_spot_based = "spot_based",
    selected_for_pilot = TRUE,
    exclusion_reason = NA_character_
  ),
  data.table(
    dataset_id = "GSE211956",
    sample_id = paste0("GSM65061", 10:17),
    platform = "10x_Visium",
    ovarian_tumor_relevant = TRUE,
    counts_available = TRUE,
    coordinates_available = TRUE,
    pathology_or_tissue_mask_available = TRUE,
    cell_resolved_or_spot_based = "spot_based",
    selected_for_pilot = c(rep(TRUE, 3), rep(FALSE, 5)),
    exclusion_reason = c(rep(NA_character_, 3),
                         rep("three-sample-per-dataset pilot cap", 5))
  ),
  data.table(
    dataset_id = "GSE189843",
    sample_id = paste0("GSM57084", 85:96),
    platform = "10x_Visium",
    ovarian_tumor_relevant = TRUE,
    counts_available = TRUE,
    coordinates_available = FALSE,
    pathology_or_tissue_mask_available = FALSE,
    cell_resolved_or_spot_based = "spot_based",
    selected_for_pilot = FALSE,
    exclusion_reason =
      "existing object build failed and usable coordinates were unavailable"
  )
), fill = TRUE)
fwrite(manifest, file.path(out_dir, "spatial_dataset_pilot_manifest.csv"),
       na = "NA")

# Derive a patient-consistent signature from the frozen calibrated dual-method
# epithelial cells; this does not alter the fixed scRNA annotation or CNV calls.
epi_counts <- readRDS(file.path(
  data_root, "research_validation_independent_cnv",
  "GSE154600_complete_final_epithelial_counts.rds"
))
cnv <- fread(file.path(out_dir, "GSE154600_calibrated_cnv_by_cell.csv.gz"))
libs <- Matrix::colSums(epi_counts)
norm <- t(t(epi_counts) / pmax(libs, 1)) * 1e4
norm@x <- log1p(norm@x)
lfc <- rbindlist(lapply(unique(cnv$patient_id), function(pid) {
  dual <- cnv[
    patient_id == pid &
      integrated_calibrated_cnv_evidence == "CALIBRATED_DUAL_METHOD_SUPPORT",
    cell_id
  ]
  other <- cnv[
    patient_id == pid &
      integrated_calibrated_cnv_evidence != "CALIBRATED_DUAL_METHOD_SUPPORT",
    cell_id
  ]
  dual <- intersect(dual, colnames(norm))
  other <- intersect(other, colnames(norm))
  if (length(dual) < 20L || length(other) < 20L) return(NULL)
  data.table(
    patient_id = pid, gene = rownames(norm),
    logfc = Matrix::rowMeans(norm[, dual, drop = FALSE]) -
      Matrix::rowMeans(norm[, other, drop = FALSE])
  )
}))
excluded <- grepl(
  "^(MT-|RPL|RPS|HLA-|IG[HKL]|PTPRC$|CD74$|LST1$|TYROBP$|FCER1G$|C1Q)",
  lfc$gene
) | lfc$gene %in% c(
  "MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "PCNA", "STMN1",
  "TUBA1B", "HMGB2", "CDK1", "CCNB1", "CCNB2"
)
signature_audit <- lfc[!excluded, .(
  n_patients_positive = sum(logfc > 0),
  n_patients_evaluated = .N,
  mean_within_patient_logfc = mean(logfc)
), by = gene][
  n_patients_positive >= 3L & mean_within_patient_logfc > .25
][order(-mean_within_patient_logfc)]
malignant_signature <- head(signature_audit$gene, 30L)
if (length(malignant_signature) < 5L) {
  stop("Fewer than five patient-consistent calibrated CNV signature genes")
}
fwrite(signature_audit,
       file.path(out_dir, "calibrated_malignant_signature_audit.csv"), na = "NA")
rm(epi_counts, norm, lfc)
gc()

extract_counts <- function(object) {
  assay <- object@assays[[object@active.assay]]
  if ("counts" %in% methods::slotNames(assay)) {
    return(as(methods::slot(assay, "counts"), "dgCMatrix"))
  }
  x <- attr(assay, "layers")$counts
  dimnames(x) <- list(
    rownames(attr(assay, "features"))[seq_len(nrow(x))],
    rownames(attr(assay, "cells"))[seq_len(ncol(x))]
  )
  as(x, "dgCMatrix")
}

read_gse211956 <- function(sample_id) {
  sid <- sub("GSM65061", "", sample_id)
  prefix <- file.path(gse211956_root, paste0(sample_id, "_SP", as.integer(sid) - 9L))
  counts <- readMM(gzfile(paste0(prefix, "_matrix.mtx.gz")))
  features <- fread(paste0(prefix, "_features.tsv.gz"), header = FALSE)
  barcodes <- fread(paste0(prefix, "_barcodes.tsv.gz"), header = FALSE)$V1
  rownames(counts) <- make.unique(features$V2)
  colnames(counts) <- barcodes
  pos <- fread(file.path(
    paste0(prefix, "_spatial"), "spatial", "tissue_positions_list.csv"
  ), header = FALSE)
  setnames(pos, c("barcode", "in_tissue", "array_row", "array_col",
                  "coord_y", "coord_x"))
  pos <- pos[in_tissue == 1L & barcode %in% colnames(counts)]
  counts <- as(counts[, pos$barcode, drop = FALSE], "dgCMatrix")
  list(counts = counts, coords = pos[, .(barcode, coord_x, coord_y)])
}

score_spots <- function(counts, coords, dataset_id, sample_id) {
  common <- intersect(coords$barcode, colnames(counts))
  counts <- counts[, common, drop = FALSE]
  coords <- coords[match(common, barcode)]
  lib <- Matrix::colSums(counts)
  norm <- t(t(counts) / pmax(lib, 1)) * 1e4
  norm@x <- log1p(norm@x)
  module <- function(gs) {
    gs <- intersect(gs, rownames(norm))
    if (!length(gs)) rep(NA_real_, ncol(norm)) else
      as.numeric(Matrix::colMeans(norm[gs, , drop = FALSE]))
  }
  gene_positive <- function(g) {
    if (g %in% rownames(counts)) as.vector(counts[g, ] > 0) else
      rep(FALSE, ncol(counts))
  }
  pct <- function(v) frank(v, ties.method = "average") / length(v)
  x <- data.table(
    dataset_id, sample_id, barcode = common,
    coord_x = coords$coord_x, coord_y = coords$coord_y,
    macrophage_score = module(c(
      "C1QA", "C1QB", "C1QC", "LST1", "FCER1G", "TYROBP"
    )),
    spp1_score = module(c("SPP1", "APOC1", "GPNMB", "TREM2", "LPL", "CTSD")),
    c1qc_score = module(c("C1QA", "C1QB", "C1QC")),
    epithelial_score = module(c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7")),
    malignant_score = module(malignant_signature),
    ITGB1_positive = gene_positive("ITGB1"),
    CD44_positive = gene_positive("CD44")
  )
  x[, `:=`(
    macrophage_pct = pct(macrophage_score),
    spp1_pct = pct(spp1_score),
    c1qc_pct = pct(c1qc_score),
    epithelial_pct = pct(epithelial_score),
    malignant_pct = pct(malignant_score)
  )]
  x[, `:=`(
    spp1_sender = macrophage_pct >= .75 & spp1_pct >= .75,
    c1qc_sender = macrophage_pct >= .75 & c1qc_pct >= .75,
    all_macrophage_sender = macrophage_pct >= .75,
    malignant_receiver = epithelial_pct >= .75 & malignant_pct >= .75,
    epithelial_receiver = epithelial_pct >= .75
  )]
  x[, `:=`(
    itgb1_receiver = malignant_receiver & ITGB1_positive,
    cd44_receiver = malignant_receiver & CD44_positive
  )]
  x
}

objects <- readRDS(gse203612_object)
spot_tables <- list()
for (sid in c("GSM6177614", "GSM6177617")) {
  object <- objects[[sid]]
  counts <- extract_counts(object)
  meta <- as.data.table(object@meta.data, keep.rownames = "barcode")
  coords <- meta[, .(barcode, coord_x, coord_y)]
  spot_tables[[sid]] <- score_spots(counts, coords, "GSE203612", sid)
}
rm(objects)
for (sid in paste0("GSM65061", 10:12)) {
  object <- read_gse211956(sid)
  spot_tables[[sid]] <- score_spots(
    object$counts, object$coords, "GSE211956", sid
  )
}
spots <- rbindlist(spot_tables, fill = TRUE)
fwrite(spots, file.path(out_dir, "spatial_spot_scores.csv.gz"),
       compress = "gzip", na = "NA")

tests <- data.table(
  spatial_test = c(
    "SPP1_macrophage_to_ITGB1_receiver",
    "SPP1_macrophage_to_CD44_receiver",
    "C1QC_macrophage_to_ITGB1_receiver",
    "all_macrophage_to_epithelial"
  ),
  sender_col = c("spp1_sender", "spp1_sender", "c1qc_sender",
                 "all_macrophage_sender"),
  receiver_col = c("itgb1_receiver", "cd44_receiver", "itgb1_receiver",
                   "epithelial_receiver")
)

set.seed(20260728)
evaluate_one <- function(x, test_row, order_name, k_neighbors) {
  sender <- x[[test_row$sender_col]]
  receiver <- x[[test_row$receiver_col]]
  n <- nrow(x)
  if (n < 20L || sum(sender) < 3L || sum(receiver) < 3L) {
    return(data.table(
      spatial_test = test_row$spatial_test, adjacency_order = order_name,
      n_spots = n, n_sender = sum(sender), n_receiver = sum(receiver),
      observed_adjacency = NA_real_, expected_adjacency = NA_real_,
      observed_expected_ratio = NA_real_, empirical_p = NA_real_,
      median_nearest_distance = NA_real_
    ))
  }
  xy <- as.matrix(x[, .(coord_x, coord_y)])
  nn <- RANN::nn2(xy, xy, k = min(k_neighbors + 1L, n))$nn.idx[, -1,
    drop = FALSE]
  neighbor_receiver <- rowMeans(matrix(receiver[nn], nrow = n))
  observed <- mean(neighbor_receiver[sender])
  null <- replicate(1000L, mean(neighbor_receiver[
    sample.int(n, sum(sender), replace = FALSE)
  ]))
  receiver_xy <- xy[receiver, , drop = FALSE]
  nearest <- RANN::nn2(receiver_xy, xy[sender, , drop = FALSE], k = 1)$nn.dists
  data.table(
    spatial_test = test_row$spatial_test, adjacency_order = order_name,
    n_spots = n, n_sender = sum(sender), n_receiver = sum(receiver),
    observed_adjacency = observed, expected_adjacency = mean(null),
    observed_expected_ratio = observed / mean(null),
    empirical_p = (1 + sum(null >= observed)) / 1001,
    median_nearest_distance = median(nearest)
  )
}

permutation <- spots[, rbindlist(lapply(seq_len(nrow(tests)), function(i) {
  rbind(
    evaluate_one(.SD, tests[i], "first_order", 6L),
    evaluate_one(.SD, tests[i], "second_order", 18L)
  )
})), by = .(dataset_id, sample_id)]
fwrite(permutation, file.path(out_dir, "spatial_permutation_results.csv"),
       na = "NA")
sample_summary <- spots[, .(
  n_spots = .N,
  n_spp1_sender = sum(spp1_sender),
  n_c1qc_sender = sum(c1qc_sender),
  n_all_macrophage_sender = sum(all_macrophage_sender),
  n_malignant_receiver = sum(malignant_receiver),
  n_itgb1_receiver = sum(itgb1_receiver),
  n_cd44_receiver = sum(cd44_receiver),
  n_epithelial_receiver = sum(epithelial_receiver)
), by = .(dataset_id, sample_id)]
fwrite(sample_summary,
       file.path(out_dir, "spatial_sender_receiver_by_sample.csv"), na = "NA")

spots[, map_class := fcase(
  spp1_sender, "SPP1 macrophage sender",
  c1qc_sender, "C1QC macrophage control",
  itgb1_receiver, "ITGB1 receiver",
  cd44_receiver, "CD44 receiver",
  malignant_receiver, "other malignant-signature receiver",
  default = "other spot"
)]
p1 <- ggplot(spots, aes(coord_x, -coord_y, color = map_class)) +
  geom_point(size = .35) +
  facet_wrap(~sample_id) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Spatial pilot sender and receiver spot map", color = NULL)
ggsave(file.path(out_dir, "spatial_sender_receiver_map.png"),
       p1, width = 13, height = 8, dpi = 180)

p2 <- ggplot(
  permutation[adjacency_order == "first_order"],
  aes(sample_id, observed_expected_ratio, color = spatial_test)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(size = 2) +
  facet_wrap(~dataset_id, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "First-order spatial adjacency effect by biological sample",
    x = NULL, y = "observed / expected", color = "test"
  )
ggsave(file.path(out_dir, "spatial_adjacency_effect_by_sample.png"),
       p2, width = 12, height = 6, dpi = 180)
message("Minimal spatial SPP1 receiver pilot complete")
