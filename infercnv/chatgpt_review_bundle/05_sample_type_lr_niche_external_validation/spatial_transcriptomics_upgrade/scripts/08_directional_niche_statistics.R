#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)
spot <- read_spot_scores(root)

coord_cols <- intersect(c("coord_x", "coord_y", "imagecol", "imagerow", "row", "col"), names(spot))
if (!all(c("coord_x", "coord_y") %in% names(spot)) && length(coord_cols) >= 2) setnames(spot, coord_cols[1:2], c("coord_x", "coord_y"))
coord <- spot[dataset == "GSE203612" & sample_id %in% c("GSM6177614", "GSM6177617") &
                is.finite(coord_x) & is.finite(coord_y)]

targets <- c("CD44", "ITGB1", "target_subclone_02_04_score", "KRAS_hypoxia_score",
             "prediction_score_CNV_Subclone_02", "prediction_score_CNV_Subclone_04")
targets <- intersect(targets, names(coord))
rows <- list()
controls <- list()

for (sid in unique(coord$sample_id)) {
  sdt <- coord[sample_id == sid]
  if (nrow(sdt) <= 7) next
  nn <- knn_index_from_xy(sdt$coord_x, sdt$coord_y, as.integer(config$spatial_statistics$k_neighbors))
  src <- high_by_fraction(sdt$SPP1_myeloid_score, config$spatial_statistics$top_fraction)
  for (target in targets) {
    tgt <- high_by_fraction(sdt[[target]], config$spatial_statistics$top_fraction)
    forward <- neighbor_test(src, tgt, nn, as.integer(config$spatial_statistics$permutations), as.integer(config$seed))
    reverse <- neighbor_test(tgt, src, nn, as.integer(config$spatial_statistics$permutations), as.integer(config$seed) + 11)
    rows[[length(rows) + 1]] <- cbind(data.table(dataset = "GSE203612", sample_id = sid,
                                                 direction = "source_to_target",
                                                 source = "SPP1_myeloid_high",
                                                 target = paste0(target, "_high")), forward)
    rows[[length(rows) + 1]] <- cbind(data.table(dataset = "GSE203612", sample_id = sid,
                                                 direction = "target_to_source",
                                                 source = paste0(target, "_high"),
                                                 target = "SPP1_myeloid_high"), reverse)
  }
  if ("KRAS_hypoxia_score" %in% names(sdt)) {
    ctrl_src <- high_by_fraction(-sdt$SPP1_myeloid_score, config$spatial_statistics$top_fraction)
    ctrl_tgt <- high_by_fraction(sdt$KRAS_hypoxia_score, config$spatial_statistics$top_fraction)
    controls[[length(controls) + 1]] <- cbind(data.table(dataset = "GSE203612", sample_id = sid,
                                                         control = "non_SPP1_high_score",
                                                         target = "KRAS_hypoxia_high"),
                                             neighbor_test(ctrl_src, ctrl_tgt, nn,
                                                           as.integer(config$spatial_statistics$permutations),
                                                           as.integer(config$seed) + 31))
  }
  set.seed(as.integer(config$seed) + 51)
  random_tgt <- sample(src)
  controls[[length(controls) + 1]] <- cbind(data.table(dataset = "GSE203612", sample_id = sid,
                                                       control = "label_permutation",
                                                       target = "random_equal_size_target"),
                                           neighbor_test(src, random_tgt, nn,
                                                         as.integer(config$spatial_statistics$permutations),
                                                         as.integer(config$seed) + 52))
}

main <- rbindlist(rows, fill = TRUE)
neg <- rbindlist(controls, fill = TRUE)
if (nrow(main) > 0) main[, fdr := p.adjust(empirical_p, "BH")]
if (nrow(neg) > 0) neg[, fdr := p.adjust(empirical_p, "BH")]

fwrite(main, file.path(root, "results", "spatial_curated", "directional_niche_statistics.csv"))
fwrite(neg, file.path(root, "results", "spatial_curated", "directional_niche_negative_controls.csv"))
