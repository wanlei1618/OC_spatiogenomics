#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)
set.seed(as.integer(config$seed))

spot <- read_spot_scores(root)
ref <- read_reference_predictions(root)
if (nrow(ref) > 0) {
  ref_cols <- intersect(c("prediction_score_CNV_Subclone_02_04", "prediction.score.CNV_Subclone_02",
                          "prediction.score.CNV_Subclone_04", "prediction_score_CNV_Subclone_02",
                          "prediction_score_CNV_Subclone_04"), names(ref))
  spot <- merge(spot, ref[, c("sample_id", "barcode", ref_cols), with = FALSE],
                by = c("sample_id", "barcode"), all.x = TRUE)
}

coord_cols <- intersect(c("coord_x", "coord_y", "imagecol", "imagerow", "row", "col"), names(spot))
if (!all(c("coord_x", "coord_y") %in% names(spot)) && length(coord_cols) >= 2) {
  setnames(spot, coord_cols[1:2], c("coord_x", "coord_y"))
}
coord <- spot[dataset == "GSE203612" & sample_id %in% c("GSM6177614", "GSM6177617") &
                is.finite(coord_x) & is.finite(coord_y)]

metrics <- c("SPP1_myeloid_score", "CD44", "ITGB1", "target_subclone_02_04_score",
             "KRAS_hypoxia_score", "prediction_score_CNV_Subclone_02",
             "prediction_score_CNV_Subclone_04")
metrics <- intersect(metrics, names(coord))
k_values <- c(4, 6, 10, 15)

auto_rows <- list()
enrich_rows <- list()
if (nrow(coord) > 0 && length(metrics) > 0) {
  for (sid in unique(coord$sample_id)) {
    sdt <- coord[sample_id == sid]
    for (k in k_values) {
      if (nrow(sdt) <= k + 1) next
      nn <- knn_index_from_xy(sdt$coord_x, sdt$coord_y, k)
      for (metric in metrics) {
        ac <- moran_geary(sdt[[metric]], nn, n_perm = as.integer(config$spatial_statistics$permutations),
                          seed = as.integer(config$seed) + k)
        ac[, `:=`(dataset = "GSE203612", sample_id = sid, metric = metric, k_neighbors = k)]
        auto_rows[[length(auto_rows) + 1]] <- ac
      }
      if (all(c("SPP1_myeloid_score", "target_subclone_02_04_score") %in% names(sdt))) {
        src <- high_by_fraction(sdt$SPP1_myeloid_score, config$spatial_statistics$top_fraction)
        tgt <- high_by_fraction(sdt$target_subclone_02_04_score, config$spatial_statistics$top_fraction)
        nt <- neighbor_test(src, tgt, nn, as.integer(config$spatial_statistics$permutations),
                            as.integer(config$seed) + k)
        nt[, `:=`(dataset = "GSE203612", sample_id = sid, source = "SPP1_myeloid_high",
                  target = "target_subclone_02_04_high", k_neighbors = k)]
        enrich_rows[[length(enrich_rows) + 1]] <- nt
      }
    }
  }
}

auto <- rbindlist(auto_rows, fill = TRUE)
enrich <- rbindlist(enrich_rows, fill = TRUE)
if (nrow(auto) > 0) auto[, `:=`(moran_q = p.adjust(moran_p, "BH"), geary_q = p.adjust(geary_p, "BH"))]
if (nrow(enrich) > 0) enrich[, empirical_q := p.adjust(empirical_p, "BH")]

fwrite(auto, file.path(root, "results", "spatial_curated", "spatial_autocorrelation.csv"))
fwrite(enrich, file.path(root, "results", "spatial_curated", "multiscale_neighborhood_enrichment.csv"))

p <- ggplot(enrich, aes(factor(k_neighbors), log2_enrichment, color = sample_id)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 2) +
  geom_line(aes(group = sample_id)) +
  labs(x = "k nearest neighbors", y = "log2 enrichment", color = "Sample") +
  theme_bw(base_size = 10)
save_plot_both(p, file.path(root, "figures", "multiscale_neighborhood_forest"), width = 7, height = 4)
