#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)
out_dir <- file.path(root, "results", "spatial_curated")

spot <- read_spot_scores(root)
neigh <- read_neighborhood(root)
qc_path <- file.path(root, "processed", "spatial_qc_filtered_summary.csv")
qc <- read_dt_if_exists(qc_path)
if (nrow(qc) == 0) {
  qc <- read_dt_if_exists(file.path(root, "results", "spatial_curated", "spatial_qc_filtered_summary.csv"))
}

grid <- CJ(
  min_features = c(100, 200, 300),
  min_counts = c(300, 500, 1000),
  max_percent_mt = c(20, 30, 40),
  top_fraction = c(0.15, 0.25, 0.35)
)

if (!all(c("nFeature_Spatial", "nCount_Spatial", "percent.mt") %in% names(spot))) {
  if (nrow(qc) > 0) {
    grid_rows <- grid[, {
      direction <- if (nrow(neigh) > 0 && "enrichment_ratio" %in% names(neigh)) {
        mean(neigh$enrichment_ratio, na.rm = TRUE) > 1
      } else {
        NA
      }
      .(sample_id = qc$sample_id,
        n_spots_retained = qc$n_spots_qc,
        score_spearman = NA_real_,
        enrichment_ratio = if (nrow(neigh) > 0) mean(neigh$enrichment_ratio, na.rm = TRUE) else NA_real_,
        empirical_p = if (nrow(neigh) > 0) min(neigh$empirical_p, na.rm = TRUE) else NA_real_,
        conclusion_direction_positive = direction,
        unstable_low_spots = qc$n_spots_qc < 50)
    }, by = .(min_features, min_counts, max_percent_mt, top_fraction)]
  } else {
    grid_rows <- grid[, .(sample_id = character(), n_spots_retained = integer(),
                          score_spearman = numeric(), enrichment_ratio = numeric(),
                          empirical_p = numeric(), conclusion_direction_positive = logical(),
                          unstable_low_spots = logical()),
                      by = .(min_features, min_counts, max_percent_mt, top_fraction)]
  }
} else {
  grid_rows <- grid[, {
    rows <- spot[nFeature_Spatial >= min_features &
                   nCount_Spatial >= min_counts &
                   percent.mt <= max_percent_mt, ]
    rows[, {
      corv <- safe_cor(SPP1_myeloid_score, target_subclone_02_04_score, "spearman")
      .(n_spots_retained = .N,
        score_spearman = unname(corv["estimate"]),
        enrichment_ratio = NA_real_,
        empirical_p = unname(corv["p"]),
        conclusion_direction_positive = is.finite(corv["estimate"]) && corv["estimate"] > 0,
        unstable_low_spots = .N < 50)
    }, by = sample_id]
  }, by = .(min_features, min_counts, max_percent_mt, top_fraction)]
}

grid_rows[, empirical_q := p.adjust(empirical_p, method = "BH")]
summary <- grid_rows[, .(
  n_parameter_sets = .N,
  median_retained_spots = median(n_spots_retained, na.rm = TRUE),
  positive_direction_rate = mean(conclusion_direction_positive, na.rm = TRUE),
  unstable_parameter_rate = mean(unstable_low_spots, na.rm = TRUE)
), by = sample_id]

fwrite(grid_rows, file.path(out_dir, "qc_sensitivity_grid.csv"))
fwrite(summary, file.path(out_dir, "qc_sensitivity_summary.csv"))

plot_dt <- copy(grid_rows)
plot_dt[, setting := paste0("F", min_features, "_C", min_counts, "_MT", max_percent_mt)]
p <- ggplot(plot_dt, aes(setting, factor(top_fraction), fill = score_spearman)) +
  geom_tile(color = "white", size = 0.2) +
  facet_wrap(~ sample_id, scales = "free_x") +
  scale_fill_gradient2(low = "#276FBF", mid = "white", high = "#C44536", na.value = "grey85") +
  labs(x = "QC setting", y = "Top fraction", fill = "Spearman") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
save_plot_both(p, file.path(root, "figures", "qc_sensitivity_heatmap"), width = 10, height = 6)
