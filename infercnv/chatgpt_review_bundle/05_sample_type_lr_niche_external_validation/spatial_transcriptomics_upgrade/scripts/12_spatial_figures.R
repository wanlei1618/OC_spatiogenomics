#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)

spot <- read_spot_scores(root)
fig_dir <- file.path(root, "figures")

if (all(c("coord_x", "coord_y") %in% names(spot))) {
  spatial <- spot[dataset == "GSE203612" & is.finite(coord_x) & is.finite(coord_y)]
  metrics <- intersect(c("SPP1_myeloid_score", "CD44", "ITGB1", "target_subclone_02_04_score",
                         "KRAS_hypoxia_score", "prediction_score_CNV_Subclone_02",
                         "prediction_score_CNV_Subclone_04"), names(spatial))
  if (nrow(spatial) > 0 && length(metrics) > 0) {
    long <- melt(spatial, id.vars = c("sample_id", "barcode", "coord_x", "coord_y"),
                 measure.vars = metrics, variable.name = "metric", value.name = "value")
    p <- ggplot(long, aes(coord_x, coord_y, color = value)) +
      geom_point(size = 0.5) +
      scale_color_viridis_c(option = "magma", na.value = "grey80") +
      coord_fixed() +
      facet_grid(sample_id ~ metric) +
      labs(x = "Spatial x", y = "Spatial y", color = "Score") +
      theme_bw(base_size = 7) +
      theme(strip.text.x = element_text(size = 6), axis.text = element_blank(), axis.ticks = element_blank())
    save_plot_both(p, file.path(fig_dir, "gse203612_spatial_score_panels"), width = 12, height = 5)
  }
}

dir_stats <- read_dt_if_exists(file.path(root, "results", "spatial_curated", "directional_niche_statistics.csv"))
if (nrow(dir_stats) > 0) {
  p <- ggplot(dir_stats[direction == "source_to_target"],
              aes(log2_enrichment, target, color = sample_id)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(size = 2) +
    labs(x = "log2 enrichment", y = "Target", color = "Sample") +
    theme_bw(base_size = 9)
  save_plot_both(p, file.path(fig_dir, "source_target_spatial_zones"), width = 8, height = 4.5)
}

grade <- read_dt_if_exists(file.path(root, "results", "meta", "evidence_grade_table.csv"))
if (nrow(grade) > 0) {
  p <- ggplot(grade, aes(evidence_layer, n_samples, fill = evidence_grade)) +
    geom_col(width = 0.6) +
    coord_flip() +
    labs(x = NULL, y = "Samples", fill = "Evidence grade") +
    theme_bw(base_size = 10)
  save_plot_both(p, file.path(fig_dir, "overall_evidence_grade"), width = 7, height = 3.5)
}

manifest <- data.table(
  figure = list.files(fig_dir, pattern = "\\.(pdf|svg)$", full.names = FALSE),
  path = file.path("figures", list.files(fig_dir, pattern = "\\.(pdf|svg)$", full.names = FALSE)),
  generated_by = "12_spatial_figures.R"
)
fwrite(manifest, file.path(root, "reports", "spatial_figure_manifest.csv"))
