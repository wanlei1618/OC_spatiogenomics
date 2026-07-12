#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)

pred <- read_reference_predictions(root)
out_dir <- file.path(root, "reference_mapping")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(pred) == 0) {
  write_empty_reason(file.path(out_dir, "reference_mapping_predictions.csv.gz"),
                     c("dataset", "sample_id", "barcode", "predicted_label", "prediction_max_score",
                       "confidence_class"),
                     "Reference mapping predictions were not available; run 04_reference_mapping_to_cnv_niches.R first.")
  fwrite(data.table(status = "not_run", reason = "missing_reference_mapping_predictions"),
         file.path(out_dir, "reference_mapping_stability.csv"))
  fwrite(data.table(status = "not_run", reason = "missing_reference_mapping_predictions"),
         file.path(out_dir, "reference_mapping_confusion_or_overlap.csv"))
} else {
  if (!"prediction_max_score" %in% names(pred) && "prediction.score.max" %in% names(pred)) {
    setnames(pred, "prediction.score.max", "prediction_max_score")
  }
  pred[, confidence_class := fifelse(prediction_max_score >= 0.6, "high",
                                     fifelse(prediction_max_score >= 0.4, "moderate", "uncertain"))]
  pred[confidence_class == "uncertain", predicted_label := "Uncertain"]
  fwrite(pred, file.path(out_dir, "reference_mapping_predictions.csv.gz"))
  stability <- pred[, .(
    n_spots = .N,
    median_prediction_score = median(prediction_max_score, na.rm = TRUE),
    low_confidence_fraction = mean(confidence_class == "uncertain", na.rm = TRUE),
    top_label = names(sort(table(predicted_label), decreasing = TRUE))[1]
  ), by = .(dataset, sample_id)]
  stability[, stability_scope := "observed_predictions; rerun with downsampling seeds when Seurat runtime is available"]
  fwrite(stability, file.path(out_dir, "reference_mapping_stability.csv"))
  overlap <- pred[, .N, by = .(sample_id, predicted_label)]
  overlap[, fraction := N / sum(N), by = sample_id]
  fwrite(overlap, file.path(out_dir, "reference_mapping_confusion_or_overlap.csv"))
}

plot_dt <- if (nrow(pred) > 0) pred else data.table(sample_id = "not_run", prediction_max_score = NA_real_)
p <- ggplot(plot_dt, aes(sample_id, prediction_max_score)) +
  geom_boxplot(outlier.shape = NA, fill = "#D9E6F2") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.25, size = 0.6) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Sample", y = "Seurat transfer max score") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_both(p, file.path(root, "figures", "reference_mapping_confidence"), width = 7, height = 4)
