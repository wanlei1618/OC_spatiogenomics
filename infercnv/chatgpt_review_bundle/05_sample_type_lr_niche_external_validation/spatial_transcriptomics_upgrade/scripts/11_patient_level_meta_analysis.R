#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)

neigh <- read_dt_if_exists(file.path(root, "results", "spatial_curated", "directional_niche_statistics.csv"))
resp <- read_dt_if_exists(file.path(root, "results", "gse189843_response", "sample_level_scores.csv"))
meta_dir <- file.path(root, "results", "meta")

spatial <- if (nrow(neigh) > 0) {
  neigh[direction == "source_to_target", .(
    evidence_layer = "coordinate-aware spatial evidence",
    effect = log2_enrichment,
    p = empirical_p,
    fdr = if ("fdr" %in% names(neigh)) fdr else NA_real_,
    note = "Within-sample coordinate-aware neighborhood statistic"
  ), by = .(sample_id, target)]
} else {
  data.table(sample_id = character(), target = character(), evidence_layer = character(),
             effect = numeric(), p = numeric(), fdr = numeric(), note = character())
}

expr <- if (nrow(resp) > 0) {
  melt(resp, id.vars = c("sample_id", "clinical_group", "n_spots"),
       variable.name = "metric", value.name = "sample_level_value")
} else {
  data.table(sample_id = character(), clinical_group = character(), metric = character(),
             sample_level_value = numeric())
}

grade <- data.table(
  evidence_layer = c("coordinate-aware spatial evidence", "expression-level replication"),
  n_samples = c(length(unique(spatial$sample_id)), length(unique(expr$sample_id))),
  scope = c("GSE203612 ovarian Visium samples only", "GSE189843 pretreatment HGSC samples only"),
  conclusion_boundary = c("Descriptive because only two coordinate-aware ovarian samples are available.",
                          "Patient-level expression replication; no coordinate neighborhood claims.")
)
grade[, evidence_grade := fifelse(n_samples >= 6, "moderate", fifelse(n_samples >= 2, "limited", "not_available"))]

fwrite(spatial, file.path(meta_dir, "spatial_evidence_sample_level.csv"))
fwrite(expr, file.path(meta_dir, "expression_evidence_sample_level.csv"))
fwrite(grade, file.path(meta_dir, "evidence_grade_table.csv"))

p <- ggplot(spatial, aes(effect, sample_id, color = target)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(size = 2) +
  labs(x = "log2 enrichment", y = "Sample", color = "Target") +
  theme_bw(base_size = 9)
save_plot_both(p, file.path(root, "figures", "sample_level_effect_forest"), width = 8, height = 4)
