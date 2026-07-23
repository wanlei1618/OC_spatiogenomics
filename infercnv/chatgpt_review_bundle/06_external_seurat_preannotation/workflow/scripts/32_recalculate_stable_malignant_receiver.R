options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v61, "stable_malignant_receiver_context_by_patient.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

lineage_names <- c("Epithelial_like", "T_NK_like", "B_Plasma_like", "Myeloid_like")
mats <- setNames(lapply(lineage_names, function(x) readRDS(file.path(
  data_root, "diagnostics_v2", "objects", "GSE154600", "lineage_inputs",
  paste0(x, "_strategy_input.rds")
))$counts), lineage_names)
manifest <- fread(file.path(v61, "GSE154600_copykat_target_manifest.csv.gz"))
stability <- fread(file.path(v61, "GSE154600_copykat_stability_by_cell.csv.gz"))

collect_from_manifest <- function(cells) {
  z <- manifest[cell_id %in% cells & available_for_copykat == TRUE]
  parts <- lapply(lineage_names, function(src) {
    ids <- z[source_lineage_input == src, cell_id]
    ids <- intersect(ids, colnames(mats[[src]]))
    if (length(ids)) mats[[src]][, ids, drop = FALSE] else NULL
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  genes <- Reduce(intersect, lapply(parts, rownames))
  ans <- do.call(cbind, lapply(parts, function(x) x[genes, , drop = FALSE]))
  ans[, !duplicated(colnames(ans)), drop = FALSE]
}

counts <- collect_from_manifest(manifest[available_for_copykat == TRUE, cell_id])
genes <- c("CD44", "ITGB1", "ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9")
present <- intersect(genes, rownames(counts))
raw <- counts[present, , drop = FALSE]
lib <- pmax(Matrix::colSums(counts), 1)
norm <- raw
if (length(norm@x)) norm@x <- log1p(norm@x * rep(1e4 / lib, diff(norm@p)))
expr <- function(g) if (g %in% rownames(norm)) as.numeric(norm[g, ]) else rep(NA_real_, ncol(norm))
det <- function(g) if (g %in% rownames(raw)) as.numeric(raw[g, ] > 0) else rep(NA_real_, ncol(raw))
cell <- data.table(cell_id = colnames(counts))
cell <- merge(
  cell,
  stability[, .(cell_id, patient_id, stability_class, malignancy_evidence)],
  by = "cell_id", all.x = TRUE
)
for (g in genes) {
  cell[[paste0(g, "_expression")]] <- expr(g)
  cell[[paste0(g, "_positive")]] <- det(g)
}

tiers <- c(
  "all_final_epithelial", "stable_malignant_supportive",
  "stable_diploid_like", "unstable_or_not_defined"
)
tier_mask <- function(d, tier) switch(
  tier,
  all_final_epithelial = rep(TRUE, nrow(d)),
  stable_malignant_supportive = d$stability_class == "STABLE_ANEUPLOID",
  stable_diploid_like = d$stability_class == "STABLE_DIPLOID",
  unstable_or_not_defined =
    d$stability_class %in% c("UNSTABLE_DISCORDANT", "MOSTLY_NOT_DEFINED", "NOT_SUBMITTED")
)
alpha <- c("ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9")

receiver <- rbindlist(lapply(unique(manifest$patient_id), function(pt) {
  d <- cell[patient_id == pt]
  final_n <- manifest[patient_id == pt, .N]
  missing_n <- manifest[patient_id == pt & available_for_copykat == FALSE, .N]
  rbindlist(lapply(tiers, function(tier) {
    z <- d[tier_mask(d, tier)]
    if (!nrow(z)) return(data.table(
      dataset_id = "GSE154600", patient_id = pt, receiver_tier = tier,
      n_final_epithelial = final_n, n_cells_without_counts = missing_n,
      n_receiver_cells = 0L, CD44_average_expression = NA_real_,
      CD44_positive_fraction = NA_real_, ITGB1_average_expression = NA_real_,
      ITGB1_positive_fraction = NA_real_,
      CD44_ITGB1_copositive_fraction = NA_real_,
      ITGA4_positive_fraction = NA_real_, ITGA5_positive_fraction = NA_real_,
      ITGAV_positive_fraction = NA_real_, ITGA8_positive_fraction = NA_real_,
      ITGA9_positive_fraction = NA_real_,
      ITGB1_any_alpha_copositive_fraction = NA_real_,
      dominant_alpha_partner = NA_character_
    ))
    alpha_cols <- paste0(alpha, "_positive")
    alpha_mat <- as.matrix(z[, ..alpha_cols])
    alpha_frac <- colMeans(alpha_mat > 0, na.rm = TRUE)
    any_alpha <- rowSums(alpha_mat > 0, na.rm = TRUE) > 0
    data.table(
      dataset_id = "GSE154600", patient_id = pt, receiver_tier = tier,
      n_final_epithelial = final_n, n_cells_without_counts = missing_n,
      n_receiver_cells = nrow(z),
      CD44_average_expression = mean(z$CD44_expression, na.rm = TRUE),
      CD44_positive_fraction = mean(z$CD44_positive, na.rm = TRUE),
      ITGB1_average_expression = mean(z$ITGB1_expression, na.rm = TRUE),
      ITGB1_positive_fraction = mean(z$ITGB1_positive, na.rm = TRUE),
      CD44_ITGB1_copositive_fraction =
        mean(z$CD44_positive > 0 & z$ITGB1_positive > 0, na.rm = TRUE),
      ITGA4_positive_fraction = alpha_frac[["ITGA4_positive"]],
      ITGA5_positive_fraction = alpha_frac[["ITGA5_positive"]],
      ITGAV_positive_fraction = alpha_frac[["ITGAV_positive"]],
      ITGA8_positive_fraction = alpha_frac[["ITGA8_positive"]],
      ITGA9_positive_fraction = alpha_frac[["ITGA9_positive"]],
      ITGB1_any_alpha_copositive_fraction =
        mean(z$ITGB1_positive > 0 & any_alpha, na.rm = TRUE),
      dominant_alpha_partner = alpha[which.max(alpha_frac)]
    )
  }), fill = TRUE)
}), fill = TRUE)

support <- function(n, fraction) fcase(
  n < 20, "NOT_EVALUABLE",
  is.na(fraction) | fraction == 0, "NOT_DETECTED",
  fraction >= .10, "SUPPORTED",
  default = "DETECTED_LOW"
)
receiver[, CD44_receiver_support := support(n_receiver_cells, CD44_positive_fraction)]
receiver[, ITGB1_receiver_support := support(n_receiver_cells, ITGB1_positive_fraction)]
receiver[, ITGB1_alpha_partner_support :=
           support(n_receiver_cells, ITGB1_any_alpha_copositive_fraction)]
receiver[, dual_CD44_ITGB1_support :=
           support(n_receiver_cells, CD44_ITGB1_copositive_fraction)]
receiver[, interpretation_note :=
           "Stable CopyKAT aneuploid is single-method supportive evidence; ITGB1 alone is not a complete receptor"]
fwrite(receiver, out, na = "NA")

status_aggregate <- function(x) {
  if (any(x == "SUPPORTED")) "SUPPORTED" else
    if (any(x == "DETECTED_LOW")) "DETECTED_LOW" else
      if (all(x == "NOT_EVALUABLE")) "NOT_EVALUABLE" else "NOT_DETECTED"
}
summary <- receiver[, .(
  n_patients = .N,
  n_evaluable_patients = sum(n_receiver_cells >= 20),
  total_receiver_cells = sum(n_receiver_cells),
  median_CD44_positive_fraction = median(CD44_positive_fraction[n_receiver_cells >= 20], na.rm = TRUE),
  median_ITGB1_positive_fraction = median(ITGB1_positive_fraction[n_receiver_cells >= 20], na.rm = TRUE),
  median_dual_positive_fraction = median(CD44_ITGB1_copositive_fraction[n_receiver_cells >= 20], na.rm = TRUE),
  median_ITGB1_any_alpha_copositive_fraction =
    median(ITGB1_any_alpha_copositive_fraction[n_receiver_cells >= 20], na.rm = TRUE),
  CD44_receiver_status = status_aggregate(CD44_receiver_support),
  ITGB1_receiver_status = status_aggregate(ITGB1_receiver_support),
  ITGB1_alpha_partner_status = status_aggregate(ITGB1_alpha_partner_support),
  dual_receiver_status = status_aggregate(dual_CD44_ITGB1_support)
), by = .(dataset_id, receiver_tier)]
fwrite(summary, file.path(v61, "stable_malignant_receiver_context_summary.csv"), na = "NA")

pdt <- melt(
  receiver[receiver_tier == "stable_malignant_supportive"],
  id.vars = c("patient_id", "n_receiver_cells"),
  measure.vars = c(
    "CD44_positive_fraction", "ITGB1_positive_fraction",
    "CD44_ITGB1_copositive_fraction", "ITGB1_any_alpha_copositive_fraction"
  ),
  variable.name = "receiver_feature", value.name = "positive_fraction"
)
p <- ggplot(pdt, aes(receiver_feature, patient_id, fill = positive_fraction)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(n_receiver_cells >= 20,
                               sprintf("%.2f", positive_fraction), "NE")), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Stable malignant-supportive receiver context",
       subtitle = "NE: fewer than 20 stable malignant-supportive cells",
       x = NULL, y = NULL)
ggsave(file.path(v61, "04_stable_malignant_receiver_context.png"),
       p, width = 10, height = 5.5, dpi = 180)
message("Stable malignant receiver context complete")
