options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v4 <- file.path(data_root, "diagnostics_v4_cross_dataset_validation")
v6 <- file.path(data_root, "diagnostics_v6_malignant_receiver_validation")
cleaned <- file.path(data_root, "diagnostics_v2_marker_ready_cleaned")
dir.create(v6, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v6, "malignant_epithelial_receiver_context.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

map_status <- function(a, c) {
  idx <- match(as.character(a$final_cluster), as.character(c$cluster))
  a[, `:=`(
    annotation_status = c$annotation_status[idx],
    canonical_support_n = as.numeric(c$canonical_support_n[idx]),
    incompatible_lineage_program = as.logical(c$incompatible_lineage_program[idx])
  )]
  a
}

consensus <- fread(file.path(v6, "malignancy_consensus_by_cell.csv.gz"))
alpha <- c("ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9")
genes <- c("CD44", "ITGB1", alpha)

summarize_cells <- function(ds, counts, assign, base_cells) {
  base_cells <- intersect(base_cells, colnames(counts))
  assign <- assign[match(base_cells, cell_id)]
  counts <- counts[, base_cells, drop = FALSE]
  lib <- pmax(Matrix::colSums(counts), 1)
  present_genes <- intersect(genes, rownames(counts))
  raw <- counts[present_genes, , drop = FALSE]
  norm <- raw
  if (length(norm@x)) norm@x <- log1p(norm@x * rep(1e4 / lib, diff(norm@p)))
  expr <- function(g) if (g %in% rownames(norm)) as.numeric(norm[g, ]) else rep(NA_real_, ncol(norm))
  det <- function(g) if (g %in% rownames(raw)) as.numeric(raw[g, ] > 0) else rep(NA_real_, ncol(raw))
  dt <- data.table(
    dataset_id = ds,
    cell_id = base_cells,
    patient_id = as.character(assign$patient_id),
    sample_id = as.character(assign$sample_id),
    CD44_expression = expr("CD44"),
    ITGB1_expression = expr("ITGB1"),
    CD44_positive = det("CD44"),
    ITGB1_positive = det("ITGB1")
  )
  for (g in alpha) {
    dt[[paste0(g, "_expression")]] <- expr(g)
    dt[[paste0(g, "_positive")]] <- det(g)
  }
  cc <- consensus[dataset_id == ds, .(cell_id, malignancy_consensus)]
  dt <- merge(dt, cc, by = "cell_id", all.x = TRUE)
  dt[is.na(malignancy_consensus), malignancy_consensus := "NOT_EVALUABLE"]
  dt
}

cell_parts <- list()
for (ds in c("GSE154600", "GSE158722")) {
  a <- fread(file.path(cleaned, ds, "cleaned_cell_assignments.csv.gz"))
  c <- fread(file.path(cleaned, ds, "cleaned_cluster_annotation_template.csv"))
  a <- map_status(a, c)
  inp <- readRDS(file.path(
    data_root, "diagnostics_v2", "objects", ds, "lineage_inputs",
    "Epithelial_like_strategy_input.rds"
  ))
  allowed <- if (ds == "GSE158722")
    c("REVIEW_PATIENT_ENRICHED", "REVIEW_AMBIGUOUS") else "REVIEW_PATIENT_ENRICHED"
  base <- a[
    final_cell_type == "Epithelial" & annotation_status %in% allowed &
      canonical_support_n >= 3 & incompatible_lineage_program != TRUE,
    cell_id
  ]
  cell_parts[[ds]] <- summarize_cells(ds, inp$counts, a, base)
  rm(a, c, inp); gc()
}

a <- fread(file.path(v6, "GSE147082_refined_cell_assignments_v3.csv.gz"))
obj <- readRDS(file.path(data_root, "GSE147082", "objects", "GSE147082_preannotation.rds"))
cnt <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
base <- unique(c(
  a[final_cell_type == "Epithelial", cell_id],
  consensus[dataset_id == "GSE147082" &
              malignancy_consensus %in% c("MALIGNANT_HIGH_CONFIDENCE", "MALIGNANT_SUPPORTIVE"),
            cell_id]
))
cell_parts[["GSE147082"]] <- summarize_cells("GSE147082", cnt, a, base)
rm(obj, cnt, a); gc()
cell_dt <- rbindlist(cell_parts, fill = TRUE)
fwrite(cell_dt, file.path(v6, "malignant_receiver_by_cell.csv.gz"),
       compress = "gzip", na = "NA")

tiers <- c(
  "all_high_confidence_epithelial", "malignant_high_confidence",
  "malignant_supportive", "diploid_epithelial"
)
tier_cells <- function(d, tier) switch(
  tier,
  all_high_confidence_epithelial = rep(TRUE, nrow(d)),
  malignant_high_confidence = d$malignancy_consensus == "MALIGNANT_HIGH_CONFIDENCE",
  malignant_supportive = d$malignancy_consensus == "MALIGNANT_SUPPORTIVE",
  diploid_epithelial = d$malignancy_consensus == "DIPLOID_SUPPORTIVE"
)

sample_keys <- unique(cell_dt[, .(dataset_id, patient_id, sample_id)])
receiver <- rbindlist(lapply(seq_len(nrow(sample_keys)), function(i) {
  key <- sample_keys[i]
  d <- cell_dt[
    dataset_id == key$dataset_id & patient_id == key$patient_id &
      sample_id == key$sample_id
  ]
  rbindlist(lapply(tiers, function(tier) {
    z <- d[tier_cells(d, tier)]
    if (!nrow(z)) return(data.table(
      dataset_id = key$dataset_id, patient_id = key$patient_id,
      sample_id = key$sample_id, receiver_tier = tier,
      n_receiver_cells = 0L, detectable_alpha_integrins = paste(
        intersect(alpha, sub("_positive$", "", grep(
          "^ITGA[0-9V]+_positive$", names(d), value = TRUE
        ))), collapse = ";"
      ),
      CD44_average_expression = NA_real_, CD44_positive_fraction = NA_real_,
      ITGB1_average_expression = NA_real_, ITGB1_positive_fraction = NA_real_,
      CD44_ITGB1_copositive_fraction = NA_real_,
      ITGB1_alpha_copositive_fraction = NA_real_,
      dominant_alpha_partner = NA_character_
    ))
    alpha_pos_cols <- paste0(alpha, "_positive")
    alpha_pos <- as.matrix(z[, ..alpha_pos_cols])
    any_alpha <- rowSums(alpha_pos > 0, na.rm = TRUE) > 0
    alpha_frac <- colMeans(alpha_pos > 0, na.rm = TRUE)
    data.table(
      dataset_id = key$dataset_id, patient_id = key$patient_id,
      sample_id = key$sample_id, receiver_tier = tier,
      n_receiver_cells = nrow(z),
      detectable_alpha_integrins = paste(
        alpha[colSums(!is.na(alpha_pos)) > 0], collapse = ";"
      ),
      CD44_average_expression = mean(z$CD44_expression, na.rm = TRUE),
      CD44_positive_fraction = mean(z$CD44_positive, na.rm = TRUE),
      ITGB1_average_expression = mean(z$ITGB1_expression, na.rm = TRUE),
      ITGB1_positive_fraction = mean(z$ITGB1_positive, na.rm = TRUE),
      CD44_ITGB1_copositive_fraction =
        mean(z$CD44_positive > 0 & z$ITGB1_positive > 0, na.rm = TRUE),
      ITGB1_alpha_copositive_fraction =
        mean(z$ITGB1_positive > 0 & any_alpha, na.rm = TRUE),
      dominant_alpha_partner = if (all(is.na(alpha_frac))) NA_character_ else
        alpha[which.max(alpha_frac)],
      ITGA4_positive_fraction = alpha_frac[["ITGA4_positive"]],
      ITGA5_positive_fraction = alpha_frac[["ITGA5_positive"]],
      ITGAV_positive_fraction = alpha_frac[["ITGAV_positive"]],
      ITGA8_positive_fraction = alpha_frac[["ITGA8_positive"]],
      ITGA9_positive_fraction = alpha_frac[["ITGA9_positive"]]
    )
  }), fill = TRUE)
}), fill = TRUE)

support <- function(n, fraction) fcase(
  n < 20, "NOT_EVALUABLE",
  is.na(fraction) | fraction == 0, "NOT_DETECTED",
  fraction >= .10, "SUPPORTED",
  default = "DETECTED_LOW"
)
receiver[, CD44_receiver_support :=
           support(n_receiver_cells, CD44_positive_fraction)]
receiver[, ITGB1_receiver_support :=
           support(n_receiver_cells, ITGB1_positive_fraction)]
receiver[, ITGB1_alpha_partner_support :=
           support(n_receiver_cells, ITGB1_alpha_copositive_fraction)]
receiver[, dual_CD44_ITGB1_support :=
           support(n_receiver_cells, CD44_ITGB1_copositive_fraction)]
receiver[, interpretation_note :=
           "ITGB1 expression alone is not evidence of a complete functional integrin receptor"]
fwrite(receiver, out, na = "NA")

pdt <- melt(
  receiver[receiver_tier %in% c("malignant_high_confidence", "malignant_supportive")],
  id.vars = c("dataset_id", "patient_id", "sample_id", "receiver_tier", "n_receiver_cells"),
  measure.vars = c(
    "CD44_positive_fraction", "ITGB1_positive_fraction",
    "CD44_ITGB1_copositive_fraction", "ITGB1_alpha_copositive_fraction"
  ),
  variable.name = "receiver_feature", value.name = "positive_fraction"
)
p <- ggplot(pdt, aes(receiver_feature, patient_id, fill = positive_fraction)) +
  geom_tile(color = "white") +
  facet_grid(receiver_tier ~ dataset_id, scales = "free", space = "free") +
  scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "Malignancy-stratified CD44, ITGB1, and alpha-integrin context",
    subtitle = "Single-method malignant-supportive cells are shown separately from unavailable high-confidence consensus",
    x = NULL, y = NULL
  )
ggsave(file.path(v6, "05_malignant_receiver_context.png"),
       p, width = 12, height = 6, dpi = 180)
message("Malignancy-stratified receiver context complete")
