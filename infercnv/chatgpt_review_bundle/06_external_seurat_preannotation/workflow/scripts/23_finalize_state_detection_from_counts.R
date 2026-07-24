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
v3 <- file.path(data_root, "diagnostics_v3_remaining_datasets")
v4 <- file.path(data_root, "diagnostics_v4_cross_dataset_validation")
v6 <- file.path(data_root, "diagnostics_v6_malignant_receiver_validation")
cleaned <- file.path(data_root, "diagnostics_v2_marker_ready_cleaned")
dir.create(v6, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v6, "state_detection_by_sample.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

roles <- c(
  GSE154600 = "primary_tumor_ecosystem",
  GSE158722 = "malignant_fluid_tumor_ecosystem",
  GSE147082 = "tumor_sensitivity_validation",
  GSE151214 = "normal_fallopian_tube_reference",
  GSE154763 = "author_annotated_myeloid_reference"
)

rank_pct <- function(x) {
  ok <- is.finite(x)
  ans <- rep(NA_real_, length(x))
  ans[ok] <- frank(x[ok], ties.method = "average") / sum(ok)
  ans
}

map_status <- function(a, c, cluster_col = "final_cluster") {
  idx <- match(as.character(a[[cluster_col]]), as.character(c$cluster))
  a[, `:=`(
    annotation_status = c$annotation_status[idx],
    canonical_support_n = as.numeric(c$canonical_support_n[idx]),
    incompatible_lineage_program = as.logical(c$incompatible_lineage_program[idx])
  )]
  a
}

raw_cell_state <- function(dataset_id, counts, assign, eligible_cells,
                           previous_scores) {
  cells <- intersect(eligible_cells, colnames(counts))
  assign <- assign[match(cells, cell_id)]
  counts <- counts[, cells, drop = FALSE]
  get_detected <- function(gene) {
    if (gene %in% rownames(counts)) as.numeric(counts[gene, ] > 0) else rep(0, length(cells))
  }
  det <- lapply(c(
    "SPP1", "APOC1", "GPNMB", "TREM2", "LPL", "CTSD",
    "C1QA", "C1QB", "C1QC", "APOE", "MRC1", "SELENOP",
    "FOLR2", "LYVE1", "CD163"
  ), get_detected)
  names(det) <- c(
    "SPP1", "APOC1", "GPNMB", "TREM2", "LPL", "CTSD",
    "C1QA", "C1QB", "C1QC", "APOE", "MRC1", "SELENOP",
    "FOLR2", "LYVE1", "CD163"
  )
  spp1_companion_n <- Reduce(`+`, det[c("APOC1", "GPNMB", "TREM2", "LPL", "CTSD")])
  c1q_detected_n <- Reduce(`+`, det[c("C1QA", "C1QB", "C1QC")])
  c1qc_companion_n <- Reduce(`+`, det[c("APOE", "MRC1", "SELENOP")])
  folr2_companion_n <- Reduce(`+`, det[c("MRC1", "SELENOP", "LYVE1", "CD163")])
  lib <- pmax(Matrix::colSums(counts), 1)
  spp1_expr <- if ("SPP1" %in% rownames(counts)) log1p(as.numeric(counts["SPP1", ]) * 1e4 / lib) else rep(NA_real_, length(cells))
  prev <- previous_scores[match(cells, cell_id)]
  data.table(
    dataset_id = dataset_id,
    dataset_role = roles[[dataset_id]],
    cell_id = cells,
    patient_id = as.character(assign$patient_id),
    sample_id = as.character(assign$sample_id),
    expression_basis = "raw_UMI_counts",
    SPP1_expression = spp1_expr,
    SPP1_detected = det$SPP1 > 0,
    SPP1_companion_detected_n = spp1_companion_n,
    SPP1_program_cell = det$SPP1 > 0 & spp1_companion_n >= 1,
    SPP1_program_strict_cell = det$SPP1 > 0 & spp1_companion_n >= 2,
    SPP1_program_percentile = prev$SPP1_dataset_percentile,
    C1Q_core_cell = c1q_detected_n >= 2,
    C1QC_program_cell = c1q_detected_n >= 2 & c1qc_companion_n >= 1,
    C1QC_program_percentile = prev$C1QC_dataset_percentile,
    FOLR2_detected = det$FOLR2 > 0,
    FOLR2_program_cell = det$FOLR2 > 0 & folr2_companion_n >= 1,
    FOLR2_program_percentile = prev$FOLR2_dataset_percentile
  )
}

previous_scores <- fread(file.path(v4, "macrophage_state_cell_scores.csv.gz"))
cell_parts <- list()

a <- fread(file.path(v4, "GSE147082_refined", "refined_cell_assignments.csv.gz"))
c <- fread(file.path(v4, "GSE147082_refined", "refined_cluster_annotation.csv"))
a <- map_status(a, c)
obj <- readRDS(file.path(data_root, "GSE147082", "objects", "GSE147082_preannotation.rds"))
cnt <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
eligible <- a[
  final_cell_type == "Macrophage" &
    annotation_status %in% c("READY_HIGH_CONFIDENCE", "READY_BROAD_TYPE_ONLY", "REVIEW_PATIENT_ENRICHED") &
    canonical_support_n >= 3 & incompatible_lineage_program != TRUE,
  cell_id
]
cell_parts[["GSE147082"]] <- raw_cell_state(
  "GSE147082", cnt, a, eligible, previous_scores[dataset_id == "GSE147082"]
)
rm(obj, cnt, a, c); gc()

a <- fread(file.path(v3, "GSE151214", "normal_reference_cell_assignments.csv.gz"))
sub <- fread(file.path(v4, "GSE151214_refined", "myeloid_cluster18_subclustering.csv"))
a[sub, on = "cell_id", final_cell_type := i.refined_type]
obj <- readRDS(file.path(data_root, "GSE151214", "objects", "GSE151214_preannotation.rds"))
cnt <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
eligible <- a[final_cell_type == "C1QC_macrophage", cell_id]
cell_parts[["GSE151214"]] <- raw_cell_state(
  "GSE151214", cnt, a, eligible, previous_scores[dataset_id == "GSE151214"]
)
rm(obj, cnt, a, sub); gc()

for (ds in c("GSE154600", "GSE158722")) {
  a <- fread(file.path(cleaned, ds, "cleaned_cell_assignments.csv.gz"))
  c <- fread(file.path(cleaned, ds, "cleaned_cluster_annotation_template.csv"))
  a <- map_status(a, c)
  inp <- readRDS(file.path(
    data_root, "diagnostics_v2", "objects", ds, "lineage_inputs",
    "Myeloid_like_strategy_input.rds"
  ))
  allowed <- if (ds == "GSE158722") "READY_HIGH_CONFIDENCE" else
    c("READY_HIGH_CONFIDENCE", "READY_BROAD_TYPE_ONLY", "REVIEW_PATIENT_ENRICHED")
  eligible <- a[
    final_cell_type == "Macrophage" & annotation_status %in% allowed &
      canonical_support_n >= 3 & incompatible_lineage_program != TRUE,
    cell_id
  ]
  cell_parts[[ds]] <- raw_cell_state(
    ds, inp$counts, a, eligible, previous_scores[dataset_id == ds]
  )
  rm(a, c, inp); gc()
}

# GSE154763 is normalized author-annotation-driven evidence and is not
# numerically pooled with raw-count datasets.
x <- fread(file.path(v4, "GSE154763_refined", "author_annotation_with_refined_state.csv.gz"))
x <- x[final_cell_type == "Macrophage"]
x[, `:=`(
  spp1_pct = rank_pct(SPP1_program),
  c1qc_pct = rank_pct(C1QC_program),
  folr2_pct = rank_pct(FOLR2_program)
)]
is_spp1_author <- grepl("SPP1", x$cell_type_original, ignore.case = TRUE)
is_c1qc_author <- grepl("C1Q", x$cell_type_original, ignore.case = TRUE)
is_folr2_author <- grepl("FOLR2", x$cell_type_original, ignore.case = TRUE)
cell_parts[["GSE154763"]] <- data.table(
  dataset_id = "GSE154763",
  dataset_role = roles[["GSE154763"]],
  cell_id = x$cell_id_harmonized,
  patient_id = as.character(x$patient),
  sample_id = as.character(x$library_id),
  expression_basis = "author_normalized_expression",
  SPP1_expression = x$SPP1_expression,
  SPP1_detected = x$SPP1_expression > 0,
  SPP1_companion_detected_n = NA_integer_,
  SPP1_program_cell = is_spp1_author | x$spp1_pct >= .75,
  SPP1_program_strict_cell = is_spp1_author & x$spp1_pct >= .75,
  SPP1_program_percentile = x$spp1_pct,
  C1Q_core_cell = is_c1qc_author | x$c1qc_pct >= .75,
  C1QC_program_cell = is_c1qc_author | x$c1qc_pct >= .75,
  C1QC_program_percentile = x$c1qc_pct,
  FOLR2_detected = is_folr2_author | x$folr2_pct >= .75,
  FOLR2_program_cell = is_folr2_author | x$folr2_pct >= .75,
  FOLR2_program_percentile = x$folr2_pct
)
rm(x); gc()

cells <- rbindlist(cell_parts, fill = TRUE)
cells[, `:=`(
  SPP1_high_macrophage = SPP1_program_percentile >= .75,
  C1QC_high_macrophage = C1QC_program_percentile >= .75,
  FOLR2_high_macrophage = FOLR2_program_percentile >= .75
)]
fwrite(cells, file.path(v6, "state_detection_by_cell.csv.gz"), compress = "gzip", na = "NA")

summarize_state <- function(d, level) {
  by_cols <- c("dataset_id", "dataset_role", "patient_id")
  if (level == "sample") by_cols <- c(by_cols, "sample_id")
  ans <- d[, .(
    n_macrophages = .N,
    expression_basis = expression_basis[1],
    SPP1_average_expression = mean(SPP1_expression, na.rm = TRUE),
    SPP1_positive_fraction = mean(SPP1_detected, na.rm = TRUE),
    SPP1_program_cell_fraction = mean(SPP1_program_cell, na.rm = TRUE),
    SPP1_program_strict_fraction = mean(SPP1_program_strict_cell, na.rm = TRUE),
    SPP1_program_median_percentile = median(SPP1_program_percentile, na.rm = TRUE),
    SPP1_high_macrophage_fraction = mean(SPP1_high_macrophage, na.rm = TRUE),
    C1QC_core_detection_fraction = mean(C1Q_core_cell, na.rm = TRUE),
    C1QC_program_cell_fraction = mean(C1QC_program_cell, na.rm = TRUE),
    C1QC_program_median_percentile = median(C1QC_program_percentile, na.rm = TRUE),
    C1QC_high_macrophage_fraction = mean(C1QC_high_macrophage, na.rm = TRUE),
    FOLR2_positive_fraction = mean(FOLR2_detected, na.rm = TRUE),
    FOLR2_program_cell_fraction = mean(FOLR2_program_cell, na.rm = TRUE),
    FOLR2_program_median_percentile = median(FOLR2_program_percentile, na.rm = TRUE),
    FOLR2_high_macrophage_fraction = mean(FOLR2_high_macrophage, na.rm = TRUE)
  ), by = by_cols]
  ans[, evaluable := n_macrophages >= 20]
  ans[, SPP1_transcript_detection := fcase(
    !evaluable, "NOT_EVALUABLE",
    SPP1_positive_fraction >= .10, "TRANSCRIPT_DETECTED",
    default = "LOW_OR_NOT_DETECTED"
  )]
  ans[, SPP1_program_support := fcase(
    !evaluable, "NOT_EVALUABLE",
    SPP1_program_cell_fraction >= .10, "PROGRAM_SUPPORTED",
    default = "PROGRAM_LOW"
  )]
  ans[, SPP1_joint_status := fcase(
    !evaluable, "NOT_EVALUABLE",
    SPP1_positive_fraction >= .10 & SPP1_program_cell_fraction >= .10,
    "TRANSCRIPT_AND_PROGRAM_PRESENT",
    SPP1_positive_fraction >= .10, "TRANSCRIPT_PRESENT_PROGRAM_LOW",
    default = "LOW_OR_NOT_DETECTED"
  )]
  ans[, SPP1_relative_enrichment := fcase(
    !evaluable, "NOT_EVALUABLE",
    SPP1_program_median_percentile >= .60, "RELATIVELY_ENRICHED",
    default = "NOT_RELATIVELY_ENRICHED"
  )]
  ans[, C1QC_core_detection := fcase(
    !evaluable, "NOT_EVALUABLE",
    C1QC_core_detection_fraction >= .10, "CORE_DETECTED",
    default = "LOW_OR_NOT_DETECTED"
  )]
  ans[, C1QC_program_support := fcase(
    !evaluable, "NOT_EVALUABLE",
    C1QC_program_cell_fraction >= .10, "PROGRAM_SUPPORTED",
    default = "PROGRAM_LOW"
  )]
  ans[, C1QC_relative_enrichment := fcase(
    !evaluable, "NOT_EVALUABLE",
    C1QC_program_median_percentile >= .60, "RELATIVELY_ENRICHED",
    default = "NOT_RELATIVELY_ENRICHED"
  )]
  ans[, FOLR2_transcript_detection := fcase(
    !evaluable, "NOT_EVALUABLE",
    FOLR2_positive_fraction >= .10, "TRANSCRIPT_DETECTED",
    default = "LOW_OR_NOT_DETECTED"
  )]
  ans[, FOLR2_program_support := fcase(
    !evaluable, "NOT_EVALUABLE",
    FOLR2_program_cell_fraction >= .10, "PROGRAM_SUPPORTED",
    default = "PROGRAM_LOW"
  )]
  ans[, FOLR2_relative_enrichment := fcase(
    !evaluable, "NOT_EVALUABLE",
    FOLR2_program_median_percentile >= .60, "RELATIVELY_ENRICHED",
    default = "NOT_RELATIVELY_ENRICHED"
  )]
  ans[]
}

by_sample <- summarize_state(cells, "sample")
by_patient <- summarize_state(cells, "patient")
fwrite(by_sample, out, na = "NA")
fwrite(by_patient, file.path(v6, "state_detection_by_patient.csv"), na = "NA")

thresholds <- CJ(
  minimum_macrophages = c(10L, 20L, 30L),
  transcript_positive_fraction = c(.05, .10, .20),
  program_cell_fraction = c(.05, .10, .20)
)
sensitivity <- rbindlist(lapply(seq_len(nrow(thresholds)), function(i) {
  th <- thresholds[i]
  p <- cells[, .(
    n_macrophages = .N,
    positive_fraction = mean(SPP1_detected),
    program_fraction = mean(SPP1_program_cell)
  ), by = .(dataset_id, dataset_role, patient_id)]
  p[, evaluable := n_macrophages >= th$minimum_macrophages]
  p[, present := evaluable &
      positive_fraction >= th$transcript_positive_fraction &
      program_fraction >= th$program_cell_fraction]
  p[, .(
    n_evaluable_patients = sum(evaluable),
    n_present_patients = sum(present),
    conclusion = if (sum(evaluable) == 0) "NOT_EVALUABLE" else
      if (sum(present) >= 2) "REPLICATED" else
        if (sum(present) == 1) "SUPPORTIVE_SINGLE_PATIENT" else "NOT_SUPPORTED"
  ), by = .(dataset_id, dataset_role)][, `:=`(
    minimum_macrophages = th$minimum_macrophages,
    transcript_positive_fraction = th$transcript_positive_fraction,
    program_cell_fraction = th$program_cell_fraction
  )]
}))
robust <- sensitivity[, {
  main <- conclusion[
    minimum_macrophages == 20 &
      transcript_positive_fraction == .10 &
      program_cell_fraction == .10
  ][1]
  cls <- if (is.na(main) || main == "NOT_EVALUABLE") "not_evaluable" else
    if (main == "REPLICATED" && all(conclusion == "REPLICATED")) "stable_replicated" else
      if (uniqueN(conclusion) > 1) "threshold_sensitive" else "not_evaluable"
  list(main_conclusion = main, SPP1_threshold_robustness = cls)
}, by = .(dataset_id, dataset_role)]
sensitivity <- merge(sensitivity, robust, by = c("dataset_id", "dataset_role"), all.x = TRUE)
fwrite(sensitivity, file.path(v6, "state_threshold_sensitivity.csv"), na = "NA")

p1 <- ggplot(
  by_sample,
  aes(SPP1_positive_fraction, SPP1_program_cell_fraction,
      color = SPP1_joint_status, size = pmin(n_macrophages, 200))
) +
  geom_vline(xintercept = .10, linetype = 2, color = "grey55") +
  geom_hline(yintercept = .10, linetype = 2, color = "grey55") +
  geom_point(alpha = .85) +
  facet_wrap(~dataset_id) +
  theme_bw() +
  labs(
    title = "SPP1 transcript detection versus companion-program support",
    x = "SPP1-positive macrophage fraction",
    y = "SPP1 + companion program-cell fraction",
    size = "macrophages"
  )
ggsave(file.path(v6, "01_spp1_transcript_vs_program.png"), p1, width = 11, height = 6, dpi = 180)

p2dt <- sensitivity[minimum_macrophages == 20]
p2dt[, threshold_pair := sprintf("%.2f / %.2f", transcript_positive_fraction, program_cell_fraction)]
p2 <- ggplot(p2dt, aes(threshold_pair, dataset_id, fill = conclusion)) +
  geom_tile(color = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(
    title = "SPP1 state threshold sensitivity (minimum 20 macrophages)",
    x = "transcript / program fraction threshold", y = NULL
  )
ggsave(file.path(v6, "02_state_threshold_sensitivity.png"), p2, width = 10, height = 5.5, dpi = 180)
message("Final count-based state detection complete")
