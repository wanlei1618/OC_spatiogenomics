#!/usr/bin/env Rscript

# Reanalyse GSE147082 from the 6,993 cells retained by repaired mitochondrial QC.

options(stringsAsFactors = FALSE, warn = 1)

required <- c("yaml", "data.table", "Matrix", "Seurat", "SeuratObject", "ggplot2")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else "."
source(file.path(script_dir, "_diagnostics_v2_common.R"))

z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
output_root <- file.path(data_root, "diagnostics_v3_remaining_datasets")
out <- file.path(output_root, "GSE147082")
if (dir.exists(out)) stop("Output already exists; refusing to overwrite: ", out)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

seed <- as.integer(cfg$project$random_seed %||% 20260718L)
set.seed(seed)
expected_cells <- 6993L

gene_upper <- function(x) toupper(as.character(x))
write_gz <- function(x, path) data.table::fwrite(x, path, na = "NA", compress = "gzip")

marker_panels <- list(
  Epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MSLN", "WFDC2", "MUC1", "MUC16", "PAX8", "CLDN3", "CLDN4"),
  T_cell = c("CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2", "CD2", "CD247", "LCK", "ITK", "BCL11B", "IL7R"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "KLRF1", "KLRC1", "PRF1", "CTSW", "XCL1", "XCL2", "FGFBP2"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22", "CD37", "CD74", "BANK1", "HLA-DRA"),
  Plasma_cell = c("MZB1", "JCHAIN", "SDC1", "DERL3", "XBP1", "SSR4", "FKBP11", "PRDX4", "IGHG1", "IGHA1"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "MRC1", "CD68", "TREM2", "SPP1", "FCGR1A", "MS4A7", "LPL", "CTSD"),
  Monocyte = c("S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "LILRB1", "LYZ", "SAT1", "LGALS3", "TYROBP"),
  pDC = c("LILRA4", "CLEC4C", "IL3RA", "GZMB", "TCF4", "SERPINF1"),
  cDC1 = c("XCR1", "CLEC9A", "BATF3", "CADM1", "IRF8", "CST3"),
  cDC2 = c("CD1C", "FCER1A", "CLEC10A", "CD1E", "HLA-DRA", "CST3"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HPGDS", "HDC"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "C7", "COL6A1", "COL6A2", "PDGFRA", "FAP", "POSTN", "CTHRC1"),
  Pericyte = c("RGS5", "MCAM", "CSPG4", "PDGFRB", "ACTA2", "MYH11", "NOTCH3", "KCNJ8", "ABCC9", "RBP1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "RAMP3", "PLVAP", "CLDN5", "ACKR1", "CA4", "RGCC", "ESAM")
)

high_specific <- list(
  Epithelial = c("EPCAM", "PAX8", "WFDC2", "MUC16"),
  T_cell = c("CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2"),
  NK_cell = c("GNLY", "KLRD1", "KLRF1", "XCL1", "XCL2"),
  B_cell = c("MS4A1", "CD19", "CD79A", "CD79B"),
  Plasma_cell = c("MZB1", "JCHAIN", "SDC1", "DERL3"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "MRC1", "TREM2"),
  Monocyte = c("S100A8", "S100A9", "FCN1", "VCAN"),
  pDC = c("LILRA4", "CLEC4C", "IL3RA", "GZMB"),
  cDC1 = c("XCR1", "CLEC9A", "BATF3", "CADM1"),
  cDC2 = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "HDC"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
  Pericyte = c("RGS5", "CSPG4", "PDGFRB", "KCNJ8"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "CLDN5")
)

state_panels <- list(
  Cycling = c("MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "STMN1", "PCNA", "CDC20", "CDK1", "CCNB1"),
  IFN_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "OASL", "IRF7", "IFI6"),
  Hypoxia = c("HIF1A", "CA9", "VEGFA", "BNIP3", "NDRG1", "EGLN3", "LDHA"),
  Stress_response = c("FOS", "JUN", "JUNB", "DDIT3", "HSPA1A", "HSPA1B", "ATF3"),
  SPP1_program = c("SPP1", "APOC1", "GPNMB", "TREM2", "LPL", "CTSD"),
  C1QC_program = c("C1QA", "C1QB", "C1QC", "APOE", "MRC1", "SELENOP"),
  FOLR2_program = c("FOLR2", "MRC1", "SELENOP", "C1QC", "LYVE1", "CD163")
)

type_family <- function(x) {
  out <- as.character(x)
  out[out %in% c("T_cell", "NK_cell")] <- "T_NK"
  out[out %in% c("B_cell", "Plasma_cell")] <- "B_Plasma"
  out[out %in% c("Macrophage", "Monocyte", "pDC", "cDC1", "cDC2")] <- "Myeloid_DC"
  out[out %in% c("Fibroblast", "Pericyte")] <- "Stromal"
  out
}

add_scores <- function(obj, panels, prefix) {
  panels <- lapply(panels, function(x) intersect(x, rownames(obj)))
  keep <- lengths(panels) >= 2L
  if (!all(keep)) stop("Insufficient genes for score panels: ", paste(names(panels)[!keep], collapse = ", "))
  obj <- AddModuleScore(obj, features = unname(panels), name = prefix,
                        nbin = 24, ctrl = 25, seed = seed, search = FALSE)
  score_cols <- paste0(prefix, seq_along(panels))
  names(score_cols) <- names(panels)
  list(object = obj, columns = score_cols)
}

score_clusters <- function(markers, clusters, score_means) {
  rows <- list()
  for (cl in clusters) {
    cm <- markers[as.character(cluster) == as.character(cl)]
    for (cell_type in names(marker_panels)) {
      hit <- cm[gene_upper %in% gene_upper(marker_panels[[cell_type]])]
      hi <- hit[gene_upper %in% gene_upper(high_specific[[cell_type]])]
      rows[[length(rows) + 1L]] <- data.table::data.table(
        cluster = as.character(cl), candidate = cell_type,
        support_n = data.table::uniqueN(hit$gene_upper),
        high_specific_n = data.table::uniqueN(hi$gene_upper),
        module_mean = score_means[cluster == as.character(cl), get(cell_type)],
        markers = if (nrow(hit)) paste(unique(hit[order(-avg_log2FC)]$gene), collapse = ";") else ""
      )
    }
  }
  scores <- data.table::rbindlist(rows)
  scores[, eligible := support_n >= 3L | high_specific_n >= 2L]
  decisions <- scores[order(cluster, -eligible, -support_n, -high_specific_n,
                           -module_mean), {
    top <- .SD[1L]
    eligible_rows <- .SD[eligible == TRUE]
    second <- if (nrow(eligible_rows) >= 2L) eligible_rows[2L] else NULL
    conflict <- !is.null(second) &&
      type_family(top$candidate) != type_family(second$candidate) &&
      second$support_n >= top$support_n - 1L &&
      second$module_mean >= top$module_mean - 0.10
    list(
      final_cell_type = if (!top$eligible || conflict) "Unresolved" else as.character(top$candidate),
      suggested_cell_type = as.character(top$candidate),
      canonical_support_n = as.integer(top$support_n),
      high_specific_support_n = as.integer(top$high_specific_n),
      canonical_markers = as.character(top$markers),
      second_candidate = if (is.null(second)) NA_character_ else as.character(second$candidate),
      incompatible_lineage_program = conflict
    )
  }, by = cluster]
  list(scores = scores, decisions = decisions)
}

object_path <- file.path(data_root, "GSE147082", "objects", "GSE147082_preannotation.rds")
qc_path <- file.path(data_root, "diagnostics_v2", "GSE147082", "01_mt_audit",
                     "qc_metadata_mt_repaired.csv.gz")
old <- readRDS(object_path)
qc <- data.table::fread(qc_path, showProgress = FALSE)
if (!"cell_id" %in% names(qc)) qc[, cell_id := paste(sample_id, original_cell_id, sep = "__")]
qc <- qc[diagnostics_v2_qc_pass == TRUE]
if (nrow(qc) != expected_cells || data.table::uniqueN(qc$cell_id) != expected_cells) {
  stop("Expected 6,993 unique repaired-QC cells")
}
if (!all(qc$cell_id %in% colnames(old))) stop("Repaired-QC cell IDs do not match count object")

old_md <- old[[]]
old_cluster_col <- if ("RNA_snn_res.0.6" %in% names(old_md)) "RNA_snn_res.0.6" else "seurat_clusters"
old_cluster <- as.character(old_md[qc$cell_id, old_cluster_col])
counts <- SeuratObject::LayerData(old, assay = "RNA", layer = "counts")[, qc$cell_id, drop = FALSE]
md <- as.data.frame(qc)
rownames(md) <- md$cell_id
md$old_cluster <- old_cluster

genes_upper <- gene_upper(rownames(counts))
detected <- matrix(0L, nrow = ncol(counts), ncol = length(high_specific),
                   dimnames = list(colnames(counts), names(high_specific)))
for (cell_type in names(high_specific)) {
  idx <- which(genes_upper %in% gene_upper(high_specific[[cell_type]]))
  if (length(idx)) detected[, cell_type] <- Matrix::colSums(counts[idx, , drop = FALSE] > 0)
}
eligible <- detected >= 2L
families_per_cell <- vapply(seq_len(nrow(eligible)), function(i) {
  data.table::uniqueN(type_family(colnames(eligible)[eligible[i, ]]))
}, FUN.VALUE = integer(1))
incompatible_cell_program <- families_per_cell >= 2L
dbl_call <- tolower(as.character(md$scDblFinder.class))
remove_doublet <- !is.na(dbl_call) & dbl_call == "doublet" & incompatible_cell_program

n_feature <- as.numeric(md$nFeature_RNA)
n_count <- as.numeric(md$nCount_RNA)
percent_mt <- as.numeric(md$percent.mt)
complexity <- log10(pmax(n_feature, 1)) / log10(pmax(n_count, 10))
stress_gene <- grepl("^MT[-.]|^RPL|^RPS|^HSPA|^FOS$|^JUN$|^JUNB$|^ATF3$", genes_upper)
stress_fraction <- Matrix::colSums(counts[stress_gene, , drop = FALSE]) /
  pmax(Matrix::colSums(counts), 1)
sample_groups <- split(seq_len(nrow(md)), as.character(md$sample_id))
low_feature <- low_count <- low_complexity <- high_mt <- rep(FALSE, nrow(md))
for (idx in sample_groups) {
  low_feature[idx] <- n_feature[idx] <= max(200, quantile(n_feature[idx], .01, na.rm = TRUE))
  low_count[idx] <- n_count[idx] <= max(500, quantile(n_count[idx], .01, na.rm = TRUE))
  low_complexity[idx] <- complexity[idx] <= quantile(complexity[idx], .01, na.rm = TRUE)
  high_mt[idx] <- percent_mt[idx] >= max(20, quantile(percent_mt[idx], .99, na.rm = TRUE))
}
no_cell_lineage <- rowSums(eligible) == 0L
remove_low_quality <- no_cell_lineage & stress_fraction >= .70 &
  (low_feature + low_count + low_complexity + high_mt >= 3L)
removal_reason <- data.table::fcase(
  remove_doublet, "REMOVE_HETEROTYPIC_DOUBLET",
  remove_low_quality, "REMOVE_LOW_QUALITY",
  default = NA_character_
)
retained <- is.na(removal_reason)

removed <- data.table::data.table(
  dataset_id = "GSE147082", cell_id = md$cell_id[!retained],
  sample_id = as.character(md$sample_id[!retained]),
  patient_id = as.character(md$patient_id[!retained]),
  old_cluster = md$old_cluster[!retained], final_cluster = NA_character_,
  final_cell_type = NA_character_, cell_subtype = NA_character_,
  cell_state = NA_character_, annotation_confidence = "REMOVED",
  patient_enriched = FALSE, doublet_call = as.character(md$scDblFinder.class[!retained]),
  removal_reason = removal_reason[!retained]
)

counts <- counts[, retained, drop = FALSE]
md <- md[retained, , drop = FALSE]
rm(old, old_md)
gc()

obj <- CreateSeuratObject(counts = counts, project = "GSE147082_mt_repaired",
                          min.cells = 0, min.features = 0, meta.data = md)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000,
                            verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, seed.use = seed, verbose = FALSE)
obj <- FindClusters(obj, resolution = c(.4, .6, .8), random.seed = seed,
                    verbose = FALSE)
cluster_col <- "RNA_snn_res.0.6"
if (!cluster_col %in% colnames(obj[[]])) stop("Resolution 0.6 cluster column missing")
Idents(obj) <- factor(obj[[cluster_col, drop = TRUE]])
obj$final_cluster <- as.character(Idents(obj))

markers <- FindAllMarkers(
  obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
  min.pct = .20, logfc.threshold = .25, return.thresh = 1, verbose = FALSE
)
markers <- data.table::as.data.table(markers, keep.rownames = FALSE)
if (!"gene" %in% names(markers)) markers[, gene := rownames(markers)]
if (!"avg_log2FC" %in% names(markers) && "avg_logFC" %in% names(markers)) {
  data.table::setnames(markers, "avg_logFC", "avg_log2FC")
}
markers[, gene_upper := gene_upper(gene)]
markers <- markers[is.finite(p_val_adj) & p_val_adj < .05 &
                     is.finite(avg_log2FC) & avg_log2FC > .25 &
                     is.finite(pct.1) & pct.1 >= .20]
markers[, dataset_id := "GSE147082"]

lineage_scored <- add_scores(obj, marker_panels, "lineage_module")
obj <- lineage_scored$object
state_scored <- add_scores(obj, state_panels, "state_module")
obj <- state_scored$object

score_md <- data.table::as.data.table(obj[[]], keep.rownames = "cell_id")
score_md[, cluster := as.character(get(cluster_col))]
lineage_means <- score_md[, lapply(.SD, mean), by = cluster,
                          .SDcols = unname(lineage_scored$columns)]
data.table::setnames(lineage_means, unname(lineage_scored$columns), names(lineage_scored$columns))
evidence <- score_clusters(markers, sort(unique(score_md$cluster)), lineage_means)

dominance <- dominance_metrics(score_md$cluster, score_md$sample_id,
                               "GSE147082", "A_uncorrected")
data.table::setnames(dominance, "normalized_shannon_entropy", "sample_entropy")
dominance[, patient_enriched := dominant_sample_fraction >= .80]

state_means <- score_md[, lapply(.SD, mean), by = cluster,
                        .SDcols = unname(state_scored$columns)]
data.table::setnames(state_means, unname(state_scored$columns), names(state_scored$columns))

data_mat <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
cluster_cells <- split(colnames(obj), as.character(obj[[cluster_col, drop = TRUE]]))
program_summary <- function(cluster, genes, score_col) {
  cells <- cluster_cells[[as.character(cluster)]]
  present <- intersect(genes, rownames(data_mat))
  mat <- data_mat[present, cells, drop = FALSE]
  spp1_present <- "SPP1" %in% rownames(data_mat)
  spp1_values <- if (spp1_present) data_mat["SPP1", cells] else rep(0, length(cells))
  sample <- as.character(obj$sample_id[match(cells, colnames(obj))])
  positive_by_sample <- data.table::data.table(sample = sample, positive = spp1_values > 0)[
    , .(n = .N, n_positive = sum(positive), positive_fraction = mean(positive)), by = sample
  ]
  list(
    module_mean = mean(obj[[score_col, drop = TRUE]][match(cells, colnames(obj))]),
    spp1_mean = mean(spp1_values), spp1_positive_fraction = mean(spp1_values > 0),
    program_genes_detected = sum(Matrix::rowSums(mat > 0) > 0),
    n_positive_samples = sum(positive_by_sample$n_positive >= 5 &
                               positive_by_sample$positive_fraction >= .05),
    positive_samples = paste(positive_by_sample[n_positive >= 5 & positive_fraction >= .05]$sample,
                             collapse = ";")
  )
}

decisions <- evidence$decisions
spp1_rows <- lapply(decisions[final_cell_type == "Macrophage"]$cluster, function(cl) {
  x <- program_summary(cl, state_panels$SPP1_program, state_scored$columns[["SPP1_program"]])
  data.table::data.table(cluster = cl, n_cells = length(cluster_cells[[cl]]),
                         n_samples = data.table::uniqueN(obj$sample_id[match(cluster_cells[[cl]], colnames(obj))]),
                         module_mean = x$module_mean, spp1_mean = x$spp1_mean,
                         spp1_positive_fraction = x$spp1_positive_fraction,
                         program_genes_detected = x$program_genes_detected,
                         n_positive_samples = x$n_positive_samples,
                         positive_samples = x$positive_samples,
                         cross_sample_replicated = x$n_positive_samples >= 2L,
                         spp1_high = x$module_mean > 0 & x$spp1_positive_fraction >= .10 &
                           x$program_genes_detected >= 3L & x$n_positive_samples >= 2L)
})
spp1_summary <- if (length(spp1_rows)) data.table::rbindlist(spp1_rows) else
  data.table::data.table(cluster = character(), n_cells = integer(), n_samples = integer(),
                         module_mean = numeric(), spp1_mean = numeric(),
                         spp1_positive_fraction = numeric(), program_genes_detected = integer(),
                         n_positive_samples = integer(), positive_samples = character(),
                         cross_sample_replicated = logical(), spp1_high = logical())

top_markers <- markers[order(cluster, -avg_log2FC),
                       .(top20_markers = paste(head(gene, 20L), collapse = ";")), by = cluster]
marker_counts <- markers[, .N, by = cluster]
cluster_table <- merge(decisions, dominance[, .(
  cluster, n_cells, n_samples, dominant_sample, dominant_sample_fraction,
  sample_entropy, patient_enriched
)], by = "cluster", all.x = TRUE)
cluster_table <- merge(cluster_table, state_means, by = "cluster", all.x = TRUE)
cluster_table <- merge(cluster_table, top_markers, by = "cluster", all.x = TRUE)
cluster_table[, n_significant_markers := marker_counts$N[match(cluster, marker_counts$cluster)]]
cluster_table[is.na(n_significant_markers), n_significant_markers := 0L]

state_marker_count <- function(cl, panel) {
  data.table::uniqueN(markers[cluster == cl & gene_upper %in% gene_upper(panel)]$gene_upper)
}
cluster_table[, cell_state := vapply(cluster, function(cl) {
  sm <- state_means[cluster == cl]
  cluster_n_samples <- cluster_table[cluster == cl, n_samples][[1L]]
  states <- character()
  if (sm$Cycling > 0 && state_marker_count(cl, state_panels$Cycling) >= 3L) states <- c(states, "Cycling")
  if (sm$IFN_response > 0 && state_marker_count(cl, state_panels$IFN_response) >= 3L) states <- c(states, "IFN_response")
  if (sm$Hypoxia > 0 && state_marker_count(cl, state_panels$Hypoxia) >= 3L) states <- c(states, "Hypoxia")
  if (sm$Stress_response > 0 && state_marker_count(cl, state_panels$Stress_response) >= 3L) states <- c(states, "Stress_response")
  if (nrow(spp1_summary[cluster == cl & spp1_high == TRUE])) states <- c(states, "SPP1_high")
  if (sm$C1QC_program > 0 && state_marker_count(cl, state_panels$C1QC_program) >= 3L && cluster_n_samples >= 2L) states <- c(states, "C1QC_high")
  if (sm$FOLR2_program > 0 && state_marker_count(cl, state_panels$FOLR2_program) >= 3L && cluster_n_samples >= 2L) states <- c(states, "FOLR2_high")
  if (length(states)) paste(states, collapse = ";") else "None"
}, FUN.VALUE = character(1))]
cluster_table[final_cell_type == "Unresolved" & grepl("(^|;)Cycling(;|$)", cell_state),
              final_cell_type := "Cycling_lineage_unknown"]
cluster_table[, annotation_status := data.table::fcase(
  incompatible_lineage_program == TRUE, "REVIEW_MIXED_OR_DOUBLET",
  final_cell_type == "Unresolved", "UNRESOLVED",
  grepl("Stress_response", cell_state) & canonical_support_n < 3L, "REVIEW_STRESS_DOMINATED",
  patient_enriched == TRUE, "REVIEW_PATIENT_ENRICHED",
  canonical_support_n >= 4L | high_specific_support_n >= 2L, "READY_HIGH_CONFIDENCE",
  canonical_support_n >= 2L, "READY_BROAD_TYPE_ONLY",
  default = "REVIEW_AMBIGUOUS"
)]
cluster_table[, annotation_confidence := data.table::fcase(
  annotation_status == "READY_HIGH_CONFIDENCE", "High",
  annotation_status == "READY_BROAD_TYPE_ONLY", "Broad_type_only",
  final_cell_type == "Unresolved", "Unresolved",
  default = "Review"
)]
cluster_table[, `:=`(
  dataset_id = "GSE147082", cell_subtype = "", notes = "New clustering from repaired mitochondrial-QC cells; old cluster IDs were not used for annotation."
)]
cluster_table <- cluster_table[, .(
  dataset_id, cluster, final_cell_type, cell_subtype, cell_state,
  annotation_status, annotation_confidence, n_cells, n_samples,
  dominant_sample, dominant_sample_fraction, sample_entropy, patient_enriched,
  canonical_markers, top20_markers, canonical_support_n,
  high_specific_support_n, second_candidate, incompatible_lineage_program,
  n_significant_markers, notes
)]

cidx <- match(as.character(obj[[cluster_col, drop = TRUE]]), cluster_table$cluster)
assignments <- data.table::data.table(
  dataset_id = "GSE147082", cell_id = colnames(obj),
  sample_id = as.character(obj$sample_id), patient_id = as.character(obj$patient_id),
  old_cluster = as.character(obj$old_cluster),
  final_cluster = as.character(obj[[cluster_col, drop = TRUE]]),
  final_cell_type = cluster_table$final_cell_type[cidx], cell_subtype = "",
  cell_state = cluster_table$cell_state[cidx],
  annotation_confidence = cluster_table$annotation_confidence[cidx],
  patient_enriched = cluster_table$patient_enriched[cidx],
  doublet_call = as.character(obj$scDblFinder.class), removal_reason = NA_character_
)
if (nrow(assignments) + nrow(removed) != expected_cells ||
    data.table::uniqueN(c(assignments$cell_id, removed$cell_id)) != expected_cells) {
  stop("Retained plus removed cells do not cover repaired-QC input")
}
if (anyNA(assignments$final_cell_type) || any(!nzchar(assignments$final_cell_type))) {
  stop("Retained cell missing final_cell_type")
}

data.table::fwrite(cluster_table, file.path(out, "cleaned_cluster_annotation_template.csv"), na = "NA")
write_gz(assignments, file.path(out, "cleaned_cell_assignments.csv.gz"))
write_gz(removed, file.path(out, "removed_cells.csv.gz"))
write_gz(markers, file.path(out, "significant_markers.csv.gz"))
data.table::fwrite(dominance, file.path(out, "cluster_sample_dominance.csv"), na = "NA")
data.table::fwrite(spp1_summary, file.path(out, "SPP1_macrophage_state_summary.csv"), na = "NA")

p_type <- DimPlot(obj, reduction = "umap", group.by = "final_cluster", label = FALSE,
                  raster = TRUE) + ggplot2::labs(title = "GSE147082 repaired-QC clusters")
obj$final_cell_type <- assignments$final_cell_type[match(colnames(obj), assignments$cell_id)]
p_type <- DimPlot(obj, reduction = "umap", group.by = "final_cell_type", label = TRUE,
                  repel = TRUE, raster = TRUE) + ggplot2::labs(title = "GSE147082 final cell type")
p_sample <- DimPlot(obj, reduction = "umap", group.by = "sample_id", raster = TRUE) +
  ggplot2::labs(title = "GSE147082 by sample")
ggplot2::ggsave(file.path(out, "UMAP_final_cell_type.png"), p_type,
                width = 9, height = 7, dpi = 180)
ggplot2::ggsave(file.path(out, "UMAP_by_sample.png"), p_sample,
                width = 9, height = 7, dpi = 180)
dot_genes <- unique(unlist(lapply(marker_panels, head, 3L), use.names = FALSE))
dot_genes <- intersect(dot_genes, rownames(obj))
p_dot <- DotPlot(obj, features = dot_genes, group.by = "final_cell_type") +
  RotatedAxis() + ggplot2::labs(title = "GSE147082 canonical marker dotplot")
ggplot2::ggsave(file.path(out, "marker_dotplot.png"), p_dot,
                width = 15, height = 8, dpi = 180)

type_counts <- assignments[, .N, by = final_cell_type][order(final_cell_type)]
removal_counts <- table(factor(removed$removal_reason,
                               levels = c("REMOVE_HETEROTYPIC_DOUBLET", "REMOVE_LOW_QUALITY")))
summary_lines <- c(
  "# GSE147082 repaired-mitochondrial-QC reanalysis", "",
  paste0("- repaired_qc_input_cells: ", expected_cells),
  paste0("- final_retained_cells: ", nrow(assignments)),
  paste0("- removed_heterotypic_doublets: ", removal_counts[["REMOVE_HETEROTYPIC_DOUBLET"]]),
  paste0("- removed_low_quality: ", removal_counts[["REMOVE_LOW_QUALITY"]]),
  "- clustering: newly recomputed from raw counts at resolutions 0.4, 0.6, and 0.8; resolution 0.6 used for the primary result.",
  "- old_7645_cell_clusters: comparison only; neither old marker calls nor old cluster IDs were used to assign final labels.",
  paste0("- patient_enriched_clusters: ", sum(cluster_table$patient_enriched)),
  paste0("- unresolved_cells: ", sum(assignments$final_cell_type %in% c("Unresolved", "Cycling_lineage_unknown"))),
  "", "## Final cell-type counts", "",
  paste0("- ", type_counts$final_cell_type, ": ", type_counts$N),
  "", "## SPP1 macrophage reference", "",
  if (nrow(spp1_summary)) paste0(
    "- cluster ", spp1_summary$cluster, ": n_cells=", spp1_summary$n_cells,
    ", module_mean=", sprintf("%.4f", spp1_summary$module_mean),
    ", SPP1_positive_fraction=", sprintf("%.4f", spp1_summary$spp1_positive_fraction),
    ", positive_samples=", ifelse(nzchar(spp1_summary$positive_samples), spp1_summary$positive_samples, "None"),
    ", cross_sample_replicated=", spp1_summary$cross_sample_replicated,
    ", SPP1_high=", spp1_summary$spp1_high
  ) else "- No macrophage cluster passed the marker-based lineage criteria.",
  "", "Patient/sample enrichment was recorded but was never a deletion criterion.",
  "Cycling, IFN, hypoxia, stress, SPP1, C1QC and FOLR2 programs are stored only in cell_state."
)
writeLines(summary_lines, file.path(out, "analysis_summary.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("GSE147082 complete: ", out)
