options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v4 <- file.path(data_root, "diagnostics_v4_cross_dataset_validation")
v5 <- file.path(data_root, "diagnostics_v5_final_calibration")
v6 <- file.path(data_root, "diagnostics_v6_malignant_receiver_validation")
dir.create(v6, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v6, "GSE147082_cluster6_final_cell_annotation.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

seed <- as.integer(cfg$project$random_seed %||% 20260718L)
set.seed(seed)
a <- fread(file.path(v4, "GSE147082_refined", "refined_cell_assignments.csv.gz"))
obj <- readRDS(file.path(data_root, "GSE147082", "objects", "GSE147082_preannotation.rds"))
counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
c6cells <- a[final_cluster == 6, cell_id]
s <- CreateSeuratObject(counts[, c6cells, drop = FALSE], min.cells = 0, min.features = 0)
s <- NormalizeData(s, verbose = FALSE)
s <- FindVariableFeatures(s, nfeatures = 2000, verbose = FALSE)
s <- ScaleData(s, features = VariableFeatures(s), verbose = FALSE)
s <- RunPCA(s, npcs = 30, verbose = FALSE)
s <- FindNeighbors(s, dims = 1:20, verbose = FALSE)
s <- FindClusters(s, resolution = .25, random.seed = seed, verbose = FALSE)
subcluster <- as.character(Idents(s))
names(subcluster) <- colnames(s)
raw <- SeuratObject::LayerData(s, assay = "RNA", layer = "counts")

any_detected <- function(genes) {
  genes <- intersect(genes, rownames(raw))
  if (!length(genes)) return(rep(FALSE, ncol(raw)))
  Matrix::colSums(raw[genes, , drop = FALSE] > 0) > 0
}
cd3 <- any_detected(c("CD3D", "CD3E"))
tcr <- any_detected(c("TRDC", "TRGC1", "TRGC2"))
nk <- any_detected(c("NCR1", "NCAM1", "FCER1G", "KLRD1"))

cell_ann <- data.table(
  dataset_id = "GSE147082",
  cell_id = colnames(s),
  parent_cluster = "6",
  subcluster = subcluster[colnames(s)],
  patient_id = a$patient_id[match(colnames(s), a$cell_id)],
  sample_id = a$sample_id[match(colnames(s), a$cell_id)],
  CD3D_CD3E_positive = cd3,
  TRDC_TRGC_positive = tcr,
  CD3_TCR_copositive = cd3 & tcr,
  NCR1_NCAM1_FCER1G_KLRD1_positive = nk
)
evidence <- cell_ann[, .(
  n_cells = .N,
  n_patients = uniqueN(patient_id),
  dominant_sample_fraction = max(table(sample_id)) / .N,
  CD3D_CD3E_positive_fraction = mean(CD3D_CD3E_positive),
  TRDC_TRGC_positive_fraction = mean(TRDC_TRGC_positive),
  CD3_TCR_copositive_fraction = mean(CD3_TCR_copositive),
  NCR1_NCAM1_FCER1G_KLRD1_positive_fraction =
    mean(NCR1_NCAM1_FCER1G_KLRD1_positive)
), by = subcluster]

evidence[, `:=`(
  final_cell_type = fcase(
    subcluster == "0", "CD8_effector_T",
    subcluster == "1", "CD8_effector_T",
    CD3_TCR_copositive_fraction >= .50, "Gamma_delta_T",
    CD3_TCR_copositive_fraction < .25 &
      NCR1_NCAM1_FCER1G_KLRD1_positive_fraction >= .50, "NK_like_unresolved",
    default = "Gamma_delta_T_NK_like"
  ),
  cell_state = fcase(
    subcluster == "1", "Cycling",
    subcluster == "2", "NK_like_cytotoxic",
    default = "Cytotoxic"
  ),
  patient_enriched = subcluster == "2",
  annotation_confidence = ifelse(subcluster == "2", "Review", "Medium")
)]
cell_ann <- merge(
  cell_ann,
  evidence[, .(subcluster, final_cell_type, cell_state,
               patient_enriched, annotation_confidence)],
  by = "subcluster", all.x = TRUE
)
setcolorder(cell_ann, c(
  "dataset_id", "cell_id", "parent_cluster", "subcluster",
  "patient_id", "sample_id", "final_cell_type", "cell_state",
  "patient_enriched", "annotation_confidence"
))
fwrite(cell_ann, out, na = "NA")
fwrite(evidence, file.path(v6, "GSE147082_cluster6_tcr_nk_evidence.csv"), na = "NA")

updated <- copy(a)
idx <- match(updated$cell_id, cell_ann$cell_id)
hit <- !is.na(idx)
updated[hit, `:=`(
  final_cell_type = cell_ann$final_cell_type[idx[hit]],
  cell_subtype = cell_ann$final_cell_type[idx[hit]],
  cell_state = cell_ann$cell_state[idx[hit]],
  annotation_confidence = cell_ann$annotation_confidence[idx[hit]],
  patient_enriched = cell_ann$patient_enriched[idx[hit]]
)]
fwrite(
  updated, file.path(v6, "GSE147082_refined_cell_assignments_v3.csv.gz"),
  compress = "gzip", na = "NA"
)

p <- melt(
  evidence,
  id.vars = c("subcluster", "final_cell_type"),
  measure.vars = c(
    "CD3D_CD3E_positive_fraction", "TRDC_TRGC_positive_fraction",
    "CD3_TCR_copositive_fraction",
    "NCR1_NCAM1_FCER1G_KLRD1_positive_fraction"
  ),
  variable.name = "evidence", value.name = "positive_fraction"
)
p <- ggplot(p, aes(evidence, subcluster, fill = positive_fraction)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", positive_fraction)), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "GSE147082 cluster 6 TCR and NK evidence", x = NULL, y = "subcluster")
ggsave(file.path(v6, "03_GSE147082_cluster6_tcr_nk_evidence.png"),
       p, width = 9, height = 4.5, dpi = 180)
message("GSE147082 cluster 6 final annotation complete")
