options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
out_dir <- file.path(data_root, "research_validation_independent_cnv")
run_root <- file.path(out_dir, "infercnv_patient_runs")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(out_dir, "GSE154600_infercnv_validation_by_patient.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

r_error <- file.path(out_dir, "infercnv_R_package_error.txt")
if (!requireNamespace("infercnv", quietly = TRUE)) {
  writeLines(c(
    "Standard inferCNV R package unavailable.",
    "BiocManager installation attempt failed before dependency resolution:",
    "Error: Bioconductor version cannot be validated; no internet connection?",
    "Warning: InternetOpenUrl failed: security channel support error.",
    "Per the task plan, infercnvpy is used as the sole fallback."
  ), r_error)
  method <- "infercnvpy_fallback"
} else {
  method <- "infercnv_R"
}
if (method == "infercnv_R")
  stop("This fixed script implements the task-authorized infercnvpy fallback only; ",
       "remove the fallback after a validated standard inferCNV environment is available")

counts <- readRDS(file.path(
  out_dir, "GSE154600_complete_final_epithelial_counts.rds"
))
stability <- fread(file.path(
  out_dir, "GSE154600_copykat_stability_by_cell_v2.csv.gz"
))
cleaned <- file.path(data_root, "diagnostics_v2_marker_ready_cleaned", "GSE154600")
assignment <- fread(file.path(cleaned, "cleaned_cell_assignments.csv.gz"))
cluster <- fread(file.path(cleaned, "cleaned_cluster_annotation_template.csv"))
idx <- match(as.character(assignment$final_cluster), as.character(cluster$cluster))
assignment[, `:=`(
  annotation_status = cluster$annotation_status[idx],
  canonical_support_n = as.numeric(cluster$canonical_support_n[idx]),
  incompatible_lineage_program =
    as.logical(cluster$incompatible_lineage_program[idx])
)]
reference_types <- c(
  "T_cell", "NK_cell", "B_cell", "Macrophage",
  "Monocyte", "cDC1", "cDC2", "pDC"
)
refs <- assignment[
  final_cell_type %in% reference_types &
    annotation_status == "READY_HIGH_CONFIDENCE" &
    canonical_support_n >= 3 &
    incompatible_lineage_program != TRUE,
  .(patient_id, cell_id, final_cell_type)
]
refs[, reference_lineage := fcase(
  final_cell_type %in% c("T_cell", "NK_cell"), "T_NK",
  final_cell_type == "B_cell", "B",
  default = "Myeloid"
)]

# Use the authoritative full object for same-patient reference counts.
full_object <- readRDS(file.path(
  data_root, "GSE154600", "objects", "GSE154600_preannotation.rds"
))
assay <- full_object@assays[[full_object@active.assay]]
full_counts <- attr(assay, "layers")$counts
dimnames(full_counts) <- list(
  rownames(attr(assay, "features"))[seq_len(nrow(full_counts))],
  rownames(attr(assay, "cells"))[seq_len(ncol(full_counts))]
)
full_counts <- as(full_counts, "dgCMatrix")
rm(full_object, assay)
gc()

genes <- rownames(full_counts)
map_first <- function(column) {
  unname(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
    column = column, multiVals = "first"
  ))
}
chromosome <- map_first("CHR")
start <- suppressWarnings(abs(as.numeric(map_first("CHRLOC"))))
end <- suppressWarnings(abs(as.numeric(map_first("CHRLOCEND"))))
gene_map <- data.table(
  gene = genes,
  chromosome = paste0("chr", chromosome),
  start = pmin(start, end, na.rm = TRUE),
  end = pmax(start, end, na.rm = TRUE)
)
gene_map[
  !chromosome %in% paste0("chr", c(1:22, "X", "Y")) |
    !is.finite(start) | !is.finite(end),
  `:=`(chromosome = NA_character_, start = NA_real_, end = NA_real_)
]
gene_map <- gene_map[!is.na(chromosome)]
gene_map <- gene_map[!duplicated(gene)]
fwrite(gene_map, file.path(run_root, "gene_order_from_org.Hs.eg.db.csv"),
       na = "NA")

python <- "C:/Users/chenfy12/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
python_lib <- file.path(out_dir, "python_lib")
if (!file.exists(python) || !dir.exists(file.path(python_lib, "infercnvpy")))
  stop("Task-authorized infercnvpy fallback environment is unavailable")
Sys.setenv(PYTHONPATH = normalizePath(
  python_lib, winslash = "/", mustWork = TRUE
))
worker <- file.path(run_root, "_infercnvpy_worker.py")
writeLines(c(
  "import sys",
  "import numpy as np",
  "import pandas as pd",
  "import scipy.io",
  "import scipy.sparse as sp",
  "import anndata as ad",
  "import scanpy as sc",
  "import infercnvpy as cnv",
  "import importlib",
  "infercnv_module = importlib.import_module('infercnvpy.tl._infercnv')",
  "infercnv_module.process_map = lambda fn, *iterables, **kwargs: [fn(*args) for args in zip(*iterables)]",
  "mtx_path, cells_path, genes_path, out_path = sys.argv[1:5]",
  "cells = pd.read_csv(cells_path)",
  "genes = pd.read_csv(genes_path)",
  "x = scipy.io.mmread(mtx_path).tocsr().transpose().tocsr()",
  "adata = ad.AnnData(X=x)",
  "adata.obs_names = cells['cell_id'].astype(str).values",
  "adata.obs['reference_group'] = cells['reference_group'].astype(str).values",
  "adata.var_names = genes['gene'].astype(str).values",
  "adata.var['chromosome'] = genes['chromosome'].astype(str).values",
  "adata.var['start'] = genes['start'].astype(float).values",
  "adata.var['end'] = genes['end'].astype(float).values",
  "sc.pp.normalize_total(adata, target_sum=1e4)",
  "sc.pp.log1p(adata)",
  "cnv.tl.infercnv(adata, reference_key='reference_group', reference_cat=['reference'], lfc_clip=3, window_size=100, step=10, dynamic_threshold=1.5, exclude_chromosomes=('chrX','chrY'), chunksize=1000, n_jobs=1, inplace=True, calculate_gene_values=False)",
  "xc = adata.obsm['X_cnv'].tocsr()",
  "scores = np.asarray(abs(xc).mean(axis=1)).ravel()",
  "refmask = (adata.obs['reference_group'].values == 'reference')",
  "ref_scores = scores[refmask]",
  "ref_median = float(np.median(ref_scores))",
  "ref_mad = float(np.median(np.abs(ref_scores - ref_median)))",
  "threshold = ref_median + 3.0 * ref_mad",
  "centromere = {'chr1':123400000,'chr2':93900000,'chr3':90900000,'chr4':50000000,'chr5':48750000,'chr6':60550000,'chr7':60100000,'chr8':45200000,'chr9':43850000,'chr10':39800000,'chr11':53400000,'chr12':35500000,'chr13':17700000,'chr14':17150000,'chr15':19000000,'chr16':36850000,'chr17':25050000,'chr18':18450000,'chr19':26150000,'chr20':28050000,'chr21':11950000,'chr22':15550000}",
  "chr_pos = adata.uns['cnv']['chr_pos']",
  "ordered = sorted(chr_pos.items(), key=lambda kv: kv[1])",
  "arm_labels = []",
  "for ii, (chrom, first) in enumerate(ordered):",
  "    last = ordered[ii+1][1] if ii+1 < len(ordered) else xc.shape[1]",
  "    nwin = last - first",
  "    gv = adata.var.loc[adata.var['chromosome'] == chrom].sort_values('start')",
  "    if len(gv) > 100:",
  "        starts = list(range(0, len(gv) - 100 + 1, 10))",
  "        centers = [float(gv.iloc[s:s+100]['start'].median()) for s in starts]",
  "    else:",
  "        centers = [float(gv['start'].median())]",
  "    if len(centers) != nwin:",
  "        centers = np.linspace(float(gv['start'].min()), float(gv['start'].max()), nwin).tolist()",
  "    arm_labels.extend([chrom + ('p' if pos < centromere[chrom] else 'q') for pos in centers])",
  "arm_labels = np.asarray(arm_labels)",
  "n_high_arms = np.zeros(adata.n_obs, dtype=int)",
  "for arm in np.unique(arm_labels):",
  "    cols = np.where(arm_labels == arm)[0]",
  "    arm_score = np.asarray(abs(xc[:, cols]).mean(axis=1)).ravel()",
  "    arm_threshold = float(np.quantile(arm_score[refmask], 0.95))",
  "    n_high_arms += arm_score > arm_threshold",
  "high = (scores > threshold) & (n_high_arms >= 2)",
  "result = pd.DataFrame({'cell_id': adata.obs_names, 'reference_group': adata.obs['reference_group'].values, 'infercnv_score': scores, 'reference_score_median': ref_median, 'reference_score_mad': ref_mad, 'infercnv_score_threshold': threshold, 'n_high_chromosome_arms': n_high_arms, 'infercnv_status': np.where(high, 'HIGH_CNV_SUPPORT', 'LOW_CNV_SUPPORT')})",
  "result.to_csv(out_path, index=False)"
), worker)

sample_references <- function(pt, seed = 20260724L) {
  x <- refs[patient_id == pt & cell_id %in% colnames(full_counts)]
  limits <- c(T_NK = 150L, B = 100L, Myeloid = 250L)
  set.seed(seed)
  rbindlist(lapply(names(limits), function(lineage) {
    ids <- x[reference_lineage == lineage, cell_id]
    if (length(ids) > limits[[lineage]]) ids <- sample(ids, limits[[lineage]])
    x[cell_id %in% ids]
  }))
}

patients <- c("T59", "T76", "T77", "T89", "T90")
per_patient <- list()
run_audit <- list()
for (pt in patients) {
  target_ids <- stability[patient_id == pt, cell_id]
  ref <- sample_references(pt)
  ids <- unique(c(target_ids, ref$cell_id))
  common_genes <- intersect(gene_map$gene, rownames(full_counts))
  input_counts <- full_counts[common_genes, ids, drop = FALSE]
  g <- gene_map[match(common_genes, gene)]
  cells <- data.table(
    cell_id = ids,
    reference_group = ifelse(ids %in% target_ids, "epithelial", "reference"),
    reference_lineage = NA_character_
  )
  cells[match(ref$cell_id, cell_id), reference_lineage := ref$reference_lineage]
  pt_dir <- file.path(run_root, pt)
  dir.create(pt_dir, recursive = TRUE, showWarnings = FALSE)
  mtx_path <- file.path(pt_dir, "counts.mtx")
  cells_path <- file.path(pt_dir, "cells.csv")
  genes_path <- file.path(pt_dir, "genes.csv")
  result_path <- file.path(pt_dir, "infercnvpy_by_cell.csv")
  Matrix::writeMM(input_counts, mtx_path)
  fwrite(cells, cells_path, na = "NA")
  fwrite(g, genes_path, na = "NA")
  err_path <- file.path(pt_dir, "infercnvpy_error.txt")
  status <- system2(
    python,
    c(shQuote(worker), shQuote(mtx_path), shQuote(cells_path),
      shQuote(genes_path), shQuote(result_path)),
    stdout = file.path(pt_dir, "infercnvpy_stdout.log"),
    stderr = err_path
  )
  if (status == 0 && file.exists(result_path)) {
    res <- fread(result_path)[reference_group == "epithelial"]
    res[, `:=`(dataset_id = "GSE154600", patient_id = pt,
               infercnv_method = method)]
    per_patient[[pt]] <- res
    run_audit[[pt]] <- data.table(
      patient_id = pt, n_targets = length(target_ids),
      n_reference_T_NK = ref[reference_lineage == "T_NK", .N],
      n_reference_B = ref[reference_lineage == "B", .N],
      n_reference_Myeloid = ref[reference_lineage == "Myeloid", .N],
      run_status = "COMPLETED", error = NA_character_
    )
  } else {
    error <- if (file.exists(err_path))
      paste(readLines(err_path, warn = FALSE), collapse = " ") else
        paste0("infercnvpy exit status ", status)
    per_patient[[pt]] <- data.table(
      cell_id = target_ids, reference_group = "epithelial",
      infercnv_score = NA_real_, reference_score_median = NA_real_,
      reference_score_mad = NA_real_, infercnv_score_threshold = NA_real_,
      n_high_chromosome_arms = NA_integer_,
      infercnv_status = "NOT_EVALUABLE", dataset_id = "GSE154600",
      patient_id = pt, infercnv_method = method
    )
    run_audit[[pt]] <- data.table(
      patient_id = pt, n_targets = length(target_ids),
      n_reference_T_NK = ref[reference_lineage == "T_NK", .N],
      n_reference_B = ref[reference_lineage == "B", .N],
      n_reference_Myeloid = ref[reference_lineage == "Myeloid", .N],
      run_status = "FAILED", error = error
    )
  }
  rm(input_counts)
  gc()
}
infer <- rbindlist(per_patient, fill = TRUE)
audit <- rbindlist(run_audit, fill = TRUE)
fwrite(audit, file.path(out_dir, "infercnv_run_audit.csv"), na = "NA")

consensus <- merge(
  stability[, .(dataset_id, patient_id, cell_id, stability_class)],
  infer[, .(
    patient_id, cell_id, infercnv_score, reference_score_median,
    reference_score_mad, infercnv_score_threshold,
    n_high_chromosome_arms, infercnv_status, infercnv_method
  )],
  by = c("patient_id", "cell_id"), all.x = TRUE
)
consensus[is.na(infercnv_status), infercnv_status := "NOT_EVALUABLE"]
consensus[, integrated_cnv_evidence := fcase(
  infercnv_status == "NOT_EVALUABLE", "NOT_EVALUABLE",
  stability_class == "STABLE_ANEUPLOID" &
    infercnv_status == "HIGH_CNV_SUPPORT", "DUAL_METHOD_MALIGNANT_SUPPORT",
  stability_class == "STABLE_ANEUPLOID", "COPYKAT_ONLY_SUPPORT",
  stability_class != "STABLE_ANEUPLOID" &
    infercnv_status == "HIGH_CNV_SUPPORT", "INFERCNV_ONLY_SUPPORT",
  default = "NO_DUAL_SUPPORT"
)]
fwrite(
  consensus,
  file.path(out_dir, "GSE154600_copykat_infercnv_consensus_by_cell.csv.gz"),
  compress = "gzip", na = "NA"
)

validation <- consensus[, .(
  n_epithelial_with_counts = .N,
  n_infercnv_high = sum(infercnv_status == "HIGH_CNV_SUPPORT"),
  n_copykat_stable_aneuploid = sum(stability_class == "STABLE_ANEUPLOID"),
  n_dual_method_support =
    sum(integrated_cnv_evidence == "DUAL_METHOD_MALIGNANT_SUPPORT"),
  n_copykat_only = sum(integrated_cnv_evidence == "COPYKAT_ONLY_SUPPORT"),
  n_infercnv_only = sum(integrated_cnv_evidence == "INFERCNV_ONLY_SUPPORT"),
  n_neither = sum(integrated_cnv_evidence == "NO_DUAL_SUPPORT"),
  n_not_evaluable = sum(integrated_cnv_evidence == "NOT_EVALUABLE")
), by = .(dataset_id, patient_id)]
validation[, copykat_infercnv_concordance := fifelse(
  n_copykat_stable_aneuploid > 0,
  n_dual_method_support / n_copykat_stable_aneuploid, NA_real_
)]
validation[, independent_cnv_method := method]
fwrite(validation, out, na = "NA")

concordance <- melt(
  validation,
  id.vars = c("dataset_id", "patient_id", "n_epithelial_with_counts"),
  measure.vars = c(
    "n_dual_method_support", "n_copykat_only", "n_infercnv_only",
    "n_neither", "n_not_evaluable"
  ),
  variable.name = "cnv_evidence_class", value.name = "n_cells"
)
concordance[, fraction_of_epithelial := n_cells / n_epithelial_with_counts]
fwrite(concordance,
       file.path(out_dir, "GSE154600_copykat_infercnv_concordance.csv"),
       na = "NA")

p <- ggplot(concordance, aes(patient_id, n_cells, fill = cnv_evidence_class)) +
  geom_col(position = "fill") +
  theme_bw() +
  labs(
    title = "GSE154600 CopyKAT and independent inferCNV concordance",
    subtitle = "infercnvpy fallback; continuous CNV signal, no HMM",
    x = NULL, y = "fraction of epithelial cells", fill = "evidence"
  )
ggsave(file.path(out_dir, "copykat_infercnv_concordance.png"),
       p, width = 10, height = 5.5, dpi = 180)
message("Independent patient-level inferCNV validation complete")
