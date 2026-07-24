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
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
previous <- file.path(data_root, "research_validation_independent_cnv")
out_dir <- file.path(data_root, "research_spatial_transition")
run_root <- file.path(out_dir, "infercnvpy_heldout_runs")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(out_dir, "infercnvpy_heldout_control_run_audit.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

targets <- readRDS(file.path(
  previous, "GSE154600_complete_final_epithelial_counts.rds"
))
stability <- fread(file.path(
  previous, "GSE154600_copykat_stability_by_cell_v2.csv.gz"
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

# Reuse the frozen gene-position table produced in the preceding independent
# CNV task instead of remapping genes with a new annotation dependency.
gene_map <- fread(file.path(
  previous, "infercnv_patient_runs", "gene_order_from_org.Hs.eg.db.csv"
))
gene_map <- gene_map[
  chromosome %chin% paste0("chr", 1:22) &
    gene %in% rownames(full_counts) & !duplicated(gene)
]
fwrite(gene_map, file.path(run_root, "gene_order_from_org.Hs.eg.db.csv"),
       na = "NA")

python <- "C:/Users/chenfy12/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
python_lib <- file.path(previous, "python_lib")
if (!file.exists(python) || !dir.exists(file.path(python_lib, "infercnvpy")))
  stop("Task-authorized infercnvpy environment is unavailable")
Sys.setenv(PYTHONPATH = normalizePath(python_lib, winslash = "/", mustWork = TRUE))
worker <- file.path(run_root, "_infercnvpy_heldout_worker.py")
writeLines(c(
  "import sys",
  "import importlib",
  "import numpy as np",
  "import pandas as pd",
  "import scipy.io",
  "import anndata as ad",
  "import scanpy as sc",
  "import infercnvpy as cnv",
  "m = importlib.import_module('infercnvpy.tl._infercnv')",
  "m.process_map = lambda fn, *it, **kw: [fn(*args) for args in zip(*it)]",
  "mtx_path, cells_path, genes_path, out_path = sys.argv[1:5]",
  "cells = pd.read_csv(cells_path)",
  "genes = pd.read_csv(genes_path)",
  "x = scipy.io.mmread(mtx_path).tocsr().transpose().tocsr()",
  "adata = ad.AnnData(X=x)",
  "adata.obs_names = cells['cell_id'].astype(str).values",
  "adata.obs['analysis_group'] = cells['analysis_group'].astype(str).values",
  "adata.var_names = genes['gene'].astype(str).values",
  "adata.var['chromosome'] = genes['chromosome'].astype(str).values",
  "adata.var['start'] = genes['start'].astype(float).values",
  "adata.var['end'] = genes['end'].astype(float).values",
  "sc.pp.normalize_total(adata, target_sum=1e4)",
  "sc.pp.log1p(adata)",
  "cnv.tl.infercnv(adata, reference_key='analysis_group', reference_cat=['baseline'], lfc_clip=3, window_size=100, step=10, dynamic_threshold=1.5, exclude_chromosomes=('chrX','chrY'), chunksize=1000, n_jobs=1, inplace=True, calculate_gene_values=False)",
  "xc = adata.obsm['X_cnv'].tocsr()",
  "scores = np.asarray(abs(xc).mean(axis=1)).ravel()",
  "grp = adata.obs['analysis_group'].values",
  "baseline = grp == 'baseline'",
  "calibration = grp == 'calibration'",
  "if calibration.sum() == 0: raise RuntimeError('No held-out calibration controls')",
  "base_median = float(np.median(scores[baseline]))",
  "base_mad = float(np.median(np.abs(scores[baseline] - base_median)))",
  "base_threshold = base_median + 3.0 * base_mad",
  "cal_threshold = float(np.quantile(scores[calibration], 0.99))",
  "global_threshold = max(base_threshold, cal_threshold)",
  "centromere = {'chr1':123400000,'chr2':93900000,'chr3':90900000,'chr4':50000000,'chr5':48750000,'chr6':60550000,'chr7':60100000,'chr8':45200000,'chr9':43850000,'chr10':39800000,'chr11':53400000,'chr12':35500000,'chr13':17700000,'chr14':17150000,'chr15':19000000,'chr16':36850000,'chr17':25050000,'chr18':18450000,'chr19':26150000,'chr20':28050000,'chr21':11950000,'chr22':15550000}",
  "ordered = sorted(adata.uns['cnv']['chr_pos'].items(), key=lambda kv: kv[1])",
  "arm_labels = []",
  "for ii, (chrom, first) in enumerate(ordered):",
  "    last = ordered[ii+1][1] if ii+1 < len(ordered) else xc.shape[1]",
  "    nwin = last - first",
  "    gv = adata.var.loc[adata.var['chromosome'] == chrom].sort_values('start')",
  "    if len(gv) > 100:",
  "        centers = [float(gv.iloc[s:s+100]['start'].median()) for s in range(0, len(gv)-100+1, 10)]",
  "    else:",
  "        centers = [float(gv['start'].median())]",
  "    if len(centers) != nwin:",
  "        centers = np.linspace(float(gv['start'].min()), float(gv['start'].max()), nwin).tolist()",
  "    arm_labels.extend([chrom + ('p' if p < centromere[chrom] else 'q') for p in centers])",
  "arm_labels = np.asarray(arm_labels)",
  "n_high_arms = np.zeros(adata.n_obs, dtype=int)",
  "for arm in np.unique(arm_labels):",
  "    cols = np.where(arm_labels == arm)[0]",
  "    arm_score = np.asarray(abs(xc[:, cols]).mean(axis=1)).ravel()",
  "    arm_threshold = float(np.quantile(arm_score[calibration], 0.99))",
  "    n_high_arms += arm_score > arm_threshold",
  "high = (scores > global_threshold) & (n_high_arms >= 2)",
  "result = pd.DataFrame({'cell_id':adata.obs_names,'analysis_group':grp,'infercnv_score':scores,'baseline_threshold':base_threshold,'calibration_threshold':cal_threshold,'global_threshold':global_threshold,'n_high_chromosome_arms':n_high_arms,'calibrated_high':high})",
  "result.to_csv(out_path,index=False)"
), worker)

split_references <- function(pt, seed) {
  x <- refs[patient_id == pt & cell_id %in% colnames(full_counts)]
  set.seed(seed)
  rbindlist(lapply(unique(x$reference_lineage), function(lineage) {
    z <- copy(x[reference_lineage == lineage])
    z <- z[sample(seq_len(nrow(z)))]
    if (nrow(z) < 10L) {
      z[, analysis_group := "baseline"]
    } else {
      nb <- floor(.60 * nrow(z))
      nc <- floor(.20 * nrow(z))
      z[, analysis_group := c(
        rep("baseline", nb), rep("calibration", nc),
        rep("test", nrow(z) - nb - nc)
      )]
    }
    z
  }))
}

patients <- c("T59", "T76", "T77", "T89", "T90")
seeds <- c(20260725L, 20260726L, 20260727L)
audits <- list()
predictions <- list()
k <- 0L
for (pt in patients) {
  target_ids <- stability[patient_id == pt, cell_id]
  for (seed in seeds) {
    k <- k + 1L
    split <- split_references(pt, seed)
    ids <- c(target_ids, split$cell_id)
    common_genes <- intersect(gene_map$gene, rownames(full_counts))
    input_counts <- full_counts[common_genes, ids, drop = FALSE]
    g <- gene_map[match(common_genes, gene)]
    cells <- rbind(
      data.table(
        cell_id = target_ids, analysis_group = "epithelial",
        reference_lineage = "Epithelial"
      ),
      split[, .(cell_id, analysis_group, reference_lineage)]
    )
    safe <- paste0(pt, "__seed", seed)
    rd <- file.path(run_root, safe)
    dir.create(rd, recursive = TRUE, showWarnings = FALSE)
    mtx <- file.path(rd, "counts.mtx")
    cell_path <- file.path(rd, "cells.csv")
    gene_path <- file.path(rd, "genes.csv")
    result_path <- file.path(rd, "heldout_by_cell.csv")
    Matrix::writeMM(input_counts, mtx)
    fwrite(cells, cell_path)
    fwrite(g, gene_path)
    err_path <- file.path(rd, "infercnvpy_error.txt")
    status <- system2(
      python,
      c(shQuote(worker), shQuote(mtx), shQuote(cell_path),
        shQuote(gene_path), shQuote(result_path)),
      stdout = file.path(rd, "infercnvpy_stdout.log"),
      stderr = err_path
    )
    if (status == 0 && file.exists(result_path)) {
      res <- fread(result_path)
      res[, `:=`(
        dataset_id = "GSE154600", patient_id = pt, split_seed = seed
      )]
      predictions[[k]] <- res
      test <- res[analysis_group == "test"]
      audits[[k]] <- data.table(
        dataset_id = "GSE154600", patient_id = pt, split_seed = seed,
        n_epithelial = length(target_ids),
        n_baseline = split[analysis_group == "baseline", .N],
        n_calibration = split[analysis_group == "calibration", .N],
        n_test = split[analysis_group == "test", .N],
        n_baseline_T_NK = split[analysis_group == "baseline" & reference_lineage == "T_NK", .N],
        n_baseline_B = split[analysis_group == "baseline" & reference_lineage == "B", .N],
        n_baseline_Myeloid = split[analysis_group == "baseline" & reference_lineage == "Myeloid", .N],
        baseline_threshold = res$baseline_threshold[1],
        calibration_control_p99 = res$calibration_threshold[1],
        global_threshold = res$global_threshold[1],
        n_test_false_positive = sum(test$calibrated_high),
        test_control_false_positive_rate = mean(test$calibrated_high),
        run_status = "COMPLETED", error = NA_character_
      )
    } else {
      error <- if (file.exists(err_path))
        paste(readLines(err_path, warn = FALSE), collapse = " ") else
          paste0("infercnvpy exit status ", status)
      audits[[k]] <- data.table(
        dataset_id = "GSE154600", patient_id = pt, split_seed = seed,
        n_epithelial = length(target_ids),
        n_baseline = split[analysis_group == "baseline", .N],
        n_calibration = split[analysis_group == "calibration", .N],
        n_test = split[analysis_group == "test", .N],
        n_baseline_T_NK = NA_integer_, n_baseline_B = NA_integer_,
        n_baseline_Myeloid = NA_integer_, baseline_threshold = NA_real_,
        calibration_control_p99 = NA_real_, global_threshold = NA_real_,
        n_test_false_positive = NA_integer_,
        test_control_false_positive_rate = NA_real_,
        run_status = "FAILED", error = error
      )
    }
    rm(input_counts)
    gc()
  }
}
audit <- rbindlist(audits, fill = TRUE)
pred <- rbindlist(predictions, fill = TRUE)
fwrite(audit, out, na = "NA")
fwrite(pred, file.path(out_dir, "infercnvpy_heldout_predictions.csv.gz"),
       compress = "gzip", na = "NA")

fpr <- audit[, .(
  n_splits_completed = sum(run_status == "COMPLETED"),
  median_test_control_fpr =
    median(test_control_false_positive_rate, na.rm = TRUE),
  n_splits_fpr_le_0_05 =
    sum(test_control_false_positive_rate <= .05, na.rm = TRUE),
  fpr_seed_20260725 =
    test_control_false_positive_rate[match(20260725L, split_seed)],
  fpr_seed_20260726 =
    test_control_false_positive_rate[match(20260726L, split_seed)],
  fpr_seed_20260727 =
    test_control_false_positive_rate[match(20260727L, split_seed)]
), by = .(dataset_id, patient_id)]
fpr[, fpr_calibration_pass :=
      n_splits_completed == 3L &
      median_test_control_fpr <= .05 &
      n_splits_fpr_le_0_05 >= 2L]
fpr[, infercnv_threshold_status := fifelse(
  fpr_calibration_pass, "INFERCNV_THRESHOLD_CALIBRATED",
  "INFERCNV_THRESHOLD_UNSTABLE"
)]
fwrite(fpr, file.path(out_dir, "infercnvpy_negative_control_fpr_by_patient.csv"),
       na = "NA")

epi <- pred[analysis_group == "epithelial", .(
  n_splits_evaluated = .N,
  n_calibrated_high = sum(calibrated_high),
  n_calibrated_low = sum(!calibrated_high),
  median_infercnv_score = median(infercnv_score),
  median_global_threshold = median(global_threshold),
  median_high_chromosome_arms = median(n_high_chromosome_arms)
), by = .(dataset_id, patient_id, cell_id)]
epi <- merge(epi, fpr[, .(
  patient_id, fpr_calibration_pass, median_test_control_fpr,
  infercnv_threshold_status
)], by = "patient_id", all.x = TRUE)
epi[, calibrated_infercnv_status := fcase(
  !fpr_calibration_pass | n_splits_evaluated < 2L, "NOT_EVALUABLE",
  n_calibrated_high >= 2L, "CALIBRATED_HIGH_CNV",
  n_calibrated_low >= 2L, "CALIBRATED_LOW_CNV",
  default = "UNSTABLE_CNV"
)]
cell <- merge(
  stability[, .(dataset_id, patient_id, cell_id, stability_class)],
  epi[, .(
    patient_id, cell_id, n_splits_evaluated, n_calibrated_high,
    n_calibrated_low, median_infercnv_score, median_global_threshold,
    median_high_chromosome_arms, median_test_control_fpr,
    infercnv_threshold_status, calibrated_infercnv_status
  )],
  by = c("patient_id", "cell_id"), all.x = TRUE
)
cell[is.na(calibrated_infercnv_status),
     calibrated_infercnv_status := "NOT_EVALUABLE"]
cell[, integrated_calibrated_cnv_evidence := fcase(
  calibrated_infercnv_status == "NOT_EVALUABLE", "NOT_EVALUABLE",
  stability_class == "STABLE_ANEUPLOID" &
    calibrated_infercnv_status == "CALIBRATED_HIGH_CNV",
  "CALIBRATED_DUAL_METHOD_SUPPORT",
  stability_class == "STABLE_ANEUPLOID", "COPYKAT_ONLY_SUPPORT",
  stability_class != "STABLE_ANEUPLOID" &
    calibrated_infercnv_status == "CALIBRATED_HIGH_CNV",
  "INFERCNV_ONLY_SUPPORT",
  default = "NO_DUAL_SUPPORT"
)]
fwrite(
  cell, file.path(out_dir, "GSE154600_calibrated_cnv_by_cell.csv.gz"),
  compress = "gzip", na = "NA"
)
by_patient <- cell[, .(
  n_final_epithelial = .N,
  n_copykat_stable_aneuploid = sum(stability_class == "STABLE_ANEUPLOID"),
  n_calibrated_high_cnv =
    sum(calibrated_infercnv_status == "CALIBRATED_HIGH_CNV"),
  n_calibrated_low_cnv =
    sum(calibrated_infercnv_status == "CALIBRATED_LOW_CNV"),
  n_unstable_cnv = sum(calibrated_infercnv_status == "UNSTABLE_CNV"),
  n_calibrated_dual_method_support =
    sum(integrated_calibrated_cnv_evidence == "CALIBRATED_DUAL_METHOD_SUPPORT"),
  n_copykat_only =
    sum(integrated_calibrated_cnv_evidence == "COPYKAT_ONLY_SUPPORT"),
  n_infercnv_only =
    sum(integrated_calibrated_cnv_evidence == "INFERCNV_ONLY_SUPPORT"),
  n_neither = sum(integrated_calibrated_cnv_evidence == "NO_DUAL_SUPPORT"),
  n_not_evaluable =
    sum(integrated_calibrated_cnv_evidence == "NOT_EVALUABLE")
), by = .(dataset_id, patient_id)]
by_patient <- merge(by_patient, fpr, by = c("dataset_id", "patient_id"))
by_patient[, copykat_calibrated_infercnv_concordance := fifelse(
  n_copykat_stable_aneuploid > 0,
  n_calibrated_dual_method_support / n_copykat_stable_aneuploid,
  NA_real_
)]
fwrite(by_patient,
       file.path(out_dir, "GSE154600_calibrated_cnv_by_patient.csv"),
       na = "NA")

fpr_plot <- audit[, .(
  patient_id, metric = "held-out test-control FPR",
  value = test_control_false_positive_rate, split_seed = factor(split_seed)
)]
composition <- melt(
  by_patient,
  id.vars = "patient_id",
  measure.vars = c(
    "n_calibrated_dual_method_support", "n_copykat_only",
    "n_infercnv_only", "n_neither", "n_not_evaluable"
  ),
  variable.name = "split_seed", value.name = "n"
)
composition[, `:=`(
  metric = "calibrated CNV evidence fraction",
  value = n / by_patient$n_final_epithelial[match(patient_id, by_patient$patient_id)]
)]
plot_dt <- rbind(
  fpr_plot,
  composition[, .(patient_id, metric, value,
                  split_seed = factor(split_seed))],
  fill = TRUE
)
p <- ggplot(plot_dt, aes(patient_id, value, fill = split_seed)) +
  geom_col(position = "dodge") +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  geom_hline(
    data = data.frame(
      metric = "held-out test-control FPR", threshold = .05
    ),
    aes(yintercept = threshold), inherit.aes = FALSE,
    linetype = 2, color = "red"
  ) +
  theme_bw() +
  labs(
    title = "Held-out inferCNV calibration and CopyKAT concordance",
    x = NULL, y = "fraction", fill = "split/evidence"
  )
ggsave(file.path(out_dir, "infercnvpy_heldout_fpr_and_concordance.png"),
       p, width = 11, height = 7, dpi = 180)
message("Held-out infercnvpy negative-control calibration complete")
