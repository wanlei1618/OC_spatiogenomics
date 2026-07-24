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
repo_root <- normalizePath(cfg$project$repo_root, winslash = "/", mustWork = TRUE)
v4 <- file.path(data_root, "diagnostics_v4_cross_dataset_validation")
v5 <- file.path(data_root, "diagnostics_v5_final_calibration")
v6 <- file.path(data_root, "diagnostics_v6_malignant_receiver_validation")
cleaned <- file.path(data_root, "diagnostics_v2_marker_ready_cleaned")
run_root <- file.path(v6, "copykat_runs")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v6, "malignancy_consensus_by_cell.csv.gz")
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

collect_columns <- function(mats, cells) {
  parts <- lapply(mats, function(m) {
    hit <- intersect(cells, colnames(m))
    if (length(hit)) m[, hit, drop = FALSE] else NULL
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  genes <- Reduce(intersect, lapply(parts, rownames))
  ans <- do.call(cbind, lapply(parts, function(m) m[genes, , drop = FALSE]))
  ans[, !duplicated(colnames(ans)), drop = FALSE]
}

# Audit target-specific prior results. The v5 chromosome-window ratio is
# retained for provenance only and is explicitly not reusable as malignancy.
candidate_paths <- unique(c(
  list.files(data_root, pattern = "(infercnv|copykat|malignan|aneuploid|diploid|cnv)",
             recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
  list.files(file.path(repo_root, "infercnv"),
             pattern = "(infercnv|copykat|malignan|aneuploid|diploid|cnv)",
             recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
))
candidate_paths <- candidate_paths[
  grepl("GSE154600|GSE147082|GSE158722", candidate_paths, ignore.case = TRUE)
]
candidate_paths <- candidate_paths[
  !startsWith(normalizePath(candidate_paths, winslash = "/", mustWork = FALSE),
              normalizePath(v6, winslash = "/", mustWork = FALSE)) &
    !grepl("/workflow/scripts/|\\.(R|png|pdf)$", candidate_paths, ignore.case = TRUE)
]
audit <- if (length(candidate_paths)) rbindlist(lapply(candidate_paths, function(p) {
  ds <- regmatches(p, regexpr("GSE(154600|147082|158722)", p, ignore.case = TRUE))
  base <- basename(p)
  method <- if (grepl("copykat", base, ignore.case = TRUE)) "CopyKAT" else
    if (grepl("infercnv", base, ignore.case = TRUE)) "inferCNV" else "CNV_like_or_label"
  formal <- method %in% c("CopyKAT", "inferCNV") &&
    !grepl("cnv_like|sensitivity", p, ignore.case = TRUE)
  data.table(
    dataset_id = toupper(ds),
    result_path = normalizePath(p, winslash = "/", mustWork = FALSE),
    method = method,
    reference_cells = NA_character_,
    n_epithelial = NA_integer_,
    n_malignant = NA_integer_,
    n_diploid = NA_integer_,
    quality_status = if (formal) "REQUIRES_TARGET_ID_AND_REFERENCE_AUDIT" else
      "EXPLORATORY_NOT_FORMAL",
    reuse_decision = "DO_NOT_REUSE"
  )
}), fill = TRUE) else data.table()
for (ds in c("GSE154600", "GSE147082", "GSE158722")) {
  if (!nrow(audit[dataset_id == ds])) audit <- rbind(
    audit,
    data.table(
      dataset_id = ds, result_path = NA_character_, method = "none_found",
      reference_cells = NA_character_, n_epithelial = NA_integer_,
      n_malignant = NA_integer_, n_diploid = NA_integer_,
      quality_status = "NO_REUSABLE_TARGET_SPECIFIC_RESULT",
      reuse_decision = "RUN_OR_MARK_NOT_EVALUABLE"
    ),
    fill = TRUE
  )
}

# CopyKAT is installed for R 4.0 on this workstation, whereas the Seurat v5
# objects require current R/Seurat. This local worker consumes version-2 RDS
# sparse matrices prepared here and writes only the prediction table.
worker <- file.path(run_root, "_copykat_worker.R")
writeLines(c(
  "options(stringsAsFactors=FALSE,warn=1)",
  ".libPaths(c('D:/Documents/R/win-library/4.0',.libPaths()))",
  "suppressPackageStartupMessages({library(copykat);library(data.table);library(Matrix)})",
  "a<-commandArgs(trailingOnly=TRUE)",
  "x<-readRDS(a[[1]]);dir.create(a[[3]],recursive=TRUE,showWarnings=FALSE);setwd(a[[3]])",
  "keep<-Matrix::rowMeans(x$counts>0)>.05;raw<-as.matrix(x$counts[keep,,drop=FALSE])",
  "ans<-tryCatch(copykat(rawmat=raw,id.type='S',ngene.chr=5,min.gene.per.cell=200,LOW.DR=.05,UP.DR=.10,win.size=25,norm.cell.names=x$normal_cells,KS.cut=.10,sam.name=a[[4]],distance='euclidean',test.emd='FALSE',output.seg='FALSE',plot.genes=FALSE,genome='hg20',n.cores=1),error=function(e)e)",
  "if(inherits(ans,'error')){writeLines(conditionMessage(ans),a[[5]]);quit(status=2)}",
  "fwrite(as.data.table(ans$prediction),a[[2]],na='NA')"
), worker)

r40 <- "D:/R/R-4.0.3/bin/x64/Rscript.exe"
if (!file.exists(r40)) stop("R 4.0 CopyKAT runtime not found: ", r40)
run_copykat <- function(dataset_id, patient_id, counts, normal_cells) {
  safe <- gsub("[^A-Za-z0-9_.-]", "_", paste(dataset_id, patient_id, sep = "__"))
  rd <- file.path(run_root, paste0(safe, "_input.rds"))
  pred <- file.path(run_root, paste0(safe, "_prediction.csv"))
  err <- file.path(run_root, paste0(safe, "_error.txt"))
  if (!file.exists(pred)) {
    saveRDS(list(counts = as(counts, "dgCMatrix"), normal_cells = normal_cells),
            rd, version = 2)
    run_dir <- file.path(run_root, safe)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    status <- system2(r40, c(
      shQuote(worker), shQuote(rd), shQuote(pred), shQuote(run_dir),
      shQuote(safe), shQuote(err)
    ))
    if (status != 0 || !file.exists(pred)) {
      msg <- if (file.exists(err)) paste(readLines(err, warn = FALSE), collapse = " ") else
        paste0("CopyKAT worker exit status ", status)
      return(list(prediction = NULL, error = msg, prediction_path = pred))
    }
  }
  list(prediction = fread(pred), error = NA_character_, prediction_path = pred)
}

predictions <- list()
run_status <- list()

# GSE154600: broad epithelial lineage is patient-enriched but has canonical
# epithelial support and no incompatible lineage program. Run each patient
# separately with same-patient immune reference cells.
a154 <- fread(file.path(cleaned, "GSE154600", "cleaned_cell_assignments.csv.gz"))
c154 <- fread(file.path(cleaned, "GSE154600", "cleaned_cluster_annotation_template.csv"))
a154 <- map_status(a154, c154)
inputs154 <- lapply(
  c("Epithelial_like", "T_NK_like", "B_Plasma_like", "Myeloid_like"),
  function(x) readRDS(file.path(
    data_root, "diagnostics_v2", "objects", "GSE154600", "lineage_inputs",
    paste0(x, "_strategy_input.rds")
  ))$counts
)
names(inputs154) <- c("epi", "tnk", "b", "my")
epi154 <- a154[
  final_cell_type == "Epithelial" &
    annotation_status == "REVIEW_PATIENT_ENRICHED" &
    canonical_support_n >= 3 & incompatible_lineage_program != TRUE,
  cell_id
]
immune_types <- c("T_cell", "NK_cell", "B_cell", "Macrophage",
                  "Monocyte", "cDC1", "cDC2", "pDC")
ref154 <- a154[
  final_cell_type %in% immune_types &
    annotation_status == "READY_HIGH_CONFIDENCE" &
    canonical_support_n >= 3 & incompatible_lineage_program != TRUE,
  cell_id
]
set.seed(as.integer(cfg$project$random_seed %||% 20260718L))
for (pt in unique(a154$patient_id)) {
  targets <- intersect(a154[patient_id == pt & cell_id %in% epi154, cell_id],
                       colnames(inputs154$epi))
  refs <- a154[patient_id == pt & cell_id %in% ref154, cell_id]
  refs <- intersect(refs, unique(unlist(lapply(inputs154[-1], colnames), use.names = FALSE)))
  if (length(refs) > 500) refs <- sample(refs, 500)
  if (length(targets) < 20 || length(refs) < 50) {
    run_status[[paste("GSE154600", pt)]] <- data.table(
      dataset_id = "GSE154600", patient_id = pt, method = "CopyKAT",
      n_target = length(targets), n_reference = length(refs),
      status = "NOT_EVALUABLE", error = "fewer than 20 targets or 50 same-patient immune references"
    )
    next
  }
  mat <- collect_columns(inputs154, c(targets, refs))
  res <- run_copykat("GSE154600", pt, mat, refs)
  run_status[[paste("GSE154600", pt)]] <- data.table(
    dataset_id = "GSE154600", patient_id = pt, method = "CopyKAT",
    n_target = length(targets), n_reference = length(refs),
    status = if (is.null(res$prediction)) "FAILED" else "COMPLETED",
    error = res$error
  )
  if (!is.null(res$prediction)) {
    p <- res$prediction
    setnames(p, names(p)[1:2], c("cell_id", "copykat_status"))
    predictions[[paste("GSE154600", pt)]] <- p[cell_id %in% targets][, `:=`(
      dataset_id = "GSE154600", patient_id = pt
    )]
    audit <- rbind(audit, data.table(
      dataset_id = "GSE154600", result_path = normalizePath(
        res$prediction_path, winslash = "/", mustWork = FALSE),
      method = "CopyKAT", reference_cells = "same-patient high-confidence immune",
      n_epithelial = length(targets),
      n_malignant = sum(tolower(p[cell_id %in% targets]$copykat_status) == "aneuploid"),
      n_diploid = sum(tolower(p[cell_id %in% targets]$copykat_status) == "diploid"),
      quality_status = "FORMAL_SINGLE_METHOD",
      reuse_decision = "REUSE_CURRENT_RUN"
    ), fill = TRUE)
  }
  rm(mat); gc()
}
rm(inputs154, c154)

# GSE147082: formal PT-2834 assessment of cluster 4, cluster 7, and the 15
# broad epithelial cells using same-patient immune references.
a147 <- fread(file.path(v6, "GSE147082_refined_cell_assignments_v3.csv.gz"))
obj147 <- readRDS(file.path(data_root, "GSE147082", "objects", "GSE147082_preannotation.rds"))
cnt147 <- SeuratObject::LayerData(obj147, assay = "RNA", layer = "counts")
pt <- "PT-2834"
targets147 <- a147[
  patient_id == pt & (final_cluster %in% c(4, 7) | final_cell_type == "Epithelial"),
  cell_id
]
refs147 <- a147[
  patient_id == pt & final_cell_type %in%
    c("T_cell", "Cytotoxic_T", "CD8_effector_T", "Gamma_delta_T",
      "B_cell", "Macrophage", "pDC", "Mast"),
  cell_id
]
targets147 <- intersect(targets147, colnames(cnt147))
refs147 <- intersect(refs147, colnames(cnt147))
if (length(refs147) > 500) refs147 <- sample(refs147, 500)
if (length(targets147) >= 20 && length(refs147) >= 10) {
  mat <- cnt147[, c(targets147, refs147), drop = FALSE]
  res <- run_copykat("GSE147082", pt, mat, refs147)
  run_status[["GSE147082 PT-2834"]] <- data.table(
    dataset_id = "GSE147082", patient_id = pt, method = "CopyKAT",
    n_target = length(targets147), n_reference = length(refs147),
    status = if (is.null(res$prediction)) "FAILED" else "COMPLETED",
    error = res$error
  )
  if (!is.null(res$prediction)) {
    p <- res$prediction
    setnames(p, names(p)[1:2], c("cell_id", "copykat_status"))
    predictions[["GSE147082 PT-2834"]] <- p[cell_id %in% targets147][, `:=`(
      dataset_id = "GSE147082", patient_id = pt
    )]
    audit <- rbind(audit, data.table(
      dataset_id = "GSE147082", result_path = normalizePath(
        res$prediction_path, winslash = "/", mustWork = FALSE),
      method = "CopyKAT", reference_cells = "same-patient high-confidence immune",
      n_epithelial = length(targets147),
      n_malignant = sum(tolower(p[cell_id %in% targets147]$copykat_status) == "aneuploid"),
      n_diploid = sum(tolower(p[cell_id %in% targets147]$copykat_status) == "diploid"),
      quality_status = "FORMAL_SINGLE_METHOD",
      reuse_decision = "REUSE_CURRENT_RUN"
    ), fill = TRUE)
  }
} else {
  run_status[["GSE147082 PT-2834"]] <- data.table(
    dataset_id = "GSE147082", patient_id = pt, method = "CopyKAT",
    n_target = length(targets147), n_reference = length(refs147),
    status = "NOT_EVALUABLE",
    error = "fewer than 20 targets or 10 same-patient immune references"
  )
}
rm(obj147, cnt147)

# GSE158722 is not formally run: platform cannot be reliably recovered from
# metadata and same-platform reference selection would therefore be fabricated.
a158 <- fread(file.path(cleaned, "GSE158722", "cleaned_cell_assignments.csv.gz"))
run_status[["GSE158722"]] <- data.table(
  dataset_id = "GSE158722", patient_id = "ALL", method = "CopyKAT",
  n_target = sum(a158$final_cell_type == "Epithelial"),
  n_reference = sum(a158$final_cell_type %in% immune_types),
  status = "NOT_EVALUABLE",
  error = "platform identity unavailable; reliable same-platform reference selection not possible"
)

pred <- rbindlist(predictions, fill = TRUE)
if (!nrow(pred)) pred <- data.table(
  cell_id = character(), copykat_status = character(),
  dataset_id = character(), patient_id = character()
)
pred[, copykat_status := tolower(copykat_status)]

targets <- rbindlist(list(
  a154[cell_id %in% epi154, .(
    dataset_id, cell_id, patient_id, sample_id, final_cluster,
    pre_cnv_cell_type = final_cell_type
  )],
  a147[patient_id == pt & (final_cluster %in% c(4, 7) | final_cell_type == "Epithelial"),
       .(dataset_id, cell_id, patient_id, sample_id, final_cluster,
         pre_cnv_cell_type = final_cell_type)],
  a158[final_cell_type == "Epithelial", .(
    dataset_id, cell_id, patient_id, sample_id, final_cluster,
    pre_cnv_cell_type = final_cell_type
  )]
), fill = TRUE)
targets <- merge(
  targets, pred[, .(dataset_id, cell_id, copykat_status)],
  by = c("dataset_id", "cell_id"), all.x = TRUE
)
targets[, `:=`(
  infercnv_status = "NOT_AVAILABLE",
  infercnv_score = NA_real_,
  copykat_aneuploid_probability = NA_real_
)]
targets[is.na(copykat_status), copykat_status := "not_evaluated"]
targets[, malignancy_consensus := fcase(
  copykat_status == "aneuploid", "MALIGNANT_SUPPORTIVE",
  default = "NOT_EVALUABLE"
)]
targets[, formal_method_note := fcase(
  dataset_id == "GSE158722",
  "CopyKAT not run because platform identity and same-platform reference were unavailable",
  copykat_status == "aneuploid",
  "CopyKAT aneuploid; inferCNV unavailable, so no two-method high-confidence consensus",
  copykat_status == "diploid",
  "CopyKAT diploid; inferCNV unavailable, so DIPLOID_SUPPORTIVE consensus not assigned",
  default = "formal method unavailable, failed, or did not return a defined call"
)]
fwrite(targets, out, compress = "gzip", na = "NA")

summary <- targets[, .(
  n_candidate_cells = .N,
  n_copykat_evaluated = sum(copykat_status %in% c("aneuploid", "diploid")),
  n_copykat_aneuploid = sum(copykat_status == "aneuploid"),
  n_copykat_diploid = sum(copykat_status == "diploid"),
  n_malignant_high_confidence = sum(malignancy_consensus == "MALIGNANT_HIGH_CONFIDENCE"),
  n_malignant_supportive = sum(malignancy_consensus == "MALIGNANT_SUPPORTIVE"),
  n_diploid_supportive = sum(malignancy_consensus == "DIPLOID_SUPPORTIVE"),
  malignancy_method = if (any(copykat_status %in% c("aneuploid", "diploid")))
    "CopyKAT_single_method" else "NOT_EVALUABLE"
), by = .(dataset_id, patient_id)]
fwrite(summary, file.path(v6, "malignancy_summary_by_patient.csv"), na = "NA")

cluster_cnv <- targets[
  dataset_id == "GSE147082" & patient_id == "PT-2834",
  .(
    n_cells = .N,
    n_copykat_evaluated = sum(copykat_status %in% c("aneuploid", "diploid")),
    n_aneuploid = sum(copykat_status == "aneuploid"),
    n_diploid = sum(copykat_status == "diploid"),
    aneuploid_fraction = if (sum(copykat_status %in% c("aneuploid", "diploid")) > 0)
      mean(copykat_status[copykat_status %in% c("aneuploid", "diploid")] == "aneuploid")
    else NA_real_
  ),
  by = .(
    target_group = fcase(
      final_cluster == 4, "cluster_4",
      final_cluster == 7, "cluster_7",
      default = "broad_epithelial"
    )
  )
]
cluster_cnv[, final_label := fcase(
  target_group == "cluster_4" & n_copykat_evaluated >= 20 & aneuploid_fraction >= .50,
  "Malignant_mesenchymal_candidate",
  target_group == "cluster_4", "Mesenchymal_stromal_candidate",
  target_group == "cluster_7" & n_copykat_evaluated >= 20 & aneuploid_fraction >= .50,
  "Malignant_cartilage_like_candidate",
  target_group == "cluster_7",
  "COL2A1_positive_chondrocyte_like_fibroblast_candidate",
  default = "Epithelial_single_cell_candidates"
)]
fwrite(cluster_cnv, file.path(v6, "GSE147082_PT2834_formal_cnv_cluster_summary.csv"), na = "NA")

# Only formal CopyKAT evidence may alter cluster 4/7 in the local full assignment.
for (cl in c(4L, 7L)) {
  label <- cluster_cnv[target_group == paste0("cluster_", cl), final_label][1]
  if (length(label) && !is.na(label)) {
    a147[patient_id == "PT-2834" & final_cluster == cl,
         `:=`(final_cell_type = label, cell_subtype = label)]
  }
}
fwrite(a147, file.path(v6, "GSE147082_refined_cell_assignments_v3.csv.gz"),
       compress = "gzip", na = "NA")
fwrite(rbindlist(run_status, fill = TRUE), file.path(v6, "malignancy_method_run_status.csv"), na = "NA")
fwrite(unique(audit), file.path(v6, "existing_malignancy_results_audit.csv"), na = "NA")

p <- melt(
  summary,
  id.vars = c("dataset_id", "patient_id"),
  measure.vars = c("n_copykat_aneuploid", "n_copykat_diploid"),
  variable.name = "formal_call", value.name = "n_cells"
)
p <- ggplot(p, aes(patient_id, n_cells, fill = formal_call)) +
  geom_col() +
  facet_wrap(~dataset_id, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Formal CopyKAT calls (single-method evidence)",
    subtitle = "No inferCNV result was available; aneuploid calls are supportive, not two-method high confidence",
    x = NULL, y = "candidate cells"
  )
ggsave(file.path(v6, "04_malignant_cnv_consensus.png"), p, width = 11, height = 5.5, dpi = 180)
message("Malignancy audit and formal CopyKAT consensus complete")
