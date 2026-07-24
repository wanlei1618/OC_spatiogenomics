options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
data_root <- normalizePath(z$cfg$project$data_root, winslash = "/", mustWork = TRUE)
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
out_dir <- file.path(data_root, "research_validation_independent_cnv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
trace_path <- file.path(out_dir, "missing_epithelial_count_trace.csv")
if (file.exists(trace_path) && !replace_generated) stop("Output exists: ", trace_path)

extract_assay5_counts <- function(object_path) {
  x <- readRDS(object_path)
  assay <- x@assays[[x@active.assay]]
  layers <- attr(assay, "layers")
  if (is.null(layers) || is.null(layers$counts))
    stop("Authoritative preannotation object does not contain an Assay5 counts layer")
  counts <- layers$counts
  feature_ids <- rownames(attr(assay, "features"))
  cell_ids <- rownames(attr(assay, "cells"))
  if (length(feature_ids) < nrow(counts) || length(cell_ids) < ncol(counts))
    stop("Assay5 feature/cell maps do not cover the counts layer")
  dimnames(counts) <- list(feature_ids[seq_len(nrow(counts))],
                           cell_ids[seq_len(ncol(counts))])
  list(counts = as(counts, "dgCMatrix"), metadata = x@meta.data)
}

object_path <- file.path(data_root, "GSE154600", "objects", "GSE154600_preannotation.rds")
if (!file.exists(object_path)) stop("Missing authoritative preannotation object: ", object_path)
full <- extract_assay5_counts(object_path)
full_counts <- full$counts
full_cells <- colnames(full_counts)

manifest <- fread(file.path(v61, "GSE154600_copykat_target_manifest.csv.gz"))
missing <- manifest[available_for_copykat == FALSE]
if (nrow(missing) != 511L)
  stop("Expected 511 missing final epithelial cells; found ", nrow(missing))

project_root <- normalizePath(file.path(data_root, "..", ".."),
                              winslash = "/", mustWork = TRUE)
raw_barcode_files <- Sys.glob(file.path(
  project_root, "*", "*", "raw_counts", "GSE154600", "extracted",
  "*_barcodes.tsv.gz"
))
read_barcode_file <- function(path) {
  patient <- sub(".*_(T[0-9]+)_barcodes.*", "\\1", basename(path))
  z <- readLines(gzfile(path), warn = FALSE)
  unique(c(z, paste0(patient, "__", z)))
}
raw_cells <- unique(unlist(lapply(raw_barcode_files, read_barcode_file),
                           use.names = FALSE))
filtered_files <- raw_barcode_files[
  grepl("filtered", raw_barcode_files, ignore.case = TRUE)
]
filtered_cells <- if (length(filtered_files))
  unique(unlist(lapply(filtered_files, read_barcode_file), use.names = FALSE)) else
    character()

trace <- missing[, .(
  patient_id,
  cell_id,
  found_in_full_object = cell_id %in% full_cells,
  found_in_raw_matrix = cell_id %in% raw_cells,
  found_in_filtered_matrix = cell_id %in% filtered_cells,
  found_in_any_lineage_input = available_for_copykat
)]
trace[, recovery_source := fcase(
  found_in_full_object & found_in_raw_matrix,
  "GSE154600_preannotation_Assay5_raw_counts;GEO_feature_barcode_matrix",
  found_in_full_object, "GSE154600_preannotation_Assay5_raw_counts",
  found_in_raw_matrix, "GEO_feature_barcode_matrix",
  default = NA_character_
)]
trace[, recoverable := found_in_full_object]
trace[, failure_reason := fcase(
  recoverable, NA_character_,
  !found_in_full_object & found_in_raw_matrix,
  "barcode_present_in_GEO_matrix_but_not_retained_in_preannotation_object",
  !found_in_full_object & !found_in_raw_matrix,
  "source_object_or_matrix_incomplete_after_exhaustive_authoritative_search",
  default = "other_confirmed_source_mismatch"
)]
fwrite(trace, trace_path, na = "NA")

target_ids <- manifest$cell_id
if (!all(target_ids %in% full_cells))
  stop("Some final epithelial cells remain absent from authoritative counts")
complete_counts <- full_counts[, target_ids, drop = FALSE]
saveRDS(
  complete_counts,
  file.path(out_dir, "GSE154600_complete_final_epithelial_counts.rds"),
  compress = TRUE, version = 2
)

coverage <- manifest[, .(n_final_epithelial = .N), by = .(dataset_id, patient_id)]
recovered <- trace[, .(n_counts_recovered = sum(recoverable)), by = patient_id]
coverage <- merge(coverage, recovered, by = "patient_id", all.x = TRUE)
coverage[is.na(n_counts_recovered), n_counts_recovered := 0L]
coverage[, n_counts_available_final := n_final_epithelial]
coverage[, final_count_coverage := n_counts_available_final / n_final_epithelial]
setcolorder(coverage, c(
  "dataset_id", "patient_id", "n_final_epithelial", "n_counts_recovered",
  "n_counts_available_final", "final_count_coverage"
))
fwrite(coverage, file.path(out_dir, "GSE154600_final_count_coverage_v2.csv"),
       na = "NA")

update_required <- coverage[, .(
  dataset_id,
  patient_id,
  n_counts_recovered,
  copykat_update_required = n_counts_recovered > 0,
  update_scope = fifelse(
    n_counts_recovered > 0,
    "rerun_existing_three_seeds_for_affected_patient",
    "no_rerun"
  ),
  parameters_changed = FALSE,
  seeds = "20260718;20260719;20260720"
)]
fwrite(update_required,
       file.path(out_dir, "copykat_update_required_by_patient.csv"), na = "NA")

# All five patients recovered cells, so rerun only those affected patients with
# the exact v6.1 seeds, CopyKAT parameters, and reference sampling logic.
reference_manifest <- fread(file.path(v61, "GSE154600_copykat_reference_manifest.csv.gz"))
reference_ids <- intersect(
  reference_manifest[available_for_copykat == TRUE, cell_id], full_cells
)
run_root <- file.path(out_dir, "copykat_recovery_runs")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
worker <- file.path(run_root, "_copykat_worker_recovery.R")
writeLines(c(
  "options(stringsAsFactors=FALSE,warn=1)",
  ".libPaths(c('D:/Documents/R/win-library/4.0',.libPaths()))",
  "suppressPackageStartupMessages({library(copykat);library(data.table);library(Matrix)})",
  "a<-commandArgs(trailingOnly=TRUE);seed<-as.integer(a[[6]]);set.seed(seed)",
  "x<-readRDS(a[[1]]);dir.create(a[[3]],recursive=TRUE,showWarnings=FALSE);setwd(a[[3]])",
  "keep<-Matrix::rowMeans(x$counts>0)>.05;raw<-as.matrix(x$counts[keep,,drop=FALSE])",
  "ans<-tryCatch(copykat(rawmat=raw,id.type='S',ngene.chr=5,min.gene.per.cell=200,LOW.DR=.05,UP.DR=.10,win.size=25,norm.cell.names=x$normal_cells,KS.cut=.10,sam.name=a[[4]],distance='euclidean',test.emd='FALSE',output.seg='FALSE',plot.genes=FALSE,genome='hg20',n.cores=1),error=function(e)e)",
  "if(inherits(ans,'error')){writeLines(conditionMessage(ans),a[[5]]);quit(status=2)}",
  "fwrite(as.data.table(ans$prediction),a[[2]],na='NA')"
), worker)
rscript <- "D:/R/R-4.0.3/bin/x64/Rscript.exe"
if (!file.exists(rscript)) stop("CopyKAT R runtime not found: ", rscript)

run_one <- function(pt, seed) {
  targets <- manifest[patient_id == pt, cell_id]
  refs_all <- reference_manifest[
    patient_id == pt &
      available_for_copykat == TRUE &
      cell_id %in% reference_ids,
    cell_id
  ]
  set.seed(seed)
  refs <- if (length(refs_all) > 500L) sample(refs_all, 500L) else refs_all
  ids <- unique(c(targets, refs))
  counts <- full_counts[, ids, drop = FALSE]
  safe <- paste0("GSE154600__", pt, "__seed", seed, "__recovered")
  input_path <- file.path(run_root, paste0(safe, "_input.rds"))
  pred_path <- file.path(run_root, paste0(safe, "_prediction.csv"))
  err_path <- file.path(run_root, paste0(safe, "_error.txt"))
  saveRDS(list(counts = counts, normal_cells = refs), input_path,
          compress = TRUE, version = 2)
  run_dir <- file.path(run_root, safe)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(rscript, c(
    shQuote(worker), shQuote(input_path), shQuote(pred_path),
    shQuote(run_dir), shQuote(safe), shQuote(err_path), as.character(seed)
  ))
  if (status != 0 || !file.exists(pred_path)) {
    error <- if (file.exists(err_path))
      paste(readLines(err_path, warn = FALSE), collapse = " ") else
        paste0("CopyKAT worker exit status ", status)
    return(list(
      run = data.table(
        dataset_id = "GSE154600", patient_id = pt, seed = seed,
        n_final_epithelial = length(targets), n_targets_submitted = length(targets),
        n_reference_used = length(refs), n_aneuploid = 0L, n_diploid = 0L,
        n_not_defined = length(targets), run_status = "FAILED", error = error
      ),
      prediction = NULL
    ))
  }
  p <- fread(pred_path)
  setnames(p, names(p)[1:2], c("cell_id", "copykat_status"))
  p[, copykat_status := tolower(copykat_status)]
  p <- p[cell_id %in% targets, .(cell_id, copykat_status)]
  p <- merge(data.table(cell_id = targets), p, by = "cell_id", all.x = TRUE)
  p[is.na(copykat_status), copykat_status := "not.defined"]
  p[, `:=`(dataset_id = "GSE154600", patient_id = pt, seed = seed)]
  list(
    run = data.table(
      dataset_id = "GSE154600", patient_id = pt, seed = seed,
      n_final_epithelial = length(targets), n_targets_submitted = length(targets),
      n_reference_used = length(refs),
      n_aneuploid = sum(p$copykat_status == "aneuploid"),
      n_diploid = sum(p$copykat_status == "diploid"),
      n_not_defined = sum(!p$copykat_status %in% c("aneuploid", "diploid")),
      run_status = "COMPLETED", error = NA_character_
    ),
    prediction = p
  )
}

patients <- c("T59", "T76", "T77", "T89", "T90")
seeds <- c(20260718L, 20260719L, 20260720L)
runs <- list()
predictions <- list()
k <- 0L
for (pt in patients) {
  for (seed in seeds) {
    k <- k + 1L
    res <- run_one(pt, seed)
    runs[[k]] <- res$run
    predictions[[k]] <- res$prediction
    gc()
  }
}
run_status <- rbindlist(runs, fill = TRUE)
pred <- rbindlist(Filter(Negate(is.null), predictions), fill = TRUE)
fwrite(run_status, file.path(out_dir, "copykat_recovery_run_status.csv"), na = "NA")
fwrite(pred, file.path(out_dir, "copykat_recovery_raw_predictions.csv.gz"),
       compress = "gzip", na = "NA")

completed <- run_status[run_status == "COMPLETED",
                        .(n_runs_attempted = .N), by = patient_id]
calls <- pred[, .(
  n_runs_defined = sum(copykat_status %in% c("aneuploid", "diploid")),
  n_aneuploid_calls = sum(copykat_status == "aneuploid"),
  n_diploid_calls = sum(copykat_status == "diploid")
), by = .(patient_id, cell_id)]
cell <- merge(
  manifest[, .(dataset_id, patient_id, cell_id)],
  completed, by = "patient_id", all.x = TRUE
)
cell <- merge(cell, calls, by = c("patient_id", "cell_id"), all.x = TRUE)
for (nm in c("n_runs_attempted", "n_runs_defined",
             "n_aneuploid_calls", "n_diploid_calls"))
  set(cell, which(is.na(cell[[nm]])), nm, 0L)
cell[, n_not_defined := pmax(n_runs_attempted - n_runs_defined, 0L)]
cell[, aneuploid_fraction_among_defined := fifelse(
  n_runs_defined > 0, n_aneuploid_calls / n_runs_defined, NA_real_
)]
cell[, diploid_fraction_among_defined := fifelse(
  n_runs_defined > 0, n_diploid_calls / n_runs_defined, NA_real_
)]
cell[, stability_class := fcase(
  n_runs_attempted == 0, "NOT_SUBMITTED",
  n_runs_defined >= 2 & aneuploid_fraction_among_defined >= .80,
  "STABLE_ANEUPLOID",
  n_runs_defined >= 2 & diploid_fraction_among_defined >= .80,
  "STABLE_DIPLOID",
  n_runs_defined >= 2, "UNSTABLE_DISCORDANT",
  default = "MOSTLY_NOT_DEFINED"
)]
fwrite(
  cell,
  file.path(out_dir, "GSE154600_copykat_stability_by_cell_v2.csv.gz"),
  compress = "gzip", na = "NA"
)
by_patient <- cell[, .(
  n_final_epithelial = .N,
  n_submitted = sum(n_runs_attempted > 0),
  n_defined_at_least_once = sum(n_runs_defined >= 1),
  n_stable_aneuploid = sum(stability_class == "STABLE_ANEUPLOID"),
  n_stable_diploid = sum(stability_class == "STABLE_DIPLOID"),
  n_unstable = sum(stability_class == "UNSTABLE_DISCORDANT"),
  n_mostly_not_defined = sum(stability_class == "MOSTLY_NOT_DEFINED"),
  n_not_submitted = sum(stability_class == "NOT_SUBMITTED")
), by = .(dataset_id, patient_id)]
by_patient[, `:=`(
  final_count_coverage = n_submitted / n_final_epithelial,
  any_defined_call_rate = n_defined_at_least_once / n_submitted
)]
fwrite(by_patient,
       file.path(out_dir, "GSE154600_copykat_stability_by_patient_v2.csv"),
       na = "NA")
message("Recovered final epithelial counts and affected-patient CopyKAT reruns complete")
