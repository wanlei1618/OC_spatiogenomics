options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
run_root <- file.path(v61, "copykat_runs_v6_1")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v61, "copykat_stability_run_status.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

lineage_names <- c("Epithelial_like", "T_NK_like", "B_Plasma_like", "Myeloid_like")
mats <- setNames(lapply(lineage_names, function(x) readRDS(file.path(
  data_root, "diagnostics_v2", "objects", "GSE154600", "lineage_inputs",
  paste0(x, "_strategy_input.rds")
))$counts), lineage_names)
target_manifest <- fread(file.path(v61, "GSE154600_copykat_target_manifest.csv.gz"))
reference_manifest <- fread(file.path(v61, "GSE154600_copykat_reference_manifest.csv.gz"))

collect_from_manifest <- function(manifest, cells) {
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

worker <- file.path(run_root, "_copykat_worker_v6_1.R")
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
r40 <- "D:/R/R-4.0.3/bin/x64/Rscript.exe"
if (!file.exists(r40)) stop("R 4.0 CopyKAT runtime not found")

run_copykat <- function(patient_id, seed, counts, normal_cells) {
  safe <- paste0("GSE154600__", patient_id, "__seed", seed)
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
      shQuote(safe), shQuote(err), as.character(seed)
    ))
    if (status != 0 || !file.exists(pred)) {
      msg <- if (file.exists(err)) paste(readLines(err, warn = FALSE), collapse = " ") else
        paste0("CopyKAT worker exit status ", status)
      return(list(prediction = NULL, error = msg, path = pred))
    }
  }
  list(prediction = fread(pred), error = NA_character_, path = pred)
}

seeds <- c(20260718L, 20260719L, 20260720L)
patients <- c("T59", "T76", "T77", "T89", "T90")
status_list <- list()
prediction_list <- list()
for (pt in patients) {
  targets <- target_manifest[patient_id == pt & available_for_copykat == TRUE, cell_id]
  refs_all <- reference_manifest[patient_id == pt & available_for_copykat == TRUE, cell_id]
  target_counts <- collect_from_manifest(target_manifest, targets)
  for (seed in seeds) {
    set.seed(seed)
    refs <- if (length(refs_all) > 500L) sample(refs_all, 500L) else refs_all
    reference_counts <- collect_from_manifest(reference_manifest, refs)
    counts <- NULL
    if (!is.null(target_counts) && !is.null(reference_counts)) {
      genes <- intersect(rownames(target_counts), rownames(reference_counts))
      counts <- cbind(
        target_counts[genes, , drop = FALSE],
        reference_counts[genes, , drop = FALSE]
      )
      counts <- counts[, !duplicated(colnames(counts)), drop = FALSE]
    }
    if (is.null(counts) || length(targets) < 20L || length(refs) < 10L) {
      status_list[[paste(pt, seed)]] <- data.table(
        dataset_id = "GSE154600", patient_id = pt, seed = seed,
        n_final_epithelial = target_manifest[patient_id == pt, .N],
        n_targets_available = length(targets), n_targets_submitted = 0L,
        n_reference_available = length(refs_all), n_reference_used = length(refs),
        n_prediction_rows = 0L, n_aneuploid = 0L, n_diploid = 0L,
        n_not_defined = length(targets), run_status = "NOT_EVALUABLE",
        error = "insufficient targets or same-patient immune references",
        prediction_path = NA_character_
      )
      next
    }
    res <- run_copykat(pt, seed, counts, refs)
    if (is.null(res$prediction)) {
      status_list[[paste(pt, seed)]] <- data.table(
        dataset_id = "GSE154600", patient_id = pt, seed = seed,
        n_final_epithelial = target_manifest[patient_id == pt, .N],
        n_targets_available = length(targets),
        n_targets_submitted = ncol(target_counts),
        n_reference_available = length(refs_all), n_reference_used = length(refs),
        n_prediction_rows = 0L, n_aneuploid = 0L, n_diploid = 0L,
        n_not_defined = ncol(target_counts), run_status = "FAILED",
        error = res$error, prediction_path = normalizePath(
          res$path, winslash = "/", mustWork = FALSE)
      )
      next
    }
    p <- res$prediction
    setnames(p, names(p)[1:2], c("cell_id", "copykat_status"))
    p[, copykat_status := tolower(copykat_status)]
    tp <- p[cell_id %in% colnames(target_counts)]
    pred_all <- data.table(cell_id = colnames(target_counts))
    pred_all <- merge(pred_all, tp[, .(cell_id, copykat_status)], by = "cell_id", all.x = TRUE)
    pred_all[is.na(copykat_status), copykat_status := "not.defined"]
    pred_all[, `:=`(
      dataset_id = "GSE154600", patient_id = pt, seed = seed,
      submitted_to_copykat = TRUE
    )]
    prediction_list[[paste(pt, seed)]] <- pred_all
    status_list[[paste(pt, seed)]] <- data.table(
      dataset_id = "GSE154600", patient_id = pt, seed = seed,
      n_final_epithelial = target_manifest[patient_id == pt, .N],
      n_targets_available = length(targets),
      n_targets_submitted = ncol(target_counts),
      n_reference_available = length(refs_all), n_reference_used = length(refs),
      n_prediction_rows = nrow(p),
      n_aneuploid = sum(pred_all$copykat_status == "aneuploid"),
      n_diploid = sum(pred_all$copykat_status == "diploid"),
      n_not_defined = sum(!pred_all$copykat_status %in% c("aneuploid", "diploid")),
      run_status = "COMPLETED", error = NA_character_,
      prediction_path = normalizePath(res$path, winslash = "/", mustWork = FALSE)
    )
    rm(counts, reference_counts); gc()
  }
  rm(target_counts); gc()
}
status <- rbindlist(status_list, fill = TRUE)
predictions <- rbindlist(prediction_list, fill = TRUE)
fwrite(status, out, na = "NA")
fwrite(predictions, file.path(v61, "copykat_stability_raw_predictions.csv.gz"),
       compress = "gzip", na = "NA")
message("GSE154600 three-seed CopyKAT stability runs complete")
