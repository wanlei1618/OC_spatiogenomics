suppressPackageStartupMessages({
  library(data.table)
})

root_dir <- "D:/OC_spatiogenomics/infercnv"
out_dir <- file.path(root_dir, "SPP1_ITGB1_CD44_hypothesis_validation_complete")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
script_dir <- file.path(out_dir, "scripts")
pkg_dir <- file.path(out_dir, "new_package_runs")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pkg_dir, recursive = TRUE, showWarnings = FALSE)

status <- data.frame(tool = character(), analysis = character(), status = character(), note = character())
add_status <- function(tool, analysis, status_value, note) {
  status <<- rbind(status, data.frame(tool = tool, analysis = analysis, status = status_value, note = note, stringsAsFactors = FALSE))
}

focus_pairs <- data.frame(
  ligand = c("SPP1", "SPP1", "APOE", "MIF", "MIF", "TGFB1", "TGFB1", "VEGFA", "VEGFA", "CSF1", "CD47"),
  receptor = c("ITGB1", "CD44", "LRP1", "CD74", "CXCR4", "TGFBR1", "TGFBR2", "KDR", "FLT1", "CSF1R", "SIRPA"),
  stringsAsFactors = FALSE
)

if (requireNamespace("liana", quietly = TRUE)) {
  liana_res <- tryCatch(liana::select_resource("Consensus")[[1]], error = function(e) e)
  if (!inherits(liana_res, "error")) {
    pair_key <- paste(liana_res$source_genesymbol, liana_res$target_genesymbol, sep = "_")
    focus_pairs$liana_consensus_present <- paste(focus_pairs$ligand, focus_pairs$receptor, sep = "_") %in% pair_key
    focus_pairs$liana_consensus_sources <- vapply(seq_len(nrow(focus_pairs)), function(i) {
      hit <- liana_res[liana_res$source_genesymbol == focus_pairs$ligand[i] & liana_res$target_genesymbol == focus_pairs$receptor[i], , drop = FALSE]
      if (nrow(hit) == 0) return("")
      paste(unique(hit$sources), collapse = ";")
    }, character(1))
    fwrite(focus_pairs, file.path(tab_dir, "LIANA_consensus_focus_pair_presence.csv"))
    add_status("LIANA", "Consensus LR resource check for predefined focus axes", "completed",
               paste(sum(focus_pairs$liana_consensus_present), "of", nrow(focus_pairs), "focus pairs present in LIANA Consensus"))
  } else {
    add_status("LIANA", "Consensus LR resource check for predefined focus axes", "failed", conditionMessage(liana_res))
  }
} else {
  add_status("LIANA", "Consensus LR resource check for predefined focus axes", "not_installed", "liana package is not available")
}

if (requireNamespace("nichenetr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(nichenetr))
  data(geneinfo_human, package = "nichenetr")
  genes_to_check <- unique(c(focus_pairs$ligand, focus_pairs$receptor, "EPCAM", "KRT8", "PAX8", "MUC16", "HIF1A", "VEGFA", "FN1", "VIM"))
  nichenet_gene_check <- data.frame(gene = genes_to_check,
                                    present_in_nichenetr_geneinfo = genes_to_check %in% geneinfo_human$symbol)
  fwrite(nichenet_gene_check, file.path(tab_dir, "NicheNet_focus_geneinfo_presence.csv"))
  add_status("NicheNet", "Local package gene universe check for focus ligands/receptors/program genes", "completed",
             paste(sum(nichenet_gene_check$present_in_nichenetr_geneinfo), "of", nrow(nichenet_gene_check), "genes present in nichenetr geneinfo_human"))
} else {
  add_status("NicheNet", "Local package gene universe check for focus ligands/receptors/program genes", "not_installed", "nichenetr package is not available")
}

if (requireNamespace("copykat", quietly = TRUE) && requireNamespace("Seurat", quietly = TRUE)) {
  copykat_note <- tryCatch({
    suppressPackageStartupMessages(library(Seurat))
    suppressPackageStartupMessages(library(copykat))
    meta_path <- file.path(root_dir, "integrated_oc_plan_analysis/tables/integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv")
    meta <- fread(meta_path, data.table = FALSE)
    meta$copykat_group <- ifelse(!is.na(meta$cnv_subclone) & meta$cnv_subclone != "", meta$cnv_subclone,
                                 ifelse(grepl("Myeloid|T/NK|B_", meta$interaction_group), "immune_reference", NA))
    set.seed(20260708)
    sampled_cells <- unlist(lapply(split(meta$cell_integrated_oc, meta$copykat_group), function(v) {
      v <- v[!is.na(v)]
      if (length(v) == 0) return(character())
      sample(v, min(length(v), ifelse(grepl("immune_reference", unique(meta$copykat_group[match(v, meta$cell_integrated_oc)])), 160, 80)))
    }), use.names = FALSE)
    sampled_cells <- unique(sampled_cells)
    obj <- tryCatch(readRDS(file.path(root_dir, "integrated_oc.RData")), error = function(e) {
      env <- new.env()
      load(file.path(root_dir, "integrated_oc.RData"), envir = env)
      env[[ls(env)[1]]]
    })
    sampled_cells <- intersect(sampled_cells, colnames(obj))
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")[, sampled_cells, drop = FALSE]
    keep_genes <- Matrix::rowSums(counts > 0) >= 5
    counts <- counts[keep_genes, , drop = FALSE]
    known <- meta[match(colnames(counts), meta$cell_integrated_oc), c("cell_integrated_oc", "copykat_group", "cnv_subclone", "interaction_group")]
    fwrite(known, file.path(tab_dir, "CopyKAT_sampled_cells_known_groups.csv"))
    oldwd <- getwd()
    setwd(pkg_dir)
    on.exit(setwd(oldwd), add = TRUE)
    ck <- copykat::copykat(rawmat = as.matrix(counts), id.type = "S", ngene.chr = 5,
                           min.gene.per.cell = 200, LOW.DR = 0.05, UP.DR = 0.1,
                           win.size = 25, norm.cell.names = known$cell_integrated_oc[known$copykat_group == "immune_reference"],
                           KS.cut = 0.1, sam.name = "integrated_oc_subclone0204_sampled",
                           distance = "euclidean", output.seg = "FALSE",
                           plot.genes = "FALSE", genome = "hg20", n.cores = 1)
    pred <- NULL
    if (is.list(ck) && "prediction" %in% names(ck)) pred <- ck$prediction
    pred_file <- file.path(pkg_dir, "integrated_oc_subclone0204_sampled_copykat_prediction.csv")
    if (!is.null(pred)) {
      fwrite(as.data.frame(pred), pred_file)
      pred_df <- as.data.frame(pred)
    } else {
      candidates <- list.files(pkg_dir, pattern = "prediction.*\\.txt$|prediction.*\\.csv$", full.names = TRUE)
      if (length(candidates) > 0) pred_df <- fread(candidates[1], data.table = FALSE) else pred_df <- data.frame()
    }
    if (nrow(pred_df) > 0) {
      cell_col <- names(pred_df)[1]
      pred_df$cell_integrated_oc <- pred_df[[cell_col]]
      merged <- merge(pred_df, known, by = "cell_integrated_oc", all.x = TRUE)
      fwrite(merged, file.path(tab_dir, "CopyKAT_sampled_prediction_with_known_subclone.csv"))
      pred_cols <- names(merged)[grepl("prediction|copykat", names(merged), ignore.case = TRUE)]
      pred_col <- if (length(pred_cols)) pred_cols[1] else names(merged)[2]
      summary_tab <- as.data.frame.matrix(table(merged$copykat_group, merged[[pred_col]], useNA = "ifany"))
      summary_tab$copykat_group <- rownames(summary_tab)
      fwrite(summary_tab, file.path(tab_dir, "CopyKAT_sampled_prediction_summary_by_known_group.csv"))
      paste("completed on", ncol(counts), "cells and", nrow(counts), "genes")
    } else {
      paste("copykat returned no parseable prediction table for", ncol(counts), "cells")
    }
  }, error = function(e) paste("failed:", conditionMessage(e)))
  add_status("CopyKAT", "Sampled malignant/CNV validation using immune cells as references", ifelse(grepl("^failed:", copykat_note), "failed", "completed"), copykat_note)
} else {
  add_status("CopyKAT", "Sampled malignant/CNV validation using immune cells as references", "not_installed", "copykat or Seurat package is not available")
}

if (requireNamespace("CaSpER", quietly = TRUE)) {
  add_status("CaSpER", "Secondary CNV validation", "package_available_not_rerun",
             "CaSpER is installed; full run was not executed because gene-position/cytoband annotation generation triggers local biomaRt SSL certificate errors in this R environment")
} else {
  add_status("CaSpER", "Secondary CNV validation", "not_installed", "CaSpER package is not available")
}

fwrite(status, file.path(tab_dir, "new_R_package_focused_validation_status.csv"))
file.copy(normalizePath("run_new_package_focused_validation_checks.R", winslash = "/", mustWork = FALSE),
          file.path(script_dir, "run_new_package_focused_validation_checks.R"), overwrite = TRUE)

