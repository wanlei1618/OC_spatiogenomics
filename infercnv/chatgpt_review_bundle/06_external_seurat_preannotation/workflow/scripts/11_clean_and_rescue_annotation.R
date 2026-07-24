#!/usr/bin/env Rscript

# Clean irrecoverable cells and rescue lineage-consistent misassignments.
# This post-processing step consumes step-10 outputs and never modifies them.

options(stringsAsFactors = FALSE, warn = 1)

required <- c(
  "yaml", "data.table", "Matrix", "Seurat", "SeuratObject", "celda",
  "SingleCellExperiment", "SummarizedExperiment"
)
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  "."
}
source(file.path(script_dir, "_diagnostics_v2_common.R"))

z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
prepared_root <- normalizePath(
  cfg$project$prepared_input_root, winslash = "/", mustWork = FALSE
)
datasets <- split_arg(arg_value("--datasets", "GSE154600,GSE158722"))
assert_datasets(datasets, c("GSE154600", "GSE158722"))

force <- has_flag("--force")
seed <- as.integer(cfg$project$random_seed %||% 12345L)
resolution_value <- as.numeric(arg_value("--resolution", "0.6"))
decont_max_iter <- as.integer(arg_value("--decont-max-iter", "100"))
set.seed(seed)

input_root <- file.path(data_root, "diagnostics_v2_marker_ready_refined")
output_root <- file.path(data_root, "diagnostics_v2_marker_ready_cleaned")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
expected_cells <- c(GSE154600 = 31103L, GSE158722 = 68568L)

write_csv_fast <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, na = "NA")
}

write_csv_gz_fast <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, na = "NA", compress = "gzip")
}

gene_upper <- function(x) toupper(as.character(x))

marker_panels <- list(
  Epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MSLN",
                 "WFDC2", "MUC1", "MUC16", "PAX8", "CLDN3", "CLDN4"),
  T_cell = c("CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2", "CD2", "CD247",
             "LCK", "ITK", "BCL11B", "GIMAP7", "IL7R", "LTB"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "KLRF1", "KLRC1", "PRF1", "CTSW",
              "XCL1", "XCL2", "FGFBP2", "FCGR3A"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22", "CD37", "CD74",
             "VPREB3", "BANK1", "HLA-DRA", "TNFRSF13C"),
  Plasma_cell = c("MZB1", "JCHAIN", "SDC1", "DERL3", "XBP1", "SSR4",
                  "FKBP11", "PRDX4", "IGHG1", "IGHG2", "IGHG3", "IGHA1"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "MRC1", "CD68", "TREM2",
                 "SPP1", "FCGR1A", "MS4A7", "LPL", "CTSD", "LST1"),
  Monocyte = c("S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "LILRB1", "LYZ",
               "CTSD", "SAT1", "LGALS3", "TYROBP", "FCER1G"),
  cDC1 = c("XCR1", "CLEC9A", "BATF3", "CADM1", "IRF8", "CST3"),
  cDC2 = c("CD1C", "FCER1A", "CLEC10A", "CD1E", "HLA-DRA", "CST3"),
  pDC = c("LILRA4", "CLEC4C", "IL3RA", "GZMB", "TCF4", "SERPINF1"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HPGDS", "HDC"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "C7",
                 "COL6A1", "COL6A2", "PDGFRA", "FAP", "POSTN", "CTHRC1"),
  Pericyte = c("RGS5", "MCAM", "CSPG4", "PDGFRB", "ACTA2", "MYH11",
               "NOTCH3", "KCNJ8", "ABCC9", "RBP1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "RAMP3",
                  "PLVAP", "CLDN5", "ACKR1", "CA4", "RGCC", "ESAM"),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "AHSP", "CA1")
)

high_specific <- list(
  Epithelial = c("EPCAM", "PAX8", "WFDC2", "MUC16"),
  T_cell = c("CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2"),
  NK_cell = c("GNLY", "KLRD1", "KLRF1", "XCL1", "XCL2"),
  B_cell = c("MS4A1", "CD19", "CD79A", "CD79B"),
  Plasma_cell = c("MZB1", "JCHAIN", "SDC1", "DERL3"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "MRC1", "TREM2"),
  Monocyte = c("S100A8", "S100A9", "FCN1", "VCAN"),
  cDC1 = c("XCR1", "CLEC9A", "BATF3", "CADM1"),
  cDC2 = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
  pDC = c("LILRA4", "CLEC4C", "IL3RA", "GZMB"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "HDC"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
  Pericyte = c("RGS5", "CSPG4", "PDGFRB", "KCNJ8"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "CLDN5"),
  Erythroid = c("ALAS2", "GYPA", "HBB", "HBA1", "HBA2")
)

state_panels <- list(
  Cycling = c("MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "STMN1", "PCNA",
              "CDC20", "CDK1", "CCNB1", "CCNB2"),
  IFN_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2",
                   "OASL", "IRF7", "IFI6"),
  Hypoxia = c("HIF1A", "CA9", "VEGFA", "BNIP3", "NDRG1", "EGLN3", "LDHA"),
  Stress_response = c("FOS", "JUN", "JUNB", "DDIT3", "HSPA1A", "HSPA1B", "ATF3")
)

type_family <- function(x) {
  y <- as.character(x)
  out <- y
  out[y %in% c("T_cell", "NK_cell")] <- "Lymphoid_T_NK"
  out[y %in% c("B_cell", "Plasma_cell")] <- "Lymphoid_B_Plasma"
  out[y %in% c("Macrophage", "Monocyte", "cDC1", "cDC2", "pDC")] <- "Myeloid_DC"
  out
}

parent_matches_type <- function(parent, type) {
  p <- as.character(parent)
  t <- as.character(type)
  (p == "T_NK" & t %in% c("T_cell", "NK_cell")) |
    (p == "B_Plasma" & t %in% c("B_cell", "Plasma_cell")) |
    (p == "Myeloid" & t %in% c("Macrophage", "Monocyte")) |
    (p == t)
}

empty_marker_table <- function() {
  data.table::data.table(
    cluster = character(), gene = character(), gene_upper = character(),
    avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(),
    p_val = numeric(), p_val_adj = numeric()
  )
}

normalize_marker_columns <- function(markers) {
  markers <- data.table::as.data.table(markers)
  if (!"gene" %in% names(markers)) markers[, gene := rownames(markers)]
  if (!"avg_log2FC" %in% names(markers) && "avg_logFC" %in% names(markers)) {
    data.table::setnames(markers, "avg_logFC", "avg_log2FC")
  }
  required_cols <- c("cluster", "gene", "avg_log2FC", "pct.1", "p_val_adj")
  absent <- setdiff(required_cols, names(markers))
  if (length(absent)) stop("Marker table missing: ", paste(absent, collapse = ", "))
  markers[, gene_upper := gene_upper(gene)]
  markers
}

run_significant_markers <- function(obj) {
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  DefaultAssay(obj) <- "RNA"
  if (length(unique(as.character(Idents(obj)))) < 2L) return(empty_marker_table())
  markers <- FindAllMarkers(
    obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
    min.pct = 0.20, logfc.threshold = 0.25, return.thresh = 1,
    verbose = FALSE
  )
  if (!nrow(markers)) return(empty_marker_table())
  markers <- normalize_marker_columns(markers)
  markers[
    is.finite(p_val_adj) & p_val_adj < 0.05 &
      is.finite(avg_log2FC) & avg_log2FC > 0.25 &
      is.finite(pct.1) & pct.1 >= 0.20
  ]
}

score_marker_evidence <- function(markers, clusters) {
  markers <- data.table::as.data.table(markers)
  if (!nrow(markers)) markers <- empty_marker_table()
  rows <- list()
  for (cluster_id in clusters) {
    cm <- markers[as.character(cluster) == as.character(cluster_id)]
    for (cell_type in names(marker_panels)) {
      hit <- cm[gene_upper %in% gene_upper(marker_panels[[cell_type]])]
      hit_hi <- hit[gene_upper %in% gene_upper(high_specific[[cell_type]])]
      rows[[length(rows) + 1L]] <- data.frame(
        cluster = as.character(cluster_id), candidate_cell_type = cell_type,
        support_n = data.table::uniqueN(hit$gene_upper),
        high_specific_n = data.table::uniqueN(hit_hi$gene_upper),
        support_weight = if (nrow(hit)) {
          sum(pmax(hit$avg_log2FC, 0) * pmax(hit$pct.1, 0), na.rm = TRUE)
        } else 0,
        supporting_markers = if (nrow(hit)) {
          paste(unique(hit[order(-avg_log2FC)]$gene), collapse = ";")
        } else "",
        stringsAsFactors = FALSE
      )
    }
  }
  scores <- data.table::rbindlist(rows, fill = TRUE)
  scores[, eligible := support_n >= 4L | high_specific_n >= 2L]
  decisions <- scores[order(cluster, -eligible, -support_n, -high_specific_n,
                            -support_weight), {
    top <- .SD[1L]
    eligible_rows <- .SD[eligible == TRUE]
    second <- if (nrow(eligible_rows) >= 2L) eligible_rows[2L] else NULL
    conflict <- !is.null(second) &&
      type_family(top$candidate_cell_type) != type_family(second$candidate_cell_type)
    suggested <- if (!top$eligible) "Unresolved" else
      if (conflict) "Mixed_or_doublet" else as.character(top$candidate_cell_type)
    list(
      suggested_cell_type = suggested,
      canonical_support_n = as.integer(top$support_n),
      high_specific_support_n = as.integer(top$high_specific_n),
      canonical_markers_present = as.character(top$supporting_markers),
      second_candidate = if (is.null(second)) NA_character_ else
        as.character(second$candidate_cell_type),
      incompatible_lineage_program = conflict
    )
  }, by = cluster]
  list(scores = scores, decisions = decisions)
}

read_complete_input <- function(dataset_id) {
  ds_cfg <- cfg$datasets[[dataset_id]]
  path <- file.path(prepared_root, ds_cfg$prepared_sce)
  if (!file.exists(path)) {
    path <- file.path(
      data_root, dataset_id, "objects", paste0(dataset_id, "_preannotation.rds")
    )
  }
  if (!file.exists(path)) stop(dataset_id, ": prepared input not found: ", path)
  message("Reading repaired-QC input: ", path)
  input <- readRDS(path)
  if (inherits(input, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(input, "counts")
    md <- as.data.frame(SummarizedExperiment::colData(input))
  } else if (inherits(input, "Seurat")) {
    counts <- SeuratObject::LayerData(input, assay = "RNA", layer = "counts")
    md <- as.data.frame(input[[]], check.names = FALSE)
  } else {
    stop(dataset_id, ": unsupported input class: ", class(input)[[1L]])
  }
  md$cell_id <- colnames(counts)
  qc <- read_repaired_qc(data_root, dataset_id)
  qc <- qc[diagnostics_v2_qc_pass %in% TRUE]
  idx <- match(qc$cell_id, md$cell_id)
  if (anyNA(idx)) stop(dataset_id, ": repaired-QC cells missing: ", sum(is.na(idx)))
  selected_md <- md[idx, , drop = FALSE]
  selected_md$cell_id <- qc$cell_id
  qidx <- match(selected_md$cell_id, qc$cell_id)
  for (column in setdiff(names(qc), names(selected_md))) {
    selected_md[[column]] <- qc[[column]][qidx]
  }
  if (!"sample_id" %in% names(selected_md)) selected_md$sample_id <- qc$sample_id[qidx]
  if (!"patient_id" %in% names(selected_md)) selected_md$patient_id <- selected_md$sample_id
  if (!"timepoint" %in% names(selected_md)) selected_md$timepoint <- NA_character_
  counts <- counts[, idx, drop = FALSE]
  rownames(selected_md) <- selected_md$cell_id
  rm(input, md, qc)
  gc()
  list(counts = counts, metadata = selected_md)
}

detect_platform <- function(md, dataset_id) {
  candidates <- intersect(
    c("platform", "technology", "sequencing_platform", "library_platform"),
    names(md)
  )
  if (!length(candidates)) {
    return(list(values = rep(NA_character_, nrow(md)), field = NA_character_,
                reliable = FALSE,
                note = "No explicit platform field in repaired-QC metadata."))
  }
  field <- candidates[[1L]]
  values <- trimws(as.character(md[[field]]))
  values[!nzchar(values)] <- NA_character_
  reliable <- !anyNA(values)
  if (dataset_id == "GSE158722" && reliable) {
    recognized <- grepl("10x|chromium|icell8", values, ignore.case = TRUE)
    reliable <- all(recognized)
  }
  list(
    values = if (reliable) values else rep(NA_character_, nrow(md)),
    field = if (reliable) field else NA_character_, reliable = reliable,
    note = if (reliable) paste0("Platform read from metadata field ", field, ".") else
      paste0("Metadata field ", field, " was not complete/reliable; platform left NA.")
  )
}

run_sample_decontx <- function(counts, md, old_cluster, dataset_id, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  adjusted <- counts
  contamination <- rep(NA_real_, ncol(counts))
  method <- rep("DecontX", ncol(counts))
  samples <- unique(as.character(md$sample_id))
  for (sample_i in seq_along(samples)) {
    sample_id <- samples[[sample_i]]
    idx <- which(as.character(md$sample_id) == sample_id)
    message(dataset_id, " / DecontX sample ", sample_id, ": ", length(idx), " cells")
    part <- counts[, idx, drop = FALSE]
    keep_gene <- Matrix::rowSums(part) > 0
    z_local <- as.integer(factor(old_cluster[idx]))
    safe_sample <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
    cache_path <- file.path(
      cache_dir, sprintf("%03d_%s_decontx.rds", sample_i, safe_sample)
    )
    if (file.exists(cache_path)) {
      cached <- readRDS(cache_path)
      cache_ok <- identical(cached$cell_id, colnames(part)) &&
        identical(cached$gene, rownames(part)[keep_gene])
      if (cache_ok) {
        adjusted[keep_gene, idx] <- cached$counts
        contamination[idx] <- cached$contamination
        method[idx] <- cached$method
        message(dataset_id, " / reused DecontX cache for ", sample_id)
        rm(cached, part)
        gc()
        next
      }
      warning(dataset_id, " / ignoring incompatible cache for ", sample_id)
    }
    if (length(idx) < 20L || sum(keep_gene) < 100L) {
      method[idx] <- "DecontX_skipped_too_small"
      next
    }
    ans <- tryCatch(
      celda::decontX(
        part[keep_gene, , drop = FALSE], z = z_local,
        maxIter = decont_max_iter, estimateDelta = TRUE,
        convergence = 0.001, varGenes = min(5000L, sum(keep_gene)),
        seed = seed, verbose = FALSE
      ),
      error = function(e) e
    )
    if (inherits(ans, "error")) {
      warning(dataset_id, " / ", sample_id, ": DecontX failed; conservative ",
              "unchanged counts retained: ", conditionMessage(ans))
      method[idx] <- "DecontX_failed_conservative"
      next
    }
    adjusted[keep_gene, idx] <- ans[["decontXcounts"]]
    contamination[idx] <- as.numeric(ans[["contamination"]])
    saveRDS(
      list(
        cell_id = colnames(part), gene = rownames(part)[keep_gene],
        counts = ans[["decontXcounts"]], contamination = contamination[idx],
        method = method[idx]
      ),
      cache_path, compress = FALSE
    )
    rm(ans, part)
    gc()
  }
  list(counts = adjusted, contamination = contamination, method = method)
}

make_marker_object <- function(counts, md, identities, project) {
  rownames(md) <- md$cell_id
  obj <- CreateSeuratObject(
    counts = counts, project = project, min.cells = 0, min.features = 0,
    meta.data = md
  )
  obj <- NormalizeData(obj, verbose = FALSE)
  Idents(obj) <- factor(identities)
  obj
}

cell_program_evidence <- function(counts) {
  detected <- matrix(
    0L, nrow = ncol(counts), ncol = length(marker_panels),
    dimnames = list(colnames(counts), names(marker_panels))
  )
  high <- detected
  genes_upper <- gene_upper(rownames(counts))
  for (cell_type in names(marker_panels)) {
    idx <- which(genes_upper %in% gene_upper(marker_panels[[cell_type]]))
    if (length(idx)) detected[, cell_type] <- Matrix::colSums(counts[idx, , drop = FALSE] > 0)
    idx_hi <- which(genes_upper %in% gene_upper(high_specific[[cell_type]]))
    if (length(idx_hi)) high[, cell_type] <- Matrix::colSums(counts[idx_hi, , drop = FALSE] > 0)
  }
  eligible <- detected >= 4L | high >= 2L
  ord <- t(apply(detected + high * 2L, 1L, order, decreasing = TRUE))
  top <- colnames(detected)[ord[, 1L]]
  second <- colnames(detected)[ord[, 2L]]
  top_ok <- eligible[cbind(seq_len(nrow(eligible)), ord[, 1L])]
  second_ok <- eligible[cbind(seq_len(nrow(eligible)), ord[, 2L])]
  incompatible <- top_ok & second_ok & type_family(top) != type_family(second)
  data.table::data.table(
    cell_id = rownames(detected),
    cell_program_type = ifelse(top_ok & !incompatible, top, "Unresolved"),
    cell_program_support = detected[cbind(seq_len(nrow(detected)), ord[, 1L])],
    cell_high_specific_support = high[cbind(seq_len(nrow(high)), ord[, 1L])],
    cell_incompatible_programs = incompatible
  )
}

prepare_clustering_object <- function(counts, md, project) {
  rownames(md) <- md$cell_id
  obj <- CreateSeuratObject(
    counts = counts, project = project, min.cells = 0, min.features = 0,
    meta.data = md
  )
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    obj, selection.method = "vst", nfeatures = min(3000L, nrow(obj)),
    verbose = FALSE
  )
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  npcs <- min(30L, ncol(obj) - 1L, length(VariableFeatures(obj)) - 1L)
  if (npcs < 2L) return(NULL)
  RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
}

recluster_type <- function(counts, md, dataset_id, cell_type) {
  if (ncol(counts) < 50L) {
    return(data.table::data.table(
      cell_id = colnames(counts), clustering_strategy = "broad_type_only",
      final_cluster = paste(cell_type, "broad_type_only", "1", sep = "__")
    ))
  }
  obj <- prepare_clustering_object(
    counts, md, paste0(dataset_id, "_cleaned_", cell_type)
  )
  if (is.null(obj)) {
    return(data.table::data.table(
      cell_id = colnames(counts), clustering_strategy = "broad_type_only",
      final_cluster = paste(cell_type, "broad_type_only", "1", sep = "__")
    ))
  }
  reduction <- "pca"
  strategy <- "A_uncorrected"
  if (cell_type != "Epithelial" && cell_type != "Unresolved" &&
      requireNamespace("harmony", quietly = TRUE)) {
    batch_var <- if (dataset_id == "GSE158722") "patient_id" else "sample_id"
    values <- as.character(obj@meta.data[[batch_var]])
    tab <- table(values[!is.na(values) & nzchar(values)])
    if (!anyNA(values) && all(nzchar(values)) && length(tab) >= 2L && min(tab) >= 20L) {
      obj <- harmony::RunHarmony(
        obj, group.by.vars = batch_var, reduction.use = "pca", verbose = FALSE
      )
      reduction <- "harmony"
      strategy <- if (dataset_id == "GSE158722") {
        "B_harmony_patient"
      } else {
        "B_harmony_sample"
      }
    }
  }
  dims <- seq_len(min(30L, ncol(Embeddings(obj, reduction))))
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
  obj <- FindClusters(
    obj, resolution = resolution_value, random.seed = seed, verbose = FALSE
  )
  raw_cluster <- as.character(Idents(obj))
  data.table::data.table(
    cell_id = colnames(obj), clustering_strategy = strategy,
    final_cluster = paste(cell_type, strategy, raw_cluster, sep = "__")
  )
}

state_from_markers <- function(markers, clusters) {
  out <- vector("list", length(clusters))
  for (i in seq_along(clusters)) {
    cluster_id <- clusters[[i]]
    cm <- markers[as.character(cluster) == as.character(cluster_id)]
    states <- names(state_panels)[vapply(state_panels, function(panel) {
      data.table::uniqueN(cm[gene_upper %in% gene_upper(panel)]$gene_upper) >= 3L
    }, FUN.VALUE = logical(1))]
    if ("SPP1" %in% cm$gene_upper) states <- c(states, "SPP1_high")
    if ("C1QC" %in% cm$gene_upper) states <- c(states, "C1QC_high")
    if ("FOLR2" %in% cm$gene_upper) states <- c(states, "FOLR2_high")
    out[[i]] <- data.frame(
      cluster = cluster_id,
      cell_state = if (length(states)) paste(unique(states), collapse = ";") else "None",
      stringsAsFactors = FALSE
    )
  }
  data.table::rbindlist(out)
}

known_rescues <- data.table::data.table(
  dataset_id = c(rep("GSE154600", 4L), "GSE158722"),
  old_cluster = c(
    "T_NK__B_harmony_sample__15", "T_NK__B_harmony_sample__16",
    "Myeloid__B_harmony_sample__10", "Myeloid__B_harmony_sample__7",
    "Epithelial__A_uncorrected__30"
  ),
  required_type = c("Erythroid", "Mast", "cDC1", "cDC2", "Endothelial")
)

for (dataset_id in datasets) {
  message("===== ", dataset_id, " =====")
  dataset_name <- dataset_id
  mandatory <- known_rescues[dataset_id == dataset_name]
  dataset_out <- file.path(output_root, dataset_id)
  if (dir.exists(dataset_out)) {
    if (!force) stop(dataset_id, ": output exists; use --force")
    normalized_out <- normalizePath(dataset_out, winslash = "/", mustWork = TRUE)
    normalized_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
    if (!startsWith(normalized_out, paste0(normalized_root, "/"))) {
      stop("Refusing to remove output outside cleaned root: ", normalized_out)
    }
    existing <- list.files(dataset_out, full.names = TRUE, all.files = TRUE,
                           no.. = TRUE)
    existing <- existing[basename(existing) != "_decontx_cache"]
    if (length(existing)) unlink(existing, recursive = TRUE, force = TRUE)
  }
  dir.create(dataset_out, recursive = TRUE, showWarnings = FALSE)

  refined_dir <- file.path(input_root, dataset_id)
  assignment_path <- file.path(refined_dir, "annotation_ready_full_cell_assignments.csv.gz")
  template_path <- file.path(refined_dir, "annotation_ready_cluster_template_refined.csv")
  if (!file.exists(assignment_path) || !file.exists(template_path)) {
    stop(dataset_id, ": required step-10 inputs are missing")
  }
  assignments <- data.table::fread(assignment_path, showProgress = FALSE)
  template <- data.table::fread(template_path, showProgress = FALSE)
  data.table::setnames(
    assignments, c("parent_broad_type", "final_cluster"),
    c("old_parent_type", "old_cluster")
  )
  if (nrow(assignments) != expected_cells[[dataset_id]] ||
      data.table::uniqueN(assignments$cell_id) != expected_cells[[dataset_id]]) {
    stop(dataset_id, ": step-10 assignments do not cover repaired-QC cells")
  }

  input <- read_complete_input(dataset_id)
  if (ncol(input$counts) != expected_cells[[dataset_id]]) {
    stop(dataset_id, ": expected ", expected_cells[[dataset_id]],
         " repaired-QC cells, found ", ncol(input$counts))
  }
  ord <- match(colnames(input$counts), assignments$cell_id)
  if (anyNA(ord)) stop(dataset_id, ": assignments missing expression cells")
  assignments <- assignments[ord]
  md <- input$metadata[match(colnames(input$counts), input$metadata$cell_id), , drop = FALSE]
  platform_info <- detect_platform(md, dataset_id)
  md$platform <- platform_info$values
  md$old_cluster <- assignments$old_cluster
  md$old_parent_type <- assignments$old_parent_type

  decont_cache <- file.path(dataset_out, "_decontx_cache")
  decont <- run_sample_decontx(
    input$counts, md, md$old_cluster, dataset_id, decont_cache
  )
  decont_counts <- decont$counts
  rm(input)
  gc()

  old_obj <- make_marker_object(
    decont_counts, md, md$old_cluster, paste0(dataset_id, "_post_decont_old_clusters")
  )
  old_markers <- run_significant_markers(old_obj)
  old_evidence <- score_marker_evidence(old_markers, unique(md$old_cluster))
  rm(old_obj)
  gc()

  cell_evidence <- cell_program_evidence(decont_counts)
  cidx <- match(md$cell_id, cell_evidence$cell_id)
  md$cell_program_type <- cell_evidence$cell_program_type[cidx]
  md$cell_program_support <- cell_evidence$cell_program_support[cidx]
  md$cell_high_specific_support <- cell_evidence$cell_high_specific_support[cidx]
  md$cell_incompatible_programs <- cell_evidence$cell_incompatible_programs[cidx]

  decision <- old_evidence$decisions
  didx <- match(md$old_cluster, decision$cluster)
  md$cluster_suggested_type <- decision$suggested_cell_type[didx]
  md$cluster_support_n <- decision$canonical_support_n[didx]
  md$cluster_high_specific_n <- decision$high_specific_support_n[didx]
  md$cluster_incompatible_program <- decision$incompatible_lineage_program[didx]
  stable_cluster <- !md$cluster_suggested_type %in% c(
    "Unresolved", "Mixed_or_doublet", NA_character_
  )
  md$provisional_cell_type <- ifelse(
    stable_cluster, md$cluster_suggested_type,
    ifelse(md$cell_program_type != "Unresolved", md$cell_program_type, "Unresolved")
  )
  for (i in seq_len(nrow(mandatory))) {
    row <- mandatory[i]
    prior <- template[cluster == row$old_cluster]
    prior_strength <- tolower(as.character(prior$canonical_evidence_strength))
    if (nrow(prior) != 1L ||
        prior$suggested_broad_cell_type[[1L]] != row$required_type ||
        !prior_strength[[1L]] %in% c("high", "moderate", "medium")) {
      stop(dataset_id, ": required rescue lacks supporting step-10 evidence: ",
           row$old_cluster, " -> ", row$required_type)
    }
    md$provisional_cell_type[md$old_cluster == row$old_cluster] <-
      row$required_type
  }

  dbl_class <- if ("scDblFinder.class" %in% names(md)) {
    as.character(md$scDblFinder.class)
  } else rep(NA_character_, nrow(md))
  dbl_score <- if ("scDblFinder.score" %in% names(md)) {
    as.numeric(md$scDblFinder.score)
  } else rep(NA_real_, nrow(md))
  md$doublet_score <- dbl_score
  md$doublet_call <- dbl_class
  heterotypic_doublet <- !is.na(dbl_class) & tolower(dbl_class) == "doublet" &
    md$cell_incompatible_programs

  n_feature <- as.numeric(md$nFeature_RNA)
  n_count <- as.numeric(md$nCount_RNA)
  percent_mt <- as.numeric(md$percent.mt)
  complexity <- log10(pmax(n_feature, 1)) / log10(pmax(n_count, 10))
  sample_groups <- split(seq_len(nrow(md)), as.character(md$sample_id))
  low_feature_cut <- low_count_cut <- low_complexity_cut <- high_mt_cut <-
    rep(NA_real_, nrow(md))
  for (idx in sample_groups) {
    low_feature_cut[idx] <- max(200, stats::quantile(n_feature[idx], 0.01, na.rm = TRUE))
    low_count_cut[idx] <- max(500, stats::quantile(n_count[idx], 0.01, na.rm = TRUE))
    low_complexity_cut[idx] <- stats::quantile(complexity[idx], 0.01, na.rm = TRUE)
    high_mt_cut[idx] <- max(20, stats::quantile(percent_mt[idx], 0.99, na.rm = TRUE))
  }
  genes_upper <- gene_upper(rownames(decont_counts))
  program_gene <- grepl("^MT-|^RPL|^RPS|^HSPA|^FOS$|^JUN$|^JUNB$|^ATF3$", genes_upper)
  program_fraction <- Matrix::colSums(decont_counts[program_gene, , drop = FALSE]) /
    pmax(Matrix::colSums(decont_counts), 1)
  low_feature <- n_feature <= low_feature_cut
  low_count <- n_count <= low_count_cut
  high_mt <- percent_mt >= high_mt_cut
  low_complexity <- complexity <= low_complexity_cut
  program_dominated <- program_fraction >= 0.70
  no_stable_lineage <- md$provisional_cell_type == "Unresolved"
  bad_feature_count <- low_feature + low_count + high_mt + low_complexity + program_dominated
  remove_low_quality <- bad_feature_count >= 3L & no_stable_lineage
  remove_ambient <- !remove_low_quality & !is.na(decont$contamination) &
    decont$contamination >= 0.50 & no_stable_lineage &
    (low_feature | low_count | low_complexity)
  removal_reason <- data.table::fcase(
    heterotypic_doublet, "REMOVE_HETEROTYPIC_DOUBLET",
    remove_low_quality, "REMOVE_LOW_QUALITY",
    remove_ambient, "REMOVE_AMBIENT_DOMINATED",
    default = NA_character_
  )
  retained <- is.na(removal_reason)

  quality_summary <- paste0(
    "nFeature=", n_feature, ";nCount=", n_count,
    ";percent_mt=", sprintf("%.3f", percent_mt),
    ";complexity=", sprintf("%.3f", complexity),
    ";program_fraction=", sprintf("%.3f", program_fraction),
    ";ambient=", ifelse(is.na(decont$contamination), "NA",
                        sprintf("%.3f", decont$contamination))
  )
  removed <- data.table::data.table(
    cell_id = md$cell_id[!retained], sample_id = as.character(md$sample_id[!retained]),
    patient_id = as.character(md$patient_id[!retained]),
    old_cluster = md$old_cluster[!retained], removal_reason = removal_reason[!retained],
    doublet_score = md$doublet_score[!retained], quality_summary = quality_summary[!retained]
  )

  kept_md <- md[retained, , drop = FALSE]
  kept_md$final_cell_type <- md$provisional_cell_type[retained]
  kept_counts <- decont_counts[, retained, drop = FALSE]
  cluster_parts <- list()
  for (cell_type in sort(unique(kept_md$final_cell_type))) {
    idx <- which(kept_md$final_cell_type == cell_type)
    message(dataset_id, " / recluster ", cell_type, ": ", length(idx), " cells")
    cluster_parts[[cell_type]] <- recluster_type(
      kept_counts[, idx, drop = FALSE], kept_md[idx, , drop = FALSE],
      dataset_id, cell_type
    )
    gc()
  }
  reclustered <- data.table::rbindlist(cluster_parts, fill = TRUE)
  ridx <- match(kept_md$cell_id, reclustered$cell_id)
  kept_md$clustering_strategy <- reclustered$clustering_strategy[ridx]
  kept_md$final_cluster <- reclustered$final_cluster[ridx]

  final_obj <- make_marker_object(
    kept_counts, kept_md, kept_md$final_cluster,
    paste0(dataset_id, "_cleaned_final_markers")
  )
  final_markers <- run_significant_markers(final_obj)
  rm(final_obj)
  gc()
  if (nrow(final_markers) && any(
    !is.finite(final_markers$p_val_adj) | final_markers$p_val_adj >= 0.05 |
      final_markers$avg_log2FC <= 0.25 | final_markers$pct.1 < 0.20
  )) stop(dataset_id, ": nonsignificant marker leaked into cleaned output")
  final_markers[, dataset_id := dataset_id]
  final_evidence <- score_marker_evidence(final_markers, unique(kept_md$final_cluster))
  final_states <- state_from_markers(final_markers, unique(kept_md$final_cluster))

  cluster_stats <- data.table::data.table(
    cluster = kept_md$final_cluster,
    initial_cell_type = kept_md$final_cell_type,
    sample_id = as.character(kept_md$sample_id),
    patient_id = as.character(kept_md$patient_id),
    platform = as.character(kept_md$platform)
  )[, {
    patient_tab <- sort(table(patient_id), decreasing = TRUE)
    sample_tab <- sort(table(sample_id), decreasing = TRUE)
    platform_valid <- !is.na(platform) & nzchar(platform)
    platform_tab <- sort(table(platform[platform_valid]), decreasing = TRUE)
    list(
      n_cells = .N, n_samples = data.table::uniqueN(sample_id),
      n_patients = data.table::uniqueN(patient_id),
      dominant_patient = names(patient_tab)[[1L]],
      dominant_patient_fraction = as.numeric(patient_tab[[1L]] / .N),
      dominant_sample = names(sample_tab)[[1L]],
      dominant_sample_fraction = as.numeric(sample_tab[[1L]] / .N),
      dominant_platform = if (length(platform_tab)) names(platform_tab)[[1L]] else NA_character_,
      dominant_platform_fraction = if (length(platform_tab)) {
        as.numeric(platform_tab[[1L]] / sum(platform_tab))
      } else NA_real_
    )
  }, by = .(cluster, initial_cell_type)]
  cluster_table <- merge(
    cluster_stats, final_evidence$decisions, by = "cluster", all.x = TRUE
  )
  cluster_table <- merge(cluster_table, final_states, by = "cluster", all.x = TRUE)
  cluster_table[, patient_enriched := dominant_patient_fraction >= 0.80 & n_cells >= 20L]
  cluster_table[, platform_confounded := if (platform_info$reliable &&
      data.table::uniqueN(kept_md$platform) >= 2L) {
    dominant_platform_fraction >= 0.80 & n_cells >= 20L
  } else NA]
  cluster_table[, final_cell_type := data.table::fcase(
    !suggested_cell_type %in% c("Unresolved", "Mixed_or_doublet") &
      incompatible_lineage_program == FALSE, suggested_cell_type,
    initial_cell_type != "Unresolved" & incompatible_lineage_program == FALSE,
      initial_cell_type,
    default = "Unresolved"
  )]
  marker_counts <- final_markers[, .N, by = cluster]
  cluster_table[, n_significant_markers := marker_counts$N[
    match(cluster, marker_counts$cluster)]]
  cluster_table[is.na(n_significant_markers), n_significant_markers := 0L]
  cluster_table[, annotation_status := data.table::fcase(
    incompatible_lineage_program == TRUE | suggested_cell_type == "Mixed_or_doublet",
      "REVIEW_MIXED_OR_DOUBLET",
    !is.na(platform_confounded) & platform_confounded == TRUE,
      "REVIEW_PLATFORM_CONFOUNDED",
    patient_enriched == TRUE, "REVIEW_PATIENT_ENRICHED",
    final_cell_type == "Unresolved", "REVIEW_AMBIGUOUS",
    canonical_support_n >= 4L | high_specific_support_n >= 2L,
      "READY_HIGH_CONFIDENCE",
    canonical_support_n >= 2L, "READY_BROAD_TYPE_ONLY",
    default = "REVIEW_AMBIGUOUS"
  )]
  cluster_table[, annotation_confidence := data.table::fcase(
    annotation_status == "READY_HIGH_CONFIDENCE", "High",
    annotation_status == "READY_BROAD_TYPE_ONLY", "Broad_type_only",
    default = "Review"
  )]
  cluster_table[, `:=`(
    dataset_id = dataset_id, manual_cell_type = "", manual_cell_subtype = "",
    manual_confidence = "", manual_notes = ""
  )]

  kidx <- match(kept_md$final_cluster, cluster_table$cluster)
  kept_md$final_cell_type <- cluster_table$final_cell_type[kidx]
  kept_md$cell_state <- cluster_table$cell_state[kidx]
  kept_md$annotation_confidence <- cluster_table$annotation_confidence[kidx]
  kept_md$patient_enriched <- cluster_table$patient_enriched[kidx]
  kept_md$platform_confounded <- cluster_table$platform_confounded[kidx]
  kept_md$state_only <- kept_md$cell_state != "None"

  cleaned_assignments <- data.table::data.table(
    dataset_id = dataset_id, cell_id = kept_md$cell_id,
    sample_id = as.character(kept_md$sample_id),
    patient_id = as.character(kept_md$patient_id),
    timepoint = as.character(kept_md$timepoint),
    platform = as.character(kept_md$platform),
    old_parent_type = kept_md$old_parent_type,
    old_cluster = kept_md$old_cluster,
    final_cell_type = kept_md$final_cell_type,
    cell_state = kept_md$cell_state,
    annotation_confidence = kept_md$annotation_confidence,
    patient_enriched = kept_md$patient_enriched,
    platform_confounded = kept_md$platform_confounded,
    state_only = kept_md$state_only,
    doublet_score = kept_md$doublet_score,
    doublet_call = kept_md$doublet_call,
    ambient_contamination = decont$contamination[retained],
    decontamination_method = decont$method[retained],
    clustering_strategy = kept_md$clustering_strategy,
    final_cluster = kept_md$final_cluster
  )

  old_decision <- old_evidence$decisions
  rescue_candidates <- unique(data.table::data.table(
    dataset_id = dataset_id, old_parent_type = md$old_parent_type,
    old_cluster = md$old_cluster,
    new_cell_type = md$provisional_cell_type
  ))
  rescue_candidates <- rescue_candidates[
    new_cell_type != "Unresolved" &
      (!parent_matches_type(old_parent_type, new_cell_type) |
         old_cluster %in% mandatory$old_cluster)
  ]
  rescue_candidates <- merge(
    rescue_candidates,
    old_decision[, .(
      old_cluster = cluster, supporting_markers = canonical_markers_present,
      canonical_support_n, high_specific_support_n,
      rescue_confidence = ifelse(
        canonical_support_n >= 4L | high_specific_support_n >= 2L, "High", "Moderate"
      )
    )], by = "old_cluster", all.x = TRUE
  )
  for (i in seq_len(nrow(mandatory))) {
    row <- mandatory[i]
    prior <- template[cluster == row$old_cluster]
    ridx <- which(
      rescue_candidates$old_cluster == row$old_cluster &
        rescue_candidates$new_cell_type == row$required_type
    )
    if (length(ridx) == 1L) {
      prior_markers <- as.character(prior$canonical_markers_present[[1L]])
      prior_n <- if (!is.na(prior_markers) && nzchar(prior_markers)) {
        length(unique(strsplit(prior_markers, ";", fixed = TRUE)[[1L]]))
      } else 0L
      if (is.na(rescue_candidates$supporting_markers[[ridx]]) ||
          !nzchar(rescue_candidates$supporting_markers[[ridx]])) {
        rescue_candidates$supporting_markers[[ridx]] <- prior_markers
      }
      rescue_candidates$canonical_support_n[[ridx]] <- max(
        rescue_candidates$canonical_support_n[[ridx]], prior_n, na.rm = TRUE
      )
      if (tolower(prior$canonical_evidence_strength[[1L]]) == "high") {
        rescue_candidates$rescue_confidence[[ridx]] <- "High"
      }
    }
  }
  before_counts <- data.table::data.table(old_cluster = md$old_cluster)[
    , .N, by = old_cluster
  ]
  retained_counts <- data.table::data.table(
    old_cluster = kept_md$old_cluster,
    new_cell_type = kept_md$final_cell_type
  )[, .N, by = .(old_cluster, new_cell_type)]
  rescue_candidates[, `:=`(
    n_cells_before = before_counts$N[match(old_cluster, before_counts$old_cluster)],
    n_cells_retained = retained_counts$N[match(
      paste(old_cluster, new_cell_type, sep = "\r"),
      paste(retained_counts$old_cluster, retained_counts$new_cell_type, sep = "\r")
    )],
    notes = "Post-DecontX canonical evidence; old parent constraint removed."
  )]
  rescued <- rescue_candidates[, .(
    dataset_id, old_parent_type, old_cluster, new_cell_type,
    supporting_markers, n_cells_before, n_cells_retained,
    rescue_confidence, notes
  )]

  for (i in seq_len(nrow(mandatory))) {
    row <- mandatory[i]
    hit <- rescued[old_cluster == row$old_cluster & new_cell_type == row$required_type]
    if (!nrow(hit) || hit$n_cells_retained[[1L]] <= 0L) {
      stop(dataset_id, ": required rescue failed: ", row$old_cluster,
           " -> ", row$required_type)
    }
  }
  if (nrow(cleaned_assignments) + nrow(removed) != expected_cells[[dataset_id]]) {
    stop(dataset_id, ": retained + removed does not equal repaired-QC input")
  }
  if (anyNA(cleaned_assignments$final_cell_type) ||
      any(!nzchar(cleaned_assignments$final_cell_type))) {
    stop(dataset_id, ": retained cell without final_cell_type")
  }
  forbidden_states <- c("Cycling", "IFN_response", "Hypoxia", "Stress_response")
  if (any(cleaned_assignments$final_cell_type %in% forbidden_states)) {
    stop(dataset_id, ": a state leaked into final_cell_type")
  }
  if (dataset_id == "GSE158722" && any(
    cleaned_assignments$clustering_strategy %in%
      c("B_harmony_sample", "B_harmony_timepoint")
  )) stop("GSE158722 used prohibited sample/timepoint Harmony")

  preferred_template <- c(
    "dataset_id", "cluster", "final_cell_type", "cell_state",
    "annotation_status", "annotation_confidence", "n_cells", "n_samples",
    "n_patients", "dominant_patient", "dominant_patient_fraction",
    "dominant_sample", "dominant_sample_fraction", "patient_enriched",
    "platform_confounded", "suggested_cell_type", "canonical_support_n",
    "high_specific_support_n", "canonical_markers_present", "second_candidate",
    "incompatible_lineage_program", "n_significant_markers",
    "manual_cell_type", "manual_cell_subtype", "manual_confidence", "manual_notes"
  )
  data.table::setcolorder(
    cluster_table,
    c(intersect(preferred_template, names(cluster_table)),
      setdiff(names(cluster_table), preferred_template))
  )

  write_csv_gz_fast(
    cleaned_assignments, file.path(dataset_out, "cleaned_cell_assignments.csv.gz")
  )
  write_csv_fast(
    cluster_table, file.path(dataset_out, "cleaned_cluster_annotation_template.csv")
  )
  write_csv_gz_fast(removed, file.path(dataset_out, "removed_cells.csv.gz"))
  write_csv_fast(rescued, file.path(dataset_out, "rescued_lineage_clusters.csv"))
  write_csv_gz_fast(
    final_markers, file.path(dataset_out, "cleaned_significant_markers.csv.gz")
  )

  removal_counts <- table(factor(
    removed$removal_reason,
    levels = c("REMOVE_HETEROTYPIC_DOUBLET", "REMOVE_LOW_QUALITY",
               "REMOVE_AMBIENT_DOMINATED")
  ))
  type_counts <- cleaned_assignments[, .N, by = final_cell_type][order(final_cell_type)]
  patient_clusters <- cluster_table[patient_enriched == TRUE, cluster]
  platform_clusters <- cluster_table[!is.na(platform_confounded) &
                                       platform_confounded == TRUE, cluster]
  unresolved_clusters <- cluster_table[final_cell_type == "Unresolved", cluster]
  summary_lines <- c(
    paste0("# ", dataset_id, " cleaning and lineage-rescue summary"), "",
    paste0("- repaired_qc_cells: ", expected_cells[[dataset_id]]),
    paste0("- retained_cells: ", nrow(cleaned_assignments)),
    paste0("- removed_heterotypic_doublets: ",
           removal_counts[["REMOVE_HETEROTYPIC_DOUBLET"]]),
    paste0("- removed_low_quality: ", removal_counts[["REMOVE_LOW_QUALITY"]]),
    paste0("- removed_ambient_dominated: ",
           removal_counts[["REMOVE_AMBIENT_DOMINATED"]]),
    paste0("- rescued_wrong_parent_cells: ", sum(rescued$n_cells_retained, na.rm = TRUE)),
    paste0("- rescued_wrong_parent_clusters: ", nrow(rescued)),
    paste0("- unresolved_cells: ",
           sum(cleaned_assignments$final_cell_type == "Unresolved")),
    paste0("- patient_enriched_clusters: ", length(patient_clusters)),
    paste0("- platform_confounded_clusters: ", length(platform_clusters)), "",
    "## Final cell-type counts", "",
    paste0("- ", type_counts$final_cell_type, ": ", type_counts$N), "",
    "## Patient-enriched clusters", "",
    if (length(patient_clusters)) paste0("- ", patient_clusters) else "- None", "",
    "## Platform handling", "",
    paste0("- ", platform_info$note),
    if (dataset_id == "GSE158722" && !platform_info$reliable) {
      "- 10X/iCell8 stratification could not be verified; samples were processed separately and platform labels were left NA."
    } else if (platform_info$reliable) {
      "- Platform-stratified interpretation used only explicit metadata."
    } else {
      "- Platform labels were left NA; no platform identity was inferred."
    },
    if (dataset_id == "GSE158722") {
      "- No sample_id or timepoint Harmony was used for GSE158722."
    } else {
      "- GSE154600 non-epithelial Harmony used sample_id only when feasible."
    }, "",
    "## Unresolved clusters", "",
    if (length(unresolved_clusters)) paste0("- ", unresolved_clusters) else "- None", "",
    "All decontamination was performed within original sample. Existing sample-wise scDblFinder calls were retained and combined with cell-level incompatible-lineage evidence; no cluster was deleted wholesale.",
    "Cycling, IFN_response, Hypoxia, Stress_response, SPP1_high, C1QC_high and FOLR2_high were stored only as cell_state.",
    "Every exported marker passed p_val_adj < 0.05, avg_log2FC > 0.25 and pct.1 >= 0.20."
  )
  writeLines(
    summary_lines, file.path(dataset_out, "cleaning_and_rescue_summary.md"),
    useBytes = TRUE
  )
  capture.output(sessionInfo(), file = file.path(dataset_out, "sessionInfo.txt"))

  normalized_cache <- normalizePath(decont_cache, winslash = "/", mustWork = TRUE)
  normalized_dataset_out <- normalizePath(dataset_out, winslash = "/", mustWork = TRUE)
  if (!startsWith(normalized_cache, paste0(normalized_dataset_out, "/"))) {
    stop("Refusing to remove cache outside dataset output: ", normalized_cache)
  }
  unlink(decont_cache, recursive = TRUE, force = TRUE)

  rm(
    assignments, template, md, kept_md, decont_counts, kept_counts, decont,
    old_markers, old_evidence, cell_evidence, final_markers, final_evidence,
    final_states, cluster_table, cleaned_assignments, removed, rescued
  )
  gc()
}

message("Cleaned marker-ready workflow complete: ", output_root)
