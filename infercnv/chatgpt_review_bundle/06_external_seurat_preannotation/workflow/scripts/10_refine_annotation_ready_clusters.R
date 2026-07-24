#!/usr/bin/env Rscript

# Refine marker-ready lineages and cover every repaired-QC cell.
# This step consumes step-09 outputs without modifying earlier diagnostics.

options(stringsAsFactors = FALSE, warn = 1)

required <- c(
  "yaml", "data.table", "Matrix", "Seurat", "SeuratObject",
  "SingleCellExperiment", "SummarizedExperiment", "ggplot2", "RANN"
)
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
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
seed <- as.integer(cfg$project$random_seed)
set.seed(seed)

dims_use <- as.integer(arg_value("--dims", "30"))
nfeatures <- as.integer(arg_value("--variable-features", "3000"))
global_resolution <- as.numeric(arg_value("--global-resolution", "0.8"))
family_resolution <- as.numeric(arg_value("--resolution", "0.6"))
knn_k <- as.integer(arg_value("--knn-k", "30"))
prediction_score_min <- as.numeric(arg_value("--prediction-score", "0.70"))
prediction_margin_min <- as.numeric(arg_value("--prediction-margin", "0.20"))
cluster_fraction_min <- as.numeric(arg_value("--cluster-fraction", "0.80"))
cluster_support_min <- as.integer(arg_value("--cluster-support", "30"))

input_root <- file.path(data_root, "diagnostics_v2_marker_ready")
output_root <- file.path(data_root, "diagnostics_v2_marker_ready_refined")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

expected_cells <- c(GSE154600 = 31103L, GSE158722 = 68568L)

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, na = "NA")
}

write_csv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, na = "NA", compress = "gzip")
}

save_png <- function(plot, path, width = 11, height = 8) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot = plot, width = width, height = height, dpi = 200,
         limitsize = FALSE)
}

marker_panels <- list(
  Epithelial = c(
    "EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT17", "KRT13",
    "KRT5", "KRT6A", "KRT6B", "MSLN", "WFDC2", "MUC1", "MUC16",
    "TACSTD2", "PAX8", "CLDN3", "CLDN4", "KRTCAP3"
  ),
  T_cell = c(
    "CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2", "CD2", "CD247",
    "LCK", "ITK", "BCL11B", "GIMAP7", "IL7R", "LTB"
  ),
  NK_cell = c(
    "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRC1", "PRF1", "CTSW",
    "XCL1", "XCL2", "FGFBP2", "FCGR3A"
  ),
  B_cell = c(
    "MS4A1", "CD79A", "CD79B", "CD19", "CD22", "CD37", "CD74",
    "VPREB3", "BANK1", "HLA-DRA", "CD83", "TNFRSF13C"
  ),
  Plasma_cell = c(
    "MZB1", "JCHAIN", "SDC1", "DERL3", "XBP1", "SSR4", "FKBP11",
    "PRDX4", "IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGKC"
  ),
  Myeloid = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "CTSD",
    "MS4A7", "LGALS3", "C1QA", "C1QB", "C1QC", "FCGR1A", "CD14",
    "APOE", "TREM2", "SPP1"
  ),
  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "BATF3",
    "LILRA4", "IL3RA", "GZMB", "CLEC4C", "CST3"
  ),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HPGDS", "HDC"),
  Fibroblast = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "C7", "COL6A1",
    "COL6A2", "PDGFRA", "FAP", "POSTN", "CTHRC1"
  ),
  Pericyte = c(
    "RGS5", "MCAM", "CSPG4", "PDGFRB", "ACTA2", "MYH11",
    "NOTCH3", "KCNJ8", "ABCC9", "RBP1"
  ),
  Endothelial = c(
    "PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "RAMP3",
    "PLVAP", "CLDN5", "ACKR1", "CA4", "RGCC", "ESAM"
  ),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "AHSP", "CA1")
)

cycling_markers <- c(
  "MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "STMN1",
  "TUBA1B", "HMGB2", "PCNA", "CDC20", "CDK1", "CCNB1", "CCNB2"
)

collapse_family <- function(x) {
  y <- as.character(x)
  out <- rep(NA_character_, length(y))
  out[y == "Epithelial"] <- "Epithelial"
  out[y %in% c("T_cell", "NK_cell")] <- "T_NK"
  out[y %in% c("B_cell", "Plasma_cell")] <- "B_Plasma"
  out[y %in% c("Myeloid", "Dendritic", "Mast")] <- "Myeloid"
  out[y == "Fibroblast"] <- "Fibroblast"
  out[y == "Pericyte"] <- "Pericyte"
  out[y == "Endothelial"] <- "Endothelial"
  out[y == "Erythroid"] <- "Erythroid"
  out
}

gene_upper <- function(x) toupper(as.character(x))

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
  if (!"avg_log2FC" %in% names(markers)) {
    if ("avg_logFC" %in% names(markers)) {
      data.table::setnames(markers, "avg_logFC", "avg_log2FC")
    } else {
      stop("Marker table has neither avg_log2FC nor avg_logFC")
    }
  }
  required_cols <- c("cluster", "gene", "avg_log2FC", "pct.1", "p_val_adj")
  absent <- setdiff(required_cols, names(markers))
  if (length(absent)) {
    stop("Marker table missing columns: ", paste(absent, collapse = ", "))
  }
  markers[, gene_upper := gene_upper(gene)]
  markers
}

run_significant_markers <- function(obj) {
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  DefaultAssay(obj) <- "RNA"
  if (length(unique(as.character(Idents(obj)))) < 2L) {
    return(empty_marker_table())
  }
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

score_cluster_marker_evidence <- function(markers, clusters) {
  markers <- data.table::as.data.table(markers)
  if (!nrow(markers)) markers <- empty_marker_table()
  if (!"gene_upper" %in% names(markers)) markers[, gene_upper := gene_upper(gene)]
  rows <- list()
  for (cluster_id in clusters) {
    cm <- markers[as.character(cluster) == as.character(cluster_id)]
    for (cell_type in names(marker_panels)) {
      hit <- cm[gene_upper %in% gene_upper(marker_panels[[cell_type]])]
      rows[[length(rows) + 1L]] <- data.frame(
        cluster = as.character(cluster_id),
        candidate_cell_type = cell_type,
        support_n = data.table::uniqueN(hit$gene_upper),
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
  decisions <- scores[order(cluster, -support_n, -support_weight), {
    top <- .SD[1L]
    second <- if (.N >= 2L) .SD[2L] else NULL
    suggested <- as.character(top$candidate_cell_type)
    evidence_status <- "SUPPORTED"
    strength <- if (top$support_n >= 4L) "high" else
      if (top$support_n >= 2L) "moderate" else "low"
    if (top$support_n < 2L) {
      suggested <- "Ambiguous"
      evidence_status <- "INSUFFICIENT_CANONICAL_SUPPORT"
    } else if (!is.null(second) && second$support_n >= 2L &&
               second$support_n >= top$support_n - 1L &&
               second$support_weight >= 0.75 * max(top$support_weight, 1e-8)) {
      if (collapse_family(top$candidate_cell_type) !=
          collapse_family(second$candidate_cell_type)) {
        suggested <- "Mixed_or_doublet"
        evidence_status <- "CONFLICTING_LINEAGE_EVIDENCE"
      } else {
        evidence_status <- "RELATED_LINEAGE_COMPETITION_REVIEW"
      }
    }
    list(
      suggested_broad_cell_type = suggested,
      canonical_evidence_strength = strength,
      evidence_status = evidence_status,
      canonical_markers_present = as.character(top$supporting_markers),
      second_candidate = if (is.null(second)) NA_character_ else
        as.character(second$candidate_cell_type),
      second_support_n = if (is.null(second)) NA_integer_ else
        as.integer(second$support_n)
    )
  }, by = cluster]
  list(scores = scores, decisions = decisions)
}

cycling_state_from_markers <- function(markers, clusters) {
  markers <- data.table::as.data.table(markers)
  if (!nrow(markers)) markers <- empty_marker_table()
  if (!"gene_upper" %in% names(markers)) markers[, gene_upper := gene_upper(gene)]
  panel <- gene_upper(cycling_markers)
  data.table::rbindlist(lapply(clusters, function(cluster_id) {
    hit <- markers[
      as.character(cluster) == as.character(cluster_id) & gene_upper %in% panel
    ]
    n_hit <- data.table::uniqueN(hit$gene_upper)
    data.frame(
      cluster = as.character(cluster_id),
      cycling_markers_present = if (nrow(hit)) {
        paste(unique(hit[order(-avg_log2FC)]$gene), collapse = ";")
      } else "",
      cycling_state = if (n_hit >= 3L) "Cycling_high" else
        if (n_hit >= 1L) "Cycling_possible" else "Non_cycling",
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)
}

prepare_object <- function(counts, metadata, project) {
  rownames(metadata) <- metadata$cell_id
  obj <- CreateSeuratObject(
    counts = counts, project = project, min.cells = 0, min.features = 0,
    meta.data = metadata
  )
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    obj, selection.method = "vst", nfeatures = min(nfeatures, nrow(obj)),
    verbose = FALSE
  )
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  npcs <- min(dims_use, ncol(obj) - 1L, length(VariableFeatures(obj)) - 1L)
  if (npcs < 2L) stop(project, ": insufficient cells/features for PCA")
  RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
}

cluster_object <- function(obj, reduction, strategy, resolution_value) {
  dims <- seq_len(min(dims_use, ncol(Embeddings(obj, reduction))))
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
  obj <- FindClusters(
    obj, resolution = resolution_value, random.seed = seed, verbose = FALSE
  )
  obj <- RunUMAP(
    obj, reduction = reduction, dims = dims, seed.use = seed,
    reduction.name = paste0("umap.", strategy), verbose = FALSE
  )
  obj
}

run_family_clustering <- function(obj, dataset_id, broad_family) {
  strategy <- "A_uncorrected"
  reduction <- "pca"
  can_correct <- !broad_family %in% c("Epithelial", "Unresolved_review") &&
    requireNamespace("harmony", quietly = TRUE)
  if (can_correct) {
    batch_var <- if (dataset_id == "GSE158722") "patient_id" else "sample_id"
    values <- as.character(obj@meta.data[[batch_var]])
    valid <- !is.na(values) & nzchar(values)
    tab <- table(values[valid])
    can_harmony <- all(valid) && length(tab) >= 2L && min(tab) >= 20L
    if (can_harmony) {
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
  obj <- cluster_object(obj, reduction, strategy, family_resolution)
  list(obj = obj, strategy = strategy)
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
  message("Reading prepared input: ", path)
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
  if (anyNA(idx)) {
    stop(dataset_id, ": repaired-QC cells missing from prepared input: ", sum(is.na(idx)))
  }
  selected_md <- md[idx, , drop = FALSE]
  selected_md$cell_id <- qc$cell_id
  qidx <- match(selected_md$cell_id, qc$cell_id)
  for (column in setdiff(names(qc), names(selected_md))) {
    selected_md[[column]] <- qc[[column]][qidx]
  }
  if (!"sample_id" %in% names(selected_md)) selected_md$sample_id <- qc$sample_id[qidx]
  if (!"patient_id" %in% names(selected_md)) selected_md$patient_id <- selected_md$sample_id
  if (!"timepoint" %in% names(selected_md)) selected_md$timepoint <- NA_character_
  cidx <- match(selected_md$cell_id, colnames(counts))
  counts <- counts[, cidx, drop = FALSE]
  rownames(selected_md) <- selected_md$cell_id
  rm(input, md, qc)
  gc()
  list(counts = counts, metadata = selected_md)
}

select_reference_clusters <- function(dataset_id) {
  template_path <- file.path(
    input_root, dataset_id, "annotation_ready_cluster_template.csv"
  )
  assignment_path <- file.path(
    input_root, dataset_id, "annotation_ready_cell_assignments_sampled.csv.gz"
  )
  if (!file.exists(template_path) || !file.exists(assignment_path)) {
    stop(dataset_id, ": step-09 marker-ready inputs are missing")
  }
  template <- data.table::fread(template_path)
  template[, reference_eligible :=
    lineage_conflict_flag == FALSE &
      !suggested_broad_cell_type %in% c("Ambiguous", "Mixed_or_doublet") &
      canonical_evidence_strength %in% c("moderate", "high") &
      !evidence_status %in% c(
        "INSUFFICIENT_CANONICAL_SUPPORT", "CONFLICTING_LINEAGE_EVIDENCE"
      ) & n_significant_markers >= 10L]
  template[, refined_broad_family := collapse_family(suggested_broad_cell_type)]
  template[is.na(refined_broad_family), reference_eligible := FALSE]
  template[, exclusion_reason := data.table::fcase(
    lineage_conflict_flag == TRUE, "lineage_conflict",
    suggested_broad_cell_type %in% c("Ambiguous", "Mixed_or_doublet"),
    "ambiguous_or_mixed",
    !canonical_evidence_strength %in% c("moderate", "high"),
    "low_canonical_evidence",
    evidence_status %in% c(
      "INSUFFICIENT_CANONICAL_SUPPORT", "CONFLICTING_LINEAGE_EVIDENCE"
    ), "excluded_evidence_status",
    n_significant_markers < 10L, "fewer_than_10_significant_markers",
    is.na(refined_broad_family), "unsupported_broad_type",
    default = ""
  )]
  selection <- template[, .(
    dataset_id,
    old_parent_broad_type = parent_broad_type,
    old_cluster = cluster,
    suggested_broad_cell_type,
    refined_broad_family,
    n_cells,
    reference_eligible,
    exclusion_reason
  )]
  assignments <- data.table::fread(assignment_path)
  refs <- merge(
    assignments[, .(cell_id, final_cluster)],
    selection[reference_eligible == TRUE,
              .(old_cluster, refined_broad_family)],
    by.x = "final_cluster", by.y = "old_cluster", all = FALSE
  )
  if (!nrow(refs)) stop(dataset_id, ": no reference-eligible cells")
  list(selection = selection, reference_cells = refs)
}

weighted_knn_predict <- function(embedding, reference_cells, k = 30L) {
  ref_idx <- match(reference_cells$cell_id, rownames(embedding))
  if (anyNA(ref_idx)) stop("Reference cells missing from full PCA: ", sum(is.na(ref_idx)))
  ref_labels <- as.character(reference_cells$refined_broad_family)
  k_use <- min(k, length(ref_idx))
  nn <- RANN::nn2(
    data = embedding[ref_idx, , drop = FALSE],
    query = embedding,
    k = k_use,
    searchtype = "standard"
  )
  families <- sort(unique(ref_labels))
  score <- matrix(0, nrow = nrow(embedding), ncol = length(families),
                  dimnames = list(rownames(embedding), families))
  for (j in seq_len(k_use)) {
    label <- ref_labels[nn$nn.idx[, j]]
    weight <- 1 / pmax(nn$nn.dists[, j], 1e-8)
    for (family in families) {
      hit <- label == family
      score[hit, family] <- score[hit, family] + weight[hit]
    }
  }
  score <- score / pmax(rowSums(score), 1e-12)
  ord <- t(apply(score, 1L, order, decreasing = TRUE))
  top_idx <- ord[, 1L]
  second_idx <- if (ncol(score) >= 2L) ord[, 2L] else top_idx
  top_score <- score[cbind(seq_len(nrow(score)), top_idx)]
  second_score <- if (ncol(score) >= 2L) {
    score[cbind(seq_len(nrow(score)), second_idx)]
  } else rep(0, nrow(score))
  predicted <- colnames(score)[top_idx]
  margin <- top_score - second_score
  status <- ifelse(
    top_score >= prediction_score_min & margin >= prediction_margin_min,
    "Confident", "Unresolved"
  )
  data.table::data.table(
    cell_id = rownames(embedding),
    predicted_broad_family = ifelse(status == "Confident", predicted, "Unresolved"),
    prediction_score = top_score,
    prediction_margin = margin,
    prediction_status = status
  )
}

global_cluster_gate <- function(predictions) {
  audit <- predictions[, {
    confident <- .SD[prediction_status == "Confident"]
    counts <- sort(table(confident$predicted_broad_family), decreasing = TRUE)
    dominant <- if (length(counts)) names(counts)[[1L]] else NA_character_
    second <- if (length(counts) >= 2L) names(counts)[[2L]] else NA_character_
    fraction <- if (length(counts)) as.numeric(counts[[1L]] / sum(counts)) else NA_real_
    second_fraction <- if (length(counts) >= 2L) {
      as.numeric(counts[[2L]] / sum(counts))
    } else 0
    refined <- if (length(counts) && counts[[1L]] >= cluster_support_min &&
                   fraction >= cluster_fraction_min) dominant else "Unresolved_review"
    list(
      n_cells = .N,
      n_confident_cells = nrow(confident),
      dominant_family = dominant,
      dominant_family_fraction = fraction,
      second_family = second,
      second_family_fraction = second_fraction,
      refined_broad_family = refined,
      cluster_assignment_status = if (refined == "Unresolved_review") {
        "UNRESOLVED_CLUSTER_REVIEW"
      } else "ASSIGNED_80PCT_CONSISTENT"
    )
  }, by = global_cluster]
  out <- merge(predictions, audit[, .(global_cluster, cluster_family = refined_broad_family)],
               by = "global_cluster", all.x = TRUE)
  out[, refined_broad_family := data.table::fcase(
    cluster_family == "Unresolved_review", "Unresolved_review",
    prediction_status == "Unresolved", cluster_family,
    predicted_broad_family == cluster_family, cluster_family,
    default = "Unresolved_review"
  )]
  out[, cluster_family := NULL]
  list(audit = audit, cells = out)
}

for (dataset_id in datasets) {
  message("===== ", dataset_id, " =====")
  dataset_out <- file.path(output_root, dataset_id)
  if (dir.exists(dataset_out)) {
    if (!force) stop(dataset_id, ": output exists; use --force")
    unlink(dataset_out, recursive = TRUE, force = TRUE)
  }
  dir.create(dataset_out, recursive = TRUE, showWarnings = FALSE)

  refs <- select_reference_clusters(dataset_id)
  write_csv(refs$selection, file.path(dataset_out, "01_reference_cluster_selection.csv"))

  input <- read_complete_input(dataset_id)
  if (ncol(input$counts) != expected_cells[[dataset_id]]) {
    stop(dataset_id, ": expected ", expected_cells[[dataset_id]],
         " repaired-QC cells, found ", ncol(input$counts))
  }
  global_obj <- prepare_object(
    input$counts, input$metadata, paste0(dataset_id, "_refined_global")
  )
  global_obj <- cluster_object(
    global_obj, "pca", "global_A_uncorrected", global_resolution
  )
  global_cluster <- as.character(Idents(global_obj))
  names(global_cluster) <- colnames(global_obj)

  embedding <- Embeddings(global_obj, "pca")[, seq_len(dims_use), drop = FALSE]
  predictions <- weighted_knn_predict(
    embedding, refs$reference_cells, knn_k
  )
  midx <- match(predictions$cell_id, rownames(global_obj@meta.data))
  predictions[, `:=`(
    dataset_id = dataset_id,
    sample_id = as.character(global_obj@meta.data$sample_id[midx]),
    patient_id = as.character(global_obj@meta.data$patient_id[midx]),
    timepoint = as.character(global_obj@meta.data$timepoint[midx]),
    global_cluster = global_cluster[cell_id]
  )]
  data.table::setcolorder(predictions, c(
    "dataset_id", "cell_id", "sample_id", "patient_id", "timepoint",
    "global_cluster", "predicted_broad_family", "prediction_score",
    "prediction_margin", "prediction_status"
  ))
  prediction_export_cols <- setdiff(names(predictions), "global_cluster")
  write_csv_gz(
    predictions[, ..prediction_export_cols],
    file.path(dataset_out, "02_full_cell_broad_family_predictions.csv.gz")
  )

  gate <- global_cluster_gate(predictions)
  write_csv(gate$audit, file.path(dataset_out, "03_global_cluster_family_audit.csv"))
  refined_cells <- gate$cells
  write_csv_gz(
    refined_cells[, .(
      dataset_id, cell_id, sample_id, patient_id, timepoint, global_cluster,
      predicted_broad_family, prediction_score, prediction_margin,
      prediction_status, refined_broad_family
    )],
    file.path(dataset_out, "04_full_cell_refined_broad_family.csv.gz")
  )

  # Keep only the global coordinates needed for plotting. Releasing the full
  # PCA/scale object prevents a second full expression object from coexisting
  # with large family objects (notably GSE158722 epithelial cells).
  global_cell_ids <- colnames(global_obj)
  global_umap <- Embeddings(global_obj, "umap.global_A_uncorrected")
  rm(global_obj, embedding, input)
  gc()

  final_marker_parts <- list()
  final_template_parts <- list()
  final_assignment_parts <- list()
  final_evidence_parts <- list()

  for (broad_family in sort(unique(refined_cells$refined_broad_family))) {
    cells <- refined_cells[refined_broad_family == broad_family, cell_id]
    if (!length(cells)) next
    message(dataset_id, " / ", broad_family, ": ", length(cells), " cells")
    family_input <- read_complete_input(dataset_id)
    cidx <- match(cells, colnames(family_input$counts))
    sub_counts <- family_input$counts[, cidx, drop = FALSE]
    sub_md <- family_input$metadata[
      match(cells, family_input$metadata$cell_id), , drop = FALSE
    ]
    rm(family_input)
    gc()
    sub_obj <- prepare_object(
      sub_counts, sub_md, paste0(dataset_id, "_refined_", broad_family)
    )
    rm(sub_counts, sub_md)
    gc()
    result <- run_family_clustering(sub_obj, dataset_id, broad_family)
    sub_obj <- result$obj
    strategy <- result$strategy
    raw_cluster <- as.character(Idents(sub_obj))
    final_cluster <- paste(broad_family, strategy, raw_cluster, sep = "__")
    Idents(sub_obj) <- factor(final_cluster)

    markers <- run_significant_markers(sub_obj)
    clusters <- sort(unique(final_cluster))
    evidence <- score_cluster_marker_evidence(markers, clusters)
    cycle <- cycling_state_from_markers(markers, clusters)
    decision <- merge(evidence$decisions, cycle, by = "cluster", all = TRUE)

    cell_stats <- data.table::data.table(
      cluster = final_cluster,
      sample_id = as.character(sub_obj$sample_id),
      patient_id = as.character(sub_obj$patient_id)
    )[, {
      sample_tab <- sort(table(sample_id), decreasing = TRUE)
      list(
        n_cells = .N,
        n_samples = data.table::uniqueN(sample_id),
        n_patients = data.table::uniqueN(patient_id),
        dominant_sample = names(sample_tab)[[1L]],
        dominant_sample_fraction = as.numeric(sample_tab[[1L]] / .N)
      )
    }, by = cluster]
    decision <- merge(decision, cell_stats, by = "cluster", all.x = TRUE)

    if (nrow(markers)) {
      markers[, `:=`(
        dataset_id = dataset_id,
        parent_broad_type = broad_family,
        clustering_strategy = strategy
      )]
      top_text <- markers[order(cluster, -avg_log2FC, p_val_adj), .(
        top20_significant_markers = paste(head(gene, 20L), collapse = ";"),
        top50_significant_markers = paste(head(gene, 50L), collapse = ";"),
        n_significant_markers = .N
      ), by = cluster]
    } else {
      top_text <- data.table::data.table(
        cluster = clusters, top20_significant_markers = "",
        top50_significant_markers = "", n_significant_markers = 0L
      )
    }
    decision <- merge(decision, top_text, by = "cluster", all.x = TRUE)
    decision[, `:=`(
      dataset_id = dataset_id,
      parent_broad_type = broad_family,
      clustering_strategy = strategy
    )]
    decision[, lineage_conflict_flag :=
      suggested_broad_cell_type == "Mixed_or_doublet" |
      (
        !suggested_broad_cell_type %in% c("Ambiguous", "Mixed_or_doublet") &
          parent_broad_type != "Unresolved_review" &
          collapse_family(suggested_broad_cell_type) != parent_broad_type &
          canonical_evidence_strength %in% c("moderate", "high")
      )]
    decision[, marker_status := data.table::fcase(
      parent_broad_type == "Unresolved_review", "REVIEW_UNRESOLVED_LINEAGE",
      suggested_broad_cell_type == "Mixed_or_doublet" |
        evidence_status == "CONFLICTING_LINEAGE_EVIDENCE",
      "REVIEW_MIXED_OR_DOUBLET",
      lineage_conflict_flag, "REVIEW_LINEAGE_CONFLICT",
      suggested_broad_cell_type == "Ambiguous" |
        evidence_status == "INSUFFICIENT_CANONICAL_SUPPORT" |
        canonical_evidence_strength == "low",
      "REVIEW_AMBIGUOUS_LINEAGE",
      n_significant_markers >= 10L, "READY_FOR_ANNOTATION",
      n_significant_markers > 0L, "FEW_SIGNIFICANT_MARKERS_REVIEW",
      default = "NO_SIGNIFICANT_MARKERS"
    )]
    decision[, `:=`(
      manual_cell_type = "",
      manual_cell_subtype = "",
      manual_confidence = "",
      manual_notes = ""
    )]

    final_marker_parts[[broad_family]] <- markers
    final_template_parts[[broad_family]] <- decision
    evidence$scores[, `:=`(
      dataset_id = dataset_id,
      parent_broad_type = broad_family,
      clustering_strategy = strategy
    )]
    final_evidence_parts[[broad_family]] <- evidence$scores
    final_assignment_parts[[broad_family]] <- data.table::data.table(
      dataset_id = dataset_id,
      cell_id = colnames(sub_obj),
      sample_id = as.character(sub_obj$sample_id),
      patient_id = as.character(sub_obj$patient_id),
      timepoint = as.character(sub_obj$timepoint),
      parent_broad_type = broad_family,
      clustering_strategy = strategy,
      final_cluster = final_cluster
    )
    rm(sub_obj, result, markers, evidence, cycle, decision)
    gc()
  }

  final_markers <- data.table::rbindlist(final_marker_parts, fill = TRUE)
  final_template <- data.table::rbindlist(final_template_parts, fill = TRUE)
  final_assignments <- data.table::rbindlist(final_assignment_parts, fill = TRUE)
  final_evidence <- data.table::rbindlist(final_evidence_parts, fill = TRUE)

  if (nrow(final_assignments) != expected_cells[[dataset_id]] ||
      data.table::uniqueN(final_assignments$cell_id) != expected_cells[[dataset_id]]) {
    stop(dataset_id, ": final assignments do not cover every repaired-QC cell")
  }
  if (nrow(final_markers) && any(
    !is.finite(final_markers$p_val_adj) | final_markers$p_val_adj >= 0.05 |
      final_markers$avg_log2FC <= 0.25 | final_markers$pct.1 < 0.20
  )) stop(dataset_id, ": nonsignificant marker leaked into refined output")

  preferred_cols <- c(
    "dataset_id", "parent_broad_type", "clustering_strategy", "cluster",
    "n_cells", "n_samples", "n_patients", "dominant_sample",
    "dominant_sample_fraction", "suggested_broad_cell_type",
    "canonical_evidence_strength", "canonical_markers_present",
    "second_candidate", "cycling_state", "cycling_markers_present",
    "n_significant_markers", "marker_status", "lineage_conflict_flag",
    "top20_significant_markers", "top50_significant_markers",
    "manual_cell_type", "manual_cell_subtype", "manual_confidence", "manual_notes"
  )
  data.table::setcolorder(
    final_template,
    c(intersect(preferred_cols, names(final_template)),
      setdiff(names(final_template), preferred_cols))
  )
  write_csv(
    final_template,
    file.path(dataset_out, "annotation_ready_cluster_template_refined.csv")
  )
  write_csv_gz(
    final_assignments,
    file.path(dataset_out, "annotation_ready_full_cell_assignments.csv.gz")
  )
  write_csv_gz(
    final_markers,
    file.path(dataset_out, "annotation_ready_significant_markers_refined.csv.gz")
  )
  write_csv(
    final_evidence,
    file.path(dataset_out, "annotation_ready_canonical_evidence_refined.csv")
  )

  refined_map <- refined_cells[
    match(global_cell_ids, cell_id), refined_broad_family
  ]
  final_map <- final_assignments[match(global_cell_ids, cell_id), final_cluster]
  plot_df <- data.frame(
    UMAP_1 = global_umap[, 1L], UMAP_2 = global_umap[, 2L],
    broad_family = refined_map, final_cluster = final_map
  )
  p_family <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = broad_family)) +
    geom_point(size = 0.12, alpha = 0.65) + theme_bw() +
    labs(title = paste(dataset_id, "refined broad families"), color = "Family")
  save_png(
    p_family, file.path(dataset_out, "UMAP_refined_broad_families.png"), 10, 7
  )
  p_cluster <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = final_cluster)) +
    geom_point(size = 0.10, alpha = 0.60) + theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste(dataset_id, "refined annotation clusters"))
  save_png(p_cluster, file.path(dataset_out, "UMAP_refined_clusters.png"), 10, 7)

  if (data.table::uniqueN(final_map) >= 2L) {
    dot_input <- read_complete_input(dataset_id)
    dot_md <- dot_input$metadata
    rownames(dot_md) <- dot_md$cell_id
    dot_obj <- CreateSeuratObject(
      counts = dot_input$counts, project = paste0(dataset_id, "_refined_dotplot"),
      min.cells = 0, min.features = 0, meta.data = dot_md
    )
    rm(dot_input)
    gc()
    dot_obj <- NormalizeData(dot_obj, verbose = FALSE)
    Idents(dot_obj) <- factor(final_map)
    panel <- unique(c(unlist(marker_panels, use.names = FALSE), cycling_markers))
    panel <- panel[panel %in% rownames(dot_obj)]
    p_dot <- DotPlot(dot_obj, features = panel) + RotatedAxis() +
      labs(
        title = paste(dataset_id, "canonical markers by refined cluster"),
        x = "Canonical marker", y = "Refined cluster"
      )
    save_png(
      p_dot, file.path(dataset_out, "canonical_marker_dotplot_refined.png"),
      18, 12
    )
    rm(dot_obj, dot_md)
    gc()
  }

  repaired_qc_cells <- expected_cells[[dataset_id]]
  reference_cells <- nrow(refs$reference_cells)
  confidently_predicted_cells <- sum(predictions$prediction_status == "Confident")
  unresolved_prediction_cells <- sum(predictions$prediction_status == "Unresolved")
  unresolved_review_cells <- sum(
    final_assignments$parent_broad_type == "Unresolved_review"
  )
  summary_lines <- c(
    paste0("# ", dataset_id, " refined marker-ready annotation summary"), "",
    paste0("- repaired_qc_cells: ", repaired_qc_cells),
    paste0("- reference_cells: ", reference_cells),
    paste0("- confidently_predicted_cells: ", confidently_predicted_cells),
    paste0("- unresolved_prediction_cells: ", unresolved_prediction_cells),
    paste0("- final_assigned_cells: ", nrow(final_assignments)),
    paste0("- unresolved_review_cells: ", unresolved_review_cells),
    paste0("- full_assignment_coverage: ",
           sprintf("%.6f", nrow(final_assignments) / repaired_qc_cells)),
    paste0("- number_of_final_clusters: ", nrow(final_template)),
    paste0("- ready_for_annotation_clusters: ",
           sum(final_template$marker_status == "READY_FOR_ANNOTATION")),
    paste0("- review_lineage_conflict_clusters: ",
           sum(final_template$marker_status == "REVIEW_LINEAGE_CONFLICT")),
    paste0("- review_ambiguous_clusters: ",
           sum(final_template$marker_status == "REVIEW_AMBIGUOUS_LINEAGE")),
    paste0("- review_mixed_or_doublet_clusters: ",
           sum(final_template$marker_status == "REVIEW_MIXED_OR_DOUBLET")),
    paste0("- few_marker_clusters: ",
           sum(final_template$marker_status == "FEW_SIGNIFICANT_MARKERS_REVIEW")),
    paste0("- no_marker_clusters: ",
           sum(final_template$marker_status == "NO_SIGNIFICANT_MARKERS")), "",
    "Epithelial remained uncorrected.",
    if (dataset_id == "GSE158722") {
      "GSE158722 used no sample/timepoint Harmony."
    } else {
      "GSE154600 non-epithelial Harmony used sample_id only when feasible."
    },
    "Cycling was treated only as a state.",
    "Every exported marker passed the significance filter."
  )
  writeLines(
    summary_lines,
    file.path(dataset_out, "annotation_ready_summary_refined.md"),
    useBytes = TRUE
  )
  capture.output(sessionInfo(), file = file.path(dataset_out, "sessionInfo.txt"))
  rm(
    global_cell_ids, global_umap, predictions, gate, refined_cells,
    final_markers, final_template, final_assignments, final_evidence
  )
  gc()
}

message("Refined marker-ready workflow complete: ", output_root)
