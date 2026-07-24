#!/usr/bin/env Rscript

# Build annotation-ready marker tables from repaired-QC cells.
# Minimal scope:
#   1) balanced global discovery subset;
#   2) cluster-level broad cell-type assignment from significant RNA markers;
#   3) Cycling is a state, never an exclusive lineage;
#   4) epithelial cells remain uncorrected;
#   5) non-epithelial cells use Harmony by sample (GSE154600/GSE147082)
#      or by patient (GSE158722);
#   6) only adjusted-P significant markers are exported for annotation.
#
# Large outputs remain under DataRoot and must not be committed.

options(stringsAsFactors = FALSE, warn = 1)

required <- c(
  "yaml", "data.table", "Matrix", "Seurat", "SingleCellExperiment",
  "SummarizedExperiment", "ggplot2"
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
assert_datasets(datasets, c("GSE147082", "GSE154600", "GSE158722"))

force <- has_flag("--force")
seed <- as.integer(cfg$project$random_seed)
set.seed(seed)

max_global_cells <- as.integer(
  arg_value("--max-global-cells",
            as.character(cfg$analysis$max_cells_per_lineage_strategy %||% 30000L))
)
resolution <- as.numeric(
  arg_value("--resolution",
            as.character(cfg$analysis$primary_resolution %||% 0.6))
)
dims_use <- as.integer(
  arg_value("--dims", as.character(cfg$analysis$dims_use %||% 30L))
)
nfeatures <- as.integer(
  arg_value("--variable-features",
            as.character(cfg$analysis$variable_features %||% 3000L))
)

output_root <- file.path(data_root, "diagnostics_v2_marker_ready")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, na = "NA")
}

write_csv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, na = "NA")
}

save_plot_pair <- function(plot, stem, width = 10, height = 7) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height,
         dpi = 200, limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height,
         limitsize = FALSE)
}

# Deliberately specific panels. These provide evidence, not irreversible labels.
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
  Mast = c(
    "TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HPGDS", "HDC"
  ),
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
  Erythroid = c(
    "HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "AHSP", "CA1"
  )
)

cycling_markers <- c(
  "MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "STMN1",
  "TUBA1B", "HMGB2", "PCNA", "CDC20", "CDK1", "CCNB1", "CCNB2"
)

cell_type_family <- function(x) {
  y <- as.character(x)
  out <- rep("other", length(y))
  out[y %in% c("Epithelial")] <- "epithelial"
  out[y %in% c("T_cell", "NK_cell", "B_cell", "Plasma_cell")] <- "lymphoid"
  out[y %in% c("Myeloid", "Dendritic", "Mast")] <- "myeloid"
  out[y %in% c("Fibroblast", "Pericyte")] <- "stromal"
  out[y %in% c("Endothelial")] <- "vascular"
  out[y %in% c("Erythroid")] <- "erythroid"
  out
}

incompatible_lineage_pair <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  cell_type_family(x) != cell_type_family(y) |
    (x %in% c("T_cell", "NK_cell") &
       y %in% c("B_cell", "Plasma_cell")) |
    (y %in% c("T_cell", "NK_cell") &
       x %in% c("B_cell", "Plasma_cell"))
}

gene_upper <- function(x) toupper(as.character(x))

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

empty_marker_table <- function() {
  data.table::data.table(
    cluster = character(), gene = character(), gene_upper = character(),
    avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(),
    p_val = numeric(), p_val_adj = numeric()
  )
}

run_significant_markers <- function(obj) {
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  DefaultAssay(obj) <- "RNA"
  if (length(unique(as.character(Idents(obj)))) < 2L) {
    return(empty_marker_table())
  }
  markers <- FindAllMarkers(
    obj,
    assay = "RNA",
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.20,
    logfc.threshold = 0.25,
    return.thresh = 1,
    verbose = FALSE
  )
  if (!nrow(markers)) return(empty_marker_table())
  markers <- normalize_marker_columns(markers)
  markers[
    is.finite(p_val_adj) &
      p_val_adj < 0.05 &
      is.finite(avg_log2FC) &
      avg_log2FC > 0.25 &
      is.finite(pct.1) &
      pct.1 >= 0.20
  ]
}

top_expressed_by_cluster <- function(obj, n = 30L) {
  DefaultAssay(obj) <- "RNA"
  avg <- tryCatch(
    AverageExpression(obj, assays = "RNA", layer = "data",
                      verbose = FALSE)$RNA,
    error = function(e) AverageExpression(obj, assays = "RNA",
                                          slot = "data",
                                          verbose = FALSE)$RNA
  )
  out <- data.table::rbindlist(lapply(seq_len(ncol(avg)), function(i) {
    values <- avg[, i]
    names(values) <- rownames(avg)
    keep <- !grepl("^(RPL|RPS|MT-|MTRNR)", toupper(names(values)))
    values <- sort(values[keep], decreasing = TRUE)
    data.frame(
      cluster = colnames(avg)[[i]],
      rank = seq_len(min(n, length(values))),
      gene = names(values)[seq_len(min(n, length(values)))],
      average_expression = unname(values[seq_len(min(n, length(values)))])
    )
  }), fill = TRUE)
  out
}

score_cluster_marker_evidence <- function(markers, clusters) {
  markers <- data.table::as.data.table(markers)
  if (!nrow(markers)) markers <- empty_marker_table()
  if (!"gene_upper" %in% names(markers)) {
    markers[, gene_upper := gene_upper(gene)]
  }

  score_rows <- list()
  for (cluster_id in clusters) {
    cluster_markers <- markers[as.character(cluster) == as.character(cluster_id)]
    for (cell_type in names(marker_panels)) {
      panel <- gene_upper(marker_panels[[cell_type]])
      hit <- cluster_markers[gene_upper %in% panel]
      score_rows[[length(score_rows) + 1L]] <- data.frame(
        cluster = as.character(cluster_id),
        candidate_cell_type = cell_type,
        support_n = data.table::uniqueN(hit$gene_upper),
        support_weight = if (nrow(hit)) {
          sum(pmax(hit$avg_log2FC, 0) * pmax(hit$pct.1, 0), na.rm = TRUE)
        } else {
          0
        },
        supporting_markers = if (nrow(hit)) {
          paste(unique(hit[order(-avg_log2FC)]$gene), collapse = ";")
        } else {
          ""
        },
        stringsAsFactors = FALSE
      )
    }
  }
  scores <- data.table::rbindlist(score_rows, fill = TRUE)

  decisions <- scores[order(cluster, -support_n, -support_weight), {
    top <- .SD[1L]
    second <- if (.N >= 2L) .SD[2L] else NULL

    suggested <- as.character(top$candidate_cell_type)
    status <- "SUPPORTED"
    strength <- if (top$support_n >= 4L) {
      "high"
    } else if (top$support_n >= 2L) {
      "moderate"
    } else {
      "low"
    }

    if (top$support_n < 2L) {
      suggested <- "Ambiguous"
      status <- "INSUFFICIENT_CANONICAL_SUPPORT"
    } else if (!is.null(second) &&
               second$support_n >= 2L &&
               second$support_n >= top$support_n - 1L &&
               second$support_weight >= 0.75 * max(top$support_weight, 1e-8)) {
      top_family <- cell_type_family(top$candidate_cell_type)
      second_family <- cell_type_family(second$candidate_cell_type)
      if (top_family != second_family) {
        suggested <- "Mixed_or_doublet"
        status <- "CONFLICTING_LINEAGE_EVIDENCE"
      } else {
        status <- "RELATED_LINEAGE_COMPETITION_REVIEW"
      }
    }

    list(
      suggested_broad_cell_type = suggested,
      suggested_family = if (suggested %in% names(marker_panels)) {
        cell_type_family(suggested)
      } else {
        NA_character_
      },
      canonical_evidence_strength = strength,
      evidence_status = status,
      top_candidate = as.character(top$candidate_cell_type),
      top_candidate_family = cell_type_family(top$candidate_cell_type),
      top_support_n = as.integer(top$support_n),
      top_support_weight = as.numeric(top$support_weight),
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
  if (!"gene_upper" %in% names(markers)) {
    markers[, gene_upper := gene_upper(gene)]
  }
  panel <- gene_upper(cycling_markers)
  data.table::rbindlist(lapply(clusters, function(cluster_id) {
    hit <- markers[
      as.character(cluster) == as.character(cluster_id) &
        gene_upper %in% panel
    ]
    n_hit <- data.table::uniqueN(hit$gene_upper)
    data.frame(
      cluster = as.character(cluster_id),
      cycling_marker_n = n_hit,
      cycling_markers_present = if (nrow(hit)) {
        paste(unique(hit[order(-avg_log2FC)]$gene), collapse = ";")
      } else {
        ""
      },
      cycling_state = if (n_hit >= 3L) {
        "Cycling_high"
      } else if (n_hit >= 1L) {
        "Cycling_possible"
      } else {
        "Non_cycling"
      },
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)
}

prepare_object <- function(counts, metadata, project, reduction_name = "pca") {
  rownames(metadata) <- metadata$cell_id
  obj <- CreateSeuratObject(
    counts = counts,
    project = project,
    min.cells = 0,
    min.features = 0,
    meta.data = metadata
  )
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = min(nfeatures, nrow(obj)),
    verbose = FALSE
  )
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  npcs <- min(max(dims_use, 20L), 50L, ncol(obj) - 1L,
              length(VariableFeatures(obj)) - 1L)
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = npcs,
    verbose = FALSE
  )
  obj
}

cluster_object <- function(obj, reduction, strategy_name) {
  dims <- seq_len(min(dims_use, ncol(Embeddings(obj, reduction))))
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
  obj <- FindClusters(
    obj,
    resolution = resolution,
    random.seed = seed,
    verbose = FALSE
  )
  obj <- RunUMAP(
    obj,
    reduction = reduction,
    dims = dims,
    seed.use = seed,
    reduction.name = paste0("umap.", strategy_name),
    verbose = FALSE
  )
  obj
}

run_lineage_clustering <- function(obj, dataset_id, broad_type) {
  strategy <- "A_uncorrected"
  reduction <- "pca"

  if (broad_type != "Epithelial" && requireNamespace("harmony", quietly = TRUE)) {
    batch_var <- if (dataset_id == "GSE158722") "patient_id" else "sample_id"
    values <- as.character(obj@meta.data[[batch_var]])
    tab <- table(values)
    can_harmony <- length(tab) >= 2L && min(tab) >= 20L

    if (can_harmony) {
      obj <- harmony::RunHarmony(
        obj,
        group.by.vars = batch_var,
        reduction.use = "pca",
        verbose = FALSE
      )
      reduction <- "harmony"
      strategy <- if (batch_var == "patient_id") {
        "B_harmony_patient"
      } else {
        "B_harmony_sample"
      }
    }
  }

  obj <- cluster_object(obj, reduction, strategy)
  list(obj = obj, strategy = strategy, reduction = reduction)
}

for (dataset_id in datasets) {
  message("===== ", dataset_id, " =====")
  dataset_out <- file.path(output_root, dataset_id)

  if (dir.exists(dataset_out)) {
    if (!force) {
      stop(
        dataset_id, ": output exists: ", dataset_out,
        ". Re-run with --force to replace only this marker-ready folder."
      )
    }
    unlink(dataset_out, recursive = TRUE, force = TRUE)
  }
  dir.create(dataset_out, recursive = TRUE, showWarnings = FALSE)

  ds_cfg <- cfg$datasets[[dataset_id]]
  prepared_path <- file.path(prepared_root, ds_cfg$prepared_sce)
  if (!file.exists(prepared_path)) {
    prepared_path <- file.path(
      data_root, dataset_id, "objects",
      paste0(dataset_id, "_preannotation.rds")
    )
  }
  if (!file.exists(prepared_path)) {
    stop(dataset_id, ": prepared input not found: ", prepared_path)
  }

  message("Reading prepared input: ", prepared_path)
  sce <- readRDS(prepared_path)
  if (inherits(sce, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(sce, "counts")
    sce_md <- as.data.frame(SummarizedExperiment::colData(sce))
  } else if (inherits(sce, "Seurat")) {
    counts <- SeuratObject::LayerData(sce, assay = "RNA", layer = "counts")
    sce_md <- as.data.frame(sce[[]], check.names = FALSE)
  } else {
    stop(dataset_id, ": unsupported prepared input class: ", class(sce)[[1L]])
  }
  sce_md$cell_id <- colnames(counts)

  qc <- read_repaired_qc(data_root, dataset_id)
  qc <- qc[diagnostics_v2_qc_pass %in% TRUE]
  n_repaired_qc_cells <- nrow(qc)
  idx <- match(qc$cell_id, sce_md$cell_id)
  if (anyNA(idx)) {
    stop(dataset_id, ": ", sum(is.na(idx)),
         " repaired-QC cells are missing from the prepared SCE")
  }

  selected_md <- sce_md[idx, , drop = FALSE]
  selected_md$cell_id <- qc$cell_id

  qidx <- match(selected_md$cell_id, qc$cell_id)
  for (column in setdiff(names(qc), names(selected_md))) {
    selected_md[[column]] <- qc[[column]][qidx]
  }

  if (!"sample_id" %in% names(selected_md)) {
    selected_md$sample_id <- qc$sample_id[qidx]
  }
  if (!"patient_id" %in% names(selected_md)) {
    selected_md$patient_id <- selected_md$sample_id
  }
  if (!"timepoint" %in% names(selected_md)) {
    selected_md$timepoint <- NA_character_
  }

  selected_md <- data.table::as.data.table(selected_md)
  chosen <- balanced_cells(
    selected_md,
    cap = max_global_cells,
    seed = seed
  )
  cidx <- match(chosen, colnames(counts))
  if (anyNA(cidx)) {
    stop(dataset_id, ": balanced cells missing from count matrix")
  }

  counts_sub <- counts[, cidx, drop = FALSE]
  md_sub <- as.data.frame(selected_md[match(chosen, cell_id)])
  rownames(md_sub) <- md_sub$cell_id

  rm(counts, sce, sce_md, qc, selected_md)
  gc()

  global_obj <- prepare_object(
    counts_sub,
    md_sub,
    project = paste0(dataset_id, "_marker_ready")
  )
  rm(counts_sub, md_sub)
  gc()

  global_obj <- cluster_object(global_obj, "pca", "global_A_uncorrected")
  global_cluster <- as.character(Idents(global_obj))
  names(global_cluster) <- colnames(global_obj)

  global_markers <- run_significant_markers(global_obj)
  global_clusters <- sort(unique(global_cluster))
  global_evidence <- score_cluster_marker_evidence(
    global_markers, global_clusters
  )
  global_cycle <- cycling_state_from_markers(
    global_markers, global_clusters
  )

  global_decisions <- merge(
    global_evidence$decisions,
    global_cycle,
    by = "cluster",
    all = TRUE
  )
  global_counts <- data.table::as.data.table(
    as.data.frame(table(cluster = global_cluster))
  )
  data.table::setnames(global_counts, "Freq", "n_cells")
  global_decisions <- merge(
    global_decisions,
    global_counts,
    by = "cluster",
    all.x = TRUE
  )
  global_decisions[, dataset_id := dataset_id]

  top_global <- if (nrow(global_markers)) {
    global_markers[
      order(cluster, -avg_log2FC, p_val_adj),
      .(
        top20_significant_markers =
          paste(head(gene, 20L), collapse = ";"),
        n_significant_markers = .N
      ),
      by = cluster
    ]
  } else {
    data.table::data.table(
      cluster = global_clusters,
      top20_significant_markers = "",
      n_significant_markers = 0L
    )
  }
  global_decisions <- merge(
    global_decisions, top_global, by = "cluster", all.x = TRUE
  )

  write_csv(
    global_decisions,
    file.path(dataset_out, "01_global_broad_cluster_review.csv")
  )
  write_csv(
    global_evidence$scores,
    file.path(dataset_out, "01_global_canonical_marker_scores.csv")
  )
  write_csv_gz(
    global_markers,
    file.path(dataset_out, "01_global_significant_markers.csv.gz")
  )

  global_meta <- data.table::data.table(
    cell_id = colnames(global_obj),
    global_cluster = global_cluster[colnames(global_obj)]
  )
  global_meta <- merge(
    global_meta,
    global_decisions[
      ,
      .(
        global_cluster = cluster,
        suggested_broad_cell_type,
        canonical_evidence_strength,
        evidence_status,
        cycling_state
      )
    ],
    by = "global_cluster",
    all.x = TRUE
  )

  midx <- match(global_meta$cell_id, rownames(global_obj@meta.data))
  global_meta[, sample_id := as.character(
    global_obj@meta.data$sample_id[midx]
  )]
  global_meta[, patient_id := as.character(
    global_obj@meta.data$patient_id[midx]
  )]
  global_meta[, timepoint := as.character(
    global_obj@meta.data$timepoint[midx]
  )]

  final_marker_parts <- list()
  final_template_parts <- list()
  final_assignment_parts <- list()
  final_evidence_parts <- list()

  resolved_types <- setdiff(
    unique(global_meta$suggested_broad_cell_type),
    c("Ambiguous", "Mixed_or_doublet", NA_character_, "")
  )

  for (broad_type in resolved_types) {
    cells <- global_meta[
      suggested_broad_cell_type == broad_type,
      cell_id
    ]
    if (length(cells) < 100L) next

    message(dataset_id, " / ", broad_type, ": ", length(cells), " cells")
    sub_obj <- subset(global_obj, cells = cells)
    sub_counts <- tryCatch(
      SeuratObject::LayerData(sub_obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(sub_obj, assay = "RNA", slot = "counts")
    )
    sub_md <- as.data.frame(sub_obj@meta.data, check.names = FALSE)
    sub_md$cell_id <- rownames(sub_md)
    sub_obj <- prepare_object(
      sub_counts,
      sub_md,
      project = paste0(dataset_id, "_", broad_type)
    )
    rm(sub_counts, sub_md)
    gc()

    result <- run_lineage_clustering(sub_obj, dataset_id, broad_type)
    sub_obj <- result$obj
    strategy <- result$strategy

    raw_cluster <- as.character(Idents(sub_obj))
    final_cluster <- paste(
      broad_type, strategy, raw_cluster, sep = "__"
    )
    Idents(sub_obj) <- factor(final_cluster)

    markers <- run_significant_markers(sub_obj)
    clusters <- sort(unique(final_cluster))
    evidence <- score_cluster_marker_evidence(markers, clusters)
    cycle <- cycling_state_from_markers(markers, clusters)
    decision <- merge(
      evidence$decisions, cycle, by = "cluster", all = TRUE
    )

    counts_table <- data.table::as.data.table(
      as.data.frame(table(cluster = final_cluster))
    )
    data.table::setnames(counts_table, "Freq", "n_cells")
    decision <- merge(decision, counts_table, by = "cluster", all.x = TRUE)

    if (nrow(markers)) {
      markers[, `:=`(
        dataset_id = dataset_id,
        parent_broad_type = broad_type,
        clustering_strategy = strategy
      )]
      top_text <- markers[
        order(cluster, -avg_log2FC, p_val_adj),
        .(
          top20_significant_markers =
            paste(head(gene, 20L), collapse = ";"),
          top50_significant_markers =
            paste(head(gene, 50L), collapse = ";"),
          n_significant_markers = .N
        ),
        by = cluster
      ]
    } else {
      top_text <- data.table::data.table(
        cluster = clusters,
        top20_significant_markers = "",
        top50_significant_markers = "",
        n_significant_markers = 0L
      )
    }

    decision <- merge(decision, top_text, by = "cluster", all.x = TRUE)
    decision[, `:=`(
      dataset_id = dataset_id,
      parent_broad_type = broad_type,
      clustering_strategy = strategy,
      marker_status = data.table::fcase(
        n_significant_markers >= 10L, "READY_FOR_ANNOTATION",
        n_significant_markers > 0L, "FEW_SIGNIFICANT_MARKERS_REVIEW",
        default = "NO_SIGNIFICANT_MARKERS"
      ),
      lineage_conflict_flag =
        suggested_broad_cell_type %in% c("Mixed_or_doublet") |
        (
          !suggested_broad_cell_type %in% c("Ambiguous") &
          !is.na(suggested_family) &
          incompatible_lineage_pair(
            suggested_broad_cell_type, broad_type
          ) &
          canonical_evidence_strength %in% c("moderate", "high")
        ),
      manual_cell_type = "",
      manual_cell_subtype = "",
      manual_confidence = "",
      manual_notes = ""
    )]

    final_marker_parts[[broad_type]] <- markers
    final_template_parts[[broad_type]] <- decision
    evidence$scores[, `:=`(
      dataset_id = dataset_id,
      parent_broad_type = broad_type,
      clustering_strategy = strategy
    )]
    final_evidence_parts[[broad_type]] <- evidence$scores

    assignment <- data.table::data.table(
      dataset_id = dataset_id,
      cell_id = colnames(sub_obj),
      parent_broad_type = broad_type,
      clustering_strategy = strategy,
      final_cluster = final_cluster,
      sample_id = as.character(sub_obj$sample_id),
      patient_id = as.character(sub_obj$patient_id),
      timepoint = as.character(sub_obj$timepoint)
    )
    final_assignment_parts[[broad_type]] <- assignment

    rm(sub_obj, result, markers, evidence, cycle, decision, assignment)
    gc()
  }

  final_markers <- data.table::rbindlist(final_marker_parts, fill = TRUE)
  final_template <- data.table::rbindlist(final_template_parts, fill = TRUE)
  final_assignments <- data.table::rbindlist(final_assignment_parts, fill = TRUE)
  final_evidence <- data.table::rbindlist(final_evidence_parts, fill = TRUE)

  if (!nrow(final_template)) {
    stop(dataset_id, ": no resolved broad cell type produced annotation clusters")
  }

  preferred_cols <- c(
    "dataset_id", "parent_broad_type", "clustering_strategy", "cluster",
    "n_cells", "suggested_broad_cell_type", "suggested_family",
    "canonical_evidence_strength",
    "canonical_markers_present", "second_candidate", "second_support_n",
    "cycling_state", "cycling_markers_present", "n_significant_markers",
    "marker_status", "lineage_conflict_flag",
    "top20_significant_markers", "top50_significant_markers",
    "manual_cell_type", "manual_cell_subtype",
    "manual_confidence", "manual_notes"
  )
  preferred_cols <- intersect(preferred_cols, names(final_template))
  data.table::setcolorder(
    final_template,
    c(preferred_cols, setdiff(names(final_template), preferred_cols))
  )

  # Hard gate: annotation-ready marker exports must contain no nonsignificant row.
  if (nrow(final_markers) &&
      all(c("p_val_adj", "avg_log2FC", "pct.1") %in% names(final_markers)) &&
      any(!is.finite(final_markers$p_val_adj) |
          final_markers$p_val_adj >= 0.05 |
          final_markers$avg_log2FC <= 0.25 |
          final_markers$pct.1 < 0.20)) {
    stop(dataset_id, ": nonsignificant marker leaked into final marker output")
  }

  if (any(final_template$parent_broad_type == "Cycling_like", na.rm = TRUE)) {
    stop(dataset_id, ": Cycling_like must not be exported as a broad lineage")
  }

  write_csv(
    final_template,
    file.path(dataset_out, "annotation_ready_cluster_template.csv")
  )
  write_csv_gz(
    final_markers,
    file.path(dataset_out, "annotation_ready_significant_markers.csv.gz")
  )
  write_csv(
    final_evidence,
    file.path(dataset_out, "annotation_ready_canonical_evidence.csv")
  )
  write_csv_gz(
    final_assignments,
    file.path(dataset_out, "annotation_ready_cell_assignments_sampled.csv.gz")
  )

  global_meta_export <- merge(
    global_meta,
    final_assignments[
      ,
      .(cell_id, final_cluster, clustering_strategy)
    ],
    by = "cell_id",
    all.x = TRUE
  )
  write_csv_gz(
    global_meta_export,
    file.path(dataset_out, "annotation_ready_global_cell_map_sampled.csv.gz")
  )

  # Plot global discovery result.
  global_obj$marker_ready_broad_type <- global_meta[
    match(colnames(global_obj), cell_id),
    suggested_broad_cell_type
  ]
  umap <- Embeddings(global_obj, "umap.global_A_uncorrected")
  plot_df <- data.frame(
    UMAP_1 = umap[, 1L],
    UMAP_2 = umap[, 2L],
    broad_type = global_obj$marker_ready_broad_type,
    sample_id = global_obj$sample_id
  )
  p_broad <- ggplot(
    plot_df,
    aes(UMAP_1, UMAP_2, color = broad_type)
  ) +
    geom_point(size = 0.18, alpha = 0.75) +
    theme_bw() +
    labs(
      title = paste(dataset_id, "marker-ready broad cell types"),
      color = "Broad type"
    )
  save_plot_pair(
    p_broad,
    file.path(dataset_out, "UMAP_marker_ready_broad_types"),
    10, 7
  )

  # Canonical marker dotplot across final clusters.
  if (nrow(final_assignments)) {
    mapped <- final_assignments[
      match(colnames(global_obj), cell_id),
      final_cluster
    ]
    keep <- !is.na(mapped) & nzchar(mapped)
    if (sum(keep) >= 100L && data.table::uniqueN(mapped[keep]) >= 2L) {
      dot_obj <- subset(global_obj, cells = colnames(global_obj)[keep])
      Idents(dot_obj) <- factor(mapped[keep])
      panel <- unique(c(unlist(marker_panels, use.names = FALSE),
                        cycling_markers))
      panel <- panel[panel %in% rownames(dot_obj)]
      if (length(panel)) {
        p_dot <- DotPlot(dot_obj, features = panel) +
          RotatedAxis() +
          labs(
            title = paste(dataset_id, "canonical markers by final cluster"),
            x = "Canonical marker",
            y = "Final cluster"
          )
        save_plot_pair(
          p_dot,
          file.path(dataset_out, "canonical_marker_dotplot"),
          18, 12
        )
      }
      rm(dot_obj)
      gc()
    }
  }

  summary_lines <- c(
    paste0("# ", dataset_id, " marker-ready annotation summary"),
    "",
    paste0("- Repaired-QC eligible cells: ", n_repaired_qc_cells),
    paste0("- Balanced discovery cells used: ", ncol(global_obj)),
    "- Global broad types were assigned from adjusted-P significant RNA markers.",
    "- Cycling is exported only as a state, never as a broad lineage.",
    "- Epithelial subclustering uses A_uncorrected.",
    if (dataset_id == "GSE158722") {
      "- Non-epithelial Harmony, when feasible, uses patient_id only; sample_id/timepoint is not corrected."
    } else {
      "- Non-epithelial Harmony, when feasible, uses sample_id."
    },
    "- All annotation-ready marker rows satisfy p_val_adj < 0.05, avg_log2FC > 0.25, and pct.1 >= 0.20.",
    paste0("- Annotation clusters exported: ", nrow(final_template)),
    paste0("- Clusters READY_FOR_ANNOTATION: ",
           sum(final_template$marker_status == "READY_FOR_ANNOTATION")),
    paste0("- Clusters needing marker review: ",
           sum(final_template$marker_status != "READY_FOR_ANNOTATION")),
    paste0("- Lineage-conflict clusters: ",
           sum(final_template$lineage_conflict_flag, na.rm = TRUE)),
    "",
    "Primary file for manual annotation:",
    "`annotation_ready_cluster_template.csv`",
    "",
    "Fill only these columns:",
    "`manual_cell_type`, `manual_cell_subtype`, `manual_confidence`, `manual_notes`."
  )
  writeLines(
    summary_lines,
    file.path(dataset_out, "annotation_ready_summary.md"),
    useBytes = TRUE
  )

  capture.output(
    sessionInfo(),
    file = file.path(dataset_out, "sessionInfo.txt")
  )

  rm(global_obj, global_markers, global_evidence, global_cycle,
     global_decisions, global_meta, final_markers, final_template,
     final_assignments, final_evidence)
  gc()
}

message("Marker-ready annotation workflow complete: ", output_root)
