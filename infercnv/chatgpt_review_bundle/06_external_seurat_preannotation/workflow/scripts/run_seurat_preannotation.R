#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1]]
}

has_flag <- function(flag) flag %in% args
config_path <- arg_value("--config", "config/five_external_datasets.yaml")
audit_only_cli <- has_flag("--audit-only")
dataset_filter_arg <- arg_value("--datasets", "")

required <- c("yaml", "Seurat", "Matrix", "dplyr", "ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                           FUN.VALUE = logical(1))]
if (length(missing)) {
  stop("Missing required R packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cfg <- yaml::read_yaml(config_path)
set.seed(as.integer(cfg$project$random_seed %||% 20260713))

safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE)
}

write_status <- function(x, path) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    dput(x, file = sub("\\.json$", ".dput", path))
  }
}

read_dataset <- function(ds) {
  type <- ds$input_type
  path <- ds$input_path

  if (type == "10x_h5") {
    if (!file.exists(path)) stop("Input not found: ", path)
    mat <- Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
    if (is.list(mat)) mat <- mat[[1]]
    return(CreateSeuratObject(mat, project = ds$dataset_id,
                              min.cells = 0, min.features = 0))
  }

  if (type == "10x_dir") {
    if (!dir.exists(path)) stop("Input directory not found: ", path)
    mat <- Read10X(path, gene.column = 2, unique.features = TRUE)
    if (is.list(mat)) mat <- mat[[1]]
    return(CreateSeuratObject(mat, project = ds$dataset_id,
                              min.cells = 0, min.features = 0))
  }

  if (type == "rds_seurat") {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("RDS is not a Seurat object")
    return(obj)
  }

  if (type == "rds_sce") {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment and SummarizedExperiment are required")
    }
    sce <- readRDS(path)
    counts <- SummarizedExperiment::assay(sce, "counts")
    obj <- CreateSeuratObject(counts, project = ds$dataset_id)
    md <- as.data.frame(SummarizedExperiment::colData(sce))
    if (!is.null(rownames(md)) && all(colnames(obj) %in% rownames(md))) {
      obj <- AddMetaData(obj, md[colnames(obj), , drop = FALSE])
    }
    return(obj)
  }

  if (type == "matrix_rds") {
    mat <- readRDS(path)
    if (!inherits(mat, "Matrix") && !is.matrix(mat)) {
      stop("matrix_rds must contain a matrix")
    }
    return(CreateSeuratObject(mat, project = ds$dataset_id))
  }

  stop("Unsupported input_type: ", type)
}

get_counts <- function(obj) {
  assay <- DefaultAssay(obj)
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
  )
}

integer_like_fraction <- function(mat, n = 100000L) {
  vals <- if (inherits(mat, "sparseMatrix")) mat@x else as.numeric(mat)
  if (!length(vals)) return(1)
  vals <- vals[sample.int(length(vals), min(n, length(vals)))]
  mean(abs(vals - round(vals)) < 1e-8)
}

mad_range <- function(x, lower_mult, upper_mult) {
  med <- median(x, na.rm = TRUE)
  md <- mad(x, center = med, constant = 1, na.rm = TRUE)
  c(lower = med - lower_mult * md, upper = med + upper_mult * md)
}

run_one <- function(ds, cfg, audit_only = FALSE) {
  id <- ds$dataset_id
  root <- file.path(cfg$project$output_root, id)
  subdirs <- c("00_input_audit", "01_qc", "02_clustering", "03_markers",
               "04_manual_annotation", "logs", "objects")
  invisible(lapply(file.path(root, subdirs), safe_dir))
  status_file <- file.path(root, "logs", "run_status.json")

  status <- list(dataset_id = id, status = "STARTED",
                 started_at = as.character(Sys.time()))
  write_status(status, status_file)

  obj <- read_dataset(ds)
  obj$dataset_id <- id
  obj$biological_role <- ds$biological_role
  obj$cell_type_manual <- NA_character_
  obj$cell_subtype_manual <- NA_character_

  counts <- get_counts(obj)
  audit <- data.frame(
    dataset_id = id,
    input_path = ds$input_path,
    input_type = ds$input_type,
    file_size_bytes = if (file.exists(ds$input_path))
      file.info(ds$input_path)$size else NA,
    n_genes = nrow(counts),
    n_cells = ncol(counts),
    duplicated_genes = sum(duplicated(rownames(counts))),
    duplicated_cells = sum(duplicated(colnames(counts))),
    integer_like_fraction = integer_like_fraction(counts),
    metadata_columns = paste(colnames(obj@meta.data), collapse = ";"),
    stringsAsFactors = FALSE
  )
  write.csv(audit, file.path(root, "00_input_audit", "input_audit.csv"),
            row.names = FALSE)

  if (nrow(counts) == 0 || ncol(counts) == 0) stop(id, ": empty matrix")
  if (audit$duplicated_cells > 0) stop(id, ": duplicated cell IDs")
  if (isTRUE(ds$declared_raw_counts) && audit$integer_like_fraction < 0.99) {
    stop(id, ": declared raw counts are not integer-like")
  }

  if (audit_only || isTRUE(cfg$project$audit_only)) {
    status$status <- "AUDIT_COMPLETE"
    status$finished_at <- as.character(Sys.time())
    write_status(status, status_file)
    return(status)
  }

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["percent.HB"]] <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGQZ]")

  sample_col <- ds$sample_column %||% ""
  if (nzchar(sample_col) && sample_col %in% colnames(obj@meta.data)) {
    obj$analysis_sample_id <- as.character(obj@meta.data[[sample_col]])
    qc_scope <- "sample_wise"
  } else {
    obj$analysis_sample_id <- id
    qc_scope <- "dataset_wise"
  }

  md <- obj@meta.data
  groups <- split(seq_len(nrow(md)), md$analysis_sample_id)
  keep <- rep(FALSE, nrow(md))
  thresholds <- list()

  for (g in names(groups)) {
    ix <- groups[[g]]
    z <- md[ix, , drop = FALSE]
    fr <- mad_range(z$nFeature_RNA, cfg$qc$mad_lower, cfg$qc$mad_upper)
    cr <- mad_range(z$nCount_RNA, cfg$qc$mad_lower, cfg$qc$mad_upper)
    mr <- mad_range(z$percent.mt, 0, cfg$qc$mad_upper)

    min_f <- max(cfg$qc$min_features, fr["lower"])
    max_f <- fr["upper"]
    min_c <- max(cfg$qc$min_counts, cr["lower"])
    max_c <- cr["upper"]
    max_mt <- min(cfg$qc$max_percent_mt, mr["upper"])

    keep[ix] <- z$nFeature_RNA >= min_f & z$nFeature_RNA <= max_f &
      z$nCount_RNA >= min_c & z$nCount_RNA <= max_c &
      z$percent.mt <= max_mt

    thresholds[[length(thresholds) + 1]] <- data.frame(
      sample_id = g, qc_scope = qc_scope,
      min_features = min_f, max_features = max_f,
      min_counts = min_c, max_counts = max_c,
      max_percent_mt = max_mt,
      n_before = nrow(z), n_after_basic_qc = sum(keep[ix])
    )
  }

  write.csv(bind_rows(thresholds),
            file.path(root, "01_qc", "qc_thresholds.csv"),
            row.names = FALSE)

  obj$basic_qc_pass <- keep
  obj2 <- subset(obj, cells = rownames(obj@meta.data)[keep])

  doublet_status <- "not_requested"
  if (isTRUE(cfg$qc$remove_doublets)) {
    available <- requireNamespace("scDblFinder", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)
    if (available) {
      sample_ids <- unique(as.character(obj2$analysis_sample_id))
      doublet_score <- setNames(rep(NA_real_, ncol(obj2)), colnames(obj2))
      doublet_class <- setNames(rep(NA_character_, ncol(obj2)), colnames(obj2))
      for (sample_id in sample_ids) {
        cells <- colnames(obj2)[obj2$analysis_sample_id == sample_id]
        message(id, ": scDblFinder sample ", sample_id,
                " (", length(cells), " cells)")
        sample_obj <- subset(obj2, cells = cells)
        sample_sce <- as.SingleCellExperiment(sample_obj)
        sample_sce <- scDblFinder::scDblFinder(sample_sce)
        sample_cd <- SummarizedExperiment::colData(sample_sce)
        doublet_score[cells] <- sample_cd$scDblFinder.score
        doublet_class[cells] <- as.character(sample_cd$scDblFinder.class)
        rm(sample_obj, sample_sce, sample_cd)
        gc()
      }
      obj2$scDblFinder.score <- unname(doublet_score[colnames(obj2)])
      obj2$scDblFinder.class <- unname(doublet_class[colnames(obj2)])
      obj2 <- subset(obj2, subset = scDblFinder.class == "singlet")
      doublet_status <- "completed"
    } else if (isTRUE(cfg$qc$allow_doublet_skip)) {
      doublet_status <- "skipped_package_unavailable"
    } else {
      stop(id, ": scDblFinder unavailable")
    }
  }

  retention <- data.frame(
    dataset_id = id, n_input = ncol(obj),
    n_after_basic_qc = sum(keep),
    n_after_doublet_filter = ncol(obj2),
    retained_fraction = ncol(obj2) / ncol(obj),
    doublet_status = doublet_status
  )
  write.csv(retention, file.path(root, "01_qc", "qc_cell_retention.csv"),
            row.names = FALSE)
  write_csv_gz(obj2@meta.data,
               file.path(root, "01_qc", "qc_metadata.csv.gz"))

  if (ncol(obj2) < cfg$analysis$min_cells_after_qc) {
    stop(id, ": fewer than minimum cells after QC")
  }

  obj2 <- NormalizeData(obj2, verbose = FALSE)
  obj2 <- FindVariableFeatures(obj2, selection.method = "vst",
                               nfeatures = cfg$analysis$variable_features,
                               verbose = FALSE)
  obj2 <- ScaleData(obj2, features = VariableFeatures(obj2), verbose = FALSE)
  obj2 <- RunPCA(obj2, npcs = cfg$analysis$n_pcs, verbose = FALSE)

  dims <- seq_len(min(cfg$analysis$dims_use,
                      ncol(Embeddings(obj2, "pca"))))
  obj2 <- FindNeighbors(obj2, dims = dims, verbose = FALSE)

  for (res in unlist(cfg$analysis$resolutions)) {
    obj2 <- FindClusters(obj2, resolution = res,
                         random.seed = cfg$project$random_seed,
                         verbose = FALSE)
  }
  obj2 <- RunUMAP(obj2, dims = dims,
                  seed.use = cfg$project$random_seed, verbose = FALSE)

  primary_res <- as.numeric(cfg$analysis$primary_resolution)
  candidates <- grep(paste0("_snn_res\\.", primary_res, "$"),
                     colnames(obj2@meta.data), value = TRUE)
  if (!length(candidates)) stop(id, ": primary resolution column not found")
  primary_col <- candidates[[1]]
  obj2$seurat_clusters <- as.character(obj2@meta.data[[primary_col]])
  Idents(obj2) <- "seurat_clusters"

  cc <- as.data.frame(table(obj2$seurat_clusters))
  colnames(cc) <- c("seurat_cluster", "n_cells")
  write.csv(cc, file.path(root, "02_clustering", "cluster_cell_counts.csv"),
            row.names = FALSE)

  cs <- as.data.frame(table(obj2$analysis_sample_id,
                            obj2$seurat_clusters))
  colnames(cs) <- c("sample_id", "seurat_cluster", "n_cells")
  write.csv(cs,
            file.path(root, "02_clustering", "cluster_by_sample_counts.csv"),
            row.names = FALSE)

  p1 <- DimPlot(obj2, group.by = "seurat_clusters",
                label = TRUE, repel = TRUE)
  ggsave(file.path(root, "02_clustering",
                   "umap_primary_resolution.pdf"), plot = p1,
         width = 8, height = 6)

  p2 <- DimPlot(obj2, group.by = "analysis_sample_id")
  ggsave(file.path(root, "02_clustering", "umap_by_sample.pdf"),
         plot = p2, width = 8, height = 6)

  for (res in unlist(cfg$analysis$resolutions)) {
    z <- grep(paste0("_snn_res\\.", res, "$"),
              colnames(obj2@meta.data), value = TRUE)
    if (length(z)) {
      p <- DimPlot(obj2, group.by = z[[1]], label = TRUE, repel = TRUE)
      ggsave(file.path(root, "02_clustering",
                       paste0("umap_resolution_", res, ".pdf")), plot = p,
             width = 8, height = 6)
    }
  }

  if (isTRUE(cfg$analysis$run_clustree) &&
      requireNamespace("clustree", quietly = TRUE)) {
    prefix <- sub(paste0(primary_res, "$"), "", primary_col)
    p <- clustree::clustree(obj2@meta.data, prefix = prefix) +
      ggplot2::guides(edge_colour = "none")
    ggsave(file.path(root, "02_clustering", "clustree.pdf"),
           plot = p, width = 10, height = 8)
  }

  DefaultAssay(obj2) <- "RNA"
  obj2 <- tryCatch(JoinLayers(obj2), error = function(e) obj2)
  saveRDS(obj2,
          file.path(root, "objects",
                    paste0(id, "_preannotation.rds")),
          compress = FALSE)

  markers <- FindAllMarkers(
    obj2, assay = "RNA",
    only.pos = isTRUE(cfg$markers$only_pos),
    test.use = cfg$markers$test_use,
    min.pct = cfg$markers$min_pct,
    logfc.threshold = cfg$markers$logfc_threshold,
    return.thresh = cfg$markers$adjusted_p_threshold,
    verbose = FALSE
  )
  if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
  write_csv_gz(markers,
               file.path(root, "03_markers",
                         "all_cluster_markers.csv.gz"))

  for (n in unlist(cfg$markers$export_top_n)) {
    top <- markers |>
      group_by(cluster) |>
      arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) |>
      slice_head(n = n) |>
      ungroup()
    write.csv(top,
              file.path(root, "03_markers",
                        paste0("top", n, "_markers_per_cluster.csv")),
              row.names = FALSE)
  }

  avg <- AverageExpression(obj2, assays = "RNA", slot = "data",
                           group.by = "seurat_clusters",
                           verbose = FALSE)$RNA
  avg <- data.frame(gene = rownames(avg), as.data.frame(avg),
                    check.names = FALSE)
  write_csv_gz(avg,
               file.path(root, "03_markers",
                         "cluster_average_expression.csv.gz"))

  panel <- c("EPCAM","KRT8","KRT18","KRT19","MSLN","WFDC2",
             "PTPRC","CD3D","CD3E","TRBC1","NKG7","GNLY",
             "CD79A","MS4A1","MZB1","JCHAIN",
             "LYZ","LST1","TYROBP","FCER1G","C1QA","C1QB","C1QC","SPP1",
             "COL1A1","COL1A2","DCN","COL3A1",
             "PECAM1","VWF","KDR","MKI67","TOP2A")
  panel <- panel[panel %in% rownames(obj2)]
  if (length(panel)) {
    p <- DotPlot(obj2, features = panel,
                 group.by = "seurat_clusters") + RotatedAxis()
    ggsave(file.path(root, "03_markers",
                     "broad_marker_dotplot.pdf"), plot = p,
           width = 16, height = 8)
  }

  top10 <- markers |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) |>
    slice_head(n = 10) |>
    pull(gene) |>
    unique()
  top10 <- top10[top10 %in% rownames(obj2)]
  if (length(top10)) {
    heatmap_cells <- WhichCells(
      obj2, downsample = 200,
      seed = cfg$project$random_seed
    )
    heatmap_obj <- subset(obj2, cells = heatmap_cells)
    p <- DoHeatmap(heatmap_obj, features = top10,
                   group.by = "seurat_clusters", raster = TRUE) +
      NoLegend()
    ggsave(file.path(root, "03_markers",
                     "top_marker_heatmap.pdf"), plot = p,
           width = 14, height = 12)
  }

  top_text <- markers |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) |>
    slice_head(n = 20) |>
    summarise(top_markers = paste(gene, collapse = ";"),
              .groups = "drop")

  manual <- cc |>
    left_join(top_text, by = c("seurat_cluster" = "cluster")) |>
    mutate(dataset_id = id, cell_type_manual = "",
           cell_subtype_manual = "", confidence = "", notes = "") |>
    select(dataset_id, seurat_cluster, n_cells, top_markers,
           cell_type_manual, cell_subtype_manual, confidence, notes)
  write.csv(manual,
            file.path(root, "04_manual_annotation",
                      "manual_annotation_template.csv"),
            row.names = FALSE)

  saveRDS(obj2,
          file.path(root, "objects",
                    paste0(id, "_preannotation.rds")),
          compress = FALSE)
  capture.output(sessionInfo(),
                 file = file.path(root, "logs", "sessionInfo.txt"))

  status$status <- "PREANNOTATION_COMPLETE_WAITING_FOR_MANUAL_CELLTYPE"
  status$finished_at <- as.character(Sys.time())
  status$n_input <- ncol(obj)
  status$n_after_qc <- ncol(obj2)
  status$primary_resolution <- primary_res
  status$n_clusters <- length(unique(obj2$seurat_clusters))
  status$n_marker_rows <- nrow(markers)
  status$doublet_status <- doublet_status
  write_status(status, status_file)
  status
}

out_root <- safe_dir(cfg$project$output_root)
audit_only <- audit_only_cli || isTRUE(cfg$project$audit_only)
datasets <- Filter(function(x) isTRUE(x$enabled), cfg$datasets)
if (nzchar(dataset_filter_arg)) {
  requested <- trimws(strsplit(dataset_filter_arg, ",", fixed = TRUE)[[1]])
  datasets <- Filter(function(x) x$dataset_id %in% requested, datasets)
  missing_requested <- setdiff(requested,
                               vapply(datasets, `[[`, character(1), "dataset_id"))
  if (length(missing_requested)) {
    stop("Requested dataset(s) are not enabled/configured: ",
         paste(missing_requested, collapse = ", "))
  }
}

results <- lapply(datasets, function(ds) {
  message("===== ", ds$dataset_id, " =====")
  result <- tryCatch(
    run_one(ds, cfg, audit_only),
    error = function(e) {
      log_dir <- safe_dir(file.path(out_root, ds$dataset_id, "logs"))
      z <- list(dataset_id = ds$dataset_id, status = "BLOCKED",
                error = conditionMessage(e),
                finished_at = as.character(Sys.time()))
      write_status(z, file.path(log_dir, "run_status.json"))
      z
    }
  )
  gc()
  result
})

summary <- bind_rows(lapply(results, as.data.frame))
write.csv(summary, file.path(out_root, "run_summary.csv"), row.names = FALSE)
if (any(summary$status == "BLOCKED")) quit(status = 2)
