#!/usr/bin/env Rscript

required <- c("yaml", "data.table", "Matrix", "Seurat", "ggplot2")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
# The largest 30k-cell epithelial RPCA reference exports ~4.46 GiB of globals
# even under a sequential plan. Six GiB admits that audited run while retaining
# a finite guard on this 8 GiB workstation.
options(future.globals.maxSize = 6 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else "."
source(file.path(script_dir, "_diagnostics_v2_common.R"))

z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
datasets <- split_arg(arg_value("--datasets", "GSE154600,GSE158722"))
assert_datasets(datasets, c("GSE154600", "GSE158722"))
skip_harmony <- has_flag("--skip-harmony")
skip_rpca <- has_flag("--skip-rpca")
resume <- has_flag("--resume")
seed <- as.integer(cfg$project$random_seed)
set.seed(seed)

cluster_and_export <- function(obj, reduction, strategy, dataset_id, lineage,
                               out, dims_use) {
  dims <- seq_len(min(dims_use, ncol(Embeddings(obj, reduction))))
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = as.numeric(cfg$analysis$primary_resolution),
                      random.seed = seed, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = reduction, dims = dims, seed.use = seed,
                 reduction.name = paste0("umap.", strategy), verbose = FALSE)
  emb <- Embeddings(obj, paste0("umap.", strategy))
  eval_emb <- Embeddings(obj, reduction)[, dims, drop = FALSE]
  md <- obj@meta.data
  clusters <- as.character(Idents(obj))
  export <- data.table::data.table(
    dataset_id = dataset_id, analysis_scope_lineage = lineage,
    provisional_broad_lineage = as.character(md$provisional_broad_lineage),
    strategy = strategy, cell_id = rownames(md),
    sample_id = as.character(md$sample_id),
    patient_id = as.character(md$patient_id),
    timepoint = as.character(md$timepoint),
    cluster = clusters, UMAP_1 = emb[, 1L], UMAP_2 = emb[, 2L]
  )
  write_csv_gz(export, file.path(out, "cell_embedding_and_clusters.csv.gz"))
  eval_export <- data.table::data.table(
    dataset_id = dataset_id, analysis_scope_lineage = lineage,
    strategy = strategy, cell_id = rownames(md),
    sample_id = as.character(md$sample_id),
    patient_id = as.character(md$patient_id),
    timepoint = as.character(md$timepoint),
    author_label_source = as.character(md$author_label_source),
    nCount_RNA = as.numeric(md$nCount_RNA),
    nFeature_RNA = as.numeric(md$nFeature_RNA),
    percent.mt = as.numeric(md$percent.mt),
    doublet_score = as.numeric(md$doublet_score),
    provisional_broad_lineage = as.character(md$provisional_broad_lineage),
    cluster = clusters
  )
  eval_export <- cbind(eval_export, as.data.frame(eval_emb))
  write_csv_gz(eval_export, file.path(out, "evaluation_reduction.csv.gz"))
  write_csv(dominance_metrics(export$cluster, export$sample_id, dataset_id,
                              strategy, lineage),
            file.path(out, "cluster_sample_metrics.csv"))
  p_cluster <- ggplot(export, aes(UMAP_1, UMAP_2, color = cluster)) +
    geom_point(size = 0.15, alpha = 0.7) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(title = paste(dataset_id, lineage, strategy, "clusters"))
  p_sample <- ggplot(export, aes(UMAP_1, UMAP_2, color = sample_id)) +
    geom_point(size = 0.15, alpha = 0.65) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(title = paste(dataset_id, lineage, strategy, "samples"))
  save_plot_pair(p_cluster, file.path(out, "UMAP_by_cluster"), 8, 6)
  save_plot_pair(p_sample, file.path(out, "UMAP_by_sample"), 9, 6)
  if (data.table::uniqueN(export$patient_id) > 1L) {
    p_patient <- ggplot(export, aes(UMAP_1, UMAP_2, color = patient_id)) +
      geom_point(size = 0.15, alpha = 0.65) + theme_bw() +
      guides(color = guide_legend(override.aes = list(size = 2))) +
      labs(title = paste(dataset_id, lineage, strategy, "patients"))
    save_plot_pair(p_patient, file.path(out, "UMAP_by_patient"), 9, 6)
  }
  invisible(obj)
}

run_baseline <- function(obj, dataset_id, lineage, out) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = as.integer(cfg$analysis$variable_features),
                              selection.method = "vst", verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj),
                npcs = as.integer(cfg$analysis$n_pcs), verbose = FALSE)
  cluster_and_export(obj, "pca", "A_uncorrected", dataset_id, lineage, out,
                     as.integer(cfg$analysis$dims_use))
}

run_harmony <- function(obj, group, strategy, dataset_id, lineage, out) {
  if (!requireNamespace("harmony", quietly = TRUE)) stop("harmony package unavailable")
  # harmony 2.0.5 names this argument reduction.use; it is the same PCA-only
  # correction requested by the workflow and does not alter RNA counts/data.
  obj <- harmony::RunHarmony(obj, group.by.vars = group,
                             reduction.use = "pca", verbose = FALSE)
  cluster_and_export(obj, "harmony", strategy, dataset_id, lineage, out,
                     as.integer(cfg$analysis$dims_use))
}

run_rpca <- function(base_obj, dataset_id, lineage, out) {
  groups <- split(colnames(base_obj), as.character(base_obj$sample_id))
  if (min(lengths(groups)) < 10L) stop("RPCA requires >=10 cells in every sample within lineage")
  objects <- SplitObject(base_obj, split.by = "sample_id")
  objects <- lapply(objects, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = min(2000L, nrow(x)), verbose = FALSE)
    x
  })
  features <- SelectIntegrationFeatures(objects, nfeatures = min(2000L, nrow(base_obj)))
  min_cells <- min(vapply(objects, ncol, numeric(1)))
  npcs <- min(as.integer(cfg$analysis$n_pcs), min_cells - 1L, length(features) - 1L)
  dims <- seq_len(min(as.integer(cfg$analysis$dims_use), npcs))
  objects <- lapply(objects, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    RunPCA(x, features = features, npcs = npcs, verbose = FALSE)
  })
  reference <- NULL
  if (length(objects) > 10L) {
    object_sizes <- vapply(objects, ncol, numeric(1))
    reference <- head(order(object_sizes, decreasing = TRUE), 5L)
  }
  anchors <- FindIntegrationAnchors(object.list = objects,
                                    anchor.features = features,
                                    reduction = "rpca", dims = dims,
                                    reference = reference,
                                    verbose = FALSE)
  # k.weight must stay below the smallest lineage-by-sample object. The fixed
  # value 30 fails for valid small Endothelial/Cycling strata (for example 26
  # cells) with a non-multiple replacement error inside IntegrateData.
  k_weight <- min(30L, min_cells - 1L)
  integrated <- IntegrateData(anchorset = anchors, dims = dims,
                              k.weight = k_weight, verbose = FALSE)
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs,
                       reduction.name = "integrated.rpca", verbose = FALSE)
  cluster_and_export(integrated, "integrated.rpca", "C_RPCA", dataset_id,
                     lineage, out, as.integer(cfg$analysis$dims_use))
}

status_path <- file.path(data_root, "diagnostics_v2", "strategy_comparison",
                         "strategy_run_status.csv")
existing_status <- if (file.exists(status_path)) {
  data.table::fread(status_path, showProgress = FALSE)[!dataset_id %in% datasets]
} else data.table::data.table()
statuses <- if (nrow(existing_status)) split(as.data.frame(existing_status),
                                             seq_len(nrow(existing_status))) else list()
for (dataset_id in datasets) {
  dataset_id_value <- dataset_id
  input_dir <- file.path(data_root, "diagnostics_v2", "objects", dataset_id,
                         "lineage_inputs")
  paths <- sort(list.files(input_dir, pattern = "_strategy_input\\.rds$", full.names = TRUE))
  if (!length(paths)) stop(dataset_id, ": no lineage strategy inputs from step 04")
  for (path in paths) {
    payload <- readRDS(path)
    lineage <- payload$lineage
    message("===== ", dataset_id, " / ", lineage, " =====")
    md <- payload$metadata
    rownames(md) <- md$cell_id
    obj <- CreateSeuratObject(payload$counts, project = dataset_id,
                              min.cells = 0, min.features = 0, meta.data = md)
    DefaultAssay(obj) <- "RNA"
    rm(payload); gc()

    base_dir <- file.path(data_root, "diagnostics_v2", dataset_id,
                          "04_strategies", lineage)
    strategy_specs <- list(
      list(id = "A_uncorrected", skip = FALSE,
           fun = function(out) run_baseline(obj, dataset_id, lineage, out)),
      list(id = "B_harmony_sample", skip = skip_harmony,
           fun = function(out) {
             baseline_path <- file.path(base_dir, "A_uncorrected", "baseline_object_temp.rds")
             base <- if (file.exists(baseline_path)) readRDS(baseline_path) else run_baseline(obj, dataset_id, lineage, file.path(base_dir, "A_uncorrected"))
             run_harmony(base, "sample_id", "B_harmony_sample", dataset_id, lineage, out)
           })
    )

    # Run baseline once and retain only a temporary normalized/PCA object until
    # Harmony completes. It is never included in the GitHub review package.
    harmony_groups <- list(B_harmony_sample = "sample_id")
    if (dataset_id == "GSE158722" &&
        data.table::uniqueN(md$patient_id) < data.table::uniqueN(md$sample_id)) {
      harmony_groups$B_harmony_patient <- "patient_id"
    }
    harmony_complete <- vapply(names(harmony_groups), function(strategy) {
      file.exists(file.path(base_dir, strategy,
                            "cell_embedding_and_clusters.csv.gz"))
    }, logical(1))
    need_baseline_obj <- !skip_harmony && !(resume && all(harmony_complete))
    a_out <- file.path(base_dir, "A_uncorrected")
    a_complete <- file.exists(file.path(a_out, "cell_embedding_and_clusters.csv.gz"))
    baseline_obj <- NULL
    if (!(resume && a_complete)) {
      result <- tryCatch({
        baseline_obj <- run_baseline(obj, dataset_id, lineage, a_out)
        list(status = "COMPLETE", message = "")
      }, error = function(e) list(status = "FAILED", message = conditionMessage(e)))
    } else {
      result <- list(status = "RESUMED_EXISTING", message = "")
      if (need_baseline_obj) {
        baseline_obj <- NormalizeData(obj, verbose = FALSE)
        baseline_obj <- FindVariableFeatures(baseline_obj, nfeatures = as.integer(cfg$analysis$variable_features), verbose = FALSE)
        baseline_obj <- ScaleData(baseline_obj, features = VariableFeatures(baseline_obj), verbose = FALSE)
        baseline_obj <- RunPCA(baseline_obj, npcs = as.integer(cfg$analysis$n_pcs), verbose = FALSE)
      }
    }
    statuses[[length(statuses) + 1L]] <- data.frame(dataset_id, lineage,
      strategy = "A_uncorrected", status = result$status, message = result$message)

    for (strategy in names(harmony_groups)) {
      out <- file.path(base_dir, strategy)
      complete <- file.exists(file.path(out, "cell_embedding_and_clusters.csv.gz"))
      if (skip_harmony) {
        result <- list(status = "SKIPPED_BY_FLAG", message = "--skip-harmony")
      } else if (resume && complete) {
        result <- list(status = "RESUMED_EXISTING", message = "")
      } else if (is.null(baseline_obj)) {
        result <- list(status = "BLOCKED_BASELINE_FAILED", message = "baseline unavailable")
      } else {
        result <- tryCatch({
          run_harmony(baseline_obj, harmony_groups[[strategy]], strategy,
                      dataset_id, lineage, out)
          list(status = "COMPLETE", message = "")
        }, error = function(e) list(status = "FAILED_CONTINUED", message = conditionMessage(e)))
      }
      statuses[[length(statuses) + 1L]] <- data.frame(dataset_id, lineage,
        strategy = strategy, status = result$status, message = result$message)
    }
    rm(baseline_obj); gc()

    out <- file.path(base_dir, "C_RPCA")
    complete <- file.exists(file.path(out, "cell_embedding_and_clusters.csv.gz"))
    if (lineage == "Combined_broad_lineages") {
      result <- list(
        status = "SKIPPED_MEMORY_GUARD_COMBINED_SCOPE",
        message = "Combined scope is an extra cross-lineage audit; RPCA remains enabled within each individual broad lineage."
      )
    } else if (skip_rpca) {
      result <- list(status = "SKIPPED_BY_FLAG", message = "--skip-rpca")
    } else if (resume && complete) {
      result <- list(status = "RESUMED_EXISTING", message = "")
    } else {
      result <- tryCatch({
        run_rpca(obj, dataset_id, lineage, out)
        list(status = "COMPLETE", message = "")
      }, error = function(e) list(status = "FAILED_CONTINUED", message = conditionMessage(e)))
    }
    statuses[[length(statuses) + 1L]] <- data.frame(dataset_id, lineage,
      strategy = "C_RPCA", status = result$status, message = result$message)
    rm(obj); gc()
    status_table <- data.table::rbindlist(statuses, fill = TRUE)
    write_csv(status_table, status_path)
    write_csv(status_table[dataset_id == dataset_id_value],
              file.path(data_root, "diagnostics_v2", "strategy_comparison",
                        paste0("strategy_run_status_", dataset_id, ".csv")))
  }
}

capture.output(sessionInfo(), file = file.path(data_root, "diagnostics_v2",
                                                "strategy_comparison", "sessionInfo_05.txt"))
message("A/B/C strategy comparison complete")
