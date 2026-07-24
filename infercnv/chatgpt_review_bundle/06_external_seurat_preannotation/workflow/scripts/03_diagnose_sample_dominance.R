#!/usr/bin/env Rscript

required <- c("yaml", "data.table", "Seurat", "ggplot2")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                            FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else "."
source(file.path(script_dir, "_diagnostics_v2_common.R"))

z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
datasets <- split_arg(arg_value("--datasets", "GSE154600,GSE158722"))
assert_datasets(datasets, c("GSE147082", "GSE154600", "GSE158722"))
forensic_root <- file.path(data_root, "diagnostics_v2", "00_forensic")
all_metrics <- data.table::fread(file.path(forensic_root, "current_cluster_sample_metrics.csv"),
                                 showProgress = FALSE)
pc_audit <- data.table::fread(file.path(forensic_root, "current_pc_sample_association.csv"),
                              showProgress = FALSE)

top_label <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(c(label = NA_character_, fraction = NA_character_))
  tab <- sort(table(x), decreasing = TRUE)
  c(label = names(tab)[[1L]], fraction = as.character(unname(tab[[1L]]) / sum(tab)))
}

for (dataset_id in datasets) {
  message("===== ", dataset_id, " =====")
  dataset_id_value <- dataset_id
  out <- file.path(data_root, "diagnostics_v2", dataset_id, "02_dominance")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  current <- all_metrics[dataset_id == dataset_id_value]
  current[, cluster := as.character(cluster)]
  marker_path <- file.path(data_root, dataset_id, "03_markers",
                           "top50_markers_per_cluster.csv")
  markers <- data.table::fread(marker_path, showProgress = FALSE)
  markers[, cluster := as.character(cluster)]
  marker_text <- markers[order(cluster, -avg_log2FC),
                         .(top_markers = paste(head(gene, 50L), collapse = ";")),
                         by = cluster]

  object_path <- file.path(data_root, dataset_id, "objects",
                           paste0(dataset_id, "_preannotation.rds"))
  size_gb <- file.info(object_path)$size / 1024^3
  can_load <- is.finite(size_gb) && size_gb <= as.numeric(cfg$project$max_object_load_gb)
  obj <- NULL
  if (can_load) obj <- readRDS(object_path)

  diagnostic <- merge(current, marker_text, by = "cluster", all.x = TRUE)
  diagnostic[, `:=`(
    top_author_celltype = NA_character_, top_author_celltype_fraction = NA_real_,
    median_nCount = NA_real_, median_nFeature = NA_real_, median_percent_mt = NA_real_,
    doublet_fraction = NA_real_, median_doublet_score = NA_real_,
    median_S_score = NA_real_, median_G2M_score = NA_real_
  )]

  plot_status <- list()
  if (!is.null(obj)) {
    md <- data.table::as.data.table(obj@meta.data, keep.rownames = "cell_id")
    cluster_col <- if ("seurat_clusters" %in% names(md)) "seurat_clusters" else
      grep("_snn_res\\.0\\.6$", names(md), value = TRUE)[[1L]]
    md[, cluster := as.character(get(cluster_col))]
    sample_col <- if ("analysis_sample_id" %in% names(md)) "analysis_sample_id" else "sample_id"
    md[, sample_id_plot := as.character(get(sample_col))]
    repaired <- read_repaired_qc(data_root, dataset_id)
    ridx <- match(md$cell_id, repaired$cell_id)
    if (any(!is.na(ridx))) md$percent.mt <- repaired$percent.mt[ridx]
    cluster_qc <- md[, {
      author_column <- intersect(c("celltype", "hpca.celltype", "encode.celltype"), names(md))
      author <- if (length(author_column)) top_label(get(author_column[[1L]])) else c(label = NA, fraction = NA)
      list(
        top_author_celltype = author[["label"]],
        top_author_celltype_fraction = as.numeric(author[["fraction"]]),
        median_nCount = as.numeric(stats::median(nCount_RNA, na.rm = TRUE)),
        median_nFeature = as.numeric(stats::median(nFeature_RNA, na.rm = TRUE)),
        median_percent_mt = if ("percent.mt" %in% names(md)) as.numeric(stats::median(percent.mt, na.rm = TRUE)) else NA_real_,
        doublet_fraction = if ("scDblFinder.class" %in% names(md)) mean(scDblFinder.class == "doublet", na.rm = TRUE) else 0,
        median_doublet_score = if ("scDblFinder.score" %in% names(md)) as.numeric(stats::median(scDblFinder.score, na.rm = TRUE)) else NA_real_,
        median_S_score = if ("S.Score" %in% names(md)) as.numeric(stats::median(S.Score, na.rm = TRUE)) else NA_real_,
        median_G2M_score = if ("G2M.Score" %in% names(md)) as.numeric(stats::median(G2M.Score, na.rm = TRUE)) else NA_real_
      )
    }, by = cluster]
    diagnostic <- merge(diagnostic[, !c("top_author_celltype", "top_author_celltype_fraction",
                                         "median_nCount", "median_nFeature", "median_percent_mt",
                                         "doublet_fraction", "median_doublet_score",
                                         "median_S_score", "median_G2M_score")],
                        cluster_qc, by = "cluster", all.x = TRUE)

    emb <- Embeddings(obj, "umap")
    plot_df <- cbind(md[match(rownames(emb), cell_id)],
                     data.table::data.table(UMAP_1 = emb[, 1L], UMAP_2 = emb[, 2L]))
    plot_umap <- function(column, filename, discrete = TRUE) {
      if (!column %in% names(plot_df) || all(is.na(plot_df[[column]]))) {
        p <- placeholder_plot(paste(dataset_id, filename), paste("Metadata unavailable:", column))
      } else {
        p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[column]])) +
          geom_point(size = 0.18, alpha = 0.7) + theme_bw() +
          labs(title = paste(dataset_id, filename), color = column)
        if (!discrete) p <- p + scale_color_viridis_c()
      }
      save_plot_pair(p, file.path(out, filename), 9, 6)
    }
    plot_umap("cluster", "UMAP_by_cluster")
    plot_umap("sample_id_plot", "UMAP_by_sample")
    patient_col <- intersect(c("patient_id", "sample_id_plot"), names(plot_df))[[1L]]
    plot_umap(patient_col, "UMAP_by_patient")
    plot_umap("timepoint", "UMAP_by_timepoint")
    plot_umap("nCount_RNA", "UMAP_by_nCount", FALSE)
    plot_umap("nFeature_RNA", "UMAP_by_nFeature", FALSE)
    plot_umap("percent.mt", "UMAP_by_percent_mt", FALSE)
    doublet_col <- intersect(c("scDblFinder.score", "doublet_score"), names(plot_df))
    plot_umap(if (length(doublet_col)) doublet_col[[1L]] else "<unavailable>",
              "UMAP_by_doublet_score", FALSE)
    plot_status[[1L]] <- data.frame(dataset_id, status = "COMPLETE_FROM_CURRENT_OBJECT",
                                    object_size_gb = size_gb, message = "")
    rm(md, plot_df, emb, obj); gc()
  } else {
    source_cluster <- file.path(data_root, dataset_id, "02_clustering",
                                "umap_primary_resolution")
    source_sample <- file.path(data_root, dataset_id, "02_clustering", "umap_by_sample")
    for (ext in c("pdf", "png")) {
      file.copy(paste0(source_cluster, ".", ext), file.path(out, paste0("UMAP_by_cluster.", ext)), overwrite = TRUE)
      file.copy(paste0(source_sample, ".", ext), file.path(out, paste0("UMAP_by_sample.", ext)), overwrite = TRUE)
    }
    reason <- paste0("Current object is ", signif(size_gb, 3),
                     " GB, above the configured ", cfg$project$max_object_load_gb,
                     " GB memory guard. Aggregate audit retained; new per-lineage UMAPs are generated in step 05.")
    for (name in c("UMAP_by_patient", "UMAP_by_timepoint", "UMAP_by_nCount",
                   "UMAP_by_nFeature", "UMAP_by_percent_mt", "UMAP_by_doublet_score")) {
      save_plot_pair(placeholder_plot(paste(dataset_id, name), reason),
                     file.path(out, name), 9, 6)
    }
    plot_status[[1L]] <- data.frame(dataset_id, status = "PARTIAL_MEMORY_GUARD",
                                    object_size_gb = size_gb, message = reason)
  }

  marker_lists <- strsplit(toupper(diagnostic$top_markers %||% ""), ";", fixed = TRUE)
  epithelial <- if (dataset_id_value == "GSE154600") {
    toupper(diagnostic$top_author_celltype) == "EPI"
  } else {
    vapply(marker_lists, function(genes) {
      any(grepl("^(EPCAM|KRT[0-9]+|MSLN|WFDC2|MUC1|MUC16|TACSTD2|PAX8)$", genes))
    }, logical(1))
  }
  epithelial[is.na(epithelial)] <- FALSE
  shared <- vapply(marker_lists, function(genes) {
    any(genes %in% c("PTPRC", "CD3D", "NKG7", "LYZ", "LST1", "COL1A1",
                     "DCN", "PECAM1", "VWF"))
  }, logical(1))
  robust_outlier <- function(x) {
    x <- as.numeric(x)
    center <- stats::median(x, na.rm = TRUE)
    spread <- stats::mad(x, center = center, na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) return(rep(FALSE, length(x)))
    is.finite(x) & abs(x - center) / spread > 3
  }
  qc_suspect <- (is.finite(diagnostic$median_percent_mt) & diagnostic$median_percent_mt > 20) |
    (is.finite(diagnostic$doublet_fraction) & diagnostic$doublet_fraction > 0.1) |
    robust_outlier(diagnostic$median_nCount) | robust_outlier(diagnostic$median_nFeature)
  diagnostic[, likely_interpretation := data.table::fcase(
    qc_suspect, "likely_technical_or_qc_effect",
    dominance_label == "strong_sample_dominance" & epithelial,
    "likely_patient_specific_tumor_state",
    shared & dominance_label != "strong_sample_dominance", "likely_shared_lineage",
    default = "mixed_or_uncertain"
  )]
  diagnostic[, interpretation_basis := paste0(
    "dominance=", dominance_label,
    "; epithelial_marker_program=", epithelial,
    "; shared_lineage_marker_program=", shared,
    "; qc_outlier_flag=", qc_suspect,
    "; cell_level_qc_available=", can_load
  )]
  counts_path <- file.path(data_root, dataset_id, "02_clustering", "cluster_by_sample_counts.csv")
  cs <- data.table::fread(counts_path, showProgress = FALSE)
  cs[, `:=`(
    patient_id = if (dataset_id_value == "GSE158722") sub("_.*$", "", sample_id) else sample_id,
    timepoint = if (dataset_id_value == "GSE158722") sub("^[^_]+_", "", sample_id) else NA_character_
  )]
  patient_composition <- cs[, .(n_cells = sum(n_cells)),
                            by = .(cluster = as.character(seurat_cluster), patient_id)]
  patient_composition[, patient_fraction := n_cells / sum(n_cells), by = cluster]
  patient_dominant <- patient_composition[order(cluster, -patient_fraction), .SD[1L], by = cluster]
  patient_dominant <- patient_dominant[, .(
    cluster, dominant_patient = patient_id,
    dominant_patient_fraction = patient_fraction
  )]
  diagnostic <- merge(diagnostic, patient_dominant, by = "cluster", all.x = TRUE)

  if (dataset_id_value == "GSE158722") {
    timepoint_composition <- cs[, .(n_cells = sum(n_cells)),
                                by = .(cluster = as.character(seurat_cluster), timepoint)]
    timepoint_composition[, timepoint_fraction := n_cells / sum(n_cells), by = cluster]
    timepoint_dominant <- timepoint_composition[order(cluster, -timepoint_fraction), .SD[1L],
                                                by = cluster]
    timepoint_dominant <- timepoint_dominant[, .(
      cluster, dominant_timepoint = timepoint,
      dominant_timepoint_fraction = timepoint_fraction
    )]
    diagnostic <- merge(diagnostic, timepoint_dominant, by = "cluster", all.x = TRUE)

    # Within-patient cluster-distribution correspondence is descriptive only:
    # timepoint is retained as biology and is never supplied as a batch variable.
    longitudinal <- cs[, .(n_cells = sum(n_cells)),
                       by = .(patient_id, timepoint,
                              cluster = as.character(seurat_cluster))]
    longitudinal[, cluster_fraction := n_cells / sum(n_cells),
                 by = .(patient_id, timepoint)]
    write_csv(longitudinal, file.path(out, "patient_timepoint_cluster_distribution.csv"))
    pair_rows <- list()
    for (patient in unique(longitudinal$patient_id)) {
      patient_dt <- longitudinal[patient_id == patient]
      tp <- sort(unique(patient_dt$timepoint))
      if (length(tp) < 2L) next
      pairs <- utils::combn(tp, 2L, simplify = FALSE)
      for (pair in pairs) {
        a <- patient_dt[timepoint == pair[[1L]], .(cluster, fraction_a = cluster_fraction)]
        b <- patient_dt[timepoint == pair[[2L]], .(cluster, fraction_b = cluster_fraction)]
        joined <- merge(a, b, by = "cluster", all = TRUE)
        joined[is.na(fraction_a), fraction_a := 0]
        joined[is.na(fraction_b), fraction_b := 0]
        pair_rows[[length(pair_rows) + 1L]] <- data.frame(
          patient_id = patient,
          timepoint_a = pair[[1L]], timepoint_b = pair[[2L]],
          cluster_distribution_similarity = 1 - sum(abs(joined$fraction_a - joined$fraction_b)) / 2,
          dominant_cluster_a = joined$cluster[[which.max(joined$fraction_a)]],
          dominant_cluster_b = joined$cluster[[which.max(joined$fraction_b)]],
          same_dominant_cluster = joined$cluster[[which.max(joined$fraction_a)]] ==
            joined$cluster[[which.max(joined$fraction_b)]],
          stringsAsFactors = FALSE
        )
      }
    }
    correspondence <- data.table::rbindlist(pair_rows, fill = TRUE)
    write_csv(correspondence,
              file.path(out, "patient_timepoint_cluster_correspondence.csv"))
  } else {
    diagnostic[, `:=`(dominant_timepoint = NA_character_,
                      dominant_timepoint_fraction = NA_real_)]
  }
  write_csv(diagnostic, file.path(out, "cluster_dominance_diagnostic_table.csv"))

  cs[, fraction := n_cells / sum(n_cells), by = seurat_cluster]
  p_heat <- ggplot(cs, aes(x = sample_id, y = factor(seurat_cluster), fill = fraction)) +
    geom_tile() + scale_fill_viridis_c() + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste(dataset_id, "cluster sample fractions"), x = "Sample", y = "Cluster")
  save_plot_pair(p_heat, file.path(out, "cluster_sample_fraction_heatmap"), 12, 8)
  p_entropy <- ggplot(diagnostic, aes(x = reorder(factor(cluster), normalized_shannon_entropy),
                                     y = normalized_shannon_entropy, fill = dominance_label)) +
    geom_col() + coord_flip() + theme_bw() +
    labs(title = paste(dataset_id, "cluster sample entropy"), x = "Cluster", y = "Normalized entropy")
  save_plot_pair(p_entropy, file.path(out, "cluster_sample_entropy"), 8, 8)
  pc <- pc_audit[dataset_id == dataset_id_value]
  if (nrow(pc) && any(is.finite(pc$sample_eta_squared))) {
    p_pc <- ggplot(pc, aes(x = PC, y = sample_eta_squared)) + geom_line() + geom_point() +
      theme_bw() + labs(title = paste(dataset_id, "PC-sample association"), y = "Sample eta squared")
  } else {
    p_pc <- placeholder_plot(paste(dataset_id, "PC-sample association"),
                             unique(pc$status %||% "PC audit unavailable"))
  }
  save_plot_pair(p_pc, file.path(out, "PC_sample_association"), 8, 6)
  write_csv(data.table::rbindlist(plot_status, fill = TRUE), file.path(out, "plot_generation_status.csv"))
}

message("Sample-dominance diagnosis complete")
