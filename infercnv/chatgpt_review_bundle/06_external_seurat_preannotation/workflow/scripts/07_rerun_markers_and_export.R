#!/usr/bin/env Rscript

required <- c("yaml", "data.table", "Matrix", "Seurat", "ggplot2")
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
assert_datasets(datasets, c("GSE154600", "GSE158722"))
recommendations <- data.table::fread(
  file.path(data_root, "diagnostics_v2", "strategy_comparison",
            "recommended_strategy_by_dataset_and_lineage.csv"), showProgress = FALSE)

run_markers <- function(obj, cluster, prefix) {
  Idents(obj) <- factor(as.character(cluster))
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  markers <- FindAllMarkers(
    obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
    min.pct = 0.20, logfc.threshold = 0.25, return.thresh = 0.05,
    verbose = FALSE
  )
  if (!"gene" %in% names(markers)) markers$gene <- rownames(markers)
  markers <- data.table::as.data.table(markers)
  if (nrow(markers)) markers[, cluster := paste(prefix, as.character(cluster), sep = "__")]
  markers
}

for (dataset_id in datasets) {
  message("===== ", dataset_id, " =====")
  dataset_id_value <- dataset_id
  ds_rec <- recommendations[dataset_id == dataset_id_value]
  out <- file.path(data_root, "diagnostics_v2", dataset_id, "05_markers")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  final_markers <- list(); final_averages <- list(); dot_rows <- list()
  cluster_rows <- list(); overlap_rows <- list(); sensitivity_rows <- list()

  for (i in seq_len(nrow(ds_rec))) {
    lineage <- ds_rec$lineage[[i]]
    selected <- ds_rec$recommended_strategy[[i]]
    input_path <- file.path(data_root, "diagnostics_v2", "objects", dataset_id,
                            "lineage_inputs", paste0(lineage, "_strategy_input.rds"))
    if (!file.exists(input_path)) next
    payload <- readRDS(input_path)
    md <- payload$metadata
    rownames(md) <- md$cell_id
    obj <- CreateSeuratObject(payload$counts, project = dataset_id,
                              min.cells = 0, min.features = 0, meta.data = md)
    obj <- NormalizeData(obj, verbose = FALSE)
    rm(payload); gc()

    strategies <- unique(c("A_uncorrected", selected))
    marker_sets <- list()
    for (strategy in strategies) {
      embedding_path <- file.path(data_root, "diagnostics_v2", dataset_id,
                                  "04_strategies", lineage, strategy,
                                  "cell_embedding_and_clusters.csv.gz")
      if (!file.exists(embedding_path)) next
      assignment <- data.table::fread(embedding_path, showProgress = FALSE)
      idx <- match(colnames(obj), assignment$cell_id)
      if (anyNA(idx)) stop(dataset_id, " / ", lineage, " / ", strategy,
                           ": cluster assignments do not match RNA object")
      prefix <- paste(lineage, strategy, sep = "__")
      markers <- run_markers(obj, assignment$cluster[idx], prefix)
      marker_sets[[strategy]] <- unique(markers[order(p_val_adj, -avg_log2FC),
                                                head(gene, 50L), by = cluster]$V1)
      sens_dir <- file.path(out, "strategy_marker_sensitivity", lineage, strategy)
      write_csv_gz(markers, file.path(sens_dir, "all_cluster_markers.csv.gz"))
      top20 <- markers[order(cluster, -avg_log2FC, p_val_adj), head(.SD, 20L), by = cluster]
      write_csv(top20, file.path(sens_dir, "top20_markers_per_cluster.csv"))
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
        dataset_id, lineage, strategy, n_marker_rows = nrow(markers),
        assay = "RNA", test_use = "wilcox", min_pct = 0.20,
        logfc_threshold = 0.25, return_threshold = 0.05
      )
      if (strategy == selected) {
        final_markers[[lineage]] <- markers
        Idents(obj) <- factor(as.character(assignment$cluster[idx]))
        avg <- AverageExpression(obj, assays = "RNA", layer = "data", verbose = FALSE)$RNA
        colnames(avg) <- paste(lineage, selected, colnames(avg), sep = "__")
        final_averages[[lineage]] <- avg
        cluster_counts <- as.data.frame(table(as.character(Idents(obj))))
        names(cluster_counts) <- c("raw_cluster", "n_cells")
        cluster_counts$dataset_id <- dataset_id
        cluster_counts$lineage <- lineage
        cluster_counts$strategy <- selected
        cluster_counts$cluster <- paste(lineage, selected, cluster_counts$raw_cluster, sep = "__")
        cluster_rows[[lineage]] <- cluster_counts
        panel <- unique(unlist(broad_marker_sets, use.names = FALSE))
        panel <- panel[panel %in% rownames(obj)]
        if (length(panel)) {
          dp <- DotPlot(obj, features = panel)$data
          dp$cluster <- paste(lineage, selected, as.character(dp$id), sep = "__")
          dp$lineage <- lineage
          dot_rows[[lineage]] <- dp
        }
      }
    }
    baseline_genes <- marker_sets$A_uncorrected %||% character()
    selected_genes <- marker_sets[[selected]] %||% character()
    overlap_rows[[lineage]] <- data.frame(
      dataset_id, lineage, selected_strategy = selected,
      n_baseline_top50_union = length(baseline_genes),
      n_selected_top50_union = length(selected_genes),
      n_intersection = length(intersect(baseline_genes, selected_genes)),
      top50_union_jaccard = if (length(union(baseline_genes, selected_genes)))
        length(intersect(baseline_genes, selected_genes)) /
          length(union(baseline_genes, selected_genes)) else NA_real_
    )
    rm(obj); gc()
  }

  markers <- data.table::rbindlist(final_markers, fill = TRUE)
  if (!nrow(markers)) stop(dataset_id, ": no final RNA marker results")
  write_csv_gz(markers, file.path(out, "all_cluster_markers.csv.gz"))
  for (n in c(20L, 50L, 100L)) {
    top <- markers[order(cluster, -avg_log2FC, p_val_adj), head(.SD, n), by = cluster]
    write_csv(top, file.path(out, paste0("top", n, "_markers_per_cluster.csv")))
  }

  all_genes <- Reduce(union, lapply(final_averages, rownames))
  avg_parts <- lapply(final_averages, function(x) {
    out_mat <- matrix(NA_real_, nrow = length(all_genes), ncol = ncol(x),
                      dimnames = list(all_genes, colnames(x)))
    out_mat[rownames(x), ] <- as.matrix(x)
    out_mat
  })
  avg <- do.call(cbind, avg_parts)
  avg_export <- data.frame(gene = rownames(avg), avg, check.names = FALSE)
  write_csv_gz(avg_export, file.path(out, "cluster_average_expression.csv.gz"))

  dots <- data.table::rbindlist(dot_rows, fill = TRUE)
  if (nrow(dots)) {
    p <- ggplot(dots, aes(x = features.plot, y = cluster, size = pct.exp,
                          color = avg.exp.scaled)) +
      geom_point() + scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = paste(dataset_id, "broad RNA marker dotplot"), x = "Gene", y = "Final cluster")
    save_plot_pair(p, file.path(out, "broad_marker_dotplot"), 16, 10)
  }

  top_genes <- unique(markers[order(cluster, -avg_log2FC, p_val_adj), head(gene, 10L), by = cluster]$V1)
  top_genes <- intersect(top_genes, rownames(avg))
  heat <- data.table::as.data.table(as.table(avg[top_genes, , drop = FALSE]))
  names(heat) <- c("gene", "cluster", "average_expression")
  heat[, scaled_expression := as.numeric(scale(average_expression)), by = gene]
  p_heat <- ggplot(heat, aes(x = cluster, y = gene, fill = scaled_expression)) +
    geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste(dataset_id, "top RNA marker heatmap"), x = "Final cluster", y = "Gene")
  save_plot_pair(p_heat, file.path(out, "top_marker_heatmap"), 16, 12)

  clusters <- data.table::rbindlist(cluster_rows, fill = TRUE)
  top_text <- markers[order(cluster, -avg_log2FC, p_val_adj),
                      .(top_markers = paste(head(gene, 20L), collapse = ";")), by = cluster]
  manual <- merge(clusters[, .(dataset_id, lineage, strategy, cluster, n_cells)],
                  top_text, by = "cluster", all.x = TRUE)
  manual[, `:=`(cell_type_manual = "", cell_subtype_manual = "",
                confidence = "", notes = "")]
  data.table::setcolorder(manual, c("dataset_id", "lineage", "strategy", "cluster",
                                    "n_cells", "top_markers", "cell_type_manual",
                                    "cell_subtype_manual", "confidence", "notes"))
  write_csv(manual, file.path(out, "manual_annotation_template.csv"))
  write_csv(data.table::rbindlist(overlap_rows, fill = TRUE),
            file.path(out, "top_marker_overlap_by_lineage.csv"))
  write_csv(data.table::rbindlist(sensitivity_rows, fill = TRUE),
            file.path(out, "marker_run_audit.csv"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
}

message("RNA-assay marker rerun complete; manual annotation fields remain blank")
