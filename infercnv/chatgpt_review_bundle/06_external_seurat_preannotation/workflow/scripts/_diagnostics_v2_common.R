options(stringsAsFactors = FALSE, warn = 1)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

arg_value <- function(flag, default = NULL, args = commandArgs(trailingOnly = TRUE)) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}

has_flag <- function(flag, args = commandArgs(trailingOnly = TRUE)) flag %in% args

split_arg <- function(x) {
  z <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  z[nzchar(z)]
}

read_diagnostics_config <- function() {
  path <- normalizePath(arg_value(
    "--config",
    "infercnv/chatgpt_review_bundle/06_external_seurat_preannotation/workflow/config/diagnostics_v2.yaml"
  ), winslash = "/", mustWork = TRUE)
  cfg <- yaml::read_yaml(path)
  cfg$project$data_root <- arg_value("--data-root", cfg$project$data_root)
  cfg$project$repo_root <- arg_value("--repo-root", cfg$project$repo_root)
  list(path = path, cfg = cfg)
}

assert_datasets <- function(requested, allowed) {
  if (!length(requested) || any(!requested %in% allowed)) {
    stop("Permitted datasets for this step: ", paste(allowed, collapse = ", "))
  }
  invisible(requested)
}

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

save_plot_pair <- function(plot, stem, width = 8, height = 6) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot,
                  width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot,
                  width = width, height = height, dpi = 180, limitsize = FALSE)
}

placeholder_plot <- function(title, reason) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = reason,
                      size = 4.2, lineheight = 1.1) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title) + ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

get_counts <- function(x) {
  if (inherits(x, "SingleCellExperiment")) {
    return(SummarizedExperiment::assay(x, "counts"))
  }
  tryCatch(
    SeuratObject::LayerData(x, assay = "RNA", layer = "counts"),
    error = function(e) Seurat::GetAssayData(x, assay = "RNA", slot = "counts")
  )
}

cell_ids_from_qc <- function(md, dataset_id) {
  if (dataset_id == "GSE154600") {
    paste(md$sample_id, md$original_cell_id, sep = "__")
  } else if (dataset_id == "GSE158722") {
    original <- as.character(md$original_cell_id)
    sample_prefix <- paste0(md$sample_id, "_")
    patient_prefix <- paste0(md$patient_id, "_")
    barcode <- ifelse(
      startsWith(original, sample_prefix),
      substring(original, nchar(sample_prefix) + 1L),
      ifelse(startsWith(original, patient_prefix),
             substring(original, nchar(patient_prefix) + 1L), original)
    )
    paste(md$patient_id, md$sample_id, barcode, sep = "__")
  } else {
    paste(md$sample_id, md$original_cell_id, sep = "__")
  }
}

read_repaired_qc <- function(data_root, dataset_id) {
  repaired <- file.path(data_root, "diagnostics_v2", dataset_id,
                        "01_mt_audit", "qc_metadata_mt_repaired.csv.gz")
  original <- file.path(data_root, dataset_id, "01_qc", "qc_metadata.csv.gz")
  path <- if (file.exists(repaired)) repaired else original
  md <- data.table::fread(path, showProgress = FALSE)
  if (!"cell_id" %in% names(md)) md$cell_id <- cell_ids_from_qc(md, dataset_id)
  if (!"diagnostics_v2_qc_pass" %in% names(md)) md$diagnostics_v2_qc_pass <- TRUE
  md
}

dominance_metrics <- function(cluster, sample, dataset_id, strategy,
                              lineage = "All") {
  tab <- data.table::data.table(cluster = as.character(cluster),
                                sample = as.character(sample))[
                                  , .N, by = .(cluster, sample)]
  total_samples <- data.table::uniqueN(sample)
  out <- tab[, {
    p <- N / sum(N)
    dominant <- which.max(N)
    entropy <- if (total_samples > 1L) -sum(p * log(p)) / log(total_samples) else NA_real_
    list(n_cells = sum(N), n_samples = .N,
         dominant_sample = sample[[dominant]],
         dominant_sample_n = N[[dominant]],
         dominant_sample_fraction = max(p),
         normalized_shannon_entropy = entropy,
         simpson_diversity = 1 - sum(p^2),
         effective_sample_number = 1 / sum(p^2))
  }, by = cluster]
  out[, `:=`(dataset_id = dataset_id, lineage = lineage, strategy = strategy)]
  out[, dominance_label := data.table::fcase(
    dominant_sample_fraction >= 0.8 & n_cells >= 100,
    "strong_sample_dominance",
    dominant_sample_fraction >= 0.6, "moderate_sample_dominance",
    default = "mixed"
  )]
  data.table::setcolorder(out, c("dataset_id", "lineage", "strategy", "cluster"))
  out[]
}

cramers_v <- function(x, y) {
  tab <- table(x, y)
  if (nrow(tab) < 2L || ncol(tab) < 2L || !sum(tab)) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tab) * min(nrow(tab) - 1L, ncol(tab) - 1L))))
}

balanced_cells <- function(md, cap, seed) {
  if (nrow(md) <= cap) return(as.character(md$cell_id))
  set.seed(seed)
  groups <- split(as.character(md$cell_id), as.character(md$sample_id))
  target <- max(1L, floor(cap / length(groups)))
  chosen <- unlist(lapply(groups, function(z) sample(z, min(length(z), target))),
                   use.names = FALSE)
  remaining <- setdiff(as.character(md$cell_id), chosen)
  if (length(chosen) < cap) {
    chosen <- c(chosen, sample(remaining, min(length(remaining), cap - length(chosen))))
  }
  chosen[seq_len(min(cap, length(chosen)))]
}

broad_marker_sets <- list(
  Epithelial_like = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MSLN", "WFDC2", "MUC1"),
  Myeloid_like = c("LYZ", "LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CTSS"),
  T_NK_like = c("CD3D", "CD3E", "TRBC1", "TRBC2", "NKG7", "GNLY", "KLRD1", "IL7R"),
  B_Plasma_like = c("CD79A", "MS4A1", "CD37", "CD74", "MZB1", "JCHAIN", "SDC1", "IGHG1"),
  Fibroblast_like = c("COL1A1", "COL1A2", "DCN", "COL3A1", "COL6A1", "PDGFRA", "LUM", "C7"),
  Endothelial_like = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "PLVAP", "CLDN5"),
  Cycling_like = c("MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "STMN1", "TUBA1B", "HMGB2")
)

strategy_cluster_column <- function(strategy) paste0("cluster_", strategy)
