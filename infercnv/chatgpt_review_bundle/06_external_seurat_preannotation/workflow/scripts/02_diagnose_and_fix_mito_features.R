#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}

has_flag <- function(flag) flag %in% args

split_arg <- function(x) {
  z <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  z[nzchar(z)]
}

required <- c("yaml", "data.table", "Matrix", "Seurat", "SingleCellExperiment",
              "SummarizedExperiment", "ggplot2")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                           FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

config_path <- normalizePath(
  arg_value("--config",
            "infercnv/chatgpt_review_bundle/06_external_seurat_preannotation/workflow/config/diagnostics_v2.yaml"),
  winslash = "/", mustWork = TRUE
)
cfg <- yaml::read_yaml(config_path)
cfg$project$data_root <- arg_value("--data-root", cfg$project$data_root)
cfg$project$repo_root <- arg_value("--repo-root", cfg$project$repo_root)
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
prepared_root <- normalizePath(cfg$project$prepared_input_root,
                               winslash = "/", mustWork = TRUE)
raw_root <- normalizePath(cfg$project$raw_count_root,
                          winslash = "/", mustWork = TRUE)
datasets <- split_arg(arg_value("--datasets", "GSE147082,GSE158722"))
allowed <- c("GSE147082", "GSE158722")
if (!length(datasets) || any(!datasets %in% allowed)) {
  stop("This script only permits: ", paste(allowed, collapse = ", "))
}

write_csv <- function(x, path) data.table::fwrite(x, path, na = "NA")

write_csv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, na = "NA")
}

pattern_definitions <- data.frame(
  pattern_name = c("uppercase_MT_hyphen", "lowercase_mt_hyphen",
                   "mixedcase_Mt_hyphen", "uppercase_MT_dot",
                   "ensembl_gene", "ensembl_transcript", "chrM_text",
                   "mitochond_text", "ribosomal_RPS_RPL", "hemoglobin"),
  regex = c("^MT-", "^mt-", "^Mt-", "^MT\\.", "^ENSG", "^ENST",
            "chrM", "mitochond", "^RP[SL]", "^HB[ABDEGQZ]"),
  ignore_case = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                  TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

pattern_audit <- function(features, source = "prepared_matrix") {
  out <- pattern_definitions
  out$source <- source
  out$n_features <- length(features)
  out$n_matches <- vapply(seq_len(nrow(out)), function(i) {
    sum(grepl(out$regex[[i]], features,
              ignore.case = out$ignore_case[[i]], perl = TRUE))
  }, numeric(1))
  out[, c("source", "pattern_name", "regex", "ignore_case",
          "n_features", "n_matches")]
}

summarize_source_file <- function(path, dataset_id) {
  header <- tryCatch(names(data.table::fread(path, nrows = 0L,
                                              showProgress = FALSE)),
                     error = identity)
  if (inherits(header, "error")) {
    return(data.frame(dataset_id = dataset_id, source_file = path,
                      status = "header_read_failed", feature_column = NA_character_,
                      n_features = NA_integer_, n_match_uppercase_MT = NA_integer_,
                      n_match_dot_MT = NA_integer_, n_match_lowercase_mt = NA_integer_,
                      n_match_ensembl = NA_integer_, message = conditionMessage(header)))
  }
  gene_column <- grep("^(Gene[ _]?Symbol|gene|genes)$", header,
                      ignore.case = TRUE)
  if (!length(gene_column)) gene_column <- 1L
  gene_column <- gene_column[[1L]]
  tab <- tryCatch(data.table::fread(path, select = gene_column,
                                    showProgress = FALSE), error = identity)
  if (inherits(tab, "error")) {
    return(data.frame(dataset_id = dataset_id, source_file = path,
                      status = "feature_read_failed",
                      feature_column = header[[gene_column]], n_features = NA_integer_,
                      n_match_uppercase_MT = NA_integer_, n_match_dot_MT = NA_integer_,
                      n_match_lowercase_mt = NA_integer_, n_match_ensembl = NA_integer_,
                      message = conditionMessage(tab)))
  }
  features <- as.character(tab[[1L]])
  data.frame(
    dataset_id = dataset_id,
    source_file = normalizePath(path, winslash = "/", mustWork = FALSE),
    status = "complete",
    feature_column = header[[gene_column]],
    n_features = length(features),
    n_match_uppercase_MT = sum(grepl("^MT-", features)),
    n_match_dot_MT = sum(grepl("^MT\\.", features)),
    n_match_lowercase_mt = sum(grepl("^mt-", features)),
    n_match_ensembl = sum(grepl("^ENSG", features)),
    message = ""
  )
}

feature_metadata_audit <- function(row_data, n_features) {
  if (is.null(row_data) || !ncol(row_data)) {
    return(data.frame(column_name = "<none>", column_class = NA_character_,
                      n_rows = n_features, n_nonmissing = 0L,
                      example_values = ""))
  }
  data.table::rbindlist(lapply(names(row_data), function(column) {
    x <- row_data[[column]]
    examples <- unique(as.character(x[!is.na(x)]))
    data.frame(column_name = column, column_class = class(x)[[1L]],
               n_rows = length(x), n_nonmissing = sum(!is.na(x)),
               example_values = paste(head(examples, 10L), collapse = ";"))
  }), fill = TRUE)
}

detect_features <- function(feature_ids, row_data) {
  result <- list(mt = character(), mt_symbols = character(),
                 mt_source = NA_character_, status = "unavailable")
  if (!is.null(row_data) && ncol(row_data)) {
    chromosome_columns <- grep("^(chromosome|chr|seqnames)$", names(row_data),
                               ignore.case = TRUE, value = TRUE)
    for (column in chromosome_columns) {
      idx <- as.character(row_data[[column]]) %in% c("MT", "chrM", "M", "chrMT")
      if (any(idx, na.rm = TRUE)) {
        result$mt <- feature_ids[idx]
        result$mt_symbols <- feature_ids[idx]
        result$mt_source <- paste0("feature_metadata_chromosome:", column)
        result$status <- "available"
        break
      }
    }
    if (!length(result$mt)) {
      symbol_columns <- grep("^(gene_?symbol|symbol|gene_?name)$", names(row_data),
                             ignore.case = TRUE, value = TRUE)
      for (column in symbol_columns) {
        symbols <- as.character(row_data[[column]])
        idx <- grepl("^MT-", symbols, ignore.case = TRUE)
        if (any(idx, na.rm = TRUE)) {
          result$mt <- feature_ids[idx]
          result$mt_symbols <- symbols[idx]
          result$mt_source <- paste0("feature_metadata_symbol:", column)
          result$status <- "available"
          break
        }
      }
    }
  }
  if (!length(result$mt)) {
    idx <- grepl("^MT-", feature_ids, ignore.case = TRUE)
    if (any(idx)) {
      result$mt <- feature_ids[idx]
      result$mt_symbols <- toupper(feature_ids[idx])
      result$mt_source <- "rownames_gene_symbol_hyphen"
      result$status <- "available"
    }
  }
  if (!length(result$mt)) {
    # R's make.names() converts canonical MT-* symbols to MT.*. This is the
    # observed GSE147082 failure mode; matching is restricted to the MT. prefix.
    idx <- grepl("^MT\\.", feature_ids, ignore.case = TRUE)
    if (any(idx)) {
      result$mt <- feature_ids[idx]
      result$mt_symbols <- sub("^MT\\.", "MT-", toupper(feature_ids[idx]))
      result$mt_source <- "rownames_gene_symbol_dot_separator"
      result$status <- "available"
    }
  }
  result$ribo <- feature_ids[grepl("^RP[SL]", feature_ids, ignore.case = TRUE)]
  result$hb <- feature_ids[grepl("^HB[ABDEGQZ]", feature_ids, ignore.case = TRUE)]
  result
}

metric_summary <- function(dataset_id, metric, before, after,
                           available, n_features_used, source) {
  stats <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(rep(NA_real_, 6L))
    c(min(x), stats::quantile(x, 0.25, names = FALSE), stats::median(x),
      mean(x), stats::quantile(x, 0.75, names = FALSE), max(x))
  }
  a <- stats(before)
  b <- stats(after)
  data.frame(
    dataset_id = dataset_id, metric = metric,
    n_cells = length(after), available = available,
    n_features_used = n_features_used, feature_source = source,
    before_min = a[[1L]], before_q25 = a[[2L]], before_median = a[[3L]],
    before_mean = a[[4L]], before_q75 = a[[5L]], before_max = a[[6L]],
    after_min = b[[1L]], after_q25 = b[[2L]], after_median = b[[3L]],
    after_mean = b[[4L]], after_q75 = b[[5L]], after_max = b[[6L]]
  )
}

save_plot_pair <- function(plot, stem, width, height) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot,
                  width = width, height = height)
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot,
                  width = width, height = height, dpi = 200)
}

mad_upper <- function(x, multiplier) {
  med <- stats::median(x, na.rm = TRUE)
  md <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
  med + multiplier * md
}

as_dgc <- function(x) {
  methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
}

table_columns_to_dgc <- function(tab, count_columns, block_size = 250L) {
  groups <- split(count_columns,
                  ceiling(seq_along(count_columns) / block_size))
  blocks <- lapply(groups, function(columns) {
    values <- as.matrix(tab[, columns, with = FALSE])
    out <- as_dgc(values)
    rm(values); gc()
    out
  })
  do.call(cbind, blocks)
}

compute_gse158722_raw_qc <- function(raw_files, out, resume = FALSE) {
  part_dir <- file.path(out, "raw_qc_parts")
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  metric_parts <- list()
  feature_parts <- list()
  for (path in raw_files) {
    patient_id <- sub("^GSE158722_", "",
                      sub("\\.counts\\.txt\\.gz$", "", basename(path)))
    metric_path <- file.path(part_dir, paste0(patient_id, "_qc_metrics.csv.gz"))
    feature_path <- file.path(part_dir, paste0(patient_id, "_features_used.csv"))
    if (resume && file.exists(metric_path) && file.exists(feature_path)) {
      message("Resuming raw QC part: ", patient_id)
      metrics <- data.table::fread(metric_path, showProgress = FALSE)
      if (!"raw_original_cell_id" %in% names(metrics)) {
        barcode <- sub("^.*__", "", metrics$cell_id)
        patient_barcode <- paste(metrics$patient_id, barcode, sep = "_")
        metrics$raw_original_cell_id <- ifelse(
          metrics$sample_id == patient_barcode,
          metrics$sample_id,
          paste(metrics$sample_id, barcode, sep = "_")
        )
      }
      metrics$raw_match_key <- paste(metrics$patient_id,
                                     metrics$raw_original_cell_id, sep = "__")
      metric_parts[[patient_id]] <- metrics
      feature_parts[[patient_id]] <- data.table::fread(feature_path,
                                                        showProgress = FALSE)
      next
    }
    message("Computing raw-source QC metrics: ", patient_id, " from ", path)
    tab <- data.table::fread(path, check.names = FALSE, showProgress = FALSE)
    gene_column <- grep("^Gene[ _]Symbol$", names(tab), ignore.case = TRUE)
    if (length(gene_column) != 1L) {
      stop(patient_id, ": cannot identify a unique Gene Symbol column")
    }
    raw_symbols <- as.character(tab[[gene_column]])
    feature_ids <- raw_symbols
    blank_symbol <- is.na(feature_ids) | !nzchar(trimws(feature_ids))
    feature_ids[blank_symbol] <- paste0("UNMAPPEDFEATURE", which(blank_symbol))
    feature_ids <- make.unique(feature_ids)
    annotation_columns <- grep("^(ENSEMBL[ _]ID|Gene[ _]ID|Gene[ _]Symbol)$",
                               names(tab), ignore.case = TRUE)
    count_columns <- setdiff(seq_len(ncol(tab)), annotation_columns)
    numeric_counts <- vapply(tab[, count_columns, with = FALSE],
                             is.numeric, logical(1))
    if (!all(numeric_counts)) {
      stop(patient_id, ": non-numeric count columns: ",
           paste(names(numeric_counts)[!numeric_counts], collapse = ", "))
    }
    original <- names(tab)[count_columns]
    parts <- strsplit(original, "_", fixed = TRUE)
    timepoint <- vapply(parts, function(x) {
      if (length(x) >= 2L) x[[2L]] else NA_character_
    }, character(1))
    sample_id <- paste(patient_id, timepoint, sep = "_")
    barcode <- vapply(parts, function(x) {
      if (length(x) >= 3L) paste(x[-(1:2)], collapse = "_") else tail(x, 1L)
    }, character(1))
    cell_id <- paste(patient_id, sample_id, barcode, sep = "__")
    mt <- feature_ids[grepl("^MT-", raw_symbols, ignore.case = TRUE)]
    ribo <- feature_ids[grepl("^RP[SL]", raw_symbols, ignore.case = TRUE)]
    hb <- feature_ids[grepl("^HB[ABDEGQZ]", raw_symbols, ignore.case = TRUE)]
    column_blocks <- split(seq_along(count_columns),
                           ceiling(seq_along(count_columns) / 200L))
    metric_blocks <- lapply(column_blocks, function(block_index) {
      columns <- count_columns[block_index]
      values <- as.matrix(tab[, columns, with = FALSE])
      mat <- as_dgc(values)
      rm(values); gc()
      rownames(mat) <- feature_ids
      colnames(mat) <- cell_id[block_index]
      obj <- CreateSeuratObject(mat, project = patient_id,
                                min.cells = 0, min.features = 0)
      if (length(mt)) {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt)
      } else {
        obj[["percent.mt"]] <- rep(NA_real_, ncol(obj))
      }
      obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo)
      obj[["percent.HB"]] <- PercentageFeatureSet(obj, features = hb)
      result <- data.frame(
        cell_id = colnames(obj), patient_id = patient_id,
        sample_id = sample_id[block_index], timepoint = timepoint[block_index],
        raw_original_cell_id = original[block_index],
        raw_match_key = paste(patient_id, original[block_index], sep = "__"),
        percent_mt_available = length(mt) > 0L,
        percent.mt = as.numeric(obj$percent.mt),
        percent.ribo = as.numeric(obj$percent.ribo),
        percent.HB = as.numeric(obj$percent.HB),
        source_file = normalizePath(path, winslash = "/", mustWork = FALSE)
      )
      rm(mat, obj); gc()
      result
    })
    metrics <- data.table::rbindlist(metric_blocks, fill = TRUE)
    feature_table <- data.table::rbindlist(list(
      data.frame(feature_type = rep("mt", length(mt)), feature_id = mt,
                 symbol = raw_symbols[match(mt, feature_ids)]),
      data.frame(feature_type = rep("ribo", length(ribo)), feature_id = ribo,
                 symbol = raw_symbols[match(ribo, feature_ids)]),
      data.frame(feature_type = rep("hb", length(hb)), feature_id = hb,
                 symbol = raw_symbols[match(hb, feature_ids)])
    ), fill = TRUE)
    feature_table$source_file <- normalizePath(path, winslash = "/",
                                                mustWork = FALSE)
    write_csv_gz(metrics, metric_path)
    write_csv(feature_table, feature_path)
    metric_parts[[patient_id]] <- metrics
    feature_parts[[patient_id]] <- feature_table
    rm(tab, metric_blocks, metrics, feature_table); gc()
  }
  list(
    metrics = data.table::rbindlist(metric_parts, fill = TRUE),
    features = data.table::rbindlist(feature_parts, fill = TRUE)
  )
}

for (dataset_id in datasets) {
  message("===== ", dataset_id, " =====")
  out <- file.path(data_root, "diagnostics_v2", dataset_id, "01_mt_audit")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ds_cfg <- cfg$datasets[[dataset_id]]
  prepared_path <- file.path(prepared_root, ds_cfg$prepared_sce)

  sce <- NULL
  if (dataset_id == "GSE147082") {
    message("Reading prepared SCE: ", prepared_path)
    sce <- readRDS(prepared_path)
    features <- rownames(sce)
    row_data <- as.data.frame(SummarizedExperiment::rowData(sce))
  } else {
    common_path <- file.path(prepared_root, ds_cfg$common_gene_file)
    features <- as.character(readRDS(common_path))
    row_data <- data.frame(row.names = features)
  }

  writeLines(head(features, 100L), file.path(out, "feature_name_head_100.txt"),
             useBytes = TRUE)
  writeLines(tail(features, 100L), file.path(out, "feature_name_tail_100.txt"),
             useBytes = TRUE)
  write_csv(pattern_audit(features),
            file.path(out, "feature_name_pattern_audit.csv"))
  write_csv(feature_metadata_audit(row_data, length(features)),
            file.path(out, "feature_metadata_columns.csv"))

  raw_dir <- file.path(raw_root, dataset_id)
  if (dataset_id == "GSE147082") raw_dir <- file.path(raw_dir, "extracted")
  raw_pattern <- if (dataset_id == "GSE147082") "\\.csv\\.gz$" else "\\.counts\\.txt\\.gz$"
  raw_files <- sort(list.files(raw_dir, pattern = raw_pattern, full.names = TRUE))
  raw_audit_path <- file.path(out, "source_file_feature_pattern_audit.csv")
  raw_audit <- if (has_flag("--resume") && file.exists(raw_audit_path)) {
    message("Resuming from completed raw feature audit: ", raw_audit_path)
    data.table::fread(raw_audit_path, showProgress = FALSE)
  } else if (length(raw_files)) {
    data.table::rbindlist(lapply(raw_files, summarize_source_file,
                                dataset_id = dataset_id), fill = TRUE)
  } else {
    data.frame(dataset_id = dataset_id, source_file = raw_dir,
               status = "no_raw_source_files_found", feature_column = NA_character_,
               n_features = NA_integer_, n_match_uppercase_MT = NA_integer_,
               n_match_dot_MT = NA_integer_, n_match_lowercase_mt = NA_integer_,
               n_match_ensembl = NA_integer_, message = "")
  }
  write_csv(raw_audit, raw_audit_path)

  detected <- detect_features(features, row_data)
  raw_mt_available <- dataset_id == "GSE158722" &&
    any(raw_audit$status == "complete" &
          as.numeric(raw_audit$n_match_uppercase_MT) > 0L, na.rm = TRUE)
  mt_available <- length(detected$mt) > 0L || raw_mt_available
  mt_source <- if (length(detected$mt)) {
    detected$mt_source
  } else if (raw_mt_available) {
    "original_raw_gene_symbol_per_patient"
  } else {
    NA_character_
  }

  feature_id_type <- if (mean(grepl("^ENSG", features)) > 0.8) {
    "ensembl_gene_id"
  } else {
    "gene_symbol_or_symbol_like"
  }
  old_md_path <- file.path(data_root, dataset_id, "01_qc", "qc_metadata.csv.gz")
  md <- data.table::fread(old_md_path, showProgress = FALSE)
  md$cell_id <- if (dataset_id == "GSE147082") {
    paste(md$sample_id, md$original_cell_id, sep = "__")
  } else {
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
  }
  before_mt <- as.numeric(md$percent.mt)
  before_ribo <- as.numeric(md$percent.ribo)
  before_hb <- as.numeric(md$percent.HB)

  if (length(detected$mt)) {
    counts <- SummarizedExperiment::assay(sce, "counts")
    obj <- CreateSeuratObject(counts = counts, project = dataset_id,
                              min.cells = 0, min.features = 0)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = detected$mt)
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = detected$ribo)
    obj[["percent.HB"]] <- PercentageFeatureSet(obj, features = detected$hb)
    idx <- match(md$cell_id, colnames(obj))
    if (anyNA(idx)) stop(dataset_id, ": failed to match ", sum(is.na(idx)),
                         " existing QC cells to prepared counts")
    md$percent.mt <- as.numeric(obj$percent.mt[idx])
    md$percent.ribo <- as.numeric(obj$percent.ribo[idx])
    md$percent.HB <- as.numeric(obj$percent.HB[idx])
    mt_table <- data.frame(feature_id = detected$mt,
                           symbol = detected$mt_symbols,
                           source = rep(detected$mt_source, length(detected$mt)))
    ribo_table <- data.frame(feature_id = detected$ribo,
                             symbol = detected$ribo,
                             source = "rownames_case_insensitive_RPS_RPL")
    hb_table <- data.frame(feature_id = detected$hb,
                           symbol = detected$hb,
                           source = "rownames_case_insensitive_HB")
    rm(obj, counts, sce); gc()
  } else if (raw_mt_available) {
    raw_qc <- compute_gse158722_raw_qc(raw_files, out,
                                       resume = has_flag("--resume"))
    md$raw_match_key <- paste(md$patient_id, md$original_cell_id, sep = "__")
    idx <- match(md$raw_match_key, raw_qc$metrics$raw_match_key)
    if (anyNA(idx)) {
      stop(dataset_id, ": failed to match ", sum(is.na(idx)),
           " existing QC cells to original raw-count QC metrics")
    }
    md$percent.mt <- as.numeric(raw_qc$metrics$percent.mt[idx])
    md$percent.ribo <- as.numeric(raw_qc$metrics$percent.ribo[idx])
    md$percent.HB <- as.numeric(raw_qc$metrics$percent.HB[idx])
    mt_table <- raw_qc$features[feature_type == "mt",
                                .(feature_id, symbol, source = source_file)]
    ribo_table <- raw_qc$features[feature_type == "ribo",
                                  .(feature_id, symbol, source = source_file)]
    hb_table <- raw_qc$features[feature_type == "hb",
                                .(feature_id, symbol, source = source_file)]
    raw_availability <- unique(raw_qc$metrics[
      , .(patient_id, source_file, percent_mt_available)
    ])
    write_csv(raw_availability,
              file.path(out, "raw_qc_patient_mt_availability.csv"))
    rm(raw_qc); gc()
  } else {
    md$percent.mt <- NA_real_
    # The ribosomal and hemoglobin patterns are valid for this matrix. The
    # existing values are retained to avoid loading the 3.8 GB prepared SCE.
    md$percent.ribo <- before_ribo
    md$percent.HB <- before_hb
    mt_table <- data.frame(feature_id = character(), symbol = character(),
                           source = character())
    ribo_table <- data.frame(feature_id = detected$ribo,
                             symbol = detected$ribo,
                             source = "rownames_case_insensitive_RPS_RPL")
    hb_table <- data.frame(feature_id = detected$hb,
                           symbol = detected$hb,
                           source = "rownames_case_insensitive_HB")
  }

  write_csv(mt_table, file.path(out, "mt_features_used.csv"))
  write_csv(ribo_table, file.path(out, "ribosomal_features_used.csv"))
  write_csv(hb_table, file.path(out, "hemoglobin_features_used.csv"))
  n_mt_features <- data.table::uniqueN(mt_table$feature_id)
  n_ribo_features <- data.table::uniqueN(ribo_table$feature_id)
  n_hb_features <- data.table::uniqueN(hb_table$feature_id)
  n_cells_mt_available <- sum(is.finite(md$percent.mt))
  fraction_cells_mt_available <- mean(is.finite(md$percent.mt))

  detection <- data.frame(
    dataset_id = dataset_id,
    n_features = length(features),
    feature_id_type = feature_id_type,
    n_match_uppercase_MT = sum(grepl("^MT-", features)),
    n_match_lowercase_mt = sum(grepl("^mt-", features)),
    n_match_dot_MT = sum(grepl("^MT\\.", features)),
    n_match_ensembl = sum(grepl("^ENSG", features)),
    n_feature_metadata_rows = nrow(row_data),
    n_feature_metadata_columns = ncol(row_data),
    candidate_mt_source = mt_source,
    n_mt_features_used = n_mt_features,
    percent_mt_available = mt_available,
    n_cells_percent_mt_available = n_cells_mt_available,
    fraction_cells_percent_mt_available = fraction_cells_mt_available,
    status = if (length(detected$mt)) {
      "available_recomputed_explicit_features"
    } else if (raw_mt_available && fraction_cells_mt_available < 1) {
      "partially_available_from_original_raw_files_some_patients_lack_mt_features"
    } else if (raw_mt_available) {
      "available_recomputed_from_original_raw_files_prepared_common_matrix_lacks_mt"
    } else {
      "unavailable_no_credible_mt_features_in_prepared_or_raw_source"
    }
  )
  write_csv(detection, file.path(out, "mt_feature_detection.csv"))

  threshold_rows <- list()
  md$mt_qc_pass <- TRUE
  for (sample_id in unique(as.character(md$analysis_sample_id))) {
    idx <- which(md$analysis_sample_id == sample_id)
    sample_mt_available <- any(is.finite(md$percent.mt[idx]))
    if (sample_mt_available) {
      max_mt <- min(as.numeric(cfg$qc$max_percent_mt),
                    mad_upper(md$percent.mt[idx], as.numeric(cfg$qc$mad_upper)))
      md$mt_qc_pass[idx] <- !is.finite(md$percent.mt[idx]) |
        md$percent.mt[idx] <= max_mt
    } else {
      max_mt <- NA_real_
      md$mt_qc_pass[idx] <- TRUE
    }
    threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, sample_id = sample_id,
      percent_mt_available = sample_mt_available,
      max_percent_mt = max_mt, n_before = length(idx),
      n_after_mt_qc = sum(md$mt_qc_pass[idx])
    )
  }
  md$diagnostics_v2_qc_pass <- md$mt_qc_pass
  write_csv(data.table::rbindlist(threshold_rows),
            file.path(out, "qc_thresholds_repaired.csv"))
  write_csv(data.frame(
    dataset_id = dataset_id, n_current_singlets = nrow(md),
    n_after_repaired_mt_qc = sum(md$diagnostics_v2_qc_pass),
    retained_fraction = mean(md$diagnostics_v2_qc_pass),
    percent_mt_available = mt_available,
    fraction_cells_percent_mt_available = fraction_cells_mt_available,
    decision = if (mt_available) "apply_explicit_feature_mt_threshold" else
      "disable_mt_filter_public_matrix_lacks_mt_features"
  ), file.path(out, "qc_cell_retention_repaired.csv"))
  write_csv_gz(md, file.path(out, "qc_metadata_mt_repaired.csv.gz"))

  summaries <- data.table::rbindlist(list(
    metric_summary(dataset_id, "percent.mt", before_mt, md$percent.mt,
                   mt_available, n_mt_features, mt_source),
    metric_summary(dataset_id, "percent.ribo", before_ribo, md$percent.ribo,
                    TRUE, n_ribo_features,
                    if (mt_available) "explicit_features_recomputed" else
                     "feature_pattern_verified_existing_metric_retained"),
    metric_summary(dataset_id, "percent.HB", before_hb, md$percent.HB,
                    TRUE, n_hb_features,
                    if (mt_available) "explicit_features_recomputed" else
                     "feature_pattern_verified_existing_metric_retained")
  ), fill = TRUE)
  write_csv(summaries, file.path(out, "qc_metric_before_after.csv"))

  if (mt_available) {
    plot_df <- md[is.finite(percent.mt)]
    p1 <- ggplot(plot_df, aes(x = percent.mt)) +
      geom_histogram(bins = 80, fill = "#2C7FB8", color = "white") +
      theme_bw() + labs(title = paste(dataset_id, "repaired percent.mt"),
                        x = "percent.mt", y = "Cells")
    p2 <- ggplot(plot_df, aes(x = analysis_sample_id, y = percent.mt)) +
      geom_violin(fill = "#7FCDBB", scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.size = 0.2) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(dataset_id, "percent.mt by sample"),
           x = "Sample", y = "percent.mt")
  } else {
    note <- paste(dataset_id,
                  "percent.mt unavailable:\nprepared public matrix contains no credible mitochondrial features")
    p1 <- ggplot() + annotate("text", x = 0, y = 0, label = note, size = 5) +
      xlim(-1, 1) + ylim(-1, 1) + theme_void()
    p2 <- p1
  }
  save_plot_pair(p1, file.path(out, "percent_mt_distribution"), 8, 6)
  save_plot_pair(p2, file.path(out, "percent_mt_by_sample"), 10, 6)

  source_note <- if (length(detected$mt)) {
    paste0("Explicit features from `", mt_source, "` were passed to ",
           "`PercentageFeatureSet(features = ...)`. The repaired metric is available.")
  } else if (raw_mt_available) {
    paste0("Per-patient original raw count files were loaded independently, and ",
           "their explicit `^MT-` features were passed to ",
           "`PercentageFeatureSet(features = ...)`. Only the resulting QC metadata ",
           "was joined by cell ID. The prepared common-gene expression matrix was ",
           "not changed and no features were added. Patients whose raw source lacks ",
           "credible mitochondrial features retain `percent.mt = NA` and do not use ",
           "an mt threshold; availability is recorded per patient.")
  } else {
    paste0("No credible mitochondrial features exist in the prepared matrix. ",
           "`percent.mt` is set to `NA`, and mitochondrial filtering is disabled. ",
           "No genes were invented or added to the matrix.")
  }
  raw_mt <- if (nrow(raw_audit)) {
    paste0("Raw source audit: ", sum(raw_audit$n_match_uppercase_MT, na.rm = TRUE),
           " `^MT-` matches and ", sum(raw_audit$n_match_dot_MT, na.rm = TRUE),
           " `^MT.` matches across ", sum(raw_audit$status == "complete"),
           " readable source files.")
  } else "Raw source audit unavailable."
  decision <- c(
    paste0("# ", dataset_id, " mitochondrial-QC decision"), "",
    paste0("- Prepared features: ", length(features)),
    paste0("- Mitochondrial features used: ", n_mt_features),
    paste0("- Feature source: ", mt_source),
    paste0("- `percent.mt` available: ", mt_available),
    paste0("- Fraction of cells with available `percent.mt`: ",
           signif(fraction_cells_mt_available, 5)),
    paste0("- Ribosomal features audited: ", n_ribo_features),
    paste0("- Hemoglobin features audited: ", n_hb_features),
    "", source_note, "", raw_mt, "",
    "The existing stage 06 result was not overwritten. Repaired metadata is stored only under `diagnostics_v2`."
  )
  writeLines(decision, file.path(out, "mt_qc_decision.md"), useBytes = TRUE)
}

message("Mitochondrial feature diagnosis complete")
