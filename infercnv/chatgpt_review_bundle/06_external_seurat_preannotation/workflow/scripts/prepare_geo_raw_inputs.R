#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required <- c("Matrix", "data.table", "Seurat", "SingleCellExperiment", "S4Vectors")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                           FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

root <- "D:/OC_spatiogenomics/公开集/单细胞/raw_counts"
prepared <- "E:/OC_spatiogenomics_archive/external_seurat_preannotation_cleanup_20260716/prepared_inputs"
dir.create(prepared, recursive = TRUE, showWarnings = FALSE)

as_dgc <- function(x) {
  methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
}

table_columns_to_dgc <- function(tab, count_columns, block_size = 500L) {
  groups <- split(count_columns,
                  ceiling(seq_along(count_columns) / block_size))
  blocks <- lapply(groups, function(cols) {
    values <- as.matrix(tab[, cols, with = FALSE])
    out <- as_dgc(values)
    rm(values); gc()
    out
  })
  do.call(cbind, blocks)
}

combine_samples <- function(mats, metadata, dataset_id) {
  if (!length(mats)) stop(dataset_id, ": no matrices")
  all_genes <- Reduce(union, lapply(mats, rownames))
  same <- vapply(mats, function(x) identical(rownames(x), all_genes), logical(1))
  if (!all(same)) {
    cat(dataset_id, ": aligning samples to union of", length(all_genes), "genes\n")
    mats <- lapply(mats, function(x) {
      row_map <- match(rownames(x), all_genes)
      Matrix::sparseMatrix(
        i = row_map[x@i + 1L], p = x@p, x = x@x,
        dims = c(length(all_genes), ncol(x)),
        dimnames = list(all_genes, colnames(x))
      )
    })
  }
  counts <- do.call(cbind, mats)
  md <- do.call(rbind, metadata)
  if (!identical(colnames(counts), rownames(md))) {
    stop(dataset_id, ": metadata row names do not match cell IDs")
  }
  if (anyDuplicated(colnames(counts))) stop(dataset_id, ": duplicate cell IDs")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(md)
  )
  out <- file.path(prepared, paste0(dataset_id, "_raw_counts_sce.rds"))
  saveRDS(sce, out, compress = FALSE)
  cat(dataset_id, nrow(sce), ncol(sce), out, "\n")
  invisible(out)
}

make_metadata <- function(cell_ids, dataset_id, sample_id,
                          patient_id = sample_id, timepoint = NA_character_,
                          original_cell_id = NA_character_) {
  data.frame(
    dataset_id = dataset_id,
    sample_id = sample_id,
    patient_id = patient_id,
    timepoint = timepoint,
    original_cell_id = original_cell_id,
    row.names = cell_ids,
    check.names = FALSE
  )
}

prepare_gse147082 <- function() {
  dataset_id <- "GSE147082"
  files <- sort(list.files(file.path(root, dataset_id, "extracted"),
                           pattern = "\\.csv\\.gz$", full.names = TRUE))
  mats <- metadata <- vector("list", length(files))
  for (i in seq_along(files)) {
    f <- files[[i]]
    sample_id <- sub("^GSM[0-9]+_", "",
                     sub("\\.csv\\.gz$", "", basename(f)))
    cat(dataset_id, sample_id, basename(f), "\n")
    tab <- data.table::fread(f, check.names = FALSE)
    genes <- make.unique(as.character(tab[[1]]))
    mat <- table_columns_to_dgc(tab, seq.int(2L, ncol(tab)))
    rownames(mat) <- genes
    original <- colnames(tab)[-1]
    cells <- paste(sample_id, original, sep = "__")
    colnames(mat) <- cells
    mats[[i]] <- mat
    metadata[[i]] <- make_metadata(cells, dataset_id, sample_id,
                                    original_cell_id = original)
    rm(tab, mat); gc()
  }
  combine_samples(mats, metadata, dataset_id)
}

prepare_gse151214 <- function() {
  dataset_id <- "GSE151214"
  files <- sort(list.files(file.path(root, dataset_id, "extracted"),
                           pattern = "\\.h5$", full.names = TRUE))
  mats <- metadata <- vector("list", length(files))
  for (i in seq_along(files)) {
    f <- files[[i]]
    sample_id <- sub("^GSM[0-9]+_", "", sub("\\.h5$", "", basename(f)))
    cat(dataset_id, sample_id, basename(f), "\n")
    mat <- Seurat::Read10X_h5(f, use.names = TRUE, unique.features = TRUE)
    if (is.list(mat)) mat <- mat[[1]]
    mat <- as_dgc(mat)
    original <- colnames(mat)
    cells <- paste(sample_id, original, sep = "__")
    colnames(mat) <- cells
    mats[[i]] <- mat
    metadata[[i]] <- make_metadata(cells, dataset_id, sample_id,
                                    original_cell_id = original)
    rm(mat); gc()
  }
  combine_samples(mats, metadata, dataset_id)
}

prepare_gse154600 <- function() {
  dataset_id <- "GSE154600"
  dir <- file.path(root, dataset_id, "extracted")
  files <- sort(list.files(dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE))
  mats <- metadata <- vector("list", length(files))
  for (i in seq_along(files)) {
    matrix_file <- files[[i]]
    stem <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file))
    sample_id <- sub("^GSM[0-9]+_", "", stem)
    gene_file <- file.path(dir, paste0(stem, "_genes.tsv.gz"))
    barcode_file <- file.path(dir, paste0(stem, "_barcodes.tsv.gz"))
    cat(dataset_id, sample_id, basename(matrix_file), "\n")
    genes <- data.table::fread(gene_file, header = FALSE)[[2]]
    barcodes <- data.table::fread(barcode_file, header = FALSE)[[1]]
    mat <- Matrix::readMM(gzfile(matrix_file))
    mat <- methods::as(mat, "dgCMatrix")
    rownames(mat) <- make.unique(as.character(genes))
    cells <- paste(sample_id, barcodes, sep = "__")
    colnames(mat) <- cells
    mats[[i]] <- mat
    md <- make_metadata(cells, dataset_id, sample_id,
                        original_cell_id = barcodes)
    author_file <- file.path(
      "D:/OC_spatiogenomics/infercnv/external_cell_annotations/raw/GSE154600",
      paste0("sample", sub("^T", "", sample_id), "_sce.rds")
    )
    if (file.exists(author_file)) {
      author <- readRDS(author_file)
      author_md <- as.data.frame(SummarizedExperiment::colData(author))
      if (!"Barcode" %in% colnames(author_md)) {
        stop(dataset_id, " ", sample_id, ": author metadata lacks Barcode")
      }
      author_index <- match(barcodes, as.character(author_md[["Barcode"]]))
      fields <- intersect(
        c("Cluster", "subtype", "margin", "hpca.celltype", "encode.celltype",
          "hpca.celltype.score", "encode.celltype.score", "celltype",
          "IMR_consensus", "DIF_consensus", "PRO_consensus", "MES_consensus"),
        colnames(author_md)
      )
      for (field in fields) md[[field]] <- author_md[[field]][author_index]
      md$author_metadata_matched <- !is.na(author_index)
      cat(dataset_id, sample_id, "author metadata matched",
          sum(md$author_metadata_matched), "of", nrow(md), "cells\n")
    }
    metadata[[i]] <- md
    rm(mat); gc()
  }
  combine_samples(mats, metadata, dataset_id)
}

prepare_gse158722 <- function() {
  dataset_id <- "GSE158722"
  common_file <- file.path(prepared, "GSE158722_common_genes.rds")
  if (!file.exists(common_file)) {
    stop(dataset_id, ": common gene list is missing: ", common_file)
  }
  common_genes <- readRDS(common_file)
  files <- sort(list.files(file.path(root, dataset_id),
                           pattern = "\\.counts\\.txt\\.gz$", full.names = TRUE))
  mats <- metadata <- vector("list", length(files))
  for (i in seq_along(files)) {
    f <- files[[i]]
    patient_id <- sub("^GSE158722_", "",
                      sub("\\.counts\\.txt\\.gz$", "", basename(f)))
    cat(dataset_id, patient_id, basename(f), "\n")
    tab <- data.table::fread(f, check.names = FALSE)
    gene_column <- grep("^Gene[ _]Symbol$", names(tab), ignore.case = TRUE)
    if (length(gene_column) != 1L) {
      stop(dataset_id, " ", patient_id, ": cannot identify Gene Symbol column")
    }
    genes <- make.unique(as.character(tab[[gene_column]]))
    row_index <- match(common_genes, genes)
    if (anyNA(row_index)) {
      stop(dataset_id, " ", patient_id, ": common gene missing unexpectedly")
    }
    tab <- tab[row_index]
    genes <- common_genes
    annotation_columns <- grep("^(ENSEMBL[ _]ID|Gene[ _]ID|Gene[ _]Symbol)$",
                               names(tab), ignore.case = TRUE)
    count_columns <- setdiff(seq_len(ncol(tab)), annotation_columns)
    numeric_counts <- vapply(tab[, count_columns, with = FALSE],
                             is.numeric, logical(1))
    if (!all(numeric_counts)) {
      stop(dataset_id, " ", patient_id,
           ": non-numeric count column(s): ",
           paste(names(numeric_counts)[!numeric_counts], collapse = ", "))
    }
    original <- colnames(tab)[count_columns]
    parts <- strsplit(original, "_", fixed = TRUE)
    timepoint <- vapply(parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_,
                        character(1))
    sample_id <- paste(patient_id, timepoint, sep = "_")
    barcode <- vapply(parts, function(x) {
      if (length(x) >= 3) paste(x[-(1:2)], collapse = "_") else tail(x, 1)
    }, character(1))
    cells <- paste(patient_id, sample_id, barcode, sep = "__")
    mat <- table_columns_to_dgc(tab, count_columns)
    rownames(mat) <- genes
    colnames(mat) <- cells
    mats[[i]] <- mat
    metadata[[i]] <- make_metadata(cells, dataset_id, sample_id,
                                    patient_id, timepoint, original)
    rm(tab, mat); gc()
  }
  combine_samples(mats, metadata, dataset_id)
}

targets <- commandArgs(trailingOnly = TRUE)
if (!length(targets)) {
  targets <- c("GSE147082", "GSE151214", "GSE154600", "GSE158722")
}
preparers <- list(
  GSE147082 = prepare_gse147082,
  GSE151214 = prepare_gse151214,
  GSE154600 = prepare_gse154600,
  GSE158722 = prepare_gse158722
)
unknown <- setdiff(targets, names(preparers))
if (length(unknown)) stop("Unknown target(s): ", paste(unknown, collapse = ", "))
for (target in targets) preparers[[target]]()
