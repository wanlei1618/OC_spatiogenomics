#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required <- c("SummarizedExperiment", "S4Vectors")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE,
                           FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

path <- "E:/OC_spatiogenomics_archive/external_seurat_preannotation_cleanup_20260716/prepared_inputs/GSE158722_raw_counts_sce.rds"
tmp <- paste0(path, ".metadata_fix_tmp")

sce <- readRDS(path)
md <- as.data.frame(SummarizedExperiment::colData(sce))
original <- as.character(md$original_cell_id)
parts <- strsplit(original, "_", fixed = TRUE)
patient_id <- vapply(parts, `[[`, character(1), 1L)
patient_number <- as.integer(sub("^P", "", patient_id))

longitudinal <- patient_number <= 9L
timepoint <- character(length(parts))
timepoint[longitudinal] <- vapply(parts[longitudinal], `[[`, character(1), 2L)
timepoint[!longitudinal & patient_number <= 17L] <- "Pre-treatment"
timepoint[!longitudinal & patient_number >= 18L] <- "Post-treatment"

barcode <- vapply(seq_along(parts), function(i) {
  x <- parts[[i]]
  if (longitudinal[[i]]) paste(x[-(1:2)], collapse = "_")
  else paste(x[-1], collapse = "_")
}, character(1))

sample_id <- paste(patient_id, timepoint, sep = "_")
cell_id <- paste(patient_id, sample_id, barcode, sep = "__")

if (anyNA(timepoint) || any(!nzchar(timepoint))) stop("Missing timepoint values")
if (anyDuplicated(cell_id)) stop("Corrected cell IDs are not unique")

colnames(sce) <- cell_id
SummarizedExperiment::colData(sce)$patient_id <- patient_id
SummarizedExperiment::colData(sce)$sample_id <- sample_id
SummarizedExperiment::colData(sce)$timepoint <- timepoint

cat("cells", ncol(sce), "patients", length(unique(patient_id)),
    "samples", length(unique(sample_id)), "\n")
print(table(timepoint))

saveRDS(sce, tmp, compress = FALSE)
cat("saved", tmp, file.info(tmp)$size, "bytes\n")
