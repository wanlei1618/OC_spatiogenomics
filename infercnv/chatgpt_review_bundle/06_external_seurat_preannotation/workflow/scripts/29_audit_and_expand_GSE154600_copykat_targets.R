options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config()
cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
cleaned <- file.path(data_root, "diagnostics_v2_marker_ready_cleaned", "GSE154600")
v61 <- file.path(data_root, "diagnostics_v6_1_copykat_stability")
dir.create(v61, recursive = TRUE, showWarnings = FALSE)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(v61, "GSE154600_copykat_target_coverage_audit.csv")
if (file.exists(out) && !replace_generated) stop("Output exists: ", out)

a <- fread(file.path(cleaned, "cleaned_cell_assignments.csv.gz"))
c <- fread(file.path(cleaned, "cleaned_cluster_annotation_template.csv"))
idx <- match(as.character(a$final_cluster), as.character(c$cluster))
a[, `:=`(
  annotation_status = c$annotation_status[idx],
  canonical_support_n = as.numeric(c$canonical_support_n[idx]),
  incompatible_lineage_program = as.logical(c$incompatible_lineage_program[idx])
)]

lineage_names <- c("Epithelial_like", "T_NK_like", "B_Plasma_like", "Myeloid_like")
lineage_paths <- file.path(
  data_root, "diagnostics_v2", "objects", "GSE154600", "lineage_inputs",
  paste0(lineage_names, "_strategy_input.rds")
)
lineage_cells <- setNames(lapply(lineage_paths, function(p) {
  x <- readRDS(p)
  colnames(x$counts)
}), lineage_names)

locate_cells <- function(cells) {
  hit <- vapply(lineage_cells, function(z) cells %in% z, logical(length(cells)))
  if (is.null(dim(hit))) hit <- matrix(hit, ncol = length(lineage_cells))
  colnames(hit) <- names(lineage_cells)
  nmatch <- rowSums(hit)
  first <- apply(hit, 1, function(x) {
    z <- which(x)
    if (length(z)) lineage_names[z[1]] else NA_character_
  })
  sources <- apply(hit, 1, function(x) paste(lineage_names[x], collapse = ";"))
  data.table(
    cell_id = cells,
    source_lineage_input = first,
    source_lineage_inputs = sources,
    n_source_matches = nmatch,
    available_for_copykat = nmatch > 0,
    missing_reason = ifelse(nmatch > 0, NA_character_, "missing_from_all_lineage_inputs")
  )
}

targets <- a[
  final_cell_type == "Epithelial" &
    annotation_status == "REVIEW_PATIENT_ENRICHED" &
    canonical_support_n >= 3 &
    incompatible_lineage_program != TRUE,
  .(dataset_id, patient_id, cell_id, final_cell_type)
]
target_manifest <- merge(targets, locate_cells(targets$cell_id), by = "cell_id", all.x = TRUE)
setcolorder(target_manifest, c(
  "dataset_id", "patient_id", "cell_id", "final_cell_type",
  "source_lineage_input", "source_lineage_inputs", "n_source_matches",
  "available_for_copykat", "missing_reason"
))
fwrite(
  target_manifest, file.path(v61, "GSE154600_copykat_target_manifest.csv.gz"),
  compress = "gzip", na = "NA"
)

reference_types <- c(
  "T_cell", "NK_cell", "B_cell", "Macrophage",
  "Monocyte", "cDC1", "cDC2", "pDC"
)
refs <- a[
  final_cell_type %in% reference_types &
    annotation_status == "READY_HIGH_CONFIDENCE" &
    canonical_support_n >= 3 &
    incompatible_lineage_program != TRUE,
  .(dataset_id, patient_id, cell_id, final_cell_type)
]
reference_manifest <- merge(refs, locate_cells(refs$cell_id), by = "cell_id", all.x = TRUE)
setcolorder(reference_manifest, c(
  "dataset_id", "patient_id", "cell_id", "final_cell_type",
  "source_lineage_input", "source_lineage_inputs", "n_source_matches",
  "available_for_copykat", "missing_reason"
))
fwrite(reference_manifest, file.path(v61, "GSE154600_copykat_reference_manifest.csv.gz"),
       compress = "gzip", na = "NA")

coverage <- target_manifest[, .(
  n_final_epithelial = .N,
  n_available_in_any_lineage_input = sum(available_for_copykat),
  n_missing_from_all_inputs = sum(!available_for_copykat),
  coverage_fraction = mean(available_for_copykat),
  n_original_epi_input = sum(grepl("(^|;)Epithelial_like($|;)", source_lineage_inputs)),
  n_rescued_from_other_inputs = sum(
    available_for_copykat &
      !grepl("(^|;)Epithelial_like($|;)", source_lineage_inputs)
  )
), by = .(dataset_id, patient_id)]
fwrite(coverage, out, na = "NA")

reference_audit <- reference_manifest[, .(
  n_final_eligible_reference = .N,
  n_available_in_any_lineage_input = sum(available_for_copykat),
  n_missing_from_all_inputs = sum(!available_for_copykat),
  coverage_fraction = mean(available_for_copykat),
  n_multiple_source_matches = sum(n_source_matches > 1)
), by = .(dataset_id, patient_id, final_cell_type)]
fwrite(reference_audit, file.path(v61, "GSE154600_copykat_reference_audit.csv"), na = "NA")
message("GSE154600 cross-lineage CopyKAT target and reference audit complete")
