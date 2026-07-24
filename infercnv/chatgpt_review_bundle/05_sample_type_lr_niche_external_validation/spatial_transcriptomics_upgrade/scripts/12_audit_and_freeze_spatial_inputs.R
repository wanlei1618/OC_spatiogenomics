options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

repo_root <- normalizePath(file.path(
  dirname(normalizePath(sub("^--file=", "", grep(
    "^--file=", commandArgs(), value = TRUE
  )[[1L]]))), "..", "..", "..", "..", ".."
), winslash = "/", mustWork = TRUE)
local_root <- "D:/OC_spatiogenomics/spatial_data/spatial_analysis_correction_v2"
dir.create(local_root, recursive = TRUE, showWarnings = FALSE)
old_upgrade <- paste0(
  "D:/OC_spatiogenomics/infercnv/",
  "spatial_transcriptomics_validation_2026-07-12/",
  "spatial_transcriptomics_upgrade"
)
old_manifest <- fread(file.path(
  old_upgrade, "metadata", "spatial_sample_manifest.csv"
))

objects_path <- paste0(
  "D:/OC_spatiogenomics/spatial_data/processed/",
  "spatial_objects_curated_scored_reference_mapped.rds"
)
objects <- readRDS(objects_path)
audit_203 <- rbindlist(lapply(
  c("GSM6177614", "GSM6177617", "GSM6177618"), function(sid) {
    old <- old_manifest[sample_id == sid]
    available <- sid %in% names(objects)
    match_fraction <- if (available) {
      meta <- objects[[sid]]@meta.data
      mean(rownames(meta) %in% colnames(objects[[sid]]))
    } else NA_real_
    data.table(
      dataset_id = "GSE203612", sample_id = sid,
      reported_tumor_type = old$disease,
      verified_tumor_type = fifelse(
        sid == "GSM6177618", "pancreatic ductal adenocarcinoma",
        "ovarian carcinoma"
      ),
      platform = "10x Genomics Visium",
      counts_available = available,
      coordinates_available = available &&
        all(c("array_row", "array_col") %in% names(objects[[sid]]@meta.data)),
      coordinate_barcode_match_fraction = match_fraction,
      include_coordinate_analysis = sid %chin% c("GSM6177614", "GSM6177617"),
      include_expression_analysis = sid %chin% c("GSM6177614", "GSM6177617"),
      analysis_role = fifelse(
        sid == "GSM6177618", "EXCLUDED_WRONG_TUMOR_TYPE",
        "PRIMARY_COORDINATE_AWARE_OVARIAN"
      ),
      exclusion_reason = fifelse(
        sid == "GSM6177618", "PDAC_NOT_OVARIAN", NA_character_
      ),
      source_path = fifelse(available, objects_path, old$source_url),
      response_group = NA_character_,
      included_by_author = TRUE,
      available_in_matrix = available,
      final_include = sid %chin% c("GSM6177614", "GSM6177617")
    )
  }
))

extract_dir <- "D:/OC_spatiogenomics/spatial_data/raw/GSE189843/extracted"
old_189 <- old_manifest[dataset == "GSE189843"]
audit_189 <- old_189[, {
  matrix_file <- list.files(
    extract_dir, pattern = paste0("^", sample_id, ".*matrix.*\\.mtx(\\.gz)?$"),
    full.names = TRUE
  )
  feature_file <- list.files(
    extract_dir, pattern = paste0("^", sample_id, ".*features.*\\.tsv(\\.gz)?$"),
    full.names = TRUE
  )
  barcode_file <- list.files(
    extract_dir, pattern = paste0("^", sample_id, ".*barcodes.*\\.tsv(\\.gz)?$"),
    full.names = TRUE
  )
  matrix_ok <- length(matrix_file) > 0L &&
    length(feature_file) > 0L && length(barcode_file) > 0L
  data.table(
    dataset_id = "GSE189843", sample_id,
    reported_tumor_type = disease,
    verified_tumor_type = "high-grade serous ovarian carcinoma",
    platform = "10x Genomics Visium FFPE",
    counts_available = matrix_ok,
    coordinates_available = FALSE,
    coordinate_barcode_match_fraction = NA_real_,
    include_coordinate_analysis = FALSE,
    include_expression_analysis = matrix_ok,
    analysis_role = "EXPRESSION_ONLY_PATIENT_LEVEL",
    exclusion_reason = fifelse(
      matrix_ok, "COORDINATES_NOT_AUTHORITATIVELY_VERIFIED",
      "EXPRESSION_MATRIX_INCOMPLETE"
    ),
    source_path = if (matrix_ok) matrix_file[[1L]] else extract_dir,
    response_group = clinical_group,
    included_by_author = TRUE,
    available_in_matrix = matrix_ok,
    final_include = matrix_ok
  )
}, by = seq_len(nrow(old_189))]
audit_189[, seq_len := NULL]
audit <- rbind(audit_203, audit_189, fill = TRUE)
fwrite(audit, file.path(local_root, "spatial_sample_audit.csv"), na = "NA")

old_tables <- file.path(
  repo_root,
  "infercnv/chatgpt_review_bundle/05_sample_type_lr_niche_external_validation",
  "sample_type_LR_niche_analysis", "tables"
)
deprecated <- data.table(
  output_or_conclusion = c(
    "all historical GSM6177618 spatial results",
    "historical three-sample ovarian pooled conclusion",
    "historical COMMOT-named expression-product outputs",
    "historical spot-level cross-sample significance"
  ),
  source_path = c(
    file.path(old_upgrade, "results"),
    file.path(old_tables, "spatial_neighborhood_enrichment.csv"),
    file.path(old_upgrade, "results"),
    file.path(old_tables, "spatial_correlation_SPP1_myeloid_Target_axis.csv")
  ),
  correction_status = c(
    "DEPRECATED_WRONG_TUMOR_TYPE",
    "DEPRECATED_SAMPLE_CONTAMINATION",
    "SIMPLIFIED_EXPRESSION_PRODUCT_NOT_COMMOT",
    "PSEUDOREPLICATION_RISK"
  ),
  replacement = c(
    "excluded from every ovarian result",
    "two verified ovarian sections analyzed separately",
    "expression-product summary only; no optimal transport claim",
    "sample/patient is the statistical unit"
  )
)
fwrite(deprecated, file.path(local_root, "deprecated_spatial_outputs.csv"),
       na = "NA")
message("Spatial input audit frozen")
