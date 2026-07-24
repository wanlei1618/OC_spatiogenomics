options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

local_root <- "D:/OC_spatiogenomics/spatial_data/spatial_analysis_correction_v2"
input_root <- Sys.glob(file.path(
  "D:/OC_spatiogenomics", "*", "ovarian_spatial_geo", "GSE211956", "suppl"
))
matrix_files <- if (length(input_root) == 1L) {
  list.files(input_root, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
} else character()
coordinate_files <- if (length(input_root) == 1L) {
  list.files(input_root, pattern = "tissue_positions_list\\.csv$",
             recursive = TRUE, full.names = TRUE)
} else character()

# Counts and coordinates are locally present, but the correction task requires
# a locally authoritative, sample-level HG-SOC identity before using this
# dataset as independent biological replication. That identity record is not
# present in the frozen local registry, so no coordinate analysis is run.
status <- data.table(
  dataset_id = "GSE211956",
  status = "NOT_RUN_INPUT_INCOMPLETE",
  counts_files_found = length(matrix_files),
  coordinate_files_found = length(coordinate_files),
  tissue_identity_verified = FALSE,
  reason = paste(
    "complete matrices and coordinates are present, but an authoritative",
    "sample-level HG-SOC identity record is absent from the frozen registry"
  ),
  action = "no download and no spatial replication performed"
)
fwrite(status, file.path(local_root, "GSE211956_status.csv"), na = "NA")
message("GSE211956 not run: authoritative HG-SOC identity not frozen")
