options(stringsAsFactors = FALSE)

read_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
  cfg <- list()
  for (line in lines) {
    key <- sub(":.*$", "", line)
    value <- trimws(sub("^[^:]+:", "", line))
    value <- gsub('^"|"$', "", value)
    cfg[[key]] <- value
  }
  cfg
}

norm_path <- function(x) gsub("/", "\\\\", x)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
config_path <- "D:/OC_spatiogenomics/infercnv/06_rigorous_validation/00_config/project_config.yaml"
cfg <- read_config(config_path)
set.seed(as.integer(cfg$random_seed))

out_root <- norm_path(cfg$output_root)
log_dir <- file.path(out_root, "logs")
table_dir <- file.path(out_root, "tables")
config_dir <- file.path(out_root, "00_config")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("01_build_checked_metadata_", timestamp, ".log"))
sink(log_file, split = TRUE)
on.exit({
  cat("\nSession info:\n")
  print(sessionInfo())
  sink()
}, add = TRUE)

cat("Starting 01_build_checked_metadata.R at", as.character(Sys.time()), "\n")
cat("Config:", config_path, "\n")

project_root <- norm_path(cfg$project_root)
input_path <- norm_path(cfg$input_metadata_path)
fallback_path <- norm_path(cfg$fallback_metadata_path)
if (!file.exists(input_path)) {
  if (file.exists(fallback_path)) {
    input_path <- fallback_path
  } else {
    stop("No configured metadata CSV exists.")
  }
}
cat("Metadata input:", input_path, "\n")

all_files <- list.files(project_root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
manifest <- data.frame(
  path = gsub("\\\\", "/", all_files),
  extension = tolower(tools::file_ext(all_files)),
  size_bytes = file.info(all_files)$size,
  modified_time = as.character(file.info(all_files)$mtime),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(table_dir, "input_file_manifest.csv"), row.names = FALSE)
cat("Input manifest files:", nrow(manifest), "\n")

md <- read.csv(input_path, check.names = FALSE)
cat("Metadata rows:", nrow(md), "columns:", ncol(md), "\n")

pick_col <- function(candidates, cols) {
  hit <- candidates[candidates %in% cols]
  if (length(hit)) hit[[1]] else NA_character_
}

cols <- colnames(md)
mapping <- data.frame(
  required_field = c("cell_id", "patient_id", "sample_id", "sample_type", "batch", "dataset",
                     "clone", "cell_type", "nCount_RNA", "nFeature_RNA", "percent.mt",
                     "doublet_score", "S.Score", "G2M.Score"),
  selected_column = NA_character_,
  status = NA_character_,
  note = NA_character_,
  stringsAsFactors = FALSE
)

candidates <- list(
  cell_id = c("cell_integrated_oc", "cell", "barcode", "cell_canonical"),
  patient_id = c("patient_id", "patient", "donor_id", "case_id"),
  sample_id = c("sample_id", "orig.ident", "batch"),
  sample_type = c("sample_type", "tissue", "site"),
  batch = c("batch", "orig.ident"),
  dataset = c("dataset", "cohort", "study"),
  clone = c("clone", "CNV_clone", "cnv_subclone", "interaction_group"),
  cell_type = c("cell_type", "annotation", "celltype"),
  nCount_RNA = c("nCount_RNA", "nCount"),
  nFeature_RNA = c("nFeature_RNA", "nFeature"),
  percent.mt = c("percent.mt", "pct_counts_mt", "mito_percent"),
  doublet_score = c("doublet_score", "DoubletFinder_score", "scDblFinder.score"),
  S.Score = c("S.Score", "S_score"),
  G2M.Score = c("G2M.Score", "G2M_score")
)

for (i in seq_len(nrow(mapping))) {
  field <- mapping$required_field[i]
  selected <- pick_col(candidates[[field]], cols)
  mapping$selected_column[i] <- selected
  mapping$status[i] <- if (is.na(selected)) "missing" else "present"
  mapping$note[i] <- if (is.na(selected)) {
    paste("No exact candidate column found among:", paste(candidates[[field]], collapse = ", "))
  } else {
    ""
  }
}

write.csv(mapping, file.path(config_dir, "metadata_field_mapping_report.csv"), row.names = FALSE)

get_or_na <- function(field) {
  selected <- mapping$selected_column[mapping$required_field == field]
  if (length(selected) == 0 || is.na(selected)) rep(NA_character_, nrow(md)) else as.character(md[[selected]])
}

clone_raw <- get_or_na("clone")
clone_label <- clone_raw
clone_label[grepl("^Subclone_", clone_label)] <- paste0("CNV_", clone_label[grepl("^Subclone_", clone_label)])
clone_label[grepl("^CNV_Subclone_", clone_label)] <- clone_label[grepl("^CNV_Subclone_", clone_label)]
clone_label[!grepl("^CNV_Subclone_", clone_label)] <- NA_character_

checked <- data.frame(
  cell_id = get_or_na("cell_id"),
  patient_id = get_or_na("patient_id"),
  sample_id = get_or_na("sample_id"),
  sample_type = get_or_na("sample_type"),
  batch = get_or_na("batch"),
  dataset = get_or_na("dataset"),
  cell_type = get_or_na("cell_type"),
  clone_raw = clone_raw,
  clone_label = clone_label,
  nCount_RNA = suppressWarnings(as.numeric(get_or_na("nCount_RNA"))),
  nFeature_RNA = suppressWarnings(as.numeric(get_or_na("nFeature_RNA"))),
  percent.mt = suppressWarnings(as.numeric(get_or_na("percent.mt"))),
  doublet_score = suppressWarnings(as.numeric(get_or_na("doublet_score"))),
  S.Score = suppressWarnings(as.numeric(get_or_na("S.Score"))),
  G2M.Score = suppressWarnings(as.numeric(get_or_na("G2M.Score"))),
  stringsAsFactors = FALSE
)

write.csv(checked, file.path(config_dir, "sample_metadata_checked.csv"), row.names = FALSE)

summary_df <- data.frame(
  metric = c("metadata_input", "n_rows", "n_columns", "n_clone_labeled_cells",
             "n_unique_samples", "n_unique_patients_nonmissing", "n_missing_patient_id",
             "n_missing_dataset", "random_seed"),
  value = c(input_path, nrow(md), ncol(md), sum(!is.na(checked$clone_label)),
            length(unique(checked$sample_id[!is.na(checked$sample_id)])),
            length(unique(checked$patient_id[!is.na(checked$patient_id)])),
            sum(is.na(checked$patient_id)), sum(is.na(checked$dataset)), cfg$random_seed),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(table_dir, "metadata_check_summary.csv"), row.names = FALSE)

cat("Checked metadata written.\n")
cat("Missing required fields:", paste(mapping$required_field[mapping$status == "missing"], collapse = ", "), "\n")
