#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)

report_dir <- file.path(root, "reports")
manifest_path <- file.path(report_dir, "spatial_validation_result_manifest.csv")
session_path <- file.path(report_dir, "spatial_validation_session_info.txt")
final_path <- file.path(report_dir, "spatial_validation_final.md")

files <- list.files(root, recursive = TRUE, full.names = TRUE)
keep <- grepl("/(results|reference_mapping|figures|reports)/", gsub("\\\\", "/", files)) &
  !grepl("spatial_validation_result_manifest.csv$", files)
manifest <- data.table(
  relative_path = gsub("\\\\", "/", substring(files[keep], nchar(root) + 2)),
  bytes = file.info(files[keep])$size,
  modified_time = as.character(file.info(files[keep])$mtime)
)
if (nrow(manifest) > 0 && requireNamespace("digest", quietly = TRUE)) {
  manifest[, sha256 := vapply(file.path(root, relative_path), digest::digest, character(1),
                              file = TRUE, algo = "sha256")]
} else {
  manifest[, sha256 := NA_character_]
}
fwrite(manifest, manifest_path)
write_session_info(session_path)

grade <- read_dt_if_exists(file.path(root, "results", "meta", "evidence_grade_table.csv"))
qc <- read_dt_if_exists(file.path(root, "results", "spatial_curated", "qc_sensitivity_summary.csv"))
response <- read_dt_if_exists(file.path(root, "results", "gse189843_response", "response_group_statistics.csv"))

lines <- c(
  "# Spatial transcriptomics validation report",
  "",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Data sources and sample audit",
  "",
  "- GSE203612 coordinate-aware ovarian Visium analysis is restricted to GSM6177614 and GSM6177617.",
  "- GSM6177618 is excluded from ovarian results because it is PDAC.",
  "- GSE189843 contributes 12 pretreatment HGSC expression-level samples only; no coordinate-neighborhood claim is made for this series.",
  "",
  "## Scoring definitions",
  "",
  "- Scores are read from `config/spatial_config.yml`, including SPP1_myeloid, target_subclone_02_04, KRAS_hypoxia, target_core_without_receptors, and SPP1_myeloid_without_SPP1.",
  "- SPP1-CD44 is treated as a candidate ligand-receptor axis. SPP1-ITGB1 is reported as an SPP1-associated ITGB1-positive adhesion/integrin program, not as a proven direct ligand-receptor interaction.",
  "",
  "## Statistical boundaries",
  "",
  "- Spot-level spatial tests describe within-sample patterns only.",
  "- Patient/sample-level summaries are the unit for cross-sample statements.",
  "- Two coordinate-aware ovarian samples are insufficient to prove a universal mechanism.",
  "- Low-confidence reference transfer spots are retained as uncertain rather than forced into a unique CNV label.",
  "",
  "## Evidence layers",
  "",
  if (nrow(grade) > 0) paste(capture.output(print(grade)), collapse = "\n") else "Evidence grade table was not available.",
  "",
  "## QC sensitivity",
  "",
  if (nrow(qc) > 0) paste(capture.output(print(qc)), collapse = "\n") else "QC sensitivity could not be summarized because upstream score tables were missing.",
  "",
  "## GSE189843 response analysis",
  "",
  if (nrow(response) > 0) paste(capture.output(print(head(response, 20))), collapse = "\n") else "Response analysis was not available.",
  "",
  "## Output manifest",
  "",
  sprintf("See `%s` for generated file sizes and SHA-256 hashes when the digest package is available.", "reports/spatial_validation_result_manifest.csv")
)
writeLines(lines, final_path)
