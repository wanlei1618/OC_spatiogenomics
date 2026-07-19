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
analysis_dir <- file.path(out_root, "01_clone_patient_confounding")
fig_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
report_dir <- file.path(out_root, "reports")
log_dir <- file.path(out_root, "logs")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("02_clone_patient_confounding_", timestamp, ".log"))
sink(log_file, split = TRUE)
on.exit({
  cat("\nSession info:\n")
  print(sessionInfo())
  sink()
}, add = TRUE)

cat("Starting 02_clone_patient_confounding.R at", as.character(Sys.time()), "\n")
md_path <- file.path(out_root, "00_config", "sample_metadata_checked.csv")
if (!file.exists(md_path)) stop("Checked metadata does not exist. Run 01_build_checked_metadata.R first.")
md <- read.csv(md_path, check.names = FALSE)
clone_md <- md[!is.na(md$clone_label) & nzchar(md$clone_label), ]
cat("Clone-labeled rows:", nrow(clone_md), "\n")

write_count_table <- function(df, group_col, file_name) {
  if (!(group_col %in% names(df)) || all(is.na(df[[group_col]]))) {
    out <- data.frame(status = paste("not_evaluable_missing", group_col), stringsAsFactors = FALSE)
    write.csv(out, file.path(analysis_dir, file_name), row.names = FALSE)
    return(out)
  }
  tab <- as.data.frame(table(df$clone_label, df[[group_col]], useNA = "ifany"))
  names(tab) <- c("clone_label", group_col, "n_cells")
  tab <- tab[tab$n_cells > 0, ]
  write.csv(tab, file.path(analysis_dir, file_name), row.names = FALSE)
  tab
}

clone_by_patient <- write_count_table(clone_md, "patient_id", "clone_by_patient_counts.csv")
clone_by_sample <- write_count_table(clone_md, "sample_id", "clone_by_sample_counts.csv")
clone_by_sample_type <- write_count_table(clone_md, "sample_type", "clone_by_sample_type_counts.csv")
clone_by_batch <- write_count_table(clone_md, "batch", "clone_by_batch_counts.csv")

diversity <- function(x) {
  p <- x / sum(x)
  shannon <- -sum(p * log(p))
  simpson <- 1 - sum(p^2)
  c(shannon = shannon, simpson = simpson)
}

metrics_for_group <- function(df, group_col) {
  clones <- sort(unique(df$clone_label))
  res <- lapply(clones, function(cl) {
    sub <- df[df$clone_label == cl, ]
    if (!(group_col %in% names(sub)) || all(is.na(sub[[group_col]]))) {
      return(data.frame(
        clone_label = cl, grouping = group_col, n_cells = nrow(sub),
        n_groups = NA_integer_, max_single_group_fraction = NA_real_,
        shannon_diversity = NA_real_, simpson_diversity = NA_real_,
        dominant_group = NA_character_, status = "not_evaluable_missing_group",
        stringsAsFactors = FALSE
      ))
    }
    counts <- sort(table(sub[[group_col]]), decreasing = TRUE)
    div <- diversity(as.numeric(counts))
    data.frame(
      clone_label = cl, grouping = group_col, n_cells = nrow(sub),
      n_groups = length(counts),
      max_single_group_fraction = as.numeric(counts[[1]]) / nrow(sub),
      shannon_diversity = div[["shannon"]],
      simpson_diversity = div[["simpson"]],
      dominant_group = names(counts)[[1]],
      status = "evaluated",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

metrics <- do.call(rbind, list(
  metrics_for_group(clone_md, "patient_id"),
  metrics_for_group(clone_md, "sample_id"),
  metrics_for_group(clone_md, "sample_type"),
  metrics_for_group(clone_md, "batch")
))
write.csv(metrics, file.path(analysis_dir, "clone_patient_confounding_metrics.csv"), row.names = FALSE)

cramers_v <- function(tab) {
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  suppressWarnings(ch <- chisq.test(tab))
  chi2 <- unname(ch$statistic)
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  sqrt(chi2 / (n * (k - 1)))
}

assoc_rows <- list()
resid_rows <- list()
for (group_col in c("patient_id", "sample_id", "sample_type", "batch")) {
  if (!(group_col %in% names(clone_md)) || all(is.na(clone_md[[group_col]]))) {
    assoc_rows[[group_col]] <- data.frame(grouping = group_col, cramers_v = NA_real_, p_value = NA_real_,
                                          status = "not_evaluable_missing_group")
    next
  }
  tab <- table(clone_md$clone_label, clone_md[[group_col]])
  suppressWarnings(ch <- chisq.test(tab))
  assoc_rows[[group_col]] <- data.frame(grouping = group_col,
                                        cramers_v = cramers_v(tab),
                                        p_value = unname(ch$p.value),
                                        status = "evaluated")
  resid <- as.data.frame(as.table(ch$stdres))
  names(resid) <- c("clone_label", "group_value", "standardized_residual")
  resid$grouping <- group_col
  resid_rows[[group_col]] <- resid
}
assoc <- do.call(rbind, assoc_rows)
assoc$p_adj_BH <- p.adjust(assoc$p_value, method = "BH")
write.csv(assoc, file.path(analysis_dir, "clone_association_cramers_v.csv"), row.names = FALSE)
if (length(resid_rows)) {
  write.csv(do.call(rbind, resid_rows), file.path(analysis_dir, "clone_group_standardized_residuals.csv"), row.names = FALSE)
}

gate <- metrics[metrics$grouping == "patient_id", ]
patient_evaluable <- any(gate$status == "evaluated")
if (patient_evaluable) {
  gate$gate_a_cross_patient_status <- ifelse(
    gate$n_groups >= as.integer(cfg$min_patients_per_state) &
      gate$max_single_group_fraction < as.numeric(cfg$max_single_patient_fraction_for_cross_patient_state),
    "passes_preliminary_patient_distribution",
    "fails_preliminary_patient_distribution"
  )
} else {
  gate$gate_a_cross_patient_status <- "not_judgable_missing_patient_id"
}
write.csv(gate, file.path(analysis_dir, "gate_a_preliminary_clone_distribution.csv"), row.names = FALSE)

plot_heatmap <- function(count_table, group_col, out_file) {
  if (!("n_cells" %in% names(count_table))) return(FALSE)
  mat <- xtabs(n_cells ~ clone_label + get(group_col), data = count_table)
  png(out_file, width = 1300, height = 850, res = 130)
  op <- par(mar = c(8, 8, 4, 2))
  image(t(log10(mat + 1))[ncol(mat):1, , drop = FALSE], axes = FALSE, col = hcl.colors(50, "YlOrRd"))
  axis(1, at = seq(0, 1, length.out = nrow(mat)), labels = rownames(mat), las = 2, cex.axis = 0.8)
  axis(2, at = seq(0, 1, length.out = ncol(mat)), labels = rev(colnames(mat)), las = 2, cex.axis = 0.7)
  title(paste("Clone by", group_col, "counts, log10(n+1)"))
  par(op)
  dev.off()
  TRUE
}

plot_heatmap(clone_by_sample, "sample_id", file.path(fig_dir, "clone_by_sample_heatmap.png"))
plot_heatmap(clone_by_sample_type, "sample_type", file.path(fig_dir, "clone_by_sample_type_heatmap.png"))
plot_heatmap(clone_by_batch, "batch", file.path(fig_dir, "clone_by_batch_heatmap.png"))
if ("n_cells" %in% names(clone_by_patient)) {
  plot_heatmap(clone_by_patient, "patient_id", file.path(fig_dir, "clone_by_patient_heatmap.png"))
}

if ("n_cells" %in% names(clone_by_sample)) {
  png(file.path(fig_dir, "clone_by_sample_mosaic.png"), width = 1200, height = 850, res = 130)
  mosaicplot(table(clone_md$clone_label, clone_md$sample_id), color = TRUE, las = 2,
             main = "Clone by sample mosaic")
  dev.off()
}

if ("n_cells" %in% names(clone_by_sample_type)) {
  png(file.path(fig_dir, "clone_by_sample_type_alluvial_like.png"), width = 1100, height = 800, res = 130)
  tab <- xtabs(n_cells ~ clone_label + sample_type, data = clone_by_sample_type)
  barplot(tab, beside = FALSE, legend.text = TRUE, las = 2,
          main = "Clone composition by sample_type", ylab = "cells")
  dev.off()
}

report_path <- file.path(report_dir, "01_clone_validity_report.md")
missing_report <- read.csv(file.path(out_root, "00_config", "metadata_field_mapping_report.csv"))
missing_fields <- missing_report$required_field[missing_report$status == "missing"]
sample_metrics <- metrics[metrics$grouping == "sample_id", ]

lines <- c(
  "# First-stage clone validity report",
  "",
  paste("Generated:", as.character(Sys.time())),
  "",
  "## Input",
  paste("- Checked metadata:", md_path),
  paste("- Clone-labeled rows:", nrow(clone_md)),
  "",
  "## Field check",
  paste("- Missing fields:", if (length(missing_fields)) paste(missing_fields, collapse = ", ") else "none"),
  "- No patient_id was inferred from sample_id; cross-patient validity remains blocked until an explicit patient mapping is supplied.",
  "",
  "## Preliminary Gate A",
  if (patient_evaluable) {
    "- patient_id is available, so clone patient distribution can be judged from `gate_a_preliminary_clone_distribution.csv`."
  } else {
    "- patient_id is missing in the selected metadata, so current CNV_Subclone_01-05 labels cannot be called cross-patient CNV states from this dataset alone."
  },
  "- sample_id-level composition was computed as a diagnostic, not as a substitute for patient-level biological replication.",
  "",
  "## Key outputs",
  "- 01_clone_patient_confounding/clone_by_sample_counts.csv",
  "- 01_clone_patient_confounding/clone_by_sample_type_counts.csv",
  "- 01_clone_patient_confounding/clone_by_batch_counts.csv",
  "- 01_clone_patient_confounding/clone_patient_confounding_metrics.csv",
  "- 01_clone_patient_confounding/clone_association_cramers_v.csv",
  "- figures/clone_by_sample_heatmap.png",
  "- figures/clone_by_sample_mosaic.png",
  "- figures/clone_by_sample_type_alluvial_like.png",
  "",
  "## Next required action",
  "Provide or construct an explicit patient_id mapping for sample_id before making cross-patient clone claims or proceeding to old-clone mechanism analyses."
)
writeLines(lines, report_path)
cat("Report written:", report_path, "\n")
