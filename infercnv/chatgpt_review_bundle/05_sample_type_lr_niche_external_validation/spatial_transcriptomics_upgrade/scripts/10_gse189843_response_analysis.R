#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]),
                                      winslash = "/", mustWork = FALSE)), "00_spatial_validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- load_config(if (length(args) >= 1) args[[1]] else NULL)
root <- path_root(config)
ensure_dirs(root)
out_dir <- file.path(root, "results", "gse189843_response")

spot <- tryCatch(read_spot_scores(root), error = function(e) data.table())
gse <- spot[dataset == "GSE189843"]
score_cols <- intersect(c("SPP1_myeloid_score", "CD44", "ITGB1", "SPP1_CD44_expr_product",
                          "SPP1_ITGB1_expr_product", "target_subclone_02_04_score",
                          "KRAS_hypoxia_score"), names(gse))

score_matrix <- function(counts, genes) {
  present <- intersect(genes, rownames(counts))
  if (length(present) == 0) return(rep(NA_real_, ncol(counts)))
  x <- as.matrix(counts[present, , drop = FALSE])
  x <- log1p(t(t(x) / pmax(colSums(counts), 1)) * 10000)
  colMeans(scale(t(x)), na.rm = TRUE)
}

gene_expr <- function(counts, gene) {
  if (!gene %in% rownames(counts)) return(rep(NA_real_, ncol(counts)))
  lib <- pmax(colSums(counts), 1)
  as.numeric(log1p(counts[gene, ] / lib * 10000))
}

read_mtx_selected <- function(matrix_file, feature_file, barcode_file, wanted_genes) {
  features <- fread(feature_file, header = FALSE)
  barcodes <- fread(barcode_file, header = FALSE)
  gene_names <- make.unique(as.character(features[[if (ncol(features) >= 2) 2 else 1]]))
  wanted <- intersect(unique(wanted_genes), gene_names)
  header <- readLines(matrix_file, n = 20)
  dim_line_index <- which(!startsWith(header, "%"))[1]
  dims <- as.integer(strsplit(header[dim_line_index], "\\s+")[[1]])
  trip <- fread(matrix_file, skip = dim_line_index, col.names = c("gene_idx", "spot_idx", "count"))
  lib <- trip[, .(library_size = sum(count)), by = spot_idx]
  selected_idx <- which(gene_names %in% wanted)
  trip <- trip[gene_idx %in% selected_idx]
  gene_map <- data.table(gene_idx = selected_idx, gene = gene_names[selected_idx])
  trip <- merge(trip, gene_map, by = "gene_idx", all.x = TRUE)
  mat <- matrix(0, nrow = length(wanted), ncol = dims[2], dimnames = list(wanted, as.character(barcodes[[1]])))
  if (nrow(trip) > 0) {
    idx <- cbind(match(trip$gene, wanted), trip$spot_idx)
    mat[idx] <- trip$count
  }
  lib_vec <- rep(1, dims[2])
  lib_vec[lib$spot_idx] <- pmax(lib$library_size, 1)
  norm <- log1p(t(t(mat) / lib_vec) * 10000)
  list(norm = norm, n_spots = dims[2])
}

score_from_norm <- function(norm, genes) {
  present <- intersect(genes, rownames(norm))
  if (length(present) == 0) return(rep(NA_real_, ncol(norm)))
  x <- norm[present, , drop = FALSE]
  if (nrow(x) == 1) return(as.numeric(scale(x[1, ])))
  colMeans(scale(t(x)), na.rm = TRUE)
}

expr_from_norm <- function(norm, gene) {
  if (!gene %in% rownames(norm)) return(rep(NA_real_, ncol(norm)))
  as.numeric(norm[gene, ])
}

if (nrow(gse) == 0) {
  manifest <- fread(file.path(script_path(), "..", "metadata", "spatial_sample_manifest.csv"))
  manifest <- manifest[dataset == "GSE189843" & toupper(include_in_ovarian_analysis) == "TRUE"]
  extract_dir <- file.path(root, "raw", "GSE189843", "extracted")
  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    sid <- manifest$sample_id[i]
    matrix_file <- list.files(extract_dir, pattern = paste0("^", sid, ".*matrix.*\\.mtx$"), full.names = TRUE)[1]
    feature_file <- list.files(extract_dir, pattern = paste0("^", sid, ".*features.*\\.tsv$"), full.names = TRUE)[1]
    barcode_file <- list.files(extract_dir, pattern = paste0("^", sid, ".*barcodes.*\\.tsv$"), full.names = TRUE)[1]
    if (any(is.na(c(matrix_file, feature_file, barcode_file)))) next
    wanted <- unique(c(unlist(config$gene_sets$spp1_myeloid),
                       unlist(config$gene_sets$target_subclone_02_04),
                       unlist(config$gene_sets$kras_hypoxia),
                       "SPP1", "CD44", "ITGB1"))
    selected <- read_mtx_selected(matrix_file, feature_file, barcode_file, wanted)
    norm <- selected$norm
    spp1 <- score_from_norm(norm, unlist(config$gene_sets$spp1_myeloid))
    target <- score_from_norm(norm, unlist(config$gene_sets$target_subclone_02_04))
    kras <- score_from_norm(norm, unlist(config$gene_sets$kras_hypoxia))
    spp1_expr <- expr_from_norm(norm, "SPP1")
    cd44 <- expr_from_norm(norm, "CD44")
    itgb1 <- expr_from_norm(norm, "ITGB1")
    rows[[length(rows) + 1]] <- data.table(
      sample_id = sid,
      clinical_group = manifest$clinical_group[i],
      n_spots = selected$n_spots,
      SPP1_myeloid_score = mean(spp1, na.rm = TRUE),
      CD44 = mean(cd44, na.rm = TRUE),
      ITGB1 = mean(itgb1, na.rm = TRUE),
      SPP1_CD44_expr_product = mean(spp1_expr * cd44, na.rm = TRUE),
      SPP1_ITGB1_expr_product = mean(spp1_expr * itgb1, na.rm = TRUE),
      target_subclone_02_04_score = mean(target, na.rm = TRUE),
      KRAS_hypoxia_score = mean(kras, na.rm = TRUE),
      SPP1_myeloid_score_high_fraction = mean(high_by_fraction(spp1, config$spatial_statistics$top_fraction)),
      target_subclone_02_04_score_high_fraction = mean(high_by_fraction(target, config$spatial_statistics$top_fraction)),
      KRAS_hypoxia_score_high_fraction = mean(high_by_fraction(kras, config$spatial_statistics$top_fraction)),
      score_score_spearman = unname(safe_cor(spp1, target, "spearman")["estimate"])
    )
  }
  sample_scores <- rbindlist(rows, fill = TRUE)
} else {
  sample_scores <- gse[, c(
    .(n_spots = .N, clinical_group = unique(clinical_group)[1]),
    lapply(.SD, mean, na.rm = TRUE),
    setNames(lapply(.SD, function(x) mean(high_by_fraction(x, config$spatial_statistics$top_fraction), na.rm = TRUE)),
             paste0(names(.SD), "_high_fraction"))
  ), by = sample_id, .SDcols = score_cols]
}

effect_rows <- list()
for (metric in setdiff(names(sample_scores), c("sample_id", "clinical_group", "n_spots"))) {
  ex <- sample_scores[tolower(clinical_group) == "excellent", get(metric)]
  poor <- sample_scores[tolower(clinical_group) == "poor", get(metric)]
  if (length(ex) > 0 && length(poor) > 0) {
    wt <- suppressWarnings(wilcox.test(ex, poor, exact = FALSE))
    all_pairs <- outer(poor, ex, "-")
    cliffs <- (sum(all_pairs > 0) - sum(all_pairs < 0)) / length(all_pairs)
    set.seed(as.integer(config$seed))
    boot <- replicate(1000, mean(sample(poor, replace = TRUE)) - mean(sample(ex, replace = TRUE)))
    effect_rows[[length(effect_rows) + 1]] <- data.table(
      metric = metric,
      excellent_median = median(ex, na.rm = TRUE),
      poor_median = median(poor, na.rm = TRUE),
      difference_poor_minus_excellent = median(poor, na.rm = TRUE) - median(ex, na.rm = TRUE),
      wilcoxon_p = wt$p.value,
      cliffs_delta_poor_vs_excellent = cliffs,
      bootstrap_ci_low = quantile(boot, 0.025, na.rm = TRUE),
      bootstrap_ci_high = quantile(boot, 0.975, na.rm = TRUE),
      n_excellent = length(ex),
      n_poor = length(poor)
    )
  }
}
stats <- rbindlist(effect_rows, fill = TRUE)
if (nrow(stats) > 0) stats[, fdr := p.adjust(wilcoxon_p, "BH")]

fwrite(sample_scores, file.path(out_dir, "sample_level_scores.csv"))
fwrite(stats, file.path(out_dir, "response_group_statistics.csv"))

plot_metric <- if ("SPP1_myeloid_score" %in% names(sample_scores)) "SPP1_myeloid_score" else score_cols[1]
p <- ggplot(sample_scores, aes(clinical_group, .data[[plot_metric]], color = clinical_group)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +
  geom_point(size = 2.2, position = position_jitter(width = 0.08, height = 0)) +
  labs(x = "Response group", y = plot_metric) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")
save_plot_both(p, file.path(root, "figures", "gse189843_response_dotbox"), width = 5, height = 4)
