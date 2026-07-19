#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
  if (is.na(file_arg)) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = FALSE))
}

load_config <- function(path = NULL) {
  if (is.null(path)) path <- file.path(script_path(), "..", "config", "spatial_config.yml")
  yaml::read_yaml(path)
}

path_root <- function(config) normalizePath(config$project_root, winslash = "/", mustWork = FALSE)

ensure_dirs <- function(root) {
  dirs <- file.path(root, c(
    "processed", "results/spatial_curated", "results/gse189843_response",
    "results/meta", "reference_mapping", "figures", "reports", "logs", "tmp"
  ))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

read_dt_if_exists <- function(path) {
  if (!file.exists(path)) return(data.table())
  data.table::fread(path)
}

write_empty_reason <- function(path, columns, reason) {
  dt <- as.data.table(setNames(rep(list(character()), length(columns)), columns))
  attr(dt, "reason") <- reason
  data.table::fwrite(dt, path)
}

read_spot_scores <- function(root) {
  path <- file.path(root, "results", "spatial_curated", "spatial_spot_scores_curated.csv.gz")
  if (!file.exists(path)) stop("Missing spot score table: ", path)
  dt <- data.table::fread(path)
  data.table::setnames(dt, old = intersect("Target_Subclone02_04_score", names(dt)),
                       new = "target_subclone_02_04_score")
  dt
}

read_neighborhood <- function(root) {
  path <- file.path(root, "results", "spatial_curated", "spatial_neighborhood_enrichment_curated.csv")
  read_dt_if_exists(path)
}

read_reference_predictions <- function(root) {
  candidates <- c(
    file.path(root, "results", "reference_mapping", "spatial_reference_mapping_predictions.csv.gz"),
    file.path(root, "reference_mapping", "reference_mapping_predictions.csv.gz")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) return(data.table())
  data.table::fread(hit)
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(c(estimate = NA_real_, p = NA_real_))
  out <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  c(estimate = unname(out$estimate), p = out$p.value)
}

high_by_fraction <- function(x, fraction) {
  if (all(!is.finite(x))) return(rep(FALSE, length(x)))
  cutoff <- as.numeric(stats::quantile(x, probs = 1 - fraction, na.rm = TRUE, names = FALSE))
  x >= cutoff
}

knn_index_from_xy <- function(x, y, k) {
  coords <- cbind(as.numeric(x), as.numeric(y))
  ok <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  if (!all(ok)) stop("Coordinates contain non-finite values")
  if (nrow(coords) <= k) stop("Too few spots for kNN")
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(coords, k = min(k + 1, nrow(coords)))$nn.idx
    return(nn[, -1, drop = FALSE])
  }
  d <- as.matrix(dist(coords))
  t(apply(d, 1, function(v) order(v)[2:(k + 1)]))
}

neighbor_test <- function(source, target, nn, n_perm = 1000, seed = 1) {
  source_idx <- which(source)
  if (length(source_idx) == 0 || sum(target) == 0) {
    return(data.table(observed = NA_real_, expected = mean(target), ratio = NA_real_,
                      log2_enrichment = NA_real_, empirical_p = NA_real_,
                      n_source = sum(source), n_target = sum(target), n_edges = 0L))
  }
  observed <- mean(target[nn[source_idx, , drop = FALSE]])
  set.seed(seed)
  perm <- replicate(n_perm, {
    perm_target <- sample(target, length(target), replace = FALSE)
    mean(perm_target[nn[source_idx, , drop = FALSE]])
  })
  expected <- mean(perm)
  ratio <- ifelse(expected > 0, observed / expected, NA_real_)
  data.table(
    observed = observed,
    expected = expected,
    ratio = ratio,
    log2_enrichment = log2(ratio),
    empirical_p = (sum(perm >= observed, na.rm = TRUE) + 1) / (length(perm) + 1),
    n_source = sum(source),
    n_target = sum(target),
    n_edges = length(source_idx) * ncol(nn)
  )
}

moran_geary <- function(values, nn, n_perm = 1000, seed = 1) {
  x <- as.numeric(values)
  ok <- is.finite(x)
  if (sum(ok) < 5) return(data.table(moran_i = NA_real_, moran_p = NA_real_, geary_c = NA_real_, geary_p = NA_real_))
  x[!ok] <- mean(x[ok])
  n <- length(x)
  x_center <- x - mean(x)
  denom <- sum(x_center^2)
  edges_i <- rep(seq_len(nrow(nn)), ncol(nn))
  edges_j <- as.vector(nn)
  w <- length(edges_i)
  calc <- function(v) {
    vc <- v - mean(v)
    d <- sum(vc^2)
    if (d == 0) return(c(I = NA_real_, C = NA_real_))
    I <- (n / w) * sum(vc[edges_i] * vc[edges_j]) / d
    C <- ((n - 1) / (2 * w)) * sum((v[edges_i] - v[edges_j])^2) / d
    c(I = I, C = C)
  }
  obs <- calc(x)
  set.seed(seed)
  perm <- replicate(n_perm, calc(sample(x, length(x), replace = FALSE)))
  data.table(
    moran_i = obs["I"],
    moran_p = (sum(perm["I", ] >= obs["I"], na.rm = TRUE) + 1) / (n_perm + 1),
    geary_c = obs["C"],
    geary_p = (sum(perm["C", ] <= obs["C"], na.rm = TRUE) + 1) / (n_perm + 1)
  )
}

save_plot_both <- function(plot, stem, width = 8, height = 5) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = width, height = height)
}

write_session_info <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(capture.output(sessionInfo())), path)
}
