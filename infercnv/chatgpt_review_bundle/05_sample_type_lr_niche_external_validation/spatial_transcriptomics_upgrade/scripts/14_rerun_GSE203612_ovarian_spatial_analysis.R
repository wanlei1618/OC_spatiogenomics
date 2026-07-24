options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(RANN)
  library(ggplot2)
})

local_root <- "D:/OC_spatiogenomics/spatial_data/spatial_analysis_correction_v2"
signatures <- fread(file.path(local_root, "updated_spatial_signatures.csv"))
objects <- readRDS(paste0(
  "D:/OC_spatiogenomics/spatial_data/processed/",
  "spatial_objects_curated_scored_reference_mapped.rds"
))
sample_ids <- c("GSM6177614", "GSM6177617")
stopifnot(all(sample_ids %in% names(objects)), !"GSM6177618" %in% sample_ids)

get_signature <- function(name) {
  signatures[signature_name == name, unique(gene)]
}
extract_counts <- function(object) {
  assay <- object@assays[[object@active.assay]]
  if ("counts" %in% methods::slotNames(assay)) {
    return(as(methods::slot(assay, "counts"), "dgCMatrix"))
  }
  x <- attr(assay, "layers")$counts
  dimnames(x) <- list(
    rownames(attr(assay, "features"))[seq_len(nrow(x))],
    rownames(attr(assay, "cells"))[seq_len(ncol(x))]
  )
  as(x, "dgCMatrix")
}
zscore <- function(v) {
  s <- sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(v)))
  (v - mean(v, na.rm = TRUE)) / s
}
score_signature <- function(norm, genes) {
  genes <- intersect(genes, rownames(norm))
  if (length(genes) < 2L) return(rep(NA_real_, ncol(norm)))
  z <- vapply(genes, function(g) zscore(as.numeric(norm[g, ])),
              numeric(ncol(norm)))
  rowMeans(z, na.rm = TRUE)
}
high_top <- function(v, fraction) {
  if (all(is.na(v))) return(rep(FALSE, length(v)))
  v >= quantile(v, 1 - fraction, na.rm = TRUE, names = FALSE)
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman",
                                 exact = FALSE))
  c(rho = unname(z$estimate), p = z$p.value, n = sum(ok))
}
hex_neighbors <- function(row, col) {
  key <- paste(row, col, sep = ":")
  offsets <- rbind(
    c(0, -2), c(0, 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1)
  )
  lapply(seq_along(row), function(i) {
    target <- paste(row[i] + offsets[, 1], col[i] + offsets[, 2], sep = ":")
    na.omit(match(target, key))
  })
}
knn_neighbors <- function(row, col, k) {
  xy <- cbind(row, col)
  idx <- RANN::nn2(xy, xy, k = min(k + 1L, nrow(xy)))$nn.idx[, -1,
    drop = FALSE]
  lapply(seq_len(nrow(idx)), function(i) idx[i, ])
}
evaluate_adjacency <- function(sender, receiver, neighbors, permutations = 1000L) {
  sender_idx <- which(sender)
  receiver <- as.logical(receiver)
  if (length(sender_idx) < 3L || sum(receiver) < 3L) {
    return(data.table(
      observed_neighbor_fraction = NA_real_,
      permuted_expected_fraction = NA_real_,
      observed_expected_ratio = NA_real_, log2_enrichment = NA_real_,
      empirical_p = NA_real_, n_sender_spots = length(sender_idx),
      n_receiver_spots = sum(receiver), n_edges = 0L
    ))
  }
  numerator <- vapply(neighbors, function(nb) sum(receiver[nb]), numeric(1))
  denominator <- lengths(neighbors)
  statistic <- function(idx) {
    d <- sum(denominator[idx])
    if (d == 0) return(NA_real_)
    sum(numerator[idx]) / d
  }
  observed <- statistic(sender_idx)
  null <- replicate(
    permutations,
    statistic(sample.int(length(sender), length(sender_idx), replace = FALSE))
  )
  expected <- mean(null, na.rm = TRUE)
  ratio <- if (is.finite(expected) && expected > 0) observed / expected else
    NA_real_
  data.table(
    observed_neighbor_fraction = observed,
    permuted_expected_fraction = expected,
    observed_expected_ratio = ratio,
    log2_enrichment = log2(ratio),
    empirical_p = (1 + sum(null >= observed, na.rm = TRUE)) /
      (1 + sum(is.finite(null))),
    n_sender_spots = length(sender_idx),
    n_receiver_spots = sum(receiver),
    n_edges = sum(denominator[sender_idx])
  )
}

thresholds <- c(top25 = .25, top15 = .15, top05 = .05)
k_values <- c(hex6 = 6L, knn4 = 4L, knn10 = 10L)
correlation_pairs <- list(
  SPP1_program_vs_ITGB1 = c("spp1_program", "ITGB1"),
  SPP1_program_vs_CD44 = c("spp1_program", "CD44"),
  SPP1_program_vs_malignant_epithelial =
    c("spp1_program", "malignant_epithelial"),
  SPP1_program_vs_Subclone02_04_KRAS_hypoxia =
    c("spp1_program", "subclone_kras_hypoxia"),
  C1QC_program_vs_ITGB1 = c("c1qc_program", "ITGB1"),
  C1QC_program_vs_CD44 = c("c1qc_program", "CD44"),
  macrophage_identity_vs_epithelial_identity =
    c("macrophage_identity", "epithelial_identity")
)

set.seed(20260729)
spot_list <- list()
correlations <- list()
neighborhood <- list()
negative <- list()
for (sid in sample_ids) {
  object <- objects[[sid]]
  counts <- extract_counts(object)
  meta <- as.data.table(object@meta.data, keep.rownames = "barcode")
  needed_coords <- c("array_row", "array_col", "coord_x", "coord_y")
  stopifnot(all(needed_coords %in% names(meta)))
  ids <- intersect(meta$barcode, colnames(counts))
  meta <- meta[match(ids, barcode)]
  counts <- counts[, ids, drop = FALSE]
  lib <- Matrix::colSums(counts)
  norm <- t(t(counts) / pmax(lib, 1)) * 1e4
  norm@x <- log1p(norm@x)
  gene_expr <- function(g) {
    if (g %in% rownames(norm)) as.numeric(norm[g, ]) else
      rep(NA_real_, ncol(norm))
  }
  score <- data.table(
    dataset_id = "GSE203612", sample_id = sid, barcode = ids,
    array_row = meta$array_row, array_col = meta$array_col,
    coord_x = meta$coord_x, coord_y = meta$coord_y,
    macrophage_identity = score_signature(
      norm, get_signature("macrophage_identity")
    ),
    spp1_program = score_signature(
      norm, get_signature("SPP1_macrophage_program")
    ),
    c1qc_program = score_signature(
      norm, get_signature("C1QC_macrophage_control")
    ),
    epithelial_identity = score_signature(
      norm, get_signature("epithelial_identity")
    ),
    malignant_epithelial = score_signature(
      norm, get_signature("CNV_supported_malignant_epithelial_signature")
    ),
    subclone_kras_hypoxia = score_signature(
      norm, get_signature("Subclone02_04_KRAS_hypoxia")
    ),
    ITGB1 = gene_expr("ITGB1"), CD44 = gene_expr("CD44")
  )
  spot_list[[sid]] <- score
  for (nm in names(correlation_pairs)) {
    pair <- correlation_pairs[[nm]]
    ct <- safe_cor(score[[pair[1]]], score[[pair[2]]])
    correlations[[length(correlations) + 1L]] <- data.table(
      dataset_id = "GSE203612", sample_id = sid,
      comparison = nm, spearman_rho = ct["rho"],
      descriptive_spot_level_p = ct["p"], n_spots = as.integer(ct["n"]),
      statistical_unit = "within-section spot; descriptive only"
    )
  }
  neighbor_sets <- list(
    hex6 = hex_neighbors(score$array_row, score$array_col),
    knn4 = knn_neighbors(score$array_row, score$array_col, 4L),
    knn10 = knn_neighbors(score$array_row, score$array_col, 10L)
  )
  for (threshold_name in names(thresholds)) {
    fraction <- thresholds[[threshold_name]]
    spp1_sender <- score$macrophage_identity >=
      median(score$macrophage_identity, na.rm = TRUE) &
      high_top(score$spp1_program, fraction)
    c1qc_sender <- score$macrophage_identity >=
      median(score$macrophage_identity, na.rm = TRUE) &
      high_top(score$c1qc_program, fraction)
    all_macro <- score$macrophage_identity >=
      median(score$macrophage_identity, na.rm = TRUE)
    epi_base <- score$epithelial_identity >=
      median(score$epithelial_identity, na.rm = TRUE)
    malignant_base <- score$malignant_epithelial >=
      median(score$malignant_epithelial, na.rm = TRUE)
    itgb1_receiver <- epi_base & malignant_base &
      high_top(score$ITGB1, fraction)
    cd44_receiver <- epi_base & malignant_base &
      high_top(score$CD44, fraction)
    malignant_high <- epi_base &
      high_top(score$malignant_epithelial, fraction)
    subclone_high <- epi_base &
      high_top(score$subclone_kras_hypoxia, fraction)
    tests <- list(
      SPP1_sender_to_ITGB1_receiver = c("spp1_sender", "itgb1_receiver"),
      SPP1_sender_to_CD44_receiver = c("spp1_sender", "cd44_receiver"),
      SPP1_sender_to_malignant_epithelial_high =
        c("spp1_sender", "malignant_high"),
      SPP1_sender_to_Subclone02_04_KRAS_hypoxia_high =
        c("spp1_sender", "subclone_high"),
      C1QC_sender_to_ITGB1_receiver = c("c1qc_sender", "itgb1_receiver"),
      C1QC_sender_to_CD44_receiver = c("c1qc_sender", "cd44_receiver"),
      all_macrophage_to_all_epithelial = c("all_macro", "epi_base")
    )
    flags <- list(
      spp1_sender = spp1_sender, c1qc_sender = c1qc_sender,
      all_macro = all_macro, itgb1_receiver = itgb1_receiver,
      cd44_receiver = cd44_receiver, malignant_high = malignant_high,
      subclone_high = subclone_high, epi_base = epi_base
    )
    for (neighbor_name in names(neighbor_sets)) {
      for (test_name in names(tests)) {
        pair <- tests[[test_name]]
        ans <- evaluate_adjacency(
          flags[[pair[1]]], flags[[pair[2]]],
          neighbor_sets[[neighbor_name]], 1000L
        )
        ans[, `:=`(
          dataset_id = "GSE203612", sample_id = sid,
          threshold = threshold_name, top_fraction = fraction,
          neighbor_definition = neighbor_name, comparison = test_name
        )]
        neighborhood[[length(neighborhood) + 1L]] <- ans
      }
    }
  }

  # One hundred average-expression-matched gene sets form a descriptive
  # sender negative-control distribution at the prespecified primary setting.
  mean_expr <- Matrix::rowMeans(norm)
  names(mean_expr) <- rownames(norm)
  candidate <- names(mean_expr)[
    is.finite(mean_expr) & !grepl("^(MT-|RPL|RPS)", names(mean_expr))
  ]
  bins <- pmin(
    10L,
    pmax(1L, ceiling(
      frank(mean_expr[candidate], ties.method = "average") /
        length(candidate) * 10
    ))
  )
  names(bins) <- candidate
  spp1_genes <- intersect(
    get_signature("SPP1_macrophage_program"), candidate
  )
  primary_neighbors <- hex_neighbors(score$array_row, score$array_col)
  for (b in seq_len(100L)) {
    random_genes <- vapply(spp1_genes, function(g) {
      pool <- names(bins)[bins == bins[g] & names(bins) != g]
      sample(pool, 1L)
    }, character(1))
    random_score <- score_signature(norm, random_genes)
    random_sender <- score$macrophage_identity >=
      median(score$macrophage_identity, na.rm = TRUE) &
      high_top(random_score, .25)
    epi_base <- score$epithelial_identity >=
      median(score$epithelial_identity, na.rm = TRUE)
    malignant_base <- score$malignant_epithelial >=
      median(score$malignant_epithelial, na.rm = TRUE)
    targets <- list(
      ITGB1_receiver = epi_base & malignant_base & high_top(score$ITGB1, .25),
      CD44_receiver = epi_base & malignant_base & high_top(score$CD44, .25)
    )
    for (target_name in names(targets)) {
      ans <- evaluate_adjacency(
        random_sender, targets[[target_name]], primary_neighbors, 1000L
      )
      ans[, `:=`(
        dataset_id = "GSE203612", sample_id = sid,
        random_set_id = b, target = target_name,
        matched_genes = paste(random_genes, collapse = ";"),
        threshold = "top25", neighbor_definition = "hex6"
      )]
      negative[[length(negative) + 1L]] <- ans
    }
  }
}

correlations <- rbindlist(correlations, fill = TRUE)
neighborhood <- rbindlist(neighborhood, fill = TRUE)
neighborhood[, BH_FDR := p.adjust(empirical_p, "BH"),
             by = .(sample_id, threshold, neighbor_definition)]
negative <- rbindlist(negative, fill = TRUE)
spots <- rbindlist(spot_list, fill = TRUE)
fwrite(correlations,
       file.path(local_root, "GSE203612_corrected_spatial_correlations.csv"),
       na = "NA")
fwrite(neighborhood,
       file.path(local_root, "GSE203612_corrected_neighborhood_results.csv"),
       na = "NA")
fwrite(negative,
       file.path(local_root, "GSE203612_negative_control_results.csv"),
       na = "NA")
fwrite(spots, file.path(local_root, "GSE203612_corrected_spot_scores.csv.gz"),
       compress = "gzip", na = "NA")

primary <- neighborhood[
  threshold == "top25" & neighbor_definition == "hex6" &
    comparison %chin% c("SPP1_sender_to_ITGB1_receiver",
                        "SPP1_sender_to_CD44_receiver")
]
conclusions <- primary[, {
  direction <- fifelse(observed_expected_ratio > 1, "ENRICHED",
                       "NOT_ENRICHED")
  .(
    comparison, observed_expected_ratio, empirical_p, BH_FDR,
    effect_direction = direction,
    within_section_support =
      observed_expected_ratio > 1 & empirical_p < .05
  )
}, by = .(dataset_id, sample_id)]
conclusions[, across_section_conclusion := {
  v <- observed_expected_ratio
  p <- empirical_p
  if (all(v > 1) && all(p < .05)) {
    "REPLICATED_IN_TWO_OVARIAN_SECTIONS"
  } else if (all(v > 1) && sum(p < .05) == 1L) {
    "DIRECTIONALLY_CONSISTENT_LIMITED_SUPPORT"
  } else if (any(v > 1) && any(v <= 1)) {
    "SPATIALLY_HETEROGENEOUS"
  } else {
    "NO_SPATIAL_ENRICHMENT"
  }
}, by = comparison]
fwrite(conclusions,
       file.path(local_root, "GSE203612_sample_level_conclusions.csv"),
       na = "NA")

spots[, map_class := fcase(
  macrophage_identity >= median(macrophage_identity) &
    high_top(spp1_program, .25), "SPP1 sender high",
  epithelial_identity >= median(epithelial_identity) &
    malignant_epithelial >= median(malignant_epithelial) &
    high_top(ITGB1, .25), "ITGB1 receiver high",
  epithelial_identity >= median(epithelial_identity) &
    malignant_epithelial >= median(malignant_epithelial) &
    high_top(CD44, .25), "CD44 receiver high",
  default = "other spot"
), by = sample_id]
p1 <- ggplot(spots, aes(coord_x, -coord_y, color = map_class)) +
  geom_point(size = .55) +
  facet_wrap(~sample_id) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Corrected ovarian spatial sender/receiver map", color = NULL)
ggsave(file.path(local_root, "GSE203612_corrected_sender_receiver_map.png"),
       p1, width = 10, height = 5, dpi = 180)

p2 <- ggplot(
  neighborhood[threshold == "top25" & neighbor_definition == "hex6"],
  aes(comparison, observed_expected_ratio, color = sample_id)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(position = position_dodge(width = .4), size = 2) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Corrected first-order adjacency effects",
    x = NULL, y = "observed / permuted expected", color = "ovarian section"
  )
ggsave(file.path(local_root, "GSE203612_corrected_adjacency_effects.png"),
       p2, width = 10, height = 6, dpi = 180)
message("Corrected GSE203612 ovarian spatial analysis complete")
