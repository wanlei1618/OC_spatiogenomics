suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

root_dir <- "D:/OC_spatiogenomics/infercnv"
script_file <- normalizePath("run_spp1_hypothesis_complete_bulk_validation.R", winslash = "/", mustWork = FALSE)
prev_dir <- file.path(root_dir, "SPP1_ITGB1_CD44_hypothesis_validation")
out_dir <- file.path(root_dir, "SPP1_ITGB1_CD44_hypothesis_validation_complete")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
script_dir <- file.path(out_dir, "scripts")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- list(
  TCGA_OV = list(
    expr = "D:/OC_Meta/TCGA_OV_expression/TCGA_counts_entrez.csv",
    clin = "D:/OC_Meta/TCGA_OV_expression/TCGA_OV_clinical.csv"
  ),
  GSE102094 = list(
    expr = "D:/OC_Meta/GSE102094/GSE102094_exp_entrez.csv",
    clin = "D:/OC_Meta/GSE102094/GSE102094_pd.csv"
  ),
  GSE32062 = list(
    expr = "D:/OC_Meta/GSE32062/GSE32062_exp_entrez.csv",
    clin = "D:/OC_Meta/GSE32062/GSE32062_GPL6480_pd.csv"
  ),
  GSE140082 = list(
    expr = "D:/OC_Meta/GSE140082/GSE140082_exp_entrez.csv",
    clin = "D:/OC_Meta/GSE140082/GSE140082_pd.csv"
  ),
  GSE49997 = list(
    expr = "D:/OC_Meta/GSE49997/GSE49997_exp_entrez.csv",
    clin = "D:/OC_Meta/GSE49997/GSE49997_pd.csv"
  )
)

signatures <- list(
  SPP1_TAM_score = c("SPP1", "CD68", "CD163", "CD14", "LST1", "TYROBP", "C1QA", "C1QB", "C1QC", "APOE", "MRC1", "MSR1", "FCGR3A", "ITGAM", "CSF1R"),
  ITGB1_CD44_tumor_score = c("ITGB1", "CD44"),
  KRAS_Hypoxia_score = c("KRAS", "ATF3", "EGR1", "FOS", "JUN", "DUSP6", "SPRY2", "MYC", "CXCL8", "PLAUR", "HIF1A", "CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "BNIP3", "NDRG1"),
  macrophage_fraction = c("CD68", "CD163", "CD14", "LST1", "TYROBP", "C1QA", "C1QB", "C1QC", "APOE", "MRC1", "MSR1", "FCGR3A", "ITGAM", "CSF1R"),
  epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "PAX8", "MUC16", "CLDN3", "CLDN4", "MSLN", "TACSTD2"),
  immune_leukocyte = c("PTPRC", "LYZ", "LST1", "TYROBP", "CD68", "C1QA", "C1QB", "C1QC")
)

clean_gene_id <- function(x) {
  x <- gsub('^"|"$', "", x)
  x <- trimws(x)
  x <- sub("\\.0$", "", x)
  x
}

normalize_id <- function(x) toupper(gsub("[[:space:]]+", "", as.character(x)))

read_expr_as_symbols <- function(path, cohort) {
  dt <- fread(path, check.names = FALSE, data.table = FALSE, showProgress = FALSE)
  header <- strsplit(readLines(path, n = 1, warn = FALSE), ",", fixed = TRUE)[[1]]
  header <- gsub('^"|"$', "", header)
  if (length(header) == ncol(dt)) names(dt) <- header
  if (names(dt)[1] == "" || is.na(names(dt)[1])) names(dt)[1] <- "sample_id"
  sample_id <- as.character(dt[[1]])
  raw_gene_ids <- clean_gene_id(names(dt)[-1])
  keep <- !is.na(raw_gene_ids) & raw_gene_ids != "" & !duplicated(seq_along(raw_gene_ids))
  expr <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(expr) <- "numeric"
  expr <- expr[, keep, drop = FALSE]
  raw_gene_ids <- raw_gene_ids[keep]
  rownames(expr) <- sample_id

  entrez <- raw_gene_ids[grepl("^[0-9]+$", raw_gene_ids)]
  mapped <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(entrez),
                                  keytype = "ENTREZID", columns = "SYMBOL")
  mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", , drop = FALSE]
  map_vec <- setNames(mapped$SYMBOL, mapped$ENTREZID)
  symbols <- unname(map_vec[raw_gene_ids])
  ensembl_like <- grepl("^ENSG", raw_gene_ids)
  fallback_idx <- is.na(symbols) & !ensembl_like & grepl("^[A-Za-z][A-Za-z0-9.-]+$", raw_gene_ids)
  if (any(fallback_idx)) symbols[fallback_idx] <- raw_gene_ids[fallback_idx]
  keep2 <- !is.na(symbols) & symbols != ""
  expr <- expr[, keep2, drop = FALSE]
  symbols <- symbols[keep2]

  finite_vals <- expr[is.finite(expr)]
  q99 <- as.numeric(quantile(finite_vals, 0.99, na.rm = TRUE))
  minv <- min(finite_vals, na.rm = TRUE)
  transform_used <- "none"
  if (minv >= 0 && q99 > 100) {
    lib <- rowSums(expr, na.rm = TRUE)
    lib[lib <= 0] <- median(lib[lib > 0], na.rm = TRUE)
    expr <- sweep(expr, 1, lib / 1e6, "/")
    expr <- log2(expr + 1)
    transform_used <- "log2(CPM+1)"
  } else if (minv >= 0 && q99 > 30) {
    expr <- log2(expr + 1)
    transform_used <- "log2(x+1)"
  }

  texpr <- t(expr)
  collapsed <- rowsum(texpr, group = symbols, reorder = FALSE)
  counts <- as.vector(table(factor(symbols, levels = rownames(collapsed))))
  collapsed <- collapsed / counts
  out <- t(collapsed)
  attr(out, "transform_used") <- transform_used
  out
}

score_signature <- function(expr, genes) {
  present <- intersect(genes, colnames(expr))
  if (length(present) == 0) return(rep(NA_real_, nrow(expr)))
  z <- scale(expr[, present, drop = FALSE])
  rowMeans(z, na.rm = TRUE)
}

prepare_clin <- function(cohort, path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (cohort == "TCGA_OV") {
    os_time <- ifelse(!is.na(x$days_to_death), x$days_to_death, x$days_to_last_followup) / 30.44
    os_event <- ifelse(tolower(x$vital_status) == "dead", 1, 0)
    out <- data.frame(sample_id = x$bcr_patient_barcode, OS_time = os_time, OS_event = os_event,
                      PFS_time = NA_real_, PFS_event = NA_real_,
                      stage = x$stage_event, grade = x$neoplasm_histologic_grade,
                      residual_disease = if ("residual_tumor" %in% names(x)) x$residual_tumor else x$tumor_residual_disease)
  } else if (cohort == "GSE102094") {
    out <- data.frame(sample_id = x$geo_accession, OS_time = as.numeric(x$os.mos), OS_event = as.numeric(x$osstatus),
                      PFS_time = as.numeric(x$pfs.mos), PFS_event = as.numeric(x$pfsstatus),
                      stage = x$stage, grade = NA_character_, residual_disease = x$surgical.outcome)
  } else if (cohort == "GSE32062") {
    out <- data.frame(sample_id = x$geo_accession, OS_time = as.numeric(x$os), OS_event = as.numeric(x$death),
                      PFS_time = as.numeric(x$pfs), PFS_event = as.numeric(x$rec),
                      stage = x$Stage, grade = x$grading, residual_disease = x$`surgery status`)
  } else if (cohort == "GSE140082") {
    os_time <- as.numeric(x$final_ostm)
    pfs_time <- as.numeric(x$final_pfstm)
    if (median(os_time, na.rm = TRUE) > 120) os_time <- os_time / 30.44
    if (median(pfs_time, na.rm = TRUE) > 120) pfs_time <- pfs_time / 30.44
    out <- data.frame(sample_id = x$geo_accession, OS_time = os_time, OS_event = as.numeric(x$final_osid),
                      PFS_time = pfs_time, PFS_event = as.numeric(x$final_pfsid),
                      stage = x$figo_stage, grade = x$newgrade, residual_disease = x$debulking_status)
  } else if (cohort == "GSE49997") {
    out <- data.frame(sample_id = x$geo_accession, OS_time = as.numeric(x$`os month`), OS_event = as.numeric(x$`os event`),
                      PFS_time = as.numeric(x$`pfs month`), PFS_event = as.numeric(x$`pfs event`),
                      stage = x$`figo grade`, grade = x$grade, residual_disease = x$`residual tumor`)
  }
  out$sample_key <- normalize_id(out$sample_id)
  out
}

scale_col <- function(x) as.numeric(scale(x))

fit_one <- function(dat, cohort, endpoint = "OS", adjusted = FALSE) {
  tcol <- paste0(endpoint, "_time")
  ecol <- paste0(endpoint, "_event")
  d <- dat[is.finite(dat[[tcol]]) & is.finite(dat[[ecol]]) & dat[[tcol]] > 0, , drop = FALSE]
  d <- d[complete.cases(d[, c("SPP1_TAM_score", "ITGB1_CD44_tumor_score", "KRAS_Hypoxia_score", "macrophage_fraction", "tumor_purity"), drop = FALSE]), , drop = FALSE]
  if (nrow(d) < 25 || sum(d[[ecol]] == 1, na.rm = TRUE) < 8) return(NULL)
  d$z_spp1_tam <- scale_col(d$SPP1_TAM_score)
  d$z_itgb1_cd44 <- scale_col(d$ITGB1_CD44_tumor_score)
  d$z_kras_hypoxia <- scale_col(d$KRAS_Hypoxia_score)
  d$z_macrophage <- scale_col(d$macrophage_fraction)
  d$z_tumor_purity <- scale_col(d$tumor_purity)

  rhs <- c("z_spp1_tam * z_itgb1_cd44", "z_kras_hypoxia", "z_macrophage", "z_tumor_purity")
  if (adjusted) {
    for (cv in c("stage", "grade", "residual_disease")) {
      vals <- d[[cv]]
      if (!all(is.na(vals)) && length(unique(vals[!is.na(vals)])) > 1) {
        d[[cv]] <- factor(vals)
        rhs <- c(rhs, cv)
      }
    }
  }
  form <- as.formula(paste0("Surv(", tcol, ", ", ecol, ") ~ ", paste(rhs, collapse = " + ")))
  fit <- tryCatch(coxph(form, data = d, singular.ok = TRUE), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(error = conditionMessage(fit), formula = deparse(form), n = nrow(d), events = sum(d[[ecol]] == 1)))
  }
  sm <- summary(fit)
  coef_df <- as.data.frame(sm$coefficients)
  ci_df <- as.data.frame(sm$conf.int)
  coef_df$term <- rownames(coef_df)
  ci_df$term <- rownames(ci_df)
  out <- merge(coef_df, ci_df[, c("term", "exp(coef)", "lower .95", "upper .95")], by = "term", all.x = TRUE)
  out$cohort <- cohort
  out$endpoint <- endpoint
  out$model <- ifelse(adjusted, "adjusted", "core")
  out$n <- nrow(d)
  out$events <- sum(d[[ecol]] == 1, na.rm = TRUE)
  form_string <- paste(deparse(form), collapse = " ")
  out$formula <- form_string
  list(fit = fit, coefficients = out, formula = form_string, n = nrow(d), events = sum(d[[ecol]] == 1, na.rm = TRUE))
}

make_scores <- function(cohort, info) {
  message("Reading ", cohort)
  expr <- read_expr_as_symbols(info$expr, cohort)
  sc <- data.frame(sample_id = rownames(expr), cohort = cohort, stringsAsFactors = FALSE)
  for (nm in names(signatures)[1:5]) sc[[nm]] <- score_signature(expr, signatures[[nm]])
  sc$immune_leukocyte <- score_signature(expr, signatures$immune_leukocyte)
  sc$tumor_purity <- sc$epithelial_tumor - sc$immune_leukocyte
  sc$sample_key <- normalize_id(sc$sample_id)
  clin <- prepare_clin(cohort, info$clin)
  merged <- merge(sc, clin, by = "sample_key", suffixes = c(".expr", ".clin"))
  merged$sample_id <- merged$sample_id.expr
  merged$transform_used <- attr(expr, "transform_used")
  present_genes <- lapply(signatures, function(g) intersect(g, colnames(expr)))
  present <- vapply(present_genes, function(g) paste(g, collapse = ";"), character(1))
  coverage <- data.frame(cohort = cohort, signature = names(present), present_genes = unname(present),
                         n_present = vapply(present_genes, length, integer(1)), stringsAsFactors = FALSE)
  fwrite(coverage, file.path(tab_dir, paste0(cohort, "_signature_gene_coverage.csv")))
  merged
}

all_scores <- rbindlist(Map(make_scores, names(cohorts), cohorts), fill = TRUE)
fwrite(all_scores, file.path(tab_dir, "bulk_signature_scores_merged.csv"))

fits <- list()
coef_list <- list()
status <- data.frame()
for (cohort_name in names(cohorts)) {
  d <- as.data.frame(all_scores)[as.data.frame(all_scores)$cohort == cohort_name, , drop = FALSE]
  for (endpoint in c("OS", "PFS")) {
    for (adjusted in c(FALSE, TRUE)) {
      nm <- paste(cohort_name, endpoint, ifelse(adjusted, "adjusted", "core"), sep = "_")
      res <- fit_one(d, cohort_name, endpoint, adjusted)
      if (is.null(res)) {
        status <- rbind(status, data.frame(cohort = cohort_name, endpoint = endpoint,
                                           model = ifelse(adjusted, "adjusted", "core"),
                                           status = "skipped_low_n_or_events", n = nrow(d), events = NA, formula = NA))
      } else if (!is.null(res$error)) {
        status <- rbind(status, data.frame(cohort = cohort_name, endpoint = endpoint,
                                           model = ifelse(adjusted, "adjusted", "core"),
                                           status = paste0("error: ", res$error), n = res$n, events = res$events, formula = res$formula))
      } else {
        fits[[nm]] <- res$fit
        coef_list[[nm]] <- res$coefficients
        status <- rbind(status, data.frame(cohort = cohort_name, endpoint = endpoint,
                                           model = ifelse(adjusted, "adjusted", "core"),
                                           status = "ok", n = res$n, events = res$events, formula = res$formula))
      }
    }
  }
}
coef_all <- rbindlist(coef_list, fill = TRUE)
fwrite(status, file.path(tab_dir, "bulk_cox_model_status.csv"))
fwrite(coef_all, file.path(tab_dir, "bulk_cox_model_coefficients.csv"))

interaction_terms <- coef_all[coef_all$term == "z_spp1_tam:z_itgb1_cd44", , drop = FALSE]
interaction_terms <- as.data.frame(interaction_terms)
hr_col <- if ("exp(coef).y" %in% names(interaction_terms)) "exp(coef).y" else if ("exp(coef)" %in% names(interaction_terms)) "exp(coef)" else "exp(coef).x"
interaction_terms$HR <- interaction_terms[[hr_col]]
interaction_terms$CI_low <- interaction_terms$`lower .95`
interaction_terms$CI_high <- interaction_terms$`upper .95`
interaction_terms$coef <- interaction_terms$coef
interaction_terms$se <- interaction_terms$`se(coef)`
interaction_terms$p_value <- interaction_terms$`Pr(>|z|)`
fwrite(interaction_terms, file.path(tab_dir, "bulk_cox_interaction_results.csv"))

meta_one <- function(x, model_name, endpoint_name) {
  x <- x[x$model == model_name & x$endpoint == endpoint_name & is.finite(x$coef) & is.finite(x$se) & x$se > 0, ]
  if (nrow(x) == 0) return(NULL)
  w <- 1 / (x$se ^ 2)
  beta <- sum(w * x$coef) / sum(w)
  se <- sqrt(1 / sum(w))
  z <- beta / se
  data.frame(endpoint = endpoint_name, model = model_name, n_cohorts = nrow(x),
             beta = beta, se = se, HR = exp(beta), CI_low = exp(beta - 1.96 * se),
             CI_high = exp(beta + 1.96 * se), p_value = 2 * pnorm(-abs(z)))
}
meta_res <- rbind(
  meta_one(interaction_terms, "core", "OS"),
  meta_one(interaction_terms, "adjusted", "OS"),
  meta_one(interaction_terms, "core", "PFS"),
  meta_one(interaction_terms, "adjusted", "PFS")
)
fwrite(meta_res, file.path(tab_dir, "bulk_meta_interaction_summary.csv"))

plot_forest <- function(dat, meta, endpoint = "OS", model = "adjusted", file_prefix = "bulk_os_adjusted") {
  x <- dat[dat$endpoint == endpoint & dat$model == model, , drop = FALSE]
  x <- x[is.finite(x$HR) & is.finite(x$CI_low) & is.finite(x$CI_high) & x$CI_low > 0 & x$CI_high > 0, , drop = FALSE]
  x <- x[order(x$cohort), ]
  if (nrow(x) == 0) return(invisible(NULL))
  m <- meta[meta$endpoint == endpoint & meta$model == model, , drop = FALSE]
  labels <- c(x$cohort, if (nrow(m)) "Fixed-effect meta" else NULL)
  hr <- c(x$HR, if (nrow(m)) m$HR else NULL)
  lo <- c(x$CI_low, if (nrow(m)) m$CI_low else NULL)
  hi <- c(x$CI_high, if (nrow(m)) m$CI_high else NULL)
  p <- c(x$p_value, if (nrow(m)) m$p_value else NULL)
  y <- rev(seq_along(labels))
  draw <- function() {
    par(mar = c(5, 8, 3, 2))
    plot(range(c(lo, hi), na.rm = TRUE), range(y) + c(-1, 1), type = "n", log = "x",
         yaxt = "n", xlab = "Hazard ratio for SPP1_TAM x ITGB1_CD44 interaction",
         ylab = "", main = paste(endpoint, model, "Cox interaction"))
    abline(v = 1, lty = 2, col = "grey50")
    segments(lo, y, hi, y, lwd = 2, col = "grey30")
    points(hr, y, pch = 19, cex = 1.1, col = ifelse(p < 0.05, "#B2182B", "#2166AC"))
    axis(2, at = y, labels = labels, las = 1)
    text(max(hi, na.rm = TRUE), y, labels = sprintf("HR %.2f [%.2f-%.2f], p=%.3g", hr, lo, hi, p),
         pos = 4, cex = 0.75, xpd = NA)
  }
  png(file.path(fig_dir, paste0(file_prefix, ".png")), width = 1800, height = 1100, res = 180)
  draw()
  dev.off()
  pdf(file.path(fig_dir, paste0(file_prefix, ".pdf")), width = 10, height = 6)
  draw()
  dev.off()
}
plot_forest(interaction_terms, meta_res, "OS", "adjusted", "bulk_OS_adjusted_interaction_forest")
plot_forest(interaction_terms, meta_res, "OS", "core", "bulk_OS_core_interaction_forest")
plot_forest(interaction_terms, meta_res, "PFS", "adjusted", "bulk_PFS_adjusted_interaction_forest")
plot_forest(interaction_terms, meta_res, "PFS", "core", "bulk_PFS_core_interaction_forest")

pkg_names <- c("survival", "edgeR", "limma", "copykat", "CaSpER", "liana", "nichenetr", "OmnipathR", "org.Hs.eg.db", "GSVA", "estimate", "MCPcounter", "survminer")
ip <- installed.packages()[, c("Package", "Version")]
pkg_status <- data.frame(package = pkg_names, installed = pkg_names %in% ip[, "Package"], version = NA_character_)
for (i in seq_along(pkg_names)) if (pkg_status$installed[i]) pkg_status$version[i] <- ip[match(pkg_names[i], ip[, "Package"]), "Version"]
fwrite(pkg_status, file.path(tab_dir, "R_package_status_for_complete_validation.csv"))

copy_previous <- function(rel) {
  src <- file.path(prev_dir, rel)
  dst <- file.path(out_dir, rel)
  if (file.exists(src)) {
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dst, overwrite = TRUE)
  }
}
for (rel in c(
  "tables/clone_malignant_identity_receptor_stability_summary.csv",
  "tables/sample_batch_stratified_receptor_program_scores.csv",
  "tables/focused_LR_predefined_axes_senders_to_Subclone02_04.csv",
  "tables/CellChat_subset_selected_sender_receiver_interactions.csv",
  "tables/pseudo_bulk_limma_focus_Subclone02_04_vs_others.csv",
  "tables/pseudo_bulk_focus_0204_GSEA_curated_pathways.csv",
  "tables/spatial_validation_marker_panel.csv",
  "tables/spatial_validation_scoring_strategy.csv",
  "figures/clone_tumor_identity_receptor_dotplot.png",
  "figures/sample_stratified_receptor_program_scores.png",
  "figures/SPP1_ITGB1_CD44_focused_LR_sender_receiver_heatmap.png",
  "figures/pseudo_bulk_focus_0204_GSEA_curated_pathways.png",
  "figures/pseudo_bulk_focus_0204_vs_others_volcano.png"
)) copy_previous(rel)

method_status <- data.frame(
  section = c("5.1", "5.2", "5.3", "5.4", "5.5"),
  analysis = c("malignant identity and receptor stability", "pseudo-bulk DEG/GSEA", "LR reproduction with CellChat/Connectome plus LIANA/NicheNet focused checks", "five-cohort Cox interaction model", "spatial validation marker/scoring plan"),
  status = c("completed from integrated_oc plus prior CNV clone metadata; CopyKAT sampled validation completed separately, CaSpER package load confirmed but full run blocked by local biomaRt SSL",
             "completed and copied from prior focused validation",
             "completed for Connectome-like and CellChat; LIANA consensus and NicheNet local gene-universe checks added separately",
             "completed in this run for TCGA_OV, GSE102094, GSE32062, GSE140082, GSE49997",
             "completed and copied from prior focused validation"),
  stringsAsFactors = FALSE
)
fwrite(method_status, file.path(tab_dir, "complete_validation_section_status.csv"))

report <- file.path(out_dir, "SPP1_ITGB1_CD44_complete_validation_report.md")
sink(report)
cat("# SPP1-ITGB1/CD44 complete validation update\n\n")
cat("Output directory: `", out_dir, "`\n\n", sep = "")
cat("## 5.4 Bulk Cox interaction model\n\n")
cat("Five cohorts were analyzed: TCGA_OV, GSE102094, GSE32062, GSE140082, and GSE49997. Expression matrices were read as sample-by-gene tables, Entrez IDs were mapped to gene symbols with org.Hs.eg.db, and cohort-wise z-score signatures were computed.\n\n")
cat("Primary interaction term: `SPP1_TAM_score x ITGB1_CD44_tumor_score`.\n\n")
cat("### Fixed-effect interaction summary\n\n")
print(meta_res)
cat("\n### Per-cohort interaction results\n\n")
print(interaction_terms[, c("cohort", "endpoint", "model", "n", "events", "HR", "CI_low", "CI_high", "p_value")])
cat("\n## 5.1-5.5 status\n\n")
print(method_status)
cat("\n## Key output tables\n\n")
cat("- `tables/bulk_signature_scores_merged.csv`\n")
cat("- `tables/bulk_cox_model_status.csv`\n")
cat("- `tables/bulk_cox_model_coefficients.csv`\n")
cat("- `tables/bulk_cox_interaction_results.csv`\n")
cat("- `tables/bulk_meta_interaction_summary.csv`\n")
cat("- `tables/R_package_status_for_complete_validation.csv`\n")
cat("- `tables/new_R_package_focused_validation_status.csv`\n")
cat("- `tables/CopyKAT_sampled_prediction_summary_by_known_group.csv`\n")
cat("- `tables/LIANA_consensus_focus_pair_presence.csv`\n")
cat("- `tables/NicheNet_focus_geneinfo_presence.csv`\n")
cat("- copied single-cell validation tables from the prior focused validation folder\n\n")
cat("## Key output figures\n\n")
cat("- `figures/bulk_OS_adjusted_interaction_forest.png`\n")
cat("- `figures/bulk_OS_core_interaction_forest.png`\n")
cat("- `figures/bulk_PFS_adjusted_interaction_forest.png`\n")
cat("- `figures/bulk_PFS_core_interaction_forest.png`\n")
cat("- copied single-cell validation figures from the prior focused validation folder\n")
sink()

file.copy(script_file, file.path(script_dir, "run_spp1_hypothesis_complete_bulk_validation.R"), overwrite = TRUE)
