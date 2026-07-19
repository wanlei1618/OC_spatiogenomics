suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(yaml)
})

Sys.setenv(TMPDIR = "D:/OC_spatiogenomics/tmp")
Sys.setenv(TEMP = "D:/OC_spatiogenomics/tmp")
Sys.setenv(TMP = "D:/OC_spatiogenomics/tmp")
dir.create("D:/OC_spatiogenomics/tmp", recursive = TRUE, showWarnings = FALSE)

config_path <- "D:/OC_spatiogenomics/infercnv/sample_type_LR_niche_analysis/config.yaml"
cfg <- yaml::read_yaml(config_path)
project_root <- cfg$project_root
tab_dir <- file.path(project_root, "tables")
fig_dir <- file.path(project_root, "figures")
obj_dir <- file.path(project_root, "objects")
log_dir <- file.path(project_root, "logs")
ext_dir <- file.path(project_root, "external_scRNA")
spatial_dir <- file.path(project_root, "spatial")
ko_dir <- file.path(project_root, "virtual_KO")
for (d in c(tab_dir, fig_dir, obj_dir, log_dir, ext_dir, spatial_dir, ko_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "run_sample_type_LR_niche_analysis.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

message("Started at ", Sys.time())

save_table <- function(x, stem) {
  fwrite(as.data.frame(x), file.path(tab_dir, paste0(stem, ".csv")))
  saveRDS(x, file.path(tab_dir, paste0(stem, ".rds")))
}

save_plot <- function(plot, stem, width = 9, height = 6, dpi = 220) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi)
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
}

read_seurat_object <- function(path) {
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(obj)) return(obj)
  env <- new.env()
  load(path, envir = env)
  env[[ls(env)[1]]]
}

std_sample_id <- function(meta) {
  if ("sample_id" %in% names(meta)) return(as.character(meta$sample_id))
  if ("sample" %in% names(meta)) return(as.character(meta$sample))
  if ("batch" %in% names(meta)) return(as.character(meta$batch))
  if ("orig.ident" %in% names(meta)) return(as.character(meta$orig.ident))
  rep("unknown_sample", nrow(meta))
}

axis_name <- function(ligand, receptor) paste(ligand, receptor, sep = "-")
axes <- rbind(
  data.frame(ligand = c("SPP1", "SPP1"), receptor = c("CD44", "ITGB1"), axis_class = "primary"),
  data.frame(ligand = c("MIF", "MIF", "APOE", "TGFB1", "TGFB1", "CXCL12"),
             receptor = c("CD74", "CXCR4", "LRP1", "TGFBR1", "TGFBR2", "CXCR4"),
             axis_class = "control")
)
axes$axis <- axis_name(axes$ligand, axes$receptor)
target_clones <- unlist(cfg$target_clones)
all_clones <- unlist(cfg$all_clones)
target_groups <- paste0("CNV_", all_clones)
focus_target_groups <- paste0("CNV_", target_clones)
min_cells <- as.integer(cfg$min_cells_group)

for (p in unlist(cfg$existing_results)) {
  if (!file.exists(p)) stop("Required input missing: ", p)
}
if (!file.exists(cfg$integrated_oc_rds)) stop("Integrated OC object missing: ", cfg$integrated_oc_rds)

metadata <- fread(cfg$existing_results$metadata, data.table = FALSE)
metadata$cell_id <- metadata$cell_integrated_oc
metadata$sample_id <- std_sample_id(metadata)
metadata$sample_type <- if ("sample_type" %in% names(metadata)) as.character(metadata$sample_type) else "unknown_sample_type"
metadata$cnv_subclone <- as.character(metadata$cnv_subclone)
metadata$target_group <- ifelse(!is.na(metadata$cnv_subclone) & metadata$cnv_subclone != "", paste0("CNV_", metadata$cnv_subclone), NA)
metadata$interaction_group <- as.character(metadata$interaction_group)
metadata$source_is_myeloid <- grepl("Myeloid|Macro|Monocyte|DC|cDC|mDC", metadata$interaction_group, ignore.case = TRUE) |
  grepl("Macrophage|Monocyte|DC", metadata$cell_type, ignore.case = TRUE)
metadata$target_is_cnv <- !is.na(metadata$target_group) & metadata$target_group %in% target_groups

save_table(metadata, "metadata_cleaned_for_sample_type_analysis")
sample_id_st <- as.data.frame.matrix(table(metadata$sample_id, metadata$sample_type))
sample_id_st$sample_id <- rownames(sample_id_st)
sample_id_st <- sample_id_st[, c("sample_id", setdiff(names(sample_id_st), "sample_id"))]
save_table(sample_id_st, "sample_id_by_sample_type_cell_counts")

sample_type_sample_id <- as.data.frame(table(sample_id = metadata$sample_id, sample_type = metadata$sample_type))
sample_type_sample_id <- sample_type_sample_id[sample_type_sample_id$Freq > 0, ]
save_table(sample_type_sample_id, "sample_id_sample_type_long_counts")
confounded <- all(rowSums(as.matrix(sample_id_st[, setdiff(names(sample_id_st), "sample_id"), drop = FALSE]) > 0) == 1)
confound_diag <- data.frame(
  diagnostic = c("sample_id_one_to_one_with_sample_type", "n_sample_ids", "n_sample_types"),
  value = c(as.character(confounded), length(unique(metadata$sample_id)), length(unique(metadata$sample_type)))
)
save_table(confound_diag, "sample_type_sample_id_confounding_diagnostic")

cnv_counts <- as.data.frame(table(sample_type = metadata$sample_type, cnv_subclone = metadata$cnv_subclone), stringsAsFactors = FALSE)
cnv_counts <- cnv_counts[!is.na(cnv_counts$cnv_subclone) & cnv_counts$cnv_subclone != "", ]
cnv_counts$total_in_sample_type <- ave(cnv_counts$Freq, cnv_counts$sample_type, FUN = sum)
cnv_counts$fraction_in_sample_type <- cnv_counts$Freq / cnv_counts$total_in_sample_type
save_table(cnv_counts, "sample_type_by_cnv_subclone_counts")

ig_counts <- as.data.frame(table(sample_type = metadata$sample_type, interaction_group = metadata$interaction_group), stringsAsFactors = FALSE)
ig_counts <- ig_counts[!is.na(ig_counts$interaction_group) & ig_counts$interaction_group != "", ]
ig_counts$total_in_sample_type <- ave(ig_counts$Freq, ig_counts$sample_type, FUN = sum)
ig_counts$fraction_in_sample_type <- ig_counts$Freq / ig_counts$total_in_sample_type
save_table(ig_counts, "sample_type_by_interaction_group_counts")

sample_type_myeloid_counts <- ig_counts[grepl("Myeloid|Macro|Monocyte|DC|cDC|mDC", ig_counts$interaction_group, ignore.case = TRUE), ]
save_table(sample_type_myeloid_counts, "sample_type_by_myeloid_source_group_counts")

p_cnv <- ggplot(cnv_counts, aes(x = sample_type, y = fraction_in_sample_type, fill = cnv_subclone)) +
  geom_col(color = "white", size = 0.1) +
  labs(x = "sample_type", y = "fraction among CNV tumor cells", fill = "CNV subclone") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_cnv, "sample_type_by_cnv_subclone_fraction_barplot", 8, 5)

ig_top <- ig_counts[order(-ig_counts$Freq), ]
top_ig <- unique(ig_top$interaction_group)[seq_len(min(20, length(unique(ig_top$interaction_group))))]
p_ig <- ggplot(ig_counts[ig_counts$interaction_group %in% top_ig, ],
               aes(x = sample_type, y = fraction_in_sample_type, fill = interaction_group)) +
  geom_col(color = "white", size = 0.05) +
  labs(x = "sample_type", y = "fraction among all cells", fill = "interaction_group") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
save_plot(p_ig, "sample_type_by_interaction_group_fraction_barplot", 10, 6)

obj <- read_seurat_object(cfg$integrated_oc_rds)
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
cells <- intersect(colnames(expr), metadata$cell_id)
metadata <- metadata[match(cells, metadata$cell_id), ]
expr <- expr[, cells, drop = FALSE]

gene_expr_stats <- function(cell_set, genes) {
  cell_set <- intersect(cell_set[!is.na(cell_set)], colnames(expr))
  present <- intersect(genes, rownames(expr))
  out <- data.frame(gene = genes, avg = NA_real_, pct = NA_real_, present = genes %in% rownames(expr))
  if (length(cell_set) == 0 || length(present) == 0) return(out)
  for (g in present) {
    v <- as.numeric(expr[g, cell_set, drop = TRUE])
    out$avg[out$gene == g] <- mean(v, na.rm = TRUE)
    out$pct[out$gene == g] <- mean(v > 0, na.rm = TRUE)
  }
  out
}

compute_lr_by_group <- function(meta, group_var) {
  group_values <- sort(unique(meta[[group_var]][!is.na(meta[[group_var]]) & meta[[group_var]] != ""]))
  rows <- list()
  idx <- 1
  for (gv in group_values) {
    m0 <- meta[meta[[group_var]] == gv, ]
    denom <- nrow(m0)
    source_groups <- sort(unique(m0$interaction_group[m0$source_is_myeloid & !is.na(m0$interaction_group)]))
    for (sg in source_groups) {
      s_idx <- !is.na(m0$interaction_group) & m0$interaction_group == sg
      s_cells <- m0$cell_id[s_idx]
      if (length(s_cells) < min_cells) next
      s_frac <- length(s_cells) / denom
      for (tg in target_groups) {
        t_idx <- !is.na(m0$target_group) & m0$target_group == tg
        t_cells <- m0$cell_id[t_idx]
        t_frac <- length(t_cells) / denom
        for (i in seq_len(nrow(axes))) {
          ligand <- axes$ligand[i]
          receptor <- axes$receptor[i]
          l_stat <- gene_expr_stats(s_cells, ligand)
          r_stat <- gene_expr_stats(t_cells, receptor)
          expr_product <- l_stat$avg[1] * r_stat$avg[1]
          pct_weighted <- l_stat$pct[1] * r_stat$pct[1] * expr_product
          abundance <- s_frac * t_frac * expr_product
          axis_score <- s_frac * t_frac * l_stat$avg[1] * r_stat$avg[1] * l_stat$pct[1] * r_stat$pct[1]
          rows[[idx]] <- data.frame(
            level = group_var, group = gv, sample_type = ifelse(group_var == "sample_type", gv, unique(m0$sample_type)[1]),
            sample_id = ifelse(group_var == "sample_id", gv, NA_character_),
            source_group = sg, target_group = tg,
            target_clone = sub("^CNV_", "", tg),
            ligand = ligand, receptor = receptor, axis = axes$axis[i], axis_class = axes$axis_class[i],
            source_n = length(s_cells), target_n = length(t_cells), total_n = denom,
            source_fraction = s_frac, target_fraction = t_frac,
            ligand_avg_source = l_stat$avg[1], receptor_avg_target = r_stat$avg[1],
            ligand_pct_source = l_stat$pct[1], receptor_pct_target = r_stat$pct[1],
            expr_product_score = expr_product,
            abundance_weighted_score = abundance,
            pct_weighted_score = pct_weighted,
            axis_score = axis_score,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        }
      }
    }
  }
  rbindlist(rows, fill = TRUE)
}

sample_type_lr <- compute_lr_by_group(metadata, "sample_type")
sample_id_lr <- compute_lr_by_group(metadata, "sample_id")
save_table(sample_type_lr, "sample_type_LR_opportunity_scores_all_axes")
save_table(sample_type_lr[sample_type_lr$axis_class == "primary", ], "sample_type_LR_opportunity_scores_primary_axes")
save_table(sample_type_lr[sample_type_lr$target_group %in% focus_target_groups, ], "sample_type_LR_opportunity_scores_target_Subclone02_04")
save_table(sample_id_lr, "sample_id_LR_opportunity_scores_all_axes")

primary_focus <- sample_type_lr[sample_type_lr$axis %in% c("SPP1-CD44", "SPP1-ITGB1") & sample_type_lr$target_group %in% focus_target_groups, ]
save_table(primary_focus, "sample_type_LR_opportunity_scores_primary_axes_target_Subclone02_04")

heat_df <- sample_type_lr[sample_type_lr$target_group %in% focus_target_groups, ]
heat_sum <- heat_df[, .(axis_score = sum(axis_score, na.rm = TRUE),
                        abundance_weighted_score = sum(abundance_weighted_score, na.rm = TRUE)),
                    by = .(sample_type, axis)]
heat_mat <- dcast(heat_sum, sample_type ~ axis, value.var = "axis_score", fill = 0)
hm <- as.matrix(heat_mat[, -1, drop = FALSE])
rownames(hm) <- heat_mat$sample_type
png(file.path(fig_dir, "heatmap_sample_type_by_axis_TargetSubclone02_04.png"), width = 1600, height = 1000, res = 180)
pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column",
         main = "Target Subclone_02/04 LR opportunity score by sample_type")
dev.off()
pdf(file.path(fig_dir, "heatmap_sample_type_by_axis_TargetSubclone02_04.pdf"), width = 9, height = 5.5)
pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column",
         main = "Target Subclone_02/04 LR opportunity score by sample_type")
dev.off()

bubble <- primary_focus
bubble$receptor <- factor(bubble$receptor, levels = c("CD44", "ITGB1"))
p_bubble <- ggplot(bubble, aes(x = target_clone, y = source_group, size = axis_score, color = ligand_avg_source)) +
  geom_point(alpha = 0.85) +
  facet_grid(sample_type ~ receptor, scales = "free_y", space = "free_y") +
  scale_size_continuous(range = c(0.2, 7)) +
  scale_color_gradient(low = "#4575B4", high = "#D73027") +
  labs(x = "target clone", y = "myeloid source group", size = "axis_score", color = "SPP1 avg") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_bubble, "bubble_source_to_target_by_sample_type_SPP1_CD44_ITGB1", 10, 9)

bar_data <- heat_sum[heat_sum$axis %in% c("SPP1-CD44", "SPP1-ITGB1"), ]
p_axis <- ggplot(bar_data, aes(x = sample_type, y = axis_score, fill = axis)) +
  geom_col(position = "dodge", color = "white") +
  labs(x = "sample_type", y = "summed axis_score to Subclone_02/04", fill = "axis") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_axis, "boxplot_axis_score_by_sample_type_external_ready", 7, 5)

clone_frac <- cnv_counts
p_stack <- ggplot(clone_frac, aes(x = sample_type, y = fraction_in_sample_type, fill = cnv_subclone)) +
  geom_col(color = "white") +
  labs(x = "sample_type", y = "target clone fraction", fill = "CNV subclone") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_stack, "stacked_bar_target_clone_fraction_by_sample_type", 8, 5)

scatter_dat <- sample_type_lr[sample_type_lr$target_group %in% focus_target_groups & sample_type_lr$axis %in% axes$axis, ]
p_scatter <- ggplot(scatter_dat, aes(x = target_fraction, y = axis_score, color = sample_type)) +
  geom_point(alpha = 0.75) +
  facet_wrap(~ axis, scales = "free_y") +
  labs(x = "target_fraction", y = "axis_score", color = "sample_type") +
  theme_bw(base_size = 10)
save_plot(p_scatter, "scatter_target_abundance_vs_axis_score_by_sample_type", 10, 6)

deg <- fread(cfg$existing_results$deg_markers, data.table = FALSE)
top_markers <- fread(cfg$existing_results$top_markers, data.table = FALSE)
name_cols <- names(deg)
gene_col <- intersect(c("gene", "Gene", "features", "feature"), name_cols)[1]
clone_col <- intersect(c("CNV_clone", "cnv_subclone", "cluster", "clone"), name_cols)[1]
fc_col <- intersect(c("avg_log2FC", "avg_logFC", "logFC"), name_cols)[1]
p_col <- intersect(c("p_val_adj", "adj.P.Val", "p_val"), name_cols)[1]
if (is.na(gene_col) || is.na(clone_col)) {
  gene_col <- names(deg)[1]
  clone_col <- names(deg)[2]
}
deg$gene_clean <- as.character(deg[[gene_col]])
deg$clone_clean <- as.character(deg[[clone_col]])
deg$fc_clean <- if (!is.na(fc_col)) as.numeric(deg[[fc_col]]) else 0
deg$p_clean <- if (!is.na(p_col)) as.numeric(deg[[p_col]]) else 1
exclude_pattern <- "^MT-|^RPL|^RPS|^MTRNR|^IG[HKL]|^TR[ABDG]"
marker_sig <- function(clone, n = 80) {
  x <- deg[deg$clone_clean %in% c(clone, paste0("CNV_", clone)) & !grepl(exclude_pattern, deg$gene_clean), ]
  x <- x[order(-x$fc_clean, x$p_clean), ]
  unique(head(x$gene_clean, n))
}
sig02 <- marker_sig("Subclone_02")
sig04 <- marker_sig("Subclone_04")
sig_common <- unique(c(intersect(sig02, sig04),
                       c("EPCAM", "KRT8", "KRT18", "KRT19", "PAX8", "MUC16", "CLDN3", "CLDN4", "KRAS", "HIF1A", "VEGFA", "CD44", "ITGB1")))
sig_common <- sig_common[sig_common %in% rownames(expr)]
sig_cd44_itgb1 <- unique(c("CD44", "ITGB1", "ITGA5", "ITGAV", "ITGA6", "ITGB4", "FN1", "COL1A1", "COL1A2",
                           "COL3A1", "COL5A1", "LAMC2", "LAMB3", "VCL", "PXN", "PTK2", "SRC", "VIM", "MMP2", "MMP14", "SERPINE1"))
sig_kras_hypoxia <- unique(c("KRAS", "ATF3", "EGR1", "FOS", "JUN", "DUSP6", "SPRY2", "MYC", "CXCL8", "PLAUR",
                             "HIF1A", "CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "BNIP3", "NDRG1"))
sig_epithelial <- unique(c("EPCAM", "KRT8", "KRT18", "KRT19", "PAX8", "MUC16", "CLDN3", "CLDN4", "MSLN", "TACSTD2"))
sig_list <- list(
  Subclone02_like = sig02,
  Subclone04_like = sig04,
  Subclone02_04_common = sig_common,
  CD44_ITGB1_target = sig_cd44_itgb1[sig_cd44_itgb1 %in% rownames(expr)],
  KRAS_hypoxia_target = sig_kras_hypoxia[sig_kras_hypoxia %in% rownames(expr)],
  tumor_epithelial = sig_epithelial[sig_epithelial %in% rownames(expr)]
)
save_signature <- function(vec, stem) save_table(data.frame(gene = vec), stem)
save_signature(sig_list$Subclone02_like, "signature_Subclone02_like")
save_signature(sig_list$Subclone04_like, "signature_Subclone04_like")
save_signature(sig_list$Subclone02_04_common, "signature_Subclone02_04_common")
save_signature(sig_list$CD44_ITGB1_target, "signature_CD44_ITGB1_target")
save_signature(sig_list$KRAS_hypoxia_target, "signature_KRAS_hypoxia_target")
save_signature(sig_list$tumor_epithelial, "signature_tumor_epithelial")
saveRDS(sig_list, file.path(obj_dir, "signature_list.rds"))

external_files <- list.files(ext_dir, recursive = TRUE, full.names = TRUE)
external_status <- data.frame(dataset = character(), status = character(), note = character())
if (length(external_files) == 0) {
  external_status <- data.frame(dataset = "none", status = "no_external_scRNA_input",
                                note = "No external scRNA files were found. Template tables were generated only.")
  save_table(data.frame(), "external_scRNA_all_datasets_sampletype_axis_scores")
  save_table(data.frame(), "external_scRNA_patient_sampletype_summary")
  save_table(data.frame(), "external_scRNA_target_signature_scores")
  save_table(data.frame(), "external_scRNA_SPP1_myeloid_CD44_ITGB1_target_scores")
  for (stem in c("external_boxplot_SPP1_CD44_by_sample_type", "external_boxplot_SPP1_ITGB1_by_sample_type",
                 "external_paired_tumor_vs_ascites_axis_score", "external_heatmap_dataset_sampletype_axis")) {
    p_empty <- ggplot() + annotate("text", x = 0, y = 0, label = "No external scRNA input provided") + theme_void()
    save_plot(p_empty, stem, 6, 4)
  }
}
save_table(external_status, "external_scRNA_validation_status")

spatial_files <- list.files(spatial_dir, recursive = TRUE, full.names = TRUE)
spatial_status <- data.frame(dataset = character(), status = character(), note = character())
if (length(spatial_files) == 0) {
  spatial_status <- data.frame(dataset = "none", status = "no_spatial_input",
                               note = "No spatial transcriptomics files were found. Template tables were generated only.")
  save_table(data.frame(), "spatial_correlation_SPP1_myeloid_Target_axis")
  save_table(data.frame(), "spatial_neighborhood_enrichment")
  save_table(data.frame(), "spatial_LR_COMMOT_scores")
  save_table(data.frame(), "spatial_LR_COMMOT_scores_by_region")
  save_table(data.frame(), "spatial_nearest_neighbor_distances")
  save_table(data.frame(), "spatial_permutation_test_SPP1_myeloid_to_target_tumor")
  for (stem in c("spatial_neighborhood_enrichment_heatmap", "spatial_COMMOT_SPP1_CD44_map",
                 "spatial_COMMOT_SPP1_ITGB1_map", "spatial_COMMOT_axis_comparison_heatmap")) {
    p_empty <- ggplot() + annotate("text", x = 0, y = 0, label = "No spatial input provided") + theme_void()
    save_plot(p_empty, stem, 6, 4)
  }
}
save_table(spatial_status, "spatial_validation_status")

ko_base <- sample_type_lr[sample_type_lr$axis %in% c("SPP1-CD44", "SPP1-ITGB1") & sample_type_lr$target_group %in% focus_target_groups, ]
ko_rows <- list()
idx <- 1
for (i in seq_len(nrow(ko_base))) {
  row <- ko_base[i, ]
  for (scenario in c("Control", "SPP1_KO", "CD44_KO", "ITGB1_KO", "CD44_ITGB1_double_KO")) {
    score <- row$axis_score
    if (scenario == "SPP1_KO") score <- 0
    if (scenario == "CD44_KO" && row$receptor == "CD44") score <- 0
    if (scenario == "ITGB1_KO" && row$receptor == "ITGB1") score <- 0
    if (scenario == "CD44_ITGB1_double_KO") score <- 0
    ko_rows[[idx]] <- cbind(row[, c("sample_type", "source_group", "target_group", "axis", "ligand", "receptor", "axis_score"), drop = FALSE],
                            scenario = scenario, perturbed_axis_score = score)
    idx <- idx + 1
  }
}
ko <- rbindlist(ko_rows, fill = TRUE)
control <- ko[ko$scenario == "Control", c("sample_type", "source_group", "target_group", "axis", "axis_score")]
names(control)[5] <- "control_axis_score"
ko <- merge(ko, control, by = c("sample_type", "source_group", "target_group", "axis"), all.x = TRUE)
ko$absolute_reduction <- ko$control_axis_score - ko$perturbed_axis_score
ko$relative_reduction <- ifelse(ko$control_axis_score > 0, ko$absolute_reduction / ko$control_axis_score, NA_real_)
save_table(ko, "virtual_KO_LR_score_reduction_by_sample_type")
save_table(data.frame(), "virtual_KO_LR_score_reduction_by_niche")
ko_sum <- ko[, .(control_axis_score = sum(control_axis_score, na.rm = TRUE),
                 perturbed_axis_score = sum(perturbed_axis_score, na.rm = TRUE),
                 absolute_reduction = sum(absolute_reduction, na.rm = TRUE)),
             by = .(sample_type, axis, scenario)]
ko_sum$relative_reduction <- ifelse(ko_sum$control_axis_score > 0, ko_sum$absolute_reduction / ko_sum$control_axis_score, NA_real_)
save_table(ko_sum, "virtual_KO_LR_score_reduction_summary_by_sample_type")
p_ko <- ggplot(ko_sum[ko_sum$scenario != "Control", ], aes(x = sample_type, y = relative_reduction, fill = scenario)) +
  geom_col(position = "dodge") +
  facet_wrap(~ axis) +
  labs(x = "sample_type", y = "relative score reduction", fill = "virtual KO") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_ko, "virtual_KO_score_reduction_barplot_by_sample_type", 9, 5)
ko_heat <- dcast(ko_sum[ko_sum$scenario != "Control", ], scenario + axis ~ sample_type, value.var = "relative_reduction", fill = 0)
kh <- as.matrix(ko_heat[, -(1:2), drop = FALSE])
rownames(kh) <- paste(ko_heat$scenario, ko_heat$axis, sep = " | ")
png(file.path(fig_dir, "virtual_KO_score_reduction_heatmap_by_axis.png"), width = 1400, height = 900, res = 180)
pheatmap(kh, cluster_rows = FALSE, cluster_cols = FALSE, main = "Virtual KO relative LR score reduction")
dev.off()
pdf(file.path(fig_dir, "virtual_KO_score_reduction_heatmap_by_axis.pdf"), width = 8, height = 5)
pheatmap(kh, cluster_rows = FALSE, cluster_cols = FALSE, main = "Virtual KO relative LR score reduction")
dev.off()

stat_summary <- rbind(
  data.frame(test = "integrated_oc_sample_type_LR_opportunity", status = "descriptive_only",
             note = "sample_type is evaluated within integrated_oc only; sample_id confounding must be considered."),
  data.frame(test = "external_scRNA_mixed_model", status = ifelse(length(external_files) == 0, "not_run_no_input", "pending"),
             note = "Requires external scRNA with patient_id and sample_type metadata."),
  data.frame(test = "spatial_permutation_or_neighborhood", status = ifelse(length(spatial_files) == 0, "not_run_no_input", "pending"),
             note = "Requires spatial transcriptomics or mIF coordinates."),
  data.frame(test = "bulk_OS_Cox", status = "not_primary_endpoint",
             note = "Bulk Cox instability does not reject a local niche-dependent communication mechanism.")
)
save_table(stat_summary, "statistical_tests_summary")
limitations <- data.frame(
  limitation = c("sample_type_sample_id_confounding", "expression_space_not_physical_space", "LR_scores_are_potential_interactions",
                 "Subclone_05_low_cell_count", "external_scRNA_missing", "spatial_data_missing", "bulk_OS_not_primary"),
  interpretation = c(
    "In integrated_oc, sample_type may be partially or fully confounded with sample_id; sample_type-specific findings are exploratory.",
    "PCA/UMAP/KNN proximity is expression-space proximity and should not be interpreted as true tissue adjacency.",
    "Ligand-receptor opportunity scores are expression-derived potentials, not proof of physical binding or signaling.",
    "Subclone_05 has low cell number and is used as a reference rather than a primary conclusion.",
    "No external scRNA input was provided in external_scRNA, so external validation tables are templates only.",
    "No spatial input was provided in spatial, so spatial co-localization and COMMOT results are templates only.",
    "Bulk OS interaction instability does not reject local sample_type/niche-dependent biology."
  )
)
save_table(limitations, "limitations_summary")

report <- file.path(project_root, "sample_type_LR_niche_analysis_report.md")
sink(report)
cat("# sample_type-dependent SPP1-CD44/ITGB1 LR niche analysis\n\n")
cat("Project root: `", project_root, "`\n\n", sep = "")
cat("## Background and rationale\n\n")
cat("This analysis tests whether SPP1+ myeloid cells show sample_type-dependent opportunity to communicate with CD44/ITGB1-expressing CNV-defined ovarian cancer subclones, especially CNV_Subclone_02 and CNV_Subclone_04. The analysis is exploratory inside integrated_oc because sample_type can be confounded with sample_id.\n\n")
cat("## Sample type diagnostics\n\n")
print(confound_diag)
cat("\nCNV clone composition and myeloid source-group composition tables are written to `tables/`.\n\n")
cat("## Integrated_oc LR opportunity score\n\n")
top_primary <- sample_type_lr[sample_type_lr$axis %in% c("SPP1-CD44", "SPP1-ITGB1") & sample_type_lr$target_group %in% focus_target_groups, ]
top_primary_sum <- top_primary[, .(axis_score = sum(axis_score, na.rm = TRUE),
                                   abundance_weighted_score = sum(abundance_weighted_score, na.rm = TRUE),
                                   source_n_total = sum(source_n, na.rm = TRUE),
                                   target_n_total = sum(target_n, na.rm = TRUE)),
                               by = .(sample_type, axis, target_group)]
print(top_primary_sum[order(sample_type, axis, target_group)])
cat("\n## Signature outputs\n\n")
cat("- Subclone_02-like genes: ", length(sig_list$Subclone02_like), "\n", sep = "")
cat("- Subclone_04-like genes: ", length(sig_list$Subclone04_like), "\n", sep = "")
cat("- Subclone_02/04 common genes: ", length(sig_list$Subclone02_04_common), "\n", sep = "")
cat("- CD44/ITGB1 target genes: ", length(sig_list$CD44_ITGB1_target), "\n", sep = "")
cat("- KRAS/hypoxia target genes: ", length(sig_list$KRAS_hypoxia_target), "\n\n", sep = "")
cat("## External scRNA validation\n\n")
print(external_status)
cat("\n## Spatial validation\n\n")
print(spatial_status)
cat("\n## Virtual KO\n\n")
print(ko_sum[order(sample_type, axis, scenario)])
cat("\n## Interpretation\n\n")
cat("SPP1-CD44/ITGB1 should be interpreted as a sample_type- and niche-dependent candidate tumor-myeloid interaction program rather than a universal bulk OS predictor. In integrated_oc, LR opportunity scores combine source abundance, target abundance, ligand/receptor expression intensity, and expression fraction. Spatial and external scRNA validation remain required before claiming true physical proximity or broad cohort generalization.\n\n")
cat("## Limitations\n\n")
print(limitations)
cat("\n## Key deliverables\n\n")
cat("- `tables/sample_type_LR_opportunity_scores_primary_axes.csv`\n")
cat("- `tables/sample_type_LR_opportunity_scores_target_Subclone02_04.csv`\n")
cat("- `figures/bubble_source_to_target_by_sample_type_SPP1_CD44_ITGB1.png`\n")
cat("- `figures/heatmap_sample_type_by_axis_TargetSubclone02_04.png`\n")
cat("- `tables/signature_Subclone02_like.csv`\n")
cat("- `tables/signature_Subclone04_like.csv`\n")
cat("- `tables/signature_Subclone02_04_common.csv`\n")
cat("- `tables/external_scRNA_all_datasets_sampletype_axis_scores.csv`\n")
cat("- `figures/external_boxplot_SPP1_CD44_by_sample_type.png`\n")
cat("- `figures/external_boxplot_SPP1_ITGB1_by_sample_type.png`\n")
cat("- `tables/spatial_correlation_SPP1_myeloid_Target_axis.csv`\n")
cat("- `tables/spatial_neighborhood_enrichment.csv`\n")
cat("- `figures/spatial_COMMOT_SPP1_CD44_map.png`\n")
cat("- `figures/spatial_COMMOT_SPP1_ITGB1_map.png`\n")
cat("- `tables/virtual_KO_LR_score_reduction_by_sample_type.csv`\n")
cat("- `figures/virtual_KO_score_reduction_barplot_by_sample_type.png`\n")
sink()

html <- file.path(project_root, "sample_type_LR_niche_analysis_report.html")
docx <- file.path(project_root, "sample_type_LR_niche_analysis_report.docx")
render_status <- data.frame(output = c("html", "docx"), status = "not_attempted", note = "")
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  rmd <- file.path(project_root, "sample_type_LR_niche_analysis_report.Rmd")
  writeLines(c("---", "title: \"sample_type LR niche analysis\"", "output:", "  html_document: default", "  word_document: default", "---", "", paste(readLines(report, warn = FALSE), collapse = "\n")), rmd)
  render_status$status[render_status$output == "html"] <- tryCatch({
    rmarkdown::render(rmd, output_format = "html_document", output_file = basename(html), output_dir = project_root, quiet = TRUE)
    "completed"
  }, error = function(e) paste0("failed: ", conditionMessage(e)))
  render_status$status[render_status$output == "docx"] <- tryCatch({
    rmarkdown::render(rmd, output_format = "word_document", output_file = basename(docx), output_dir = project_root, quiet = TRUE)
    "completed"
  }, error = function(e) paste0("failed: ", conditionMessage(e)))
}
save_table(render_status, "final_report_render_status")

sink(file.path(log_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
message("Finished at ", Sys.time())
