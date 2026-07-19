suppressPackageStartupMessages({
  library(Matrix)
  library(limma)
  library(ggplot2)
  library(fgsea)
})

base_dir <- "D:/OC_spatiogenomics/infercnv"
out_dir <- file.path(base_dir, "SPP1_ITGB1_CD44_hypothesis_validation")
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
scripts_dir <- file.path(out_dir, "scripts")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading integrated_oc object and previous annotations...")
obj <- readRDS(file.path(base_dir, "integrated_oc.RData"))
expr <- obj@assays[["RNA"]]@data
meta <- read.csv(file.path(base_dir, "integrated_oc_plan_analysis/tables/integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("cell_integrated_oc" %in% colnames(meta))

common_cells <- intersect(colnames(expr), meta$cell_integrated_oc)
expr <- expr[, common_cells, drop = FALSE]
meta <- meta[match(common_cells, meta$cell_integrated_oc), , drop = FALSE]
rownames(meta) <- meta$cell_integrated_oc

cell_scores <- read.csv(file.path(base_dir, "integrated_oc_plan_analysis/tables/cell_level_function_module_scores.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
cell_scores <- cell_scores[match(common_cells, cell_scores$cell_integrated_oc), , drop = FALSE]

tumor_idx <- which(!is.na(meta$cnv_subclone) & meta$cnv_subclone != "")
tumor_meta <- meta[tumor_idx, , drop = FALSE]
tumor_expr <- expr[, tumor_idx, drop = FALSE]
tumor_scores <- cell_scores[tumor_idx, , drop = FALSE]
tumor_meta$focus_clone <- ifelse(tumor_meta$cnv_subclone %in% c("Subclone_02", "Subclone_04"),
                                 "Subclone_02_04", "Other_subclones")

present <- function(genes) intersect(genes, rownames(expr))
epithelial_genes <- c("PAX8", "MUC16", "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4")
immune_genes <- c("PTPRC", "LYZ", "LST1", "TYROBP", "CD68", "C1QA", "C1QB", "C1QC")
receptor_genes <- c("ITGB1", "CD44")
sender_ligands <- c("SPP1", "APOE", "MIF", "TGFB1", "VEGFA", "CSF1", "CD47")
all_marker_genes <- unique(c(epithelial_genes, immune_genes, receptor_genes, sender_ligands))
marker_present <- present(all_marker_genes)
missing_markers <- setdiff(all_marker_genes, marker_present)
write.csv(data.frame(gene = all_marker_genes, present = all_marker_genes %in% marker_present),
          file.path(tables_dir, "marker_gene_presence_in_integrated_oc.csv"), row.names = FALSE)

mean_sparse <- function(x) {
  if (length(x) == 0) return(rep(NA_real_, nrow(expr)))
  Matrix::rowMeans(expr[, x, drop = FALSE])
}

gene_group_stats <- function(mat, md, genes, group_col) {
  genes <- intersect(genes, rownames(mat))
  groups <- sort(unique(md[[group_col]]))
  out <- list()
  k <- 1
  for (g in groups) {
    cells <- rownames(md)[md[[group_col]] == g]
    sub <- mat[genes, cells, drop = FALSE]
    for (gene in genes) {
      vals <- as.numeric(sub[gene, ])
      out[[k]] <- data.frame(
        group = g,
        gene = gene,
        mean_expr = mean(vals),
        median_expr = median(vals),
        pct_expr = mean(vals > 0),
        n_cells = length(vals),
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  do.call(rbind, out)
}

message("5.1 Marker identity and receptor stability...")
clone_gene_stats <- gene_group_stats(tumor_expr, tumor_meta, all_marker_genes, "cnv_subclone")
clone_gene_stats$gene_class <- ifelse(clone_gene_stats$gene %in% epithelial_genes, "Tumor epithelial",
                               ifelse(clone_gene_stats$gene %in% immune_genes, "Immune contamination",
                               ifelse(clone_gene_stats$gene %in% receptor_genes, "Receptor", "Ligand/control")))
write.csv(clone_gene_stats, file.path(tables_dir, "clone_marker_receptor_expression_summary.csv"), row.names = FALSE)

plot_stats <- clone_gene_stats[clone_gene_stats$gene %in% c(epithelial_genes, immune_genes, receptor_genes), ]
p1 <- ggplot(plot_stats, aes(x = group, y = gene)) +
  geom_point(aes(size = pct_expr, color = mean_expr)) +
  scale_color_gradient(low = "grey90", high = "#B2182B") +
  scale_size(range = c(0.5, 7)) +
  facet_grid(gene_class ~ ., scales = "free_y", space = "free_y") +
  labs(x = "CNV subclone", y = NULL, color = "Mean expression", size = "Pct > 0") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95"))
ggsave(file.path(figures_dir, "clone_tumor_identity_receptor_dotplot.png"), p1, width = 8.5, height = 7, dpi = 300)
ggsave(file.path(figures_dir, "clone_tumor_identity_receptor_dotplot.pdf"), p1, width = 8.5, height = 7)

score_genes <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[genes, , drop = FALSE])
}
tumor_cell_validation <- data.frame(
  cell_integrated_oc = colnames(tumor_expr),
  cnv_subclone = tumor_meta$cnv_subclone,
  batch = tumor_meta$batch,
  sample_type = tumor_meta$sample_type,
  focus_clone = tumor_meta$focus_clone,
  epithelial_score = score_genes(tumor_expr, epithelial_genes),
  immune_contamination_score = score_genes(tumor_expr, immune_genes),
  ITGB1 = if ("ITGB1" %in% rownames(tumor_expr)) as.numeric(tumor_expr["ITGB1", ]) else NA_real_,
  CD44 = if ("CD44" %in% rownames(tumor_expr)) as.numeric(tumor_expr["CD44", ]) else NA_real_,
  KRAS_activation = tumor_scores$KRAS_activation,
  Hypoxia = tumor_scores$Hypoxia,
  Stemness_epithelial = tumor_scores$Stemness_epithelial,
  stringsAsFactors = FALSE
)
write.csv(tumor_cell_validation, file.path(tables_dir, "tumor_cell_marker_receptor_scores.csv"), row.names = FALSE)

clone_validation <- aggregate(cbind(epithelial_score, immune_contamination_score, ITGB1, CD44,
                                    KRAS_activation, Hypoxia, Stemness_epithelial) ~ cnv_subclone,
                              data = tumor_cell_validation, FUN = function(x) mean(x, na.rm = TRUE))
clone_validation$n_cells <- as.integer(table(tumor_cell_validation$cnv_subclone)[clone_validation$cnv_subclone])
write.csv(clone_validation, file.path(tables_dir, "clone_malignant_identity_receptor_stability_summary.csv"), row.names = FALSE)

sample_clone_summary <- aggregate(cbind(ITGB1, CD44, KRAS_activation, Hypoxia, Stemness_epithelial,
                                        epithelial_score, immune_contamination_score) ~ batch + sample_type + cnv_subclone + focus_clone,
                                  data = tumor_cell_validation,
                                  FUN = function(x) mean(x, na.rm = TRUE))
sample_clone_summary$n_cells <- as.integer(aggregate(cell_integrated_oc ~ batch + sample_type + cnv_subclone + focus_clone,
                                                     data = tumor_cell_validation, FUN = length)$cell_integrated_oc)
write.csv(sample_clone_summary, file.path(tables_dir, "sample_batch_stratified_receptor_program_scores.csv"), row.names = FALSE)

long_sample <- reshape(sample_clone_summary,
                       varying = c("ITGB1", "CD44", "KRAS_activation", "Hypoxia", "Stemness_epithelial"),
                       v.names = "mean_score", timevar = "feature",
                       times = c("ITGB1", "CD44", "KRAS_activation", "Hypoxia", "Stemness_epithelial"),
                       direction = "long")
long_sample$feature <- factor(long_sample$feature, levels = c("ITGB1", "CD44", "KRAS_activation", "Hypoxia", "Stemness_epithelial"))
p2 <- ggplot(long_sample, aes(x = cnv_subclone, y = mean_score, color = batch, group = batch)) +
  geom_point(aes(size = n_cells), alpha = 0.85, position = position_jitter(width = 0.08, height = 0)) +
  geom_line(alpha = 0.35) +
  facet_wrap(~feature, scales = "free_y", ncol = 3) +
  labs(x = "CNV subclone", y = "Sample x clone mean score", color = "Batch/sample", size = "Cells") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figures_dir, "sample_stratified_receptor_program_scores.png"), p2, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(figures_dir, "sample_stratified_receptor_program_scores.pdf"), p2, width = 10, height = 6.5)

feature_tests <- lapply(c("ITGB1", "CD44", "KRAS_activation", "Hypoxia", "Stemness_epithelial",
                          "epithelial_score", "immune_contamination_score"), function(feat) {
  d <- sample_clone_summary[!is.na(sample_clone_summary[[feat]]), ]
  wt <- tryCatch(wilcox.test(d[[feat]] ~ d$focus_clone), error = function(e) NULL)
  data.frame(
    feature = feat,
    mean_focus_02_04 = mean(d[[feat]][d$focus_clone == "Subclone_02_04"], na.rm = TRUE),
    mean_other = mean(d[[feat]][d$focus_clone == "Other_subclones"], na.rm = TRUE),
    delta_focus_minus_other = mean(d[[feat]][d$focus_clone == "Subclone_02_04"], na.rm = TRUE) -
      mean(d[[feat]][d$focus_clone == "Other_subclones"], na.rm = TRUE),
    wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value,
    stringsAsFactors = FALSE
  )
})
feature_tests <- do.call(rbind, feature_tests)
feature_tests$wilcox_fdr <- p.adjust(feature_tests$wilcox_p, method = "BH")
write.csv(feature_tests, file.path(tables_dir, "focus_clone_02_04_vs_others_sample_stratified_tests.csv"), row.names = FALSE)

message("5.2 Pseudo-bulk expression program validation...")
pb_groups <- paste(tumor_meta$batch, tumor_meta$cnv_subclone, sep = "|")
pb_fac <- factor(pb_groups)
pb_mm <- sparseMatrix(i = seq_along(pb_fac), j = as.integer(pb_fac), x = 1,
                      dims = c(length(pb_fac), nlevels(pb_fac)),
                      dimnames = list(NULL, levels(pb_fac)))
n_by_group <- as.numeric(Matrix::colSums(pb_mm))
pb_avg <- tumor_expr %*% pb_mm
pb_avg <- t(t(pb_avg) / pmax(n_by_group, 1))
pb_meta <- do.call(rbind, strsplit(colnames(pb_avg), "\\|", fixed = FALSE))
pb_meta <- data.frame(pb_group = colnames(pb_avg), batch = pb_meta[, 1], cnv_subclone = pb_meta[, 2],
                      n_cells = n_by_group, stringsAsFactors = FALSE)
pb_meta$focus_clone <- ifelse(pb_meta$cnv_subclone %in% c("Subclone_02", "Subclone_04"),
                              "Subclone_02_04", "Other_subclones")
write.csv(pb_meta, file.path(tables_dir, "pseudo_bulk_sample_clone_metadata_focus_0204.csv"), row.names = FALSE)

keep <- pb_meta$n_cells >= 10
pb_use <- as.matrix(pb_avg[, keep, drop = FALSE])
pbm <- pb_meta[keep, , drop = FALSE]
design <- model.matrix(~ focus_clone + batch, data = pbm)
fit <- lmFit(pb_use, design)
fit <- eBayes(fit)
coef_name <- grep("^focus_cloneSubclone_02_04$", colnames(design), value = TRUE)
if (length(coef_name) == 0) {
  coef_name <- grep("focus_clone", colnames(design), value = TRUE)[1]
}
deg <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
deg$gene <- rownames(deg)
deg <- deg[, c("gene", setdiff(colnames(deg), "gene"))]
write.csv(deg, file.path(tables_dir, "pseudo_bulk_limma_focus_Subclone02_04_vs_others.csv"), row.names = FALSE)

pathways <- list(
  HALLMARK_HYPOXIA = c("CA9","VEGFA","SLC2A1","ENO1","LDHA","PGK1","ALDOA","BNIP3","NDRG1","ADM","P4HA1","PLOD2","ANGPTL4","HILPDA","LOX","SERPINE1"),
  KRAS_SIGNALING_UP = c("KRAS","DUSP6","FOSL1","JUN","FOS","EGR1","EREG","AREG","MYC","SPRY2","ATF3","IER3","PLAUR","CXCL8","IL6"),
  FAK_INTEGRIN_AKT_ERK = c("PTK2","SRC","PXN","VCL","TLN1","ITGB1","ITGA5","ITGAV","CD44","FN1","COL1A1","COL1A2","MAPK1","MAPK3","AKT1","PIK3CA"),
  ECM_REMODELING_MIGRATION = c("SPP1","FN1","COL1A1","COL1A2","COL3A1","COL5A1","MMP2","MMP9","MMP14","THBS1","VCAN","TNC","LOX","SERPINE1","VIM"),
  INTEGRIN_CELL_ADHESION = c("ITGB1","ITGB5","ITGA2","ITGA3","ITGA5","ITGAV","CD44","ICAM1","VCAM1","LAMC1","LAMB1","FN1","PXN","VCL","TLN1"),
  STEMNESS_EPITHELIAL = c("EPCAM","KRT8","KRT18","KRT19","CLDN3","CLDN4","SOX9","PROM1","ALDH1A1","MUC16","PAX8","MSLN","TACSTD2")
)
pathways <- lapply(pathways, function(g) intersect(g, deg$gene))
ranks <- deg$t
names(ranks) <- deg$gene
ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
fg <- fgsea(pathways = pathways, stats = ranks, minSize = 5, maxSize = 500, nperm = 10000)
fg <- fg[order(fg$padj, -abs(fg$NES)), ]
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
write.csv(fg, file.path(tables_dir, "pseudo_bulk_focus_0204_GSEA_curated_pathways.csv"), row.names = FALSE)

top_deg <- head(deg[order(deg$P.Value), ], 30)
p3 <- ggplot(deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.35, size = 0.7, color = "grey45") +
  geom_point(data = subset(deg, gene %in% unique(unlist(pathways))), color = "#2166AC", alpha = 0.8, size = 1.2) +
  geom_text(data = top_deg[1:min(12, nrow(top_deg)), ], aes(label = gene), size = 2.5, vjust = -0.5) +
  labs(x = "Focus Subclone_02/04 vs others logFC", y = "-log10(P value)") +
  theme_bw(base_size = 10)
ggsave(file.path(figures_dir, "pseudo_bulk_focus_0204_vs_others_volcano.png"), p3, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "pseudo_bulk_focus_0204_vs_others_volcano.pdf"), p3, width = 7.5, height = 6)

p4 <- ggplot(fg, aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.1)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey65")) +
  labs(x = NULL, y = "NES", fill = "FDR < 0.1") +
  theme_bw(base_size = 10)
ggsave(file.path(figures_dir, "pseudo_bulk_focus_0204_GSEA_curated_pathways.png"), p4, width = 7.5, height = 4.5, dpi = 300)
ggsave(file.path(figures_dir, "pseudo_bulk_focus_0204_GSEA_curated_pathways.pdf"), p4, width = 7.5, height = 4.5)

message("5.3 Focused LR replication and target-program proxy...")
focused_lr <- read.csv(file.path(base_dir, "integrated_oc_plan_analysis/tables/focused_LR_axes_involving_cnv_subclones.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
senders <- c("Myeloid_Macro-Inflammatory_TNF", "Myeloid_Interferon-Responsive Myeloid", "Myeloid_Macro-C3/CX3CR1")
receivers <- c("CNV_Subclone_02", "CNV_Subclone_04")
axes <- data.frame(
  ligand = c("SPP1", "SPP1", "APOE", "MIF", "TGFB1", "VEGFA", "VEGFA", "CSF1", "CD47"),
  receptor = c("ITGB1", "CD44", "LRP1", "CD74", "TGFBR1", "KDR", "FLT1", "CSF1R", "SIRPA"),
  stringsAsFactors = FALSE
)
axis_lr <- focused_lr[focused_lr$source_group %in% senders &
                        focused_lr$target_group %in% receivers &
                        paste(focused_lr$ligand, focused_lr$receptor) %in% paste(axes$ligand, axes$receptor), ]
axis_lr <- axis_lr[order(axis_lr$ligand, axis_lr$receptor, -axis_lr$score), ]
write.csv(axis_lr, file.path(tables_dir, "focused_LR_predefined_axes_senders_to_Subclone02_04.csv"), row.names = FALSE)

lr_rank <- focused_lr[focused_lr$source_group %in% senders & focused_lr$target_group %in% receivers, ]
lr_rank <- lr_rank[order(lr_rank$source_group, lr_rank$target_group, -lr_rank$score), ]
lr_rank$rank_in_source_target <- ave(-lr_rank$score, lr_rank$source_group, lr_rank$target_group,
                                     FUN = function(x) rank(x, ties.method = "min"))
write.csv(lr_rank, file.path(tables_dir, "focused_LR_ranked_all_pairs_selected_senders_receivers.csv"), row.names = FALSE)

target_sets <- list(
  FAK_AKT_ERK = c("PTK2","SRC","PXN","VCL","TLN1","ITGB1","CD44","MAPK1","MAPK3","AKT1","PIK3CA"),
  ECM_REMODELING = c("FN1","COL1A1","COL1A2","COL3A1","MMP2","MMP9","MMP14","THBS1","VCAN","TNC","LOX","SERPINE1","VIM"),
  MIGRATION_INVASION = c("VIM","SNAI1","SNAI2","ZEB1","MMP2","MMP9","ITGB1","CD44","CXCL8","PLAUR","SERPINE1"),
  HYPOXIA_RESPONSE = pathways$HALLMARK_HYPOXIA
)
target_enrich <- lapply(names(target_sets), function(nm) {
  genes <- intersect(target_sets[[nm]], deg$gene)
  top_genes <- deg$gene[deg$adj.P.Val < 0.1 & deg$logFC > 0]
  universe <- deg$gene
  mat <- matrix(c(
    sum(universe %in% genes & universe %in% top_genes),
    sum(universe %in% genes & !(universe %in% top_genes)),
    sum(!(universe %in% genes) & universe %in% top_genes),
    sum(!(universe %in% genes) & !(universe %in% top_genes))
  ), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  data.frame(program = nm, genes_in_set = length(genes), up_DEG_overlap = mat[1, 1],
             fisher_p = ft$p.value, odds_ratio = unname(ft$estimate),
             overlap_genes = paste(intersect(genes, top_genes), collapse = ";"),
             stringsAsFactors = FALSE)
})
target_enrich <- do.call(rbind, target_enrich)
target_enrich$fdr <- p.adjust(target_enrich$fisher_p, method = "BH")
write.csv(target_enrich, file.path(tables_dir, "SPP1_axis_target_program_proxy_enrichment.csv"), row.names = FALSE)

cellchat_status <- data.frame(step = "CellChat", status = "not_run", note = "", stringsAsFactors = FALSE)
cellchat_file <- file.path(tables_dir, "CellChat_subset_selected_sender_receiver_interactions.csv")
try({
  suppressPackageStartupMessages(library(CellChat))
  data(CellChatDB.human)
  selected_groups <- unique(c(senders, receivers))
  cc_cells <- meta$cell_integrated_oc[meta$interaction_group %in% selected_groups]
  cc_groups <- meta$interaction_group[match(cc_cells, meta$cell_integrated_oc)]
  cc_data <- expr[, cc_cells, drop = FALSE]
  cc_meta <- data.frame(labels = cc_groups, row.names = cc_cells, stringsAsFactors = FALSE)
  cellchat <- createCellChat(object = cc_data, meta = cc_meta, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  df.net <- df.net[df.net$source %in% senders & df.net$target %in% receivers, ]
  write.csv(df.net, cellchat_file, row.names = FALSE)
  cellchat_status <- data.frame(step = "CellChat", status = "completed",
                                note = paste("Rows:", nrow(df.net)), stringsAsFactors = FALSE)
}, silent = TRUE)
write.csv(cellchat_status, file.path(tables_dir, "CellChat_run_status.csv"), row.names = FALSE)

p5_data <- axis_lr[axis_lr$ligand == "SPP1" & axis_lr$receptor %in% c("ITGB1", "CD44"), ]
if (nrow(p5_data) > 0) {
  p5 <- ggplot(p5_data, aes(x = target_group, y = source_group, fill = score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste(ligand, receptor, sprintf("%.2f", score), sep = "\n")), size = 2.5) +
    scale_fill_gradient(low = "grey92", high = "#B2182B") +
    labs(x = "Receiver clone", y = "Sender myeloid subtype", fill = "LR score") +
    theme_bw(base_size = 10)
  ggsave(file.path(figures_dir, "SPP1_ITGB1_CD44_focused_LR_sender_receiver_heatmap.png"), p5, width = 8.5, height = 4.8, dpi = 300)
  ggsave(file.path(figures_dir, "SPP1_ITGB1_CD44_focused_LR_sender_receiver_heatmap.pdf"), p5, width = 8.5, height = 4.8)
}

message("5.4 Bulk model template and 5.5 spatial panel...")
bulk_inputs <- data.frame(
  variable = c("OS_time", "OS_event", "SPP1_TAM_score", "ITGB1_CD44_tumor_score", "KRAS_Hypoxia_score",
               "macrophage_fraction", "tumor_purity", "stage", "grade", "residual_disease"),
  definition = c("overall survival time", "overall survival event, 1=event",
                 "bulk score for SPP1+ TAM/macrophage program",
                 "bulk tumor receptor score, mean/z-score of ITGB1 and CD44 with tumor epithelial deconvolution",
                 "combined KRAS activation and hypoxia score",
                 "estimated macrophage fraction, e.g. CIBERSORT/xCell/EPIC",
                 "tumor purity estimate, e.g. ABSOLUTE/ESTIMATE",
                 "clinical stage", "tumor grade", "residual disease status"),
  required = TRUE,
  stringsAsFactors = FALSE
)
write.csv(bulk_inputs, file.path(tables_dir, "bulk_interaction_model_required_inputs.csv"), row.names = FALSE)

bulk_script <- c(
  "library(survival)",
  "bulk_meta <- read.csv('bulk_meta_with_scores.csv')",
  "fit <- coxph(Surv(OS_time, OS_event) ~ SPP1_TAM_score * ITGB1_CD44_tumor_score +",
  "               KRAS_Hypoxia_score + macrophage_fraction + tumor_purity +",
  "               stage + grade + residual_disease, data = bulk_meta)",
  "summary(fit)"
)
writeLines(bulk_script, file.path(scripts_dir, "bulk_interaction_cox_model_template.R"))

spatial_panel <- data.frame(
  class = c("Tumor", "Tumor", "Receptor", "Receptor", "Macrophage", "Macrophage", "Macrophage",
            "Ligand", "State", "State", "State"),
  marker = c("PAX8", "EPCAM/KRT8", "ITGB1", "CD44", "CD68", "CD163", "CD14",
             "SPP1", "HIF1A/CAIX", "VEGFA", "pFAK"),
  use = c("define ovarian epithelial tumor nuclei/cells", "backup tumor epithelial gate",
          "receiver receptor on tumor clone", "receiver receptor on tumor clone",
          "pan-macrophage gate", "TAM/M2-like enrichment", "monocyte/macrophage gate",
          "ligand-producing macrophage subset", "hypoxia state near interface",
          "angiogenic/hypoxia-associated state", "downstream integrin/FAK activation"),
  stringsAsFactors = FALSE
)
write.csv(spatial_panel, file.path(tables_dir, "spatial_validation_marker_panel.csv"), row.names = FALSE)

spatial_metrics <- data.frame(
  metric = c("nearest_neighbor_distance", "interface_density", "state_intensity_near_interface", "clinical_association"),
  definition = c("distance from SPP1+CD68+ macrophage to nearest PAX8+ITGB1+/CD44+ tumor cell",
                 "density of SPP1+ macrophages within a fixed radius of receptor-positive tumor cells",
                 "pFAK/HIF1A/CAIX/VEGFA intensity in tumor cells near SPP1+ macrophages",
                 "association of high-interface regions with ascites/progression/clinical endpoint"),
  output = c("per-cell and per-region distance distribution", "region-level interaction score",
             "state activation score", "statistical model with clinical covariates"),
  stringsAsFactors = FALSE
)
write.csv(spatial_metrics, file.path(tables_dir, "spatial_validation_scoring_strategy.csv"), row.names = FALSE)

copykat_status <- data.frame(
  method = c("inferCNV_HMM_existing", "CopyKAT", "CaSpER"),
  status = c("completed_from_existing_inferCNV_subclones_and_CNV_burden",
             ifelse("copykat" %in% rownames(installed.packages()), "available_not_run", "not_installed"),
             ifelse("CaSpER" %in% rownames(installed.packages()), "available_not_run", "not_installed")),
  note = c("Subclone_02/04 are cells with inferCNV-derived CNV labels; use CNV burden and event tables for malignant validation.",
           "Install/run separately to provide orthogonal aneuploidy calls.",
           "Install/run separately to provide RNA-CNV orthogonal validation."),
  stringsAsFactors = FALSE
)
write.csv(copykat_status, file.path(tables_dir, "orthogonal_CNV_method_status_CopyKAT_CaSpER.csv"), row.names = FALSE)

report <- c(
  "# SPP1-ITGB1/CD44 hypothesis validation",
  "",
  "## Output scope",
  "This folder validates the hypothesis that myeloid/macrophage-derived SPP1 signals to CNV_Subclone_02/04 through ITGB1/CD44.",
  "",
  "## 5.1 Target clone malignant identity and receptor robustness",
  paste0("- Tumor/receptor/immune marker summary: `tables/clone_marker_receptor_expression_summary.csv`."),
  paste0("- Sample/batch-stratified receptor and state scores: `tables/sample_batch_stratified_receptor_program_scores.csv`."),
  paste0("- Focus Subclone_02/04 vs other clone tests: `tables/focus_clone_02_04_vs_others_sample_stratified_tests.csv`."),
  paste0("- Main figures: `figures/clone_tumor_identity_receptor_dotplot.png`, `figures/sample_stratified_receptor_program_scores.png`."),
  "",
  "## 5.2 Pseudo-bulk expression program validation",
  "- Pseudo-bulk groups were sample/batch x cnv_subclone, using normalized RNA expression averages because raw counts are not fully exposed in the object.",
  "- DEG table: `tables/pseudo_bulk_limma_focus_Subclone02_04_vs_others.csv`.",
  "- Curated GSEA table: `tables/pseudo_bulk_focus_0204_GSEA_curated_pathways.csv`.",
  "",
  "## 5.3 LR replication",
  "- Existing Connectome-like focused LR results directly support SPP1->ITGB1/CD44 from myeloid/macrophage senders to CNV_Subclone_02/04.",
  "- Focused predefined axis table: `tables/focused_LR_predefined_axes_senders_to_Subclone02_04.csv`.",
  "- Ranked selected sender/receiver table: `tables/focused_LR_ranked_all_pairs_selected_senders_receivers.csv`.",
  "- CellChat run status: `tables/CellChat_run_status.csv`.",
  "- NicheNet package was not installed; a target-program enrichment proxy is provided in `tables/SPP1_axis_target_program_proxy_enrichment.csv`.",
  "",
  "## 5.4 Bulk interaction model",
  "- Bulk cohorts were not present locally, so the Cox interaction model was not executed.",
  "- Required inputs: `tables/bulk_interaction_model_required_inputs.csv`.",
  "- Template script: `scripts/bulk_interaction_cox_model_template.R`.",
  "",
  "## 5.5 Spatial validation pre-analysis",
  "- Marker panel: `tables/spatial_validation_marker_panel.csv`.",
  "- Scoring strategy: `tables/spatial_validation_scoring_strategy.csv`.",
  "",
  "## Environment caveats",
  paste0("- Missing marker genes in integrated_oc RNA matrix: ", ifelse(length(missing_markers) == 0, "none", paste(missing_markers, collapse = ", ")), "."),
  "- CopyKAT/CaSpER/LIANA/NicheNet/DESeq2/edgeR were not installed in this local R 4.0.3 environment.",
  "- The Seurat package itself is not loadable here, but the Seurat object's metadata and RNA data matrix were accessible through object slots."
)
writeLines(report, file.path(out_dir, "SPP1_ITGB1_CD44_hypothesis_validation_report.md"))
writeLines(report, file.path(out_dir, "SPP1_ITGB1_CD44_hypothesis_validation_report.txt"))

message("Done: ", out_dir)
