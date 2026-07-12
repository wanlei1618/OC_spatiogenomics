suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

setClass("KeyMixin", contains = "VIRTUAL", slots = list(key = "character"))
setClass("LogMap", contains = "matrix")
setClass(
  "StdAssay",
  contains = c("VIRTUAL", "KeyMixin"),
  slots = c(
    layers = "list",
    cells = "LogMap",
    features = "LogMap",
    default = "integer",
    assay.orig = "character",
    meta.data = "data.frame",
    misc = "list"
  )
)
setClass("Assay5", contains = "StdAssay")
setClass(
  "Assay",
  contains = "KeyMixin",
  slots = c(
    counts = "AnyMatrix",
    data = "AnyMatrix",
    scale.data = "matrix",
    assay.orig = "OptionalCharacter",
    var.features = "vector",
    meta.features = "data.frame",
    misc = "OptionalList"
  )
)
setClass(
  "Seurat",
  slots = c(
    assays = "list",
    meta.data = "data.frame",
    active.assay = "character",
    active.ident = "factor",
    graphs = "list",
    neighbors = "list",
    reductions = "list",
    images = "list",
    project.name = "character",
    misc = "list",
    version = "package_version",
    commands = "list",
    tools = "list"
  )
)

out_dir <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

integrated_file <- "D:/OC_spatiogenomics/infercnv/integrated_oc.RData"
t_file <- "D:/OC_spatiogenomics/infercnv/integratedocTcells.RData"
mye_file <- "D:/OC_spatiogenomics/infercnv/integratedocMyecells.RData"
b_file_candidates <- c(
  "D:/OC_spatiogenomics/infercnv/integratedocBcells.RData",
  "D:/OC_spatiogenomics/integratedocBcells.RData"
)
b_file <- b_file_candidates[file.exists(b_file_candidates)][1]
subclone_file <- file.path(out_dir, "infercnv_cell_to_subclone_k5.csv")
lr_file <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/work/omnipath_ligrecextra.tsv"

canonical_cell <- function(x) {
  sub("^([^_]+)_\\1_", "\\1_", x)
}

message("Loading integrated_oc...")
integrated_oc <- readRDS(integrated_file)
meta <- integrated_oc@meta.data
meta$cell_integrated_oc <- rownames(meta)
meta$cell_canonical <- canonical_cell(rownames(meta))

message("Mapping inferCNV subclones...")
subclone <- read.csv(subclone_file, stringsAsFactors = FALSE)
subclone$cell_canonical <- canonical_cell(subclone$cell)
subclone$cell_style2_canonical <- canonical_cell(subclone$cell_for_integratedocTcells_style)
meta$cnv_subclone <- NA_character_
idx <- match(meta$cell_canonical, subclone$cell_canonical)
meta$cnv_subclone[!is.na(idx)] <- subclone$cnv_subclone[idx[!is.na(idx)]]
idx2 <- is.na(meta$cnv_subclone)
idx_alt <- match(meta$cell_canonical[idx2], subclone$cell_style2_canonical)
meta$cnv_subclone[idx2][!is.na(idx_alt)] <- subclone$cnv_subclone[idx_alt[!is.na(idx_alt)]]

message("Mapping T/NK subtypes...")
t_obj <- readRDS(t_file)
t_meta <- t_obj@meta.data
t_meta$cell_canonical <- canonical_cell(rownames(t_meta))
meta$t_nk_subtype <- NA_character_
t_idx <- match(meta$cell_canonical, t_meta$cell_canonical)
meta$t_nk_subtype[!is.na(t_idx)] <- as.character(t_meta$cell_type[t_idx[!is.na(t_idx)]])

message("Mapping myeloid subtypes...")
mye_obj <- readRDS(mye_file)
mye_meta <- mye_obj@meta.data
mye_meta$cell_canonical <- canonical_cell(rownames(mye_meta))
meta$myeloid_subtype <- NA_character_
m_idx <- match(meta$cell_canonical, mye_meta$cell_canonical)
meta$myeloid_subtype[!is.na(m_idx)] <- as.character(mye_meta$cell_type[m_idx[!is.na(m_idx)]])

meta$b_subtype <- NA_character_
b_source <- "not found"
if (!is.na(b_file)) {
  message("Mapping B-cell subtypes from: ", b_file)
  b_obj <- readRDS(b_file)
  b_meta <- b_obj@meta.data
  b_meta$cell_canonical <- canonical_cell(rownames(b_meta))
  b_idx <- match(meta$cell_canonical, b_meta$cell_canonical)
  meta$b_subtype[!is.na(b_idx)] <- as.character(b_meta$cell_type[b_idx[!is.na(b_idx)]])
  b_source <- b_file
  write.csv(
    data.frame(
      cell_integrated_oc = rownames(meta)[!is.na(b_idx)],
      b_subtype = meta$b_subtype[!is.na(b_idx)]
    ),
    file.path(out_dir, "integrated_oc_Bcell_subtypes_mapped_from_integratedocBcells.csv"),
    row.names = FALSE
  )
} else {
  message("integratedocBcells.RData not found; falling back to B marker signatures.")
}
b_cells <- meta$cell_type == "B cells" & is.na(meta$b_subtype)
if (is.na(b_file) && any(b_cells)) {
  expr_for_b <- integrated_oc@assays$RNA@data
  b_cell_names <- rownames(meta)[b_cells]
  b_sigs <- list(
    PC_IGHG = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "MZB1", "JCHAIN", "XBP1", "TNFRSF17", "PRDM1"),
    `Bm_stress-response` = c("HSPA1B", "HSPA1A", "HSP90AA1", "DNAJB1", "HSPE1", "NR4A1", "FOSB", "DUSP2", "CD83"),
    `Classical-Bm_TXNIP` = c("TXNIP", "LAPTM5", "BTG1", "YBX3", "RCSD1", "SMAP2", "AFF3", "LINC00926"),
    `Early-PC_MS4A1low` = c("MZB1", "JCHAIN", "XBP1", "TXNDC5", "PPIB", "SEC11C", "IGHM", "IGHA1", "IGHA2"),
    Bn_TCL1A = c("TCL1A", "FCER2", "MS4A1", "IL4R", "SELL", "CD37", "BACH2", "NIBAN3")
  )
  b_scores <- matrix(
    NA_real_,
    nrow = length(b_cell_names),
    ncol = length(b_sigs),
    dimnames = list(b_cell_names, names(b_sigs))
  )
  for (nm in names(b_sigs)) {
    genes <- intersect(b_sigs[[nm]], rownames(expr_for_b))
    if (length(genes) == 0) next
    b_scores[, nm] <- Matrix::colMeans(expr_for_b[genes, b_cell_names, drop = FALSE])
  }
  b_scores[is.na(b_scores)] <- -Inf
  raw_b_label <- colnames(b_scores)[max.col(b_scores, ties.method = "first")]

  ms4a1 <- if ("MS4A1" %in% rownames(expr_for_b)) as.numeric(expr_for_b["MS4A1", b_cell_names]) else rep(0, length(b_cell_names))
  pc_core <- intersect(c("MZB1", "JCHAIN", "XBP1", "TXNDC5"), rownames(expr_for_b))
  pc_core_score <- if (length(pc_core) > 0) Matrix::colMeans(expr_for_b[pc_core, b_cell_names, drop = FALSE]) else rep(0, length(b_cell_names))
  ighg_core <- intersect(c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), rownames(expr_for_b))
  ighg_score <- if (length(ighg_core) > 0) Matrix::colMeans(expr_for_b[ighg_core, b_cell_names, drop = FALSE]) else rep(0, length(b_cell_names))

  ## Keep a separate early-PC state for plasma-program cells with relatively low IGHG.
  early_pc <- pc_core_score >= stats::quantile(pc_core_score, 0.75, na.rm = TRUE) &
    ighg_score < stats::quantile(ighg_score, 0.75, na.rm = TRUE) &
    ms4a1 < stats::quantile(ms4a1, 0.60, na.rm = TRUE)
  raw_b_label[early_pc] <- "Early-PC_MS4A1low"

  meta$b_subtype[b_cells] <- raw_b_label
  write.csv(
    data.frame(cell_integrated_oc = b_cell_names, b_subtype = raw_b_label, b_scores),
    file.path(out_dir, "integrated_oc_Bcell_subtypes_inferred_from_markers.csv"),
    row.names = FALSE
  )
}

meta$interaction_group <- as.character(meta$cell_type)
meta$interaction_group[!is.na(meta$cnv_subclone)] <- paste0("CNV_", meta$cnv_subclone[!is.na(meta$cnv_subclone)])
meta$interaction_group[is.na(meta$cnv_subclone) & !is.na(meta$t_nk_subtype)] <- paste0("T_NK_", meta$t_nk_subtype[is.na(meta$cnv_subclone) & !is.na(meta$t_nk_subtype)])
meta$interaction_group[is.na(meta$cnv_subclone) & is.na(meta$t_nk_subtype) & !is.na(meta$myeloid_subtype)] <- paste0("Myeloid_", meta$myeloid_subtype[is.na(meta$cnv_subclone) & is.na(meta$t_nk_subtype) & !is.na(meta$myeloid_subtype)])
meta$interaction_group[is.na(meta$cnv_subclone) & is.na(meta$t_nk_subtype) & is.na(meta$myeloid_subtype) & !is.na(meta$b_subtype)] <- paste0("B_", meta$b_subtype[is.na(meta$cnv_subclone) & is.na(meta$t_nk_subtype) & is.na(meta$myeloid_subtype) & !is.na(meta$b_subtype)])
meta$interaction_group[is.na(meta$interaction_group) | meta$interaction_group == ""] <- "Unassigned"

group_counts <- sort(table(meta$interaction_group), decreasing = TRUE)
keep_groups <- names(group_counts[group_counts >= 20])
keep_cells <- meta$interaction_group %in% keep_groups

write.csv(meta, file.path(out_dir, "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv"), row.names = FALSE)
write.csv(data.frame(interaction_group = names(group_counts), n_cells = as.integer(group_counts)),
          file.path(out_dir, "integrated_oc_interaction_group_cell_counts.csv"), row.names = FALSE)

mapping_summary <- data.frame(
  item = c("integrated_oc cells", "cnv subclone mapped", "T/NK subtype mapped", "myeloid subtype mapped", "B subtype source", "B subtype mapped", "groups >=20 cells"),
  value = c(nrow(meta), sum(!is.na(meta$cnv_subclone)), sum(!is.na(meta$t_nk_subtype)),
            sum(!is.na(meta$myeloid_subtype)), b_source, sum(!is.na(meta$b_subtype)), length(keep_groups))
)
write.csv(mapping_summary, file.path(out_dir, "integrated_oc_mapping_summary.csv"), row.names = FALSE)

message("Preparing expression matrix...")
expr <- integrated_oc@assays$RNA@data
if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("RNA@data lacks dimnames; gene or cell names are unavailable.")
}
cell_order <- match(meta$cell_integrated_oc[keep_cells], colnames(expr))
ok <- !is.na(cell_order)
expr <- expr[, cell_order[ok], drop = FALSE]
groups <- meta$interaction_group[keep_cells][ok]

message("Loading ligand-receptor database...")
lr <- read.delim(lr_file, stringsAsFactors = FALSE, quote = "", comment.char = "")
lr <- unique(lr[, c("source_genesymbol", "target_genesymbol")])
colnames(lr) <- c("ligand", "receptor")
lr <- lr[lr$ligand != "" & lr$receptor != "" & !is.na(lr$ligand) & !is.na(lr$receptor), ]
lr <- lr[lr$ligand %in% rownames(expr) & lr$receptor %in% rownames(expr), ]
lr <- unique(lr)
write.csv(lr, file.path(out_dir, "omnipath_ligrecextra_LR_pairs_used.csv"), row.names = FALSE)

genes <- sort(unique(c(lr$ligand, lr$receptor)))
expr_lr <- expr[genes, , drop = FALSE]
group_levels <- names(sort(table(groups), decreasing = TRUE))

message("Computing group average expression and detection fraction...")
avg_expr <- matrix(0, nrow = length(genes), ncol = length(group_levels), dimnames = list(genes, group_levels))
pct_expr <- matrix(0, nrow = length(genes), ncol = length(group_levels), dimnames = list(genes, group_levels))
for (g in group_levels) {
  cols <- which(groups == g)
  sub <- expr_lr[, cols, drop = FALSE]
  avg_expr[, g] <- Matrix::rowMeans(sub)
  pct_expr[, g] <- Matrix::rowMeans(sub > 0)
}
write.csv(avg_expr, file.path(out_dir, "lr_gene_average_expression_by_group.csv"))
write.csv(pct_expr, file.path(out_dir, "lr_gene_detection_fraction_by_group.csv"))

message("Scoring ligand-receptor interactions...")
lig_idx <- match(lr$ligand, rownames(avg_expr))
rec_idx <- match(lr$receptor, rownames(avg_expr))

summary_list <- vector("list", length(group_levels) * length(group_levels))
top_list <- vector("list", length(group_levels) * length(group_levels))
n <- 0L
min_pct <- 0.10
min_avg <- 0.05
top_n_per_pair <- 50

for (src in group_levels) {
  lig_avg <- avg_expr[lig_idx, src]
  lig_pct <- pct_expr[lig_idx, src]
  for (tgt in group_levels) {
    rec_avg <- avg_expr[rec_idx, tgt]
    rec_pct <- pct_expr[rec_idx, tgt]
    pass <- lig_avg >= min_avg & rec_avg >= min_avg & lig_pct >= min_pct & rec_pct >= min_pct
    score <- lig_avg * rec_avg * lig_pct * rec_pct
    score[!pass] <- 0
    positive <- which(score > 0)
    n <- n + 1L
    if (length(positive) > 0) {
      ord_all <- positive[order(score[positive], decreasing = TRUE)]
      best <- ord_all[1]
      summary_list[[n]] <- data.frame(
        source_group = src,
        target_group = tgt,
        n_lr_pairs = length(positive),
        total_score = sum(score[positive]),
        mean_score = mean(score[positive]),
        max_score = score[best],
        top_ligand = lr$ligand[best],
        top_receptor = lr$receptor[best],
        stringsAsFactors = FALSE
      )
      keep <- head(ord_all, top_n_per_pair)
      top_list[[n]] <- data.frame(
        source_group = src,
        target_group = tgt,
        ligand = lr$ligand[keep],
        receptor = lr$receptor[keep],
        score = score[keep],
        ligand_avg = lig_avg[keep],
        receptor_avg = rec_avg[keep],
        ligand_pct = lig_pct[keep],
        receptor_pct = rec_pct[keep],
        stringsAsFactors = FALSE
      )
    } else {
      summary_list[[n]] <- data.frame(
        source_group = src, target_group = tgt, n_lr_pairs = 0,
        total_score = 0, mean_score = 0, max_score = 0,
        top_ligand = NA_character_, top_receptor = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
}

pair_summary <- do.call(rbind, summary_list)
top_interactions <- do.call(rbind, top_list[!vapply(top_list, is.null, logical(1))])
top_interactions <- top_interactions[order(top_interactions$score, decreasing = TRUE), ]

write.csv(pair_summary[order(pair_summary$total_score, decreasing = TRUE), ],
          file.path(out_dir, "lr_group_pair_summary_connectome_like.csv"), row.names = FALSE)
write.csv(top_interactions,
          file.path(out_dir, "lr_top_interactions_top50_per_group_pair.csv"), row.names = FALSE)
write.csv(head(top_interactions, 5000),
          file.path(out_dir, "lr_top5000_interactions_global.csv"), row.names = FALSE)

subclone_involved <- grepl("^CNV_Subclone", top_interactions$source_group) | grepl("^CNV_Subclone", top_interactions$target_group)
write.csv(top_interactions[subclone_involved, ],
          file.path(out_dir, "lr_top_interactions_involving_cnv_subclones.csv"), row.names = FALSE)
write.csv(top_interactions[grepl("^CNV_Subclone", top_interactions$source_group) & !grepl("^CNV_Subclone", top_interactions$target_group), ],
          file.path(out_dir, "lr_top_interactions_cnv_subclone_to_celltype.csv"), row.names = FALSE)
write.csv(top_interactions[!grepl("^CNV_Subclone", top_interactions$source_group) & grepl("^CNV_Subclone", top_interactions$target_group), ],
          file.path(out_dir, "lr_top_interactions_celltype_to_cnv_subclone.csv"), row.names = FALSE)

score_mat <- matrix(0, nrow = length(group_levels), ncol = length(group_levels),
                    dimnames = list(group_levels, group_levels))
score_mat[cbind(pair_summary$source_group, pair_summary$target_group)] <- pair_summary$total_score
log_score_mat <- log10(score_mat + 1)

group_type <- ifelse(grepl("^CNV_Subclone", group_levels), "CNV subclone",
                     ifelse(grepl("^T_NK_", group_levels), "T/NK subtype",
                            ifelse(grepl("^Myeloid_", group_levels), "Myeloid subtype",
                                   ifelse(grepl("^B_", group_levels), "B subtype", "Other broad type"))))
type_cols <- c(
  "CNV subclone" = "#4C78A8",
  "T/NK subtype" = "#59A14F",
  "Myeloid subtype" = "#F28E2B",
  "B subtype" = "#B07AA1",
  "Other broad type" = "#BAB0AC"
)
ha_col <- HeatmapAnnotation(Target = group_type, col = list(Target = type_cols), simple_anno_size = unit(3, "mm"))
ha_row <- rowAnnotation(Source = group_type, col = list(Source = type_cols), simple_anno_size = unit(3, "mm"))

ht <- Heatmap(
  log_score_mat,
  name = "log10(score+1)",
  col = colorRamp2(c(min(log_score_mat), median(log_score_mat), max(log_score_mat)), c("white", "#F4A582", "#B2182B")),
  top_annotation = ha_col,
  left_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title = "Target group",
  row_title = "Source group",
  use_raster = TRUE
)
pdf(file.path(out_dir, "lr_group_pair_total_score_heatmap.pdf"), width = 12, height = 11)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
png(file.path(out_dir, "lr_group_pair_total_score_heatmap.png"), width = 2400, height = 2200, res = 200)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cnv_groups <- group_levels[grepl("^CNV_Subclone", group_levels)]
non_cnv_groups <- setdiff(group_levels, cnv_groups)
if (length(cnv_groups) > 0 && length(non_cnv_groups) > 0) {
  cnv_to_cell <- log10(score_mat[cnv_groups, non_cnv_groups, drop = FALSE] + 1)
  cell_to_cnv <- log10(score_mat[non_cnv_groups, cnv_groups, drop = FALSE] + 1)
  pdf(file.path(out_dir, "lr_cnv_subclone_celltype_directional_heatmaps.pdf"), width = 13, height = 8)
  ht1 <- Heatmap(cnv_to_cell, name = "CNV->cell", col = colorRamp2(c(min(cnv_to_cell), median(cnv_to_cell), max(cnv_to_cell)), c("white", "#F4A582", "#B2182B")),
                 row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 6), column_title = "Subclone ligands -> cell-type receptors")
  ht2 <- Heatmap(t(cell_to_cnv), name = "cell->CNV", col = colorRamp2(c(min(cell_to_cnv), median(cell_to_cnv), max(cell_to_cnv)), c("white", "#D7B5D8", "#762A83")),
                 row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 6), column_title = "Cell-type ligands -> subclone receptors")
  draw(ht1 %v% ht2, heatmap_legend_side = "right")
  dev.off()
}

top_subclone <- top_interactions[subclone_involved, ]
top_subclone <- head(top_subclone[order(top_subclone$score, decreasing = TRUE), ], 40)
if (nrow(top_subclone) > 0) {
  top_subclone$pair <- paste(top_subclone$ligand, top_subclone$receptor, sep = " -> ")
  top_subclone$direction <- paste(top_subclone$source_group, top_subclone$target_group, sep = " => ")
  top_subclone$label <- paste(top_subclone$direction, top_subclone$pair, sep = "\n")
  p <- ggplot(top_subclone, aes(x = reorder(label, score), y = score, fill = grepl("^CNV_Subclone", source_group))) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#4C78A8", "FALSE" = "#F28E2B"), guide = "none") +
    labs(x = NULL, y = "LR score", title = "Top ligand-receptor interactions involving CNV subclones") +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 5))
  ggsave(file.path(out_dir, "lr_top40_subclone_involving_pairs.png"), p, width = 12, height = 10, dpi = 220)
  ggsave(file.path(out_dir, "lr_top40_subclone_involving_pairs.pdf"), p, width = 12, height = 10)
}

message("Done.")
message("Groups: ", length(group_levels), "; LR pairs used: ", nrow(lr), "; top interactions: ", nrow(top_interactions))
