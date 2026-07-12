suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(future)
})

set.seed(1)
future::plan("sequential")
options(stringsAsFactors = FALSE)

setClass("KeyMixin", contains = "VIRTUAL", slots = list(key = "character"))
setClass("LogMap", contains = "matrix")
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

out_dir <- "D:/spatiogenomics_new"
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
report_dir <- file.path(out_dir, "reports")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

candidate_genes <- c("CXCL12", "CXCR4", "MDK", "SDC4", "CPT1A")
signature_list <- list(
  CXCL12_CXCR4 = c("CXCL12", "CXCR4"),
  MDK_SDC4 = c("MDK", "SDC4"),
  FAO_CPT1A = c("CPT1A", "CPT1B", "CPT1C", "CPT2", "SLC25A20", "ACSL1", "ACSL3",
                "ACSL4", "ACSL5", "CD36", "ACADM", "ACADS", "ACADSB", "ACADVL",
                "HADHA", "HADHB", "ACOX1", "ETFA", "ETFB", "ETFDH", "PPARA", "PPARGC1A"),
  TLS = c("CXCL13", "CCL19", "CCL21", "CCR7", "LTA", "LTB", "MS4A1",
          "CD19", "CD79A", "BANK1", "CD3D", "CD3E", "CD4", "CD8A",
          "CD8B", "CXCR5", "ICOS", "PDCD1", "BCL6", "IL21", "CD40LG",
          "LAMP3", "POU2AF1"),
  CA_MSC = c("ACTA2", "FAP", "PDGFRA", "PDGFRB", "THY1", "PDPN", "POSTN",
             "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "VIM",
             "TAGLN", "MYL9", "CXCL12", "IL6", "TGFB1", "S100A4", "MMP2",
             "VCAN", "SPARC")
)

canonical_cell <- function(x) {
  sub("^([^_]+)_\\1_", "\\1_", x)
}

mean_expr <- function(mat, groups, genes) {
  genes <- intersect(genes, rownames(mat))
  out <- do.call(rbind, lapply(genes, function(g) {
    x <- as.numeric(mat[g, ])
    do.call(rbind, lapply(names(split(seq_along(x), groups)), function(grp) {
      idx <- which(groups == grp)
      data.frame(
        group = grp,
        gene = g,
        mean_expr = mean(x[idx]),
        median_expr = median(x[idx]),
        pct_expr = mean(x[idx] > 0),
        n_cells = length(idx)
      )
    }))
  }))
  rownames(out) <- NULL
  out[, c("group", "gene", "n_cells", "mean_expr", "median_expr", "pct_expr")]
}

score_matrix <- function(mat, signatures) {
  rows <- list()
  for (nm in names(signatures)) {
    genes <- intersect(signatures[[nm]], rownames(mat))
    if (length(genes) == 0) next
    rows[[nm]] <- Matrix::colMeans(mat[genes, , drop = FALSE])
  }
  as.data.frame(rows)
}

message("Loading integrated_oc.RData")
obj <- readRDS("D:/OC_spatiogenomics/infercnv/integrated_oc.RData")
meta <- obj@meta.data
meta$cell_integrated_oc <- rownames(meta)
meta$cell_canonical <- canonical_cell(rownames(meta))
meta$cell_type <- as.character(meta$cell_type)
data.input <- obj@assays$RNA@data

common <- intersect(colnames(data.input), rownames(meta))
data.input <- data.input[, common]
meta <- meta[common, , drop = FALSE]
write.csv(as.data.frame(table(meta$cell_type)), file.path(table_dir, "integrated_oc_cell_type_counts_from_RData.csv"), row.names = FALSE)

message("Summarizing cell-type gene sources")
cell_source <- mean_expr(data.input, meta$cell_type, candidate_genes)
names(cell_source)[names(cell_source) == "group"] <- "cell_type"
cell_source <- cell_source[order(cell_source$gene, -cell_source$mean_expr, -cell_source$pct_expr), ]
write.csv(cell_source, file.path(table_dir, "integrated_oc_candidate_gene_cell_type_expression.csv"), row.names = FALSE)

top_sources <- do.call(rbind, lapply(split(cell_source, cell_source$gene), function(d) {
  d <- d[order(-d$mean_expr, -d$pct_expr), ]
  data.frame(
    gene = d$gene[1],
    top_cell_types = paste(sprintf("%s(mean=%.3f,pct=%.2f,n=%d)", d$cell_type, d$mean_expr, d$pct_expr, d$n_cells), collapse = "; ")
  )
}))
write.csv(top_sources, file.path(table_dir, "integrated_oc_candidate_gene_cell_sources.csv"), row.names = FALSE)

message("Summarizing inferCNV clone expression")
clone_file <- "D:/OC_spatiogenomics/infercnv/integrated_oc_plan_analysis/source_outputs/infercnv_cell_to_subclone_k5.csv"
clone_map <- read.csv(clone_file)
meta_clone <- merge(
  meta,
  clone_map[, c("cell_for_integratedocTcells_style", "cnv_subclone")],
  by.x = "cell_integrated_oc",
  by.y = "cell_for_integratedocTcells_style",
  all.x = FALSE,
  all.y = FALSE
)
clone_cells <- intersect(meta_clone$cell_integrated_oc, colnames(data.input))
meta_clone <- meta_clone[match(clone_cells, meta_clone$cell_integrated_oc), ]
clone_mat <- data.input[, clone_cells, drop = FALSE]
clone_gene <- mean_expr(clone_mat, meta_clone$cnv_subclone, candidate_genes)
names(clone_gene)[names(clone_gene) == "group"] <- "cnv_subclone"
write.csv(clone_gene, file.path(table_dir, "infercnv_clone_candidate_gene_expression_summary.csv"), row.names = FALSE)

clone_scores <- score_matrix(clone_mat, signature_list)
clone_scores$cnv_subclone <- meta_clone$cnv_subclone
clone_score_summary <- do.call(rbind, lapply(setdiff(colnames(clone_scores), "cnv_subclone"), function(sig) {
  do.call(rbind, lapply(split(clone_scores[[sig]], clone_scores$cnv_subclone), function(x) {
    data.frame(signature = sig, mean_score = mean(x), median_score = median(x), n_cells = length(x))
  }))
}))
clone_score_summary$cnv_subclone <- rownames(clone_score_summary)
rownames(clone_score_summary) <- NULL
clone_score_summary <- clone_score_summary[, c("cnv_subclone", "signature", "n_cells", "mean_score", "median_score")]
write.csv(clone_score_summary, file.path(table_dir, "infercnv_clone_candidate_signature_score_summary.csv"), row.names = FALSE)

clone_stats <- do.call(rbind, lapply(candidate_genes, function(g) {
  if (!g %in% rownames(clone_mat)) return(NULL)
  x <- as.numeric(clone_mat[g, ])
  groups <- meta_clone$cnv_subclone
  kw <- kruskal.test(split(x, groups))
  means <- tapply(x, groups, mean)
  data.frame(
    gene = g,
    n_clones = length(means),
    max_clone = names(which.max(means)),
    min_clone = names(which.min(means)),
    max_mean = max(means),
    min_mean = min(means),
    delta_max_min = max(means) - min(means),
    kruskal_p = kw$p.value
  )
}))
clone_stats$fdr <- p.adjust(clone_stats$kruskal_p, method = "BH")
write.csv(clone_stats, file.path(table_dir, "infercnv_clone_candidate_gene_specificity_tests.csv"), row.names = FALSE)

message("Running CellChat by metadata cell_type")
cellchat_rds <- file.path(table_dir, "cellchat_integrated_oc_cell_type.rds")
if (file.exists(cellchat_rds)) {
  message("Using existing CellChat RDS: ", cellchat_rds)
  cellchat <- readRDS(cellchat_rds)
} else {
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type", do.sparse = TRUE)
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  message("CellChat signaling genes retained: ", nrow(cellchat@data.signaling))
  message("Candidate genes in data.signaling: ", paste(intersect(candidate_genes, rownames(cellchat@data.signaling)), collapse = ", "))
  message("Candidate genes in RNA matrix: ", paste(intersect(candidate_genes, rownames(data.input)), collapse = ", "))
  lr_all <- CellChatDB$interaction
  lr_keep <- vapply(seq_len(nrow(lr_all)), function(i) {
    genes <- extractGeneSubsetFromPair(lr_all[i, , drop = FALSE], object = cellchat, combined = TRUE)
    all(genes %in% rownames(cellchat@data.signaling))
  }, logical(1))
  lr_use <- lr_all[lr_keep, , drop = FALSE]
  rownames(lr_use) <- lr_use$interaction_name
  message("CellChat LR pairs retained after gene availability filter: ", nrow(lr_use), " / ", nrow(lr_all))
  message("Candidate LR pairs retained: ", paste(intersect(c("CXCL12_CXCR4", "MDK_SDC4"), rownames(lr_use)), collapse = ", "))
  cellchat@LR$LRsig <- lr_use
  cellchat <- computeCommunProb(cellchat, LR.use = lr_use, raw.use = TRUE, population.size = TRUE, nboot = 100, seed.use = 1)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, cellchat_rds)
}

comm <- subsetCommunication(cellchat, thresh = 1)
write.csv(comm, file.path(table_dir, "cellchat_cell_type_all_communications.csv"), row.names = FALSE)
candidate_names <- intersect(c("CXCL12_CXCR4", "MDK_SDC4"), dimnames(cellchat@net$prob)[[3]])
candidate_axes <- do.call(rbind, lapply(candidate_names, function(ax) {
  prob <- cellchat@net$prob[, , ax]
  pval <- cellchat@net$pval[, , ax]
  grid <- expand.grid(source = rownames(prob), target = colnames(prob), stringsAsFactors = FALSE)
  grid$ligand <- sub("_.*$", "", ax)
  grid$receptor <- sub("^.*_", "", ax)
  grid$prob <- as.numeric(prob[cbind(grid$source, grid$target)])
  grid$pval <- as.numeric(pval[cbind(grid$source, grid$target)])
  grid$interaction_name <- ax
  grid$pathway_name <- ifelse(ax == "CXCL12_CXCR4", "CXCL", "MK")
  grid$axis <- ax
  grid$significant <- grid$pval < 0.05 & grid$prob > 0
  grid
}))
candidate_axes <- candidate_axes[order(candidate_axes$axis, -candidate_axes$prob, candidate_axes$pval), ]
candidate_axes_sig <- candidate_axes[candidate_axes$significant, , drop = FALSE]
write.csv(candidate_axes, file.path(table_dir, "singlecell_LR_axis_hits_CXCL12_CXCR4_MDK_SDC4.csv"), row.names = FALSE)
write.csv(candidate_axes, file.path(table_dir, "cellchat_cell_type_candidate_axes.csv"), row.names = FALSE)

png(file.path(figure_dir, "cellchat_cell_type_candidate_axis_prob.png"), width = 2200, height = 1400, res = 180)
plot_axes <- candidate_axes[candidate_axes$prob > 0, , drop = FALSE]
if (nrow(plot_axes) > 0) {
  op <- par(mar = c(6, 12, 4, 2))
  plot_axes <- head(plot_axes, 30)
  labels <- paste(plot_axes$source, "->", plot_axes$target, plot_axes$axis, ifelse(plot_axes$significant, "*", ""))
  barplot(plot_axes$prob, names.arg = labels, horiz = TRUE, las = 1,
          col = ifelse(plot_axes$axis == "MDK_SDC4", "#2F6C8F", "#A4514F"),
          xlab = "CellChat probability", main = "Candidate CellChat axes by metadata cell_type")
  par(op)
} else {
  plot.new()
  text(0.5, 0.5, "No CXCL12-CXCR4 or MDK-SDC4 interactions passed CellChat threshold")
}
dev.off()

report <- c(
  "# CellChat rerun with integrated_oc metadata cell_type",
  "",
  sprintf("Cells used: %d; genes in RNA data: %d; grouping field: metadata$cell_type.", ncol(data.input), nrow(data.input)),
  "",
  "Cell-type counts:",
  paste(capture.output(print(table(meta$cell_type))), collapse = "\n"),
  "",
  sprintf("CellChat communications passing threshold: %d.", nrow(comm)),
  sprintf("Candidate CXCL12-CXCR4 / MDK-SDC4 significant directed pairs: %d / %d total directed pairs.", nrow(candidate_axes_sig), nrow(candidate_axes)),
  "",
  "Top candidate axes:",
  paste(capture.output(print(head(candidate_axes[, intersect(c('source','target','ligand','receptor','prob','pval','pathway_name','axis','significant'), colnames(candidate_axes))], 20))), collapse = "\n"),
  "",
  "Clone specificity tests:",
  paste(capture.output(print(clone_stats)), collapse = "\n")
)
writeLines(report, file.path(report_dir, "cellchat_celltype_rerun_report.md"), useBytes = TRUE)
message("Done CellChat rerun")
