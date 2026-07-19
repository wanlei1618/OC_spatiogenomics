# Usage:
# 1. Set full_seurat_rds to the full OO/OC Seurat object that contains both malignant
#    cells and the T/NK cells in integratedocTcells.
# 2. Run with an R environment where Seurat and CellChat are installed.

library(Seurat)
library(CellChat)
library(patchwork)

full_seurat_rds <- "D:/OC_spatiogenomics/path_to_full_OO_or_OC_seurat_object.rds"
subclone_csv <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs/infercnv_cell_to_subclone_k5.csv"
out_dir <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs/cellchat_subclone_vs_tcelltype"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(full_seurat_rds)
subclone <- read.csv(subclone_csv, stringsAsFactors = FALSE)

obj$cnv_subclone <- NA_character_
obj$cnv_subclone[match(subclone$cell, colnames(obj))] <- subclone$cnv_subclone
obj$cnv_subclone[match(subclone$cell_for_integratedocTcells_style, colnames(obj))] <- subclone$cnv_subclone

stopifnot("cell_type" %in% colnames(obj@meta.data))

obj$cellchat_group <- ifelse(
  !is.na(obj$cnv_subclone),
  obj$cnv_subclone,
  as.character(obj$cell_type)
)

keep_groups <- c(unique(subclone$cnv_subclone), unique(as.character(obj$cell_type)))
obj_use <- subset(obj, subset = cellchat_group %in% keep_groups)

data.input <- GetAssayData(obj_use, assay = "RNA", slot = "data")
meta <- data.frame(group = obj_use$cellchat_group, row.names = colnames(obj_use))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
CellChatDB <- CellChatDB.human
cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path(out_dir, "cellchat_subclone_vs_tcelltype.rds"))
write.csv(subsetCommunication(cellchat), file.path(out_dir, "cellchat_subclone_vs_tcelltype_communications.csv"), row.names = FALSE)

pdf(file.path(out_dir, "cellchat_interaction_count.pdf"), width = 8, height = 7)
netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)), weight.scale = TRUE, label.edge = FALSE)
dev.off()

pdf(file.path(out_dir, "cellchat_interaction_weight.pdf"), width = 8, height = 7)
netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)), weight.scale = TRUE, label.edge = FALSE)
dev.off()
