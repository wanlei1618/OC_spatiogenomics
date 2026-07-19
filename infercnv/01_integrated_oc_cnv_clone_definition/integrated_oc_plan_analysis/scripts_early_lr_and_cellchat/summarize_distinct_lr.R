out_dir <- "C:/Users/chenfy12/Documents/Codex/2026-07-05/list-files-home-shpc-006-oc/outputs"
x <- read.csv(file.path(out_dir, "lr_top_interactions_involving_cnv_subclones.csv"), stringsAsFactors = FALSE)
x$lr_pair <- paste(x$ligand, x$receptor, sep = " -> ")
x <- x[order(x$score, decreasing = TRUE), ]
distinct_lr <- x[!duplicated(x$lr_pair), ]
write.csv(
  head(distinct_lr, 200),
  file.path(out_dir, "lr_top_distinct_LR_pairs_involving_cnv_subclones.csv"),
  row.names = FALSE
)
print(head(distinct_lr[, c("source_group", "target_group", "ligand", "receptor", "score", "ligand_pct", "receptor_pct")], 30))
