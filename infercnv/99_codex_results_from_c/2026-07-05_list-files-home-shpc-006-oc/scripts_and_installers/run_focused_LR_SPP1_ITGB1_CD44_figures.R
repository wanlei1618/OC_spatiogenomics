############################################################
## focused LR axes visualization: SPP1 -> ITGB1/CD44
## Outputs 4 publication-ready figures under:
## D:/OC_spatiogenomics/infercnv/LR_SPP1_ITGB1_CD44_figures
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tibble)
  library(grid)
})

base_dir <- "D:/OC_spatiogenomics/infercnv"
lr_file <- file.path(base_dir, "integrated_oc_plan_analysis/tables/focused_LR_axes_involving_cnv_subclones.csv")
outdir <- file.path(base_dir, "LR_SPP1_ITGB1_CD44_figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_lr_table <- function(file) {
  x0 <- tryCatch(
    read.csv(file, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (!is.null(x0)) return(x0)
  encodings <- c("UTF-8", "UTF-8-BOM", "GB18030", "GBK", "latin1")
  for (enc in encodings) {
    x <- tryCatch(
      read.csv(file, fileEncoding = enc, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(x)) return(x)
  }
  stop("Failed to read CSV: ", file)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

lr_raw <- read_lr_table(lr_file)
required_cols <- c("source_group", "target_group", "ligand", "receptor", "score",
                   "ligand_avg", "receptor_avg", "ligand_pct", "receptor_pct")
missing_cols <- setdiff(required_cols, colnames(lr_raw))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

lr <- lr_raw %>%
  mutate(
    across(c(source_group, target_group, ligand, receptor), as.character),
    score = safe_num(score),
    ligand_avg = safe_num(ligand_avg),
    receptor_avg = safe_num(receptor_avg),
    ligand_pct = safe_num(ligand_pct),
    receptor_pct = safe_num(receptor_pct),
    focus_axis = if ("focus_axis" %in% colnames(.)) as.character(focus_axis) else paste(ligand, receptor, sep = "_")
  ) %>%
  filter(!is.na(score))

target_clones <- sprintf("CNV_Subclone_%02d", 1:5)
priority_clones <- c("CNV_Subclone_02", "CNV_Subclone_04")
myeloid_pattern <- "Myeloid|Macro|Macrophage|Monocyte|DC"
axis_of_interest <- c("SPP1_ITGB1", "SPP1_CD44", "APOE_LRP1", "MIF_CD74", "MIF_CD74_CXCR4",
                      "TGFB1_TGFBR1", "TGFB1_TGFBR2", "VEGFA_KDR", "VEGFA_FLT1",
                      "CSF1_CSF1R", "CD47_SIRPA")

lr_myeloid_to_cnv <- lr %>%
  filter(grepl(myeloid_pattern, source_group, ignore.case = TRUE), target_group %in% target_clones)

spp1_lr <- lr_myeloid_to_cnv %>%
  filter(ligand == "SPP1", receptor %in% c("ITGB1", "CD44"))

if (nrow(spp1_lr) == 0) stop("No SPP1 -> ITGB1/CD44 interaction rows found.")

top_n_sources <- 10
source_rank <- spp1_lr %>%
  group_by(source_group) %>%
  summarise(
    max_score = max(score, na.rm = TRUE),
    mean_score = mean(score, na.rm = TRUE),
    max_ligand_avg = max(ligand_avg, na.rm = TRUE),
    max_ligand_pct = max(ligand_pct, na.rm = TRUE),
    n_targets = n_distinct(target_group),
    .groups = "drop"
  ) %>%
  arrange(desc(max_score), desc(max_ligand_pct))

source_order <- source_rank %>% slice_head(n = top_n_sources) %>% pull(source_group)
target_order <- target_clones[target_clones %in% unique(lr_myeloid_to_cnv$target_group)]
receptor_order <- c("ITGB1", "CD44")

############################
## Figure 1: bubble heatmap
############################

plot_df <- spp1_lr %>%
  filter(source_group %in% source_order) %>%
  mutate(
    source_group = factor(source_group, levels = rev(source_order)),
    target_group = factor(target_group, levels = target_order),
    receptor = factor(receptor, levels = receptor_order),
    axis_label = factor(paste0("SPP1-", receptor), levels = c("SPP1-ITGB1", "SPP1-CD44"))
  )

bg_df <- tidyr::expand_grid(
  source_group = factor(source_order, levels = rev(source_order)),
  target_group = factor(priority_clones, levels = target_order),
  axis_label = factor(levels(plot_df$axis_label), levels = levels(plot_df$axis_label))
)

p1 <- ggplot() +
  geom_tile(data = bg_df, aes(x = target_group, y = source_group),
            fill = "grey92", color = NA, width = 0.95, height = 0.95) +
  geom_point(data = plot_df,
             aes(x = target_group, y = source_group, size = score, fill = score),
             shape = 21, color = "black", stroke = 0.25, alpha = 0.96) +
  facet_wrap(~ axis_label, nrow = 1) +
  scale_size_continuous(range = c(1.8, 10), name = "LR score") +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "LR score") +
  labs(
    title = "SPP1-myeloid/macrophage signaling to CNV subclones",
    subtitle = "Bubble size/color = LR score; grey background highlights CNV_Subclone_02/04",
    x = "Target CNV subclone", y = "Source myeloid/macrophage group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey88", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    panel.grid = element_line(color = "grey91")
  )

ggsave(file.path(outdir, "Fig1_SPP1_ITGB1_CD44_bubble_heatmap.pdf"), p1, width = 12, height = 6.5)
ggsave(file.path(outdir, "Fig1_SPP1_ITGB1_CD44_bubble_heatmap.png"), p1, width = 12, height = 6.5, dpi = 300)

############################
## Figure 2: expression support
############################

ligand_expr_df <- spp1_lr %>%
  filter(source_group %in% source_order) %>%
  group_by(source_group, ligand) %>%
  summarise(
    ligand_avg = max(ligand_avg, na.rm = TRUE),
    ligand_pct = max(ligand_pct, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source_group = factor(source_group, levels = rev(source_order)),
         ligand = factor(ligand, levels = "SPP1"))

p2_left <- ggplot(ligand_expr_df, aes(x = ligand, y = source_group)) +
  geom_point(aes(size = ligand_pct, fill = ligand_avg),
             shape = 21, color = "black", stroke = 0.25) +
  scale_size_continuous(range = c(2.2, 9), labels = percent_format(accuracy = 1), name = "Ligand pct") +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Ligand avg") +
  labs(title = "Source-side SPP1 expression", x = NULL, y = "Source group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold"))

receptor_expr_df <- spp1_lr %>%
  group_by(target_group, receptor) %>%
  summarise(
    receptor_avg = max(receptor_avg, na.rm = TRUE),
    receptor_pct = max(receptor_pct, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(target_group = factor(target_group, levels = target_order),
         receptor = factor(receptor, levels = receptor_order))

bg_receptor_df <- tidyr::expand_grid(
  target_group = factor(priority_clones, levels = target_order),
  receptor = factor(receptor_order, levels = receptor_order)
)

p2_right <- ggplot() +
  geom_tile(data = bg_receptor_df, aes(x = target_group, y = receptor),
            fill = "grey92", color = NA, width = 0.95, height = 0.95) +
  geom_point(data = receptor_expr_df,
             aes(x = target_group, y = receptor, size = receptor_pct, fill = receptor_avg),
             shape = 21, color = "black", stroke = 0.25) +
  scale_size_continuous(range = c(2.2, 9), labels = percent_format(accuracy = 1), name = "Receptor pct") +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Receptor avg") +
  labs(title = "Target-side ITGB1/CD44 expression",
       subtitle = "Grey background highlights CNV_Subclone_02/04",
       x = "Target CNV subclone", y = "Receptor") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

p2 <- p2_left + p2_right + plot_layout(widths = c(1.15, 1.35))
ggsave(file.path(outdir, "Fig2_SPP1_ligand_receptor_expression_dotplot.pdf"), p2, width = 13, height = 6.5)
ggsave(file.path(outdir, "Fig2_SPP1_ligand_receptor_expression_dotplot.png"), p2, width = 13, height = 6.5, dpi = 300)

############################
## Figure 3: candidate axis priority
############################

axis_rank_df <- lr_myeloid_to_cnv %>%
  mutate(is_priority_target = target_group %in% priority_clones) %>%
  group_by(focus_axis, ligand, receptor) %>%
  summarise(
    n_interactions = n(),
    n_sources = n_distinct(source_group),
    n_targets = n_distinct(target_group),
    total_score = sum(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    mean_score = mean(score, na.rm = TRUE),
    top5_mean_score = mean(sort(score, decreasing = TRUE)[seq_len(min(5, n()))], na.rm = TRUE),
    priority_score_sum = sum(score[is_priority_target], na.rm = TRUE),
    priority_fraction = priority_score_sum / total_score,
    mean_ligand_pct = mean(ligand_pct, na.rm = TRUE),
    mean_receptor_pct = mean(receptor_pct, na.rm = TRUE),
    expr_support = sqrt(pmax(mean_ligand_pct, 0) * pmax(mean_receptor_pct, 0)),
    .groups = "drop"
  ) %>%
  mutate(
    validation_priority_index = top5_mean_score * priority_fraction * expr_support * log1p(n_sources + n_targets),
    axis_class = case_when(
      ligand == "SPP1" & receptor %in% c("ITGB1", "CD44") ~ "SPP1 axis",
      focus_axis %in% axis_of_interest ~ "Other candidate axis",
      TRUE ~ "Other LR axis"
    )
  ) %>%
  arrange(desc(validation_priority_index))

priority_plot_df <- axis_rank_df %>%
  filter(row_number() <= 18 | focus_axis %in% c("SPP1_ITGB1", "SPP1_CD44", "APOE_LRP1", "MIF_CD74_CXCR4")) %>%
  distinct(focus_axis, .keep_all = TRUE) %>%
  arrange(validation_priority_index) %>%
  mutate(focus_axis = factor(focus_axis, levels = focus_axis))

p3 <- ggplot(priority_plot_df, aes(x = validation_priority_index, y = focus_axis)) +
  geom_col(aes(fill = axis_class), width = 0.65, alpha = 0.9) +
  geom_point(aes(size = max_score, color = priority_fraction), stroke = 0.25) +
  scale_fill_manual(values = c("SPP1 axis" = "#B2182B", "Other candidate axis" = "#4D9221", "Other LR axis" = "grey70"),
                    name = "Axis class") +
  scale_color_gradient(low = "#92C5DE", high = "#B2182B",
                       labels = percent_format(accuracy = 1),
                       name = "Score fraction\nto 02/04") +
  scale_size_continuous(range = c(2.5, 8), name = "Max score") +
  labs(
    title = "Prioritization of candidate ligand-receptor axes",
    subtitle = "Index combines strength, CNV_Subclone_02/04 specificity, expression support,\nand source/target coverage",
    x = "Validation priority index", y = "Focused LR axis"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "right")

ggsave(file.path(outdir, "Fig3_LR_axis_priority_plot.pdf"), p3, width = 10, height = 6.2)
ggsave(file.path(outdir, "Fig3_LR_axis_priority_plot.png"), p3, width = 10, height = 6.2, dpi = 300)

############################
## Figure 4: simplified network
############################

network_sources <- spp1_lr %>%
  filter(target_group %in% priority_clones) %>%
  group_by(source_group) %>%
  summarise(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_score)) %>%
  slice_head(n = 7) %>%
  pull(source_group)

network_edges <- spp1_lr %>%
  filter(source_group %in% network_sources, target_group %in% priority_clones) %>%
  mutate(edge_label = paste0("SPP1-", receptor))

source_nodes <- tibble(
  node = network_sources,
  node_type = "SPP1+ myeloid/macrophage",
  x = 1,
  y = seq_along(network_sources)
)
target_nodes <- tibble(
  node = rev(priority_clones),
  node_type = "CNV target clone",
  x = 3,
  y = seq(from = 2, to = length(network_sources) - 1, length.out = length(priority_clones))
)
nodes <- bind_rows(source_nodes, target_nodes)

edges_plot <- network_edges %>%
  left_join(nodes %>% select(source_group = node, x_from = x, y_from = y), by = "source_group") %>%
  left_join(nodes %>% select(target_group = node, x_to = x, y_to = y), by = "target_group")

p4 <- ggplot() +
  geom_curve(
    data = edges_plot,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
        size = score, linetype = receptor, alpha = score),
    curvature = 0.18, color = "grey25"
  ) +
  geom_point(data = nodes, aes(x = x, y = y, shape = node_type, fill = node_type),
             size = 7.2, color = "black") +
  geom_text(data = source_nodes, aes(x = x - 0.08, y = y, label = node),
            hjust = 1, size = 3.5) +
  geom_text(data = target_nodes, aes(x = x + 0.08, y = y, label = node),
            hjust = 0, size = 4.1, fontface = "bold") +
  annotate("text", x = 2, y = length(network_sources) + 0.35,
           label = "SPP1 -> ITGB1/CD44", fontface = "bold", size = 4.2) +
  scale_size_continuous(range = c(0.35, 2.9), name = "LR score") +
  scale_alpha_continuous(range = c(0.25, 0.95), guide = "none") +
  scale_shape_manual(values = c("SPP1+ myeloid/macrophage" = 21, "CNV target clone" = 24), name = NULL) +
  scale_fill_manual(values = c("SPP1+ myeloid/macrophage" = "#92C5DE", "CNV target clone" = "#F4A582"), name = NULL) +
  labs(
    title = "Simplified SPP1-myeloid/macrophage to CNV_Subclone_02/04 network",
    subtitle = "Edges are restricted to SPP1-ITGB1/CD44 interactions",
    x = NULL, y = NULL
  ) +
  coord_cartesian(xlim = c(0.05, 3.95), ylim = c(0.5, length(network_sources) + 0.65), clip = "off") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 95, 10, 205))

ggsave(file.path(outdir, "Fig4_SPP1_myeloid_to_CNV_network.pdf"), p4, width = 12, height = 6.5)
ggsave(file.path(outdir, "Fig4_SPP1_myeloid_to_CNV_network.png"), p4, width = 12, height = 6.5, dpi = 300)

############################
## Export plotting tables and run summary
############################

write.csv(plot_df, file.path(outdir, "Fig1_SPP1_bubble_plot_data.csv"), row.names = FALSE)
write.csv(ligand_expr_df, file.path(outdir, "Fig2_source_SPP1_expression_data.csv"), row.names = FALSE)
write.csv(receptor_expr_df, file.path(outdir, "Fig2_target_ITGB1_CD44_expression_data.csv"), row.names = FALSE)
write.csv(axis_rank_df, file.path(outdir, "Fig3_LR_axis_priority_summary.csv"), row.names = FALSE)
write.csv(network_edges, file.path(outdir, "Fig4_network_edge_data.csv"), row.names = FALSE)

summary_lines <- c(
  "# Focused LR SPP1-ITGB1/CD44 visualization",
  "",
  paste0("Input: ", normalizePath(lr_file, winslash = "/")),
  paste0("Output: ", normalizePath(outdir, winslash = "/")),
  "",
  paste0("Total LR rows: ", nrow(lr)),
  paste0("Myeloid/macrophage -> CNV rows: ", nrow(lr_myeloid_to_cnv)),
  paste0("SPP1 -> ITGB1/CD44 rows: ", nrow(spp1_lr)),
  "",
  "Figures:",
  "- Fig1_SPP1_ITGB1_CD44_bubble_heatmap.pdf/png",
  "- Fig2_SPP1_ligand_receptor_expression_dotplot.pdf/png",
  "- Fig3_LR_axis_priority_plot.pdf/png",
  "- Fig4_SPP1_myeloid_to_CNV_network.pdf/png"
)
writeLines(summary_lines, file.path(outdir, "README_SPP1_ITGB1_CD44_figures.md"))
cat(paste(summary_lines, collapse = "\n"), "\n")
