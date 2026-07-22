options(stringsAsFactors = FALSE, warn = 1)

required <- c("yaml", "data.table", "Matrix", "Seurat", "SeuratObject", "ggplot2")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else "."
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config(); cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash = "/", mustWork = TRUE)
out <- file.path(data_root, "diagnostics_v3_remaining_datasets", "GSE151214")
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly = TRUE)
if (dir.exists(out) && !replace_generated) stop("Output already exists; refusing to overwrite: ", out)
if (!dir.exists(out)) dir.create(out, recursive = TRUE, showWarnings = FALSE)
seed <- as.integer(cfg$project$random_seed %||% 20260718L); set.seed(seed)
expected_cells <- 48837L
write_gz <- function(x, path) data.table::fwrite(x, path, na = "NA", compress = "gzip")
upper <- function(x) toupper(as.character(x))

panels <- list(
  Secretory_epithelial = c("EPCAM","PAX8","OVGP1","KRT8","KRT18","KRT19","KRT7","WFDC2","MUC1","SLPI"),
  Ciliated_epithelial = c("FOXJ1","PIFO","TPPP3","CAPS","C20ORF85","RSPH1","CFAP126","CFAP157"),
  Transitional_epithelial = c("KRT13","KRT4","KRT17","KRT5","KRT14","KRT15","CLDN4","SFN"),
  T_cell = c("CD3D","CD3E","CD3G","TRBC1","TRBC2","CD2","IL7R","LCK"),
  NK_cell = c("NKG7","GNLY","KLRD1","KLRF1","PRF1","XCL1","XCL2","FGFBP2"),
  Gamma_delta_T = c("TRDC","TRGC1","TRGC2","CD3D","CD3E","KLRB1"),
  B_cell = c("MS4A1","CD79A","CD79B","CD19","CD74","BANK1"),
  Plasma_cell = c("MZB1","JCHAIN","SDC1","DERL3","XBP1","IGHG1"),
  Monocyte = c("S100A8","S100A9","FCN1","VCAN","CTSS","LYZ"),
  Macrophage = c("C1QA","C1QB","C1QC","APOE","MRC1","CD68","TREM2","SPP1"),
  Dendritic = c("FCER1A","CD1C","CLEC10A","CLEC9A","XCR1","LILRA4","CLEC4C"),
  Mast = c("TPSAB1","TPSB2","CPA3","KIT","MS4A2","HDC"),
  Fibroblast = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
  Smooth_muscle = c("ACTA2","MYH11","TAGLN","CNN1","MYL9","DES"),
  Pericyte = c("RGS5","CSPG4","PDGFRB","KCNJ8","ABCC9","NOTCH3"),
  Blood_endothelial = c("PECAM1","VWF","KDR","EMCN","RAMP2","PLVAP","CLDN5"),
  Lymphatic_endothelial = c("PDPN","LYVE1","FLT4","CCL21","PROX1","CLDN5")
)
high_specific <- list(
  Secretory_epithelial=c("OVGP1","PAX8","WFDC2","EPCAM"), Ciliated_epithelial=c("FOXJ1","PIFO","TPPP3","CAPS"),
  Transitional_epithelial=c("KRT13","KRT4","KRT17","KRT5"), T_cell=c("CD3D","CD3E","TRBC1","TRBC2"),
  NK_cell=c("GNLY","KLRD1","XCL1","XCL2"), Gamma_delta_T=c("TRDC","TRGC1","TRGC2"),
  B_cell=c("MS4A1","CD79A","CD79B","CD19"), Plasma_cell=c("MZB1","JCHAIN","SDC1","DERL3"),
  Monocyte=c("S100A8","S100A9","FCN1","VCAN"), Macrophage=c("C1QA","C1QB","C1QC","MRC1"),
  Dendritic=c("FCER1A","CD1C","CLEC9A","XCR1","LILRA4"), Mast=c("TPSAB1","TPSB2","CPA3","KIT"),
  Fibroblast=c("COL1A1","COL1A2","DCN","LUM"), Smooth_muscle=c("MYH11","TAGLN","CNN1","ACTA2"),
  Pericyte=c("RGS5","CSPG4","PDGFRB","KCNJ8"), Blood_endothelial=c("PECAM1","VWF","KDR","EMCN"),
  Lymphatic_endothelial=c("PDPN","LYVE1","FLT4","CCL21")
)
state_panels <- list(
  IFN_response=c("ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","OAS2","IRF7"),
  Stress_response=c("FOS","JUN","JUNB","ATF3","HSPA1A","HSPA1B","DDIT3"),
  Activation=c("CD69","FOS","JUNB","TNFRSF9","IL2RA","NFKBIA")
)
type_family <- function(x) {
  y <- as.character(x)
  y[y %in% c("Secretory_epithelial","Ciliated_epithelial","Transitional_epithelial")] <- "Epithelial"
  y[y %in% c("T_cell","NK_cell","Gamma_delta_T")] <- "T_NK"
  y[y %in% c("B_cell","Plasma_cell")] <- "B_Plasma"
  y[y %in% c("Monocyte","Macrophage","Dendritic")] <- "Myeloid"
  y[y %in% c("Fibroblast","Smooth_muscle","Pericyte")] <- "Stromal"
  y[y %in% c("Blood_endothelial","Lymphatic_endothelial")] <- "Endothelial"
  y
}

object_path <- file.path(data_root, "GSE151214", "objects", "GSE151214_preannotation.rds")
marker_path <- file.path(data_root, "GSE151214", "03_markers", "all_cluster_markers.csv.gz")
obj <- readRDS(object_path)
if (ncol(obj) != expected_cells) stop("Expected 48,837 post-QC/doublet cells")
md <- obj[[]]
sample_col <- if ("sample_id" %in% names(md)) "sample_id" else if ("orig.ident" %in% names(md)) "orig.ident" else stop("No sample identifier")
cluster_col <- if ("RNA_snn_res.0.6" %in% names(md)) "RNA_snn_res.0.6" else "seurat_clusters"
if (!cluster_col %in% names(md)) stop("No cluster field")
obj$sample_id_v3 <- as.character(md[[sample_col]])
obj$cluster_v3 <- as.character(md[[cluster_col]])

mt_genes <- grep("^MT-", rownames(obj), value = TRUE, ignore.case = TRUE)
pct_col <- grep("^percent[._]?mt$", names(md), value = TRUE, ignore.case = TRUE)
if (!length(pct_col)) stop("No mitochondrial percentage column for audit")
pct <- as.numeric(md[[pct_col[[1L]]]])
mt_audit <- data.table::rbindlist(list(
  data.table::data.table(scope="global", sample_id="ALL", n_cells=ncol(obj), n_mt_features=length(mt_genes),
                         mt_min=min(pct,na.rm=TRUE), mt_median=median(pct,na.rm=TRUE), mt_mean=mean(pct,na.rm=TRUE),
                         mt_p95=quantile(pct,.95,na.rm=TRUE), mt_max=max(pct,na.rm=TRUE)),
  data.table::data.table(sample_id=obj$sample_id_v3, pct=pct)[, .(scope="sample", n_cells=.N,
    n_mt_features=length(mt_genes), mt_min=min(pct,na.rm=TRUE), mt_median=median(pct,na.rm=TRUE),
    mt_mean=mean(pct,na.rm=TRUE), mt_p95=quantile(pct,.95,na.rm=TRUE), mt_max=max(pct,na.rm=TRUE)), by=sample_id]
), fill=TRUE)
mt_pass <- length(mt_genes) >= 10L && all(is.finite(pct)) && max(pct) < 50 && median(pct) > 0
mt_audit[, `:=`(audit_status=if (mt_pass) "PASS_EXISTING_QC_VALID" else "FAIL_REVIEW_REQUIRED",
                action=if (mt_pass) "No mitochondrial QC rerun" else "Review required")]
data.table::fwrite(mt_audit, file.path(out,"mt_feature_audit.csv"), na="NA")
if (!mt_pass) stop("Mitochondrial feature audit failed")

markers <- data.table::fread(marker_path, showProgress=FALSE)
markers <- markers[is.finite(p_val_adj) & p_val_adj < .05 & is.finite(avg_log2FC) & avg_log2FC > .25 &
                     is.finite(pct.1) & pct.1 >= .2]
markers[, `:=`(cluster=as.character(cluster), gene_upper=upper(gene))]
score_features <- lapply(panels, function(x) intersect(x, rownames(obj)))
if (any(lengths(score_features) < 2L)) stop("Insufficient genes in lineage panel")
obj <- AddModuleScore(obj, features=unname(score_features), name="lineage_v3", nbin=24, ctrl=25, seed=seed, search=FALSE)
score_cols <- paste0("lineage_v3", seq_along(panels)); names(score_cols) <- names(panels)
scores <- data.table::as.data.table(obj[[]], keep.rownames="cell_id")
scores[, cluster := as.character(cluster_v3)]
means <- scores[, lapply(.SD, mean), by=cluster, .SDcols=unname(score_cols)]
data.table::setnames(means, unname(score_cols), names(score_cols))

decision_rows <- lapply(sort(unique(scores$cluster)), function(cl) {
  cm <- markers[cluster==cl]
  cand <- data.table::rbindlist(lapply(names(panels), function(tp) {
    hit <- cm[gene_upper %in% upper(panels[[tp]])]
    hi <- hit[gene_upper %in% upper(high_specific[[tp]])]
    data.table::data.table(candidate=tp, support_n=data.table::uniqueN(hit$gene_upper),
      high_specific_n=data.table::uniqueN(hi$gene_upper), module_mean=means[cluster==cl, get(tp)],
      canonical_markers=paste(unique(hit[order(-avg_log2FC)]$gene),collapse=";"))
  }))
  cand[, eligible := support_n>=3L | high_specific_n>=2L]
  data.table::setorder(cand, -eligible, -support_n, -high_specific_n, -module_mean)
  top <- cand[1L]; second <- cand[2L]
  conflict <- top$eligible && second$eligible && type_family(top$candidate) != type_family(second$candidate) &&
    second$support_n >= top$support_n-1L && second$module_mean >= top$module_mean-.1
  data.table::data.table(cluster=cl, final_cell_type=if (!top$eligible || conflict) "Unresolved" else top$candidate,
    canonical_support_n=top$support_n, high_specific_support_n=top$high_specific_n,
    canonical_markers=top$canonical_markers, second_candidate=second$candidate, incompatible_lineage_program=conflict)
})
cluster_table <- data.table::rbindlist(decision_rows)
# Resolve normal-reference sublineages only when their defining marker combinations
# are significant; no cluster number is used in these rules.
gamma_clusters <- markers[, .(marker_rule = all(c("TRDC","TRGC1","TRGC2") %in% gene_upper) &&
  any(c("CD3D","CD3E") %in% gene_upper)), by=cluster][marker_rule == TRUE, cluster]
pericyte_clusters <- markers[, .(marker_rule = "RGS5" %in% gene_upper &&
  sum(c("RGS5","NOTCH3","PDGFRB","CSPG4","KCNJ8","ABCC9") %in% gene_upper) >= 3L), by=cluster][marker_rule == TRUE, cluster]
cluster_table[cluster %in% gamma_clusters, `:=`(final_cell_type="Gamma_delta_T", incompatible_lineage_program=FALSE)]
cluster_table[cluster %in% pericyte_clusters, `:=`(final_cell_type="Pericyte", incompatible_lineage_program=FALSE)]
dominance <- dominance_metrics(obj$cluster_v3, obj$sample_id_v3, "GSE151214", "normal_reference_uncorrected")
data.table::setnames(dominance,"normalized_shannon_entropy","sample_entropy")
dominance[, patient_enriched := dominant_sample_fraction >= .8]
cluster_table <- merge(cluster_table, dominance[,.(cluster,n_cells,n_samples,dominant_sample,dominant_sample_fraction,sample_entropy,patient_enriched)], by="cluster", all.x=TRUE)

state_features <- lapply(state_panels, function(x) intersect(x, rownames(obj)))
obj <- AddModuleScore(obj, features=unname(state_features), name="state_v3", nbin=24, ctrl=25, seed=seed, search=FALSE)
state_cols <- paste0("state_v3",seq_along(state_panels)); names(state_cols)<-names(state_panels)
state_md <- data.table::as.data.table(obj[[]], keep.rownames="cell_id"); state_md[,cluster:=as.character(cluster_v3)]
state_means <- state_md[,lapply(.SD,mean),by=cluster,.SDcols=unname(state_cols)]
data.table::setnames(state_means,unname(state_cols),names(state_cols))
cluster_table <- merge(cluster_table,state_means,by="cluster",all.x=TRUE)
cluster_table[, cell_state := vapply(cluster,function(cl) {
  sm<-state_means[cluster==cl]; z<-character()
  for (nm in names(state_panels)) if (sm[[nm]]>0 && data.table::uniqueN(markers[cluster==cl & gene_upper %in% upper(state_panels[[nm]])]$gene_upper)>=3L) z<-c(z,nm)
  if(length(z)) paste(z,collapse=";") else "None"
},character(1))]
marker_counts <- markers[,.N,by=cluster]
cluster_table[, n_significant_markers:=marker_counts$N[match(cluster,marker_counts$cluster)]]
cluster_table[is.na(n_significant_markers),n_significant_markers:=0L]
top_markers <- markers[order(cluster,-avg_log2FC),.(top20_markers=paste(head(gene,20L),collapse=";")),by=cluster]
cluster_table <- merge(cluster_table,top_markers,by="cluster",all.x=TRUE)
cluster_table[, `:=`(
  annotation_status=data.table::fcase(incompatible_lineage_program,"REVIEW_MIXED_OR_DOUBLET",final_cell_type=="Unresolved","UNRESOLVED",
    patient_enriched,"REVIEW_PATIENT_ENRICHED",canonical_support_n>=4L|high_specific_support_n>=2L,"READY_HIGH_CONFIDENCE",default="READY_BROAD_TYPE_ONLY"),
  annotation_confidence=data.table::fcase(final_cell_type=="Unresolved","Unresolved",canonical_support_n>=4L|high_specific_support_n>=2L,"High",default="Broad_type_only"),
  dataset_id="GSE151214", analysis_role="normal_fallopian_tube_reference",
  notes="Existing post-QC/doublet object; old cluster IDs group cells only and were not hardcoded to labels. Excluded from tumor TME, ligand-receptor meta-analysis, malignant-state inference and CNV analysis."
)]

idx <- match(obj$cluster_v3,cluster_table$cluster)
patient_value <- if("patient_id" %in% names(md)) as.character(md$patient_id) else obj$sample_id_v3
dbl_col <- grep("scDblFinder.class|doublet_call",names(md),value=TRUE)[1L]
doublet_value <- if(!is.na(dbl_col)) as.character(md[[dbl_col]]) else "post_filter_no_call"
assignments <- data.table::data.table(dataset_id="GSE151214",cell_id=colnames(obj),sample_id=obj$sample_id_v3,
  patient_id=patient_value,old_cluster=obj$cluster_v3,final_cluster=obj$cluster_v3,
  final_cell_type=cluster_table$final_cell_type[idx],cell_subtype="",cell_state=cluster_table$cell_state[idx],
  annotation_confidence=cluster_table$annotation_confidence[idx],patient_enriched=cluster_table$patient_enriched[idx],removal_reason=NA_character_)
assignments[, doublet_call:=doublet_value]
if(nrow(assignments)!=expected_cells || anyNA(assignments$final_cell_type)) stop("Assignment coverage failed")
removed <- assignments[0]

data_mat <- SeuratObject::LayerData(obj,assay="RNA",layer="data")
background_genes <- intersect(c("CD44","ITGB1","SPP1","C1QA","C1QB","C1QC","FOLR2","TREM2","EPCAM","PAX8","OVGP1","FOXJ1"),rownames(data_mat))
bg_long <- data.table::rbindlist(lapply(background_genes,function(g) {
  v<-as.numeric(data_mat[g,]); data.table::data.table(final_cell_type=assignments$final_cell_type,sample_id=assignments$sample_id,gene=g,value=v)[,
    .(n_cells=.N,average_expression=mean(value),positive_fraction=mean(value>0)),by=.(final_cell_type,sample_id,gene)]
}))
ci_genes <- intersect(c("CD44","ITGB1"),rownames(obj))
obj <- AddModuleScore(obj,features=list(ci_genes),name="CD44_ITGB1_program",nbin=24,ctrl=25,seed=seed,search=FALSE)
ci <- data.table::data.table(final_cell_type=assignments$final_cell_type,sample_id=assignments$sample_id,score=obj$CD44_ITGB1_program1)[,
  .(CD44_ITGB1_module_mean=mean(score)),by=.(final_cell_type,sample_id)]
background <- merge(bg_long,ci,by=c("final_cell_type","sample_id"),all.x=TRUE)
sp_genes <- intersect(c("SPP1","APOC1","GPNMB","TREM2","LPL","CTSD"),rownames(obj))
obj <- AddModuleScore(obj,features=list(sp_genes),name="SPP1_normal_program",nbin=24,ctrl=25,seed=seed,search=FALSE)
sp_module <- data.table::data.table(final_cell_type=assignments$final_cell_type,sample_id=assignments$sample_id,
  SPP1_module_mean=obj$SPP1_normal_program1)[,.(SPP1_module_mean=mean(SPP1_module_mean)),by=.(final_cell_type,sample_id)]
spp1_background <- background[gene %in% c("SPP1","C1QA","C1QB","C1QC","FOLR2","TREM2") & final_cell_type %in% c("Monocyte","Macrophage","Dendritic")]
spp1_background <- merge(spp1_background,sp_module,by=c("final_cell_type","sample_id"),all.x=TRUE)

keep_cols <- c("dataset_id","analysis_role","cluster","final_cell_type","cell_state","annotation_status","annotation_confidence","n_cells","n_samples","dominant_sample","dominant_sample_fraction","sample_entropy","patient_enriched","canonical_markers","top20_markers","canonical_support_n","high_specific_support_n","second_candidate","incompatible_lineage_program","n_significant_markers","notes")
cluster_table[,cell_subtype:=""]
keep_cols <- append(keep_cols,"cell_subtype",after=4L)
data.table::fwrite(cluster_table[,..keep_cols],file.path(out,"normal_reference_cluster_annotation.csv"),na="NA")
write_gz(assignments,file.path(out,"normal_reference_cell_assignments.csv.gz"))
write_gz(removed,file.path(out,"removed_cells.csv.gz"))
write_gz(markers,file.path(out,"significant_markers.csv.gz"))
data.table::fwrite(dominance,file.path(out,"cluster_sample_dominance.csv"),na="NA")
data.table::fwrite(background,file.path(out,"normal_CD44_ITGB1_background.csv"),na="NA")
data.table::fwrite(spp1_background,file.path(out,"normal_SPP1_myeloid_background.csv"),na="NA")
obj$final_cell_type <- assignments$final_cell_type
p1<-DimPlot(obj,reduction="umap",group.by="final_cell_type",label=TRUE,repel=TRUE,raster=TRUE)+ggplot2::labs(title="GSE151214 normal fallopian-tube reference")
ggplot2::ggsave(file.path(out,"UMAP_normal_cell_types.png"),p1,width=9,height=7,dpi=180)
dot_genes<-intersect(unique(unlist(lapply(panels,head,3L),use.names=FALSE)),rownames(obj))
p2<-DotPlot(obj,features=dot_genes,group.by="final_cell_type")+RotatedAxis()+ggplot2::labs(title="GSE151214 canonical markers")
ggplot2::ggsave(file.path(out,"marker_dotplot.png"),p2,width=16,height=8,dpi=180)

counts<-assignments[,.N,by=final_cell_type][order(final_cell_type)]
summary_lines<-c("# GSE151214 normal fallopian-tube reference","",
  paste0("- input_post_qc_doublet_cells: ",expected_cells),paste0("- retained_cells: ",nrow(assignments)),
  paste0("- mitochondrial_audit: ",if(mt_pass) "PASS_EXISTING_QC_VALID" else "FAIL"),
  paste0("- mitochondrial_features: ",length(mt_genes)),"- mitochondrial_qc_rerun: no",
  "- analysis_role: normal_fallopian_tube_reference",
  "- exclusions: tumor TME meta-analysis; ligand-receptor meta-analysis; malignant-state inference; CNV analysis.",
  "- per-sample artifact review: existing post-QC/doublet-filtered object was retained; no additional cell met a high-confidence doublet, severe-low-quality, or ambient-RNA-only removal rule.",
  "- sample dominance is annotation metadata only; epithelial cells remain uncorrected.","","## Final cell-type counts","",
  paste0("- ",counts$final_cell_type,": ",counts$N),"","CD44/ITGB1 and SPP1/C1QC/FOLR2/TREM2 normal-background summaries include average expression, positive fraction and a CD44/ITGB1 module score.")
writeLines(summary_lines,file.path(out,"analysis_summary.md"),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
message("GSE151214 complete: ",out)
