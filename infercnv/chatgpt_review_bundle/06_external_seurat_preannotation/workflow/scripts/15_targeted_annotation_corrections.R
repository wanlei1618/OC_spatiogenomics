options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(Seurat); library(ggplot2)})
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "_diagnostics_v2_common.R"))
z <- read_diagnostics_config(); cfg <- z$cfg
data_root <- normalizePath(cfg$project$data_root, winslash="/", mustWork=TRUE)
v3 <- file.path(data_root,"diagnostics_v3_remaining_datasets")
out_root <- file.path(data_root,"diagnostics_v4_cross_dataset_validation")
dir.create(out_root,recursive=TRUE,showWarnings=FALSE)
seed <- as.integer(cfg$project$random_seed %||% 20260718L); set.seed(seed)
replace_generated <- "--replace-generated-output" %in% commandArgs(trailingOnly=TRUE)
write_gz <- function(x,p) fwrite(x,p,na="NA",compress="gzip")
norm_target <- function(counts, genes) {
  genes <- intersect(genes,rownames(counts)); lib <- pmax(Matrix::colSums(counts),1)
  x <- counts[genes,,drop=FALSE]; x@x <- log1p(x@x * rep.int(1e4/lib,diff(x@p)))
  x
}
row_module <- function(mat, genes) {
  use <- intersect(genes,rownames(mat)); if(!length(use)) return(rep(NA_real_,ncol(mat)))
  Matrix::colMeans(mat[use,,drop=FALSE])
}
positive_fraction <- function(mat, genes) {
  use <- intersect(genes,rownames(mat)); if(!length(use)) return(rep(0,ncol(mat)))
  Matrix::colMeans(mat[use,,drop=FALSE] > 0)
}

## GSE147082 targeted review -------------------------------------------------
o147 <- file.path(out_root,"GSE147082_refined")
if(dir.exists(o147) && !replace_generated) stop("Output exists: ",o147)
if(!dir.exists(o147)) dir.create(o147,recursive=TRUE)
a147 <- fread(file.path(v3,"GSE147082","cleaned_cell_assignments.csv.gz"))
c147 <- fread(file.path(v3,"GSE147082","cleaned_cluster_annotation_template.csv"))
m147 <- fread(file.path(v3,"GSE147082","significant_markers.csv.gz"))
obj147 <- readRDS(file.path(data_root,"GSE147082","objects","GSE147082_preannotation.rds"))
counts147 <- SeuratObject::LayerData(obj147,assay="RNA",layer="counts")[,a147$cell_id,drop=FALSE]
all_cnv_genes <- rownames(counts147)
tmp_ref <- "D:/TEMP/ucsc_hg38_refGene.txt.gz"
if(!file.exists(tmp_ref)) utils::download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",tmp_ref,mode="wb",quiet=TRUE,timeout=300)
refgene <- fread(tmp_ref,header=FALSE,select=c(3,5,6,13),showProgress=FALSE)
setnames(refgene,c("chrom","start","end","gene")); refgene <- refgene[chrom %chin% paste0("chr",c(1:22,"X"))]
refgene <- refgene[,.(chrom=chrom[1L],position=as.integer(median((start+end)/2))),by=gene]
gene_order <- merge(data.table(gene=all_cnv_genes,gene_upper=toupper(all_cnv_genes)),
                    refgene[,.(gene_upper=toupper(gene),chrom,position)],by="gene_upper")
gene_order[,chrom_order:=match(chrom,paste0("chr",c(1:22,"X")))]
setorder(gene_order,chrom_order,position)
obj147s <- CreateSeuratObject(counts147,min.cells=0,min.features=0)
obj147s <- NormalizeData(obj147s,verbose=FALSE)
d147 <- SeuratObject::LayerData(obj147s,assay="RNA",layer="data")
ref_types <- c("T_cell","B_cell","NK_cell","Macrophage","Monocyte","pDC","cDC1","cDC2","Endothelial","Fibroblast","Pericyte")
ref_cells <- a147[final_cell_type %in% ref_types & annotation_confidence %in% c("High","Broad_type_only"),cell_id]
if(length(ref_cells)<100) stop("Insufficient normal reference cells for GSE147082 CNV-equivalent score")
ref_by_type <- lapply(intersect(ref_types,unique(a147$final_cell_type)),function(tp) {
  zc <- intersect(a147[final_cell_type==tp & annotation_confidence %in% c("High","Broad_type_only"),cell_id],colnames(d147))
  if(!length(zc)) return(NULL); Matrix::rowMeans(d147[,zc,drop=FALSE])
}); ref_by_type <- Filter(Negate(is.null),ref_by_type)
ref_mean <- Reduce("+",ref_by_type)/length(ref_by_type)
group_cells <- list(
  cluster_4=a147[final_cluster==4,cell_id], cluster_7=a147[final_cluster==7,cell_id],
  tumor_epithelial=a147[final_cell_type=="Epithelial" & annotation_confidence %in% c("High","Broad_type_only","Review"),cell_id],
  normal_fibro_pericyte=a147[final_cell_type %in% c("Fibroblast","Pericyte") & annotation_confidence %in% c("High","Broad_type_only","Review"),cell_id]
)
profiles <- lapply(group_cells,function(cells) Matrix::rowMeans(d147[,intersect(cells,colnames(d147)),drop=FALSE])-ref_mean)
smooth_one <- function(v) {
  dt <- copy(gene_order); dt[,value:=v[match(gene,rownames(d147))]]; dt <- dt[is.finite(value)]
  dt[,smooth:=if(.N>=51) stats::runmed(value,k=min(101L,if(.N%%2L) .N else .N-1L),endrule="median") else value,by=chrom]
  dt
}
smoothed <- lapply(profiles,smooth_one)
epi_profile <- smoothed$tumor_epithelial$smooth
str_profile <- smoothed$normal_fibro_pericyte$smooth
cnv_rows <- rbindlist(lapply(names(smoothed),function(nm){x<-smoothed[[nm]]; data.table(
  comparison_group=nm,n_cells=length(group_cells[[nm]]),n_ordered_genes=nrow(x),
  cnv_intensity=mean(abs(x$smooth),na.rm=TRUE),
  correlation_to_tumor_epithelial=cor(x$smooth,epi_profile,use="complete.obs"),
  correlation_to_normal_stroma=cor(x$smooth,str_profile,use="complete.obs"))}))
epi_genes <- c("EPCAM","KRT7","KRT8","KRT18","KRT19","PAX8","WFDC2","MSLN","MUC1")
stroma_genes <- c("COL1A1","COL1A2","COL3A1","DCN","LUM","FBN1","FBN2","COL6A1","COL6A2","PDGFRA")
cart_genes <- c("COL2A1","CHAD","MATN3","COL9A1","HAPLN1","ACAN","COL11A2","MATN1","MATN4","ITGA10")
mal_mes_genes <- c("IL13RA2","SIX1","LGR6","CLU","FBN2","MATN2","VIM","SNAI2","ZEB1")
module_summary <- rbindlist(lapply(c("4","7"),function(cl){cells<-a147[final_cluster==cl,cell_id]; data.table(
  cluster=cl,n_cells=length(cells),epithelial_module=mean(row_module(d147[,cells,drop=FALSE],epi_genes)),
  stromal_module=mean(row_module(d147[,cells,drop=FALSE],stroma_genes)),
  cartilage_module=mean(row_module(d147[,cells,drop=FALSE],cart_genes)),
  malignant_mesenchymal_module=mean(row_module(d147[,cells,drop=FALSE],mal_mes_genes)),
  cartilage_marker_support=uniqueN(m147[cluster==cl & toupper(gene)%in%cart_genes,toupper(gene)]),
  epithelial_marker_support=uniqueN(m147[cluster==cl & toupper(gene)%in%epi_genes,toupper(gene)]),
  principal_markers=paste(head(m147[cluster==cl][order(-avg_log2FC),gene],20),collapse=";"))}))
target_cnv <- cnv_rows[comparison_group %chin% c("cluster_4","cluster_7")]
target_cnv[,cluster:=sub("cluster_","",comparison_group)]
target_cnv <- merge(target_cnv,module_summary,by="cluster")
target_cnv[,n_cells:=n_cells.x]
target_cnv[,c("n_cells.x","n_cells.y"):=NULL]
epi_int <- cnv_rows[comparison_group=="tumor_epithelial",cnv_intensity]
str_int <- cnv_rows[comparison_group=="normal_fibro_pericyte",cnv_intensity]
mid <- mean(c(epi_int,str_int))
target_cnv[,cnv_normal_like:=correlation_to_normal_stroma>=correlation_to_tumor_epithelial]
target_cnv[,final_judgment:=fcase(
  !cnv_normal_like & cnv_intensity>=mid & epithelial_module>=stromal_module,"Malignant_epithelial_mesenchymal_like",
  cnv_normal_like & cartilage_marker_support>=5,"COL2A1_positive_chondrocyte_like_fibroblast",
  cnv_normal_like & stromal_module>epithelial_module,"Mesenchymal_stromal",
  default="Unresolved_mesenchymal")]
fwrite(target_cnv,file.path(o147,"target_cluster_cnv_summary.csv"),na="NA")

# Cluster 6: targeted low-resolution reclustering and lineage programs.
c6cells <- a147[final_cluster==6,cell_id]
c6 <- CreateSeuratObject(counts147[,c6cells,drop=FALSE],min.cells=0,min.features=0)
c6 <- NormalizeData(c6,verbose=FALSE); c6 <- FindVariableFeatures(c6,nfeatures=2000,verbose=FALSE)
c6 <- ScaleData(c6,features=VariableFeatures(c6),verbose=FALSE); c6 <- RunPCA(c6,npcs=30,verbose=FALSE)
c6 <- FindNeighbors(c6,dims=1:20,verbose=FALSE); c6 <- RunUMAP(c6,dims=1:20,seed.use=seed,verbose=FALSE)
c6 <- FindClusters(c6,resolution=.25,random.seed=seed,verbose=FALSE)
d6 <- SeuratObject::LayerData(c6,assay="RNA",layer="data")
tgenes<-c("CD3D","CD3E","TRBC1","TRBC2","CD2"); gdgenes<-c("TRDC","TRGC1","TRGC2"); nkgenes<-c("NKG7","GNLY","KLRD1","XCL1","XCL2","PRF1")
c6dt <- data.table(cell_id=colnames(c6),subcluster=as.character(Idents(c6)),T_score=row_module(d6,tgenes),gamma_delta_score=row_module(d6,gdgenes),NK_score=row_module(d6,nkgenes),CD3_positive_fraction=positive_fraction(d6,c("CD3D","CD3E")),TRG_positive_fraction=positive_fraction(d6,gdgenes),NK_positive_fraction=positive_fraction(d6,nkgenes))
c6sum <- c6dt[,.(n_cells=.N,T_score=mean(T_score),gamma_delta_score=mean(gamma_delta_score),NK_score=mean(NK_score),CD3_positive_fraction=mean(CD3_positive_fraction),TRG_positive_fraction=mean(TRG_positive_fraction),NK_positive_fraction=mean(NK_positive_fraction)),by=subcluster]
c6sum[,refined_type:=fcase(
  TRG_positive_fraction>=.25 & gamma_delta_score>=.35,"Gamma_delta_T",
  CD3_positive_fraction<.25 & NK_positive_fraction>=.35,"NK_cell",
  CD3_positive_fraction>=.4 & T_score>=.3 & TRG_positive_fraction<.2 & NK_score<T_score*1.5,"Cytotoxic_T",
  CD3_positive_fraction>=.3 & NK_positive_fraction>=.3,"NKT_like_unresolved_cytotoxic",
  default="Unresolved_cytotoxic_lymphocyte")]
c6dt[,refined_type:=c6sum$refined_type[match(subcluster,c6sum$subcluster)]]

review <- rbindlist(list(
  target_cnv[,.(target=paste0("cluster_",cluster),n_cells,decision=final_judgment,cnv_intensity,epithelial_module,stromal_module,cartilage_module,principal_markers)],
  c6sum[,.(target=paste0("cluster_6_subcluster_",subcluster),n_cells,decision=refined_type,cnv_intensity=NA_real_,epithelial_module=NA_real_,stromal_module=NA_real_,cartilage_module=NA_real_,principal_markers=paste(c(tgenes,gdgenes,nkgenes),collapse=";"))],
  data.table(target="cluster_19",n_cells=a147[final_cluster==19,.N],decision="Unresolved",cnv_intensity=NA_real_,epithelial_module=NA_real_,stromal_module=NA_real_,cartilage_module=NA_real_,principal_markers=paste(head(m147[cluster==19][order(-avg_log2FC),gene],20),collapse=";"))
),fill=TRUE)
fwrite(review,file.path(o147,"target_cluster_review.csv"),na="NA")
refined147 <- copy(a147)
refined147[,cell_subtype:=as.character(cell_subtype)]
for(cl in c("4","7")) refined147[final_cluster==cl,`:=`(final_cell_type=target_cnv[cluster==cl,final_judgment],cell_subtype=target_cnv[cluster==cl,final_judgment],annotation_confidence=ifelse(grepl("Unresolved",target_cnv[cluster==cl,final_judgment]),"Unresolved","Review"))]
refined147[final_cluster==6,`:=`(final_cell_type=c6dt$refined_type[match(cell_id,c6dt$cell_id)],cell_subtype=c6dt$refined_type[match(cell_id,c6dt$cell_id)],annotation_confidence="Review")]
refined147[final_cluster==19,`:=`(final_cell_type="Unresolved",cell_subtype="",annotation_confidence="Unresolved")]
write_gz(refined147,file.path(o147,"refined_cell_assignments.csv.gz"))
refined_cluster <- copy(c147)
refined_cluster[,cell_subtype:=as.character(cell_subtype)]
for(cl in c("4","7")) refined_cluster[cluster==cl,`:=`(final_cell_type=target_cnv[cluster==cl,final_judgment],cell_subtype=target_cnv[cluster==cl,final_judgment],annotation_status=ifelse(grepl("Unresolved",target_cnv[cluster==cl,final_judgment]),"UNRESOLVED","REVIEW_PATIENT_ENRICHED"),annotation_confidence="Review",notes=paste0(notes," Targeted CNV-equivalent and lineage-module review."))]
refined_cluster[cluster==6,`:=`(final_cell_type="Cytotoxic_lymphocyte_subclustered",cell_subtype=paste(unique(c6sum$refined_type),collapse=";"),annotation_status="REVIEW_AMBIGUOUS",annotation_confidence="Review",notes=paste0(notes," Cluster 6 was subclustered; use refined cell assignments."))]
fwrite(refined_cluster,file.path(o147,"refined_cluster_annotation.csv"),na="NA")
emb <- as.data.table(Embeddings(obj147,"umap"),keep.rownames="cell_id"); emb_dims<-setdiff(names(emb),"cell_id")[1:2]; setnames(emb,emb_dims,c("dim1","dim2")); emb <- merge(emb,refined147[,.(cell_id,final_cluster,final_cell_type)],by="cell_id")
pt <- emb[final_cluster %in% c(4,6,7,19)]
p <- ggplot(pt,aes(dim1,dim2,color=final_cell_type))+geom_point(size=.7,alpha=.8)+theme_bw()+labs(title="GSE147082 targeted clusters on existing embedding",color="refined type")
ggsave(file.path(o147,"UMAP_refined_target_clusters.png"),p,width=9,height=7,dpi=180)
writeLines(c("# GSE147082 targeted refinement","",paste0("- cluster 4: ",target_cnv[cluster=="4",final_judgment]),paste0("- cluster 7: ",target_cnv[cluster=="7",final_judgment]),paste0("- cluster 6 subtypes: ",paste(c6sum$refined_type,c6sum$n_cells,sep="=",collapse="; ")),"- cluster 19: Unresolved and excluded from downstream proportion/state/LR analyses.","- CNV-equivalent method: chromosome-ordered 101-gene running-median deviation from balanced high-confidence normal lineages using UCSC hg38 refGene order.","- Patient enrichment was not a removal criterion."),file.path(o147,"analysis_summary.md"))
rm(obj147,obj147s,c6,d147,d6,counts147); gc()

## GSE151214 cluster 18 and normal reference -------------------------------
o151 <- file.path(out_root,"GSE151214_refined")
if(dir.exists(o151) && !replace_generated) stop("Output exists: ",o151)
if(!dir.exists(o151)) dir.create(o151,recursive=TRUE)
a151 <- fread(file.path(v3,"GSE151214","normal_reference_cell_assignments.csv.gz"))
c151 <- fread(file.path(v3,"GSE151214","normal_reference_cluster_annotation.csv"))
obj151 <- readRDS(file.path(data_root,"GSE151214","objects","GSE151214_preannotation.rds"))
counts151 <- SeuratObject::LayerData(obj151,assay="RNA",layer="counts")[,a151$cell_id,drop=FALSE]
c18cells <- a151[final_cluster==18,cell_id]
s18 <- CreateSeuratObject(counts151[,c18cells,drop=FALSE],min.cells=0,min.features=0)
s18 <- NormalizeData(s18,verbose=FALSE); s18 <- FindVariableFeatures(s18,nfeatures=2000,verbose=FALSE); s18 <- ScaleData(s18,features=VariableFeatures(s18),verbose=FALSE); s18 <- RunPCA(s18,npcs=25,verbose=FALSE); s18 <- FindNeighbors(s18,dims=1:15,verbose=FALSE); s18 <- FindClusters(s18,resolution=.2,random.seed=seed,verbose=FALSE)
d18 <- SeuratObject::LayerData(s18,assay="RNA",layer="data")
macro_genes<-c("C1QA","C1QB","C1QC","APOE","MRC1","CD68","LST1"); dc_genes<-c("CD1C","FCER1A","CLEC10A","CD1E","HLA-DRA")
s18dt<-data.table(cell_id=colnames(s18),subcluster=as.character(Idents(s18)),macrophage_score=row_module(d18,macro_genes),dc_score=row_module(d18,dc_genes),macrophage_positive_fraction=positive_fraction(d18,macro_genes),dc_positive_fraction=positive_fraction(d18,dc_genes))
s18sum<-s18dt[,.(n_cells=.N,macrophage_score=mean(macrophage_score),dc_score=mean(dc_score),macrophage_positive_fraction=mean(macrophage_positive_fraction),dc_positive_fraction=mean(dc_positive_fraction)),by=subcluster]
s18sum[,refined_type:=fcase(macrophage_positive_fraction>=.35 & macrophage_score>dc_score*1.15,"C1QC_macrophage",dc_positive_fraction>=.3 & dc_score>macrophage_score*1.15,"CD1C_CLEC10A_dendritic_like",default="Macrophage_DC_mixed")]
s18dt[,refined_type:=s18sum$refined_type[match(subcluster,s18sum$subcluster)]]
fwrite(merge(s18dt,s18sum[,.(subcluster,n_cells,refined_type)],by=c("subcluster","refined_type")),file.path(o151,"myeloid_cluster18_subclustering.csv"),na="NA")
a151ref<-copy(a151); a151ref[,cell_subtype:=as.character(cell_subtype)]; a151ref[final_cluster==18,`:=`(final_cell_type=s18dt$refined_type[match(cell_id,s18dt$cell_id)],cell_subtype=s18dt$refined_type[match(cell_id,s18dt$cell_id)],annotation_confidence="Review")]
c151[annotation_status=="REVIEW_PATIENT_ENRICHED",annotation_confidence:="Review"]
c18rows<-s18sum[,.(dataset_id="GSE151214",analysis_role="normal_fallopian_tube_reference",cluster=paste0("18_sub",subcluster),final_cell_type=refined_type,cell_subtype=refined_type,cell_state="None",annotation_status="REVIEW_AMBIGUOUS",annotation_confidence="Review",n_cells,n_samples=NA_integer_,dominant_sample=NA_character_,dominant_sample_fraction=NA_real_,sample_entropy=NA_real_,patient_enriched=FALSE,canonical_markers=paste(c(macro_genes,dc_genes),collapse=";"),top20_markers="targeted subclustering",canonical_support_n=NA_integer_,high_specific_support_n=NA_integer_,second_candidate=NA_character_,incompatible_lineage_program=refined_type=="Macrophage_DC_mixed",n_significant_markers=NA_integer_,notes="Targeted low-resolution subclustering of original cluster 18.")]
c151ref<-rbindlist(list(c151[cluster!=18],c18rows),fill=TRUE)
fwrite(c151ref,file.path(o151,"normal_reference_annotation_refined.csv"),na="NA")

features <- c("CD44","ITGB1","SPP1","C1QA","C1QB","C1QC","FOLR2","TREM2")
nx <- norm_target(counts151,features); lib151<-pmax(Matrix::colSums(counts151),1)
groups<-split(seq_len(nrow(a151ref)),paste(a151ref$sample_id,a151ref$final_cell_type,sep="||"))
bg<-rbindlist(lapply(groups,function(ii){cells<-a151ref$cell_id[ii]; use<-match(cells,colnames(nx)); use<-use[!is.na(use)]; if(!length(use)) return(NULL); rawidx<-match(cells,colnames(counts151)); rawidx<-rawidx[!is.na(rawidx)]; rbindlist(lapply(features,function(g){v<-if(g%in%rownames(nx)) as.numeric(nx[g,use]) else rep(NA_real_,length(use)); rawsum<-if(g%in%rownames(counts151)) sum(counts151[g,rawidx]) else NA_real_; data.table(dataset_id="GSE151214",sample_id=a151ref$sample_id[ii[1]],final_cell_type=a151ref$final_cell_type[ii[1]],feature=g,feature_type="gene",n_cells=length(use),average_expression=mean(v,na.rm=TRUE),positive_fraction=mean(v>0,na.rm=TRUE),pseudobulk_expression=log1p(1e6*rawsum/sum(lib151[rawidx]))) }))}))
bg[,evaluable:=n_cells>=20]
bg[,within_dataset_percentile:=ifelse(evaluable,frank(average_expression,ties.method="average")/sum(evaluable),NA_real_),by=feature]
excluded<-bg[evaluable==FALSE]; fwrite(excluded,file.path(o151,"excluded_small_groups.csv"),na="NA")
cohort<-bg[evaluable==TRUE,.(median_across_samples=median(average_expression,na.rm=TRUE),IQR_across_samples=IQR(average_expression,na.rm=TRUE),median_positive_fraction=median(positive_fraction,na.rm=TRUE),n_evaluable_samples=.N),by=.(dataset_id,final_cell_type,feature,feature_type)]
bg<-merge(bg,cohort[,.(final_cell_type,feature,median_across_samples,IQR_across_samples,n_evaluable_samples)],by=c("final_cell_type","feature"),all.x=TRUE)
fwrite(bg,file.path(o151,"normal_background_by_sample.csv"),na="NA"); fwrite(cohort,file.path(o151,"normal_background_cohort_summary.csv"),na="NA")
writeLines(c("# GSE151214 targeted normal-reference correction","",paste0("- cluster 18 subclusters: ",paste(s18sum$refined_type,s18sum$n_cells,sep="=",collapse="; ")),"- Original 48,837-cell QC and clustering were not rerun.","- Groups with fewer than 20 cells are retained in excluded_small_groups.csv and excluded from quantitative baseline summaries.","- REVIEW_PATIENT_ENRICHED rows now have annotation_confidence=Review.","- This dataset remains a normal fallopian-tube reference and is excluded from tumor effect aggregation."),file.path(o151,"analysis_summary.md"))
rm(obj151,s18,d18,counts151,nx); gc()

## GSE154763 state thresholds ------------------------------------------------
o154 <- file.path(out_root,"GSE154763_refined")
if(dir.exists(o154) && !replace_generated) stop("Output exists: ",o154)
if(!dir.exists(o154)) dir.create(o154,recursive=TRUE)
a154<-fread(file.path(v3,"GSE154763","author_annotation_harmonized.csv.gz"))
s154<-fread(file.path(v3,"GSE154763","myeloid_state_scores.csv.gz"))
key<-"cell_id_harmonized"; setkeyv(a154,key); setkeyv(s154,key); x154<-merge(a154,s154,by=key,suffixes=c("","_score"))
programs<-c("SPP1_program","C1QC_program","FOLR2_program","lipid_associated_program","GPNMB_hypoxia_program","inflammatory_monocyte_program")
lineage_programs<-list(Macrophage=programs[1:5],Monocyte=c("inflammatory_monocyte_program","SPP1_program"))
audit<-rbindlist(lapply(names(lineage_programs),function(lin) rbindlist(lapply(lineage_programs[[lin]],function(pg)data.table(lineage=lin,program=pg,q75=quantile(x154[final_cell_type==lin,get(pg)],.75,na.rm=TRUE),threshold_rule="score > max(lineage_q75,0); top-second >= 0.15")))))
x154[,`:=`(module_state="None",module_state_top_score=NA_real_,module_state_second_score=NA_real_,module_state_margin=NA_real_)]
for(lin in names(lineage_programs)){
  ii<-which(x154$final_cell_type==lin); pgs<-lineage_programs[[lin]]; zmat<-as.matrix(x154[ii,..pgs]); ord<-t(apply(zmat,1,order,decreasing=TRUE)); top<-zmat[cbind(seq_len(nrow(zmat)),ord[,1])]; second<-zmat[cbind(seq_len(nrow(zmat)),ord[,2])]; topname<-pgs[ord[,1]]; th<-audit$q75[match(paste(lin,topname),paste(audit$lineage,audit$program))]; state<-ifelse(top<=pmax(th,0),"None",ifelse(top-second<.15,"Mixed",sub("_program$","",topname))); x154$module_state[ii]<-state; x154$module_state_top_score[ii]<-top; x154$module_state_second_score[ii]<-second; x154$module_state_margin[ii]<-top-second
}
x154[!final_cell_type%in%c("Macrophage","Monocyte"),module_state:="None"]
fwrite(audit,file.path(o154,"state_threshold_audit.csv"),na="NA")
write_gz(x154,file.path(o154,"author_annotation_with_refined_state.csv.gz"))
sample_col<-if("library_id"%in%names(x154))"library_id" else "sample_id"; patient_col<-if("patient"%in%names(x154))"patient" else "patient_id"
byps<-x154[,c(list(n_cells=.N,subtype_fraction=.N/.N[1L],SPP1_positive_fraction=mean(SPP1_positive,na.rm=TRUE)),lapply(.SD,median,na.rm=TRUE)),by=c(patient_col,sample_col,"cell_type_original","final_cell_type"),.SDcols=programs]
# Correct subtype fractions within patient/sample.
byps[,subtype_fraction:=n_cells/sum(n_cells),by=c(patient_col,sample_col)]
fwrite(byps,file.path(o154,"state_by_patient_sample.csv"),na="NA")
expected<-fcase(grepl("Macro_SPP1",x154$cell_type_original),"SPP1",grepl("Macro_C1QC",x154$cell_type_original),"C1QC",grepl("Mono_CD14",x154$cell_type_original),"inflammatory_monocyte",default=NA_character_)
x154[,expected_state:=expected]
conc<-x154[!is.na(expected_state),.(n_cells=.N,n_concordant=sum(module_state==expected_state),concordance_fraction=mean(module_state==expected_state),none_fraction=mean(module_state=="None"),mixed_fraction=mean(module_state=="Mixed")),by=cell_type_original]
fwrite(conc,file.path(o154,"author_subtype_state_concordance.csv"),na="NA")
pair<-x154[get(patient_col)=="P20190304",.(n_cells=.N,SPP1_positive_fraction=mean(SPP1_positive,na.rm=TRUE),SPP1_program_median=median(SPP1_program),C1QC_program_median=median(C1QC_program),FOLR2_program_median=median(FOLR2_program)),by=c(patient_col,sample_col,"tissue")]
pair[,interpretation:="descriptive paired comparison; no statistical inference"]
fwrite(pair,file.path(o154,"paired_tumor_normal_description.csv"),na="NA")
writeLines(c("# GSE154763 refined module-state validation","",paste0("- cells: ",nrow(x154)),paste0("- None: ",sum(x154$module_state=="None")),paste0("- Mixed: ",sum(x154$module_state=="Mixed")),paste0("- assigned lineage-restricted states: ",sum(!x154$module_state%in%c("None","Mixed"))),"- Author subtype remains the primary annotation; module_state is validation only.","- DC and Mast cells are never forced into macrophage states.","- Threshold: lineage-specific 75th percentile, positive score, and top-minus-second >= 0.15."),file.path(o154,"analysis_summary.md"))
message("Targeted corrections complete: ",out_root)
