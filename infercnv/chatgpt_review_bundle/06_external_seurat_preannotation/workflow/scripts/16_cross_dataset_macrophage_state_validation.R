options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(Matrix);library(SeuratObject);library(ggplot2)})
script_arg<-grep("^--file=",commandArgs(),value=TRUE); script_dir<-dirname(normalizePath(sub("^--file=","",script_arg[[1L]])))
source(file.path(script_dir,"_diagnostics_v2_common.R")); z<-read_diagnostics_config(); cfg<-z$cfg
data_root<-normalizePath(cfg$project$data_root,winslash="/",mustWork=TRUE); v4<-file.path(data_root,"diagnostics_v4_cross_dataset_validation"); dir.create(v4,recursive=TRUE,showWarnings=FALSE)
replace_generated<-"--replace-generated-output"%in%commandArgs(trailingOnly=TRUE)
outfile<-file.path(v4,"macrophage_state_by_sample.csv"); if(file.exists(outfile)&&!replace_generated)stop("Output exists: ",outfile)
v3<-file.path(data_root,"diagnostics_v3_remaining_datasets"); cleaned<-file.path(data_root,"diagnostics_v2_marker_ready_cleaned")
programs<-list(SPP1=c("SPP1","APOC1","GPNMB","TREM2","LPL","CTSD"),C1QC=c("C1QA","C1QB","C1QC","APOE","MRC1","SELENOP"),FOLR2=c("FOLR2","MRC1","SELENOP","C1QC","LYVE1","CD163"),lipid=c("TREM2","APOC1","GPNMB","LPL","CTSD","LGALS3"),inflammatory=c("S100A8","S100A9","FCN1","VCAN","CTSS","IL1B"))
roles<-c(GSE154600="primary_tumor_ecosystem",GSE158722="malignant_fluid_tumor_ecosystem",GSE147082="tumor_sensitivity_validation",GSE151214="normal_fallopian_tube_reference",GSE154763="author_annotated_myeloid_reference")
rank_pct<-function(x)frank(x,ties.method="average")/sum(is.finite(x))
score_counts<-function(dataset_id,counts,assign,eligible_cells){
  cells<-intersect(eligible_cells,colnames(counts)); assign<-assign[match(cells,cell_id)]; counts<-counts[,cells,drop=FALSE]
  lib<-pmax(Matrix::colSums(counts),1); needed<-intersect(unique(unlist(programs)),rownames(counts)); x<-counts[needed,,drop=FALSE]; x@x<-log1p(x@x*rep(1e4/lib,diff(x@p)))
  ans<-data.table(dataset_id=dataset_id,dataset_role=roles[[dataset_id]],cell_id=cells,patient_id=as.character(assign$patient_id),sample_id=as.character(assign$sample_id),SPP1_expression=if("SPP1"%in%rownames(x))as.numeric(x["SPP1",])else NA_real_,SPP1_positive=if("SPP1"%in%rownames(counts))as.numeric(counts["SPP1",]>0)else NA_real_)
  for(nm in names(programs)){g<-intersect(programs[[nm]],rownames(x)); ans[[paste0(nm,"_raw_score")]]<-Matrix::colMeans(x[g,,drop=FALSE]); ans[[paste0(nm,"_gene_positive_fraction")]]<-Matrix::colMeans(counts[g,,drop=FALSE]>0)}
  for(nm in names(programs)){raw<-paste0(nm,"_raw_score"); ans[[paste0(nm,"_dataset_percentile")]]<-rank_pct(ans[[raw]]); ans[,paste0(nm,"_sample_percentile"):=rank_pct(get(raw)),by=sample_id]}
  ans
}
map_status<-function(a,c,cluster_col="final_cluster"){
  ckey<-as.character(c$cluster); idx<-match(as.character(a[[cluster_col]]),ckey); a[,`:=`(annotation_status=c$annotation_status[idx],canonical_support_n=as.numeric(c$canonical_support_n[idx]),incompatible_lineage_program=as.logical(c$incompatible_lineage_program[idx]))]; a
}
cell_scores<-list()

# GSE147082: targeted refined assignments, original counts only.
a<-fread(file.path(v4,"GSE147082_refined","refined_cell_assignments.csv.gz")); c<-fread(file.path(v4,"GSE147082_refined","refined_cluster_annotation.csv")); a<-map_status(a,c)
o<-readRDS(file.path(data_root,"GSE147082","objects","GSE147082_preannotation.rds")); cnt<-SeuratObject::LayerData(o,assay="RNA",layer="counts")
elig<-a[final_cell_type=="Macrophage" & annotation_status%in%c("READY_HIGH_CONFIDENCE","READY_BROAD_TYPE_ONLY","REVIEW_PATIENT_ENRICHED") & canonical_support_n>=3 & incompatible_lineage_program!=TRUE,cell_id]
cell_scores[["GSE147082"]]<-score_counts("GSE147082",cnt,a,elig); rm(o,cnt,a,c);gc()

# GSE151214: only targeted C1QC macrophage subclusters, never mixed/DC.
a<-fread(file.path(v3,"GSE151214","normal_reference_cell_assignments.csv.gz")); sub<-fread(file.path(v4,"GSE151214_refined","myeloid_cluster18_subclustering.csv")); a[sub,on="cell_id",`:=`(final_cell_type=i.refined_type)]
o<-readRDS(file.path(data_root,"GSE151214","objects","GSE151214_preannotation.rds")); cnt<-SeuratObject::LayerData(o,assay="RNA",layer="counts")
elig<-a[final_cell_type=="C1QC_macrophage",cell_id]; cell_scores[["GSE151214"]]<-score_counts("GSE151214",cnt,a,elig); rm(o,cnt,a,sub);gc()

# GSE154600 and GSE158722 use compact prepared myeloid count inputs.
for(ds in c("GSE154600","GSE158722")){
  a<-fread(file.path(cleaned,ds,"cleaned_cell_assignments.csv.gz")); c<-fread(file.path(cleaned,ds,"cleaned_cluster_annotation_template.csv")); a<-map_status(a,c)
  inp<-readRDS(file.path(data_root,"diagnostics_v2","objects",ds,"lineage_inputs","Myeloid_like_strategy_input.rds")); cnt<-inp$counts
  allowed<-if(ds=="GSE158722")"READY_HIGH_CONFIDENCE" else c("READY_HIGH_CONFIDENCE","READY_BROAD_TYPE_ONLY","REVIEW_PATIENT_ENRICHED")
  elig<-a[final_cell_type=="Macrophage" & annotation_status%in%allowed & canonical_support_n>=3 & incompatible_lineage_program!=TRUE,cell_id]
  cell_scores[[ds]]<-score_counts(ds,cnt,a,elig); rm(a,c,inp,cnt);gc()
}

# GSE154763 uses normalized-expression z scores, converted only to within-reference percentiles.
x<-fread(file.path(v4,"GSE154763_refined","author_annotation_with_refined_state.csv.gz")); x<-x[final_cell_type=="Macrophage"]
s<-data.table(dataset_id="GSE154763",dataset_role=roles[["GSE154763"]],cell_id=x$cell_id_harmonized,patient_id=as.character(x$patient),sample_id=as.character(x$library_id),SPP1_expression=x$SPP1_expression,SPP1_positive=as.numeric(x$SPP1_positive))
mapcols<-c(SPP1="SPP1_program",C1QC="C1QC_program",FOLR2="FOLR2_program",lipid="lipid_associated_program",inflammatory="inflammatory_monocyte_program")
for(nm in names(mapcols)){s[[paste0(nm,"_raw_score")]]<-x[[mapcols[[nm]]]]; s[[paste0(nm,"_gene_positive_fraction")]]<-NA_real_; s[[paste0(nm,"_dataset_percentile")]]<-rank_pct(s[[paste0(nm,"_raw_score")]]); s[,paste0(nm,"_sample_percentile"):=rank_pct(get(paste0(nm,"_raw_score"))),by=sample_id]}
cell_scores[["GSE154763"]]<-s; rm(x,s);gc()
cell_dt<-rbindlist(cell_scores,fill=TRUE); fwrite(cell_dt,file.path(v4,"macrophage_state_cell_scores.csv.gz"),compress="gzip",na="NA")

for(nm in c("SPP1","C1QC","FOLR2"))cell_dt[,paste0(nm,"_high"):=get(paste0(nm,"_dataset_percentile"))>=.75]
sample_state<-cell_dt[,.(n_macrophages=.N,SPP1_average_expression=mean(SPP1_expression,na.rm=TRUE),SPP1_positive_fraction=mean(SPP1_positive,na.rm=TRUE),SPP1_program_median_percentile=median(SPP1_dataset_percentile,na.rm=TRUE),C1QC_program_median_percentile=median(C1QC_dataset_percentile,na.rm=TRUE),FOLR2_program_median_percentile=median(FOLR2_dataset_percentile,na.rm=TRUE),lipid_program_median_percentile=median(lipid_dataset_percentile,na.rm=TRUE),SPP1_high_macrophage_fraction=mean(SPP1_high),C1QC_high_macrophage_fraction=mean(C1QC_high),FOLR2_high_macrophage_fraction=mean(FOLR2_high)),by=.(dataset_id,dataset_role,patient_id,sample_id)]
sample_state[,evaluable:=n_macrophages>=20]; fwrite(sample_state,outfile,na="NA")

patient_state<-cell_dt[,.(n_macrophages=.N,SPP1_positive_fraction=mean(SPP1_positive,na.rm=TRUE),SPP1_median=median(SPP1_dataset_percentile),C1QC_median=median(C1QC_dataset_percentile),FOLR2_median=median(FOLR2_dataset_percentile)),by=.(dataset_id,dataset_role,patient_id)]
patient_state[,evaluable:=n_macrophages>=20]
repro<-rbindlist(lapply(c("SPP1","C1QC","FOLR2"),function(st)patient_state[,{
  med<-get(paste0(st,"_median")); pos<-if(st=="SPP1")SPP1_positive_fraction else med
  ev<-evaluable; positive<-ev & med>=.60 & pos>=.10; nep<-uniqueN(patient_id[ev]); npos<-uniqueN(patient_id[positive])
  status<-if(nep==0)"NOT_EVALUABLE" else if(npos>=2)"REPLICATED" else if(nep==1&&npos==1)"SINGLE_PATIENT_ONLY" else if(npos==1)"SUPPORTIVE" else "NOT_SUPPORTED"
  list(n_high_confidence_macrophages=sum(n_macrophages),n_evaluable_patients=nep,n_positive_patients=npos,positive_patients=paste(patient_id[positive],collapse=";"),state_status=status)
},by=.(dataset_id,dataset_role)][,state:=st]))
setcolorder(repro,c("dataset_id","dataset_role","state")); fwrite(repro,file.path(v4,"macrophage_state_reproducibility_matrix.csv"),na="NA")

heat<-melt(sample_state,id.vars=c("dataset_id","sample_id","evaluable"),measure.vars=c("SPP1_program_median_percentile","C1QC_program_median_percentile","FOLR2_program_median_percentile"),variable.name="state",value.name="median_percentile")[evaluable==TRUE,.(median_percentile=median(median_percentile,na.rm=TRUE)),by=.(dataset_id,state)]
p1<-ggplot(heat,aes(state,dataset_id,fill=median_percentile))+geom_tile(color="white")+scale_fill_viridis_c(limits=c(0,1))+theme_bw()+labs(title="Macrophage state percentiles by dataset",x=NULL,y=NULL)
ggsave(file.path(v4,"01_macrophage_state_heatmap.png"),p1,width=8,height=4.5,dpi=180)
p2<-ggplot(sample_state[dataset_id%in%c("GSE154600","GSE158722","GSE147082")&evaluable==TRUE],aes(dataset_id,SPP1_high_macrophage_fraction,color=patient_id))+geom_jitter(width=.15,height=0,size=2)+theme_bw()+labs(title="Tumor-dataset SPP1-high macrophage prevalence",x=NULL,y="fraction")
ggsave(file.path(v4,"02_tumor_spp1_macrophage_fraction.png"),p2,width=8,height=5,dpi=180)
message("Cross-dataset macrophage state validation complete")
