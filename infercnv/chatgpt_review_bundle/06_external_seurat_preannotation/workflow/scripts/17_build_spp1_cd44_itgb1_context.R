options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(Matrix);library(SeuratObject);library(ggplot2)})
script_arg<-grep("^--file=",commandArgs(),value=TRUE);script_dir<-dirname(normalizePath(sub("^--file=","",script_arg[[1L]])))
source(file.path(script_dir,"_diagnostics_v2_common.R"));z<-read_diagnostics_config();cfg<-z$cfg
data_root<-normalizePath(cfg$project$data_root,winslash="/",mustWork=TRUE);v4<-file.path(data_root,"diagnostics_v4_cross_dataset_validation");v3<-file.path(data_root,"diagnostics_v3_remaining_datasets");cleaned<-file.path(data_root,"diagnostics_v2_marker_ready_cleaned")
replace_generated<-"--replace-generated-output"%in%commandArgs(trailingOnly=TRUE);out<-file.path(v4,"spp1_cd44_itgb1_context_by_sample.csv");if(file.exists(out)&&!replace_generated)stop("Output exists: ",out)
sender<-fread(file.path(v4,"macrophage_state_by_sample.csv"))[dataset_id%in%c("GSE147082","GSE154600","GSE158722")]
receiver_counts<-function(ds,counts,a,c,strict=FALSE){
  idx<-match(as.character(a$final_cluster),as.character(c$cluster));a[,`:=`(annotation_status=c$annotation_status[idx],canonical_support_n=as.numeric(c$canonical_support_n[idx]),incompatible_lineage_program=as.logical(c$incompatible_lineage_program[idx]))]
  allowed<-if(strict)"READY_HIGH_CONFIDENCE"else c("READY_HIGH_CONFIDENCE","READY_BROAD_TYPE_ONLY","REVIEW_PATIENT_ENRICHED")
  a<-a[final_cell_type%in%c("Epithelial","Malignant_epithelial_mesenchymal_like")&annotation_status%in%allowed&canonical_support_n>=3&incompatible_lineage_program!=TRUE]
  cells<-intersect(a$cell_id,colnames(counts));a<-a[match(cells,cell_id)];counts<-counts[,cells,drop=FALSE];lib<-pmax(Matrix::colSums(counts),1)
  g<-intersect(c("CD44","ITGB1"),rownames(counts));x<-counts[g,,drop=FALSE];x@x<-log1p(x@x*rep(1e4/lib,diff(x@p)))
  dt<-data.table(dataset_id=ds,patient_id=as.character(a$patient_id),sample_id=as.character(a$sample_id),cell_id=cells,malignant_epithelial=a$final_cell_type=="Malignant_epithelial_mesenchymal_like",CD44_expression=if("CD44"%in%rownames(x))as.numeric(x["CD44",])else NA_real_,ITGB1_expression=if("ITGB1"%in%rownames(x))as.numeric(x["ITGB1",])else NA_real_,CD44_positive=if("CD44"%in%rownames(counts))as.numeric(counts["CD44",]>0)else NA_real_,ITGB1_positive=if("ITGB1"%in%rownames(counts))as.numeric(counts["ITGB1",]>0)else NA_real_)
  dt[,.(n_high_confidence_epithelial=.N,CD44_average_expression=mean(CD44_expression,na.rm=TRUE),CD44_positive_fraction=mean(CD44_positive,na.rm=TRUE),ITGB1_average_expression=mean(ITGB1_expression,na.rm=TRUE),ITGB1_positive_fraction=mean(ITGB1_positive,na.rm=TRUE),CD44_ITGB1_copositive_fraction=mean(CD44_positive>0&ITGB1_positive>0,na.rm=TRUE),n_malignant_epithelial=sum(malignant_epithelial),malignant_CD44_positive_fraction=if(any(malignant_epithelial))mean(CD44_positive[malignant_epithelial],na.rm=TRUE)else NA_real_,malignant_ITGB1_positive_fraction=if(any(malignant_epithelial))mean(ITGB1_positive[malignant_epithelial],na.rm=TRUE)else NA_real_),by=.(dataset_id,patient_id,sample_id)]
}
receivers<-list()
a<-fread(file.path(v4,"GSE147082_refined","refined_cell_assignments.csv.gz"));c<-fread(file.path(v4,"GSE147082_refined","refined_cluster_annotation.csv"));o<-readRDS(file.path(data_root,"GSE147082","objects","GSE147082_preannotation.rds"));cnt<-SeuratObject::LayerData(o,assay="RNA",layer="counts");receivers[["GSE147082"]]<-receiver_counts("GSE147082",cnt,a,c,FALSE);rm(a,c,o,cnt);gc()
for(ds in c("GSE154600","GSE158722")){a<-fread(file.path(cleaned,ds,"cleaned_cell_assignments.csv.gz"));c<-fread(file.path(cleaned,ds,"cleaned_cluster_annotation_template.csv"));inp<-readRDS(file.path(data_root,"diagnostics_v2","objects",ds,"lineage_inputs","Epithelial_like_strategy_input.rds"));receivers[[ds]]<-receiver_counts(ds,inp$counts,a,c,ds=="GSE158722");rm(a,c,inp);gc()}
recv<-rbindlist(receivers,fill=TRUE);ctx<-merge(sender,recv,by=c("dataset_id","patient_id","sample_id"),all=FALSE)
ctx[,receiver_evaluable:=n_high_confidence_epithelial>=20];ctx[,context_evaluable:=evaluable&receiver_evaluable];ctx[,context_name:="SPP1-associated CD44/ITGB1-positive adhesion context"]
corr<-ctx[context_evaluable==TRUE,{n=.N;if(n>=6){z<-suppressWarnings(cor.test(SPP1_high_macrophage_fraction,CD44_ITGB1_copositive_fraction,method="spearman",exact=FALSE));list(n_evaluable_samples=n,spearman_rho=unname(z$estimate),p_value=z$p.value,analysis_type="Spearman sample-level")}else list(n_evaluable_samples=n,spearman_rho=NA_real_,p_value=NA_real_,analysis_type="descriptive only; n < 6")},by=dataset_id]
ctx<-merge(ctx,corr,by="dataset_id",all.x=TRUE);fwrite(ctx,out,na="NA")
normal<-fread(file.path(v4,"GSE151214_refined","normal_background_cohort_summary.csv"));normal[,dataset_role:="normal_fallopian_tube_reference"];normal[,interpretation:="normal background only; excluded from tumor effect aggregation"]
fwrite(normal,file.path(v4,"normal_background_summary.csv"),na="NA")
p3<-ggplot(ctx[context_evaluable==TRUE],aes(SPP1_high_macrophage_fraction,CD44_ITGB1_copositive_fraction,color=patient_id))+geom_point(size=2)+geom_smooth(method="lm",se=FALSE,color="grey40",linewidth=.6)+facet_wrap(~dataset_id,scales="free")+theme_bw()+labs(title="SPP1-associated CD44/ITGB1-positive adhesion context",x="SPP1-high macrophage fraction",y="epithelial CD44/ITGB1 co-positive fraction")
ggsave(file.path(v4,"03_spp1_epithelial_context_scatter.png"),p3,width=11,height=5,dpi=180)
np<-normal[final_cell_type%in%c("Secretory_epithelial","Ciliated_epithelial","C1QC_macrophage","Macrophage_DC_mixed")&feature%in%c("CD44","ITGB1","SPP1")]
p4<-ggplot(np,aes(final_cell_type,median_positive_fraction,fill=feature))+geom_col(position="dodge")+theme_bw()+theme(axis.text.x=element_text(angle=25,hjust=1))+labs(title="GSE151214 normal-reference positive fractions",x=NULL,y="median across evaluable samples")
ggsave(file.path(v4,"04_normal_reference_background.png"),p4,width=9,height=5,dpi=180)
message("SPP1-CD44/ITGB1 context analysis complete")
