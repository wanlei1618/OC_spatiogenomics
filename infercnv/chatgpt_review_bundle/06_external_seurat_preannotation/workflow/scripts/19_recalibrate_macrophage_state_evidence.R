options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
script_arg<-grep("^--file=",commandArgs(),value=TRUE);script_dir<-dirname(normalizePath(sub("^--file=","",script_arg[[1L]])))
source(file.path(script_dir,"_diagnostics_v2_common.R"));z<-read_diagnostics_config();cfg<-z$cfg
data_root<-normalizePath(cfg$project$data_root,winslash="/",mustWork=TRUE);v4<-file.path(data_root,"diagnostics_v4_cross_dataset_validation");v5<-file.path(data_root,"diagnostics_v5_final_calibration");dir.create(v5,recursive=TRUE,showWarnings=FALSE)
replace_generated<-"--replace-generated-output"%in%commandArgs(trailingOnly=TRUE);out<-file.path(v5,"macrophage_spp1_presence_by_sample.csv");if(file.exists(out)&&!replace_generated)stop("Output exists: ",out)
s<-fread(file.path(v4,"macrophage_state_by_sample.csv"))
s[,SPP1_presence:=fcase(n_macrophages<20,"NOT_EVALUABLE_FOR_REPLICATION",SPP1_positive_fraction>=.10&SPP1_high_macrophage_fraction>=.10,"PRESENT",default="ABSENT_OR_LOW")]
s[,SPP1_relative_enrichment:=fcase(n_macrophages<20,"NOT_EVALUABLE",SPP1_program_median_percentile>=.60,"RELATIVELY_ENRICHED",default="NOT_RELATIVELY_ENRICHED")]
s[,evaluable:=n_macrophages>=20]
spp1sample<-s[,.(dataset_id,dataset_role,patient_id,sample_id,n_macrophages,SPP1_average_expression,SPP1_positive_fraction,SPP1_high_macrophage_fraction,SPP1_program_median_percentile,SPP1_presence,SPP1_relative_enrichment,evaluable)]
fwrite(spp1sample,out,na="NA")
patient<-spp1sample[,.(any_sample_present=any(SPP1_presence=="PRESENT"),n_evaluable_samples=sum(evaluable),n_present_samples=sum(SPP1_presence=="PRESENT"),median_SPP1_positive_fraction=median(SPP1_positive_fraction[evaluable],na.rm=TRUE),median_SPP1_high_fraction=median(SPP1_high_macrophage_fraction[evaluable],na.rm=TRUE),any_relative_enrichment=any(SPP1_relative_enrichment=="RELATIVELY_ENRICHED"),evaluable_patient=any(evaluable),sample_timepoint_summary=paste(paste0(sample_id,":",SPP1_presence),collapse=";")),by=.(dataset_id,dataset_role,patient_id)]
fwrite(patient,file.path(v5,"macrophage_spp1_presence_by_patient.csv"),na="NA")
repro<-patient[,{
  nep=sum(evaluable_patient);np=sum(any_sample_present&evaluable_patient)
  status<-if(nep==0)"NOT_EVALUABLE"else if(np>=2)"REPLICATED"else if(np==1)"SUPPORTIVE_SINGLE_PATIENT"else if(nep>=2)"NOT_REPLICATED"else"NOT_EVALUABLE"
  list(n_evaluable_patients=nep,n_present_patients=np,present_patients=paste(patient_id[any_sample_present&evaluable_patient],collapse=";"),SPP1_reproducibility_status=status)
},by=.(dataset_id,dataset_role)]
fwrite(repro,file.path(v5,"macrophage_spp1_reproducibility_by_dataset.csv"),na="NA")
enrich<-spp1sample[,.(dataset_id,dataset_role,patient_id,sample_id,n_macrophages,evaluable,SPP1_program_median_percentile,SPP1_relative_enrichment)]
fwrite(enrich,file.path(v5,"macrophage_spp1_relative_enrichment.csv"),na="NA")

presence<-rbindlist(lapply(c("SPP1","C1QC","FOLR2"),function(st){high<-paste0(st,"_high_macrophage_fraction");med<-paste0(st,"_program_median_percentile");z<-s[,.(dataset_id,dataset_role,patient_id,sample_id,n_macrophages,state=st,state_high_macrophage_fraction=get(high),state_program_median_percentile=get(med))];z[,evaluable:=n_macrophages>=20];z[,state_presence:=fcase(!evaluable,"NOT_EVALUABLE_FOR_REPLICATION",state_high_macrophage_fraction>=.10,"PRESENT",default="ABSENT_OR_LOW")];z[,relative_enrichment:=fcase(!evaluable,"NOT_EVALUABLE",state_program_median_percentile>=.60,"RELATIVELY_ENRICHED",default="NOT_RELATIVELY_ENRICHED")];z}))
fwrite(presence,file.path(v5,"macrophage_state_presence_matrix.csv"),na="NA")
fwrite(presence[,.(dataset_id,dataset_role,patient_id,sample_id,state,n_macrophages,evaluable,state_program_median_percentile,relative_enrichment)],file.path(v5,"macrophage_state_relative_enrichment_matrix.csv"),na="NA")
state_repro<-presence[,.(any_present=any(state_presence=="PRESENT"),evaluable_patient=any(evaluable)),by=.(dataset_id,dataset_role,state,patient_id)][,{
  nep=sum(evaluable_patient);np=sum(any_present&evaluable_patient);status<-if(nep==0)"NOT_EVALUABLE"else if(np>=2)"REPLICATED"else if(np==1)"SUPPORTIVE_SINGLE_PATIENT"else if(nep>=2)"NOT_REPLICATED"else"NOT_EVALUABLE";list(n_evaluable_patients=nep,n_present_patients=np,present_patients=paste(patient_id[any_present&evaluable_patient],collapse=";"),reproducibility_status=status)
},by=.(dataset_id,dataset_role,state)]
fwrite(state_repro,file.path(v5,"macrophage_state_reproducibility_calibrated.csv"),na="NA")

# Author subtypes are the primary GSE154763 validation unit.
x154<-fread(file.path(v4,"GSE154763_refined","author_annotation_with_refined_state.csv.gz"))
author<-x154[,.(n_cells=.N,n_patients=uniqueN(patient),SPP1_program_median=median(SPP1_program),C1QC_program_median=median(C1QC_program),FOLR2_program_median=median(FOLR2_program),lipid_program_median=median(lipid_associated_program),SPP1_positive_fraction=mean(SPP1_positive)),by=cell_type_original]
fwrite(author,file.path(v5,"author_subtype_program_distribution.csv"),na="NA")

p1dt<-merge(repro,spp1sample[evaluable==TRUE,.(n_enriched=sum(SPP1_relative_enrichment=="RELATIVELY_ENRICHED"),n_samples=.N),by=dataset_id],by="dataset_id",all.x=TRUE)
p1<-ggplot(p1dt,aes(dataset_id,n_present_patients,fill=SPP1_reproducibility_status))+geom_col()+geom_point(aes(y=n_enriched),shape=21,size=3,fill="white")+theme_bw()+labs(title="SPP1 presence (bars) and relative enrichment (white points)",x=NULL,y="patient/sample count")
ggsave(file.path(v5,"01_spp1_presence_vs_relative_enrichment.png"),p1,width=9,height=5,dpi=180)
p2<-ggplot(spp1sample,aes(SPP1_positive_fraction,SPP1_high_macrophage_fraction,color=SPP1_presence,shape=evaluable))+geom_point(size=2.5)+facet_wrap(~dataset_id,scales="free")+theme_bw()+labs(title="Patient/sample SPP1 expression and high-state fractions",x="SPP1 positive fraction",y="SPP1-high macrophage fraction")
ggsave(file.path(v5,"02_patient_spp1_fraction.png"),p2,width=11,height=6,dpi=180)
p3dt<-presence[evaluable==TRUE,.(presence_fraction=mean(state_presence=="PRESENT")),by=.(dataset_id,state)]
p3<-ggplot(p3dt,aes(state,dataset_id,fill=presence_fraction))+geom_tile(color="white")+scale_fill_viridis_c(limits=c(0,1))+theme_bw()+labs(title="Macrophage state presence across evaluable samples",x=NULL,y=NULL)
ggsave(file.path(v5,"03_state_presence_matrix.png"),p3,width=7,height=4.5,dpi=180)
message("Macrophage evidence recalibration complete")
