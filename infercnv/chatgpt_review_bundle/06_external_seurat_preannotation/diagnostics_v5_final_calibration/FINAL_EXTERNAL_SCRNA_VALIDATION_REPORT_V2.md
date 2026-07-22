# Final external scRNA validation report v2

## Statistical definitions

SPP1/state presence, cross-patient reproducibility and relative enrichment are distinct outputs. Presence uses evaluable samples and high-state fractions; reproducibility counts patients with at least one evaluable PRESENT sample; within-dataset median percentile >=0.60 is used only for relative enrichment. Cells are never treated as independent statistical replicates.

## 1. SPP1 presence, reproducibility and relative enrichment

- GSE154600: 5/5 evaluable patients PRESENT; reproducibility=REPLICATED; relatively enriched samples=none at percentile >=0.60. This dataset supports cross-patient SPP1 macrophage presence even though no patient median exceeds the relative-enrichment percentile threshold.
- GSE147082: 4/4 evaluable patients PRESENT; reproducibility=REPLICATED; relatively enriched samples=PT-3232/PT-3232.
- GSE158722: 1/2 evaluable patients PRESENT; reproducibility=SUPPORTIVE_SINGLE_PATIENT; relatively enriched samples=none at percentile >=0.60. P07_Time1 is evaluable and PRESENT; P04_Time3 is evaluable but ABSENT_OR_LOW under the joint positive/high-fraction rule. All other small timepoints remain descriptive and are NOT_EVALUABLE_FOR_REPLICATION, not negative.
- GSE151214: 4/5 evaluable patients PRESENT; reproducibility=REPLICATED; relatively enriched samples=FT01/FT01; FT-SA23510/FT-SA23510; FT-SA23515/FT-SA23515; normal-reference evidence only.
- GSE154763: 2/2 evaluable patients PRESENT; reproducibility=REPLICATED; relatively enriched samples=P20190304/OV-P20190304-T; author-annotated myeloid-reference evidence only.

## 2. C1QC and FOLR2 calibration

- GSE154600: C1QC REPLICATED (4/5); FOLR2 REPLICATED (5/5).
- GSE147082: C1QC REPLICATED (3/4); FOLR2 REPLICATED (4/4).
- GSE158722: C1QC SUPPORTIVE_SINGLE_PATIENT (1/2); FOLR2 SUPPORTIVE_SINGLE_PATIENT (1/2); small-sample support only.
SPP1, C1QC and FOLR2 can coexist at the presence level. Relative enrichment identifies which patient/sample is shifted within its own dataset and must not be used to negate presence.

## 3. Why GSE154600 is not SPP1-negative

All five evaluable patients are PRESENT: T59 positive=0.947, high=0.301; T76 positive=0.842, high=0.132; T77 positive=0.408, high=0.192; T89 positive=0.599, high=0.288; T90 positive=0.305, high=0.337. Therefore SPP1 cross-patient presence is REPLICATED. The absence of a sample median percentile >=0.60 only means no patient is relatively enriched by that separate criterion.

## 4. GSE147082 cluster 6

True FindAllMarkers results support: subcluster 0 CD8_effector_T n=298 markers=INPP4B;CD8A;AC015849.1;CCL5; subcluster 1 Cycling_CD8_effector_T n=99 markers=RRM2;FBXO43;CDCA2;UBE2C;SHCBP1;BUB1B;HJURP;CCNB1;KIF23;CDCA5;CCNA2;TOP2A;KIF18B;CDK1;DLGAP5;BIRC5;FAM111B;PCLAF;GTSE1;AURKA; subcluster 2 Gamma_delta_T n=73 markers=SH2D1B;FCER1G;NCAM1;NCR1;TYROBP;KLRB1;TXK;KIR2DL4;GNLY;KLRC1;LINC00299;ITGAX;AC092821.3;MATK;ATP8B4;KRT86;CLNK;LAT2;PLCG2;RIN3. The former two identically named Cytotoxic_T groups are now CD8_effector_T and Cycling_CD8_effector_T; cycling is interpreted as a state. Gamma_delta_T remains distinct.

## 5. PT-2834 patient-internal CNV-like sensitivity

Cluster 4: Mesenchymal_stromal_candidate; epithelial-profile correlation=0.070. Cluster 7: COL2A1_positive_chondrocyte_like_fibroblast_candidate. Method=exploratory_patient_internal_cnv_like. Only 15 broad-lineage epithelial reference cells were available, so this is a candidate-level sensitivity result, not definitive malignant-cell classification.

## 6. Epithelial CD44/ITGB1 receiver context

- GSE154600: CD44=DETECTED_LOW, ITGB1=SUPPORTED, dual-positive=DETECTED_LOW, context=DESCRIPTIVE_COEXISTENCE.
- GSE147082: CD44=DETECTED_LOW, ITGB1=SUPPORTED, dual-positive=DETECTED_LOW, context=DESCRIPTIVE_COEXISTENCE.
- GSE158722: CD44=NOT_AVAILABLE, ITGB1=NOT_AVAILABLE, dual-positive=NOT_AVAILABLE, context=NOT_EVALUABLE.
ITGB1 receiver support is stronger than CD44 in the evaluable tumor datasets; dual-positive fractions are detected at low levels. Receivers are broad-lineage Epithelial and are not described as uniformly confirmed malignant cells. The result is an SPP1-associated epithelial CD44/ITGB1 context, not a proven direct ligand-receptor interaction.

## 7. GSE154763 author-subtype validation

Macro_SPP1: SPP1 program median=0.277. Macro_C1QC: SPP1=0.006, C1QC=0.585, FOLR2=0.430. Primary validation is based on whole author-subtype distributions; single-cell module_state remains sensitivity analysis only.

## 8. Evidence freeze and limitations

The v2 matrix removes eligible_for_primary_conclusion and instead records tumor-context eligibility and explicit support flags. Dataset roles remain fixed. GSE158722 is INCONCLUSIVE rather than negative. No result proves direct binding, causal regulation or directionality; PT-2834 CNV-like classification and all sender/receiver results remain descriptive and require orthogonal or experimental validation.
