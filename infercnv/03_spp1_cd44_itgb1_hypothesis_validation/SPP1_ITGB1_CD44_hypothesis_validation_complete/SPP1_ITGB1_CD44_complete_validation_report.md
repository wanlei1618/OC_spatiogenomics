# SPP1-ITGB1/CD44 complete validation update

Output directory: `D:/OC_spatiogenomics/infercnv/SPP1_ITGB1_CD44_hypothesis_validation_complete`

## 5.4 Bulk Cox interaction model

Five cohorts were analyzed: TCGA_OV, GSE102094, GSE32062, GSE140082, and GSE49997. Expression matrices were read as sample-by-gene tables, Entrez IDs were mapped to gene symbols with org.Hs.eg.db, and cohort-wise z-score signatures were computed.

Primary interaction term: `SPP1_TAM_score x ITGB1_CD44_tumor_score`.

### Fixed-effect interaction summary

  endpoint    model n_cohorts         beta         se        HR    CI_low
1       OS     core         5 -0.001272303 0.03979452 0.9987285 0.9237909
2       OS adjusted         5  0.010869080 0.04123243 1.0109284 0.9324437
3      PFS     core         4 -0.040939185 0.03570651 0.9598875 0.8950068
4      PFS adjusted         4 -0.020158084 0.03610334 0.9800437 0.9130902
   CI_high   p_value
1 1.079745 0.9744945
2 1.096019 0.7920842
3 1.029472 0.2515689
4 1.051907 0.5766095

### Per-cohort interaction results

      cohort endpoint    model   n events        HR    CI_low  CI_high
1    TCGA_OV       OS     core 377    205 1.0769304 0.9532076 1.216712
2    TCGA_OV       OS adjusted 377    205 1.0889161 0.9593113 1.236031
3  GSE102094       OS     core  84     19 1.0358291 0.6821084 1.572979
4  GSE102094       OS adjusted  84     19 0.9871011 0.6077389 1.603269
5  GSE102094      PFS     core  85     57 0.8195370 0.6015846 1.116453
6  GSE102094      PFS adjusted  85     57 0.8483453 0.6116524 1.176632
7   GSE32062       OS     core 260    121 0.9459843 0.8156019 1.097210
8   GSE32062       OS adjusted 260    121 0.9394997 0.8123819 1.086508
9   GSE32062      PFS     core 260    193 1.0290705 0.9222768 1.148230
10  GSE32062      PFS adjusted 260    193 1.0265705 0.9231370 1.141593
11 GSE140082       OS     core 380     96 0.9412275 0.7909790 1.120016
12 GSE140082       OS adjusted 380     96 0.9791449 0.8055728 1.190115
13 GSE140082      PFS     core 380    235 0.9043852 0.8058831 1.014927
14 GSE140082      PFS adjusted 380    235 0.9353698 0.8278484 1.056856
15  GSE49997       OS     core 194     57 0.9358890 0.7090281 1.235336
16  GSE49997       OS adjusted 194     57 0.9939759 0.7533403 1.311476
17  GSE49997      PFS     core 194    124 0.9690718 0.8187892 1.146938
18  GSE49997      PFS adjusted 194    124 0.9904262 0.8356725 1.173838
     p_value
1  0.2339244
2  0.1876753
3  0.8688250
4  0.9581601
5  0.2070790
6  0.3244237
7  0.4630162
8  0.4001366
9  0.6082238
10 0.6284132
11 0.4948502
12 0.8323435
13 0.0876116
14 0.2835447
15 0.6399223
16 0.9659223
17 0.7148081
18 0.9116354

## 5.1-5.5 status

  section
1     5.1
2     5.2
3     5.3
4     5.4
5     5.5
                                                                    analysis
1                                  malignant identity and receptor stability
2                                                       pseudo-bulk DEG/GSEA
3 LR reproduction with CellChat/Connectome plus LIANA/NicheNet focused checks
4                                          five-cohort Cox interaction model
5                                     spatial validation marker/scoring plan
                                                                                                              status
1 completed from integrated_oc plus prior CNV clone metadata; CopyKAT sampled validation completed, CaSpER package load confirmed but full run blocked by local biomaRt SSL
2                                                                 completed and copied from prior focused validation
3          completed for Connectome-like and CellChat; LIANA consensus and NicheNet local gene-universe checks added
4                                        completed in this run for TCGA_OV, GSE102094, GSE32062, GSE140082, GSE49997
5                                                                 completed and copied from prior focused validation

## Key output tables

- `tables/bulk_signature_scores_merged.csv`
- `tables/bulk_cox_model_status.csv`
- `tables/bulk_cox_model_coefficients.csv`
- `tables/bulk_cox_interaction_results.csv`
- `tables/bulk_meta_interaction_summary.csv`
- `tables/R_package_status_for_complete_validation.csv`
- `tables/new_R_package_focused_validation_status.csv`
- `tables/CopyKAT_sampled_prediction_summary_by_known_group.csv`
- `tables/LIANA_consensus_focus_pair_presence.csv`
- `tables/NicheNet_focus_geneinfo_presence.csv`
- copied single-cell validation tables from the prior focused validation folder

## Key output figures

- `figures/bulk_OS_adjusted_interaction_forest.png`
- `figures/bulk_OS_core_interaction_forest.png`
- `figures/bulk_PFS_adjusted_interaction_forest.png`
- `figures/bulk_PFS_core_interaction_forest.png`
- copied single-cell validation figures from the prior focused validation folder

## New package focused checks

CopyKAT was run on a sampled integrated_oc subset with immune cells as normal references. The summary supported malignant/aneuploid status for the target clones: Subclone_02 was 80/80 aneuploid, Subclone_04 was 79/80 aneuploid, while immune references were mostly diploid (132 diploid, 0 aneuploid, 28 not defined).

LIANA Consensus contained 3 of the predefined focus pairs: SPP1-CD44, APOE-LRP1, and CD47-SIRPA. SPP1-ITGB1 was not present in LIANA Consensus, so the SPP1-ITGB1 evidence remains from the expression-based Connectome-like focused LR table rather than LIANA consensus curation.

NicheNet local package geneinfo contained all 25 checked focus genes. A full NicheNet ligand-target model was not rerun because the package installation does not include the large ligand-target prior matrix locally; the result table records the local gene-universe check.

CaSpER is installed, but full rerun was not executed because its annotation generation triggered a local biomaRt SSL certificate error in this R environment.
