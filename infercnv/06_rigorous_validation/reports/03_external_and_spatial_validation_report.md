# P1 execution report from available tables

Generated: 2026-07-11T17:36:30

## Scope

This run executes P1 analyses that are possible from existing local outputs. It does not overwrite the P0 warning: integrated_oc lacks `patient_id`, and current CNV clone labels are sample/sample_type/batch confounded.

## P1 status

| Task | Status | Main output | Limitation |
|---|---|---|---|
| P1-1 external scRNA projection | executed_from_existing_tables_partial | `06_external_scrna_projection/external_scrna_patient_level_associations.csv` | expression/signature projection; not validated CNV state transfer after patient-wise P0 |
| P1-2 LR competition | executed_from_existing_tables_partial | `07_ligand_receptor_competition/external_lr_competition_axis_ranking.csv` | SPP1 axes are candidate axes; integrated clone specificity remains sample-confounded |
| P1-3 full NicheNet ligand-target | blocked_not_completed | `08_nichenet_ligand_target/nichenet_execution_blocker.csv` | requires patient-level stable 02/04 DEG from P0/P0-5 |
| P1-4 spatial deconvolution/neighborhood | executed_sample_level_meta_partial | `09_spatial_deconvolution_neighborhood/spatial_random_effects_meta_analysis.csv` | uses existing scores and neighborhood table; no full deconvolution method found |
| P1-5 virtual perturbation | executed_score_dependency_summary | `10_virtual_perturbation_causal/virtual_perturbation_score_dependency_summary.csv` | score arithmetic dependency, not causal KO proof |

## P1-1 External scRNA patient/sample-level association

External rows used:

- patient/sample summary: 22
- signature rows: 132
- axis-score rows: 132

Primary association, all patient-samples:

- `source_myeloid_fraction` vs `Subclone02_04_common_mean_target_tumor`: Spearman r = 0.0514, n = 22, BH q = 0.9956

Output: `06_external_scrna_projection/external_scrna_patient_level_associations.csv`

## P1-2 LR competition

Top external LR axes by mean axis score:

| sample_type | axis | n | mean score | rank |
|---|---|---:|---:|---:|
| tumor | MIF-CD74 | 1 | 0.0204 | 1 |
| peritoneal_implant | SPP1-ITGB1 | 21 | 0.0202 | 1 |
| peritoneal_implant | MIF-CD74 | 21 | 0.0152 | 2 |
| peritoneal_implant | SPP1-CD44 | 21 | 0.0057 | 3 |
| peritoneal_implant | APOE-LRP1 | 21 | 9.390e-04 | 4 |
| tumor | TGFB1-TGFBR1 | 1 | 1.500e-04 | 2 |
| peritoneal_implant | TGFB1-TGFBR1 | 21 | 3.279e-05 | 5 |
| tumor | APOE-LRP1 | 1 | 2.988e-05 | 3 |
| tumor | SPP1-ITGB1 | 1 | 1.864e-06 | 4 |
| tumor | SPP1-CD44 | 1 | 1.149e-06 | 5 |

Interpretation: SPP1-CD44 and SPP1-ITGB1 remain candidate axes, but the competition table also keeps MIF/APOE/TGFB/CXCL12 controls. If controls outrank SPP1 in a stratum, the SPP1 claim should be downgraded for that stratum.

## P1-3 NicheNet

Full NicheNet was not run because the required patient-level stable 02/04 DEG target set is not available after Gate A/B. A blocker table was written to:

`08_nichenet_ligand_target/nichenet_execution_blocker.csv`

## P1-4 Spatial sample-level meta-analysis

All spatial samples random-effects meta-analysis:

- k = 15
- pooled Spearman r = 0.1356
- 95% CI = [0.0364, 0.2321]
- BH q = 0.0113
- I2 = 0.9593

GSE203612 coordinate-available subset:

- k = 3
- pooled Spearman r = 0.1068
- 95% CI = [-0.1953, 0.3903]
- BH q = 0.4910

Neighborhood enrichment meta-analysis:

- pooled enrichment ratio = 1.0329
- 95% CI = [0.8944, 1.1929]

Output: `09_spatial_deconvolution_neighborhood/spatial_random_effects_meta_analysis.csv`

## P1-5 Virtual perturbation

Existing virtual KO tables were summarized as expression-score dependency, not causal perturbation. Output:

`10_virtual_perturbation_causal/virtual_perturbation_score_dependency_summary.csv`

## Conclusion Classification

- Partially supported: spatial co-expression gradients across existing spatial samples; peritoneal-implant external LR ranking where SPP1-ITGB1 ranks first and SPP1-CD44 remains in the candidate set.
- Not supported in this execution: direct external patient-sample association between SPP1/myeloid fraction and 02/04-like tumor signature (`source_myeloid_fraction` vs `Subclone02_04_common_mean_target_tumor`, Spearman r = 0.0514, BH q = 0.9956).
- Not supported as complete: cross-patient CNV-state mechanism, full NicheNet mechanism, causal KO directionality.
- Cannot judge from current data structure: whether integrated_oc CNV_Subclone_02/04 are true cross-patient programs.
