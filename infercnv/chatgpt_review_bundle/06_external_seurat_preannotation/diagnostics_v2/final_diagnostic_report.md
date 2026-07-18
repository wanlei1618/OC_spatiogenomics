# Final diagnostic report

## Scope and safeguards

Only GSE147082, GSE154600, and GSE158722 were processed. GSE154763 was never run. Existing stage 06 `results/` and local result folders were not overwritten. No final `cell_type` or `cell_subtype` was assigned. All marker tests used the RNA assay.

## 1-3. Why percent.mt was zero, whether it is usable, and feature sources

- GSE147082: the prepared row names use the R `make.names` dot form (`MT.*`), while the old workflow searched only `^MT-`. The repaired calculation explicitly used 33 mitochondrial features. Availability: available_recomputed_explicit_features. Retained after repaired mt QC: 6993/7645.
- GSE158722: the prepared common-gene matrix contains no mitochondrial features, but original per-patient raw files were audited. Patients with `MT-` features were explicitly recalculated from their own raw counts; patients lacking credible mt features retain `NA` and skip mt filtering. No genes were added to the prepared matrix. Availability: partially_available_from_original_raw_files_some_patients_lack_mt_features; cell-level available fraction: 0.9889; unavailable raw-source patients: P01, P02, P03.

## 4. Strong sample-dominant clusters in GSE154600

- Cluster 0: T76 (99.8%, n=4169)
- Cluster 1: T59 (97.0%, n=3186)
- Cluster 10: T59 (95.0%, n=843)
- Cluster 14: T59 (100.0%, n=432)
- Cluster 16: T77 (99.7%, n=386)
- Cluster 17: T76 (99.7%, n=371)
- Cluster 23: T76 (100.0%, n=258)
- Cluster 25: T76 (99.5%, n=204)
- Cluster 26: T77 (99.4%, n=175)
- Cluster 3: T76 (98.1%, n=3102)
- Cluster 4: T59 (99.3%, n=2600)
- Cluster 5: T77 (99.8%, n=1980)
- Cluster 6: T77 (99.7%, n=1706)
- Cluster 7: T76 (100.0%, n=1511)
- Cluster 8: T59 (99.9%, n=1510)
- Cluster 9: T90 (99.6%, n=977)

## 5. Strong patient/timepoint-dominant clusters in GSE158722

- Cluster 11: P11_Pre-treatment (98.4%, n=3359)
- Cluster 12: P23_Post-treatment (99.4%, n=2789)
- Cluster 16: P14_Pre-treatment (99.8%, n=2347)
- Cluster 17: P13_Pre-treatment (98.8%, n=2107)
- Cluster 19: P15_Pre-treatment (90.2%, n=1771)
- Cluster 20: P08_Time3 (93.7%, n=1726)
- Cluster 22: P18_Post-treatment (98.8%, n=1652)
- Cluster 24: P19_Post-treatment (96.5%, n=1249)
- Cluster 26: P24_Post-treatment (100.0%, n=619)
- Cluster 28: P16_Pre-treatment (99.0%, n=298)
- Cluster 29: P01_Time1 (100.0%, n=215)
- Cluster 3: P21_Post-treatment (86.4%, n=6466)
- Cluster 30: P01_Time2 (100.0%, n=162)
- Cluster 31: P03_Time3 (97.9%, n=146)
- Cluster 32: P03_Time2 (100.0%, n=121)
- Cluster 4: P07_Time2 (87.7%, n=6365)
- Cluster 5: P12_Pre-treatment (99.2%, n=5954)
- Cluster 6: P08_Time1 (82.9%, n=4977)

## 6. Biological versus technical interpretation

Dominance is an audit flag, not a deletion rule. Strong clusters carrying epithelial/tumor programs (for example EPCAM/KRT/WFDC2/MSLN) are retained as candidate patient-specific malignant states. Clusters whose dominance coincides with abnormal nCount, nFeature, percent.mt, or doublet score are flagged as technical/QC suspects. Shared immune and stromal lineages are evaluated with correction sensitivity analyses.

- GSE154600 / `mixed_or_uncertain` (n=22): clusters 0, 1, 10, 13, 14, 15, 16, 17, 18, 2, 20, 21, 22, 23, 25, 27, 28, 29, 3, 4, 6, 7
- GSE154600 / `likely_shared_lineage` (n=4): clusters 11, 12, 19, 24
- GSE154600 / `likely_patient_specific_tumor_state` (n=4): clusters 26, 5, 8, 9
- GSE158722 / `mixed_or_uncertain` (n=24): clusters 0, 1, 10, 11, 14, 16, 17, 18, 19, 2, 21, 23, 24, 26, 27, 29, 30, 31, 32, 33, 4, 5, 6, 9
- GSE158722 / `likely_patient_specific_tumor_state` (n=5): clusters 12, 20, 22, 28, 3
- GSE158722 / `likely_shared_lineage` (n=5): clusters 13, 15, 25, 7, 8

## 7-8. Harmony/RPCA effects and lineage-preserving choices

- Strategy status counts: FAILED_CONTINUED 8; RESUMED_EXISTING 42; SKIPPED_MEMORY_GUARD_COMBINED_SCOPE 2
- GSE154600 / B_Plasma_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.8558.
- GSE154600 / Cycling_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.7617.
- GSE154600 / Endothelial_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.9615.
- GSE154600 / Epithelial_like: `A_uncorrected` - Preserve the uncorrected patient-specific epithelial structure as the primary result. Harmony/RPCA remain parallel shared-state sensitivity outputs when completed.
- GSE154600 / Fibroblast_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.196.
- GSE154600 / Myeloid_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.3755.
- GSE154600 / T_NK_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.4825.
- GSE158722 / B_Plasma_like: `B_harmony_patient` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.5241.
- GSE158722 / Cycling_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.3964.
- GSE158722 / Epithelial_like: `A_uncorrected` - Preserve the uncorrected patient-specific epithelial structure as the primary result. Harmony/RPCA remain parallel shared-state sensitivity outputs when completed.
- GSE158722 / Fibroblast_like: `A_uncorrected` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.8188.
- GSE158722 / Myeloid_like: `A_uncorrected` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=-0.008601.
- GSE158722 / T_NK_like: `B_harmony_patient` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.1891.

- GSE154600 / B_Plasma_like / `B_harmony_sample`: median dominant-sample fraction 0.85 -> 0.436; entropy 0.327 -> 0.809; ARI vs A=0.406; rare-population retention=0.541.
- GSE154600 / Cycling_like / `B_harmony_sample`: median dominant-sample fraction 0.806 -> 0.531; entropy 0.378 -> 0.745; ARI vs A=0.455; rare-population retention=NA.
- GSE154600 / Endothelial_like / `B_harmony_sample`: median dominant-sample fraction 0.853 -> 0.376; entropy 0.301 -> 0.748; ARI vs A=0.446; rare-population retention=NA.
- GSE154600 / Epithelial_like / `A_uncorrected`: median dominant-sample fraction 0.997 -> 0.997; entropy 0.0457 -> 0.0457; ARI vs A=1; rare-population retention=1.
- GSE154600 / Fibroblast_like / `B_harmony_sample`: median dominant-sample fraction 0.958 -> 0.818; entropy 0.218 -> 0.365; ARI vs A=0.781; rare-population retention=0.725.
- GSE154600 / Myeloid_like / `B_harmony_sample`: median dominant-sample fraction 0.966 -> 0.677; entropy 0.194 -> 0.554; ARI vs A=0.48; rare-population retention=0.429.
- GSE154600 / T_NK_like / `B_harmony_sample`: median dominant-sample fraction 0.994 -> 0.561; entropy 0.0922 -> 0.567; ARI vs A=0.588; rare-population retention=0.273.
- GSE158722 / B_Plasma_like / `B_harmony_patient`: median dominant-sample fraction 0.787 -> 0.631; entropy 0.254 -> 0.304; ARI vs A=0.815; rare-population retention=NA.
- GSE158722 / Cycling_like / `B_harmony_sample`: median dominant-sample fraction 0.875 -> 0.416; entropy 0.247 -> 0.548; ARI vs A=0.198; rare-population retention=0.221.
- GSE158722 / Epithelial_like / `A_uncorrected`: median dominant-sample fraction 0.956 -> 0.956; entropy 0.171 -> 0.171; ARI vs A=1; rare-population retention=1.
- GSE158722 / Fibroblast_like / `A_uncorrected`: median dominant-sample fraction 0.455 -> 0.455; entropy 0.435 -> 0.435; ARI vs A=1; rare-population retention=1.
- GSE158722 / Myeloid_like / `A_uncorrected`: median dominant-sample fraction 0.892 -> 0.892; entropy 0.167 -> 0.167; ARI vs A=1; rare-population retention=1.
- GSE158722 / T_NK_like / `B_harmony_patient`: median dominant-sample fraction 0.985 -> 0.753; entropy 0.0981 -> 0.191; ARI vs A=0.907; rare-population retention=0.6.

### RPCA completed comparisons and blockers

- Completed RPCA / GSE154600 / B_Plasma_like: median dominant-sample fraction=0.624; entropy=0.539; ARI vs A=0.524; rare-population retention=0.518.
- Completed RPCA / GSE154600 / Epithelial_like: median dominant-sample fraction=0.901; entropy=0.375; ARI vs A=0.527; rare-population retention=0.364.
- Completed RPCA / GSE154600 / Fibroblast_like: median dominant-sample fraction=0.78; entropy=0.42; ARI vs A=0.649; rare-population retention=0.514.
- Completed RPCA / GSE154600 / Myeloid_like: median dominant-sample fraction=0.655; entropy=0.619; ARI vs A=0.275; rare-population retention=0.282.
- Completed RPCA / GSE154600 / T_NK_like: median dominant-sample fraction=0.809; entropy=0.47; ARI vs A=0.436; rare-population retention=0.371.
- RPCA blocker / GSE154600 / Cycling_like: number of items to replace is not a multiple of replacement length
- RPCA blocker / GSE154600 / Endothelial_like: number of items to replace is not a multiple of replacement length
- RPCA blocker / GSE158722 / B_Plasma_like: RPCA requires >=10 cells in every sample within lineage
- RPCA blocker / GSE158722 / Cycling_like: RPCA requires >=10 cells in every sample within lineage
- RPCA blocker / GSE158722 / Epithelial_like: The total size of the 10 globals exported for future expression ('FUN()') is 6.34 GiB. This exceeds the maximum allowed size 6.00 GiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help(""future.globals.maxSize"", package = ""future"") for how to adjust or remove the default threshold via an R option The three largest globals are 'FUN' (4.46 GiB of class 'function'), 'object.list' (1.89 GiB of class 'list') and 'NNHelper' (9.67 KiB of class 'function')
- RPCA blocker / GSE158722 / Fibroblast_like: RPCA requires >=10 cells in every sample within lineage
- RPCA blocker / GSE158722 / Myeloid_like: RPCA requires >=10 cells in every sample within lineage
- RPCA blocker / GSE158722 / T_NK_like: RPCA requires >=10 cells in every sample within lineage

## 9. Patient-specific tumor clusters to retain

- GSE154600 cluster 26: dominant sample T77 (99.4%); retained as uncorrected epithelial/tumor-state candidate.
- GSE154600 cluster 5: dominant sample T77 (99.8%); retained as uncorrected epithelial/tumor-state candidate.
- GSE154600 cluster 8: dominant sample T59 (99.9%); retained as uncorrected epithelial/tumor-state candidate.
- GSE154600 cluster 9: dominant sample T90 (99.6%); retained as uncorrected epithelial/tumor-state candidate.
- GSE158722 cluster 12: dominant sample P23_Post-treatment (99.4%); retained as uncorrected epithelial/tumor-state candidate.
- GSE158722 cluster 20: dominant sample P08_Time3 (93.7%); retained as uncorrected epithelial/tumor-state candidate.
- GSE158722 cluster 22: dominant sample P18_Post-treatment (98.8%); retained as uncorrected epithelial/tumor-state candidate.
- GSE158722 cluster 28: dominant sample P16_Pre-treatment (99.0%); retained as uncorrected epithelial/tumor-state candidate.
- GSE158722 cluster 3: dominant sample P21_Post-treatment (86.4%); retained as uncorrected epithelial/tumor-state candidate.
All epithelial-like uncorrected clusters remain available as the primary view. Corrected epithelial results, when successful, are parallel shared-state sensitivity outputs and never replace the uncorrected result.

## 10. Version for subsequent manual cell-type review

Use the strategy listed per dataset and non-epithelial lineage in `strategy_comparison/recommended_strategy_by_dataset_and_lineage.csv`. For epithelial-like cells, use `A_uncorrected`. The blank `manual_annotation_template.csv` files are the handoff point.

## 11. Remaining limitations

- iLISI and kBET remain skipped when their packages are unavailable; their status is non-blocking and recorded.
- GSE158722 current-object PC/UMAP metadata audit may be memory-guarded; newly generated lineage strategy embeddings provide the sensitivity views.
- Broad lineages are provisional diagnostic strata, not final annotations.
- Marker overlap is a lineage-level gene-set sensitivity metric; cluster identities are not assumed to map one-to-one after correction.

- The task description flags a possible T61/T77 identity conflict, but T61 is absent from the supplied GSE154600 result metadata. This package therefore retains source sample IDs and makes no resolved patient-identity claim for T77.

## Reproducibility

Run `workflow/scripts/run_diagnostics_v2.ps1`. Large RDS/count payloads remain local under `diagnostics_v2/objects` and are excluded from GitHub.
