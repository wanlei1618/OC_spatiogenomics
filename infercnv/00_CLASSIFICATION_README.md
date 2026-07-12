# inferCNV project reclassification and ChatGPT validation input

Generated: 2026-07-10 12:05:53

Updated with 2026-07-05 Codex results: 2026-07-10 12:30

Project root: `D:\OC_spatiogenomics\infercnv`

## Hypothesis to verify

SPP1+ myeloid/macrophage cells may shape a sample-type/niche-specific invasive ecology for KRAS/hypoxia-high ovarian cancer CNV subclones through SPP1-CD44 and SPP1-ITGB1 axes.

## New directory classes

| Directory | Role |
|---|---|
| `00_raw_objects_and_infercnv` | Raw objects and inferCNV calling layer: inferCNV outputs plus integrated_oc and immune subtype RData objects. |
| `01_integrated_oc_cnv_clone_definition` | integrated_oc metadata, CNV subclone definition, TME/LR/kNN exploratory layer. |
| `02_cnv_expression_state_analysis` | CNV burden and expression-state interpretation for KRAS/hypoxia/TF proxy programs. |
| `03_spp1_cd44_itgb1_hypothesis_validation` | Focused SPP1-CD44/ITGB1 validation, CellChat/LIANA/CopyKAT/bulk Cox, and dedicated LR figures. |
| `04_sample_type_scores` | sample_type-stratified SPP1_TAM and ITGB1/CD44 tumor score analysis. |
| `05_sample_type_lr_niche_external_validation` | Main sample_type LR niche analysis with Zhang2022 scRNA and GSE203612/GSE189843 spatial validation. |
| `90_environment_install_logs` | Package install and environment logs. |
| `99_codex_results_from_c` | Source archive and migration manifests for inferCNV-related Codex outputs originally found under C drive. |

## 2026-07-05 Codex results added

The newly added 2026-07-05 Codex results were reviewed against the classified main directories.

- 58 unique files were moved into the main classified directories.
- 36 files were kept only in `99_codex_results_from_c` because same-name/same-size copies already exist in the main classified directories.
- The added unique files mainly supplement early inferCNV subclone outputs, early integrated_oc immune subtype mapping, broad LR interaction tables, object-inspection scripts, external-data preparation scripts, and OmniPath/R package logs.

Detailed records:

- `D:\OC_spatiogenomics\infercnv\99_codex_results_from_c\2026-07-05_reclassification_manifest.csv`
- `D:\OC_spatiogenomics\infercnv\99_codex_results_from_c\2026-07-05_file_inventory.csv`
- `D:\OC_spatiogenomics\infercnv\file_reclassification_manifest_2026-07-10_with_7-5.csv`

## Locked source folders still present

Two old root-level source folders were copied into classified destinations. After the 2026-07-05 update, their files are verified as covered by the classified destinations, but Windows still did not allow deleting the old root-level copies because files or folders were in use:

- `D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis`
- `D:\OC_spatiogenomics\infercnv\CNV_expression_joint_analysis_integrated_oc`

The classified copies are:

- `D:\OC_spatiogenomics\infercnv\01_integrated_oc_cnv_clone_definition\integrated_oc_plan_analysis`
- `D:\OC_spatiogenomics\infercnv\02_cnv_expression_state_analysis\CNV_expression_joint_analysis_integrated_oc`

## Key reports to read in order

1. `D:\OC_spatiogenomics\infercnv\01_integrated_oc_cnv_clone_definition\integrated_oc_plan_analysis\integrated_oc_infercnv_plan_task_breakdown_and_results.md`
2. `D:\OC_spatiogenomics\infercnv\02_cnv_expression_state_analysis\CNV_expression_joint_analysis_integrated_oc\CNV_expression_joint_analysis_task_breakdown_and_results.md`
3. `D:\OC_spatiogenomics\infercnv\03_spp1_cd44_itgb1_hypothesis_validation\SPP1_ITGB1_CD44_hypothesis_validation\SPP1_ITGB1_CD44_hypothesis_validation_report.md`
4. `D:\OC_spatiogenomics\infercnv\03_spp1_cd44_itgb1_hypothesis_validation\SPP1_ITGB1_CD44_hypothesis_validation_complete\SPP1_ITGB1_CD44_complete_validation_report.md`
5. `D:\OC_spatiogenomics\infercnv\03_spp1_cd44_itgb1_hypothesis_validation\LR_SPP1_ITGB1_CD44_figures\README_SPP1_ITGB1_CD44_figures.md`
6. `D:\OC_spatiogenomics\infercnv\04_sample_type_scores\sample_type_SPP1_TAM_ITGB1_CD44_scores\sample_type_SPP1_TAM_ITGB1_CD44_score_report.md`
7. `D:\OC_spatiogenomics\infercnv\05_sample_type_lr_niche_external_validation\sample_type_LR_niche_analysis\sample_type_LR_niche_analysis_report.md`

## Core evidence files

- `01_integrated_oc_cnv_clone_definition\integrated_oc_plan_analysis\tables\integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv`
- `01_integrated_oc_cnv_clone_definition\integrated_oc_plan_analysis\tables\focused_LR_axes_involving_cnv_subclones.csv`
- `02_cnv_expression_state_analysis\CNV_expression_joint_analysis_integrated_oc\tables\CNV_clone_burden_summary.csv`
- `02_cnv_expression_state_analysis\CNV_expression_joint_analysis_integrated_oc\tables\tumor_cell_CNV_burden_signature_TF_scores.csv`
- `03_spp1_cd44_itgb1_hypothesis_validation\SPP1_ITGB1_CD44_hypothesis_validation_complete\tables\focused_LR_predefined_axes_senders_to_Subclone02_04.csv`
- `03_spp1_cd44_itgb1_hypothesis_validation\SPP1_ITGB1_CD44_hypothesis_validation_complete\tables\LIANA_consensus_focus_pair_presence.csv`
- `03_spp1_cd44_itgb1_hypothesis_validation\SPP1_ITGB1_CD44_hypothesis_validation_complete\tables\bulk_meta_interaction_summary.csv`
- `04_sample_type_scores\sample_type_SPP1_TAM_ITGB1_CD44_scores\tables\sample_type_score_summary.csv`
- `05_sample_type_lr_niche_external_validation\sample_type_LR_niche_analysis\tables\sample_type_LR_opportunity_scores_primary_axes_target_Subclone02_04.csv`
- `05_sample_type_lr_niche_external_validation\sample_type_LR_niche_analysis\tables\external_scRNA_all_datasets_sampletype_axis_scores.csv`
- `05_sample_type_lr_niche_external_validation\sample_type_LR_niche_analysis\tables\spatial_correlation_SPP1_myeloid_Target_axis.csv`
- `05_sample_type_lr_niche_external_validation\sample_type_LR_niche_analysis\tables\spatial_neighborhood_enrichment.csv`
- `05_sample_type_lr_niche_external_validation\sample_type_LR_niche_analysis\tables\limitations_summary.csv`

## Suggested ChatGPT validation prompt

Please verify whether the reorganized files support this logic:

1. inferCNV/integrated_oc define five ovarian cancer CNV subclones.
2. CNV_Subclone_02 and CNV_Subclone_04 are the main target clones because they show KRAS-active and hypoxia/high-CNV/immune-modulatory features.
3. SPP1+ myeloid/macrophage groups show candidate LR communication toward CNV_Subclone_02/04 through SPP1-CD44 and SPP1-ITGB1.
4. SPP1-CD44 is better supported by CellChat/LIANA-style evidence; SPP1-ITGB1 is mainly supported by expression-product/Connectome-like evidence and should be interpreted as an adhesion/integrin candidate axis.
5. integrated_oc sample_type findings are exploratory because sample_id is one-to-one with sample_type.
6. Zhang2022 ovarian scRNA and GSE203612/GSE189843 spatial data provide external expression-level validation for tumor/peritoneal-implant-like niche specificity.
7. Bulk Cox interaction results are not stable and should not be the primary endpoint.
8. Final claim should be: SPP1-CD44/ITGB1 is a candidate sample_type/niche-dependent tumor-myeloid communication program, not proof of physical ligand-receptor binding or universal OS prediction.

## Generated files

- `D:\OC_spatiogenomics\infercnv\file_reclassification_manifest.csv`
- `D:\OC_spatiogenomics\infercnv\00_CLASSIFICATION_README.md`
- `D:\OC_spatiogenomics\infercnv\CHATGPT_VALIDATION_INPUT.md`
