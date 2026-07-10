# inferCNV ChatGPT review bundle

This package is a compact review set extracted from:

`D:\OC_spatiogenomics\infercnv`

It excludes large binary/intermediate files such as `.RData`, `.rds`, `.infercnv_obj`, raw inferCNV matrices, install logs, and duplicated source folders.

## Recommended reading order

1. `CHATGPT_VALIDATION_INPUT.md`
2. `01_integrated_oc_cnv_clone_definition/integrated_oc_plan_analysis/integrated_oc_infercnv_plan_task_breakdown_and_results.md`
3. `02_cnv_expression_state_analysis/CNV_expression_joint_analysis_integrated_oc/CNV_expression_joint_analysis_task_breakdown_and_results.md`
4. `03_spp1_cd44_itgb1_hypothesis_validation/SPP1_ITGB1_CD44_hypothesis_validation_complete/SPP1_ITGB1_CD44_complete_validation_report.md`
5. `04_sample_type_scores/sample_type_SPP1_TAM_ITGB1_CD44_scores/sample_type_SPP1_TAM_ITGB1_CD44_score_report.md`
6. `05_sample_type_lr_niche_external_validation/sample_type_LR_niche_analysis/sample_type_LR_niche_analysis_report.md`

## What to review

Please evaluate whether the provided reports and tables support the central claim:

SPP1-positive myeloid/macrophage cells may shape a sample-type/niche-specific invasive ecology for KRAS/hypoxia-high ovarian cancer CNV subclones through SPP1-CD44 and SPP1-ITGB1 axes.

Focus especially on:

- Whether inferCNV/integrated_oc evidence supports the five CNV subclones.
- Whether CNV_Subclone_02 and CNV_Subclone_04 are justified as target clones.
- Whether SPP1-CD44 is better supported than SPP1-ITGB1 by ligand-receptor evidence.
- Whether sample_type effects are exploratory because sample_id is one-to-one with sample_type.
- Whether external scRNA/spatial validation supports expression-level niche specificity.
- Whether bulk Cox results should be treated as secondary and unstable.

## Suggested final output

Please return:

1. A verdict: supported, partially supported, or not supported.
2. Strongest evidence.
3. Weakest links or confounders.
4. Claims that should be softened.
5. Additional analyses needed before manuscript-level confidence.
