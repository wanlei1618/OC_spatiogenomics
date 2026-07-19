# FINAL CODEX EXECUTION SUMMARY

## Scope completed

First-stage P0 entry work was completed under:

`D:/OC_spatiogenomics/infercnv/06_rigorous_validation`

Completed actions:

- Created the requested validation directory structure.
- Created `00_config/project_config.yaml`.
- Created `CODEX_EXECUTION_PLAN.md`.
- Implemented and ran `scripts/R/01_build_checked_metadata.R`.
- Implemented and ran `scripts/R/02_clone_patient_confounding.R`.
- Built an input manifest with 650 project files.
- Built checked metadata from the selected integrated metadata CSV.
- Generated clone by sample, sample_type and batch count tables.
- Generated clone confounding metrics, Cramer's V association statistics, standardized residuals and first-stage figures.
- Generated `reports/01_clone_validity_report.md` and `reports/01_clone_validity_report.html`.

## Main result

The current metadata contains 49,326 rows and 6,507 clone-labeled cells, but it does not contain `patient_id`.

Because `patient_id` is missing, Gate A cannot judge whether CNV_Subclone_01-05 are true cross-patient CNV states. No patient mapping was inferred from sample_id.

## Sample-level warning signal

At sample level, all clone labels are dominated by sample `31V1`:

| clone | cells | detected samples | max single-sample fraction | dominant sample |
|---|---:|---:|---:|---|
| CNV_Subclone_01 | 3318 | 4 | 0.805 | 31V1 |
| CNV_Subclone_02 | 1071 | 3 | 0.992 | 31V1 |
| CNV_Subclone_03 | 1021 | 2 | 0.976 | 31V1 |
| CNV_Subclone_04 | 1001 | 2 | 0.943 | 31V1 |
| CNV_Subclone_05 | 96 | 2 | 0.990 | 31V1 |

Clone association with sample_id/sample_type/batch was detectable, with Cramer's V = 0.156 and BH-adjusted P = 5.07e-94. Because sample_id, sample_type and batch are highly overlapping in the selected table, this is a diagnostic warning rather than independent biological evidence.

## Blocker before continuing old-clone mechanism analysis

An explicit `sample_id` to `patient_id` mapping is required before current clone labels can be described as cross-patient states or used for downstream mechanism claims.

Recommended next step: add a reviewed patient mapping table to `00_config/` and rerun scripts 01 and 02 before proceeding to patient-wise CNV reconstruction or pseudobulk functional validation.

## P1 execution update

P1 was executed where possible from existing local tables, with outputs written under `06_external_scrna_projection`, `07_ligand_receptor_competition`, `08_nichenet_ligand_target`, `09_spatial_deconvolution_neighborhood`, and `10_virtual_perturbation_causal`.

Key results:

- External patient-sample association did not support direct SPP1/myeloid fraction vs 02/04-like tumor signature association: Spearman r = 0.0514, BH q = 0.9956.
- External peritoneal-implant LR ranking placed SPP1-ITGB1 first and SPP1-CD44 third among evaluated axes; tumor ranking is not stable because only one tumor sample is present.
- Spatial random-effects meta-analysis across 15 samples gave pooled Spearman r = 0.1356, 95% CI 0.0364 to 0.2321, BH q = 0.0113, with high heterogeneity.
- GSE203612 coordinate-available neighborhood enrichment was weak: pooled enrichment ratio = 1.0329, 95% CI 0.8944 to 1.1929.
- Full NicheNet remains blocked because patient-level stable 02/04 meta-DEG targets are not available after Gate A/B.

Current P1 conclusion: candidate niche evidence is partial and mixed; direct external SPP1-myeloid to 02/04-like association is not supported by the available patient-sample table.

## P2 execution update

# P2 bulk, CNV-expression and mediation report

Generated: 2026-07-11T20:12:53

## P2-1 Bulk endpoint repositioning

The five-cohort Cox interaction results were retained as negative or non-significant and are not used as the primary endpoint.

Output: `11_bulk_negative_validation/bulk_survival_interaction_negative_results_retained.csv`

Bulk signature covariation was re-positioned as a co-variation check. Main meta result for `SPP1_TAM_score` vs `KRAS_Hypoxia_score`:

- pooled Spearman r = 0.1042
- 95% CI = [0.0134, 0.1933]
- BH q = 0.0306
- I2 = 0.6098

Output: `11_bulk_negative_validation/bulk_signature_covariation_random_effects_meta.csv`

## P2-2 CNV-expression dosage independent validation status

Local outputs contain 776 positive RNA-derived CNV-expression coupling candidates at the existing threshold. These were relabelled as:

`RNA-derived CNV-expression coupling candidates`

They were not upgraded to DNA dosage drivers because no independent DNA CNV plus RNA validation table was found locally.

Output: `12_integrated_evidence_scoring/RNA_derived_CNV_expression_coupling_candidates.csv`

## P2-3 Statistical mediation exploration

Bulk-level exploratory mediation was run per cohort:

`SPP1_TAM_score -> ITGB1_CD44_tumor_score -> KRAS_Hypoxia_score`

This is explicitly statistical mediation, not causal proof and not local spatial niche evidence. Cohorts with bootstrap CI excluding 0: 5.

Output: `12_integrated_evidence_scoring/bulk_exploratory_mediation_bootstrap.csv`

## Integrated conclusion

The current evidence supports cautious, mixed candidate biology rather than a validated mechanism. The strongest limitation remains unresolved P0/Gate A: integrated_oc lacks `patient_id`, so current CNV_Subclone_02/04 cannot yet be promoted to cross-patient CNV programs.
