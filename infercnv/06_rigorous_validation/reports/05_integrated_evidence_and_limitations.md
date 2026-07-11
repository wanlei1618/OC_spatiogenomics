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
