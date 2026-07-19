# Patient-level functional validation report

Generated: 2026-07-11T20:38:38

## Status

P0-5 / Gate B is not fully completed because the required `patient_id` field is missing from the integrated_oc metadata and current CNV clone labels failed preliminary sample-confounding checks.

## What is available

- Existing sample-level CNV-expression outputs are present in `CNV_expression_joint_analysis_integrated_oc`.
- Existing functional score summaries and pseudobulk-style tables can be used only as exploratory sample-level evidence.
- These results must not be reported as patient-level biological replication or as primary P values from independent patients.

## Why this is not marked complete

The task book requires patient/sample-level pseudobulk and leave-one-patient-out/meta-analysis. Because `patient_id` is absent and only 4 sample IDs are available in the integrated metadata, Gate B cannot be judged as a validated patient-level 02/04 functional state.

## Required to complete

1. Provide a reviewed sample-to-patient mapping.
2. Rebuild patient-wise CNV programs or verify existing clone labels after Gate A.
3. Re-run pseudobulk DEG/GSEA/PROGENy/TF activity with patient/sample as the replicate unit.
4. Run leave-one-patient-out sensitivity analysis.

## Conclusion

Classification: 因数据结构无法判断 / blocked by missing patient-level replication.
