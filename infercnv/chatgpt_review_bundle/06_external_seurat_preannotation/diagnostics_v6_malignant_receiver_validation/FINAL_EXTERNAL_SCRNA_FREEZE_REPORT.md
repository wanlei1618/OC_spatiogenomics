# Final external scRNA evidence freeze report

## Frozen scope

This report freezes the external scRNA preprocessing and evidence definitions after targeted count-based state detection, GSE147082 cluster 6 correction, formal CopyKAT assessment, and malignancy-stratified receiver analysis. No dataset-wide QC, clustering, or broad annotation was rerun.

## Dataset-level evidence

- GSE147082: SPP1 transcript=DETECTED_REPLICATED; program=SUPPORTED_REPLICATED; reproducibility=REPLICATED; malignancy=NOT_AVAILABLE; receiver context=NOT_EVALUABLE_OR_NOT_ESTABLISHED.
- GSE151214: SPP1 transcript=DETECTED_REPLICATED; program=SUPPORTED_REPLICATED; reproducibility=REPLICATED; malignancy=NOT_APPLICABLE_REFERENCE; receiver context=REFERENCE_BACKGROUND_ONLY.
- GSE154600: SPP1 transcript=DETECTED_REPLICATED; program=SUPPORTED_REPLICATED; reproducibility=REPLICATED; malignancy=MALIGNANT_SUPPORTIVE_SINGLE_METHOD; receiver context=NOT_EVALUABLE_OR_NOT_ESTABLISHED.
- GSE154763: SPP1 transcript=DETECTED_REPLICATED; program=SUPPORTED_REPLICATED; reproducibility=REPLICATED; malignancy=NOT_APPLICABLE_REFERENCE; receiver context=REFERENCE_BACKGROUND_ONLY.
- GSE158722: SPP1 transcript=DETECTED_REPLICATED; program=SUPPORTED_REPLICATED; reproducibility=REPLICATED; malignancy=NOT_AVAILABLE; receiver context=NOT_EVALUABLE_OR_NOT_ESTABLISHED.

SPP1 transcript detection, companion-program support, and within-dataset relative enrichment are separate endpoints. The top-quartile SPP1-high fraction is retained only as a relative indicator and never defines transcript presence.

## GSE158722 P04_Time3

P04_Time3 contains 36 high-confidence macrophages; 80.6% detect SPP1 transcript and 52.8% meet the SPP1-plus-companion program definition. Its final status is `TRANSCRIPT_AND_PROGRAM_PRESENT`, not SPP1-negative.

## Threshold sensitivity

- GSE147082: REPLICATED; robustness=stable_replicated.
- GSE151214: REPLICATED; robustness=stable_replicated.
- GSE154600: REPLICATED; robustness=stable_replicated.
- GSE154763: REPLICATED; robustness=stable_replicated.
- GSE158722: REPLICATED; robustness=stable_replicated.

The formal analysis uses at least 20 macrophages, transcript-positive fraction at least 0.10, and program-cell fraction at least 0.10. Sensitivity spans 10/20/30 cells and 0.05/0.10/0.20 transcript and program thresholds.

## GSE147082 cluster 6

Cycling is stored only in `cell_state`. Subcluster 2 has CD3/TCR co-positive fraction 0.164 and NK-marker-positive fraction 0.986; it is frozen as `NK_like_unresolved` with `cell_state=NK_like_cytotoxic`, `patient_enriched=TRUE`, and confidence `Review`.

## Formal malignancy assessment

No reusable target-specific inferCNV or CopyKAT result was found. New patient-internal CopyKAT runs evaluated 1671 GSE154600 epithelial candidates and 819 GSE147082 PT-2834 candidates. GSE154600 has 1659 single-method malignant-supportive cells and 0 two-method high-confidence cells. GSE147082 has 0 malignant-supportive and 0 high-confidence cells.

PT-2834 cluster 4: `Mesenchymal_stromal_candidate` (CopyKAT aneuploid fraction 0).
PT-2834 cluster 7: `COL2A1_positive_chondrocyte_like_fibroblast_candidate` (CopyKAT aneuploid fraction 0).

The previous CNV-like intensity ratio is audit-only and has no role in malignancy classification. CopyKAT diploid calls are not promoted to `DIPLOID_SUPPORTIVE` because inferCNV is unavailable. Likewise, CopyKAT aneuploid calls remain `MALIGNANT_SUPPORTIVE`, not two-method `MALIGNANT_HIGH_CONFIDENCE`.

GSE158722 malignancy is `NOT_EVALUABLE`: platform identity cannot be reliably recovered, so a same-platform reference cannot be selected without fabrication.

## Malignant receiver context

Among evaluable GSE154600 malignant-supportive samples (n=5), median positive fractions are CD44=0.1, ITGB1=0.235, CD44/ITGB1 dual=0.042, and ITGB1/alpha-integrin co-positive=0.01.

CD44, ITGB1, ITGB1-alpha partner, and dual support are reported separately. ITGB1 expression alone is not evidence of a complete functional integrin receptor, and expression co-occurrence does not establish direct SPP1 binding.

## Interpretation boundary

- SPP1 transcript and companion-program support recur across multiple patients, with patient heterogeneity.
- GSE151214 and GSE154763 remain reference-only and are excluded from tumor effect aggregation.
- `tumor_specificity_status` is `NOT_ESTABLISHED` for every dataset.
- The data do not establish tumor-specific SPP1 macrophages, direct SPP1-ITGB1 binding, receptor activation, causality, or spatial contact.
- Spatial assays and/or wet-lab perturbation are required for mechanism validation.

## Frozen outputs

- `state_detection_by_sample.csv` and `state_detection_by_patient.csv`
- `state_threshold_sensitivity.csv`
- `GSE147082_cluster6_final_cell_annotation.csv` and `GSE147082_cluster6_tcr_nk_evidence.csv`
- `existing_malignancy_results_audit.csv`, `malignancy_summary_by_patient.csv`, and local per-cell consensus
- `malignant_epithelial_receiver_context.csv`
- `external_scrna_evidence_matrix_v3.csv`

External scRNA preprocessing and evidence definitions are frozen at v6; subsequent work should proceed to spatial and mechanistic validation rather than another ordinary threshold revision.
