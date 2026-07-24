# Corrected spatial validation report

## 1. Sample audit and GSM6177618 exclusion

Only GSM6177614 and GSM6177617 are coordinate-aware ovarian sections. GSM6177618 is verified PDAC and is excluded from coordinate analysis, expression analysis, pooled ovarian conclusions, figures and evidence matrices (`PDAC_NOT_OVARIAN`).

## 2. GSE203612 SPP1 to ITGB1

GSM6177614 OE=0.618, empirical p=1.0000; GSM6177617 OE=1.206, empirical p=0.0080. Dataset conclusion: **SPATIALLY_HETEROGENEOUS**.

## 3. GSE203612 SPP1 to CD44

GSM6177614 OE=1.026, empirical p=0.3956; GSM6177617 OE=1.273, empirical p=0.0010. Dataset conclusion: **LIMITED_SPATIAL_SUPPORT**.

## 4. C1QC and general macrophage controls

C1QC to ITGB1: GSM6177614 OE=0.535, p=1.0000; GSM6177617 OE=1.171, p=0.0230

C1QC to CD44: GSM6177614 OE=0.833, p=0.9650; GSM6177617 OE=1.048, p=0.2557

General macrophage to epithelial: GSM6177614 OE=0.765, p=1.0000; GSM6177617 OE=1.000, p=1.0000

## 5. GSE189843 expression-only evidence

12 author-included samples with available matrices were retained (6 Excellent, 6 Poor). Coordinates remain unverified; these results are sample/patient-level expression summaries and must not be called spatial colocalization, proximity or neighborhood enrichment.

- median_SPP1_program: poor-minus-excellent median -0.0317, Cliff delta -0.444, exact p 0.2403
- median_ITGB1: poor-minus-excellent median 0.0306, Cliff delta 0.056, exact p 0.8939
- median_CD44: poor-minus-excellent median -0.5740, Cliff delta -0.611, exact p 0.0606
- SPP1_ITGB1_sample_internal_spearman: poor-minus-excellent median -0.0336, Cliff delta -0.222, exact p 0.5887
- SPP1_CD44_sample_internal_spearman: poor-minus-excellent median -0.0232, Cliff delta -0.278, exact p 0.4848

## 6. Optional GSE211956 replication

**NOT_RUN_INPUT_INCOMPLETE**: complete matrices and coordinates are present, but an authoritative sample-level HG-SOC identity record is absent from the frozen registry. No new download or coordinate analysis was performed.

## 7. Wet-lab priority

Priority: **CD44**, based on the corrected two-section specificity comparison rather than expression-product maps.

## 8. Interpretation boundary

The analysis can describe enrichment, depletion or heterogeneity of SPP1 macrophage programs near ITGB1/CD44-positive CNV-supported epithelial expression regions. It does not prove direct SPP1 binding, receptor activation, causal tumor progression or tumor-specificity of SPP1 macrophages. No output is represented as COMMOT or optimal transport.
