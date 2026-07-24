# Final calibrated CNV and spatial pilot decision

## 1. Held-out immune negative-control FPR

- T59: median FPR 0.0082 (3/3 splits <=0.05), INFERCNV_THRESHOLD_CALIBRATED
- T76: median FPR 0.0062 (3/3 splits <=0.05), INFERCNV_THRESHOLD_CALIBRATED
- T77: median FPR 0.0116 (3/3 splits <=0.05), INFERCNV_THRESHOLD_CALIBRATED
- T89: median FPR 0.0076 (3/3 splits <=0.05), INFERCNV_THRESHOLD_CALIBRATED
- T90: median FPR 0.0050 (3/3 splits <=0.05), INFERCNV_THRESHOLD_CALIBRATED

## 2. T77 CNV method discordance

- Epithelial__A_uncorrected__4: 245 CopyKAT-only cells (54.8% of cluster; 45.5% of T77 CopyKAT-only)
- Epithelial__A_uncorrected__5: 166 CopyKAT-only cells (51.9% of cluster; 30.8% of T77 CopyKAT-only)
- Epithelial__A_uncorrected__6: 99 CopyKAT-only cells (43.0% of cluster; 18.4% of T77 CopyKAT-only)
- Epithelial__A_uncorrected__7: 29 CopyKAT-only cells (16.4% of cluster; 5.4% of T77 CopyKAT-only)
- Epithelial__A_uncorrected__9: 0 CopyKAT-only cells (0.0% of cluster; 0.0% of T77 CopyKAT-only)

The cluster table separates sequencing-depth and epithelial-marker summaries from method calls; it does not reinterpret fixed clusters. Clusters 4, 5 and 6 account for 94.6% of T77 CopyKAT-only cells. Clusters 4 and 5 retain higher median depth and epithelial scores than clusters 6 and 7, so the discordance is cluster-associated and cannot be attributed uniformly to low depth.

## 3. Receiver pseudobulk and depth strata

- CD44: REPLICATED_ROBUST
- CD44_ITGB1_dual: REPLICATED_ROBUST
- ITGB1: REPLICATED_ROBUST
- ITGB1_any_alpha: SINGLE_PATIENT_SUPPORT

The primary tier uses only CALIBRATED_DUAL_METHOD_SUPPORT cells. CopyKAT-stable cells are reported separately as sensitivity evidence.

## 4. Spatial datasets entering the pilot

- GSE203612/GSM6177614
- GSE203612/GSM6177617
- GSE211956/GSM6506110
- GSE211956/GSM6506111
- GSE211956/GSM6506112

## 5. SPP1 macrophage proximity to ITGB1 receiver

Status: **REPLICATED_SPATIAL_SUPPORT**. GSE203612/GSM6177614 OE=0.646, p=1.0000; GSE203612/GSM6177617 OE=1.453, p=0.0110; GSE211956/GSM6506110 OE=1.377, p=0.0010; GSE211956/GSM6506111 OE=0.494, p=1.0000; GSE211956/GSM6506112 OE=1.862, p=0.0010

## 6. SPP1 macrophage proximity to CD44 receiver

Status: **REPLICATED_SPATIAL_SUPPORT**. GSE203612/GSM6177614 OE=0.854, p=0.8921; GSE203612/GSM6177617 OE=1.410, p=0.0250; GSE211956/GSM6506110 OE=1.689, p=0.0010; GSE211956/GSM6506111 OE=0.458, p=1.0000; GSE211956/GSM6506112 OE=2.060, p=0.0010

## 7. Comparison with C1QC macrophage control

GSE203612/GSM6177614 OE=0.503, p=1.0000; GSE203612/GSM6177617 OE=1.660, p=0.0010; GSE211956/GSM6506110 OE=1.018, p=0.3826; GSE211956/GSM6506111 OE=0.557, p=1.0000; GSE211956/GSM6506112 OE=1.527, p=0.0010

SPP1-to-ITGB1 is stronger than the matched C1QC control in two of the three SPP1-supporting samples, but C1QC is stronger in GSM6177617. The decision status explicitly downgrades an SPP1 result if it is not stronger than the matched C1QC or general-macrophage control.

## 8. Wet-lab decision

Proceed to a focused wet-lab follow-up: **yes**. Priority receiver: **ITGB1**.

## 9. Interpretation boundary

The supported wording is an SPP1-macrophage spatial proximity to a CNV-supported epithelial ITGB1/CD44 expression context. The analysis does not prove direct SPP1 binding, receptor activation, causal tumor progression, or tumor specificity of SPP1 macrophages.

## 10. Remaining limitations

Visium spots are mixtures, inferred signature scores are not cell identity calls, samples rather than spots are biological replicates, and the pilot contains only the locally available technically usable datasets. Spatial association requires orthogonal validation.
