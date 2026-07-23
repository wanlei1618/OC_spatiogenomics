# Final external scRNA freeze report v3.1

## Scope

v6.1 audits GSE154600 CopyKAT target coverage, repeats same-patient immune-reference sampling with seeds 20260718/20260719/20260720, builds a stable single-method malignancy layer, recalculates receiver expression, and corrects GSE154763 reference terminology. v6 outputs remain unchanged; no full QC, clustering, broad annotation, or SPP1 threshold was rerun.

## CopyKAT coverage and stability

- T59: final=1513; available/submitted=1422/1422 (coverage 0.940); any defined=860 (rate 0.605); stable aneuploid=843; stable diploid=1; unstable=15; mostly not defined=563; not submitted=91.
- T76: final=125; available/submitted=37/37 (coverage 0.296); any defined=21 (rate 0.568); stable aneuploid=17; stable diploid=1; unstable=3; mostly not defined=16; not submitted=88.
- T77: final=1175; available/submitted=1023/1023 (coverage 0.871); any defined=563 (rate 0.550); stable aneuploid=535; stable diploid=3; unstable=22; mostly not defined=463; not submitted=152.
- T89: final=181; available/submitted=173/173 (coverage 0.956); any defined=49 (rate 0.283); stable aneuploid=45; stable diploid=0; unstable=2; mostly not defined=126; not submitted=8.
- T90: final=975; available/submitted=803/803 (coverage 0.824); any defined=263 (rate 0.328); stable aneuploid=256; stable diploid=1; unstable=0; mostly not defined=546; not submitted=172.

Across GSE154600, 3969 final epithelial cells were audited; 3458 were found in at least one of the four specified lineage count inputs and submitted; 1756 received a defined call at least once; 1696 were stable aneuploid; 6 were stable diploid-like; 42 were discordant; 1714 were mostly not defined; and 511 were not submitted because counts were absent.

## T76 coverage correction

v6 submitted 37 of 125 T76 final epithelial cells. Cross-lineage collection still finds 37 of 125 (coverage 0.296); the remaining 88 are absent from all four task-specified lineage count matrices. Coverage therefore did not numerically increase, but the omission is now explicit and classified as `NOT_SUBMITTED`, rather than silently excluded.

## Three-seed reference stability

- T59: aneuploid 850-853; diploid 4-9; not defined 562-565.
- T76: aneuploid 17-20; diploid 1-4; not defined 16-16.
- T77: aneuploid 540-544; diploid 7-17; not defined 465-476.
- T89: aneuploid 46-47; diploid 0-2; not defined 125-127.
- T90: aneuploid 254-262; diploid 1-1; not defined 540-548.

Reference resampling produced 1696 stable aneuploid and 6 stable diploid-like calls, with 42 discordant cells. Seed-to-seed totals varied modestly but did not overturn the dominant patient-level pattern.

Stable CopyKAT aneuploid is `MALIGNANT_SUPPORTIVE_STABLE`, which remains single-method supportive evidence. It is not double-method high-confidence malignancy because inferCNV is unavailable.

## Stable malignant receiver

- T59: n=843; CD44=0.163 (SUPPORTED); ITGB1=0.361 (SUPPORTED); ITGB1/any-alpha=0.058 (DETECTED_LOW); dual=0.083 (DETECTED_LOW); dominant alpha=ITGAV.
- T76: n=17; CD44=0.059 (NOT_EVALUABLE); ITGB1=0.235 (NOT_EVALUABLE); ITGB1/any-alpha=0.000 (NOT_EVALUABLE); dual=0.059 (NOT_EVALUABLE); dominant alpha=ITGAV.
- T77: n=535; CD44=0.043 (DETECTED_LOW); ITGB1=0.135 (SUPPORTED); ITGB1/any-alpha=0.009 (DETECTED_LOW); dual=0.006 (DETECTED_LOW); dominant alpha=ITGAV.
- T89: n=45; CD44=0.000 (NOT_DETECTED); ITGB1=0.089 (DETECTED_LOW); ITGB1/any-alpha=0.000 (NOT_DETECTED); dual=0.000 (NOT_DETECTED); dominant alpha=ITGAV.
- T90: n=256; CD44=0.027 (DETECTED_LOW); ITGB1=0.094 (DETECTED_LOW); ITGB1/any-alpha=0.000 (NOT_DETECTED); dual=0.004 (DETECTED_LOW); dominant alpha=ITGA4.

Across evaluable patients, stable malignant-supportive receiver status is CD44=SUPPORTED, ITGB1=SUPPORTED, ITGB1-alpha partner=DETECTED_LOW, and CD44/ITGB1 dual=DETECTED_LOW.

ITGB1 expression alone does not establish a complete functional receptor, and co-expression does not establish direct SPP1 binding.

## GSE147082 retained interpretation

- Cluster 4 remains `Mesenchymal_stromal_candidate`.
- Cluster 7 remains `COL2A1_positive_chondrocyte_like_fibroblast_candidate`.

The prior single CopyKAT diploid result is retained as single-method evidence and does not over-upgrade or force a malignant/normal label.

## GSE154763 terminology correction

- SPP1: true transcript detection plus separate reference program support; expression columns available=TRUE.
- C1QC: true core transcript detection plus separate reference program support; expression columns available=TRUE.
- FOLR2: true transcript detection plus separate reference program support; expression columns available=TRUE.

SPP1, C1QC-core, and FOLR2 transcript detection now use their explicit normalized expression columns (>0). Reference program support is reported separately. Threshold robustness is termed `reference_consistent`, `reference_threshold_sensitive`, or `reference_not_evaluable`, rather than raw-count `stable_replicated`.

## Frozen evidence

- SPP1/C1QC/FOLR2 v6 raw-count thresholds and conclusions remain frozen.
- GSE154600 CopyKAT coverage, defined-call rate, three-seed stability, and stable receiver context are frozen at v3.1.
- GSE151214 remains normal-reference only; GSE154763 remains author-normalized myeloid-reference only.
- `tumor_specificity_status` remains `NOT_ESTABLISHED` for every dataset.

## Remaining validation

inferCNV or another independent formal CNV method is required for double-method high-confidence malignancy. Spatial assays are required to establish sender-receiver proximity. Wet-lab assays are required to establish direct binding, receptor activation, and causality. The available data do not establish tumor-specific SPP1 macrophages or SPP1-driven progression.

External scRNA evidence is frozen at v6.1; subsequent work should move to independent CNV, spatial, and mechanistic validation rather than another ordinary threshold revision.
