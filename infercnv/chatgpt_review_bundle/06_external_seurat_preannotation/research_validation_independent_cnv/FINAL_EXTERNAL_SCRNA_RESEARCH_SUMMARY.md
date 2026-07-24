# Final external scRNA research transition summary

## 1. Recovery of 511 final epithelial cells

- T59: traced=91; recovered=91; found in full object=91; found in GEO feature-barcode matrix=91.
- T76: traced=88; recovered=88; found in full object=88; found in GEO feature-barcode matrix=88.
- T77: traced=152; recovered=152; found in full object=152; found in GEO feature-barcode matrix=152.
- T89: traced=8; recovered=8; found in full object=8; found in GEO feature-barcode matrix=8.
- T90: traced=172; recovered=172; found in full object=172; found in GEO feature-barcode matrix=172.

All 511 previously unsubmitted final epithelial cells were recovered from the authoritative GSE154600 preannotation Assay5 raw-count layer. Their GEO feature-barcode source was also traced. The prior omission was caused by incomplete lineage-strategy matrices, not by absent raw counts.

## 2. Final count coverage

- T59: 1513/1513 (1); recovered=91.
- T76: 125/125 (1); recovered=88.
- T77: 1175/1175 (1); recovered=152.
- T89: 181/181 (1); recovered=8.
- T90: 975/975 (1); recovered=172.

Final GSE154600 epithelial count coverage is 3969/3969 (1). T76 is 125/125 (1).

## 3. CopyKAT defined-cell technical bias

Patient-stratified descriptive audit status: `HIGHER_DEPTH_IN_DEFINED_CALLS`. Cell-level tests are treated only as technical diagnostics, not as patient-level biological replication.

## 4. Independent CNV support

- T59: high CNV=1111; stable CopyKAT aneuploid=870; dual-method=799; CopyKAT-only=71; inferCNV-only=312; concordance among stable CopyKAT=0.918.
- T76: high CNV=87; stable CopyKAT aneuploid=58; dual-method=38; CopyKAT-only=20; inferCNV-only=49; concordance among stable CopyKAT=0.655.
- T77: high CNV=227; stable CopyKAT aneuploid=610; dual-method=118; CopyKAT-only=492; inferCNV-only=109; concordance among stable CopyKAT=0.193.
- T89: high CNV=119; stable CopyKAT aneuploid=44; dual-method=33; CopyKAT-only=11; inferCNV-only=86; concordance among stable CopyKAT=0.750.
- T90: high CNV=289; stable CopyKAT aneuploid=243; dual-method=162; CopyKAT-only=81; inferCNV-only=127; concordance among stable CopyKAT=0.667.

Across patients, inferCNV-high cells=1833; dual-method supportive cells=1150.

The standard inferCNV R package was unavailable because the local Bioconductor installation could not validate through the Windows security channel. The task-authorized infercnvpy fallback was therefore used for continuous CNV signal validation. No HMM subclone calls and no `confirmed malignant` wording are used.

## 5. Patient-replicated receiver expression

- T59: tier=DUAL_METHOD_MALIGNANT_SUPPORT; n=799; CD44=0.217 (SUPPORTED); ITGB1=0.538 (SUPPORTED); ITGB1-alpha=0.088 (DETECTED_LOW); CD44/ITGB1=0.134 (SUPPORTED).
- T76: tier=DUAL_METHOD_MALIGNANT_SUPPORT; n=38; CD44=0.184 (SUPPORTED); ITGB1=0.500 (SUPPORTED); ITGB1-alpha=0.079 (DETECTED_LOW); CD44/ITGB1=0.158 (SUPPORTED).
- T77: tier=DUAL_METHOD_MALIGNANT_SUPPORT; n=118; CD44=0.034 (DETECTED_LOW); ITGB1=0.246 (SUPPORTED); ITGB1-alpha=0.017 (DETECTED_LOW); CD44/ITGB1=0.017 (DETECTED_LOW).
- T89: tier=DUAL_METHOD_MALIGNANT_SUPPORT; n=33; CD44=0.152 (SUPPORTED); ITGB1=0.121 (SUPPORTED); ITGB1-alpha=0.000 (NOT_DETECTED); CD44/ITGB1=0.061 (DETECTED_LOW).
- T90: tier=DUAL_METHOD_MALIGNANT_SUPPORT; n=162; CD44=0.031 (DETECTED_LOW); ITGB1=0.148 (SUPPORTED); ITGB1-alpha=0.000 (NOT_DETECTED); CD44/ITGB1=0.006 (DETECTED_LOW).

- CD44: REPLICATED_SUPPORTED; evaluable patients=5; supported=3; low=2; not detected=0.
- ITGB1: REPLICATED_SUPPORTED; evaluable patients=5; supported=5; low=0; not detected=0.
- ITGB1_alpha: REPLICATED_LOW; evaluable patients=5; supported=0; low=3; not detected=2.
- dual_CD44_ITGB1: REPLICATED_SUPPORTED; evaluable patients=5; supported=2; low=3; not detected=0.

`malignant_receiver_expression_context` is `DESCRIPTIVE_ONLY`: there is no sample-level association statistic, spatial proximity evidence, direct ligand-receptor evidence, or causal evidence.

## 6. Research transition

GSE154600 has `spatial_validation_priority = HIGH`. The external scRNA evidence is sufficient to prioritize a spatial validation study of SPP1 macrophages and CNV-supported epithelial ITGB1/CD44 context, but it does not itself establish spatial interaction, direct receptor engagement, tumor specificity, or causality.

External scRNA preprocessing stops here: no v6.2/v6.3 threshold versions, additional CopyKAT seeds, SPP1 threshold changes, or new cell-type cleanup are introduced.
