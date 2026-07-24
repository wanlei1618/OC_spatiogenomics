# Final external scRNA validation report

## Scope and statistical boundary

No whole-dataset reclustering or QC was repeated. All cross-dataset comparisons use within-dataset percentiles, positive fractions and patient/sample-level prevalence; raw module scores from different methods are never compared. GSE151214 and GSE154763 are references only and are excluded from tumor-effect aggregation.

## 1. GSE147082 targeted annotation

- Cluster 4: **Mesenchymal_stromal**. CNV-equivalent intensity=0.0059, normal-stroma similarity=0.499, tumor-epithelial similarity=0.268.
- Cluster 7: **COL2A1_positive_chondrocyte_like_fibroblast** with 10 cartilage-program markers; its CNV-equivalent profile is more similar to normal stroma than tumor epithelium.
- Cluster 6: targeted subclustering yielded Cytotoxic_T n=298; Cytotoxic_T n=99; Gamma_delta_T n=73. It is no longer treated as ordinary T cells.
- Cluster 19 remains Unresolved and is excluded from downstream proportion, state, LR and differential analyses.

## 2. GSE151214 normal myeloid reference

Cluster 18 split into C1QC_macrophage n=420; CD1C_CLEC10A_dendritic_like n=151; Macrophage_DC_mixed n=90. Macrophage_DC_mixed cells are not counted as pure macrophage background. Groups with fewer than 20 cells remain audited but are excluded from quantitative summaries.
Normal C1QC macrophage SPP1 background across evaluable sample groups: median average expression=0.8157, median positive fraction=0.270, n=5 samples.

## 3. GSE154763 author-driven state validation

Strict lineage-specific thresholds assigned None to 1658 cells and Mixed to 787 cells; DC/Mast cells remain None. Author subtype remains primary.
Macro_SPP1 median SPP1 program=0.277 versus Macro_C1QC=0.006; Macro_C1QC median C1QC/FOLR2 programs=0.585/0.430. Cohort-level ordering supports the author labels, while strict cell-level concordance is partial.
The P20190304 tumor/normal comparison is descriptive only; no inference is made from one paired patient.

## 4. Macrophage-state reproducibility

- GSE154600: SPP1 NOT_SUPPORTED (0/5 evaluable patients; none); C1QC REPLICATED (3/5 evaluable patients; T77;T89;T90); FOLR2 REPLICATED (3/5 evaluable patients; T77;T89;T90).
- GSE158722: SPP1 NOT_SUPPORTED (0/2 evaluable patients; none); C1QC NOT_SUPPORTED (0/2 evaluable patients; none); FOLR2 NOT_SUPPORTED (0/2 evaluable patients; none).
- GSE147082: SPP1 SUPPORTIVE (1/4 evaluable patients; PT-3232); C1QC REPLICATED (2/4 evaluable patients; PT-3232;PT-3401); FOLR2 SUPPORTIVE (1/4 evaluable patients; PT-3401).
C1QC/FOLR2 frequently co-vary, whereas SPP1 shows dataset- and patient-dependent overlap, consistent with partially separated states along a continuum rather than mutually exclusive universal classes.

## 5. SPP1-associated CD44/ITGB1-positive adhesion context

- GSE154600: 5 samples; descriptive only; n < 6.
- GSE147082: 4 samples; descriptive only; n < 6.
- GSE158722: not evaluable because no matched high-confidence sender/receiver sample set passed both cell-count gates.
These are sample-level co-occurrence descriptions. They are not evidence of a proven direct SPP1-ITGB1 ligand-receptor interaction or causal regulation.

## 6. Evidence roles

GSE154600 is the primary tumor support dataset; GSE147082 provides sensitivity support; GSE158722 is inconclusive under strict high-confidence gates; GSE151214 is a normal reference; GSE154763 is an author-annotated myeloid reference. Only the three tumor datasets are eligible for tumor conclusions, and eligibility still depends on their explicit state/context gates.

## 7. Remaining limitations

The analyses cannot prove direct receptor binding, directionality or causality; cannot treat cells as independent replicates; cannot infer GSE154763 raw-count QC/CNV; cannot recover GSE158722 unresolved cells into the validation; and cannot claim correlation significance when fewer than six evaluable samples are available.
