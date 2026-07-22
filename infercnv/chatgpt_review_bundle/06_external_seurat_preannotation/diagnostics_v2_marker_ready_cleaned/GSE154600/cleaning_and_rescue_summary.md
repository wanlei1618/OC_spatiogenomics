# GSE154600 cleaning and lineage-rescue summary

- repaired_qc_cells: 31103
- retained_cells: 31101
- removed_heterotypic_doublets: 0
- removed_low_quality: 0
- removed_ambient_dominated: 2
- rescued_wrong_parent_cells: 3911
- rescued_wrong_parent_clusters: 38
- unresolved_cells: 1110
- patient_enriched_clusters: 31
- platform_confounded_clusters: 0

## Final cell-type counts

- B_cell: 748
- Endothelial: 553
- Epithelial: 4943
- Erythroid: 72
- Fibroblast: 3738
- Macrophage: 4754
- Mast: 38
- Monocyte: 638
- NK_cell: 966
- Pericyte: 11
- Plasma_cell: 1036
- T_cell: 12161
- Unresolved: 1110
- cDC1: 34
- cDC2: 282
- pDC: 17

## Patient-enriched clusters

- Endothelial__B_harmony_sample__1
- Epithelial__A_uncorrected__0
- Epithelial__A_uncorrected__1
- Epithelial__A_uncorrected__2
- Epithelial__A_uncorrected__3
- Epithelial__A_uncorrected__4
- Epithelial__A_uncorrected__5
- Epithelial__A_uncorrected__6
- Epithelial__A_uncorrected__7
- Epithelial__A_uncorrected__8
- Epithelial__A_uncorrected__9
- Fibroblast__A_uncorrected__0
- Fibroblast__A_uncorrected__1
- Fibroblast__A_uncorrected__2
- Fibroblast__A_uncorrected__3
- Macrophage__B_harmony_sample__2
- Macrophage__B_harmony_sample__5
- Macrophage__B_harmony_sample__6
- NK_cell__B_harmony_sample__0
- NK_cell__B_harmony_sample__3
- Plasma_cell__B_harmony_sample__2
- Plasma_cell__B_harmony_sample__3
- Plasma_cell__B_harmony_sample__4
- T_cell__B_harmony_sample__0
- T_cell__B_harmony_sample__3
- Unresolved__A_uncorrected__0
- Unresolved__A_uncorrected__1
- Unresolved__A_uncorrected__5
- Unresolved__A_uncorrected__6
- cDC1__broad_type_only__1
- cDC2__A_uncorrected__0

## Platform handling

- No explicit platform field in repaired-QC metadata.
- Platform labels were left NA; no platform identity was inferred.
- GSE154600 non-epithelial Harmony used sample_id only when feasible.

## Unresolved clusters

- Unresolved__A_uncorrected__0
- Unresolved__A_uncorrected__1
- Unresolved__A_uncorrected__2
- Unresolved__A_uncorrected__3
- Unresolved__A_uncorrected__4
- Unresolved__A_uncorrected__5
- Unresolved__A_uncorrected__6

All decontamination was performed within original sample. Existing sample-wise scDblFinder calls were retained and combined with cell-level incompatible-lineage evidence; no cluster was deleted wholesale.
Cycling, IFN_response, Hypoxia, Stress_response, SPP1_high, C1QC_high and FOLR2_high were stored only as cell_state.
Every exported marker passed p_val_adj < 0.05, avg_log2FC > 0.25 and pct.1 >= 0.20.
