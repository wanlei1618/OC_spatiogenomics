# GSE158722 cleaning and lineage-rescue summary

- repaired_qc_cells: 68568
- retained_cells: 68540
- removed_heterotypic_doublets: 0
- removed_low_quality: 1
- removed_ambient_dominated: 27
- rescued_wrong_parent_cells: 1846
- rescued_wrong_parent_clusters: 43
- unresolved_cells: 25403
- patient_enriched_clusters: 74
- platform_confounded_clusters: 0

## Final cell-type counts

- B_cell: 170
- Endothelial: 38
- Epithelial: 34453
- Erythroid: 6
- Fibroblast: 984
- Macrophage: 1242
- Monocyte: 3329
- NK_cell: 1041
- Plasma_cell: 11
- T_cell: 1687
- Unresolved: 25403
- cDC1: 160
- pDC: 16

## Patient-enriched clusters

- B_cell__A_uncorrected__0
- B_cell__A_uncorrected__1
- Epithelial__A_uncorrected__0
- Epithelial__A_uncorrected__1
- Epithelial__A_uncorrected__10
- Epithelial__A_uncorrected__11
- Epithelial__A_uncorrected__13
- Epithelial__A_uncorrected__14
- Epithelial__A_uncorrected__16
- Epithelial__A_uncorrected__17
- Epithelial__A_uncorrected__18
- Epithelial__A_uncorrected__19
- Epithelial__A_uncorrected__2
- Epithelial__A_uncorrected__20
- Epithelial__A_uncorrected__21
- Epithelial__A_uncorrected__22
- Epithelial__A_uncorrected__23
- Epithelial__A_uncorrected__24
- Epithelial__A_uncorrected__26
- Epithelial__A_uncorrected__27
- Epithelial__A_uncorrected__3
- Epithelial__A_uncorrected__4
- Epithelial__A_uncorrected__5
- Epithelial__A_uncorrected__6
- Epithelial__A_uncorrected__7
- Epithelial__A_uncorrected__8
- Epithelial__A_uncorrected__9
- Fibroblast__A_uncorrected__5
- Macrophage__A_uncorrected__0
- Macrophage__A_uncorrected__1
- Macrophage__A_uncorrected__4
- Monocyte__A_uncorrected__0
- Monocyte__A_uncorrected__1
- Monocyte__A_uncorrected__2
- Monocyte__A_uncorrected__3
- Monocyte__A_uncorrected__4
- Monocyte__A_uncorrected__6
- Monocyte__A_uncorrected__7
- NK_cell__A_uncorrected__0
- NK_cell__A_uncorrected__1
- NK_cell__A_uncorrected__2
- NK_cell__A_uncorrected__3
- T_cell__A_uncorrected__0
- T_cell__A_uncorrected__1
- T_cell__A_uncorrected__2
- T_cell__A_uncorrected__3
- T_cell__A_uncorrected__4
- T_cell__A_uncorrected__5
- Unresolved__A_uncorrected__0
- Unresolved__A_uncorrected__10
- Unresolved__A_uncorrected__11
- Unresolved__A_uncorrected__12
- Unresolved__A_uncorrected__13
- Unresolved__A_uncorrected__15
- Unresolved__A_uncorrected__18
- Unresolved__A_uncorrected__19
- Unresolved__A_uncorrected__2
- Unresolved__A_uncorrected__20
- Unresolved__A_uncorrected__21
- Unresolved__A_uncorrected__22
- Unresolved__A_uncorrected__23
- Unresolved__A_uncorrected__24
- Unresolved__A_uncorrected__25
- Unresolved__A_uncorrected__26
- Unresolved__A_uncorrected__27
- Unresolved__A_uncorrected__28
- Unresolved__A_uncorrected__29
- Unresolved__A_uncorrected__3
- Unresolved__A_uncorrected__31
- Unresolved__A_uncorrected__4
- Unresolved__A_uncorrected__6
- Unresolved__A_uncorrected__8
- cDC1__A_uncorrected__0
- cDC1__A_uncorrected__1

## Platform handling

- No explicit platform field in repaired-QC metadata.
- 10X/iCell8 stratification could not be verified; samples were processed separately and platform labels were left NA.
- No sample_id or timepoint Harmony was used for GSE158722.

## Unresolved clusters

- Unresolved__A_uncorrected__0
- Unresolved__A_uncorrected__1
- Unresolved__A_uncorrected__10
- Unresolved__A_uncorrected__11
- Unresolved__A_uncorrected__12
- Unresolved__A_uncorrected__13
- Unresolved__A_uncorrected__14
- Unresolved__A_uncorrected__15
- Unresolved__A_uncorrected__16
- Unresolved__A_uncorrected__17
- Unresolved__A_uncorrected__18
- Unresolved__A_uncorrected__19
- Unresolved__A_uncorrected__2
- Unresolved__A_uncorrected__20
- Unresolved__A_uncorrected__21
- Unresolved__A_uncorrected__22
- Unresolved__A_uncorrected__23
- Unresolved__A_uncorrected__24
- Unresolved__A_uncorrected__25
- Unresolved__A_uncorrected__26
- Unresolved__A_uncorrected__27
- Unresolved__A_uncorrected__28
- Unresolved__A_uncorrected__29
- Unresolved__A_uncorrected__3
- Unresolved__A_uncorrected__30
- Unresolved__A_uncorrected__31
- Unresolved__A_uncorrected__4
- Unresolved__A_uncorrected__5
- Unresolved__A_uncorrected__6
- Unresolved__A_uncorrected__7
- Unresolved__A_uncorrected__8
- Unresolved__A_uncorrected__9

All decontamination was performed within original sample. Existing sample-wise scDblFinder calls were retained and combined with cell-level incompatible-lineage evidence; no cluster was deleted wholesale.
Cycling, IFN_response, Hypoxia, Stress_response, SPP1_high, C1QC_high and FOLR2_high were stored only as cell_state.
Every exported marker passed p_val_adj < 0.05, avg_log2FC > 0.25 and pct.1 >= 0.20.
