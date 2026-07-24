# GSE147082 repaired-mitochondrial-QC reanalysis

- repaired_qc_input_cells: 6993
- final_retained_cells: 6993
- removed_heterotypic_doublets: 0
- removed_low_quality: 0
- clustering: newly recomputed from raw counts at resolutions 0.4, 0.6, and 0.8; resolution 0.6 used for the primary result.
- old_7645_cell_clusters: comparison only; neither old marker calls nor old cluster IDs were used to assign final labels.
- patient_enriched_clusters: 8
- unresolved_cells: 959

## Final cell-type counts

- B_cell: 352
- Endothelial: 144
- Epithelial: 1885
- Fibroblast: 1668
- Macrophage: 595
- Mast: 33
- Pericyte: 135
- Plasma_cell: 201
- T_cell: 976
- Unresolved: 959
- pDC: 45

## SPP1 macrophage reference

- cluster 2: n_cells=595, module_mean=0.5358, SPP1_positive_fraction=0.3647, positive_samples=PT-3232;PT-5150;PT-6885;PT-3401;PT-2834, cross_sample_replicated=TRUE, SPP1_high=TRUE

Patient/sample enrichment was recorded but was never a deletion criterion.
Cycling, IFN, hypoxia, stress, SPP1, C1QC and FOLR2 programs are stored only in cell_state.
