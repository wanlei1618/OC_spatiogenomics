# GSE158722 marker-ready annotation summary

- Repaired-QC eligible cells: 68568
- Balanced discovery cells used: 30000
- Global broad types were assigned from adjusted-P significant RNA markers.
- Cycling is exported only as a state, never as a broad lineage.
- Epithelial subclustering uses A_uncorrected.
- Non-epithelial Harmony, when feasible, uses patient_id only; sample_id/timepoint is not corrected.
- All annotation-ready marker rows satisfy p_val_adj < 0.05, avg_log2FC > 0.25, and pct.1 >= 0.20.
- Annotation clusters exported: 48
- Clusters READY_FOR_ANNOTATION: 48
- Clusters needing marker review: 0
- Lineage-conflict clusters: 14

Primary file for manual annotation:
`annotation_ready_cluster_template.csv`

Fill only these columns:
`manual_cell_type`, `manual_cell_subtype`, `manual_confidence`, `manual_notes`.
