# GSE154600 marker-ready annotation summary

- Repaired-QC eligible cells: 31103
- Balanced discovery cells used: 30000
- Global broad types were assigned from adjusted-P significant RNA markers.
- Cycling is exported only as a state, never as a broad lineage.
- Epithelial subclustering uses A_uncorrected.
- Non-epithelial Harmony, when feasible, uses sample_id.
- All annotation-ready marker rows satisfy p_val_adj < 0.05, avg_log2FC > 0.25, and pct.1 >= 0.20.
- Annotation clusters exported: 58
- Clusters READY_FOR_ANNOTATION: 54
- Clusters needing marker review: 4
- Lineage-conflict clusters: 23

Primary file for manual annotation:
`annotation_ready_cluster_template.csv`

Fill only these columns:
`manual_cell_type`, `manual_cell_subtype`, `manual_confidence`, `manual_notes`.
