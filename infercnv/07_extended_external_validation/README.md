# Extended external scRNA annotation upgrade

This directory separates external scRNA validation by dataset role and annotation provenance. The earlier marker-score-only outputs are archived under `tables/archive_marker_score_v1/` and are no longer treated as author-original annotation results.

## Dataset Roles

- `GSE154600`: tumor ecosystem, author SCE annotations, primary only when sample identity and annotation audit pass.
- `GSE158722`: malignant-fluid tumor ecosystem, author cell annotations, primary only when barcode/composite-key audit passes.
- `GSE147082`: tumor ecosystem, secondary reannotation only, sensitivity analysis only.
- `GSE154763`: myeloid reference only; no tumor target and no cohort tumor LR analysis.
- `GSE151214`: normal fallopian tube reference; no ovarian tumor TME meta-analysis.

## Key Outputs

- `config/external_scrna_dataset_registry.csv`
- `annotations/manifests/original_annotation_download_manifest.csv`
- `annotations/harmonized/*annotations.csv.gz`
- `tables/external_scrna_annotation_match_audit.csv`
- `tables/external_scrna_dataset_reclassification.csv`
- `tables/curated_external_scRNA_primary_meta_summary.csv`
- `tables/curated_external_scRNA_sensitivity_meta_summary.csv`
- `reports/external_scrna_annotation_upgrade_report.md`

## Interpretation Boundaries

SPP1-CD44 may be treated as a candidate ligand-receptor axis. SPP1-ITGB1 should be described as an SPP1-associated ITGB1-positive adhesion/integrin program unless independently confirmed by an authoritative LR resource. Expression products are reported as expression potential, not complete cell-cell communication inference.

## Legacy Marker-Score Outputs

Legacy files such as `extended_external_scRNA_dataset_level_summary.csv`, `extended_external_scRNA_meta_summary.csv`, `extended_external_scRNA_celltype_composition.csv`, and `extended_external_scRNA_status.csv` are retained only for comparison after archival.
