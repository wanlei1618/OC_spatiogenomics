# First-stage clone validity report

Generated: 2026-07-11 12:08:00

## Input
- Checked metadata: D:\OC_spatiogenomics\infercnv\06_rigorous_validation/00_config/sample_metadata_checked.csv
- Clone-labeled rows: 6507

## Field check
- Missing fields: patient_id, dataset, doublet_score, S.Score, G2M.Score
- No patient_id was inferred from sample_id; cross-patient validity remains blocked until an explicit patient mapping is supplied.

## Preliminary Gate A
- patient_id is missing in the selected metadata, so current CNV_Subclone_01-05 labels cannot be called cross-patient CNV states from this dataset alone.
- sample_id-level composition was computed as a diagnostic, not as a substitute for patient-level biological replication.

## Key outputs
- 01_clone_patient_confounding/clone_by_sample_counts.csv
- 01_clone_patient_confounding/clone_by_sample_type_counts.csv
- 01_clone_patient_confounding/clone_by_batch_counts.csv
- 01_clone_patient_confounding/clone_patient_confounding_metrics.csv
- 01_clone_patient_confounding/clone_association_cramers_v.csv
- figures/clone_by_sample_heatmap.png
- figures/clone_by_sample_mosaic.png
- figures/clone_by_sample_type_alluvial_like.png

## Next required action
Provide or construct an explicit patient_id mapping for sample_id before making cross-patient clone claims or proceeding to old-clone mechanism analyses.
