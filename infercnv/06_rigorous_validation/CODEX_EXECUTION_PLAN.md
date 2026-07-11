# CODEX EXECUTION PLAN: 06 rigorous validation

This plan follows the supplied task book and limits the first execution phase to P0 entry analyses.

## Phase 1 scope

1. Scan existing project files and build an input manifest.
2. Check whether available metadata contains patient_id, sample_id, sample_type, batch, dataset, clone and QC fields.
3. Generate a missing-field and potential-mapping report without guessing patient mappings.
4. Create the 06_rigorous_validation directory structure and project_config.yaml.
5. Implement and run:
   - scripts/R/01_build_checked_metadata.R
   - scripts/R/02_clone_patient_confounding.R
6. Write a first-stage clone validity report and a Gate A preliminary decision.

## Guardrails

- Large outputs and temporary files stay under D:/OC_spatiogenomics.
- The main biological replicate is patient_id when available, otherwise sample_id can be reported only as sample-level evidence.
- Cell-level counts are treated as observations/composition summaries, not as independent biological replicates.
- Missing patient_id or dataset fields are reported as blocking limitations for cross-patient claims.
- P1/P2 analyses, NicheNet, spatial deconvolution and virtual perturbation are not run in this phase.
