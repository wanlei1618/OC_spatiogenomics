# GSE158722 mitochondrial-QC decision

- Prepared features: 18392
- Mitochondrial features used: 13
- Feature source: original_raw_gene_symbol_per_patient
- `percent.mt` available: TRUE
- Fraction of cells with available `percent.mt`: 0.98886
- Ribosomal features audited: 105
- Hemoglobin features audited: 10

Per-patient original raw count files were loaded independently, and their explicit `^MT-` features were passed to `PercentageFeatureSet(features = ...)`. Only the resulting QC metadata was joined by cell ID. The prepared common-gene expression matrix was not changed and no features were added. Patients whose raw source lacks credible mitochondrial features retain `percent.mt = NA` and do not use an mt threshold; availability is recorded per patient.

Raw source audit: 260 `^MT-` matches and 0 `^MT.` matches across 23 readable source files.

The existing stage 06 result was not overwritten. Repaired metadata is stored only under `diagnostics_v2`.
