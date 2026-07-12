# External scRNA original annotation upgrade report

## Scope

This upgrade reclassifies the five external scRNA datasets by biological role and annotation provenance. It archives the marker-score v1 summaries and separates primary tumor-ecosystem, sensitivity, myeloid-only, and normal-reference evidence layers.

## Dataset reclassification

| dataset_id | existing_repo_id | biological_role | annotation_status | analysis_role | primary_tumor_lr | sensitivity_tumor_lr | myeloid_source_reference | normal_reference | main_action |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE147082 | OV_GSE147082 | tumor_ecosystem | secondary_reannotation_only | sensitivity_secondary | no | yes | yes | no | sensitivity only; audit cell-count difference |
| GSE151214 | OV_GSE151214 | normal_fallopian_tube_reference | original_barcode_annotation_not_public | normal_reference | no | no | no | yes | remove from ovarian tumor TME meta |
| GSE154600 | OV_GSE154600 | tumor_ecosystem | author_original_SCE | tumor_ecosystem_primary | yes_if_audit_passes | yes | yes | no | preserve author labels; resolve T61/T77 |
| GSE154763 | OV_GSE154763 | myeloid_reference_only | author_original_GEO_metadata | myeloid_reference_only | no | no | yes | no | source-state validation only |
| GSE158722 | OV_GSE158722 | tumor_ecosystem_malignant_fluid | author_original_GEO_annotations | tumor_ecosystem_primary | yes_if_audit_passes | yes | yes | no | patient/timepoint/barcode matching |

## Download manifest and SHA-256

| dataset_id | file_id | path | bytes | sha256 | status |
| --- | --- | --- | --- | --- | --- |
| GSE154763 | GSE154763_OV-FTC_metadata.csv.gz | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154763\GSE154763_OV-FTC_metadata.csv.gz | 152616 | 7623fb51dead57d695951d32ded0c5a7cd7d63c26a4489800bca93c2f3b70d20 | downloaded |
| GSE154763 | GSE154763_OV-FTC_normalized_expression.csv.gz | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154763\GSE154763_OV-FTC_normalized_expression.csv.gz | 16563141 | 77864de964af093f1dfa1831bd3477a4f9190fdecafb85c5e87b4b43c47faafa | downloaded |
| GSE158722 | GSE158722.cell_annotations.txt.gz | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE158722\GSE158722.cell_annotations.txt.gz | 1308677 | f640e9278be56982b7625f11e5bad2b427319ea11c713acd9f28358089bb5051 | downloaded |
| GSE154600 | sample59_sce.rds | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154600\sample59_sce.rds | 113345251 | e055a319e1ba03dcf2b162ba23c92f24b219c6e9e5357a070fc88e26bb922734 | downloaded |
| GSE154600 | sample76_sce.rds | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154600\sample76_sce.rds | 99310448 | 48a7ec84c9016687daa3e884c871f9e0417b2142d71d9224cd09cfc82635c734 | downloaded |
| GSE154600 | sample77_sce.rds | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154600\sample77_sce.rds | 71853370 | c6c4198cfa08277f949b2257ae4286d725fd1cb19d9eee046e115b84b921e9a0 | downloaded |
| GSE154600 | sample89_sce.rds | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154600\sample89_sce.rds | 36849742 | c4fdd67a461dc8f7f2d9c2670a6caa1d3af6e9efad8f9173afad7b8a385786c4 | downloaded |
| GSE154600 | sample90_sce.rds | D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154600\sample90_sce.rds | 28452145 | dd0737d26f27761f9d4e93d412940e3b58994088eaf6a5ba8005564462926d1f | downloaded |

## Barcode/composite-key audit

| dataset_id | n_expression_cells | n_annotation_cells | n_exact_matches | n_normalized_matches | n_unmatched_expression | n_unmatched_annotation | match_fraction_expression | match_fraction_annotation | duplicate_raw_barcode | duplicate_normalized_barcode | duplicate_cell_key | selected_match_method | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE147082 | 9796 | 9796 | 9796 | 9796 | 0 | 0 | 1.0 | 1.0 | 0 | 0 | 0 | expression_h5_barcode_inventory | pass |
| GSE151214 | 59446 | 59446 | 59446 | 59446 | 0 | 0 | 1.0 | 1.0 | 0 | 0 | 0 | expression_h5_barcode_inventory | reference_only |
| GSE154600 | 42253 | 42253 | 42253 | 42253 | 0 | 0 | 1.0 | 1.0 | 0 | 0 | 0 | author_sce_self_contained | pass |
| GSE154763 | 3888 | 3888 | 0 | 0 | 3888 | 3888 | 0.0 | 0.0 | 1 | 0 | 0 | no_match | reference_only |
| GSE158722 | 96846 | 63793 | 0 | 0 | 96846 | 63793 | 0.0 | 0.0 | 0 | 0 | 0 | no_match | fail |

## Annotation status

| dataset_id | analysis_role | annotation_status | n_annotation_cells | primary_tumor_lr | sensitivity_tumor_lr | myeloid_source_reference | normal_reference | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE147082 | sensitivity_secondary | secondary_reannotation_only | 9796 | no | yes | yes | no | available |
| GSE151214 | normal_reference | original_barcode_annotation_not_public | 59446 | no | no | no | yes | available |
| GSE154600 | tumor_ecosystem_primary | author_original_SCE | 42253 | yes_if_audit_passes | yes | yes | no | available |
| GSE154763 | myeloid_reference_only | author_original_GEO_metadata | 3888 | no | no | yes | no | available |
| GSE158722 | tumor_ecosystem_primary | author_original_GEO_annotations | 63793 | yes_if_audit_passes | yes | yes | no | available |

## T61/T77 boundary

GSE154600 T77 is not assumed to equal GEO T61. If optional raw GSE154600 matrices are not downloaded, T77 remains unresolved and is excluded from primary-analysis eligibility.

## Interpretation boundaries

- GSE154763 is myeloid_reference_only and does not create tumor target cells or cohort LR scores.
- GSE151214 is a normal fallopian tube reference and does not enter ovarian tumor TME meta-analysis.
- GSE147082 is secondary_reannotation only and is restricted to sensitivity analysis.
- Expression potential is not named as a complete ligand-receptor inference.
- SPP1-ITGB1 is described as an SPP1-associated ITGB1-positive adhesion/integrin program unless an authoritative LR resource confirms a direct pair.

## Unresolved issues

- Some public expression H5 files may be absent from the lightweight repository archive; when absent, match status is reported as fail rather than filled by row-order merging.
- GSE154600 T61/T77 requires optional raw matrix audit for a unique resolution.
