# sample_type-dependent SPP1-CD44/ITGB1 LR niche analysis

Project root: `D:/OC_spatiogenomics/infercnv/sample_type_LR_niche_analysis`

Spatial-data curation update: 2026-07-12

## 1. Scope

This report integrates:

1. exploratory sample-type analysis in `integrated_oc`;
2. external ovarian single-cell validation;
3. curated spatial-transcriptomics validation;
4. explicit separation of coordinate-aware and expression-only evidence.

The central hypothesis remains that SPP1-positive myeloid/macrophage cells may participate in a niche-dependent program involving CD44-positive or ITGB1-positive KRAS/hypoxia-high ovarian cancer CNV subclones. Ligand-receptor and spatial scores are interpreted as candidate interaction opportunities, not proof of physical binding or causal signaling.

## 2. Integrated_oc exploratory analysis

`integrated_oc` remains exploratory because `sample_type` is one-to-one with `sample_id`. Cell-level P values cannot substitute for patient-level replication.

The primary-axis opportunity scores are strongest in solid-tumor-like contexts and are focused on `CNV_Subclone_02/04`. Full source-group and sample-type results are retained in:

- `tables/sample_type_LR_opportunity_scores_primary_axes_target_Subclone02_04.csv`
- `tables/sample_type_score_summary.csv`

These findings generate a niche hypothesis but do not establish a general sample-type effect.

## 3. External scRNA validation: Zhang2022 ovarian

Macrophage/dendritic cells were treated as candidate SPP1-source cells and malignant cells as target cells. Peritoneum, omentum and mesentery were grouped as implant-like niches.

Aggregated axis scores were higher in peritoneal-implant-like samples than in the primary-tumor category:

| sample_type | SPP1-CD44 | SPP1-ITGB1 | MIF-CD74 |
|---|---:|---:|---:|
| peritoneal_implant | 0.005744 | 0.020217 | 0.015222 |
| tumor | 0.00000115 | 0.00000186 | 0.020403 |

This result is expression-level support only. Anatomical site, patient and treatment phase are not fully balanced.

Key outputs:

- `tables/external_scRNA_all_datasets_sampletype_axis_scores.csv`
- `tables/external_scRNA_patient_sampletype_summary.csv`
- `tables/external_scRNA_target_signature_scores.csv`
- `tables/external_scRNA_SPP1_myeloid_CD44_ITGB1_target_scores.csv`

## 4. Spatial-transcriptomics data curation

### 4.1 GSE203612

The GSE203612 series is pan-cancer. Only the following Visium samples are ovarian carcinoma:

| sample_id | GEO title | inclusion |
|---|---|---|
| GSM6177614 | NYU_OVCA1_Vis | included |
| GSM6177617 | NYU_OVCA3_Vis | included |
| GSM6177618 | NYU_PDAC1_Vis | excluded |

`GSM6177618` is titled and sourced as primary pancreatic ductal adenocarcinoma. GEO contains internally inconsistent characteristic fields that say ovarian carcinoma/ovary, but the specific sample title and source identify PDAC. It is excluded from all ovarian summaries and neighborhood analyses.

Both valid ovarian samples release filtered matrices, tissue positions, scalefactors and tissue images, and therefore support coordinate-aware exploratory analysis.

### 4.2 GSE189843

GSE189843 contains 12 pretreatment high-grade serous ovarian carcinoma samples:

- `GSM5708485`–`GSM5708490`: excellent response to neoadjuvant chemotherapy;
- `GSM5708491`–`GSM5708496`: poor response.

The GEO supplementary archive provides count matrices, barcodes, features and tissue images, but no spot-position or scalefactor files. Consequently, these samples are used only for expression-level spot scoring and sample-level summaries. They are not used for neighborhood, distance or coordinate-dependent ligand-receptor claims.

The curated sample metadata are stored in:

- `spatial_transcriptomics_upgrade/metadata/spatial_sample_manifest.csv`
- `spatial_transcriptomics_upgrade/metadata/download_manifest.csv`

## 5. Corrected spatial results

### 5.1 Coordinate-aware GSE203612 ovarian samples

| sample_id | n_spots | Spearman r | direction | neighborhood enrichment ratio |
|---|---:|---:|---|---:|
| GSM6177614 | 1762 | -0.197 | negative | 0.893 |
| GSM6177617 | 1661 | 0.221 | positive | 1.100 |

The two valid ovarian samples show discordant directions. One sample has lower-than-global target-neighbor frequency around SPP1-myeloid-high spots, whereas the other shows only mild enrichment.

Therefore GSE203612 does not support a universal spatial colocalization claim. It supports only heterogeneous, sample-dependent spatial association.

### 5.2 Expression-only GSE189843 samples

Across the six excellent responders, the median sample-level Spearman correlation was 0.131; five of six samples were positive and three remained significant after within-dataset BH correction.

Across the six poor responders, the median sample-level Spearman correlation was 0.160; five of six samples were positive and four remained significant after within-dataset BH correction.

These summaries should not be interpreted as a validated response-group difference because:

- there are only six patients per group;
- spots within a sample are not independent patient replicates;
- tissue coordinates are unavailable;
- no formal patient-level differential model has yet been applied.

Corrected tables:

- `tables/spatial_correlation_SPP1_myeloid_Target_axis.csv`
- `tables/spatial_neighborhood_enrichment.csv`
- `tables/limitations_summary.csv`

## 6. Improved reproducible workflow

A reproducible spatial data layer has been added under:

`spatial_transcriptomics_upgrade/`

The workflow:

1. downloads only curated GEO inputs to the D drive;
2. creates SpaceRanger-compatible GSE203612 directories;
3. builds coordinate-aware GSE203612 and expression-only GSE189843 Seurat objects;
4. performs QC and gene-program scoring;
5. limits neighborhood tests to samples with released coordinates;
6. transfers myeloid and CNV-subclone reference states from `integrated_oc`;
7. audits outputs to prevent excluded samples or coordinate-unavailable samples from entering invalid analyses.

See `spatial_transcriptomics_upgrade/README.md` for commands.

## 7. Interpretation

The curated public spatial data support the following cautious model:

> SPP1-myeloid and CD44/ITGB1/KRAS-hypoxia-related tumor programs can coexist in ovarian cancer tissue, but the strength and direction of spatial association vary between samples. Current data support a candidate niche-dependent program rather than a universal spatial interaction.

The evidence hierarchy is:

1. **SPP1-CD44 candidate communication:** supported by single-cell expression and ligand-receptor database evidence;
2. **SPP1-associated ITGB1-positive adhesion program:** supported by expression co-occurrence, but not proof of direct SPP1-ITGB1 binding;
3. **spatial colocalization:** heterogeneous and exploratory;
4. **patient-level clinical association:** not established.

## 8. Limitations and next steps

| limitation | required action |
|---|---|
| only two valid coordinate-aware GSE203612 ovarian samples | add independent ovarian Visium/CosMx/Xenium cohorts |
| GSE189843 coordinates absent from GEO supplement | request original SpaceRanger spatial files from the authors |
| spot-level P values affected by spatial autocorrelation | use spatially constrained null models and patient-level meta-analysis |
| Visium spots contain mixed cells | run cell2location/RCTD with `integrated_oc` as reference |
| expression-product scores are not physical interaction evidence | run formal COMMOT/NicheNet and validate by RNAscope or multiplex IF |
| CNV labels are inferred from scRNA reference | validate mapped niches with DNA/CNV-aware pathology or targeted assays |
