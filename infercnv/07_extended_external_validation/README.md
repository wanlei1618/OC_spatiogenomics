# Extended external validation for SPP1 macrophage to ITGB1/CD44 tumor axis

This directory extends the previous validation by adding five ovarian scRNA h5 datasets and two spatial transcriptomics GEO datasets that were not included in the first pass.

## Added external scRNA datasets

- OV_GSE147082: `D:\OC_spatiogenomics\公开集\单细胞\OV_GSE147082_expression.h5`
- OV_GSE151214: `D:\OC_spatiogenomics\公开集\单细胞\OV_GSE151214_expression.h5`
- OV_GSE154600: `D:\OC_spatiogenomics\公开集\单细胞\OV_GSE154600_expression.h5`
- OV_GSE154763: `D:\OC_spatiogenomics\公开集\单细胞\OV_GSE154763_expression.h5`
- OV_GSE158722: `D:\OC_spatiogenomics\公开集\单细胞\OV_GSE158722_expression.h5`

The h5 files did not contain usable cell metadata in the 10x matrix groups, so broad cell classes were inferred from marker scores. This is a validation-by-expression-potential analysis, not a curated annotation analysis.

## Added spatial datasets

- GSE211956: `D:\OC_spatiogenomics\公开集\ovarian_spatial_geo\GSE211956\seurat_visium`
- GSE227019: `D:\OC_spatiogenomics\公开集\ovarian_spatial_geo\GSE227019\seurat_visium`

## Key outputs

- `tables/extended_external_scRNA_dataset_level_summary.csv`
- `tables/extended_external_scRNA_meta_summary.csv`
- `tables/extended_external_scRNA_celltype_composition.csv`
- `tables/extended_spatial_sample_summary.csv`
- `tables/extended_spatial_neighborhood_enrichment.csv`
- `tables/extended_spatial_LR_scores.csv`
- `tables/extended_spatial_virtual_KO.csv`
- `tables/extended_validation_evidence_matrix.csv`

## Methods summary

Selected genes from the mechanism-axis gene sets were streamed from 10x h5/mtx matrices, normalized as log1p(CP10K), and summarized into module scores. scRNA source cells were marker-inferred myeloid/macrophage cells; target cells were marker-inferred tumor/epithelial cells in the top quartile of the target-axis score. Spatial spots were scored with the same signatures. Neighborhood enrichment used array coordinates and a radius of 2.01 times the median nearest-neighbor distance. Spatial virtual KO set SPP1 or receptor expression to zero in the focused source/target compartments and recalculated focused LR scores.

## scRNA meta summary

| axis | n_dataset | mean_lr_score | median_lr_score | mean_source_fraction | mean_target_fraction |
| --- | --- | --- | --- | --- | --- |
| SPP1-CD44 | 5 | 0.5356942613008834 | 0.5520318707383952 | 0.2687069985326319 | 0.12995464967555426 |
| SPP1-ITGB1 | 5 | 1.4606735365141588 | 0.9107579524179952 | 0.2687069985326319 | 0.12995464967555426 |

## spatial LR meta summary

| axis | n_spatial_samples | mean_focused_lr_score | median_focused_lr_score |
| --- | --- | --- | --- |
| SPP1-CD44 | 14 | 0.9231581189147146 | 0.2721289566804439 |
| SPP1-ITGB1 | 14 | 1.573933305493485 | 0.9848603970404799 |
