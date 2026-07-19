# External validation datasets

This directory contains the public external validation datasets used by the
inferCNV sample-type ligand-receptor niche analysis.

## Single-cell RNA-seq

- `Zhang2022_Ovarian_scRNA/`
  - Source used in the analysis: Zhang2022 ovarian single-cell RNA-seq.
  - Files used: `Cells.csv`, `Samples.csv`, `Genes.txt`,
    `Exp_data_UMIcounts.mtx`.
  - Analysis role: external scRNA validation of SPP1 myeloid source cells and
    CD44/ITGB1 target tumor programs.

- `OV_GSE147082_scRNA/`
- `OV_GSE151214_scRNA/`
- `OV_GSE154600_scRNA/`
- `OV_GSE154763_scRNA/`
- `OV_GSE158722_scRNA/`
  - Source used in the extended validation: ovarian cancer 10x-style h5
    expression matrices.
  - Files used: one `OV_GSE*_expression.h5` file per dataset. The
    `OV_GSE158722_expression.h5` source file is stored as ordered split parts
    because the complete file is larger than GitHub's single LFS object limit.
  - Analysis role: extended external scRNA validation of marker-inferred
    SPP1+ myeloid/macrophage source cells and ITGB1/CD44 target-high tumor
    programs. These h5 files do not include curated cell metadata in the
    matrix group, so the extended analysis uses marker-score cell-type
    inference.

## Spatial transcriptomics

- `spatial/GSE203612/seurat_visium/GSM6177614/`
- `spatial/GSE203612/seurat_visium/GSM6177617/`
- `spatial/GSE203612/seurat_visium/GSM6177618/`
  - Source used in the analysis: GSE203612 ovarian Visium samples.
  - Files used: `filtered_feature_bc_matrix.h5` plus the `spatial/` coordinate
    and image files.
  - Analysis role: spot-level scoring, correlation analysis, and coordinate-
    based neighborhood enrichment.

- `spatial/GSE189843/suppl/`
  - Source used in the analysis: GSE189843 HGSC Visium supplement files.
  - Files used: GSM5708485 through GSM5708496 matrix/barcode/feature files.
  - Analysis role: spot-level scoring and ligand-receptor expression-product
    summaries. The local copy does not include usable spot coordinate files for
    neighborhood enrichment.

- `spatial/GSE211956/seurat_visium/`
  - Source used in the extended validation: GSE211956 Visium samples
    GSM6506110 through GSM6506117.
  - Files used: 10x matrix triplets under `filtered_feature_bc_matrix/` plus
    the `spatial/` coordinate and image files.
  - Analysis role: spot-level signature scoring, SPP1 macrophage to target-axis
    correlation, coordinate neighborhood enrichment, focused spatial LR scores,
    and spatial virtual KO.

- `spatial/GSE227019/seurat_visium/`
  - Source used in the extended validation: GSE227019 Visium samples
    GSM7090083 through GSM7090088.
  - Files used: `filtered_feature_bc_matrix.h5` plus the `spatial/` coordinate
    and image files.
  - Analysis role: spot-level signature scoring, SPP1 macrophage to target-axis
    correlation, coordinate neighborhood enrichment, focused spatial LR scores,
    and spatial virtual KO.

## Analysis outputs

- Initial external validation outputs are under
  `infercnv/05_sample_type_lr_niche_external_validation/`.
- Extended validation outputs for the added five scRNA datasets and two spatial
  datasets are under `infercnv/07_extended_external_validation/`.

Large matrix, archive, HDF5, compressed, and image files are tracked with Git
LFS.
