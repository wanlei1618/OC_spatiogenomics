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

Large matrix, archive, HDF5, compressed, and image files are tracked with Git
LFS.
