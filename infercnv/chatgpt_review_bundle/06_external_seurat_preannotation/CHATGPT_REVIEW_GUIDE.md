# ChatGPT review guide

This document is an index for reviewing the external Seurat pre-annotation
without access to the large RDS objects.

## Recommended review order

For each completed dataset:

1. Read `logs/run_status.json` for final cell, cluster, marker, and doublet
   counts.
2. Read `00_input_audit/input_audit.csv` to confirm the input was count-like.
3. Read `01_qc/qc_thresholds.csv` and `01_qc/qc_cell_retention.csv` to assess
   filtering severity.
4. Review `02_clustering/umap_primary_resolution.pdf`,
   `umap_by_sample.pdf`, and `clustree.pdf` for cluster/sample structure and
   resolution stability.
5. Read `03_markers/top20_markers_per_cluster.csv` and inspect
   `broad_marker_dotplot.pdf` and `top_marker_heatmap.pdf`.
6. Propose cell types and subtypes in
   `04_manual_annotation/manual_annotation_template.csv` while preserving the
   existing cluster IDs.

## Suggested annotation hierarchy

Start with broad compartments before assigning subtypes:

- epithelial/tumor: EPCAM, KRT8, KRT18, KRT19, MSLN, WFDC2
- T/NK: PTPRC, CD3D, CD3E, TRBC1, NKG7, GNLY
- B/plasma: CD79A, MS4A1, MZB1, JCHAIN
- myeloid: LYZ, LST1, TYROBP, FCER1G, C1QA/B/C, SPP1
- fibroblast/stromal: COL1A1, COL1A2, DCN, COL3A1
- endothelial: PECAM1, VWF, KDR
- cycling: MKI67, TOP2A

Use multiple concordant markers and inspect sample distribution before
assigning rare or disease-specific subtypes. Do not infer malignancy from one
marker alone.

## Dataset-specific interpretation

- **GSE151214** is the normal fallopian tube reference.
- **GSE154763** must not be treated as a raw-count analysis result. It is
  intentionally blocked because the available matrix is normalized and
  non-integer.
- **GSE154600** preserves available author metadata fields by barcode for
  comparison, but the uploaded `cell_type_manual` fields remain blank.
- **GSE158722** contains 23 patients and 39 biological sample/timepoint groups;
  P01-P09 have longitudinal timepoints, P10-P17 are pre-treatment, and P18-P24
  are post-treatment.

## Review questions

1. Are any clusters dominated by one sample or patient?
2. Do adjacent resolutions split biologically coherent lineages or technical
   states?
3. Are low-cell clusters supported by multiple markers?
4. Do epithelial clusters show plausible lineage/tumor programs rather than
   ambient RNA alone?
5. Are myeloid subtypes supported by coherent C1QC/SPP1/inflammatory programs?
6. Are cycling clusters better represented as states layered onto a parent
   lineage?
7. Which proposed labels should be marked low confidence and revisited with
   additional references?

## Files not directly previewable on GitHub

The full marker and average-expression files are gzip-compressed CSVs. They can
be downloaded and read with R, Python, or standard gzip-aware tools. The
uncompressed Top 20/50/100 tables are included for browser and ChatGPT review.

