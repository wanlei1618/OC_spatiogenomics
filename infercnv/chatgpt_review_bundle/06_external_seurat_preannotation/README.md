# External scRNA-seq Seurat pre-annotation

This directory contains the reproducible workflow and reviewable outputs from
independent Seurat preprocessing of five external ovarian-cancer-related
single-cell datasets.

Large Seurat objects and raw count matrices are intentionally excluded from
GitHub. The uploaded material contains the scripts, configuration, input
audits, QC summaries, multi-resolution clustering outputs, marker tables,
figures, run logs, and blank manual-annotation templates.

## Dataset summary

| Dataset | Biological role | Input cells | Cells after QC/doublet filtering | Clusters at resolution 0.6 | Marker rows | Status |
|---|---|---:|---:|---:|---:|---|
| GSE147082 | Ovarian tumor ecosystem | 9,885 | 7,645 | 21 | 15,774 | Complete, waiting for manual annotation |
| GSE151214 | Normal fallopian tube reference | 72,132 | 48,837 | 26 | 15,159 | Complete, waiting for manual annotation |
| GSE154600 | Ovarian tumor ecosystem | 52,121 | 31,103 | 30 | 14,961 | Complete, waiting for manual annotation |
| GSE154763 | Myeloid-enriched reference | — | — | — | — | Blocked: public matrix is normalized/non-integer and no raw counts were available |
| GSE158722 | Malignant-fluid tumor ecosystem | 145,515 | 98,832 | 34 | 52,380 | Complete, waiting for manual annotation |

## Quick review links

Every PDF figure also has a 200-DPI PNG with the same basename in the same
directory. Use the PNG links for direct GitHub and ChatGPT visual review; the
PDF files remain the source-quality versions.

### GSE147082

- [Run status](results/GSE147082/logs/run_status.json)
- [QC retention](results/GSE147082/01_qc/qc_cell_retention.csv)
- [Cluster counts](results/GSE147082/02_clustering/cluster_cell_counts.csv)
- [Primary UMAP (PNG)](results/GSE147082/02_clustering/umap_primary_resolution.png) ([PDF](results/GSE147082/02_clustering/umap_primary_resolution.pdf))
- [Top 20 markers](results/GSE147082/03_markers/top20_markers_per_cluster.csv)
- [Marker dot plot (PNG)](results/GSE147082/03_markers/broad_marker_dotplot.png) ([PDF](results/GSE147082/03_markers/broad_marker_dotplot.pdf))
- [Manual annotation template](results/GSE147082/04_manual_annotation/manual_annotation_template.csv)

### GSE151214

- [Run status](results/GSE151214/logs/run_status.json)
- [QC retention](results/GSE151214/01_qc/qc_cell_retention.csv)
- [Cluster counts](results/GSE151214/02_clustering/cluster_cell_counts.csv)
- [Primary UMAP (PNG)](results/GSE151214/02_clustering/umap_primary_resolution.png) ([PDF](results/GSE151214/02_clustering/umap_primary_resolution.pdf))
- [Top 20 markers](results/GSE151214/03_markers/top20_markers_per_cluster.csv)
- [Marker dot plot (PNG)](results/GSE151214/03_markers/broad_marker_dotplot.png) ([PDF](results/GSE151214/03_markers/broad_marker_dotplot.pdf))
- [Manual annotation template](results/GSE151214/04_manual_annotation/manual_annotation_template.csv)

### GSE154600

- [Run status](results/GSE154600/logs/run_status.json)
- [QC retention](results/GSE154600/01_qc/qc_cell_retention.csv)
- [Cluster counts](results/GSE154600/02_clustering/cluster_cell_counts.csv)
- [Primary UMAP (PNG)](results/GSE154600/02_clustering/umap_primary_resolution.png) ([PDF](results/GSE154600/02_clustering/umap_primary_resolution.pdf))
- [Top 20 markers](results/GSE154600/03_markers/top20_markers_per_cluster.csv)
- [Marker dot plot (PNG)](results/GSE154600/03_markers/broad_marker_dotplot.png) ([PDF](results/GSE154600/03_markers/broad_marker_dotplot.pdf))
- [Manual annotation template](results/GSE154600/04_manual_annotation/manual_annotation_template.csv)

### GSE154763

- [Input audit](results/GSE154763/00_input_audit/input_audit.csv)
- [Blocked run status](results/GSE154763/logs/run_status.json)

### GSE158722

- [Run status](results/GSE158722/logs/run_status.json)
- [QC retention](results/GSE158722/01_qc/qc_cell_retention.csv)
- [Cluster counts](results/GSE158722/02_clustering/cluster_cell_counts.csv)
- [Primary UMAP (PNG)](results/GSE158722/02_clustering/umap_primary_resolution.png) ([PDF](results/GSE158722/02_clustering/umap_primary_resolution.pdf))
- [Top 20 markers](results/GSE158722/03_markers/top20_markers_per_cluster.csv)
- [Marker dot plot (PNG)](results/GSE158722/03_markers/broad_marker_dotplot.png) ([PDF](results/GSE158722/03_markers/broad_marker_dotplot.pdf))
- [Manual annotation template](results/GSE158722/04_manual_annotation/manual_annotation_template.csv)

The complete positive-marker tables are stored as
`results/<dataset>/03_markers/all_cluster_markers.csv.gz`. Average-expression
matrices are stored as `cluster_average_expression.csv.gz`.

## Workflow

- [Main Seurat workflow](workflow/scripts/run_seurat_preannotation.R)
- [PowerShell launcher](workflow/scripts/run_skill.ps1)
- [GEO raw-input preparation](workflow/scripts/prepare_geo_raw_inputs.R)
- [GSE158722 metadata correction](workflow/scripts/fix_gse158722_metadata.R)
- [Marker-stage checkpoint recovery](workflow/scripts/resume_markers.R)
- [Clean artifacts and rescue misassigned lineages](workflow/scripts/11_clean_and_rescue_annotation.R)
- [Exact local configuration used](workflow/config/five_external_datasets.local-used.yaml)
- [Portable configuration template](workflow/config/five_external_datasets.example.yaml)
- [Original workflow instructions](workflow/SKILL.md)
- [ChatGPT review guide](CHATGPT_REVIEW_GUIDE.md)

## Analysis outline

1. Read raw integer count matrices and audit their count-like properties.
2. Apply per-sample QC thresholds based on fixed minimums and MAD bounds.
3. Run `scDblFinder` sequentially per biological sample to control memory.
4. Apply `LogNormalize`, select 3,000 variable genes, scale, and run 50 PCs.
5. Use the first 30 PCs for neighbors and UMAP.
6. Cluster independently at resolutions 0.2, 0.4, 0.6, 0.8, 1.0, and 1.2.
7. Use resolution 0.6 as the primary clustering and run positive Wilcoxon
   `FindAllMarkers` with `min.pct=0.20`, `logfc.threshold=0.25`, and adjusted
   p-value threshold 0.05.
8. Leave `cell_type_manual` and `cell_subtype_manual` empty for expert review.

## Step 11: clean and rescue annotation

Run the two datasets sequentially, reviewing GSE154600 before starting
GSE158722:

```text
Rscript workflow/scripts/11_clean_and_rescue_annotation.R --datasets GSE154600 --force
Rscript workflow/scripts/11_clean_and_rescue_annotation.R --datasets GSE158722 --force
```

The local outputs are written below
`diagnostics_v2_marker_ready_cleaned/<dataset>/`. Only the small cluster-level
annotation, rescue, and summary tables are kept in this review bundle; full
cell assignments, removed-cell tables, marker tables, caches, and logs remain
local.

## Steps 12-14: remaining scRNA datasets

Run these scripts in order:

```text
Rscript workflow/scripts/12_reanalyze_GSE147082_after_mt_qc.R
Rscript workflow/scripts/13_clean_normal_reference_GSE151214.R
python workflow/scripts/14_audit_and_score_GSE154763_reference.py
```

Outputs are written locally below `diagnostics_v3_remaining_datasets/`.
GSE147082 starts from the 6,993 repaired-mitochondrial-QC cells. GSE151214 is
restricted to its normal fallopian-tube reference role. GSE154763 uses only
author annotations and normalized-expression module scores after a unique ID
match; it never treats normalized expression as raw counts. Large RDS,
expression, marker, and full per-cell files remain local and are not committed.

## Steps 15-18: final external scRNA validation

Run the targeted corrections and sample/patient-level validation in order:

```text
Rscript workflow/scripts/15_targeted_annotation_corrections.R
Rscript workflow/scripts/16_cross_dataset_macrophage_state_validation.R
Rscript workflow/scripts/17_build_spp1_cd44_itgb1_context.R
Rscript workflow/scripts/18_generate_external_scrna_final_report.R
```

These steps reuse existing results and lineage count inputs; they do not repeat
whole-dataset QC or clustering. Cross-dataset conclusions use within-dataset
percentiles, positive fractions, and sample/patient-level prevalence. Local
outputs are written under `diagnostics_v4_cross_dataset_validation/`; full
per-cell refinements and score tables stay on the D drive.

## Steps 19-22: evidence calibration and freeze

Run the final calibration after steps 15-18:

```text
Rscript workflow/scripts/19_recalibrate_macrophage_state_evidence.R
Rscript workflow/scripts/20_refine_GSE147082_cluster6_and_cnv.R
Rscript workflow/scripts/21_rebuild_external_evidence_matrix_v2.R
Rscript workflow/scripts/22_generate_final_external_report_v2.R
```

These steps separate state presence, cross-patient reproducibility, and
within-dataset relative enrichment; run real targeted markers for GSE147082
cluster 6; and add the PT-2834 patient-internal CNV-like sensitivity analysis.
Outputs are written under `diagnostics_v5_final_calibration/`. The v2 evidence
matrix and report are the frozen external scRNA interpretation.

## Excluded large files

The following local files were excluded because they total approximately
10.3 GB and exceed normal GitHub file limits:

- `GSE147082_preannotation.rds`
- `GSE151214_preannotation.rds`
- `GSE154600_preannotation.rds`
- `GSE158722_preannotation.rds`

Raw/prepared count objects are also excluded. The uploaded tables and figures
are sufficient for QC, clustering, marker, and manual-annotation review.

