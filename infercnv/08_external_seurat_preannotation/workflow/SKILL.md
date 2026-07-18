---
name: single-cell-seurat-clustering-markers
description: >
  Analyze scRNA-seq datasets independently with Seurat from count-matrix
  loading through QC, optional doublet removal, normalization, PCA, UMAP,
  multi-resolution clustering and FindAllMarkers. Export cluster marker
  evidence for manual annotation and stop before final cell-type assignment.
---

# Single-cell Seurat clustering and marker discovery

## Purpose

Run a reproducible pre-annotation workflow for every dataset independently:

```text
counts → QC → doublet check → normalization → PCA → clustering
→ UMAP → FindAllMarkers → manual annotation template → STOP
```

The user assigns final cell types.

## Mandatory boundaries

The agent MUST:

- analyze each dataset separately;
- preserve raw counts in the RNA assay;
- keep `seurat_clusters`;
- create empty `cell_type_manual` and `cell_subtype_manual` fields;
- export all positive markers and top 20/50/100 markers per cluster;
- save a pre-annotation Seurat object;
- stop after marker discovery.

The agent MUST NOT:

- assign final cell types automatically;
- merge the five datasets before independent review;
- overwrite user annotations;
- use a single marker to name a cluster;
- run inferCNV, CopyKAT, CellChat, LIANA, NicheNet or survival analysis;
- treat GSE151214 as ovarian tumor;
- treat GSE154763 as a complete tumor ecosystem;
- join metadata by row order;
- use placeholder labels as completed annotations.

Automated reference labels, when explicitly requested, must be stored only as
`provisional_reference_label`.

## Default datasets

- GSE147082: ovarian tumor ecosystem.
- GSE151214: normal fallopian-tube reference.
- GSE154600: ovarian tumor ecosystem.
- GSE154763: myeloid-enriched reference.
- GSE158722: malignant-fluid ovarian tumor ecosystem.

## Supported inputs

- `10x_h5`
- `10x_dir`
- `rds_seurat`
- `rds_sce`
- `matrix_rds`

## Required workflow

### 1. Input audit

Report:

- path and file size;
- matrix dimensions;
- integer-like count fraction;
- duplicated genes and cell IDs;
- metadata columns;
- sample fields when present.

Stop on empty matrices, duplicated cell IDs, inaccessible count layers, or
non-integer matrices declared as raw counts.

### 2. QC

Calculate:

- `nFeature_RNA`
- `nCount_RNA`
- `percent.mt`
- `percent.ribo`
- `percent.HB`

Use sample-wise median/MAD thresholds when a valid sample column exists.
Otherwise use dataset-wise thresholds and record the limitation.

Defaults:

- minimum genes: 200;
- minimum UMIs: 500;
- maximum mitochondrial fraction: 25%;
- MAD multiplier: 3.

### 3. Doublets

Use `scDblFinder` when installed and enabled. Preserve its score and class.
If unavailable, continue only when `allow_doublet_skip: true`.

### 4. Seurat workflow

Default:

```text
NormalizeData(LogNormalize)
FindVariableFeatures(vst, 3000)
ScaleData
RunPCA(50 PCs)
FindNeighbors
FindClusters at 0.2–1.2
RunUMAP
```

Do not integrate datasets during this stage.

### 5. Markers

Use the RNA assay and run:

```r
FindAllMarkers(
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.20,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)
```

Export complete and top-ranked marker tables. Do not remove mitochondrial,
ribosomal, immunoglobulin or cell-cycle genes from the complete table.

### 6. Manual annotation handoff

Create:

```text
manual_annotation_template.csv
```

with:

```text
dataset_id,seurat_cluster,n_cells,top_markers,
cell_type_manual,cell_subtype_manual,confidence,notes
```

The last four fields must remain empty.

Final status:

```text
PREANNOTATION_COMPLETE_WAITING_FOR_MANUAL_CELLTYPE
```

## Dataset-specific constraints

### GSE147082

Run independent clustering and marker discovery. No final annotation.

### GSE151214

Treat as normal fallopian tube. Secretory, ciliated and other lineages are
hypotheses for manual review, not automatic labels.

### GSE154600

Preserve author metadata including `celltype`, `subtype`, `Cluster`,
`hpca.celltype` and `encode.celltype` when present. Never use `subtype` as a
cell-type field.

### GSE154763

Treat as myeloid enriched. Do not require epithelial/tumor populations.

### GSE158722

Preserve patient, sample and timepoint metadata. Cell IDs must be unique across
samples. Prefer:

```text
patient_id__sample_id__barcode
```

Never join samples with naked 10x barcodes.

## Output structure

```text
<output_root>/<dataset_id>/
├── 00_input_audit/
├── 01_qc/
├── 02_clustering/
├── 03_markers/
├── 04_manual_annotation/
├── logs/
└── objects/
```

Required outputs:

```text
00_input_audit/input_audit.csv
01_qc/qc_thresholds.csv
01_qc/qc_cell_retention.csv
01_qc/qc_metadata.csv.gz
02_clustering/cluster_cell_counts.csv
02_clustering/cluster_by_sample_counts.csv
02_clustering/umap_primary_resolution.pdf
02_clustering/umap_by_sample.pdf
03_markers/all_cluster_markers.csv.gz
03_markers/top20_markers_per_cluster.csv
03_markers/top50_markers_per_cluster.csv
03_markers/top100_markers_per_cluster.csv
03_markers/cluster_average_expression.csv.gz
03_markers/broad_marker_dotplot.pdf
03_markers/top_marker_heatmap.pdf
04_manual_annotation/manual_annotation_template.csv
objects/<dataset_id>_preannotation.rds
logs/sessionInfo.txt
logs/run_status.json
```

## Acceptance gates

Pass only if:

- at least 200 cells remain;
- the primary clustering contains at least 2 clusters;
- every cluster has exported markers;
- marker columns include gene, cluster, avg_log2FC, pct.1, pct.2 and p_val_adj;
- `cell_type_manual` remains empty;
- no downstream mechanism analysis has run.

## Execution

Audit first:

```powershell
powershell -ExecutionPolicy Bypass `
  -File scripts\run_skill.ps1 `
  -Config config\five_external_datasets.yaml `
  -AuditOnly
```

Then run the full workflow after reviewing the audit:

```powershell
powershell -ExecutionPolicy Bypass `
  -File scripts\run_skill.ps1 `
  -Config config\five_external_datasets.yaml
```
