# single-cell-seurat-clustering-markers

将本文件夹复制到：

```text
D:\OC_spatiogenomics\.agents\skills\single-cell-seurat-clustering-markers
```

先运行输入审计：

```powershell
cd D:\OC_spatiogenomics\.agents\skills\single-cell-seurat-clustering-markers

powershell -ExecutionPolicy Bypass `
  -File scripts\run_skill.ps1 `
  -Config config\five_external_datasets.yaml `
  -AuditOnly
```

确认路径和矩阵正确后再完整运行：

```powershell
powershell -ExecutionPolicy Bypass `
  -File scripts\run_skill.ps1 `
  -Config config\five_external_datasets.yaml
```

主要R包：

```r
install.packages(c("Seurat","yaml","dplyr","ggplot2","patchwork","jsonlite","clustree"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment","SummarizedExperiment","scDblFinder"))
```

完成后填写每个数据集的：

```text
04_manual_annotation/manual_annotation_template.csv
```

Skill不会自动写入最终cell type。
