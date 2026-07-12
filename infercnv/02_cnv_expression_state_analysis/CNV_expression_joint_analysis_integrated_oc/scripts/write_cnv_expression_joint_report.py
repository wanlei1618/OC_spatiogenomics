# -*- coding: utf-8 -*-
from pathlib import Path
import pandas as pd

base = Path(r"D:\OC_spatiogenomics\infercnv\CNV_expression_joint_analysis_integrated_oc")
tables = base / "tables"
figures = base / "figures"

burden = pd.read_csv(tables / "CNV_clone_burden_summary.csv")
func = pd.read_csv(tables / "CNV_clone_functional_and_TF_activity_summary.csv")
dosage = pd.read_csv(tables / "CNV_expression_dosage_genes_positive_rho_padj005.csv")
cor_pairs = pd.read_csv(tables / "CNV_transcriptome_coupling_spearman_pairs.csv")
markers = pd.read_csv(tables / "CNV_clone_top20_markers_per_clone.csv")
pb_meta = pd.read_csv(tables / "pseudo_bulk_sample_clone_metadata.csv")

def md_table(df, cols, n=20, float_cols=None):
    float_cols = float_cols or []
    lines = []
    for _, row in df.head(n).iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if c in float_cols:
                try:
                    v = f"{float(v):.4g}"
                except Exception:
                    pass
            vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

def top_feature(feature):
    x = func[func["feature"] == feature].sort_values("mean_score", ascending=False)
    return "; ".join([f"{r.CNV_clone}={r.mean_score:.3f}" for r in x.head(5).itertuples()])

features = ["EMT", "Hypoxia", "KRAS_UP", "KRAS_DN", "LAPTM5_axis", "Proliferation",
            "Immune_modulatory", "TF_SNAI2", "TF_ATF3", "TF_HIF1A", "TF_STAT3", "TF_MYC", "TF_FOXM1"]

marker_summary = []
for cl, x in markers.groupby("cluster"):
    genes = ", ".join(x.sort_values(["p_val_adj", "avg_log2FC"], ascending=[True, False]).head(10)["gene"])
    marker_summary.append({"CNV_clone": cl, "top_marker_genes": genes})
marker_summary = pd.DataFrame(marker_summary)

report = f"""# CNV + 表达联合分析：integrated_oc

结果目录：`{base}`

## 1. 按文档计划拆解的任务与完成情况

| 模块 | 文档计划 | 本次实现 | 主要输出 |
|---|---|---|---|
| 肿瘤细胞准备 | 使用已回填的 CNV clone，提取肿瘤细胞 | 已完成；6507 个 CNV observation 细胞，5 个 CNV subclone | `tables/integrated_oc_metadata_with_CNV_clone_and_subtypes.csv` |
| 单细胞 clone DEG | FindAllMarkers / MAST | 用 RNA@data 做 one-vs-rest Wilcoxon marker 分析 | `tables/CNV_clone_DEG_single_cell_wilcoxon_one_vs_rest.csv` |
| pseudo-bulk DEG | sample × CNV_clone 聚合后 edgeR/limma | 用 sample_id × CNV_clone 聚合 counts，并用 limma-trend 比较各 clone vs Subclone_01 | `tables/pseudo_bulk_DEG_limma_clone_vs_Subclone_01.csv` |
| CNV burden | inferCNV gene matrix 计算每细胞 CNV 强度 | 用 pred_cnv_genes.dat 构建 HMM state-3 gene × subcluster 矩阵，映射到细胞 | `tables/CNV_clone_burden_summary.csv`, `figures/CNV_burden_by_clone_violin.png` |
| chromosome-level CNV score | 按染色体聚合 CNV signal | 已完成，每个 clone 的 chromosome-level CNV score | `tables/chromosome_level_CNV_score_by_clone.csv`, `figures/chromosome_level_CNV_score_by_clone_heatmap.png` |
| CNV-expression dosage | CNV signal 与表达相关 | 在 52 个 inferCNV subcluster 水平计算 gene-wise Spearman correlation | `tables/CNV_expression_dosage_correlation_by_infercnv_subcluster.csv` |
| EMT/Hypoxia/KRAS/LAPTM5 | signature score 比较 | 用 curated gene sets 计算表达 signature score | `tables/tumor_cell_CNV_burden_signature_TF_scores.csv` |
| TF activity | DoRothEA/decoupleR | 当前环境无 dorothea/decoupleR，改用 TF target signature proxy score | `tables/CNV_clone_functional_and_TF_activity_summary.csv` |
| 耦合模型 | CNV burden 与 EMT/Hypoxia/KRAS/TF 相关 | 已完成 Spearman 相关矩阵和配对显著性 | `tables/CNV_transcriptome_coupling_spearman_*` |

## 2. CNV burden 排名

| CNV_clone | mean | median | sd |
|---|---:|---:|---:|
{md_table(burden.sort_values("mean", ascending=False), ["CNV_clone", "mean", "median", "sd"], 10, ["mean", "median", "sd"])}

结论：`Subclone_04` 的 CNV burden 最高，其次为 `Subclone_02/05/03`；`Subclone_01` 最低。

## 3. 功能 signature 与 TF target proxy activity

| feature | clone ranking by mean score |
|---|---|
""" + "\n".join([f"| {f} | {top_feature(f)} |" for f in features]) + f"""

解释要点：
- `Subclone_04`：CNV burden 最高，同时 Hypoxia、Immune-modulatory、TF_HIF1A、TF_STAT3 较高。
- `Subclone_02`：KRAS_UP、TF_ATF3、TF_STAT3、TF_MYC 较高，呈 KRAS/ATF3/STAT3-active 表型。
- `Subclone_01`：CNV burden 最低，但 EMT、LAPTM5_axis、TF_SNAI2、TF_FOXM1 相对更高，提示它更像 EMT/LAPTM5-high transcriptional-state clone，而不是 high-CNV clone。
- `Subclone_05`：细胞数少，仅 96 个，所有解释应谨慎。

## 4. clone marker genes

| CNV_clone | top marker genes |
|---|---|
{md_table(marker_summary, ["CNV_clone", "top_marker_genes"], 10)}

完整 marker 表见：

- `tables/CNV_clone_DEG_single_cell_wilcoxon_one_vs_rest.csv`
- `tables/CNV_clone_top20_markers_per_clone.csv`
- `figures/CNV_clone_top_marker_average_expression_heatmap.png`

## 5. CNV-expression dosage genes

基于 inferCNV subcluster 层面的 CNV HMM state 与平均表达做 Spearman correlation。正相关且 BH 校正显著的前 30 个候选：

| gene | rho | padj | n_groups |
|---|---:|---:|---:|
{md_table(dosage.sort_values("rho", ascending=False), ["gene", "rho", "padj", "n_groups"], 30, ["rho", "padj"])}

解释要点：
- 高相关 dosage candidates 包括 `MSLN`, `WFDC2`, `TACSTD2`, `ELF3`, `KRT8`, `SLPI` 等卵巢癌/上皮相关基因。
- 这些基因更适合作为 “CNV-dosage coupled expression genes” 候选，而不是直接等同于 driver gene。

## 6. CNV-transcriptome coupling

最强相关 pair：

| var1 | var2 | rho | padj |
|---|---|---:|---:|
{md_table(cor_pairs, ["var1", "var2", "rho", "padj"], 30, ["rho", "padj"])}

和 CNV_burden 直接相关的 pair 可在 `tables/CNV_transcriptome_coupling_spearman_pairs.csv` 中筛选 `var1 == CNV_burden` 或 `var2 == CNV_burden`。

## 7. pseudo-bulk 设计

sample × CNV_clone pseudo-bulk 样本数：{pb_meta.shape[0]}

| pb_group | sample_id | CNV_clone |
|---|---|---|
{md_table(pb_meta, ["pb_group", "sample_id", "CNV_clone"], 30)}

注意：当前只有 4 个样本来源，pseudo-bulk DEG 主要用于探索，正式统计解释仍需更多 patient/sample 复制。

## 8. 主要图件

- `figures/CNV_burden_by_clone_violin.png`
- `figures/chromosome_level_CNV_score_by_clone_heatmap.png`
- `figures/CNV_clone_signature_and_TF_activity_heatmap.png`
- `figures/CNV_clone_top_marker_average_expression_heatmap.png`
- `figures/top_CNV_expression_dosage_gene_expression_by_subcluster.png`
- `figures/CNV_transcriptome_coupling_correlation_heatmap.png`

## 9. 推荐结果表述

基于 inferCNV HMM gene-level CNV state，我们将 integrated_oc 中 6507 个肿瘤 observation 细胞划分为 5 个 CNV-defined subclone，并计算每个细胞/clone 的 CNV burden。`Subclone_04` 和 `Subclone_02` 具有较高 CNV burden，并分别呈现 hypoxia/immune-modulatory 与 KRAS/ATF3/STAT3-active 表达状态。`Subclone_01` CNV burden 较低但 EMT、LAPTM5 和 SNAI2 target score 较高，提示 CNV clone 与转录状态并非一一对应，而是存在“CNV burden high clone”和“EMT/LAPTM5 transcriptional-state clone”的分化。CNV-expression dosage 分析进一步识别出 MSLN、WFDC2、TACSTD2、ELF3、KRT8 等与 CNV state 显著正相关的候选 dosage genes，为构建 CNV-driven transcriptional evolution model 提供依据。

## 10. 方法限制

- `TF activity` 使用 TF target signature proxy score，不是 DoRothEA/VIPER 的正式 regulon activity。
- `CNV_burden` 来自 inferCNV HMM predicted gene state，而不是直接读取 inferCNV `expr.data` RDS；这是因为当前 R 环境未安装 inferCNV 包。
- pseudo-bulk 复制数有限，结果作为探索性 DEG 使用。
- `Subclone_05` 细胞数少，统计稳定性弱。
"""

(base / "CNV_expression_joint_analysis_task_breakdown_and_results.md").write_text(report, encoding="utf-8")
(base / "CNV_expression_joint_analysis_task_breakdown_and_results.txt").write_text(report, encoding="utf-8")
print(base / "CNV_expression_joint_analysis_task_breakdown_and_results.md")
