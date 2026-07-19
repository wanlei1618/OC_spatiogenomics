# -*- coding: utf-8 -*-
from pathlib import Path
import pandas as pd

plan_dir = Path(r"D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis")
tables = plan_dir / "tables"
figures = plan_dir / "figures"
scripts = plan_dir / "scripts"
source = plan_dir / "source_outputs"

clone_counts = pd.read_csv(tables / "cnv_subclone_cell_counts.csv")
mapping = pd.read_csv(source / "integrated_oc_mapping_summary.csv")
group_counts = pd.read_csv(source / "integrated_oc_interaction_group_cell_counts.csv")
func = pd.read_csv(tables / "cnv_subclone_function_module_score_summary.csv")
tme_top = pd.read_csv(tables / "knn_top_TME_enriched_neighbor_labels.csv")
lr_focus = pd.read_csv(tables / "focused_LR_axes_involving_cnv_subclones.csv")
distinct_lr = pd.read_csv(source / "lr_top_distinct_LR_pairs_involving_cnv_subclones.csv")

def top_module(module, n=3):
    x = func[func["module"] == module].sort_values("mean_score", ascending=False).head(n)
    return "; ".join([f"{r.cnv_subclone}={r.mean_score:.3f}" for r in x.itertuples()])

def table_lines(df, cols, n=10):
    out = []
    for r in df.head(n).itertuples(index=False):
        vals = [getattr(r, c) for c in cols]
        out.append("| " + " | ".join(str(v) for v in vals) + " |")
    return "\n".join(out)

lr_focus_small = lr_focus.sort_values("score", ascending=False).head(20).copy()
lr_focus_small["score"] = lr_focus_small["score"].map(lambda x: f"{x:.3f}")
tme_small = tme_top.head(20).copy()
for c in ["observed_prop", "expected_prop_perm_mean", "enrichment_ratio", "z_score", "empirical_p_greater"]:
    tme_small[c] = tme_small[c].map(lambda x: f"{x:.4g}")

report = f"""# integrated_oc 基于 inferCNV 的 CNV 亚克隆与 TME 互作分析结果

生成目录：`{plan_dir}`

## 1. 按文档计划拆解的任务与完成状态

| 模块 | 文档要求 | 本次完成情况 | 主要输出 |
|---|---|---|---|
| 数据准备 | 使用 integrated_oc RNA assay；整理 cell_type/sample metadata | 已完成，保留 CNV clone、T/NK、myeloid、B 亚型和 interaction_group | `tables/integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv` |
| CNV clone 回填 | 将 inferCNV subclone 定位回 integrated_oc | 已完成，6507 个 CNV observation 细胞映射为 5 个 subclone | `tables/cnv_subclone_cell_counts.csv` |
| B 细胞亚型定位 | 使用 integratedocBcells.RData 的 cell_type | 已完成，2274 个 B 子集细胞映射；未进入 B 子集的 659 个保留为 broad B cells | `source_outputs/integrated_oc_Bcell_subtypes_mapped_from_integratedocBcells.csv` |
| CNV clone 组成 | 比较各 clone 的组成和 sample_type 分布 | 已完成 | `tables/cnv_subclone_composition_by_sample_type.csv`, `figures/cnv_subclone_composition_by_sample_type.png` |
| 功能注释 | EMT、hypoxia、KRAS、LAPTM5、proliferation、immune-modulatory 等评分 | 已完成，基于 RNA@data 计算模块均值 | `tables/cnv_subclone_function_module_score_summary.csv`, `figures/cnv_subclone_function_module_score_heatmap.png` |
| kNN 邻近分析 | PCA/kNN 邻居组成、permutation test | 已完成，PCA 1:30、k=30、1000 次置换 | `tables/knn_neighbor_enrichment_*`, `figures/knn_*` |
| 配体-受体互作 | CellChat/Connectome/NicheNet 比较 clone 与 TME 通讯 | 已完成 Connectome-like OmniPath LR scoring | `source_outputs/lr_*`, `tables/focused_LR_axes_involving_cnv_subclones.csv` |
| 机制模型 | 总结 CNV clone - TME 关系 | 已在本报告下方形成候选模型 | 本报告 |

## 2. CNV subclone 规模

| cnv_subclone | n_cells | fraction_of_cnv_cells |
|---|---:|---:|
{table_lines(clone_counts.assign(fraction_of_cnv_cells=clone_counts['fraction_of_cnv_cells'].map(lambda x: f'{x:.3f}')), ['cnv_subclone','n_cells','fraction_of_cnv_cells'], 10)}

## 3. integrated_oc 映射概况

| item | value |
|---|---|
{table_lines(mapping, ['item','value'], 20)}

## 4. B 细胞亚型映射

`integratedocBcells.RData` 中真实 B 亚型已映射回 `integrated_oc`。当前 interaction_group 中 B 相关细胞数：

| interaction_group | n_cells |
|---|---:|
{table_lines(group_counts[group_counts['interaction_group'].astype(str).str.startswith('B')], ['interaction_group','n_cells'], 20)}

## 5. 功能注释主要结论

- EMT 模块最高：{top_module('EMT')}
- Hypoxia 模块最高：{top_module('Hypoxia')}
- KRAS activation 模块最高：{top_module('KRAS_activation')}
- LAPTM5 axis 最高：{top_module('LAPTM5_axis')}
- Immune-modulatory 模块最高：{top_module('Immune_modulatory')}
- Proliferation 模块最高：{top_module('Proliferation')}
- Stemness/epithelial 模块最高：{top_module('Stemness_epithelial')}

解释要点：
- `Subclone_02` 表现出较高 KRAS activation 和 epithelial/stemness 模块。
- `Subclone_04` 表现出较高 hypoxia、immune-modulatory、KRAS activation。
- `Subclone_01` 在 EMT、LAPTM5 axis、M2-related、proliferation 上相对更高，但其自身细胞数最大，后续需结合样本来源和 CNV 事件谨慎解释。

## 6. PCA/kNN 邻近富集

方法：PCA 前 30 维，k=30；对 neighbor label 进行 1000 次随机置换，计算 observed proportion、expected proportion、enrichment ratio、z-score 和 empirical P。

TME-only 邻近富集最靠前结果：

| cnv_subclone | neighbor_label | observed_prop | expected_prop | enrichment_ratio | z_score | p_greater |
|---|---|---:|---:|---:|---:|---:|
{table_lines(tme_small.rename(columns={'expected_prop_perm_mean':'expected_prop','empirical_p_greater':'p_greater'}), ['cnv_subclone','neighbor_label','observed_prop','expected_prop','enrichment_ratio','z_score','p_greater'], 20)}

解释要点：
- PCA/kNN 邻近首先显示 CNV clone 自身强聚集，这是 clone transcriptional-state coherence 的证据。
- 排除 CNV clone 自身后，`Subclone_01` 和 `Subclone_03` 对 broad `T cells` 有轻度富集。
- 本数据中 M2/CAF-like 邻近信号在 PCA 空间中没有成为最强富集项；因此不能写作“某 clone 空间邻近 M2/CAF”，更准确应写作“在表达空间 kNN 中未见强 M2/CAF 邻近富集，通讯层面存在 myeloid/B 相关 LR 信号”。

## 7. 配体-受体互作重点

按照文档列出的重点通路筛选，CNV subclone 相关强轴包括：

| focus_axis | source_group | target_group | ligand | receptor | score |
|---|---|---|---|---|---:|
{table_lines(lr_focus_small, ['focus_axis','source_group','target_group','ligand','receptor','score'], 20)}

去重后的 CNV subclone 相关强 LR 候选见：

- `source_outputs/lr_top_distinct_LR_pairs_involving_cnv_subclones.csv`
- `tables/focused_LR_axes_involving_cnv_subclones.csv`

主要模式：
- Tumor clone -> myeloid/DC/B：`MIF -> CD74` 是最强、最稳定的通讯轴。
- Myeloid/macrophage -> tumor clone：`TIMP1 -> CD63`、`C3 -> IFITM3`、`APOE -> LRP1` 是较突出的候选。
- Smooth muscle/CAF-like -> tumor clone：`CXCL12 -> CXCR4` 可作为 CAF-like 调控候选，但得分低于 MIF-CD74 主轴。

## 8. 可写入文章/汇报的候选模型

当前结果更支持如下表达：

> 基于 inferCNV，我们在 integrated_oc 的 malignant/Other 细胞中定义了 5 个 CNV subclone。不同 subclone 呈现不同的功能状态，其中 Subclone_02/04 更偏 KRAS/hypoxia/immune-modulatory，Subclone_01 更偏 EMT/LAPTM5/M2-related。PCA-kNN 分析显示 CNV subclone 在表达空间中存在明显自聚集，而 TME-only 邻近富集较弱，仅 Subclone_01/03 对 T cells 有轻度富集。Connectome-like ligand-receptor 分析提示不同 CNV subclone 与 myeloid/DC/B 细胞之间存在显著的 MIF-CD74 通讯轴，myeloid 细胞还可能通过 TIMP1-CD63、C3-IFITM3、APOE-LRP1 等轴影响肿瘤 CNV subclone。

## 9. 结果目录结构

- `tables/`：整理后的分析表格，包括 clone 组成、功能评分、kNN 邻近、重点 LR 轴。
- `figures/`：功能评分热图、clone 组成图、kNN 邻近热图。
- `source_outputs/`：此前 inferCNV subclone、CNV heatmap、LR 互作输出的原始结果副本。
- `scripts/`：可复现脚本。

## 10. 注意事项

- kNN 邻近是表达空间邻近，不是组织空间邻近。
- LR score 是表达层面的潜在通讯，不证明真实物理互作。
- `Subclone_05` 只有 96 个细胞，统计稳定性弱，建议作为低频 clone 谨慎解释。
- 如果后续有 patient/condition metadata，应进一步按 patient 内部复核 clone 和互作，避免样本差异误判为 clone 差异。
"""

(plan_dir / "integrated_oc_infercnv_plan_task_breakdown_and_results.md").write_text(report, encoding="utf-8")
(plan_dir / "integrated_oc_infercnv_plan_task_breakdown_and_results.txt").write_text(report, encoding="utf-8")
print(plan_dir / "integrated_oc_infercnv_plan_task_breakdown_and_results.md")
