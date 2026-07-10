# integrated_oc 基于 inferCNV 的 CNV 亚克隆与 TME 互作分析结果

生成目录：`D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis`

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
| Subclone_01 | 3318 | 0.510 |
| Subclone_02 | 1071 | 0.165 |
| Subclone_03 | 1021 | 0.157 |
| Subclone_04 | 1001 | 0.154 |
| Subclone_05 | 96 | 0.015 |

## 3. integrated_oc 映射概况

| item | value |
|---|---|
| integrated_oc cells | 49326 |
| cnv subclone mapped | 6507 |
| T/NK subtype mapped | 15174 |
| myeloid subtype mapped | 13148 |
| B subtype source | D:/OC_spatiogenomics/infercnv/integratedocBcells.RData |
| B subtype mapped | 2274 |
| groups >=20 cells | 44 |

## 4. B 细胞亚型映射

`integratedocBcells.RData` 中真实 B 亚型已映射回 `integrated_oc`。当前 interaction_group 中 B 相关细胞数：

| interaction_group | n_cells |
|---|---:|
| B_Bn_TCL1A | 1183 |
| B cells | 659 |
| B_Classical-Bm_TXNIP | 633 |
| B_PC_IGHG | 281 |
| B_Bm_stress-response | 110 |
| B_Early-PC_MS4A1low | 67 |

## 5. 功能注释主要结论

- EMT 模块最高：Subclone_01=0.409; Subclone_05=0.305; Subclone_02=0.292
- Hypoxia 模块最高：Subclone_04=0.494; Subclone_03=0.473; Subclone_02=0.467
- KRAS activation 模块最高：Subclone_02=1.034; Subclone_04=0.921; Subclone_05=0.617
- LAPTM5 axis 最高：Subclone_01=0.316; Subclone_04=0.107; Subclone_05=0.105
- Immune-modulatory 模块最高：Subclone_04=0.984; Subclone_01=0.943; Subclone_02=0.875
- Proliferation 模块最高：Subclone_01=0.161; Subclone_03=0.151; Subclone_04=0.115
- Stemness/epithelial 模块最高：Subclone_02=0.990; Subclone_04=0.914; Subclone_05=0.800

解释要点：
- `Subclone_02` 表现出较高 KRAS activation 和 epithelial/stemness 模块。
- `Subclone_04` 表现出较高 hypoxia、immune-modulatory、KRAS activation。
- `Subclone_01` 在 EMT、LAPTM5 axis、M2-related、proliferation 上相对更高，但其自身细胞数最大，后续需结合样本来源和 CNV 事件谨慎解释。

## 6. PCA/kNN 邻近富集

方法：PCA 前 30 维，k=30；对 neighbor label 进行 1000 次随机置换，计算 observed proportion、expected proportion、enrichment ratio、z-score 和 empirical P。

TME-only 邻近富集最靠前结果：

| cnv_subclone | neighbor_label | observed_prop | expected_prop | enrichment_ratio | z_score | p_greater |
|---|---|---:|---:|---:|---:|---:|
| Subclone_01 | T cells | 0.03793 | 0.02408 | 1.575 | 5.673 | 0.000999 |
| Subclone_03 | T cells | 0.03399 | 0.024 | 1.416 | 3.148 | 0.001998 |
| Subclone_01 | T_NK_CD4+ T naive | 0.03029 | 0.02427 | 1.248 | 2.355 | 0.008991 |
| Subclone_04 | T cells | 0.02531 | 0.02405 | 1.052 | 0.4243 | 0.3167 |
| Subclone_03 | T_NK_CD4+ T naive | 0.01867 | 0.02431 | 0.7682 | -1.815 | 0.968 |
| Subclone_05 | B_PC_IGHG | 0.002083 | 0.005606 | 0.3716 | -1.569 | 0.975 |
| Subclone_04 | T_NK_CD4+ T regulatory | 0.005761 | 0.009175 | 0.6279 | -1.808 | 0.98 |
| Subclone_05 | T cells | 0.01389 | 0.02408 | 0.5768 | -2.24 | 0.996 |
| Subclone_04 | T_NK_CD4+ T naive | 0.01745 | 0.02435 | 0.7167 | -2.282 | 0.996 |
| Subclone_05 | B_Early-PC_MS4A1low | 0 | 0.001318 | 0 | -1.183 | 1 |
| Subclone_05 | NK cells | 0 | 0.001541 | 0 | -1.243 | 1 |
| Subclone_05 | Myeloid_mDC_LAMP3 | 0 | 0.001933 | 0 | -1.46 | 1 |
| Subclone_02 | B_Early-PC_MS4A1low | 0 | 0.001329 | 0 | -1.562 | 1 |
| Subclone_05 | B_Bm_stress-response | 0 | 0.002206 | 0 | -1.636 | 1 |
| Subclone_02 | NK cells | 0 | 0.001566 | 0 | -1.64 | 1 |
| Subclone_02 | Myeloid_mDC_LAMP3 | 0 | 0.001936 | 0 | -1.855 | 1 |
| Subclone_03 | B_Early-PC_MS4A1low | 0 | 0.001353 | 0 | -1.88 | 1 |
| Subclone_04 | B_Early-PC_MS4A1low | 0 | 0.001359 | 0 | -1.885 | 1 |
| Subclone_03 | NK cells | 0 | 0.001538 | 0 | -1.941 | 1 |
| Subclone_04 | NK cells | 0 | 0.001564 | 0 | -1.985 | 1 |

解释要点：
- PCA/kNN 邻近首先显示 CNV clone 自身强聚集，这是 clone transcriptional-state coherence 的证据。
- 排除 CNV clone 自身后，`Subclone_01` 和 `Subclone_03` 对 broad `T cells` 有轻度富集。
- 本数据中 M2/CAF-like 邻近信号在 PCA 空间中没有成为最强富集项；因此不能写作“某 clone 空间邻近 M2/CAF”，更准确应写作“在表达空间 kNN 中未见强 M2/CAF 邻近富集，通讯层面存在 myeloid/B 相关 LR 信号”。

## 7. 配体-受体互作重点

按照文档列出的重点通路筛选，CNV subclone 相关强轴包括：

| focus_axis | source_group | target_group | ligand | receptor | score |
|---|---|---|---|---|---:|
| MIF_CD74_CXCR4 | CNV_Subclone_03 | Myeloid_cDC1 | MIF | CD74 | 19.718 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | Myeloid_cDC2-like | MIF | CD74 | 18.151 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | Myeloid_cDC2 | MIF | CD74 | 17.817 |
| MIF_CD74_CXCR4 | CNV_Subclone_04 | Myeloid_cDC1 | MIF | CD74 | 17.518 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | Myeloid_mDC_LAMP3 | MIF | CD74 | 17.386 |
| MIF_CD74_CXCR4 | CNV_Subclone_02 | Myeloid_cDC1 | MIF | CD74 | 16.863 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | DC | MIF | CD74 | 16.666 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | B_Bn_TCL1A | MIF | CD74 | 16.369 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | Myeloid_Macro-C3/CX3CR1 | MIF | CD74 | 16.334 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | B_Classical-Bm_TXNIP | MIF | CD74 | 16.277 |
| MIF_CD74_CXCR4 | CNV_Subclone_04 | Myeloid_cDC2-like | MIF | CD74 | 16.126 |
| MIF_CD74_CXCR4 | CNV_Subclone_04 | Myeloid_cDC2 | MIF | CD74 | 15.828 |
| MIF_CD74_CXCR4 | CNV_Subclone_05 | Myeloid_cDC1 | MIF | CD74 | 15.758 |
| MIF_CD74_CXCR4 | CNV_Subclone_02 | Myeloid_cDC2-like | MIF | CD74 | 15.523 |
| MIF_CD74_CXCR4 | CNV_Subclone_04 | Myeloid_mDC_LAMP3 | MIF | CD74 | 15.445 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | B_Bm_stress-response | MIF | CD74 | 15.307 |
| MIF_CD74_CXCR4 | CNV_Subclone_02 | Myeloid_cDC2 | MIF | CD74 | 15.237 |
| MIF_CD74_CXCR4 | CNV_Subclone_02 | Myeloid_mDC_LAMP3 | MIF | CD74 | 14.868 |
| MIF_CD74_CXCR4 | CNV_Subclone_04 | DC | MIF | CD74 | 14.806 |
| MIF_CD74_CXCR4 | CNV_Subclone_03 | Myeloid_Interferon-Responsive Myeloid | MIF | CD74 | 14.774 |

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
