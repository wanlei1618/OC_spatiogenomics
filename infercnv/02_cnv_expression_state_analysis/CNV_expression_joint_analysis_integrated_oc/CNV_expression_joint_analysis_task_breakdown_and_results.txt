# CNV + 表达联合分析：integrated_oc

结果目录：`D:\OC_spatiogenomics\infercnv\CNV_expression_joint_analysis_integrated_oc`

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
| Subclone_04 | 0.584 | 0.5864 | 0.04111 |
| Subclone_02 | 0.535 | 0.5679 | 0.04295 |
| Subclone_05 | 0.5014 | 0.4854 | 0.02784 |
| Subclone_03 | 0.4816 | 0.509 | 0.0672 |
| Subclone_01 | 0.2564 | 0.2755 | 0.1454 |

结论：`Subclone_04` 的 CNV burden 最高，其次为 `Subclone_02/05/03`；`Subclone_01` 最低。

## 3. 功能 signature 与 TF target proxy activity

| feature | clone ranking by mean score |
|---|---|
| EMT | Subclone_01=0.409; Subclone_05=0.305; Subclone_02=0.292; Subclone_04=0.267; Subclone_03=0.192 |
| Hypoxia | Subclone_04=0.494; Subclone_03=0.473; Subclone_02=0.467; Subclone_01=0.344; Subclone_05=0.294 |
| KRAS_UP | Subclone_02=0.867; Subclone_04=0.773; Subclone_05=0.514; Subclone_01=0.454; Subclone_03=0.425 |
| KRAS_DN | Subclone_02=1.349; Subclone_04=1.267; Subclone_01=0.725; Subclone_05=0.676; Subclone_03=0.505 |
| LAPTM5_axis | Subclone_01=0.316; Subclone_04=0.107; Subclone_05=0.105; Subclone_02=0.053; Subclone_03=0.043 |
| Proliferation | Subclone_01=0.161; Subclone_03=0.151; Subclone_04=0.115; Subclone_02=0.110; Subclone_05=0.083 |
| Immune_modulatory | Subclone_04=0.984; Subclone_01=0.943; Subclone_02=0.875; Subclone_05=0.798; Subclone_03=0.684 |
| TF_SNAI2 | Subclone_01=0.415; Subclone_04=0.160; Subclone_02=0.144; Subclone_05=0.133; Subclone_03=0.088 |
| TF_ATF3 | Subclone_02=2.468; Subclone_04=2.219; Subclone_05=1.382; Subclone_01=1.220; Subclone_03=1.018 |
| TF_HIF1A | Subclone_04=0.351; Subclone_02=0.331; Subclone_03=0.297; Subclone_01=0.222; Subclone_05=0.182 |
| TF_STAT3 | Subclone_02=0.742; Subclone_04=0.721; Subclone_01=0.472; Subclone_05=0.414; Subclone_03=0.317 |
| TF_MYC | Subclone_02=0.245; Subclone_03=0.228; Subclone_04=0.221; Subclone_01=0.181; Subclone_05=0.136 |
| TF_FOXM1 | Subclone_01=0.139; Subclone_03=0.097; Subclone_04=0.073; Subclone_05=0.060; Subclone_02=0.050 |

解释要点：
- `Subclone_04`：CNV burden 最高，同时 Hypoxia、Immune-modulatory、TF_HIF1A、TF_STAT3 较高。
- `Subclone_02`：KRAS_UP、TF_ATF3、TF_STAT3、TF_MYC 较高，呈 KRAS/ATF3/STAT3-active 表型。
- `Subclone_01`：CNV burden 最低，但 EMT、LAPTM5_axis、TF_SNAI2、TF_FOXM1 相对更高，提示它更像 EMT/LAPTM5-high transcriptional-state clone，而不是 high-CNV clone。
- `Subclone_05`：细胞数少，仅 96 个，所有解释应谨慎。

## 4. clone marker genes

| CNV_clone | top marker genes |
|---|---|
| Subclone_01 | VIM, TMSB4X, GMFG, COTL1, ETS1, VAMP5, IL32, LGALS1, COL4A1, SRGN |
| Subclone_02 | MSLN, BCAM, CLDN4, CD9, SOX9, FOLR1, RARRES1, COL9A3, KLK6, ELF3 |
| Subclone_03 | GSTP1, RPL7, RPL10A, RPL9, RPS2, ATP5MC3, UBA52, RPS5, HSPE1, RPS16 |
| Subclone_04 | CYP4B1, CLDN7, EPCAM, SLC7A2, PERP, SPINT2, CLDN3, USP53, DSC2, CHI3L1 |
| Subclone_05 | MSLN, SCGB3A1, KRTCAP2, CLU, MT-ND4L, FKBP2, MT-ND4, MT-ND3, MT-ND5, MT-CO2 |

完整 marker 表见：

- `tables/CNV_clone_DEG_single_cell_wilcoxon_one_vs_rest.csv`
- `tables/CNV_clone_top20_markers_per_clone.csv`
- `figures/CNV_clone_top_marker_average_expression_heatmap.png`

## 5. CNV-expression dosage genes

基于 inferCNV subcluster 层面的 CNV HMM state 与平均表达做 Spearman correlation。正相关且 BH 校正显著的前 30 个候选：

| gene | rho | padj | n_groups |
|---|---:|---:|---:|
| MSLN | 0.8229 | 3.838e-10 | 52 |
| COX20 | 0.7976 | 2.792e-09 | 52 |
| S100A13 | 0.797 | 2.792e-09 | 52 |
| NDUFB9 | 0.7803 | 1.229e-08 | 52 |
| TSTD1 | 0.7769 | 1.373e-08 | 52 |
| WFDC2 | 0.7667 | 3.01e-08 | 52 |
| ADSS | 0.7653 | 3.01e-08 | 52 |
| C1orf35 | 0.7614 | 3.782e-08 | 52 |
| RPS5 | 0.7548 | 5.558e-08 | 52 |
| SLC39A4 | 0.754 | 5.558e-08 | 52 |
| NDUFS8 | 0.7536 | 5.558e-08 | 52 |
| COX6C | 0.7505 | 6.695e-08 | 52 |
| TACSTD2 | 0.7458 | 8.903e-08 | 52 |
| ELF3 | 0.7454 | 8.903e-08 | 52 |
| KCNK15 | 0.7398 | 1.323e-07 | 52 |
| AAAS | 0.7097 | 1.292e-06 | 52 |
| SCX | 0.7046 | 1.702e-06 | 52 |
| PRKRA | 0.7041 | 1.702e-06 | 52 |
| MVB12A | 0.7034 | 1.702e-06 | 52 |
| SDHC | 0.7005 | 1.969e-06 | 52 |
| BTG2 | 0.6997 | 1.969e-06 | 52 |
| ATPAF1 | 0.6993 | 1.969e-06 | 52 |
| COX4I1 | 0.698 | 2.063e-06 | 52 |
| OPN3 | 0.6969 | 2.127e-06 | 52 |
| LY6E | 0.6948 | 2.368e-06 | 52 |
| SLPI | 0.6917 | 2.803e-06 | 52 |
| HMGA1 | 0.6911 | 2.806e-06 | 52 |
| GNG5 | 0.6901 | 2.902e-06 | 52 |
| KRT8 | 0.6862 | 3.627e-06 | 52 |
| PPCS | 0.6833 | 4.249e-06 | 52 |

解释要点：
- 高相关 dosage candidates 包括 `MSLN`, `WFDC2`, `TACSTD2`, `ELF3`, `KRT8`, `SLPI` 等卵巢癌/上皮相关基因。
- 这些基因更适合作为 “CNV-dosage coupled expression genes” 候选，而不是直接等同于 driver gene。

## 6. CNV-transcriptome coupling

最强相关 pair：

| var1 | var2 | rho | padj |
|---|---|---:|---:|
| LAPTM5_axis | LAPTM5_expr | 1 | 0 |
| KRAS_UP | TF_ATF3 | 0.9287 | 0 |
| KRAS_DN | TF_ATF3 | 0.8376 | 0 |
| KRAS_DN | TF_STAT3 | 0.8222 | 0 |
| Hypoxia | TF_HIF1A | 0.8194 | 0 |
| KRAS_UP | KRAS_DN | 0.7628 | 0 |
| Proliferation | TF_MYC | 0.7328 | 0 |
| EMT | TF_SNAI2 | 0.6984 | 0 |
| TF_ATF3 | TF_STAT3 | 0.6792 | 0 |
| KRAS_UP | TF_STAT3 | 0.6625 | 0 |
| Proliferation | TF_FOXM1 | 0.5998 | 0 |
| Immune_modulatory | TF_STAT3 | 0.4971 | 0 |
| KRAS_DN | TF_RELA | 0.4887 | 0 |
| TF_STAT3 | TF_RELA | 0.4845 | 0 |
| Immune_modulatory | TF_SNAI2 | 0.4815 | 0 |
| Immune_modulatory | TF_RELA | 0.4745 | 0 |
| KRAS_DN | Immune_modulatory | 0.453 | 0 |
| LAPTM5_axis | Immune_modulatory | 0.4343 | 9.383e-297 |
| LAPTM5_expr | Immune_modulatory | 0.4343 | 9.383e-297 |
| TF_ATF3 | TF_RELA | 0.3833 | 5.982e-226 |
| LAPTM5_axis | TF_SNAI2 | 0.3676 | 1.209e-206 |
| LAPTM5_expr | TF_SNAI2 | 0.3676 | 1.209e-206 |
| CNV_burden | TF_SNAI2 | -0.3604 | 4.143e-198 |
| LAPTM5_axis | TF_RELA | 0.3588 | 2.642e-196 |
| LAPTM5_expr | TF_RELA | 0.3588 | 2.642e-196 |
| KRAS_UP | TF_RELA | 0.3574 | 1.031e-194 |
| EMT | Immune_modulatory | 0.3519 | 2.024e-188 |
| Hypoxia | TF_MYC | 0.3333 | 5.551e-168 |
| Immune_modulatory | TF_ATF3 | 0.3273 | 1.179e-161 |
| TF_SNAI2 | TF_RELA | 0.327 | 2.09e-161 |

和 CNV_burden 直接相关的 pair 可在 `tables/CNV_transcriptome_coupling_spearman_pairs.csv` 中筛选 `var1 == CNV_burden` 或 `var2 == CNV_burden`。

## 7. pseudo-bulk 设计

sample × CNV_clone pseudo-bulk 样本数：13

| pb_group | sample_id | CNV_clone |
|---|---|---|
| 31V1|Subclone_01 | 31V1 | Subclone_01 |
| 31V1|Subclone_03 | 31V1 | Subclone_03 |
| 31V1|Subclone_02 | 31V1 | Subclone_02 |
| 31V1|Subclone_04 | 31V1 | Subclone_04 |
| 31V1|Subclone_05 | 31V1 | Subclone_05 |
| 30V1|Subclone_01 | 30V1 | Subclone_01 |
| 30V1|Subclone_04 | 30V1 | Subclone_04 |
| 30V1|Subclone_03 | 30V1 | Subclone_03 |
| 30V1|Subclone_02 | 30V1 | Subclone_02 |
| 30V1|Subclone_05 | 30V1 | Subclone_05 |
| 30V2|Subclone_01 | 30V2 | Subclone_01 |
| 30V2|Subclone_02 | 30V2 | Subclone_02 |
| 30V3|Subclone_01 | 30V3 | Subclone_01 |

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
