# C 盘 Codex infercnv 相关结果迁移说明

本目录用于集中保存从 `C:\Users\chenfy12\Documents\Codex` 中识别出的 infercnv / OC spatiogenomics 相关结果与辅助脚本。迁移后，后续梳理、验证和提交给 ChatGPT 的材料均以 `D:\OC_spatiogenomics\infercnv` 为统一工作位置。

## 已迁移来源

1. `2026-07-05_list-files-home-shpc-006-oc`
   - 迁移了该 Codex 会话下的 `outputs`。
   - 迁移了与 infercnv / OC spatiogenomics 分析有关的顶层脚本、安装脚本和工作文件。
   - 保留未迁移：`work\rpkgs` 等 R 包缓存目录，这些属于运行环境依赖，不是结果分析文件。
   - 2026-07-10 复核后，58 个主线目录缺失的 2026-07-05 文件已移动补入 `00`-`05` 和 `90` 分类目录；36 个与主线目录同名同大小的文件保留在本来源归档中，避免重复投放。

2. `2026-07-07_6-bulk-cxcl12-cxcr4-score-mdk`
   - 迁移了引用 `D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis` 的 CXCL12/CXCR4/MDK/SDC4 相关分析脚本。
   - 该会话的实际结果输出位置在 `D:\spatiogenomics_new`，不是 C 盘 Codex 结果，因此此处只归档 C 盘中的脚本/说明文件。
   - 保留未迁移：`work\r_pack*` 等 R 包缓存目录。

3. `2026-07-10_d-oc-spatiogenomics-infercnv-md`
   - 迁移了本次 infercnv 文件归类过程中曾生成在 C 盘 Codex 会话下的 `outputs`。
   - 其中包括前一版用于 ChatGPT 验证的梳理文件。

## 清单文件

详细迁移记录见：

`D:\OC_spatiogenomics\infercnv\99_codex_results_from_c\codex_c_drive_infercnv_migration_manifest.csv`

2026-07-05 结果复核与二次归类记录见：

`D:\OC_spatiogenomics\infercnv\99_codex_results_from_c\2026-07-05_reclassification_manifest.csv`

该表记录了原始位置、目标位置、迁移类别与状态，可用于追踪每个文件的来源。

## 迁移后复查

已对 `C:\Users\chenfy12\Documents\Codex` 进行两类复查：

1. 文件名/路径复查：未发现仍留在 C 盘 Codex 目录中的 infercnv、OC_spatiogenomics、sample_type、SPP1、CD44、ITGB1、CNV、spatiogenomics、ovarian 相关命中文件。
2. 文本内容复查：排除 `.git` 和 R 包缓存目录后，未发现仍留在 C 盘 Codex 目录中的上述关键词命中文件。

因此，当前与 infercnv 课题直接相关的 C 盘 Codex 结果已经统一归档到本目录。
