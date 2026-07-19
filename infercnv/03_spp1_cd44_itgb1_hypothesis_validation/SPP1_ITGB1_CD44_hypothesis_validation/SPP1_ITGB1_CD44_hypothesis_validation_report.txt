# SPP1-ITGB1/CD44 hypothesis validation

## Main conclusion
The current data support the malignant/tumor identity of CNV_Subclone_02/04 and provide focused LR evidence for a myeloid/macrophage SPP1 axis toward these clones. The strongest support is:

- CNV_Subclone_02/04 have higher epithelial/tumor scores and lower immune-contamination scores than the other CNV clones.
- ITGB1 is higher in the Subclone_02/04 focus group at the sample x clone level; CD44 is strongest in Subclone_04 rather than uniformly high in both Subclone_02 and Subclone_04.
- Existing Connectome-like focused LR results support SPP1->ITGB1 and SPP1->CD44 from Myeloid_Macro-Inflammatory_TNF, Myeloid_Interferon-Responsive Myeloid, and Myeloid_Macro-C3/CX3CR1 to CNV_Subclone_02/04.
- CellChat successfully reran on the selected sender/receiver groups and reproduced SPP1->CD44 to CNV_Subclone_04. It did not return a single SPP1->ITGB1 interaction in the exported CellChat table, likely because CellChat's SPP1 pathway encodes a narrower receptor set than the Connectome-like database.
- Pseudo-bulk DEG/GSEA supports KRAS signaling enrichment in Subclone_02/04 vs other clones, while curated ECM remodeling was negatively enriched. Therefore the SPP1-ITGB1/CD44 axis is a strong LR candidate, but the downstream ECM/FAK target-program evidence is not yet fully confirmed in this dataset.

## Output scope
This folder validates the hypothesis that myeloid/macrophage-derived SPP1 signals to CNV_Subclone_02/04 through ITGB1/CD44.

## 5.1 Target clone malignant identity and receptor robustness
- Tumor/receptor/immune marker summary: `tables/clone_marker_receptor_expression_summary.csv`.
- Sample/batch-stratified receptor and state scores: `tables/sample_batch_stratified_receptor_program_scores.csv`.
- Focus Subclone_02/04 vs other clone tests: `tables/focus_clone_02_04_vs_others_sample_stratified_tests.csv`.
- Main figures: `figures/clone_tumor_identity_receptor_dotplot.png`, `figures/sample_stratified_receptor_program_scores.png`.

## 5.2 Pseudo-bulk expression program validation
- Pseudo-bulk groups were sample/batch x cnv_subclone, using normalized RNA expression averages because raw counts are not fully exposed in the object.
- DEG table: `tables/pseudo_bulk_limma_focus_Subclone02_04_vs_others.csv`.
- Curated GSEA table: `tables/pseudo_bulk_focus_0204_GSEA_curated_pathways.csv`.

## 5.3 LR replication
- Existing Connectome-like focused LR results directly support SPP1->ITGB1/CD44 from myeloid/macrophage senders to CNV_Subclone_02/04.
- Focused predefined axis table: `tables/focused_LR_predefined_axes_senders_to_Subclone02_04.csv`.
- Ranked selected sender/receiver table: `tables/focused_LR_ranked_all_pairs_selected_senders_receivers.csv`.
- CellChat run status: `tables/CellChat_run_status.csv`.
- NicheNet package was not installed; a target-program enrichment proxy is provided in `tables/SPP1_axis_target_program_proxy_enrichment.csv`.

## 5.4 Bulk interaction model
- Bulk cohorts were not present locally, so the Cox interaction model was not executed.
- Required inputs: `tables/bulk_interaction_model_required_inputs.csv`.
- Template script: `scripts/bulk_interaction_cox_model_template.R`.

## 5.5 Spatial validation pre-analysis
- Marker panel: `tables/spatial_validation_marker_panel.csv`.
- Scoring strategy: `tables/spatial_validation_scoring_strategy.csv`.

## Environment caveats
- Missing marker genes in integrated_oc RNA matrix: none.
- CopyKAT/CaSpER/LIANA/NicheNet/DESeq2/edgeR were not installed in this local R 4.0.3 environment.
- The Seurat package itself is not loadable here, but the Seurat object's metadata and RNA data matrix were accessible through object slots.
