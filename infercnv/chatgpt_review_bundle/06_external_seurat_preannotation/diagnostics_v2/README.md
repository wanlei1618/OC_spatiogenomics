# diagnostics_v2 review package

This small package contains forensic audits, mitochondrial-QC decisions, sample-dominance tables, A/B/C strategy summaries, RNA marker tables, plots, logs, and blank manual annotation templates.

Start with [`final_diagnostic_report.md`](final_diagnostic_report.md), then inspect:

- `00_forensic/current_result_snapshot.md`
- `GSE147082/01_mt_audit/mt_qc_decision.md`
- `GSE158722/01_mt_audit/mt_qc_decision.md`
- `GSE154600/02_dominance/cluster_dominance_diagnostic_table.csv`
- `GSE158722/02_dominance/cluster_dominance_diagnostic_table.csv`
- `GSE158722/02_dominance/patient_timepoint_cluster_correspondence.csv`
- `strategy_comparison/batch_strategy_decision_report.md`
- each dataset's `05_markers/top20_markers_per_cluster.csv` and blank `manual_annotation_template.csv`.

Large RDS files, count matrices, per-cell embeddings, and per-cell QC/lineage assignments are intentionally local-only.
