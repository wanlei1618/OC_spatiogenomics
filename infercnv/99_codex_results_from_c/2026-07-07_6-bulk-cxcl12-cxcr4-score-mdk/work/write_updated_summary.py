from pathlib import Path
import pandas as pd

OUT = Path(r"D:\spatiogenomics_new")
T = OUT / "tables"
R = OUT / "reports"
R.mkdir(parents=True, exist_ok=True)


def csv_block(df, cols=None, n=12):
    if cols:
        df = df[cols]
    return df.head(n).to_csv(index=False).replace(",", " | ")


bulk_meta = pd.read_csv(T / "bulk_meta_cox_results.csv")
resist = pd.read_csv(T / "bulk_platinum_resistance_cluster_results.csv")
axes = pd.read_csv(T / "singlecell_LR_axis_hits_CXCL12_CXCR4_MDK_SDC4.csv")
sources = pd.read_csv(T / "integrated_oc_candidate_gene_cell_sources.csv")
clone = pd.read_csv(T / "infercnv_clone_candidate_gene_specificity_tests.csv")
external = pd.read_csv(T / "external_CPTAC_PDC_DepMap_HPA_candidate_filter.csv")

axis_summary = (
    axes.groupby("axis")
    .agg(
        directed_pairs=("axis", "size"),
        nonzero_pairs=("prob", lambda x: int((x > 0).sum())),
        significant_pairs=("significant", lambda x: int(x.astype(bool).sum())),
        max_prob=("prob", "max"),
    )
    .reset_index()
)

lines = []
lines.append("# Spatiogenomics candidate-axis analysis")
lines.append("")
lines.append("This report replaces the previous single-cell communication layer. CellChat 1.6.1 was installed locally after installing Rtools40, and the integrated OC single-cell object was rerun using the metadata `cell_type` field.")
lines.append("")
lines.append("Bulk data: six cohorts from the local `D:/OC_found/...` dataset folder. Single-cell data: `D:/OC_spatiogenomics/infercnv/integrated_oc.RData`, RNA assay normalized `data`, grouped by metadata `cell_type`.")
lines.append("")
lines.append("PFS uses the available `DFS/dfsstatus` columns. Platinum-resistance analysis uses the available clinical `cluster` field; interpret as cluster-associated platinum-resistance signal unless the original coding map is confirmed.")
lines.append("")
lines.append("## Bulk Meta Cox Highlights")
top = bulk_meta[bulk_meta["model"] == "univariate"].sort_values("p")
lines.append(csv_block(top, ["endpoint", "signature", "n_cohorts", "meta_hr", "meta_ci_low", "meta_ci_high", "p", "fdr"], 10))
lines.append("")
lines.append("## Platinum Cluster Highlights")
lines.append(csv_block(resist.sort_values("p"), ["cohort", "signature", "test", "effect", "p", "fdr"], 12))
lines.append("")
lines.append("## Integrated OC Cell Sources")
lines.append(csv_block(sources, ["gene", "top_cell_types"], 10))
lines.append("")
lines.append("## CellChat Cell-Type Communication")
lines.append("CellChat was run on 49,326 cells and 8 metadata cell_type groups. The human CellChatDB LR table was filtered to LR pairs with all ligand/receptor subunits present in RNA signaling genes, retaining 1,720 / 1,939 LR pairs. `computeCommunProb` used raw RNA signaling expression, population-size weighting, and 100 bootstraps.")
lines.append("")
lines.append(csv_block(axis_summary, ["axis", "directed_pairs", "nonzero_pairs", "significant_pairs", "max_prob"], 10))
lines.append("")
sig_axes = axes[axes["significant"].astype(bool)].sort_values("prob", ascending=False)
lines.append("Significant candidate-axis directed pairs:")
lines.append(csv_block(sig_axes, ["source", "target", "ligand", "receptor", "prob", "pval", "axis"], 20))
lines.append("")
lines.append("Interpretation: CXCL12-CXCR4 is supported from Smooth_muscle_cells to immune/Other targets. MDK-SDC4 was present in the CellChatDB and retained for calculation, but no cell_type-level directed pair had nonzero/significant CellChat probability under the triMean-based CellChat model; this is consistent with sparse SDC4 expression by broad cell_type.")
lines.append("")
lines.append("## InferCNV Clone Specificity")
lines.append(csv_block(clone.sort_values("kruskal_p"), ["gene", "max_clone", "min_clone", "max_mean", "min_mean", "delta_max_min", "kruskal_p", "fdr"], 10))
lines.append("")
lines.append("Clone-level expression remains strongly clone-specific for MDK, SDC4, CPT1A, CXCR4, and CXCL12. MDK is highest in Subclone_04; SDC4 is highest in Subclone_02; CXCR4 and CPT1A are highest in Subclone_04.")
lines.append("")
lines.append("## External Filter")
lines.append(csv_block(external, ["Gene", "priority", "filter_score", "filter_rationale", "singlecell_lr_support", "max_singlecell_lr_score"], 10))
lines.append("")
lines.append("HPA was queried from the official downloadable protein atlas table. PDC/CPTAC ovarian study and DepMap are retained as orthogonal filtering sources; DepMap numeric matrix download remained blocked by browser verification in this environment.")
lines.append("")
lines.append("## Main Outputs")
for p in sorted(T.glob("*.csv")):
    lines.append(f"- `{p}`")
for p in sorted(T.glob("*.rds")):
    lines.append(f"- `{p}`")
for p in sorted((OUT / "figures").glob("*.png")):
    lines.append(f"- `{p}`")
for p in sorted(R.glob("*.md")):
    if p.name != "summary_report.md":
        lines.append(f"- `{p}`")

(R / "summary_report.md").write_text("\n".join(lines), encoding="utf-8")
print("Updated summary_report.md")
