import os
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.api import Logit, add_constant
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


DATA_DIR = Path(r"D:\OC_found\数据集")
SC_DIR = Path(r"D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis")
OUT_DIR = Path(r"D:\spatiogenomics_new")
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"
REPORT_DIR = OUT_DIR / "reports"

for d in [OUT_DIR, FIG_DIR, TABLE_DIR, REPORT_DIR]:
    d.mkdir(parents=True, exist_ok=True)


SIGNATURES = {
    "CXCL12_CXCR4": ["CXCL12", "CXCR4"],
    "MDK_SDC4": ["MDK", "SDC4"],
    "FAO_CPT1A": [
        "CPT1A", "CPT1B", "CPT1C", "CPT2", "SLC25A20", "ACSL1", "ACSL3",
        "ACSL4", "ACSL5", "CD36", "ACADM", "ACADS", "ACADSB", "ACADVL",
        "HADHA", "HADHB", "ACOX1", "ETFA", "ETFB", "ETFDH", "PPARA",
        "PPARGC1A",
    ],
    "TLS": [
        "CXCL13", "CCL19", "CCL21", "CCR7", "LTA", "LTB", "MS4A1",
        "CD19", "CD79A", "BANK1", "CD3D", "CD3E", "CD4", "CD8A",
        "CD8B", "CXCR5", "ICOS", "PDCD1", "BCL6", "IL21", "CD40LG",
        "LAMP3", "POU2AF1",
    ],
    "CA_MSC": [
        "ACTA2", "FAP", "PDGFRA", "PDGFRB", "THY1", "PDPN", "POSTN",
        "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "VIM",
        "TAGLN", "MYL9", "CXCL12", "IL6", "TGFB1", "S100A4", "MMP2",
        "VCAN", "SPARC",
    ],
}


COHORTS = {
    "GSE10": ("GSE10_inter.txt", "GSE10_OS.csv"),
    "GSE14": ("GSE14_inter.txt", "GSE14_OS.csv"),
    "GSE32": ("GSE32_inter.txt", "GSE32_OS.csv"),
    "GSE49": ("GSE49_inter.txt", "GSE49_OS.csv"),
    "ICGC": ("ICGC_micro_inter.txt", "ICGC_OS.csv"),
    "TCGA": ("TCGA_inter.txt", "TCGA_OS.csv"),
}


def bh(df, p_col="p"):
    out = df.copy()
    mask = out[p_col].notna()
    out["fdr"] = np.nan
    if mask.any():
        out.loc[mask, "fdr"] = multipletests(out.loc[mask, p_col], method="fdr_bh")[1]
    return out


def read_expr(path):
    expr = pd.read_csv(path, sep="\t")
    gene_col = expr.columns[0]
    expr[gene_col] = expr[gene_col].astype(str).str.upper()
    expr = expr.groupby(gene_col).mean(numeric_only=True)
    return expr


def zscore_rows(mat):
    vals = mat.astype(float)
    mean = vals.mean(axis=1)
    sd = vals.std(axis=1, ddof=0).replace(0, np.nan)
    return vals.sub(mean, axis=0).div(sd, axis=0)


def score_signatures(expr):
    log_expr = np.log1p(expr)
    z = zscore_rows(log_expr)
    scores = {}
    coverage = []
    for sig, genes in SIGNATURES.items():
        present = [g for g in genes if g in z.index]
        missing = [g for g in genes if g not in z.index]
        if present:
            scores[sig] = z.loc[present].mean(axis=0)
        else:
            scores[sig] = pd.Series(np.nan, index=z.columns)
        coverage.append({
            "signature": sig,
            "n_genes": len(genes),
            "n_present": len(present),
            "present_genes": ";".join(present),
            "missing_genes": ";".join(missing),
        })
    return pd.DataFrame(scores), pd.DataFrame(coverage)


def clean_clinical(clin):
    clin = clin.copy()
    clin.index = clin.index.astype(str)
    for col in ["age", "OS", "osstatus", "DFS", "dfsstatus", "cluster"]:
        if col in clin.columns:
            clin[col] = pd.to_numeric(clin[col], errors="coerce")
    if "stage" in clin.columns:
        clin["stage"] = clin["stage"].astype(str).str.upper().str.extract(r"(II|III|IV|I)", expand=False)
    return clin


def cox_one(df, time_col, event_col, score_col, adjust=False):
    cols = [time_col, event_col, score_col]
    tmp = df[cols + (["age", "stage"] if adjust and "age" in df.columns and "stage" in df.columns else [])].copy()
    tmp = tmp.rename(columns={time_col: "time", event_col: "event", score_col: "score"})
    tmp = tmp.replace([np.inf, -np.inf], np.nan).dropna(subset=["time", "event", "score"])
    tmp = tmp[(tmp["time"] > 0) & (tmp["event"].isin([0, 1]))]
    if tmp.shape[0] < 20 or tmp["event"].sum() < 5:
        return None
    tmp["score"] = (tmp["score"] - tmp["score"].mean()) / tmp["score"].std(ddof=0)
    model_cols = ["time", "event", "score"]
    if adjust and "age" in tmp.columns and "stage" in tmp.columns:
        tmp["age"] = pd.to_numeric(tmp["age"], errors="coerce")
        stage_dum = pd.get_dummies(tmp["stage"], prefix="stage", drop_first=True)
        tmp = pd.concat([tmp.drop(columns=["stage"]), stage_dum], axis=1)
        model_cols = ["time", "event", "score", "age"] + list(stage_dum.columns)
        tmp = tmp.dropna(subset=model_cols)
    try:
        cph = CoxPHFitter()
        cph.fit(tmp[model_cols], duration_col="time", event_col="event")
        row = cph.summary.loc["score"]
        return {
            "n": int(tmp.shape[0]),
            "events": int(tmp["event"].sum()),
            "hr": float(row["exp(coef)"]),
            "ci_low": float(row["exp(coef) lower 95%"]),
            "ci_high": float(row["exp(coef) upper 95%"]),
            "coef": float(row["coef"]),
            "se": float(row["se(coef)"]),
            "p": float(row["p"]),
        }
    except Exception as exc:
        return {"n": int(tmp.shape[0]), "events": int(tmp["event"].sum()), "error": str(exc), "p": np.nan}


def logrank_one(df, time_col, event_col, score_col):
    tmp = df[[time_col, event_col, score_col]].copy()
    tmp = tmp.rename(columns={time_col: "time", event_col: "event", score_col: "score"})
    tmp = tmp.replace([np.inf, -np.inf], np.nan).dropna()
    tmp = tmp[(tmp["time"] > 0) & (tmp["event"].isin([0, 1]))]
    if tmp.shape[0] < 20:
        return None
    med = tmp["score"].median()
    hi = tmp["score"] >= med
    if hi.sum() < 5 or (~hi).sum() < 5:
        return None
    res = logrank_test(tmp.loc[hi, "time"], tmp.loc[~hi, "time"], tmp.loc[hi, "event"], tmp.loc[~hi, "event"])
    return {
        "n_high": int(hi.sum()),
        "n_low": int((~hi).sum()),
        "median_cut": float(med),
        "p": float(res.p_value),
    }


def resistance_tests(df, score_col):
    out = []
    tmp = df[["cluster", score_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if tmp["cluster"].nunique() >= 2:
        groups = [g[score_col].values for _, g in tmp.groupby("cluster")]
        if all(len(g) > 1 for g in groups):
            kw = stats.kruskal(*groups)
            out.append({
                "test": "Kruskal_cluster_all",
                "n": int(tmp.shape[0]),
                "effect": np.nan,
                "p": float(kw.pvalue),
            })
    bin_tmp = tmp[tmp["cluster"].isin([0, 2])].copy()
    if bin_tmp["cluster"].nunique() == 2 and min(bin_tmp["cluster"].value_counts()) >= 10:
        x0 = bin_tmp.loc[bin_tmp["cluster"] == 0, score_col]
        x2 = bin_tmp.loc[bin_tmp["cluster"] == 2, score_col]
        mw = stats.mannwhitneyu(x2, x0, alternative="two-sided")
        out.append({
            "test": "MannWhitney_cluster2_vs_0",
            "n": int(bin_tmp.shape[0]),
            "effect": float(x2.median() - x0.median()),
            "p": float(mw.pvalue),
        })
        try:
            design = add_constant(((bin_tmp[score_col] - bin_tmp[score_col].mean()) / bin_tmp[score_col].std(ddof=0)).rename("score"))
            y = (bin_tmp["cluster"] == 2).astype(int)
            fit = Logit(y, design).fit(disp=False)
            out.append({
                "test": "Logit_cluster2_vs_0",
                "n": int(bin_tmp.shape[0]),
                "effect": float(np.exp(fit.params["score"])),
                "p": float(fit.pvalues["score"]),
            })
        except Exception:
            pass
    return out


def meta_analyze(cox_df):
    rows = []
    for (endpoint, signature, model), sub in cox_df.dropna(subset=["coef", "se"]).groupby(["endpoint", "signature", "model"]):
        sub = sub[sub["se"] > 0]
        if sub.empty:
            continue
        w = 1 / (sub["se"] ** 2)
        coef = float((w * sub["coef"]).sum() / w.sum())
        se = float(math.sqrt(1 / w.sum()))
        z = coef / se
        p = float(2 * stats.norm.sf(abs(z)))
        rows.append({
            "endpoint": endpoint,
            "signature": signature,
            "model": model,
            "n_cohorts": int(sub.shape[0]),
            "meta_hr": float(np.exp(coef)),
            "meta_ci_low": float(np.exp(coef - 1.96 * se)),
            "meta_ci_high": float(np.exp(coef + 1.96 * se)),
            "meta_coef": coef,
            "meta_se": se,
            "p": p,
        })
    return bh(pd.DataFrame(rows), "p") if rows else pd.DataFrame()


def run_bulk():
    all_scores = []
    all_cov = []
    cox_rows = []
    lr_rows = []
    resist_rows = []
    for cohort, (expr_file, clin_file) in COHORTS.items():
        expr = read_expr(DATA_DIR / expr_file)
        clin = clean_clinical(pd.read_csv(DATA_DIR / clin_file, index_col=0))
        common = [s for s in expr.columns if s in clin.index]
        expr = expr[common]
        clin = clin.loc[common]
        scores, coverage = score_signatures(expr)
        coverage.insert(0, "cohort", cohort)
        coverage.insert(1, "n_samples", len(common))
        all_cov.append(coverage)
        merged = clin.join(scores)
        merged.insert(0, "sample_id", merged.index)
        merged.insert(0, "cohort", cohort)
        all_scores.append(merged)

        for sig in SIGNATURES:
            for endpoint, time_col, event_col in [("OS", "OS", "osstatus"), ("PFS_DFS", "DFS", "dfsstatus")]:
                for adjust in [False, True]:
                    res = cox_one(merged, time_col, event_col, sig, adjust=adjust)
                    if res:
                        res.update({"cohort": cohort, "endpoint": endpoint, "signature": sig, "model": "adjusted_age_stage" if adjust else "univariate"})
                        cox_rows.append(res)
                res = logrank_one(merged, time_col, event_col, sig)
                if res:
                    res.update({"cohort": cohort, "endpoint": endpoint, "signature": sig})
                    lr_rows.append(res)
            for rr in resistance_tests(merged, sig):
                rr.update({"cohort": cohort, "signature": sig})
                resist_rows.append(rr)

    score_df = pd.concat(all_scores, ignore_index=True)
    cov_df = pd.concat(all_cov, ignore_index=True)
    cox_df = bh(pd.DataFrame(cox_rows), "p")
    lr_df = bh(pd.DataFrame(lr_rows), "p")
    resist_df = bh(pd.DataFrame(resist_rows), "p")
    meta_df = meta_analyze(cox_df)

    score_df.to_csv(TABLE_DIR / "bulk_signature_scores_by_sample.csv", index=False, encoding="utf-8-sig")
    cov_df.to_csv(TABLE_DIR / "bulk_signature_gene_coverage.csv", index=False, encoding="utf-8-sig")
    cox_df.to_csv(TABLE_DIR / "bulk_cox_OS_PFS_results.csv", index=False, encoding="utf-8-sig")
    lr_df.to_csv(TABLE_DIR / "bulk_logrank_median_split_results.csv", index=False, encoding="utf-8-sig")
    resist_df.to_csv(TABLE_DIR / "bulk_platinum_resistance_cluster_results.csv", index=False, encoding="utf-8-sig")
    meta_df.to_csv(TABLE_DIR / "bulk_meta_cox_results.csv", index=False, encoding="utf-8-sig")
    return score_df, cov_df, cox_df, lr_df, resist_df, meta_df


def summarize_single_cell():
    meta_path = SC_DIR / "tables" / "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv"
    lr_path = SC_DIR / "source_outputs" / "lr_top_interactions_involving_cnv_subclones.csv"
    focused_path = SC_DIR / "tables" / "focused_LR_axes_involving_cnv_subclones.csv"
    clone_path = SC_DIR / "tables" / "cnv_subclone_key_gene_expression_summary.csv"
    module_clone_path = SC_DIR / "tables" / "cnv_subclone_function_module_score_summary.csv"

    meta = pd.read_csv(meta_path)
    cell_counts = meta["interaction_group"].value_counts(dropna=False).rename_axis("interaction_group").reset_index(name="n_cells")
    cell_counts.to_csv(TABLE_DIR / "integrated_oc_interaction_group_cell_counts.csv", index=False, encoding="utf-8-sig")

    lr = pd.read_csv(lr_path)
    lr["axis"] = lr["ligand"].astype(str).str.upper() + "_" + lr["receptor"].astype(str).str.upper()
    axes = ["CXCL12_CXCR4", "MDK_SDC4"]
    axis_hits = lr[lr["axis"].isin(axes)].sort_values(["axis", "score"], ascending=[True, False])
    axis_hits.to_csv(TABLE_DIR / "singlecell_LR_axis_hits_CXCL12_CXCR4_MDK_SDC4.csv", index=False, encoding="utf-8-sig")

    focused = pd.read_csv(focused_path) if focused_path.exists() else pd.DataFrame()
    if not focused.empty:
        focused.to_csv(TABLE_DIR / "CellChat_LIANA_like_focused_LR_axes.csv", index=False, encoding="utf-8-sig")

    clone_key = pd.read_csv(clone_path)
    module_clone = pd.read_csv(module_clone_path)
    clone_key.to_csv(TABLE_DIR / "infercnv_clone_key_gene_expression_summary_existing.csv", index=False, encoding="utf-8-sig")
    module_clone.to_csv(TABLE_DIR / "infercnv_clone_module_score_summary_existing.csv", index=False, encoding="utf-8-sig")

    source_rows = []
    for gene, role in [("CXCL12", "ligand"), ("CXCR4", "receptor"), ("MDK", "ligand"), ("SDC4", "receptor")]:
        rows = lr[(lr["ligand"].str.upper() == gene) | (lr["receptor"].str.upper() == gene)].copy()
        if rows.empty:
            source_rows.append({"gene": gene, "role": role, "top_groups": "", "note": "not detected in exported LR table"})
            continue
        group_col = "source_group" if role == "ligand" else "target_group"
        top = (
            rows.groupby(group_col)
            .agg(max_score=("score", "max"), mean_expr=("ligand_avg" if role == "ligand" else "receptor_avg", "max"), pct=("ligand_pct" if role == "ligand" else "receptor_pct", "max"))
            .sort_values(["max_score", "mean_expr"], ascending=False)
            .head(8)
            .reset_index()
        )
        source_rows.append({
            "gene": gene,
            "role": role,
            "top_groups": "; ".join([f"{r[group_col]}(score={r.max_score:.3g},expr={r.mean_expr:.3g},pct={r.pct:.2g})" for _, r in top.iterrows()]),
            "note": "inferred from exported ligand-receptor averages",
        })
    pd.DataFrame(source_rows).to_csv(TABLE_DIR / "integrated_oc_candidate_gene_cell_sources.csv", index=False, encoding="utf-8-sig")
    return axis_hits, clone_key, module_clone


def make_plots(cox_df, resist_df, meta_df, axis_hits):
    sns.set(style="whitegrid", font_scale=0.9)
    uni = cox_df[(cox_df["model"] == "univariate") & cox_df["hr"].notna()].copy()
    if not uni.empty:
        uni["log10_fdr"] = -np.log10(uni["fdr"].clip(lower=1e-300))
        for endpoint in ["OS", "PFS_DFS"]:
            sub = uni[uni["endpoint"] == endpoint]
            if sub.empty:
                continue
            pivot = sub.pivot_table(index="signature", columns="cohort", values="hr")
            plt.figure(figsize=(7.2, 3.8))
            sns.heatmap(np.log2(pivot), center=0, cmap="vlag", annot=pivot.round(2), fmt="", linewidths=.5)
            plt.title(f"{endpoint} Cox HR per SD score (log2 color, HR labels)")
            plt.tight_layout()
            plt.savefig(FIG_DIR / f"bulk_{endpoint}_cox_hr_heatmap.png", dpi=220)
            plt.close()

    if not meta_df.empty:
        sub = meta_df[meta_df["model"] == "univariate"].sort_values(["endpoint", "meta_hr"])
        plt.figure(figsize=(7, 4.8))
        y = np.arange(sub.shape[0])
        plt.errorbar(sub["meta_hr"], y, xerr=[sub["meta_hr"] - sub["meta_ci_low"], sub["meta_ci_high"] - sub["meta_hr"]], fmt="o", color="#2f5d7c")
        plt.axvline(1, color="grey", linestyle="--", linewidth=1)
        plt.yticks(y, sub["endpoint"] + " | " + sub["signature"])
        plt.xscale("log")
        plt.xlabel("Meta-analysis HR per SD")
        plt.tight_layout()
        plt.savefig(FIG_DIR / "bulk_meta_cox_forest.png", dpi=220)
        plt.close()

    if not resist_df.empty:
        sub = resist_df[resist_df["test"] == "MannWhitney_cluster2_vs_0"]
        if not sub.empty:
            pivot = sub.pivot_table(index="signature", columns="cohort", values="effect")
            plt.figure(figsize=(7.2, 3.8))
            sns.heatmap(pivot, center=0, cmap="vlag", annot=True, fmt=".2f", linewidths=.5)
            plt.title("Score median difference: cluster 2 - cluster 0")
            plt.tight_layout()
            plt.savefig(FIG_DIR / "bulk_platinum_cluster_score_difference_heatmap.png", dpi=220)
            plt.close()

    if axis_hits is not None and not axis_hits.empty:
        top = axis_hits.sort_values("score", ascending=False).head(25).copy()
        top["pair"] = top["source_group"] + " -> " + top["target_group"] + " | " + top["axis"]
        plt.figure(figsize=(8.5, max(4, 0.25 * len(top))))
        sns.barplot(data=top, y="pair", x="score", hue="axis", dodge=False)
        plt.xlabel("LR score")
        plt.ylabel("")
        plt.tight_layout()
        plt.savefig(FIG_DIR / "singlecell_candidate_LR_axis_top_hits.png", dpi=220)
        plt.close()


def write_report(cox_df, lr_df, resist_df, meta_df, axis_hits):
    def md_table(df):
        return df.to_csv(index=False).replace(",", " | ")

    lines = []
    lines.append("# Spatiogenomics candidate-axis analysis\n")
    lines.append("Input data: six bulk cohorts from the local `D:/OC_found/...` dataset folder; integrated OC / inferCNV exported tables from `D:/OC_spatiogenomics/infercnv/integrated_oc_plan_analysis`.\n")
    lines.append("Scores: ligand-receptor axes are mean z-scores of the two genes; FAO/CPT1A, TLS, and CA-MSC are mean z-scores of curated gene sets. PFS uses the available `DFS/dfsstatus` columns.\n")
    lines.append("Platinum-resistance readout: the available clinical field is `cluster`; analyses include all-cluster Kruskal-Wallis and `cluster=2` vs `cluster=0` binary tests. Interpret this as cluster-associated platinum-resistance signal unless the original coding map is confirmed.\n")
    lines.append("Communication-axis validation: CellChat/LIANA packages were not installed in the local R 4.0.3 environment, so this run uses the exported integrated-OC ligand-receptor tables from the existing analysis directory as a CellChat/LIANA-like validation layer.\n")
    if not meta_df.empty:
        lines.append("## Bulk meta Cox highlights\n")
        top = meta_df[meta_df["model"] == "univariate"].sort_values("p").head(12)
        lines.append(md_table(top[["endpoint", "signature", "n_cohorts", "meta_hr", "meta_ci_low", "meta_ci_high", "p", "fdr"]]))
        lines.append("")
    if not resist_df.empty:
        lines.append("## Platinum cluster highlights\n")
        top = resist_df.sort_values("p").head(12)
        lines.append(md_table(top[["cohort", "signature", "test", "effect", "p", "fdr"]]))
        lines.append("")
    if axis_hits is not None and not axis_hits.empty:
        lines.append("## Integrated OC communication-axis hits\n")
        lines.append(md_table(axis_hits.sort_values("score", ascending=False).head(15)))
        lines.append("")
    lines.append("## Main outputs\n")
    for p in sorted(TABLE_DIR.glob("*.csv")):
        lines.append(f"- `{p}`")
    for p in sorted(FIG_DIR.glob("*.png")):
        lines.append(f"- `{p}`")
    lines.append("- External filter report: `D:/spatiogenomics_new/reports/external_candidate_filter.md`")
    (REPORT_DIR / "summary_report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    (TABLE_DIR / "signature_definitions.json").write_text(json.dumps(SIGNATURES, indent=2), encoding="utf-8")
    score_df, cov_df, cox_df, lr_df, resist_df, meta_df = run_bulk()
    axis_hits, clone_key, module_clone = summarize_single_cell()
    make_plots(cox_df, resist_df, meta_df, axis_hits)
    write_report(cox_df, lr_df, resist_df, meta_df, axis_hits)
    print(f"Finished. Outputs: {OUT_DIR}")


if __name__ == "__main__":
    main()
