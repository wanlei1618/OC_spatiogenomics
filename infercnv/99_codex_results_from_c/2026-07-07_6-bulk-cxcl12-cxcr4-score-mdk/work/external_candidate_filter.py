import io
import zipfile
import urllib.request
from pathlib import Path

import pandas as pd


OUT_DIR = Path(r"D:\spatiogenomics_new")
TABLE_DIR = OUT_DIR / "tables"
REPORT_DIR = OUT_DIR / "reports"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

CANDIDATES = ["CXCL12", "CXCR4", "MDK", "SDC4", "CPT1A"]
HPA_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
PDC_STUDY_URL = "https://pdc.cancer.gov/pdc/study/PDC000110"
DEPMAP_DOWNLOAD_URL = "https://depmap.org/portal/download/all/"


def load_hpa_candidates():
    data = urllib.request.urlopen(HPA_URL, timeout=60).read()
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        with z.open("proteinatlas.tsv") as f:
            cols = [
                "Gene", "Ensembl", "Gene description", "Uniprot", "Protein class",
                "Evidence", "HPA evidence", "RNA tissue specificity",
                "RNA tissue distribution", "RNA cancer specificity",
                "RNA cancer distribution", "Protein tissue specificity",
                "Protein tissue distribution", "Secretome location",
                "Secretome function", "Subcellular location",
                "Cancer prognostics - Ovary Serous Cystadenocarcinoma (TCGA)",
                "Cancer prognostics - Ovary Serous Cystadenocarcinoma (validation)",
            ]
            df = pd.read_csv(f, sep="\t", usecols=lambda c: c in cols)
    return df[df["Gene"].isin(CANDIDATES)].copy()


def local_evidence():
    meta = pd.read_csv(TABLE_DIR / "bulk_meta_cox_results.csv")
    lr = pd.read_csv(TABLE_DIR / "singlecell_LR_axis_hits_CXCL12_CXCR4_MDK_SDC4.csv")
    sources = pd.read_csv(TABLE_DIR / "integrated_oc_candidate_gene_cell_sources.csv")
    resist = pd.read_csv(TABLE_DIR / "bulk_platinum_resistance_cluster_results.csv")

    rows = []
    for gene in CANDIDATES:
        sigs = []
        if gene in ["CXCL12", "CXCR4"]:
            sigs.append("CXCL12_CXCR4")
        if gene in ["MDK", "SDC4"]:
            sigs.append("MDK_SDC4")
        if gene == "CPT1A":
            sigs.append("FAO_CPT1A")
        meta_sub = meta[(meta["signature"].isin(sigs)) & (meta["model"] == "univariate")]
        res_sub = resist[resist["signature"].isin(sigs)].sort_values("p").head(3)
        source_sub = sources[sources["gene"] == gene]
        lr_sub = lr[(lr.get("ligand", pd.Series(dtype=str)).astype(str).str.upper() == gene) |
                    (lr.get("receptor", pd.Series(dtype=str)).astype(str).str.upper() == gene)]
        source_field = ""
        if not source_sub.empty:
            source_field = source_sub["top_cell_types"].iloc[0] if "top_cell_types" in source_sub.columns else source_sub.iloc[0].get("top_groups", "")
        score_col = "prob" if "prob" in lr_sub.columns else "score"
        if "significant" in lr_sub.columns:
            lr_support = int(((lr_sub[score_col] > 0) & (lr_sub["significant"].astype(bool))).sum())
        else:
            lr_support = int((lr_sub[score_col] > 0).sum()) if score_col in lr_sub.columns else int(lr_sub.shape[0])
        rows.append({
            "Gene": gene,
            "bulk_meta_signal": "; ".join(
                f"{r.endpoint}:{r.signature} HR={r.meta_hr:.2f},FDR={r.fdr:.3g}"
                for _, r in meta_sub.sort_values("p").iterrows()
            ),
            "best_platinum_cluster_signal": "; ".join(
                f"{r.cohort}:{r.test} p={r.p:.2g},effect={r.effect if pd.notna(r.effect) else 'NA'}"
                for _, r in res_sub.iterrows()
            ),
            "integrated_oc_cell_source": source_field,
            "singlecell_lr_support": lr_support,
            "max_singlecell_lr_score": float(lr_sub[score_col].max()) if not lr_sub.empty and score_col in lr_sub.columns else 0.0,
        })
    return pd.DataFrame(rows)


def score_candidate(row):
    score = 0
    reasons = []
    if "Evidence at protein level" in str(row.get("Evidence", "")):
        score += 1
        reasons.append("HPA protein evidence")
    text = " ".join(str(row.get(c, "")) for c in ["Protein class", "Secretome location", "Subcellular location"])
    if any(key.lower() in text.lower() for key in ["secreted", "membrane", "plasma proteins", "potential drug targets", "transporters"]):
        score += 1
        reasons.append("secreted/membrane/druggable class")
    if row.get("singlecell_lr_support", 0) > 0:
        score += 2
        reasons.append("CellChat cell_type LR axis detected")
    if "FDR=0.000" in str(row.get("bulk_meta_signal", "")) or "FDR=0.001" in str(row.get("bulk_meta_signal", "")):
        score += 1
        reasons.append("bulk meta survival signal")
    if str(row.get("best_platinum_cluster_signal", "")):
        score += 1
        reasons.append("platinum cluster association")
    if row["Gene"] in ["CXCR4", "CPT1A", "SDC4"]:
        score += 1
        reasons.append("cell-surface/metabolic intervention plausibility")
    priority = "High" if score >= 5 else "Medium" if score >= 3 else "Low"
    return pd.Series({"filter_score": score, "priority": priority, "filter_rationale": "; ".join(reasons)})


def main():
    hpa = load_hpa_candidates()
    loc = local_evidence()
    combined = loc.merge(hpa, on="Gene", how="left")
    scored = combined.join(combined.apply(score_candidate, axis=1))
    scored["PDC_CPTAC_filter"] = (
        "Use PDC000110 CPTAC ovarian proteome as orthogonal protein-level verification; "
        "machine-readable gene quantification was not pulled in this run."
    )
    scored["DepMap_filter"] = (
        "Use DepMap CRISPR/RNAi dependency and expression matrices to deprioritize pan-essential or non-expressed targets; "
        "portal download was blocked by browser verification in this environment."
    )
    scored["source_links"] = f"HPA={HPA_URL}; PDC_CPTAC={PDC_STUDY_URL}; DepMap={DEPMAP_DOWNLOAD_URL}"
    scored = scored.sort_values(["priority", "filter_score", "max_singlecell_lr_score"], ascending=[True, False, False])
    scored.to_csv(TABLE_DIR / "external_CPTAC_PDC_DepMap_HPA_candidate_filter.csv", index=False, encoding="utf-8-sig")

    lines = [
        "# External candidate filter",
        "",
        "HPA fields were downloaded from the official `proteinatlas.tsv.zip` table. PDC/CPTAC and DepMap were used as official filtering sources; DepMap matrix download was blocked by browser verification in this environment, so DepMap is recorded as a pending quantitative filter rather than a numeric pass/fail.",
        "",
        scored[[
            "Gene", "priority", "filter_score", "filter_rationale",
            "bulk_meta_signal", "singlecell_lr_support", "max_singlecell_lr_score",
            "Evidence", "Protein class",
            "Cancer prognostics - Ovary Serous Cystadenocarcinoma (TCGA)",
        ]].to_csv(index=False).replace(",", " | "),
        "",
        f"HPA download: {HPA_URL}",
        f"PDC/CPTAC ovarian study: {PDC_STUDY_URL}",
        f"DepMap download portal: {DEPMAP_DOWNLOAD_URL}",
    ]
    (REPORT_DIR / "external_candidate_filter.md").write_text("\n".join(lines), encoding="utf-8")
    print("Wrote external candidate filter")


if __name__ == "__main__":
    main()
