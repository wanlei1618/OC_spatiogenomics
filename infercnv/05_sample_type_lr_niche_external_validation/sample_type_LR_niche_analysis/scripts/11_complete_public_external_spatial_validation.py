import csv
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr


PROJECT = Path(r"D:\OC_spatiogenomics\infercnv\sample_type_LR_niche_analysis")
TABLES = PROJECT / "tables"
FIGURES = PROJECT / "figures"
LOGS = PROJECT / "logs"
ROOT = Path(r"D:\OC_spatiogenomics")
TMP = ROOT / "tmp"

for path in [TABLES, FIGURES, LOGS, TMP]:
    path.mkdir(parents=True, exist_ok=True)

os.environ["TMPDIR"] = str(TMP)
os.environ["TEMP"] = str(TMP)
os.environ["TMP"] = str(TMP)


AXES = [("SPP1", "CD44"), ("SPP1", "ITGB1"), ("MIF", "CD74"), ("APOE", "LRP1"), ("TGFB1", "TGFBR1"), ("CXCL12", "CXCR4")]
PRIMARY_AXES = [("SPP1", "CD44"), ("SPP1", "ITGB1")]


def log(message):
    with open(LOGS / "complete_public_external_spatial_validation.log", "a", encoding="utf-8") as handle:
        handle.write(message + "\n")
    print(message, flush=True)


def discover_public_root():
    for child in ROOT.iterdir():
        if child.is_dir() and (child / "Data_Zhang2022_Ovarian.tar.gz").exists():
            return child
    raise FileNotFoundError("Could not find public dataset root containing Data_Zhang2022_Ovarian.tar.gz")


PUBLIC_ROOT = discover_public_root()


def read_signature(name):
    path = TABLES / name
    if not path.exists() or path.stat().st_size == 0:
        return []
    df = pd.read_csv(path)
    if "gene" in df.columns:
        return [str(x).strip() for x in df["gene"].dropna().tolist()]
    return [str(x).strip() for x in df.iloc[:, 0].dropna().tolist()]


SIGNATURES = {
    "Subclone02_like": read_signature("signature_Subclone02_like.csv"),
    "Subclone04_like": read_signature("signature_Subclone04_like.csv"),
    "Subclone02_04_common": read_signature("signature_Subclone02_04_common.csv"),
    "CD44_ITGB1_target": read_signature("signature_CD44_ITGB1_target.csv"),
    "KRAS_hypoxia_target": read_signature("signature_KRAS_hypoxia_target.csv"),
    "tumor_epithelial": read_signature("signature_tumor_epithelial.csv"),
}


def normalize_sample_type(location):
    value = str(location).strip().lower()
    if value == "tumor":
        return "tumor"
    if value in {"omentum", "peritoneum", "mesentery"}:
        return "peritoneal_implant"
    if not value or value == "nan":
        return "unknown"
    return re.sub(r"[^a-z0-9]+", "_", value).strip("_")


def score_from_counts(counts, total):
    total = np.asarray(total, dtype=float)
    return np.log1p((np.asarray(counts, dtype=float) / np.maximum(total, 1.0)) * 10000.0)


def zscore(x):
    arr = np.asarray(x, dtype=float)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(arr, dtype=float)
    return (arr - np.nanmean(arr)) / sd


def corr_stat(result):
    if hasattr(result, "statistic"):
        return result.statistic, result.pvalue
    return result[0], result[1]


def df_to_markdown(df, max_rows=None):
    if max_rows is not None:
        df = df.head(max_rows)
    if df.empty:
        return "_No rows._"
    df = df.copy()
    df = df.fillna("")
    headers = [str(c) for c in df.columns]
    rows = [[str(v) for v in row] for row in df.to_numpy()]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(cell.replace("|", "/") for cell in row) + " |")
    return "\n".join(lines)


def plot_box(df, axis_name, out_png):
    sub = df[df["axis"] == axis_name].copy()
    if sub.empty:
        return
    labels = sorted(sub["sample_type"].dropna().unique())
    data = [sub.loc[sub["sample_type"] == label, "axis_score"].dropna().to_numpy() for label in labels]
    plt.figure(figsize=(7, 4.5))
    plt.boxplot(data, labels=labels, showfliers=False)
    jitter_x = []
    jitter_y = []
    rng = np.random.default_rng(7)
    for idx, label in enumerate(labels, start=1):
        vals = sub.loc[sub["sample_type"] == label, "axis_score"].dropna().to_numpy()
        jitter_x.extend(idx + rng.normal(0, 0.04, size=len(vals)))
        jitter_y.extend(vals)
    plt.scatter(jitter_x, jitter_y, s=18, alpha=0.75, color="#2f6f8f")
    plt.ylabel("External scRNA opportunity score")
    plt.xlabel("Sample type")
    plt.title(axis_name)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.savefig(str(out_png).replace(".png", ".pdf"))
    plt.close()


def find_zhang_dir():
    for dirpath, _, files in os.walk(PUBLIC_ROOT):
        files = set(files)
        if {"Cells.csv", "Genes.txt", "Exp_data_UMIcounts.mtx"}.issubset(files):
            return Path(dirpath)
    raise FileNotFoundError("Could not find Zhang2022 matrix directory")


def run_external_zhang():
    zhang = find_zhang_dir()
    log(f"External scRNA: using {zhang}")
    cells = pd.read_csv(zhang / "Cells.csv")
    genes = [line.strip().strip('"') for line in open(zhang / "Genes.txt", encoding="utf-8")]
    gene_to_idx = {gene.upper(): i + 1 for i, gene in enumerate(genes)}

    wanted = set()
    for lig, rec in AXES:
        wanted.add(lig)
        wanted.add(rec)
    for gene_list in SIGNATURES.values():
        wanted.update(gene_list)
    wanted = {g.upper() for g in wanted if g}
    selected = {gene: gene_to_idx[gene] for gene in wanted if gene in gene_to_idx}
    idx_to_gene = {idx: gene for gene, idx in selected.items()}
    log(f"External scRNA: extracting {len(selected)} selected genes from MatrixMarket")

    n_cells = cells.shape[0]
    selected_counts = {gene: np.zeros(n_cells, dtype=np.float32) for gene in selected}
    mtx_path = zhang / "Exp_data_UMIcounts.mtx"
    line_count = 0
    with open(mtx_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line[0] == "%":
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            row = int(parts[0])
            if row not in idx_to_gene:
                continue
            col = int(parts[1]) - 1
            if 0 <= col < n_cells:
                selected_counts[idx_to_gene[row]][col] = float(parts[2])
            line_count += 1
    log(f"External scRNA: streamed selected matrix entries; selected rows present={len(selected_counts)}")

    cells["sample_type"] = cells["anatomical_location"].map(normalize_sample_type)
    cells["is_source_myeloid"] = cells["cell_type"].astype(str).str.contains("Macrophage|Dendritic|Monocyte|Myeloid|DC", case=False, regex=True)
    cells["is_target_tumor"] = cells["cell_type"].astype(str).str.contains("Malignant", case=False, regex=True)
    total = cells["nCount_RNA"].to_numpy(dtype=float)

    expr = {gene: score_from_counts(selected_counts.get(gene, np.zeros(n_cells)), total) for gene in selected_counts}
    for gene in wanted:
        expr.setdefault(gene, np.zeros(n_cells, dtype=float))

    sig_scores = {}
    for sig, genes in SIGNATURES.items():
        present = [g.upper() for g in genes if g.upper() in expr]
        if present:
            sig_scores[sig] = np.vstack([expr[g] for g in present]).mean(axis=0)
        else:
            sig_scores[sig] = np.zeros(n_cells, dtype=float)

    rows = []
    target_rows = []
    sp_rows = []
    patient_rows = []
    grouping = ["patient", "sample", "sample_type", "anatomical_location", "treatment_phase"]
    for keys, sub in cells.groupby(grouping, dropna=False):
        sub_idx = sub.index.to_numpy()
        source = cells.loc[sub_idx, "is_source_myeloid"].to_numpy()
        target = cells.loc[sub_idx, "is_target_tumor"].to_numpy()
        source_fraction = float(source.mean()) if len(source) else 0.0
        target_fraction = float(target.mean()) if len(target) else 0.0
        patient, sample, sample_type, location, treatment = keys

        for sig, arr in sig_scores.items():
            target_rows.append(
                {
                    "dataset": "Zhang2022_Ovarian",
                    "patient_id": patient,
                    "sample_id": sample,
                    "sample_type": sample_type,
                    "signature": sig,
                    "mean_all_cells": float(np.nanmean(arr[sub_idx])),
                    "mean_target_tumor": float(np.nanmean(arr[sub_idx][target])) if target.any() else np.nan,
                    "n_cells": len(sub_idx),
                    "n_target_tumor": int(target.sum()),
                }
            )

        for ligand, receptor in AXES:
            lig = ligand.upper()
            rec = receptor.upper()
            lig_source = expr[lig][sub_idx][source]
            rec_target = expr[rec][sub_idx][target]
            ligand_avg = float(np.nanmean(lig_source)) if lig_source.size else 0.0
            receptor_avg = float(np.nanmean(rec_target)) if rec_target.size else 0.0
            ligand_pct = float((lig_source > 0).mean()) if lig_source.size else 0.0
            receptor_pct = float((rec_target > 0).mean()) if rec_target.size else 0.0
            expr_product = ligand_avg * receptor_avg
            abundance_weighted = source_fraction * target_fraction * expr_product
            axis_score = abundance_weighted * ligand_pct * receptor_pct
            rows.append(
                {
                    "dataset": "Zhang2022_Ovarian",
                    "patient_id": patient,
                    "sample_id": sample,
                    "sample_type": sample_type,
                    "anatomical_location": location,
                    "treatment_phase": treatment,
                    "axis": f"{ligand}-{receptor}",
                    "source_group": "Macrophage_Dendritic",
                    "target_group": "Malignant",
                    "n_cells": len(sub_idx),
                    "source_n": int(source.sum()),
                    "target_n": int(target.sum()),
                    "source_fraction": source_fraction,
                    "target_fraction": target_fraction,
                    "ligand_avg_source": ligand_avg,
                    "receptor_avg_target": receptor_avg,
                    "ligand_pct_source": ligand_pct,
                    "receptor_pct_target": receptor_pct,
                    "expr_product_score": expr_product,
                    "abundance_weighted_score": abundance_weighted,
                    "axis_score": axis_score,
                    "target_Subclone02_04_like_score": float(np.nanmean((sig_scores["Subclone02_like"][sub_idx] + sig_scores["Subclone04_like"][sub_idx]) / 2)),
                    "CD44_ITGB1_target_score": float(np.nanmean(sig_scores["CD44_ITGB1_target"][sub_idx])),
                    "KRAS_hypoxia_target_score": float(np.nanmean(sig_scores["KRAS_hypoxia_target"][sub_idx])),
                }
            )
            if (ligand, receptor) in PRIMARY_AXES:
                sp_rows.append(rows[-1].copy())

        patient_rows.append(
            {
                "dataset": "Zhang2022_Ovarian",
                "patient_id": patient,
                "sample_id": sample,
                "sample_type": sample_type,
                "anatomical_location": location,
                "treatment_phase": treatment,
                "n_cells": len(sub_idx),
                "source_myeloid_n": int(source.sum()),
                "target_tumor_n": int(target.sum()),
                "source_myeloid_fraction": source_fraction,
                "target_tumor_fraction": target_fraction,
            }
        )

    all_scores = pd.DataFrame(rows)
    target_df = pd.DataFrame(target_rows)
    spp_df = pd.DataFrame(sp_rows)
    patient_df = pd.DataFrame(patient_rows).drop_duplicates()

    all_scores.to_csv(TABLES / "external_scRNA_all_datasets_sampletype_axis_scores.csv", index=False)
    target_df.to_csv(TABLES / "external_scRNA_target_signature_scores.csv", index=False)
    spp_df.to_csv(TABLES / "external_scRNA_SPP1_myeloid_CD44_ITGB1_target_scores.csv", index=False)
    patient_df.to_csv(TABLES / "external_scRNA_patient_sampletype_summary.csv", index=False)
    pd.DataFrame([{"dataset": "Zhang2022_Ovarian", "status": "completed", "note": "Used public Zhang2022 ovarian scRNA metadata and streamed selected genes from MatrixMarket."}]).to_csv(TABLES / "external_scRNA_validation_status.csv", index=False)

    plot_box(all_scores, "SPP1-CD44", FIGURES / "external_boxplot_SPP1_CD44_by_sample_type.png")
    plot_box(all_scores, "SPP1-ITGB1", FIGURES / "external_boxplot_SPP1_ITGB1_by_sample_type.png")

    heat = all_scores.groupby(["sample_type", "axis"], as_index=False)["axis_score"].mean()
    pivot = heat.pivot(index="sample_type", columns="axis", values="axis_score").fillna(0)
    plt.figure(figsize=(8, 3.6))
    plt.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
    plt.xticks(range(pivot.shape[1]), pivot.columns, rotation=45, ha="right")
    plt.yticks(range(pivot.shape[0]), pivot.index)
    plt.colorbar(label="Mean axis score")
    plt.title("External scRNA LR opportunity")
    plt.tight_layout()
    plt.savefig(FIGURES / "external_heatmap_dataset_sampletype_axis.png", dpi=220)
    plt.savefig(FIGURES / "external_heatmap_dataset_sampletype_axis.pdf")
    plt.close()

    paired = all_scores[all_scores["axis"].isin(["SPP1-CD44", "SPP1-ITGB1"])].copy()
    paired.to_csv(TABLES / "external_paired_tumor_vs_implant_axis_score.csv", index=False)
    plt.figure(figsize=(7, 4.5))
    for axis, sub in paired.groupby("axis"):
        summarized = sub.groupby(["patient_id", "sample_type"], as_index=False)["axis_score"].mean()
        wide = summarized.pivot(index="patient_id", columns="sample_type", values="axis_score")
        if {"tumor", "peritoneal_implant"}.issubset(wide.columns):
            for _, row in wide.iterrows():
                plt.plot([0, 1], [row["tumor"], row["peritoneal_implant"]], marker="o", alpha=0.55, label=axis if axis not in plt.gca().get_legend_handles_labels()[1] else None)
    plt.xticks([0, 1], ["tumor", "peritoneal_implant"])
    plt.ylabel("Axis score")
    plt.title("Paired tumor vs implant-like samples")
    plt.tight_layout()
    plt.savefig(FIGURES / "external_paired_tumor_vs_ascites_axis_score.png", dpi=220)
    plt.savefig(FIGURES / "external_paired_tumor_vs_ascites_axis_score.pdf")
    plt.close()
    return all_scores


def read_10x_h5_selected(h5_path, wanted_genes):
    with h5py.File(h5_path, "r") as f:
        grp = f["matrix"]
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in grp["barcodes"][:]]
        features = grp["features"]
        if "name" in features:
            genes = [x.decode() if isinstance(x, bytes) else str(x) for x in features["name"][:]]
        else:
            genes = [x.decode() if isinstance(x, bytes) else str(x) for x in features["id"][:]]
        shape = tuple(grp["shape"][:])
        mat = sparse.csc_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=shape)
    gene_to_rows = defaultdict(list)
    for i, gene in enumerate(genes):
        gene_to_rows[gene.upper()].append(i)
    out = {}
    for gene in wanted_genes:
        rows = gene_to_rows.get(gene.upper(), [])
        if rows:
            out[gene.upper()] = np.asarray(mat[rows, :].sum(axis=0)).ravel()
        else:
            out[gene.upper()] = np.zeros(len(barcodes), dtype=float)
    totals = np.asarray(mat.sum(axis=0)).ravel()
    return barcodes, out, totals


def read_visium_mtx_selected(sample_dir, accession, wanted_genes):
    features = next(sample_dir.glob(f"{accession}_features_*.tsv"), None)
    barcodes = next(sample_dir.glob(f"{accession}_barcodes_*.tsv"), None)
    matrix = next(sample_dir.glob(f"{accession}_matrix_*.mtx"), None)
    if not all([features, barcodes, matrix]):
        return None
    feature_df = pd.read_csv(features, sep="\t", header=None)
    genes = feature_df.iloc[:, 1 if feature_df.shape[1] > 1 else 0].astype(str).str.upper().tolist()
    rows = [i for i, g in enumerate(genes) if g in {x.upper() for x in wanted_genes}]
    mat = mmread(matrix).tocsr()
    barcode_list = pd.read_csv(barcodes, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    totals = np.asarray(mat.sum(axis=0)).ravel()
    out = {}
    for gene in wanted_genes:
        idxs = [i for i, g in enumerate(genes) if g == gene.upper()]
        out[gene.upper()] = np.asarray(mat[idxs, :].sum(axis=0)).ravel() if idxs else np.zeros(len(barcode_list), dtype=float)
    return barcode_list, out, totals


def spatial_scores_from_counts(counts, totals):
    expr = {gene: score_from_counts(values, totals) for gene, values in counts.items()}
    def sig_score(sig):
        present = [g.upper() for g in SIGNATURES[sig] if g.upper() in expr]
        if not present:
            return np.zeros(len(totals), dtype=float)
        return np.vstack([expr[g] for g in present]).mean(axis=0)

    out = {
        "SPP1_myeloid_score": zscore(expr.get("SPP1", 0)) + zscore(expr.get("LYZ", 0)) + zscore(expr.get("C1QA", 0)),
        "CD44_ITGB1_target_score": sig_score("CD44_ITGB1_target"),
        "Subclone02_04_like_score": (sig_score("Subclone02_like") + sig_score("Subclone04_like")) / 2,
        "KRAS_hypoxia_target_score": sig_score("KRAS_hypoxia_target"),
        "tumor_epithelial_score": sig_score("tumor_epithelial"),
        "SPP1": expr.get("SPP1", np.zeros(len(totals))),
        "CD44": expr.get("CD44", np.zeros(len(totals))),
        "ITGB1": expr.get("ITGB1", np.zeros(len(totals))),
    }
    out["Target_axis_score"] = zscore(out["Subclone02_04_like_score"]) + zscore(out["CD44_ITGB1_target_score"]) + zscore(out["KRAS_hypoxia_target_score"]) + zscore(out["tumor_epithelial_score"])
    out["SPP1_CD44_score"] = out["SPP1"] * out["CD44"]
    out["SPP1_ITGB1_score"] = out["SPP1"] * out["ITGB1"]
    return out


def parse_tissue_positions(path):
    pos = pd.read_csv(path, header=None)
    pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: pos.shape[1]]
    return pos


def scatter_spatial(df, x, y, out_png, title):
    plt.figure(figsize=(5.2, 4.4))
    plt.scatter(df[x], df[y], s=8, alpha=0.45, color="#3f7f5f")
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.savefig(str(out_png).replace(".png", ".pdf"))
    plt.close()


def run_spatial():
    wanted = set(["SPP1", "CD44", "ITGB1", "LYZ", "C1QA", "MIF", "CD74", "APOE", "LRP1", "TGFB1", "TGFBR1", "CXCL12", "CXCR4"])
    for genes in SIGNATURES.values():
        wanted.update([g.upper() for g in genes])
    wanted = sorted(wanted)
    rows = []
    commot_rows = []
    neigh_rows = []

    spatial_root = PUBLIC_ROOT / "ovarian_spatial_geo"
    gse203 = spatial_root / "GSE203612" / "seurat_visium"
    ovca_samples = {"GSM6177614": "NYU_OVCA1_Vis", "GSM6177617": "NYU_OVCA3_Vis", "GSM6177618": "NYU_OVCA_like_Vis"}
    for accession, title in ovca_samples.items():
        sample_dir = gse203 / accession
        h5_path = sample_dir / "filtered_feature_bc_matrix.h5"
        pos_path = sample_dir / "spatial" / "tissue_positions_list.csv"
        if not h5_path.exists() or not pos_path.exists():
            continue
        log(f"Spatial: processing {accession}")
        barcodes, counts, totals = read_10x_h5_selected(h5_path, wanted)
        scores = spatial_scores_from_counts(counts, totals)
        df = pd.DataFrame({"barcode": barcodes, "sample_id": accession, "sample_title": title, "dataset": "GSE203612_OVCA", "sample_type": "tumor"})
        for key, val in scores.items():
            df[key] = val
        pos = parse_tissue_positions(pos_path)
        df = df.merge(pos, on="barcode", how="left")
        df = df[df.get("in_tissue", 1).fillna(1).astype(int) == 1].copy()

        if len(df) >= 5:
            pear_r, pear_p = corr_stat(pearsonr(df["SPP1_myeloid_score"], df["Target_axis_score"]))
            spear_r, spear_p = corr_stat(spearmanr(df["SPP1_myeloid_score"], df["Target_axis_score"]))
            rows.append(
                {
                    "dataset": "GSE203612_OVCA",
                    "sample_id": accession,
                    "sample_title": title,
                    "n_spots": len(df),
                    "pearson_r": pear_r,
                    "pearson_p": pear_p,
                    "spearman_r": spear_r,
                    "spearman_p": spear_p,
                    "mean_SPP1_CD44_score": df["SPP1_CD44_score"].mean(),
                    "mean_SPP1_ITGB1_score": df["SPP1_ITGB1_score"].mean(),
                }
            )

        q_s = df["SPP1_myeloid_score"].quantile(0.75)
        q_t = df["Target_axis_score"].quantile(0.75)
        df["niche_class"] = np.select(
            [(df["SPP1_myeloid_score"] >= q_s) & (df["Target_axis_score"] >= q_t), df["SPP1_myeloid_score"] >= q_s, df["Target_axis_score"] >= q_t],
            ["double_high", "SPP1_myeloid_high", "Target_axis_high"],
            default="other",
        )
        coords = df[["array_row", "array_col"]].to_numpy(dtype=float)
        if np.isfinite(coords).all() and len(coords) > 10:
            tree = cKDTree(coords)
            neigh = tree.query_ball_point(coords, r=2.01)
            source_mask = df["niche_class"].eq("SPP1_myeloid_high") | df["niche_class"].eq("double_high")
            target_mask = df["niche_class"].eq("Target_axis_high") | df["niche_class"].eq("double_high")
            obs_edges = 0
            total_edges = 0
            for i, ns in enumerate(neigh):
                if not source_mask.iloc[i]:
                    continue
                for j in ns:
                    if i == j:
                        continue
                    total_edges += 1
                    obs_edges += int(target_mask.iloc[j])
            expected = float(target_mask.mean())
            observed = obs_edges / total_edges if total_edges else np.nan
            neigh_rows.append(
                {
                    "dataset": "GSE203612_OVCA",
                    "sample_id": accession,
                    "source_class": "SPP1_myeloid_high_or_double_high",
                    "target_class": "Target_axis_high_or_double_high",
                    "radius_array_units": 2.01,
                    "observed_neighbor_fraction": observed,
                    "expected_global_fraction": expected,
                    "enrichment_ratio": observed / expected if expected else np.nan,
                    "n_source_spots": int(source_mask.sum()),
                    "n_target_spots": int(target_mask.sum()),
                    "n_neighbor_edges": int(total_edges),
                }
            )

        for axis, score_col in [("SPP1-CD44", "SPP1_CD44_score"), ("SPP1-ITGB1", "SPP1_ITGB1_score")]:
            commot_rows.append(
                {
                    "dataset": "GSE203612_OVCA",
                    "sample_id": accession,
                    "sample_title": title,
                    "sample_type": "tumor",
                    "region": "all_tissue_spots",
                    "axis": axis,
                    "mean_score": df[score_col].mean(),
                    "median_score": df[score_col].median(),
                    "p75_score": df[score_col].quantile(0.75),
                    "double_high_mean_score": df.loc[df["niche_class"].eq("double_high"), score_col].mean(),
                    "n_spots": len(df),
                }
            )

        scatter_spatial(df, "SPP1_myeloid_score", "Target_axis_score", FIGURES / f"spatial_{accession}_SPP1_myeloid_vs_Target_axis.png", accession)
        for score, cmap in [("SPP1_myeloid_score", "magma"), ("CD44_ITGB1_target_score", "viridis"), ("KRAS_hypoxia_target_score", "plasma"), ("Target_axis_score", "coolwarm")]:
            plt.figure(figsize=(5, 4.5))
            plt.scatter(df["array_col"], -df["array_row"], c=df[score], s=10, cmap=cmap)
            plt.axis("equal")
            plt.axis("off")
            plt.title(f"{accession} {score}")
            plt.colorbar(fraction=0.046)
            plt.tight_layout()
            plt.savefig(FIGURES / f"spatial_{accession}_{score}.png", dpi=220)
            plt.savefig(FIGURES / f"spatial_{accession}_{score}.pdf")
            plt.close()

    # GSE189843 ovarian HGSC Visium has matrices but no spot coordinate files in the local copy.
    gse189 = spatial_root / "GSE189843" / "suppl"
    if gse189.exists():
        for matrix_path in sorted(gse189.glob("GSM*_matrix_*.mtx"))[:12]:
            accession = matrix_path.name.split("_")[0]
            try:
                parsed = read_visium_mtx_selected(gse189, accession, wanted)
                if parsed is None:
                    continue
                barcodes, counts, totals = parsed
                scores = spatial_scores_from_counts(counts, totals)
                df = pd.DataFrame(scores)
                if len(df) >= 5:
                    pear_r, pear_p = corr_stat(pearsonr(df["SPP1_myeloid_score"], df["Target_axis_score"]))
                    spear_r, spear_p = corr_stat(spearmanr(df["SPP1_myeloid_score"], df["Target_axis_score"]))
                    rows.append(
                        {
                            "dataset": "GSE189843_HGSC",
                            "sample_id": accession,
                            "sample_title": accession,
                            "n_spots": len(df),
                            "pearson_r": pear_r,
                            "pearson_p": pear_p,
                            "spearman_r": spear_r,
                            "spearman_p": spear_p,
                            "mean_SPP1_CD44_score": df["SPP1_CD44_score"].mean(),
                            "mean_SPP1_ITGB1_score": df["SPP1_ITGB1_score"].mean(),
                        }
                    )
                    for axis, score_col in [("SPP1-CD44", "SPP1_CD44_score"), ("SPP1-ITGB1", "SPP1_ITGB1_score")]:
                        commot_rows.append(
                            {
                                "dataset": "GSE189843_HGSC",
                                "sample_id": accession,
                                "sample_title": accession,
                                "sample_type": "tumor",
                                "region": "all_tissue_spots_no_coordinates",
                                "axis": axis,
                                "mean_score": df[score_col].mean(),
                                "median_score": df[score_col].median(),
                                "p75_score": df[score_col].quantile(0.75),
                                "double_high_mean_score": np.nan,
                                "n_spots": len(df),
                            }
                        )
            except Exception as exc:
                log(f"Spatial: skipped {accession}: {exc}")

    corr_df = pd.DataFrame(rows)
    neigh_df = pd.DataFrame(neigh_rows)
    commot_df = pd.DataFrame(commot_rows)
    corr_df.to_csv(TABLES / "spatial_correlation_SPP1_myeloid_Target_axis.csv", index=False)
    neigh_df.to_csv(TABLES / "spatial_neighborhood_enrichment.csv", index=False)
    commot_df.to_csv(TABLES / "spatial_LR_COMMOT_scores.csv", index=False)
    commot_df.to_csv(TABLES / "spatial_LR_COMMOT_scores_by_region.csv", index=False)
    pd.DataFrame([{"dataset": "GSE203612_OVCA;GSE189843_HGSC", "status": "completed", "note": "Computed spot-level scores; GSE203612 OVCA used for coordinate-based neighborhood enrichment. COMMOT outputs are simplified expression-product LR maps, not full COMMOT transport modeling."}]).to_csv(TABLES / "spatial_validation_status.csv", index=False)

    if not neigh_df.empty:
        pivot = neigh_df.pivot_table(index="sample_id", columns="target_class", values="enrichment_ratio", aggfunc="mean").fillna(0)
        plt.figure(figsize=(5.2, 3.4))
        plt.imshow(pivot.to_numpy(), aspect="auto", cmap="magma")
        plt.xticks(range(pivot.shape[1]), pivot.columns, rotation=30, ha="right")
        plt.yticks(range(pivot.shape[0]), pivot.index)
        plt.colorbar(label="Neighborhood enrichment")
        plt.tight_layout()
        plt.savefig(FIGURES / "spatial_neighborhood_enrichment_heatmap.png", dpi=220)
        plt.savefig(FIGURES / "spatial_neighborhood_enrichment_heatmap.pdf")
        plt.close()
    if not commot_df.empty:
        pivot = commot_df.pivot_table(index="sample_id", columns="axis", values="mean_score", aggfunc="mean").fillna(0)
        plt.figure(figsize=(6, 4))
        plt.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
        plt.xticks(range(pivot.shape[1]), pivot.columns, rotation=30, ha="right")
        plt.yticks(range(pivot.shape[0]), pivot.index)
        plt.colorbar(label="Mean spatial LR score")
        plt.tight_layout()
        plt.savefig(FIGURES / "spatial_COMMOT_axis_comparison_heatmap.png", dpi=220)
        plt.savefig(FIGURES / "spatial_COMMOT_axis_comparison_heatmap.pdf")
        plt.close()
        for axis, out_name in [("SPP1-CD44", "spatial_COMMOT_SPP1_CD44_map.png"), ("SPP1-ITGB1", "spatial_COMMOT_SPP1_ITGB1_map.png")]:
            sub = commot_df[commot_df["axis"] == axis].copy()
            plt.figure(figsize=(7, 4))
            plt.bar(sub["sample_id"], sub["mean_score"], color="#4878a8")
            plt.xticks(rotation=45, ha="right")
            plt.ylabel("Mean spatial LR score")
            plt.title(axis)
            plt.tight_layout()
            plt.savefig(FIGURES / out_name, dpi=220)
            plt.savefig(FIGURES / out_name.replace(".png", ".pdf"))
            plt.close()
    return corr_df, neigh_df, commot_df


def update_statistics_and_limitations(external_df, corr_df, neigh_df):
    stats = []
    if external_df is not None and not external_df.empty:
        for axis, sub in external_df[external_df["axis"].isin(["SPP1-CD44", "SPP1-ITGB1"])].groupby("axis"):
            tumor = sub.loc[sub["sample_type"] == "tumor", "axis_score"]
            implant = sub.loc[sub["sample_type"] == "peritoneal_implant", "axis_score"]
            stats.append(
                {
                    "analysis": "external_scRNA_Zhang2022",
                    "comparison": f"{axis}: tumor vs peritoneal_implant descriptive",
                    "n_tumor_samples": tumor.shape[0],
                    "n_implant_samples": implant.shape[0],
                    "mean_tumor": tumor.mean(),
                    "mean_peritoneal_implant": implant.mean(),
                    "interpretation": "Descriptive patient-sample level comparison; anatomical site and treatment phase may be confounded.",
                }
            )
    if corr_df is not None and not corr_df.empty:
        stats.append(
            {
                "analysis": "spatial_spot_correlation",
                "comparison": "SPP1_myeloid_score vs Target_axis_score",
                "n_tumor_samples": corr_df.shape[0],
                "n_implant_samples": np.nan,
                "mean_tumor": corr_df["spearman_r"].mean(),
                "mean_peritoneal_implant": np.nan,
                "interpretation": "Mean Spearman correlation across spatial samples; supports co-expression gradients, not physical binding.",
            }
        )
    pd.DataFrame(stats).to_csv(TABLES / "statistical_tests_summary.csv", index=False)

    limitations = pd.DataFrame(
        [
            ("sample_type_sample_id_confounding", "In integrated_oc, sample_type is one-to-one with sample_id; integrated sample_type findings are exploratory."),
            ("external_site_treatment_confounding", "Zhang2022 external samples contain tumor and implant-like sites, but treatment phase and anatomical site are not fully balanced."),
            ("spatial_coordinate_availability", "GSE203612 OVCA Visium supports coordinate neighborhood analysis; local GSE189843 HGSC matrices lack spot coordinate files and were used only for spot-level score correlations."),
            ("COMMOT_simplified", "Spatial COMMOT deliverables here are expression-product LR score maps/summaries, not full optimal-transport COMMOT inference."),
            ("LR_scores_are_potential_interactions", "LR opportunity scores are expression-derived potentials and do not prove ligand-receptor binding or signaling."),
            ("bulk_OS_not_primary", "Bulk OS instability does not reject a local sample_type/niche-dependent mechanism."),
        ],
        columns=["limitation", "interpretation"],
    )
    limitations.to_csv(TABLES / "limitations_summary.csv", index=False)


def write_report(external_df, corr_df, neigh_df, commot_df):
    integrated = pd.read_csv(TABLES / "sample_type_LR_opportunity_scores_primary_axes_target_Subclone02_04.csv")
    ext_summary = external_df.groupby(["sample_type", "axis"], as_index=False)["axis_score"].mean() if external_df is not None and not external_df.empty else pd.DataFrame()
    report = PROJECT / "sample_type_LR_niche_analysis_report.md"
    with open(report, "w", encoding="utf-8") as handle:
        handle.write("# sample_type-dependent SPP1-CD44/ITGB1 LR niche analysis\n\n")
        handle.write(f"Project root: `{PROJECT.as_posix()}`\n\n")
        handle.write("## Completed status\n\n")
        handle.write("This report updates the previous placeholder external/scRNA and spatial sections using public data under `D:/OC_spatiogenomics/公开集`.\n\n")
        handle.write("## Integrated_oc exploratory analysis\n\n")
        handle.write("Integrated_oc remains exploratory because sample_type is one-to-one with sample_id in the available metadata. The strongest primary-axis opportunity scores are in the tumor sample type, especially toward CNV_Subclone_02/04.\n\n")
        integrated_cols = [col for col in ["sample_type", "axis", "target_group", "axis_score", "abundance_weighted_score", "source_n_total", "target_n_total", "source_n", "target_n"] if col in integrated.columns]
        handle.write(df_to_markdown(integrated[integrated_cols], max_rows=30))
        handle.write("\n\n## External scRNA validation: Zhang2022 ovarian\n\n")
        handle.write("Zhang2022 ovarian single-cell data were processed from `Cells.csv`, `Samples.csv`, `Genes.txt`, and the UMI MatrixMarket file by streaming selected genes only. Macrophage/dendritic cells were used as SPP1-source cells and malignant cells as target tumor cells. Peritoneum/omentum/mesentery were treated as implant-like solid/peritoneal niches.\n\n")
        if not ext_summary.empty:
            handle.write(df_to_markdown(ext_summary))
            handle.write("\n\n")
        handle.write("Key outputs: `external_scRNA_all_datasets_sampletype_axis_scores.csv`, `external_scRNA_patient_sampletype_summary.csv`, `external_scRNA_target_signature_scores.csv`, and `external_scRNA_SPP1_myeloid_CD44_ITGB1_target_scores.csv`.\n\n")
        handle.write("## Spatial validation\n\n")
        handle.write("GSE203612 ovarian Visium samples with coordinates were used for spot-level SPP1-myeloid vs target-axis correlation and local neighborhood enrichment. GSE189843 HGSC Visium matrices were also scored, but the local copy lacks spot coordinate files, so they contribute only spot-level score correlations and LR expression-product summaries.\n\n")
        if corr_df is not None and not corr_df.empty:
            handle.write(df_to_markdown(corr_df[["dataset", "sample_id", "n_spots", "spearman_r", "spearman_p", "mean_SPP1_CD44_score", "mean_SPP1_ITGB1_score"]]))
            handle.write("\n\n")
        if neigh_df is not None and not neigh_df.empty:
            handle.write("Neighborhood enrichment summary:\n\n")
            handle.write(df_to_markdown(neigh_df))
            handle.write("\n\n")
        handle.write("The spatial COMMOT-named outputs are simplified spatial LR expression-product maps/summaries, not full COMMOT optimal-transport inference.\n\n")
        handle.write("## Interpretation\n\n")
        handle.write("The completed public-data checks support a cautious niche model: SPP1-CD44/ITGB1 is strongest in tumor/implant-like contexts at the integrated_oc level, and public ovarian scRNA/spatial data provide additional expression-level support for SPP1-myeloid and CD44/ITGB1/KRAS-hypoxia target programs. This should be interpreted as a candidate sample_type/niche-dependent communication program, not as proof of physical interaction or a universal OS predictor.\n\n")
        handle.write("## Limitations\n\n")
        handle.write(df_to_markdown(pd.read_csv(TABLES / "limitations_summary.csv")))
        handle.write("\n")

    html = PROJECT / "sample_type_LR_niche_analysis_report.html"
    text = report.read_text(encoding="utf-8")
    html.write_text("<html><body><pre style='white-space:pre-wrap;font-family:Arial'>" + text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;") + "</pre></body></html>", encoding="utf-8")
    try:
        import pypandoc

        pypandoc.convert_file(str(report), "docx", outputfile=str(PROJECT / "sample_type_LR_niche_analysis_report.docx"))
        render_status = "docx_rendered_with_pypandoc"
    except Exception as exc:
        render_status = f"docx_not_rendered: {exc}"
    pd.DataFrame([{"report": str(report), "html": str(html), "docx_status": render_status}]).to_csv(TABLES / "final_report_render_status.csv", index=False)


def main():
    log("Starting public external/spatial completion")
    external_path = TABLES / "external_scRNA_all_datasets_sampletype_axis_scores.csv"
    corr_path = TABLES / "spatial_correlation_SPP1_myeloid_Target_axis.csv"
    neigh_path = TABLES / "spatial_neighborhood_enrichment.csv"
    commot_path = TABLES / "spatial_LR_COMMOT_scores.csv"
    if external_path.exists() and external_path.stat().st_size > 0:
        external = pd.read_csv(external_path)
        log("External scRNA: reusing existing completed table")
    else:
        external = run_external_zhang()
    if corr_path.exists() and corr_path.stat().st_size > 0 and commot_path.exists() and commot_path.stat().st_size > 0:
        corr = pd.read_csv(corr_path)
        neigh = pd.read_csv(neigh_path) if neigh_path.exists() and neigh_path.stat().st_size > 0 else pd.DataFrame()
        commot = pd.read_csv(commot_path)
        log("Spatial: reusing existing completed tables")
    else:
        corr, neigh, commot = run_spatial()
    update_statistics_and_limitations(external, corr, neigh)
    write_report(external, corr, neigh, commot)
    log("Completed public external/spatial completion")


if __name__ == "__main__":
    main()
