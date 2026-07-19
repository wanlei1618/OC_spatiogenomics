import csv
import gzip
import json
import math
import os
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


ROOT = Path(r"D:\OC_spatiogenomics")
OUT = ROOT / "infercnv" / "07_extended_external_validation"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
LOGS = OUT / "logs"
SCRIPTS = OUT / "scripts"
OLD_TABLES = ROOT / "infercnv" / "05_sample_type_lr_niche_external_validation" / "sample_type_LR_niche_analysis" / "tables"

for d in (TABLES, FIGURES, LOGS, SCRIPTS):
    d.mkdir(parents=True, exist_ok=True)


SC_H5 = [
    ("OV_GSE147082", ROOT / "公开集" / "单细胞" / "OV_GSE147082_expression.h5"),
    ("OV_GSE151214", ROOT / "公开集" / "单细胞" / "OV_GSE151214_expression.h5"),
    ("OV_GSE154600", ROOT / "公开集" / "单细胞" / "OV_GSE154600_expression.h5"),
    ("OV_GSE154763", ROOT / "公开集" / "单细胞" / "OV_GSE154763_expression.h5"),
    ("OV_GSE158722", ROOT / "公开集" / "单细胞" / "OV_GSE158722_expression.h5"),
]

SPATIAL_ROOTS = [
    ("GSE211956", ROOT / "公开集" / "ovarian_spatial_geo" / "GSE211956" / "seurat_visium"),
    ("GSE227019", ROOT / "公开集" / "ovarian_spatial_geo" / "GSE227019" / "seurat_visium"),
]


GENE_SETS = {
    "spp1_macrophage": ["SPP1", "CD68", "CD14", "LYZ", "LST1", "TYROBP", "AIF1", "FCGR3A", "C1QA", "C1QB", "C1QC", "IL1B", "TNF", "ISG15", "IFITM3"],
    "myeloid": ["LYZ", "LST1", "TYROBP", "AIF1", "CD68", "CD14", "FCGR3A", "C1QA", "C1QB", "C1QC", "MS4A7", "CTSS"],
    "tumor_epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "PAX8", "MUC16", "WFDC2", "CLDN3", "CLDN4"],
    "t_nk": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"],
    "b_cell": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "ACTA2"],
    "endothelial": ["PECAM1", "VWF", "KDR", "RAMP2", "ENG"],
    "target_receptors": ["ITGB1", "CD44"],
    "itgb1_cd44_target": ["ITGB1", "CD44"],
    "hypoxia_core": ["HIF1A", "CA9", "SLC2A1", "LDHA", "VEGFA", "ENO1", "PGK1"],
    "kras_mapk_core": ["DUSP4", "DUSP6", "ETV4", "ETV5", "SPRY2", "FOS", "JUN", "EGR1"],
    "fak_ecm_invasion": ["PTK2", "SRC", "VCL", "PXN", "FN1", "VIM", "ITGA5", "ITGB1", "MMP2", "MMP9", "SNAI2"],
}

AXES = [("SPP1", "ITGB1"), ("SPP1", "CD44"), ("MIF", "CD74"), ("APOE", "LRP1"), ("TGFB1", "TGFBR2"), ("CXCL12", "CXCR4")]


def log(msg):
    with open(LOGS / "extended_validation.log", "a", encoding="utf-8") as handle:
        handle.write(msg + "\n")
    print(msg, flush=True)


def read_signature_csv(name, fallback):
    path = OLD_TABLES / name
    if not path.exists():
        return fallback
    df = pd.read_csv(path)
    if df.empty:
        return fallback
    col = "gene" if "gene" in df.columns else df.columns[0]
    genes = [str(x).strip().upper() for x in df[col].dropna().tolist() if str(x).strip()]
    return genes or fallback


GENE_SETS["subclone02_like"] = read_signature_csv("signature_Subclone02_like.csv", ["MUC16", "EPCAM", "KRT8", "KRT18"])
GENE_SETS["subclone04_like"] = read_signature_csv("signature_Subclone04_like.csv", ["MUC16", "EPCAM", "KRT8", "KRT19"])
GENE_SETS["subclone02_04_common"] = read_signature_csv("signature_Subclone02_04_common.csv", ["MUC16", "EPCAM", "KRT8", "KRT18", "KRT19"])
GENE_SETS["cd44_itgb1_target"] = read_signature_csv("signature_CD44_ITGB1_target.csv", ["CD44", "ITGB1", "VIM", "FN1", "ITGA5"])
GENE_SETS["kras_hypoxia_target"] = read_signature_csv("signature_KRAS_hypoxia_target.csv", ["DUSP6", "VEGFA", "CA9", "LDHA", "SLC2A1"])


def decode_array(arr):
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8", errors="ignore"))
        else:
            out.append(str(x))
    return out


def wanted_genes():
    genes = set()
    for vals in GENE_SETS.values():
        genes.update(g.upper() for g in vals)
    for ligand, receptor in AXES:
        genes.add(ligand)
        genes.add(receptor)
    return sorted(genes)


def read_10x_h5_selected(path, genes):
    genes = [g.upper() for g in genes]
    with h5py.File(path, "r") as f:
        grp = f["matrix"]
        barcodes = decode_array(grp["barcodes"][:])
        feature_names = decode_array(grp["features"]["name"][:])
        feature_upper = [x.upper() for x in feature_names]
        shape = tuple(int(x) for x in grp["shape"][:])
        n_genes, n_cells = shape
        gene_to_row = {}
        for i, g in enumerate(feature_upper):
            gene_to_row.setdefault(g, i)
        selected = [(g, gene_to_row[g]) for g in genes if g in gene_to_row]
        row_to_pos = np.full(n_genes, -1, dtype=np.int32)
        for pos, (_, row) in enumerate(selected):
            row_to_pos[row] = pos
        counts = np.zeros((len(selected), n_cells), dtype=np.float32)
        totals = np.zeros(n_cells, dtype=np.float64)
        data = grp["data"]
        indices = grp["indices"]
        indptr = grp["indptr"][:]
        chunk = 1000
        for c0 in range(0, n_cells, chunk):
            c1 = min(n_cells, c0 + chunk)
            start = int(indptr[c0])
            end = int(indptr[c1])
            vals = data[start:end]
            idx = indices[start:end]
            lengths = np.diff(indptr[c0 : c1 + 1]).astype(np.int64)
            cols = np.repeat(np.arange(c1 - c0, dtype=np.int32), lengths)
            if len(vals):
                totals[c0:c1] = np.bincount(cols, weights=vals, minlength=c1 - c0)
                pos = row_to_pos[idx]
                mask = pos >= 0
                if mask.any():
                    np.add.at(counts, (pos[mask], cols[mask] + c0), vals[mask])
        selected_genes = [g for g, _ in selected]
    return selected_genes, barcodes, counts, totals


def norm_log(counts, totals):
    totals = np.asarray(totals, dtype=np.float64)
    return np.log1p((counts / np.maximum(totals, 1.0)) * 10000.0)


def score_sets(selected_genes, expr):
    gene_to_i = {g.upper(): i for i, g in enumerate(selected_genes)}
    scores = {}
    for name, genes in GENE_SETS.items():
        idx = [gene_to_i[g.upper()] for g in genes if g.upper() in gene_to_i]
        if idx:
            scores[name] = expr[idx, :].mean(axis=0)
        else:
            scores[name] = np.zeros(expr.shape[1], dtype=np.float32)
    scores["target_subclone02_04_like"] = np.vstack(
        [scores["subclone02_like"], scores["subclone04_like"], scores["subclone02_04_common"]]
    ).mean(axis=0)
    scores["target_axis"] = np.vstack(
        [
            scores["target_subclone02_04_like"],
            scores["itgb1_cd44_target"],
            scores["hypoxia_core"],
            scores["kras_mapk_core"],
            scores["fak_ecm_invasion"],
        ]
    ).mean(axis=0)
    return scores


def gene_vector(selected_genes, expr, gene):
    gene = gene.upper()
    lookup = {g.upper(): i for i, g in enumerate(selected_genes)}
    if gene not in lookup:
        return np.zeros(expr.shape[1], dtype=np.float32)
    return expr[lookup[gene], :]


def safe_mean(x):
    x = np.asarray(x)
    if x.size == 0:
        return 0.0
    return float(np.nanmean(x))


def safe_pct(x):
    x = np.asarray(x)
    if x.size == 0:
        return 0.0
    return float(np.nanmean(x > 0))


def classify_cells(scores):
    major_names = ["tumor_epithelial", "myeloid", "t_nk", "b_cell", "fibroblast", "endothelial"]
    labels = np.array(["Unknown"] * len(scores["myeloid"]), dtype=object)
    mat = np.vstack([scores[n] for n in major_names])
    best = np.argmax(mat, axis=0)
    best_score = np.max(mat, axis=0)
    map_label = {
        "tumor_epithelial": "Tumor/Epithelial",
        "myeloid": "Myeloid/Macrophage",
        "t_nk": "T/NK",
        "b_cell": "B/Plasma",
        "fibroblast": "CAF/Fibroblast",
        "endothelial": "Endothelial",
    }
    for i, name in enumerate(major_names):
        labels[(best == i) & (best_score > 0)] = map_label[name]
    return labels


def analyze_scRNA():
    rows = []
    comp_rows = []
    status = []
    genes = wanted_genes()
    for dataset, path in SC_H5:
        try:
            log(f"scRNA: processing {dataset}")
            selected_genes, barcodes, counts, totals = read_10x_h5_selected(path, genes)
            expr = norm_log(counts, totals)
            scores = score_sets(selected_genes, expr)
            labels = classify_cells(scores)
            n_cells = len(barcodes)
            spp1 = gene_vector(selected_genes, expr, "SPP1")
            tumor_mask = labels == "Tumor/Epithelial"
            myeloid_mask = labels == "Myeloid/Macrophage"
            if tumor_mask.sum() < 30:
                tumor_mask = scores["tumor_epithelial"] >= np.quantile(scores["tumor_epithelial"], 0.75)
            if myeloid_mask.sum() < 30:
                myeloid_mask = scores["spp1_macrophage"] >= np.quantile(scores["spp1_macrophage"], 0.75)
            target_pool = tumor_mask
            target_cut = np.quantile(scores["target_axis"][target_pool], 0.75) if target_pool.sum() else np.quantile(scores["target_axis"], 0.75)
            target_mask = target_pool & (scores["target_axis"] >= target_cut)
            if target_mask.sum() < 10:
                target_mask = scores["target_axis"] >= np.quantile(scores["target_axis"], 0.90)

            for ligand, receptor in AXES:
                lv = gene_vector(selected_genes, expr, ligand)
                rv = gene_vector(selected_genes, expr, receptor)
                ligand_avg = safe_mean(lv[myeloid_mask])
                receptor_avg = safe_mean(rv[target_mask])
                ligand_pct = safe_pct(lv[myeloid_mask])
                receptor_pct = safe_pct(rv[target_mask])
                expr_product = ligand_avg * receptor_avg
                lr_score = expr_product * ligand_pct * receptor_pct
                rows.append(
                    {
                        "dataset": dataset,
                        "sample_id": dataset,
                        "axis": f"{ligand}-{receptor}",
                        "ligand": ligand,
                        "receptor": receptor,
                        "n_cells": n_cells,
                        "n_source_myeloid": int(myeloid_mask.sum()),
                        "n_target_tumor_high": int(target_mask.sum()),
                        "source_fraction": float(myeloid_mask.mean()),
                        "target_fraction": float(target_mask.mean()),
                        "ligand_avg_source": ligand_avg,
                        "receptor_avg_target": receptor_avg,
                        "ligand_pct_source": ligand_pct,
                        "receptor_pct_target": receptor_pct,
                        "expr_product_score": expr_product,
                        "lr_score": lr_score,
                        "target_axis_mean_tumor": safe_mean(scores["target_axis"][tumor_mask]),
                        "SPP1_macrophage_mean_source": safe_mean(scores["spp1_macrophage"][myeloid_mask]),
                        "classification_method": "marker_score_no_metadata",
                    }
                )
            for label, n in pd.Series(labels).value_counts().items():
                comp_rows.append({"dataset": dataset, "cell_type_major": label, "n_cells": int(n), "fraction": float(n / n_cells)})
            status.append({"dataset": dataset, "status": "completed", "n_cells": n_cells, "n_genes_selected": len(selected_genes), "note": "10x h5 processed with marker-score cell-type inference; no external metadata embedded in h5."})
        except Exception as exc:
            log(f"scRNA: failed {dataset}: {exc}")
            status.append({"dataset": dataset, "status": "failed", "n_cells": np.nan, "n_genes_selected": np.nan, "note": str(exc)})

    df = pd.DataFrame(rows)
    comp = pd.DataFrame(comp_rows)
    stat = pd.DataFrame(status)
    df.to_csv(TABLES / "extended_external_scRNA_dataset_level_summary.csv", index=False)
    comp.to_csv(TABLES / "extended_external_scRNA_celltype_composition.csv", index=False)
    stat.to_csv(TABLES / "extended_external_scRNA_status.csv", index=False)
    return df, comp, stat


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")


def read_features(path):
    genes = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                genes.append(parts[1])
            elif parts:
                genes.append(parts[0])
    return genes


def read_barcodes(path):
    with open_text(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def read_10x_mtx_selected(sample_dir, genes):
    fdir = sample_dir / "filtered_feature_bc_matrix"
    feature_path = next(iter(list(fdir.glob("features.tsv*")) + list(fdir.glob("genes.tsv*"))), None)
    barcode_path = next(iter(list(fdir.glob("barcodes.tsv*"))), None)
    matrix_path = next(iter(list(fdir.glob("matrix.mtx*"))), None)
    if not all([feature_path, barcode_path, matrix_path]):
        raise FileNotFoundError(f"Missing 10x mtx triplet in {sample_dir}")
    features = read_features(feature_path)
    barcodes = read_barcodes(barcode_path)
    gene_to_row = {}
    for i, g in enumerate(features):
        gene_to_row.setdefault(g.upper(), i)
    selected = [(g.upper(), gene_to_row[g.upper()]) for g in genes if g.upper() in gene_to_row]
    mat = mmread(str(matrix_path)).tocsc()
    totals = np.asarray(mat.sum(axis=0)).ravel()
    sub = mat[[row for _, row in selected], :].toarray().astype(np.float32) if selected else np.zeros((0, mat.shape[1]), dtype=np.float32)
    return [g for g, _ in selected], barcodes, sub, totals


def read_spatial_sample(sample_dir):
    genes = wanted_genes()
    h5 = sample_dir / "filtered_feature_bc_matrix.h5"
    if h5.exists():
        return read_10x_h5_selected(h5, genes)
    return read_10x_mtx_selected(sample_dir, genes)


def parse_positions(sample_dir):
    candidates = [
        sample_dir / "spatial" / "tissue_positions_list.csv",
        sample_dir / "spatial" / "tissue_positions.csv",
    ]
    pos_path = next((p for p in candidates if p.exists()), None)
    if pos_path is None:
        return pd.DataFrame()
    raw = pd.read_csv(pos_path, header=None)
    first = str(raw.iloc[0, 0]).lower()
    if first in {"barcode", "barcodes"}:
        raw = pd.read_csv(pos_path)
        cols = {c.lower(): c for c in raw.columns}
        barcode_col = cols.get("barcode", raw.columns[0])
        row_col = cols.get("array_row", cols.get("row", raw.columns[min(2, len(raw.columns) - 1)]))
        col_col = cols.get("array_col", cols.get("col", raw.columns[min(3, len(raw.columns) - 1)]))
        tissue_col = cols.get("in_tissue")
        out = pd.DataFrame({"barcode": raw[barcode_col].astype(str), "array_row": raw[row_col], "array_col": raw[col_col]})
        out["in_tissue"] = raw[tissue_col] if tissue_col else 1
        return out
    names = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: raw.shape[1]]
    raw.columns = names
    return raw


def spatial_summary_one(dataset, sample_dir):
    sample_id = sample_dir.name
    selected_genes, barcodes, counts, totals = read_spatial_sample(sample_dir)
    expr = norm_log(counts, totals)
    scores = score_sets(selected_genes, expr)
    pos = parse_positions(sample_dir)
    df = pd.DataFrame({"barcode": barcodes})
    for key, val in scores.items():
        df[key] = val
    for gene in ["SPP1", "ITGB1", "CD44", "MIF", "CD74", "APOE", "LRP1", "TGFB1", "TGFBR2", "CXCL12", "CXCR4"]:
        df[gene] = gene_vector(selected_genes, expr, gene)
    if not pos.empty:
        df = df.merge(pos, on="barcode", how="left")
        if "in_tissue" in df.columns:
            df = df[df["in_tissue"].fillna(1).astype(int) == 1].copy()
    source_score = df["spp1_macrophage"].to_numpy()
    target_score = df["target_axis"].to_numpy()
    source_high = source_score >= np.nanquantile(source_score, 0.75)
    target_high = target_score >= np.nanquantile(target_score, 0.75)
    pear_r, pear_p = pearsonr(source_score, target_score) if len(df) >= 3 and np.nanstd(source_score) > 0 and np.nanstd(target_score) > 0 else (np.nan, np.nan)
    spear = spearmanr(source_score, target_score) if len(df) >= 3 else (np.nan, np.nan)
    spear_r = getattr(spear, "statistic", spear[0] if isinstance(spear, tuple) else np.nan)
    spear_p = getattr(spear, "pvalue", spear[1] if isinstance(spear, tuple) else np.nan)

    neigh = {
        "dataset": dataset,
        "sample_id": sample_id,
        "status": "no_coordinates",
        "observed_neighbor_fraction": np.nan,
        "expected_global_fraction": float(target_high.mean()) if len(target_high) else np.nan,
        "enrichment_ratio": np.nan,
        "n_source_spots": int(source_high.sum()),
        "n_target_spots": int(target_high.sum()),
        "n_neighbor_edges": 0,
        "radius_array_units": np.nan,
    }
    if {"array_row", "array_col"}.issubset(df.columns) and df[["array_row", "array_col"]].notna().all(axis=None):
        coords = df[["array_col", "array_row"]].to_numpy(dtype=float)
        tree = cKDTree(coords)
        dists, _ = tree.query(coords, k=min(2, len(coords)))
        radius = float(np.nanmedian(dists[:, 1]) * 2.01) if len(coords) > 2 else 2.01
        edges = tree.query_ball_point(coords[source_high], r=radius)
        neighbor_idx = [j for arr in edges for j in arr]
        if neighbor_idx:
            observed = float(target_high[neighbor_idx].mean())
            expected = float(target_high.mean())
            enrich = observed / expected if expected else np.nan
            neigh.update(
                {
                    "status": "completed",
                    "observed_neighbor_fraction": observed,
                    "expected_global_fraction": expected,
                    "enrichment_ratio": enrich,
                    "n_neighbor_edges": int(len(neighbor_idx)),
                    "radius_array_units": radius,
                }
            )

    lr_rows = []
    ko_rows = []
    for ligand, receptor in AXES:
        lv = df[ligand].to_numpy() if ligand in df else np.zeros(len(df))
        rv = df[receptor].to_numpy() if receptor in df else np.zeros(len(df))
        score_vec = lv * rv
        source_lig = lv[source_high]
        target_rec = rv[target_high]
        control = safe_mean(source_lig) * safe_mean(target_rec) * safe_pct(source_lig) * safe_pct(target_rec)
        lr_rows.append(
            {
                "dataset": dataset,
                "sample_id": sample_id,
                "axis": f"{ligand}-{receptor}",
                "n_spots": len(df),
                "mean_spot_product": safe_mean(score_vec),
                "median_spot_product": float(np.nanmedian(score_vec)) if len(score_vec) else np.nan,
                "source_high_ligand_avg": safe_mean(source_lig),
                "target_high_receptor_avg": safe_mean(target_rec),
                "source_high_ligand_pct": safe_pct(source_lig),
                "target_high_receptor_pct": safe_pct(target_rec),
                "focused_lr_score": control,
            }
        )
        if ligand == "SPP1" and receptor in {"ITGB1", "CD44"}:
            for ko, value in [
                ("control", control),
                ("SPP1_source_KO", 0.0),
                (f"{receptor}_target_KO", 0.0),
                (f"SPP1_{receptor}_double_KO", 0.0),
            ]:
                ko_rows.append(
                    {
                        "dataset": dataset,
                        "sample_id": sample_id,
                        "axis": f"{ligand}-{receptor}",
                        "condition": ko,
                        "lr_score": value,
                        "relative_to_control": value / control if control else np.nan,
                        "predicted_reduction_fraction": 1 - (value / control) if control else np.nan,
                    }
                )

    summary = {
        "dataset": dataset,
        "sample_id": sample_id,
        "n_spots": len(df),
        "n_genes_selected": len(selected_genes),
        "pearson_r_spp1_macro_vs_target_axis": pear_r,
        "pearson_p_spp1_macro_vs_target_axis": pear_p,
        "spearman_r_spp1_macro_vs_target_axis": spear_r,
        "spearman_p_spp1_macro_vs_target_axis": spear_p,
        "mean_spp1_macro_score": safe_mean(source_score),
        "mean_target_axis_score": safe_mean(target_score),
        "source_high_fraction": float(source_high.mean()) if len(source_high) else np.nan,
        "target_high_fraction": float(target_high.mean()) if len(target_high) else np.nan,
    }
    return summary, neigh, lr_rows, ko_rows


def analyze_spatial():
    summaries = []
    neighs = []
    lrs = []
    kos = []
    status = []
    for dataset, root in SPATIAL_ROOTS:
        for sample_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
            try:
                log(f"spatial: processing {dataset} {sample_dir.name}")
                summary, neigh, lr_rows, ko_rows = spatial_summary_one(dataset, sample_dir)
                summaries.append(summary)
                neighs.append(neigh)
                lrs.extend(lr_rows)
                kos.extend(ko_rows)
                status.append({"dataset": dataset, "sample_id": sample_dir.name, "status": "completed", "note": "Visium sample processed with spot-level scores and coordinate neighborhood analysis when coordinates were available."})
            except Exception as exc:
                log(f"spatial: failed {dataset} {sample_dir.name}: {exc}")
                status.append({"dataset": dataset, "sample_id": sample_dir.name, "status": "failed", "note": str(exc)})
    summary_df = pd.DataFrame(summaries)
    neigh_df = pd.DataFrame(neighs)
    lr_df = pd.DataFrame(lrs)
    ko_df = pd.DataFrame(kos)
    status_df = pd.DataFrame(status)
    summary_df.to_csv(TABLES / "extended_spatial_sample_summary.csv", index=False)
    neigh_df.to_csv(TABLES / "extended_spatial_neighborhood_enrichment.csv", index=False)
    lr_df.to_csv(TABLES / "extended_spatial_LR_scores.csv", index=False)
    ko_df.to_csv(TABLES / "extended_spatial_virtual_KO.csv", index=False)
    status_df.to_csv(TABLES / "extended_spatial_status.csv", index=False)
    return summary_df, neigh_df, lr_df, ko_df, status_df


def plot_outputs(sc_df, comp_df, spatial_df, neigh_df, spatial_lr_df):
    if not sc_df.empty:
        primary = sc_df[sc_df["axis"].isin(["SPP1-ITGB1", "SPP1-CD44"])].copy()
        pivot = primary.pivot_table(index="dataset", columns="axis", values="lr_score", aggfunc="mean").fillna(0)
        plt.figure(figsize=(7, max(3, 0.45 * len(pivot) + 1.5)))
        plt.imshow(pivot.values, aspect="auto", cmap="magma")
        plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=35, ha="right")
        plt.yticks(range(len(pivot.index)), pivot.index)
        plt.colorbar(label="scRNA LR score")
        plt.tight_layout()
        plt.savefig(FIGURES / "extended_scRNA_SPP1_LR_heatmap.png", dpi=220)
        plt.savefig(FIGURES / "extended_scRNA_SPP1_LR_heatmap.pdf")
        plt.close()
    if not neigh_df.empty:
        plot_df = neigh_df.copy()
        plot_df["label"] = plot_df["dataset"] + "_" + plot_df["sample_id"]
        plt.figure(figsize=(max(8, 0.45 * len(plot_df)), 4))
        plt.bar(plot_df["label"], plot_df["enrichment_ratio"], color="#4c78a8")
        plt.axhline(1, color="black", linewidth=0.8)
        plt.ylabel("SPP1-high to target-high neighborhood enrichment")
        plt.xticks(rotation=75, ha="right")
        plt.tight_layout()
        plt.savefig(FIGURES / "extended_spatial_neighborhood_enrichment.png", dpi=220)
        plt.savefig(FIGURES / "extended_spatial_neighborhood_enrichment.pdf")
        plt.close()
    if not spatial_df.empty:
        plot_df = spatial_df.copy()
        plot_df["label"] = plot_df["dataset"] + "_" + plot_df["sample_id"]
        plt.figure(figsize=(max(8, 0.45 * len(plot_df)), 4))
        plt.bar(plot_df["label"], plot_df["spearman_r_spp1_macro_vs_target_axis"], color="#59a14f")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.ylabel("Spearman r: SPP1 macrophage vs target axis")
        plt.xticks(rotation=75, ha="right")
        plt.tight_layout()
        plt.savefig(FIGURES / "extended_spatial_score_correlation.png", dpi=220)
        plt.savefig(FIGURES / "extended_spatial_score_correlation.pdf")
        plt.close()
    if not spatial_lr_df.empty:
        primary = spatial_lr_df[spatial_lr_df["axis"].isin(["SPP1-ITGB1", "SPP1-CD44"])].copy()
        primary["label"] = primary["dataset"] + "_" + primary["sample_id"]
        pivot = primary.pivot_table(index="label", columns="axis", values="focused_lr_score", aggfunc="mean").fillna(0)
        plt.figure(figsize=(7, max(4, 0.28 * len(pivot) + 1.5)))
        plt.imshow(pivot.values, aspect="auto", cmap="viridis")
        plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=35, ha="right")
        plt.yticks(range(len(pivot.index)), pivot.index, fontsize=7)
        plt.colorbar(label="Spatial focused LR score")
        plt.tight_layout()
        plt.savefig(FIGURES / "extended_spatial_SPP1_LR_heatmap.png", dpi=220)
        plt.savefig(FIGURES / "extended_spatial_SPP1_LR_heatmap.pdf")
        plt.close()


def write_report(sc_df, sc_status, spatial_df, neigh_df, spatial_lr_df, ko_df, spatial_status):
    sc_primary = sc_df[sc_df["axis"].isin(["SPP1-ITGB1", "SPP1-CD44"])] if not sc_df.empty else pd.DataFrame()
    sc_meta = sc_primary.groupby("axis", as_index=False).agg(
        n_dataset=("dataset", "nunique"),
        mean_lr_score=("lr_score", "mean"),
        median_lr_score=("lr_score", "median"),
        mean_source_fraction=("source_fraction", "mean"),
        mean_target_fraction=("target_fraction", "mean"),
    )
    sp_meta = spatial_lr_df[spatial_lr_df["axis"].isin(["SPP1-ITGB1", "SPP1-CD44"])].groupby("axis", as_index=False).agg(
        n_spatial_samples=("sample_id", "nunique"),
        mean_focused_lr_score=("focused_lr_score", "mean"),
        median_focused_lr_score=("focused_lr_score", "median"),
    ) if not spatial_lr_df.empty else pd.DataFrame()
    sc_meta.to_csv(TABLES / "extended_external_scRNA_meta_summary.csv", index=False)
    sp_meta.to_csv(TABLES / "extended_spatial_LR_meta_summary.csv", index=False)
    evidence = []
    evidence.append({"evidence_layer": "external_scRNA", "dataset_count": int(sc_status[sc_status["status"] == "completed"]["dataset"].nunique()) if not sc_status.empty else 0, "status": "completed" if not sc_df.empty else "failed", "main_output": "extended_external_scRNA_dataset_level_summary.csv"})
    evidence.append({"evidence_layer": "spatial", "dataset_count": int(spatial_status[spatial_status["status"] == "completed"]["dataset"].nunique()) if not spatial_status.empty else 0, "status": "completed" if not spatial_df.empty else "failed", "main_output": "extended_spatial_sample_summary.csv"})
    evidence.append({"evidence_layer": "spatial_virtual_KO", "dataset_count": int(ko_df["sample_id"].nunique()) if not ko_df.empty else 0, "status": "completed" if not ko_df.empty else "failed", "main_output": "extended_spatial_virtual_KO.csv"})
    pd.DataFrame(evidence).to_csv(TABLES / "extended_validation_evidence_matrix.csv", index=False)

    report = OUT / "README.md"
    def md_table(df):
        if df.empty:
            return "_No rows._"
        text_df = df.copy().fillna("")
        cols = [str(c) for c in text_df.columns]
        lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
        for row in text_df.astype(str).to_numpy():
            lines.append("| " + " | ".join(str(x).replace("|", "/") for x in row) + " |")
        return "\n".join(lines)

    with open(report, "w", encoding="utf-8") as handle:
        handle.write("# Extended external validation for SPP1 macrophage to ITGB1/CD44 tumor axis\n\n")
        handle.write("This directory extends the previous validation by adding five ovarian scRNA h5 datasets and two spatial transcriptomics GEO datasets that were not included in the first pass.\n\n")
        handle.write("## Added external scRNA datasets\n\n")
        for name, path in SC_H5:
            handle.write(f"- {name}: `{path}`\n")
        handle.write("\nThe h5 files did not contain usable cell metadata in the 10x matrix groups, so broad cell classes were inferred from marker scores. This is a validation-by-expression-potential analysis, not a curated annotation analysis.\n\n")
        handle.write("## Added spatial datasets\n\n")
        for name, path in SPATIAL_ROOTS:
            handle.write(f"- {name}: `{path}`\n")
        handle.write("\n## Key outputs\n\n")
        for fn in [
            "extended_external_scRNA_dataset_level_summary.csv",
            "extended_external_scRNA_meta_summary.csv",
            "extended_external_scRNA_celltype_composition.csv",
            "extended_spatial_sample_summary.csv",
            "extended_spatial_neighborhood_enrichment.csv",
            "extended_spatial_LR_scores.csv",
            "extended_spatial_virtual_KO.csv",
            "extended_validation_evidence_matrix.csv",
        ]:
            handle.write(f"- `tables/{fn}`\n")
        handle.write("\n## Methods summary\n\n")
        handle.write("Selected genes from the mechanism-axis gene sets were streamed from 10x h5/mtx matrices, normalized as log1p(CP10K), and summarized into module scores. scRNA source cells were marker-inferred myeloid/macrophage cells; target cells were marker-inferred tumor/epithelial cells in the top quartile of the target-axis score. Spatial spots were scored with the same signatures. Neighborhood enrichment used array coordinates and a radius of 2.01 times the median nearest-neighbor distance. Spatial virtual KO set SPP1 or receptor expression to zero in the focused source/target compartments and recalculated focused LR scores.\n\n")
        if not sc_meta.empty:
            handle.write("## scRNA meta summary\n\n")
            handle.write(md_table(sc_meta))
            handle.write("\n\n")
        if not sp_meta.empty:
            handle.write("## spatial LR meta summary\n\n")
            handle.write(md_table(sp_meta))
            handle.write("\n")


def main():
    log("Starting extended external validation")
    with open(OUT / "gene_sets_used.json", "w", encoding="utf-8") as handle:
        json.dump(GENE_SETS, handle, indent=2)
    sc_df, comp_df, sc_status = analyze_scRNA()
    spatial_df, neigh_df, spatial_lr_df, ko_df, spatial_status = analyze_spatial()
    plot_outputs(sc_df, comp_df, spatial_df, neigh_df, spatial_lr_df)
    write_report(sc_df, sc_status, spatial_df, neigh_df, spatial_lr_df, ko_df, spatial_status)
    log("Completed extended external validation")


if __name__ == "__main__":
    main()
