from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plan_dir = Path(r"D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis")
tables_dir = plan_dir / "tables"
figures_dir = plan_dir / "figures"

k = 30
n_perm = 1000
rng = np.random.default_rng(20260707)

print("Loading PCA and metadata...")
pca = pd.read_csv(tables_dir / "integrated_oc_pca30_for_knn.csv")
meta = pd.read_csv(tables_dir / "integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv")
meta = meta.set_index("cell_integrated_oc").loc[pca["cell"]].reset_index()

X = pca.drop(columns=["cell"]).to_numpy(dtype=np.float32)
interaction = meta["interaction_group"].fillna("Unassigned").astype(str).to_numpy()
broad = meta["cell_type"].fillna("Unassigned").astype(str).to_numpy()
clone = meta["cnv_subclone"].fillna("").astype(str).to_numpy()

print("Computing nearest neighbors...")
nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", algorithm="auto", n_jobs=-1)
nn.fit(X)
indices = nn.kneighbors(X, return_distance=False)
indices = indices[:, 1:(k + 1)]

def enrichment_for_label(label_array, label_name):
    levels = pd.Index(pd.Series(label_array).unique()).sort_values().tolist()
    code_map = {v: i for i, v in enumerate(levels)}
    codes = np.array([code_map[v] for v in label_array], dtype=np.int32)
    global_counts = np.bincount(codes, minlength=len(levels))
    global_prop = global_counts / global_counts.sum()

    rows = []
    clone_levels = sorted([x for x in pd.Series(clone).unique() if x])
    for cl in clone_levels:
      clone_idx = np.where(clone == cl)[0]
      flat = indices[clone_idx, :].reshape(-1)
      obs_counts = np.bincount(codes[flat], minlength=len(levels))
      obs_prop = obs_counts / obs_counts.sum()

      perm_props = np.zeros((n_perm, len(levels)), dtype=np.float32)
      for p in range(n_perm):
        shuffled = rng.permutation(codes)
        pc = np.bincount(shuffled[flat], minlength=len(levels))
        perm_props[p, :] = pc / pc.sum()

      perm_mean = perm_props.mean(axis=0)
      perm_sd = perm_props.std(axis=0, ddof=1)
      z = (obs_prop - perm_mean) / np.where(perm_sd == 0, np.nan, perm_sd)
      p_greater = (np.sum(perm_props >= obs_prop[None, :], axis=0) + 1) / (n_perm + 1)
      p_less = (np.sum(perm_props <= obs_prop[None, :], axis=0) + 1) / (n_perm + 1)
      enrich = obs_prop / np.where(perm_mean == 0, np.nan, perm_mean)

      for i, lab in enumerate(levels):
        rows.append({
          "label_type": label_name,
          "cnv_subclone": cl,
          "neighbor_label": lab,
          "clone_n_cells": len(clone_idx),
          "neighbor_n": int(obs_counts[i]),
          "observed_prop": float(obs_prop[i]),
          "expected_prop_perm_mean": float(perm_mean[i]),
          "enrichment_ratio": float(enrich[i]) if np.isfinite(enrich[i]) else np.nan,
          "z_score": float(z[i]) if np.isfinite(z[i]) else np.nan,
          "empirical_p_greater": float(p_greater[i]),
          "empirical_p_less": float(p_less[i]),
          "global_prop": float(global_prop[i]),
        })
    return pd.DataFrame(rows)

print("Running enrichment tests for interaction_group...")
fine = enrichment_for_label(interaction, "interaction_group")
fine.to_csv(tables_dir / "knn_neighbor_enrichment_interaction_group_k30_perm1000.csv", index=False)

print("Running enrichment tests for broad cell_type...")
broad_df = enrichment_for_label(broad, "cell_type")
broad_df.to_csv(tables_dir / "knn_neighbor_enrichment_cell_type_k30_perm1000.csv", index=False)

focus_patterns = [
    "Myeloid_Macro-M2-like", "Myeloid_Macro-LYVE1", "Myeloid_Macro-THBS1",
    "Myeloid_Macro-C3/CX3CR1", "Myeloid_Macro-Inflammatory_M1",
    "Myeloid_Macro-Inflammatory_TNF", "Myeloid_Macro-TIMD4",
    "Smooth_muscle_cells", "B_", "T_NK_CD8", "T_NK_CD4", "T_NK_CD56", "T_NK_",
    "DC", "Myeloid_cDC", "CNV_Subclone"
]
focus_mask = fine["neighbor_label"].apply(lambda x: any(str(x).startswith(p) or str(x) == p for p in focus_patterns))
focus = fine[focus_mask].sort_values(["cnv_subclone", "empirical_p_greater", "enrichment_ratio"], ascending=[True, True, False])
focus.to_csv(tables_dir / "knn_neighbor_enrichment_focus_tme_groups.csv", index=False)

def plot_heatmap(df, value, out_prefix, top_n=None):
    mat_df = df.copy()
    if top_n is not None:
        max_abs = mat_df.groupby("neighbor_label")[value].apply(lambda s: np.nanmax(np.abs(s)))
        keep = max_abs.sort_values(ascending=False).head(top_n).index
        mat_df = mat_df[mat_df["neighbor_label"].isin(keep)]
    mat = mat_df.pivot(index="neighbor_label", columns="cnv_subclone", values=value)
    mat = mat.replace([np.inf, -np.inf], np.nan).fillna(0)
    plt.figure(figsize=(7, max(4, 0.22 * mat.shape[0] + 1.8)))
    sns.heatmap(mat, cmap="vlag", center=0 if value == "z_score" else 1,
                linewidths=0.1, linecolor="white", cbar_kws={"label": value})
    plt.title(out_prefix.replace("_", " "))
    plt.xlabel("CNV subclone")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(figures_dir / f"{out_prefix}.png", dpi=220)
    plt.savefig(figures_dir / f"{out_prefix}.pdf")
    plt.close()

print("Plotting heatmaps...")
plot_heatmap(fine, "z_score", "knn_interaction_group_neighbor_zscore_top35", top_n=35)
plot_heatmap(fine, "enrichment_ratio", "knn_interaction_group_neighbor_enrichment_top35", top_n=35)
plot_heatmap(broad_df, "z_score", "knn_broad_cell_type_neighbor_zscore")
plot_heatmap(broad_df, "enrichment_ratio", "knn_broad_cell_type_neighbor_enrichment")

top_enriched = fine.sort_values(["empirical_p_greater", "z_score"], ascending=[True, False]).head(200)
top_enriched.to_csv(tables_dir / "knn_top_enriched_neighbor_labels.csv", index=False)

print("Done.")
