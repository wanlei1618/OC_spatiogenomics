from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plan_dir = Path(r"D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis")
tables_dir = plan_dir / "tables"
figures_dir = plan_dir / "figures"

fine = pd.read_csv(tables_dir / "knn_neighbor_enrichment_interaction_group_k30_perm1000.csv")
tme = fine[~fine["neighbor_label"].astype(str).str.startswith("CNV_Subclone")].copy()
tme = tme.sort_values(["empirical_p_greater", "z_score"], ascending=[True, False])
tme.to_csv(tables_dir / "knn_neighbor_enrichment_TME_only_interaction_group.csv", index=False)
tme.head(200).to_csv(tables_dir / "knn_top_TME_enriched_neighbor_labels.csv", index=False)

focus_terms = [
    "Myeloid_Macro-M2-like", "Myeloid_Macro-LYVE1", "Myeloid_Macro-THBS1",
    "Myeloid_Macro-C3/CX3CR1", "Myeloid_Macro-Inflammatory_M1",
    "Myeloid_Macro-Inflammatory_TNF", "Myeloid_Macro-TIMD4", "Macrophages",
    "Smooth_muscle_cells", "B_Bn_TCL1A", "B_Classical-Bm_TXNIP", "B_PC_IGHG",
    "B_Bm_stress-response", "B_Early-PC_MS4A1low", "T_NK_CD8+ T cytotoxic",
    "T_NK_CD8+ T effector", "T_NK_CD8+ T resident memory", "T_NK_CD8+ T CXCL13",
    "T_NK_CD4+ T regulatory", "T_NK_CD4+ T effector", "T_NK_CD4+ T naive",
    "T_NK_CD56dim cytotoxic NK", "T_NK_CD56bright regulatory", "DC"
]
focus = tme[tme["neighbor_label"].isin(focus_terms)].copy()
focus.to_csv(tables_dir / "knn_neighbor_enrichment_TME_focus_groups.csv", index=False)

def plot(df, value, outfile, title):
    mat = df.pivot(index="neighbor_label", columns="cnv_subclone", values=value)
    mat = mat.replace([np.inf, -np.inf], np.nan).fillna(0)
    plt.figure(figsize=(7, max(5, 0.3 * mat.shape[0] + 2)))
    sns.heatmap(mat, cmap="vlag", center=0 if value == "z_score" else 1,
                linewidths=0.1, linecolor="white", cbar_kws={"label": value})
    plt.title(title)
    plt.xlabel("CNV subclone")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(figures_dir / f"{outfile}.png", dpi=220)
    plt.savefig(figures_dir / f"{outfile}.pdf")
    plt.close()

plot(focus, "z_score", "knn_TME_focus_neighbor_zscore", "TME-focused kNN neighbor enrichment z-score")
plot(focus, "enrichment_ratio", "knn_TME_focus_neighbor_enrichment", "TME-focused kNN neighbor enrichment ratio")
print("Done TME-only outputs")
