import csv
import math
import os
from collections import defaultdict
from datetime import datetime


PROJECT_ROOT = r"D:\OC_spatiogenomics\infercnv"
OUT_ROOT = os.path.join(PROJECT_ROOT, "06_rigorous_validation")
SRC_ROOT = os.path.join(PROJECT_ROOT, "05_sample_type_lr_niche_external_validation", "sample_type_LR_niche_analysis")
TABLE_SRC = os.path.join(SRC_ROOT, "tables")

DIRS = {
    "external": os.path.join(OUT_ROOT, "06_external_scrna_projection"),
    "lr": os.path.join(OUT_ROOT, "07_ligand_receptor_competition"),
    "nichenet": os.path.join(OUT_ROOT, "08_nichenet_ligand_target"),
    "spatial": os.path.join(OUT_ROOT, "09_spatial_deconvolution_neighborhood"),
    "perturb": os.path.join(OUT_ROOT, "10_virtual_perturbation_causal"),
    "tables": os.path.join(OUT_ROOT, "tables"),
    "reports": os.path.join(OUT_ROOT, "reports"),
    "logs": os.path.join(OUT_ROOT, "logs"),
}
for path in DIRS.values():
    os.makedirs(path, exist_ok=True)

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_PATH = os.path.join(DIRS["logs"], f"11_execute_p1_from_existing_tables_{TIMESTAMP}.log")


def read_csv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def fnum(value):
    try:
        if value in ("", None):
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def mean(values):
    values = [v for v in values if not math.isnan(v)]
    return sum(values) / len(values) if values else math.nan


def median(values):
    values = sorted(v for v in values if not math.isnan(v))
    if not values:
        return math.nan
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2


def rankdata(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg_rank = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg_rank
        i = j + 1
    return ranks


def pearson(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if not math.isnan(a) and not math.isnan(b)]
    n = len(pairs)
    if n < 3:
        return math.nan, n
    xs, ys = zip(*pairs)
    mx, my = mean(xs), mean(ys)
    sx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    sy = math.sqrt(sum((b - my) ** 2 for b in ys))
    if sx == 0 or sy == 0:
        return math.nan, n
    return sum((a - mx) * (b - my) for a, b in pairs) / (sx * sy), n


def normal_p_from_z(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def spearman(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if not math.isnan(a) and not math.isnan(b)]
    n = len(pairs)
    if n < 3:
        return math.nan, math.nan, n
    xs, ys = zip(*pairs)
    rx = rankdata(list(xs))
    ry = rankdata(list(ys))
    r, _ = pearson(rx, ry)
    if math.isnan(r) or abs(r) >= 1 or n <= 3:
        return r, math.nan, n
    # Fisher-z normal approximation; used only as a compact screening statistic.
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(max(n - 3, 1))
    return r, normal_p_from_z(z), n


def bh_adjust(rows, p_col="p_value", out_col="p_adj_BH"):
    indexed = [(i, fnum(row.get(p_col))) for i, row in enumerate(rows)]
    valid = sorted([(i, p) for i, p in indexed if not math.isnan(p)], key=lambda x: x[1])
    m = len(valid)
    adjusted = [math.nan] * len(rows)
    running = 1.0
    for rank, (i, p) in reversed(list(enumerate(valid, start=1))):
        running = min(running, p * m / rank)
        adjusted[i] = running
    for i, row in enumerate(rows):
        row[out_col] = adjusted[i]


def random_effects_meta(effects, variances):
    data = [(e, v) for e, v in zip(effects, variances) if not math.isnan(e) and not math.isnan(v) and v > 0]
    k = len(data)
    if k == 0:
        return {"k": 0}
    eff = [e for e, _ in data]
    var = [v for _, v in data]
    w = [1 / v for v in var]
    fixed = sum(wi * ei for wi, ei in zip(w, eff)) / sum(w)
    q = sum(wi * (ei - fixed) ** 2 for wi, ei in zip(w, eff))
    c = sum(w) - sum(wi ** 2 for wi in w) / sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    wr = [1 / (v + tau2) for v in var]
    pooled = sum(wi * ei for wi, ei in zip(wr, eff)) / sum(wr)
    se = math.sqrt(1 / sum(wr))
    z = pooled / se if se > 0 else math.nan
    return {
        "k": k,
        "pooled_effect": pooled,
        "se": se,
        "ci_low": pooled - 1.96 * se,
        "ci_high": pooled + 1.96 * se,
        "p_value": normal_p_from_z(z) if not math.isnan(z) else math.nan,
        "tau2": tau2,
        "I2": max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0,
    }


def fmt(x, digits=4):
    if isinstance(x, str):
        return x
    if x is None or math.isnan(x):
        return "NA"
    if abs(x) >= 1000 or (abs(x) < 0.001 and x != 0):
        return f"{x:.3e}"
    return f"{x:.{digits}f}"


with open(LOG_PATH, "w", encoding="utf-8") as log:
    log.write(f"Started {datetime.now().isoformat()}\n")
    log.write(f"Source tables: {TABLE_SRC}\n")

patient_summary = read_csv(os.path.join(TABLE_SRC, "external_scRNA_patient_sampletype_summary.csv"))
signature_rows = read_csv(os.path.join(TABLE_SRC, "external_scRNA_target_signature_scores.csv"))
external_axis = read_csv(os.path.join(TABLE_SRC, "external_scRNA_all_datasets_sampletype_axis_scores.csv"))
integrated_lr = read_csv(os.path.join(TABLE_SRC, "sample_id_LR_opportunity_scores_all_axes.csv"))
spatial_corr = read_csv(os.path.join(TABLE_SRC, "spatial_correlation_SPP1_myeloid_Target_axis.csv"))
spatial_neigh = read_csv(os.path.join(TABLE_SRC, "spatial_neighborhood_enrichment.csv"))
virtual_ko = read_csv(os.path.join(TABLE_SRC, "virtual_KO_LR_score_reduction_summary_by_sample_type.csv"))

# P1-1: external patient-level association.
sample_records = {}
for row in patient_summary:
    key = row["sample_id"]
    sample_records[key] = {
        "dataset": row.get("dataset", ""),
        "patient_id": row.get("patient_id", ""),
        "sample_id": key,
        "sample_type": row.get("sample_type", ""),
        "anatomical_location": row.get("anatomical_location", ""),
        "treatment_phase": row.get("treatment_phase", ""),
        "n_cells": fnum(row.get("n_cells")),
        "source_myeloid_fraction": fnum(row.get("source_myeloid_fraction")),
        "target_tumor_fraction": fnum(row.get("target_tumor_fraction")),
    }
for row in signature_rows:
    key = row["sample_id"]
    if key not in sample_records:
        continue
    sig = row.get("signature", "")
    sample_records[key][f"{sig}_mean_all_cells"] = fnum(row.get("mean_all_cells"))
    sample_records[key][f"{sig}_mean_target_tumor"] = fnum(row.get("mean_target_tumor"))

sample_table = list(sample_records.values())
external_sample_fields = [
    "dataset", "patient_id", "sample_id", "sample_type", "anatomical_location", "treatment_phase",
    "n_cells", "source_myeloid_fraction", "target_tumor_fraction",
    "Subclone02_like_mean_target_tumor", "Subclone04_like_mean_target_tumor",
    "Subclone02_04_common_mean_target_tumor", "CD44_ITGB1_target_mean_target_tumor",
    "KRAS_hypoxia_target_mean_target_tumor",
]
write_csv(os.path.join(DIRS["external"], "external_scrna_patient_sample_level_scores.csv"), sample_table, external_sample_fields)

assoc_pairs = [
    ("source_myeloid_fraction", "Subclone02_like_mean_target_tumor"),
    ("source_myeloid_fraction", "Subclone04_like_mean_target_tumor"),
    ("source_myeloid_fraction", "Subclone02_04_common_mean_target_tumor"),
    ("source_myeloid_fraction", "CD44_ITGB1_target_mean_target_tumor"),
    ("source_myeloid_fraction", "KRAS_hypoxia_target_mean_target_tumor"),
    ("target_tumor_fraction", "Subclone02_04_common_mean_target_tumor"),
]
external_assoc = []
for xcol, ycol in assoc_pairs:
    for stratum_name, rows in {
        "all_samples": sample_table,
        "tumor": [r for r in sample_table if r.get("sample_type") == "tumor"],
        "peritoneal_implant": [r for r in sample_table if r.get("sample_type") == "peritoneal_implant"],
    }.items():
        r, p, n = spearman([fnum(v.get(xcol)) for v in rows], [fnum(v.get(ycol)) for v in rows])
        external_assoc.append({
            "stratum": stratum_name,
            "x": xcol,
            "y": ycol,
            "n_samples": n,
            "spearman_r": r,
            "p_value": p,
            "note": "patient-sample level; site and treatment may be confounded",
        })
bh_adjust(external_assoc)
write_csv(os.path.join(DIRS["external"], "external_scrna_patient_level_associations.csv"), external_assoc,
          ["stratum", "x", "y", "n_samples", "spearman_r", "p_value", "p_adj_BH", "note"])

patient_agg = defaultdict(lambda: defaultdict(list))
for row in sample_table:
    patient = row.get("patient_id", "")
    if not patient:
        continue
    for col in ["source_myeloid_fraction", "Subclone02_04_common_mean_target_tumor",
                "CD44_ITGB1_target_mean_target_tumor", "KRAS_hypoxia_target_mean_target_tumor"]:
        patient_agg[patient][col].append(fnum(row.get(col)))
patient_rows = []
for patient, vals in sorted(patient_agg.items()):
    out = {"patient_id": patient}
    for col, values in vals.items():
        out[col] = mean(values)
    patient_rows.append(out)
write_csv(os.path.join(DIRS["external"], "external_scrna_patient_aggregated_scores.csv"), patient_rows,
          ["patient_id", "source_myeloid_fraction", "Subclone02_04_common_mean_target_tumor",
           "CD44_ITGB1_target_mean_target_tumor", "KRAS_hypoxia_target_mean_target_tumor"])

loocv_rows = []
for ycol in ["Subclone02_04_common_mean_target_tumor", "CD44_ITGB1_target_mean_target_tumor", "KRAS_hypoxia_target_mean_target_tumor"]:
    full_r, full_p, full_n = spearman([fnum(r.get("source_myeloid_fraction")) for r in patient_rows],
                                      [fnum(r.get(ycol)) for r in patient_rows])
    for left_out in [r["patient_id"] for r in patient_rows]:
        rows = [r for r in patient_rows if r["patient_id"] != left_out]
        r, p, n = spearman([fnum(v.get("source_myeloid_fraction")) for v in rows], [fnum(v.get(ycol)) for v in rows])
        loocv_rows.append({
            "outcome": ycol,
            "full_spearman_r": full_r,
            "full_n_patients": full_n,
            "left_out_patient_id": left_out,
            "loo_spearman_r": r,
            "loo_n_patients": n,
            "direction_preserved": "" if math.isnan(full_r) or math.isnan(r) else (full_r == 0 or full_r * r > 0),
        })
write_csv(os.path.join(DIRS["external"], "external_scrna_leave_one_patient_out.csv"), loocv_rows,
          ["outcome", "full_spearman_r", "full_n_patients", "left_out_patient_id", "loo_spearman_r", "loo_n_patients", "direction_preserved"])

# P1-2: LR competition.
axes_required = {"SPP1-CD44", "SPP1-ITGB1", "MIF-CD74", "APOE-LRP1", "TIMP1-CD63", "TGFB1-TGFBR1", "TGFB1-TGFBR2", "CXCL12-CXCR4"}
ext_axis_summary = []
grouped = defaultdict(list)
for row in external_axis:
    axis = row.get("axis", "")
    if axis in axes_required:
        grouped[(row.get("dataset", ""), row.get("sample_type", ""), axis)].append(fnum(row.get("axis_score")))
for (dataset, sample_type, axis), values in grouped.items():
    ext_axis_summary.append({
        "source": "external_scRNA",
        "dataset": dataset,
        "sample_type": sample_type,
        "axis": axis,
        "n_patient_samples": len([v for v in values if not math.isnan(v)]),
        "mean_axis_score": mean(values),
        "median_axis_score": median(values),
    })
for sample_type in sorted({r["sample_type"] for r in ext_axis_summary}):
    rows = [r for r in ext_axis_summary if r["sample_type"] == sample_type]
    rows.sort(key=lambda r: fnum(r["mean_axis_score"]), reverse=True)
    for idx, row in enumerate(rows, 1):
        row["rank_within_sample_type"] = idx
write_csv(os.path.join(DIRS["lr"], "external_lr_competition_axis_ranking.csv"), ext_axis_summary,
          ["source", "dataset", "sample_type", "axis", "n_patient_samples", "mean_axis_score", "median_axis_score", "rank_within_sample_type"])

integrated_axis_summary = []
int_grouped = defaultdict(list)
for row in integrated_lr:
    axis = row.get("axis", "")
    target = row.get("target_clone", "")
    if axis in axes_required and target in {"Subclone_02", "Subclone_04"}:
        int_grouped[(row.get("sample_type", ""), row.get("sample_id", ""), target, axis)].append(fnum(row.get("axis_score")))
for (sample_type, sample_id, target_clone, axis), values in int_grouped.items():
    integrated_axis_summary.append({
        "source": "integrated_oc_exploratory",
        "sample_type": sample_type,
        "sample_id": sample_id,
        "target_clone": target_clone,
        "axis": axis,
        "mean_axis_score": mean(values),
        "median_axis_score": median(values),
        "n_source_target_pairs": len([v for v in values if not math.isnan(v)]),
        "interpretation": "exploratory because integrated_oc sample_id/sample_type/batch are confounded",
    })
for key in sorted({(r["sample_id"], r["target_clone"]) for r in integrated_axis_summary}):
    rows = [r for r in integrated_axis_summary if (r["sample_id"], r["target_clone"]) == key]
    rows.sort(key=lambda r: fnum(r["mean_axis_score"]), reverse=True)
    for idx, row in enumerate(rows, 1):
        row["rank_within_sample_target"] = idx
write_csv(os.path.join(DIRS["lr"], "integrated_oc_lr_competition_axis_ranking_exploratory.csv"), integrated_axis_summary,
          ["source", "sample_type", "sample_id", "target_clone", "axis", "mean_axis_score", "median_axis_score",
           "n_source_target_pairs", "rank_within_sample_target", "interpretation"])

# P1-4: spatial sample-level meta-analysis.
spatial_rows = []
for row in spatial_corr:
    n = fnum(row.get("n_spots"))
    r = fnum(row.get("spearman_r"))
    if math.isnan(n) or math.isnan(r) or n <= 3 or abs(r) >= 1:
        z = math.nan
        var = math.nan
    else:
        z = 0.5 * math.log((1 + r) / (1 - r))
        var = 1 / (n - 3)
    spatial_rows.append({
        "dataset": row.get("dataset", ""),
        "sample_id": row.get("sample_id", ""),
        "n_spots": int(n) if not math.isnan(n) else "",
        "spearman_r": r,
        "fisher_z": z,
        "variance": var,
        "spot_level_p_original": row.get("spearman_p", ""),
        "meta_note": "effect size only; spot-level p not treated as biological replicate",
    })
write_csv(os.path.join(DIRS["spatial"], "spatial_sample_level_effects_for_meta.csv"), spatial_rows,
          ["dataset", "sample_id", "n_spots", "spearman_r", "fisher_z", "variance", "spot_level_p_original", "meta_note"])

meta_rows = []
for dataset in sorted(set(r["dataset"] for r in spatial_rows) | {"ALL"}):
    rows = spatial_rows if dataset == "ALL" else [r for r in spatial_rows if r["dataset"] == dataset]
    res = random_effects_meta([fnum(r["fisher_z"]) for r in rows], [fnum(r["variance"]) for r in rows])
    if res.get("k", 0):
        pooled_z = res["pooled_effect"]
        pooled_r = math.tanh(pooled_z)
        ci_low_r = math.tanh(res["ci_low"])
        ci_high_r = math.tanh(res["ci_high"])
    else:
        pooled_r = ci_low_r = ci_high_r = math.nan
    meta_rows.append({
        "dataset": dataset,
        "k_samples": res.get("k", 0),
        "pooled_fisher_z": res.get("pooled_effect", math.nan),
        "pooled_spearman_r": pooled_r,
        "ci_low_r": ci_low_r,
        "ci_high_r": ci_high_r,
        "p_value": res.get("p_value", math.nan),
        "tau2": res.get("tau2", math.nan),
        "I2": res.get("I2", math.nan),
        "method": "DerSimonian-Laird random-effects meta-analysis on sample-level Fisher z",
    })
bh_adjust(meta_rows)
write_csv(os.path.join(DIRS["spatial"], "spatial_random_effects_meta_analysis.csv"), meta_rows,
          ["dataset", "k_samples", "pooled_fisher_z", "pooled_spearman_r", "ci_low_r", "ci_high_r",
           "p_value", "p_adj_BH", "tau2", "I2", "method"])

neigh_rows = []
for row in spatial_neigh:
    er = fnum(row.get("enrichment_ratio"))
    source_n = fnum(row.get("n_source_spots"))
    target_n = fnum(row.get("n_target_spots"))
    if math.isnan(er) or er <= 0:
        log_er = math.nan
    else:
        log_er = math.log(er)
    variance = (1 / max(source_n, 1) + 1 / max(target_n, 1)) if not math.isnan(source_n) and not math.isnan(target_n) else math.nan
    neigh_rows.append({
        "dataset": row.get("dataset", ""),
        "sample_id": row.get("sample_id", ""),
        "enrichment_ratio": er,
        "log_enrichment_ratio": log_er,
        "variance_approx": variance,
        "n_source_spots": row.get("n_source_spots", ""),
        "n_target_spots": row.get("n_target_spots", ""),
    })
neigh_meta = random_effects_meta([fnum(r["log_enrichment_ratio"]) for r in neigh_rows], [fnum(r["variance_approx"]) for r in neigh_rows])
neigh_meta_row = {
    "dataset": "GSE203612_OVCA_coordinate_neighborhood",
    "k_samples": neigh_meta.get("k", 0),
    "pooled_log_enrichment_ratio": neigh_meta.get("pooled_effect", math.nan),
    "pooled_enrichment_ratio": math.exp(neigh_meta.get("pooled_effect", math.nan)) if neigh_meta.get("k", 0) else math.nan,
    "ci_low_enrichment_ratio": math.exp(neigh_meta.get("ci_low", math.nan)) if neigh_meta.get("k", 0) else math.nan,
    "ci_high_enrichment_ratio": math.exp(neigh_meta.get("ci_high", math.nan)) if neigh_meta.get("k", 0) else math.nan,
    "p_value": neigh_meta.get("p_value", math.nan),
    "I2": neigh_meta.get("I2", math.nan),
    "method": "random-effects meta-analysis on log neighborhood enrichment ratios; variance approximation",
}
write_csv(os.path.join(DIRS["spatial"], "spatial_neighborhood_random_effects_meta.csv"), [neigh_meta_row],
          ["dataset", "k_samples", "pooled_log_enrichment_ratio", "pooled_enrichment_ratio",
           "ci_low_enrichment_ratio", "ci_high_enrichment_ratio", "p_value", "I2", "method"])

# P1-5: virtual perturbation summary.
ko_out = []
for row in virtual_ko:
    axis = row.get("axis", "")
    scenario = row.get("scenario", "")
    if axis in {"SPP1-CD44", "SPP1-ITGB1"}:
        ko_out.append({
            "sample_type": row.get("sample_type", ""),
            "axis": axis,
            "scenario": scenario,
            "control_axis_score": row.get("control_axis_score", ""),
            "perturbed_axis_score": row.get("perturbed_axis_score", ""),
            "absolute_reduction": row.get("absolute_reduction", ""),
            "relative_reduction": row.get("relative_reduction", ""),
            "conclusion_scope": "computed expression-score dependency only; not causal perturbation evidence",
        })
write_csv(os.path.join(DIRS["perturb"], "virtual_perturbation_score_dependency_summary.csv"), ko_out,
          ["sample_type", "axis", "scenario", "control_axis_score", "perturbed_axis_score",
           "absolute_reduction", "relative_reduction", "conclusion_scope"])

# P1-3: NicheNet blocker.
nichenet_rows = [{
    "required_output": "ligand activity ranking / SPP1 regulatory potential / top SPP1 target genes",
    "status": "blocked_not_found",
    "reason": "No complete NicheNet result table was found in project outputs.",
    "required_input_before_run": "patient-level stable 02/04 meta-DEG target gene set after P0 Gate A/B",
}]
write_csv(os.path.join(DIRS["nichenet"], "nichenet_execution_blocker.csv"), nichenet_rows,
          ["required_output", "status", "reason", "required_input_before_run"])

status_rows = [
    {"task": "P1-1 external scRNA projection", "status": "executed_from_existing_tables_partial",
     "main_output": "06_external_scrna_projection/external_scrna_patient_level_associations.csv",
     "limitation": "expression/signature projection; not validated CNV state transfer after patient-wise P0"},
    {"task": "P1-2 LR competition", "status": "executed_from_existing_tables_partial",
     "main_output": "07_ligand_receptor_competition/external_lr_competition_axis_ranking.csv",
     "limitation": "SPP1 axes are candidate axes; integrated clone specificity remains sample-confounded"},
    {"task": "P1-3 full NicheNet ligand-target", "status": "blocked_not_completed",
     "main_output": "08_nichenet_ligand_target/nichenet_execution_blocker.csv",
     "limitation": "requires patient-level stable 02/04 DEG from P0/P0-5"},
    {"task": "P1-4 spatial deconvolution/neighborhood", "status": "executed_sample_level_meta_partial",
     "main_output": "09_spatial_deconvolution_neighborhood/spatial_random_effects_meta_analysis.csv",
     "limitation": "uses existing scores and neighborhood table; no full deconvolution method found"},
    {"task": "P1-5 virtual perturbation", "status": "executed_score_dependency_summary",
     "main_output": "10_virtual_perturbation_causal/virtual_perturbation_score_dependency_summary.csv",
     "limitation": "score arithmetic dependency, not causal KO proof"},
]
write_csv(os.path.join(DIRS["tables"], "p1_execution_status.csv"), status_rows,
          ["task", "status", "main_output", "limitation"])

top_external = sorted(ext_axis_summary, key=lambda r: fnum(r.get("mean_axis_score")), reverse=True)[:10]
spatial_all = next((r for r in meta_rows if r["dataset"] == "ALL"), {})
spatial_gse203612 = next((r for r in meta_rows if r["dataset"] == "GSE203612_OVCA"), {})
assoc_0204 = next((r for r in external_assoc if r["stratum"] == "all_samples" and r["y"] == "Subclone02_04_common_mean_target_tumor" and r["x"] == "source_myeloid_fraction"), {})

report_md = f"""# P1 execution report from available tables

Generated: {datetime.now().isoformat(timespec="seconds")}

## Scope

This run executes P1 analyses that are possible from existing local outputs. It does not overwrite the P0 warning: integrated_oc lacks `patient_id`, and current CNV clone labels are sample/sample_type/batch confounded.

## P1 status

| Task | Status | Main output | Limitation |
|---|---|---|---|
"""
for row in status_rows:
    report_md += f"| {row['task']} | {row['status']} | `{row['main_output']}` | {row['limitation']} |\n"

report_md += f"""
## P1-1 External scRNA patient/sample-level association

External rows used:

- patient/sample summary: {len(patient_summary)}
- signature rows: {len(signature_rows)}
- axis-score rows: {len(external_axis)}

Primary association, all patient-samples:

- `source_myeloid_fraction` vs `Subclone02_04_common_mean_target_tumor`: Spearman r = {fmt(fnum(assoc_0204.get('spearman_r')))}, n = {assoc_0204.get('n_samples', 'NA')}, BH q = {fmt(fnum(assoc_0204.get('p_adj_BH')))}

Output: `06_external_scrna_projection/external_scrna_patient_level_associations.csv`

## P1-2 LR competition

Top external LR axes by mean axis score:

| sample_type | axis | n | mean score | rank |
|---|---|---:|---:|---:|
"""
for row in top_external:
    report_md += f"| {row.get('sample_type')} | {row.get('axis')} | {row.get('n_patient_samples')} | {fmt(fnum(row.get('mean_axis_score')))} | {row.get('rank_within_sample_type', '')} |\n"

report_md += f"""
Interpretation: SPP1-CD44 and SPP1-ITGB1 remain candidate axes, but the competition table also keeps MIF/APOE/TGFB/CXCL12 controls. If controls outrank SPP1 in a stratum, the SPP1 claim should be downgraded for that stratum.

## P1-3 NicheNet

Full NicheNet was not run because the required patient-level stable 02/04 DEG target set is not available after Gate A/B. A blocker table was written to:

`08_nichenet_ligand_target/nichenet_execution_blocker.csv`

## P1-4 Spatial sample-level meta-analysis

All spatial samples random-effects meta-analysis:

- k = {spatial_all.get('k_samples', 'NA')}
- pooled Spearman r = {fmt(fnum(spatial_all.get('pooled_spearman_r')))}
- 95% CI = [{fmt(fnum(spatial_all.get('ci_low_r')))}, {fmt(fnum(spatial_all.get('ci_high_r')))}]
- BH q = {fmt(fnum(spatial_all.get('p_adj_BH')))}
- I2 = {fmt(fnum(spatial_all.get('I2')))}

GSE203612 coordinate-available subset:

- k = {spatial_gse203612.get('k_samples', 'NA')}
- pooled Spearman r = {fmt(fnum(spatial_gse203612.get('pooled_spearman_r')))}
- 95% CI = [{fmt(fnum(spatial_gse203612.get('ci_low_r')))}, {fmt(fnum(spatial_gse203612.get('ci_high_r')))}]
- BH q = {fmt(fnum(spatial_gse203612.get('p_adj_BH')))}

Neighborhood enrichment meta-analysis:

- pooled enrichment ratio = {fmt(fnum(neigh_meta_row.get('pooled_enrichment_ratio')))}
- 95% CI = [{fmt(fnum(neigh_meta_row.get('ci_low_enrichment_ratio')))}, {fmt(fnum(neigh_meta_row.get('ci_high_enrichment_ratio')))}]

Output: `09_spatial_deconvolution_neighborhood/spatial_random_effects_meta_analysis.csv`

## P1-5 Virtual perturbation

Existing virtual KO tables were summarized as expression-score dependency, not causal perturbation. Output:

`10_virtual_perturbation_causal/virtual_perturbation_score_dependency_summary.csv`

## Conclusion Classification

- Partially supported: external expression-level SPP1/myeloid and 02/04-like target association; spatial co-expression gradients; candidate LR axes.
- Not supported as complete: cross-patient CNV-state mechanism, full NicheNet mechanism, causal KO directionality.
- Cannot judge from current data structure: whether integrated_oc CNV_Subclone_02/04 are true cross-patient programs.
"""

report_md_path = os.path.join(DIRS["reports"], "03_external_and_spatial_validation_report.md")
with open(report_md_path, "w", encoding="utf-8") as handle:
    handle.write(report_md)

ligand_report = report_md.replace("# P1 execution report from available tables", "# P1 ligand-target and virtual perturbation execution report")
with open(os.path.join(DIRS["reports"], "04_ligand_target_and_virtual_KO_report.md"), "w", encoding="utf-8") as handle:
    handle.write(ligand_report)

html = report_md
html = html.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
html = html.replace("\n", "<br>\n")
with open(os.path.join(DIRS["reports"], "03_external_and_spatial_validation_report.html"), "w", encoding="utf-8") as handle:
    handle.write("<!doctype html><meta charset='utf-8'><title>P1 execution report</title>"
                 "<body style='font-family:Arial,Microsoft YaHei,sans-serif;line-height:1.55;margin:32px'>"
                 f"{html}</body>")

with open(LOG_PATH, "a", encoding="utf-8") as log:
    log.write(f"patient_summary rows: {len(patient_summary)}\n")
    log.write(f"signature rows: {len(signature_rows)}\n")
    log.write(f"external axis rows: {len(external_axis)}\n")
    log.write(f"integrated LR rows: {len(integrated_lr)}\n")
    log.write(f"spatial rows: {len(spatial_corr)}\n")
    log.write(f"Finished {datetime.now().isoformat()}\n")
