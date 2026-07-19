import csv
import math
import os
from collections import defaultdict
from datetime import datetime


ROOT = r"D:\OC_spatiogenomics\infercnv"
OUT = os.path.join(ROOT, "06_rigorous_validation")
SRC = os.path.join(ROOT, "05_sample_type_lr_niche_external_validation", "sample_type_LR_niche_analysis")
TABLES = os.path.join(OUT, "tables")
REPORTS = os.path.join(OUT, "reports")
LOGS = os.path.join(OUT, "logs")
P1_DIRS = {
    "external": os.path.join(OUT, "06_external_scrna_projection"),
    "lr": os.path.join(OUT, "07_ligand_receptor_competition"),
    "nichenet": os.path.join(OUT, "08_nichenet_ligand_target"),
    "spatial": os.path.join(OUT, "09_spatial_deconvolution_neighborhood"),
    "perturb": os.path.join(OUT, "10_virtual_perturbation_causal"),
}

for path in [TABLES, REPORTS, LOGS, *P1_DIRS.values()]:
    os.makedirs(path, exist_ok=True)


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def to_float(value):
    try:
        return float(value)
    except Exception:
        return math.nan


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


log_path = os.path.join(LOGS, f"10_p1_existing_evidence_audit_{datetime.now():%Y%m%d_%H%M%S}.log")
with open(log_path, "w", encoding="utf-8") as log:
    log.write(f"Started {datetime.now().isoformat()}\n")
    log.write(f"Source: {SRC}\n")

external_path = os.path.join(SRC, "tables", "external_scRNA_all_datasets_sampletype_axis_scores.csv")
spatial_path = os.path.join(SRC, "tables", "spatial_correlation_SPP1_myeloid_Target_axis.csv")
neigh_path = os.path.join(SRC, "tables", "spatial_neighborhood_enrichment.csv")
lr_path = os.path.join(SRC, "tables", "sample_id_LR_opportunity_scores_all_axes.csv")
ko_path = os.path.join(SRC, "tables", "virtual_KO_LR_score_reduction_summary_by_sample_type.csv")
limitations_path = os.path.join(SRC, "tables", "limitations_summary.csv")

external = read_csv(external_path) if os.path.exists(external_path) else []
spatial = read_csv(spatial_path) if os.path.exists(spatial_path) else []
neigh = read_csv(neigh_path) if os.path.exists(neigh_path) else []
lr = read_csv(lr_path) if os.path.exists(lr_path) else []
ko = read_csv(ko_path) if os.path.exists(ko_path) else []
limitations = read_csv(limitations_path) if os.path.exists(limitations_path) else []

external_summary = []
by_axis_type = defaultdict(list)
for row in external:
    by_axis_type[(row.get("dataset", ""), row.get("sample_type", ""), row.get("axis", ""))].append(to_float(row.get("axis_score", "")))
for (dataset, sample_type, axis), values in sorted(by_axis_type.items()):
    values = [v for v in values if not math.isnan(v)]
    if not values:
        continue
    external_summary.append({
        "dataset": dataset,
        "sample_type": sample_type,
        "axis": axis,
        "n_patient_samples": len(values),
        "mean_axis_score": sum(values) / len(values),
        "median_axis_score": sorted(values)[len(values) // 2],
    })
write_csv(
    os.path.join(P1_DIRS["external"], "external_scrna_axis_summary_from_existing_outputs.csv"),
    external_summary,
    ["dataset", "sample_type", "axis", "n_patient_samples", "mean_axis_score", "median_axis_score"],
)

spatial_summary = []
for row in spatial:
    r = to_float(row.get("spearman_r", ""))
    spatial_summary.append({
        "dataset": row.get("dataset", ""),
        "sample_id": row.get("sample_id", ""),
        "n_spots": row.get("n_spots", ""),
        "spearman_r": r,
        "spearman_p": row.get("spearman_p", ""),
        "direction": "positive" if r > 0 else "negative" if r < 0 else "zero",
        "has_coordinates_for_neighborhood": "yes" if row.get("dataset", "") == "GSE203612_OVCA" else "not_in_existing_neighborhood_table",
    })
write_csv(
    os.path.join(P1_DIRS["spatial"], "spatial_correlation_summary_from_existing_outputs.csv"),
    spatial_summary,
    ["dataset", "sample_id", "n_spots", "spearman_r", "spearman_p", "direction", "has_coordinates_for_neighborhood"],
)

neigh_summary = []
for row in neigh:
    neigh_summary.append({
        "dataset": row.get("dataset", ""),
        "sample_id": row.get("sample_id", ""),
        "source_class": row.get("source_class", ""),
        "target_class": row.get("target_class", ""),
        "enrichment_ratio": row.get("enrichment_ratio", ""),
        "n_source_spots": row.get("n_source_spots", ""),
        "n_target_spots": row.get("n_target_spots", ""),
    })
write_csv(
    os.path.join(P1_DIRS["spatial"], "spatial_neighborhood_summary_from_existing_outputs.csv"),
    neigh_summary,
    ["dataset", "sample_id", "source_class", "target_class", "enrichment_ratio", "n_source_spots", "n_target_spots"],
)

axes_of_interest = {"SPP1-CD44", "SPP1-ITGB1", "MIF-CD74", "APOE-LRP1", "TIMP1-CD63", "TGFB1-TGFBR1", "TGFB1-TGFBR2", "CXCL12-CXCR4"}
lr_summary = []
lr_grouped = defaultdict(list)
for row in lr:
    axis = row.get("axis", "")
    if axis in axes_of_interest and row.get("target_clone", "") in {"Subclone_02", "Subclone_04"}:
        lr_grouped[(row.get("level", ""), row.get("group", ""), row.get("target_clone", ""), axis)].append(to_float(row.get("axis_score", "")))
for (level, group, target_clone, axis), values in sorted(lr_grouped.items()):
    values = [v for v in values if not math.isnan(v)]
    if values:
        lr_summary.append({
            "level": level,
            "group": group,
            "target_clone": target_clone,
            "axis": axis,
            "n_rows": len(values),
            "mean_axis_score": sum(values) / len(values),
            "max_axis_score": max(values),
        })
write_csv(
    os.path.join(P1_DIRS["lr"], "lr_competition_summary_from_existing_outputs.csv"),
    lr_summary,
    ["level", "group", "target_clone", "axis", "n_rows", "mean_axis_score", "max_axis_score"],
)

ko_summary = []
for row in ko:
    if row.get("axis", "") in {"SPP1-CD44", "SPP1-ITGB1"}:
        ko_summary.append({
            "sample_type": row.get("sample_type", ""),
            "axis": row.get("axis", ""),
            "scenario": row.get("scenario", ""),
            "control_axis_score": row.get("control_axis_score", ""),
            "perturbed_axis_score": row.get("perturbed_axis_score", ""),
            "relative_reduction": row.get("relative_reduction", ""),
            "interpretation": "expression-score arithmetic perturbation, not causal KO evidence",
        })
write_csv(
    os.path.join(P1_DIRS["perturb"], "virtual_ko_summary_from_existing_outputs.csv"),
    ko_summary,
    ["sample_type", "axis", "scenario", "control_axis_score", "perturbed_axis_score", "relative_reduction", "interpretation"],
)

status_rows = [
    {
        "task": "P1-1 external scRNA projection",
        "status": "partially_supported_existing_output",
        "reason": "External Zhang2022 patient_id/sample_id tables exist, but they are expression/signature projections rather than validated CNV program transfer after Gate A.",
    },
    {
        "task": "P1-2 LR competition",
        "status": "partially_supported_existing_output",
        "reason": "Multiple LR axes are available, but integrated_oc target clone interpretation remains sample_type/sample_id confounded.",
    },
    {
        "task": "P1-3 full NicheNet ligand-target",
        "status": "not_completed",
        "reason": "No complete NicheNet ligand activity ranking, SPP1 regulatory potential, or SPP1 target gene table was found.",
    },
    {
        "task": "P1-4 spatial deconvolution/neighborhood",
        "status": "partial_spatial_expression_neighborhood_only",
        "reason": "Existing outputs include spot-level correlations and GSE203612 neighborhood enrichment; full deconvolution plus sample-level random-effects meta-analysis was not found.",
    },
    {
        "task": "P1-5 virtual perturbation",
        "status": "exploratory_score_perturbation_only",
        "reason": "Existing virtual KO outputs are expression-score reductions, not CellOracle/scTenifold/NicheNet causal perturbation.",
    },
]
write_csv(os.path.join(TABLES, "p1_task_status_audit.csv"), status_rows, ["task", "status", "reason"])

positive_spatial = sum(1 for row in spatial_summary if row["direction"] == "positive")
negative_spatial = sum(1 for row in spatial_summary if row["direction"] == "negative")

report = f"""# P1 evidence audit under rigorous validation

Generated: {datetime.now().isoformat(timespec="seconds")}

## Overall decision

P1 cannot be marked complete as mechanistic validation, because P0/Gate A remains blocked by missing `patient_id` in the integrated_oc metadata and the current clone labels are strongly sample/sample_type/batch confounded.

Existing P1-like outputs were found and audited. They can support a cautious candidate niche model, but not a definitive SPP1-driven, cross-patient CNV state mechanism.

## P1 task status

| Task | Status | Reason |
|---|---|---|
"""
for row in status_rows:
    report += f"| {row['task']} | {row['status']} | {row['reason']} |\n"

report += f"""
## External scRNA

External scRNA tables contain {len(external)} axis-score rows and include explicit `patient_id`, `sample_id`, `sample_type`, anatomical location and treatment phase fields. The summarized output is:

`06_external_scrna_projection/external_scrna_axis_summary_from_existing_outputs.csv`

Interpretation: useful for candidate expression-level replication, but not a replacement for patient-wise CNV reconstruction or validated 02/04-like state transfer.

## LR competition

LR competition summaries were generated for SPP1-CD44, SPP1-ITGB1, MIF-CD74, APOE-LRP1, TIMP1-CD63, TGFB1-TGFBR1/2 and CXCL12-CXCR4 where present.

`07_ligand_receptor_competition/lr_competition_summary_from_existing_outputs.csv`

Interpretation: SPP1-CD44 and SPP1-ITGB1 remain candidate axes. Because integrated_oc sample_id, sample_type and batch are overlapping, this cannot establish clone-specific mechanism.

## NicheNet

No complete NicheNet ligand-target output was found. This means P1-3 is not completed. Checking package installation or gene universe is not enough for this task; it requires ligand activity ranking, SPP1 regulatory potential and target-gene enrichment based on patient-level stable DEG.

## Spatial evidence

Spatial correlation rows audited: {len(spatial)}. Direction count: {positive_spatial} positive, {negative_spatial} negative.

Neighborhood enrichment rows audited: {len(neigh)}.

Outputs:

- `09_spatial_deconvolution_neighborhood/spatial_correlation_summary_from_existing_outputs.csv`
- `09_spatial_deconvolution_neighborhood/spatial_neighborhood_summary_from_existing_outputs.csv`

Interpretation: GSE203612 has coordinate-based neighborhood evidence, while other local spatial tables are mainly spot-level score correlations. This remains partial because full deconvolution and sample-level random-effects meta-analysis were not found.

## Virtual perturbation

Virtual KO summaries were found and copied into:

`10_virtual_perturbation_causal/virtual_ko_summary_from_existing_outputs.csv`

Interpretation: these are arithmetic expression-score perturbations and must be labelled computational predictions, not causal proof.

## Required before true P1 completion

1. Add explicit patient_id mapping for integrated_oc samples.
2. Finish P0 patient-wise CNV reconstruction/stability/QC deconfounding.
3. Define patient-level stable 02/04 DEG or validated 02/04-like target gene set.
4. Run full NicheNet using SPP1+ myeloid sender and validated 02/04-like tumor receiver.
5. Run at least one spatial deconvolution method plus a sensitivity method, then aggregate per section/sample.
"""

with open(os.path.join(REPORTS, "03_external_and_spatial_validation_report.md"), "w", encoding="utf-8") as handle:
    handle.write(report)

with open(os.path.join(REPORTS, "04_ligand_target_and_virtual_KO_report.md"), "w", encoding="utf-8") as handle:
    handle.write(report.replace("# P1 evidence audit under rigorous validation", "# Ligand-target and virtual KO audit"))

with open(log_path, "a", encoding="utf-8") as log:
    log.write(f"External rows: {len(external)}\n")
    log.write(f"Spatial rows: {len(spatial)}\n")
    log.write(f"LR rows: {len(lr)}\n")
    log.write(f"KO rows: {len(ko)}\n")
    log.write(f"Finished {datetime.now().isoformat()}\n")
