import csv
import os
from datetime import datetime


ROOT = r"D:\OC_spatiogenomics\infercnv\06_rigorous_validation"
REPORTS = os.path.join(ROOT, "reports")
TABLES = os.path.join(ROOT, "tables")
LOGS = os.path.join(ROOT, "logs")
os.makedirs(REPORTS, exist_ok=True)
os.makedirs(TABLES, exist_ok=True)
os.makedirs(LOGS, exist_ok=True)


def exists(rel):
    return os.path.exists(os.path.join(ROOT, rel.replace("/", os.sep)))


def read_text(rel):
    path = os.path.join(ROOT, rel.replace("/", os.sep))
    if not os.path.exists(path):
        return ""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def write_text(rel, text):
    path = os.path.join(ROOT, rel.replace("/", os.sep))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def md_to_html(md, title):
    body = (
        md.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace("\n", "<br>\n")
    )
    return (
        "<!doctype html><meta charset='utf-8'>"
        f"<title>{title}</title>"
        "<body style='font-family:Arial,Microsoft YaHei,sans-serif;line-height:1.55;margin:32px'>"
        f"{body}</body>"
    )


generated = datetime.now().isoformat(timespec="seconds")

functional_md = f"""# Patient-level functional validation report

Generated: {generated}

## Status

P0-5 / Gate B is not fully completed because the required `patient_id` field is missing from the integrated_oc metadata and current CNV clone labels failed preliminary sample-confounding checks.

## What is available

- Existing sample-level CNV-expression outputs are present in `CNV_expression_joint_analysis_integrated_oc`.
- Existing functional score summaries and pseudobulk-style tables can be used only as exploratory sample-level evidence.
- These results must not be reported as patient-level biological replication or as primary P values from independent patients.

## Why this is not marked complete

The task book requires patient/sample-level pseudobulk and leave-one-patient-out/meta-analysis. Because `patient_id` is absent and only 4 sample IDs are available in the integrated metadata, Gate B cannot be judged as a validated patient-level 02/04 functional state.

## Required to complete

1. Provide a reviewed sample-to-patient mapping.
2. Rebuild patient-wise CNV programs or verify existing clone labels after Gate A.
3. Re-run pseudobulk DEG/GSEA/PROGENy/TF activity with patient/sample as the replicate unit.
4. Run leave-one-patient-out sensitivity analysis.

## Conclusion

Classification: 因数据结构无法判断 / blocked by missing patient-level replication.
"""
write_text("reports/02_patient_level_functional_report.md", functional_md)
write_text("reports/02_patient_level_functional_report.html", md_to_html(functional_md, "Patient-level functional validation report"))

ligand_md = read_text("reports/04_ligand_target_and_virtual_KO_report.md")
if ligand_md:
    write_text("reports/04_ligand_target_and_virtual_KO_report.html", md_to_html(ligand_md, "Ligand-target and virtual KO report"))

items = [
    {
        "section": "Directory/config",
        "item": "Create 06_rigorous_validation directory structure and project_config.yaml",
        "status": "completed",
        "evidence": "00_config/project_config.yaml; required directories created",
        "limitation": "",
    },
    {
        "section": "P0-1",
        "item": "Clone x patient/sample confounding diagnosis",
        "status": "completed_with_patient_blocker",
        "evidence": "01_clone_patient_confounding/*.csv; reports/01_clone_validity_report.html",
        "limitation": "patient_id missing; patient-level diagnosis not evaluable",
    },
    {
        "section": "P0-2",
        "item": "Patient-wise CNV reconstruction",
        "status": "blocked_not_completed",
        "evidence": "reports/01_clone_validity_report.html",
        "limitation": "requires explicit patient_id and patient-wise malignant/reference cell grouping",
    },
    {
        "section": "P0-3",
        "item": "Clone stability and multi-algorithm consistency",
        "status": "blocked_not_completed",
        "evidence": "No stability/ARI/bootstrap outputs in 03_clone_stability_consensus",
        "limitation": "should follow patient-wise CNV reconstruction; current clone labels sample-confounded",
    },
    {
        "section": "P0-4",
        "item": "QC, contamination and transcription-state deconfounding",
        "status": "partial_field_check_only",
        "evidence": "00_config/metadata_field_mapping_report.csv",
        "limitation": "doublet_score, S.Score, G2M.Score and patient_id missing; no mixed model deconfounding completed",
    },
    {
        "section": "P0-5",
        "item": "Patient/sample-level functional validation",
        "status": "blocked_exploratory_only",
        "evidence": "reports/02_patient_level_functional_report.html",
        "limitation": "patient-level pseudobulk and leave-one-patient-out cannot be validated without patient_id",
    },
    {
        "section": "P1-1",
        "item": "External scRNA projection of 02/04-like state",
        "status": "partial_completed_from_existing_tables",
        "evidence": "06_external_scrna_projection/external_scrna_patient_level_associations.csv",
        "limitation": "expression/signature projection; not validated CNV program transfer",
    },
    {
        "section": "P1-2",
        "item": "LR competition axis analysis",
        "status": "partial_completed_from_existing_tables",
        "evidence": "07_ligand_receptor_competition/external_lr_competition_axis_ranking.csv",
        "limitation": "SPP1 axis is candidate only; integrated clone specificity sample-confounded",
    },
    {
        "section": "P1-3",
        "item": "Full NicheNet ligand-target analysis",
        "status": "blocked_not_completed",
        "evidence": "08_nichenet_ligand_target/nichenet_execution_blocker.csv",
        "limitation": "requires patient-level stable 02/04 meta-DEG target set",
    },
    {
        "section": "P1-4",
        "item": "Spatial deconvolution and neighborhood analysis",
        "status": "partial_completed",
        "evidence": "09_spatial_deconvolution_neighborhood/spatial_random_effects_meta_analysis.csv",
        "limitation": "sample-level meta and existing neighborhood only; no full deconvolution method completed",
    },
    {
        "section": "P1-5",
        "item": "Virtual perturbation and directionality",
        "status": "partial_score_dependency_only",
        "evidence": "10_virtual_perturbation_causal/virtual_perturbation_score_dependency_summary.csv",
        "limitation": "score arithmetic; not causal KO evidence",
    },
    {
        "section": "P2-1",
        "item": "Bulk result repositioning",
        "status": "completed",
        "evidence": "11_bulk_negative_validation/bulk_survival_interaction_negative_results_retained.csv",
        "limitation": "OS/PFS retained as negative/not primary",
    },
    {
        "section": "P2-2",
        "item": "CNV-expression dosage independent validation",
        "status": "partial_relabel_completed",
        "evidence": "12_integrated_evidence_scoring/RNA_derived_CNV_expression_coupling_candidates.csv",
        "limitation": "no independent DNA CNV + RNA validation found locally",
    },
    {
        "section": "P2-3",
        "item": "Statistical mediation / SEM exploration",
        "status": "partial_bulk_exploratory_completed",
        "evidence": "12_integrated_evidence_scoring/bulk_exploratory_mediation_bootstrap.csv",
        "limitation": "bulk statistical mediation only; not causal proof",
    },
    {
        "section": "Final reports",
        "item": "Required final report file set",
        "status": "completed_as_reports_with_limitations",
        "evidence": "reports/01..05 plus FINAL_CODEX_EXECUTION_SUMMARY.md",
        "limitation": "reports include blocked/partial classifications where data requirements are unmet",
    },
]

fields = ["section", "item", "status", "evidence", "limitation"]
with open(os.path.join(TABLES, "taskbook_completion_audit.csv"), "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(items)

summary_lines = [
    "# Taskbook completion audit",
    "",
    f"Generated: {generated}",
    "",
    "## Executive status",
    "",
    "All taskbook parts have now been either executed where possible or explicitly accounted for with a blocker/limitation. The full taskbook is not scientifically complete because key prerequisites are missing, most importantly `patient_id` for integrated_oc and patient-wise CNV reconstruction.",
    "",
    "## Audit table",
    "",
    "| Section | Item | Status | Evidence | Limitation |",
    "|---|---|---|---|---|",
]
for row in items:
    summary_lines.append(
        f"| {row['section']} | {row['item']} | {row['status']} | `{row['evidence']}` | {row['limitation']} |"
    )

summary_lines.extend([
    "",
    "## Final conclusion classes",
    "",
    "- Supported: directory/config setup, P0-1 sample-level confounding diagnosis, P2 bulk repositioning.",
    "- Partially supported: P1 external/LR/spatial/virtual-perturbation evidence, P2 bulk covariation, P2 exploratory mediation, RNA-derived CNV-expression coupling candidates.",
    "- Not supported: direct external SPP1-myeloid to 02/04-like association, OS/PFS interaction, complete NicheNet mechanism.",
    "- Not judgeable from current data structure: cross-patient CNV clone validity and patient-level 02/04 functional state, due to missing patient_id.",
])

audit_md = "\n".join(summary_lines) + "\n"
write_text("reports/TASKBOOK_COMPLETION_AUDIT.md", audit_md)
write_text("reports/TASKBOOK_COMPLETION_AUDIT.html", md_to_html(audit_md, "Taskbook completion audit"))

final_required = [
    "reports/01_clone_validity_report.html",
    "reports/02_patient_level_functional_report.html",
    "reports/03_external_and_spatial_validation_report.html",
    "reports/04_ligand_target_and_virtual_KO_report.html",
    "reports/05_integrated_evidence_and_limitations.docx",
    "reports/FINAL_CODEX_EXECUTION_SUMMARY.md",
]
final_rows = []
for rel in final_required:
    final_rows.append({
        "required_file": rel,
        "exists": exists(rel),
        "status": "present" if exists(rel) else "missing",
    })
write_csv_path = os.path.join(TABLES, "final_required_outputs_audit.csv")
with open(write_csv_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["required_file", "exists", "status"])
    writer.writeheader()
    writer.writerows(final_rows)

with open(os.path.join(LOGS, f"20_generate_taskbook_completion_audit_{datetime.now():%Y%m%d_%H%M%S}.log"), "w", encoding="utf-8") as handle:
    handle.write(f"Generated completion audit at {generated}\n")
    handle.write(f"Rows: {len(items)}\n")
