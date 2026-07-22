#!/usr/bin/env python3
"""Audit GSE154763 IDs and score author-defined myeloid states on normalized expression."""
import argparse
import gzip
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd


DEFAULT_ROOT = r"D:\OC_spatiogenomics\infercnv\external_seurat_preannotation"
DEFAULT_META = r"D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154763\GSE154763_OV-FTC_metadata.csv.gz"
DEFAULT_EXPR = r"D:\OC_spatiogenomics\infercnv\external_cell_annotations\raw\GSE154763\GSE154763_OV-FTC_normalized_expression.csv.gz"


def clean_id(x, remove_suffix=False, uppercase=False):
    s = pd.Series(x, dtype="object").fillna("").astype(str).str.strip()
    s = s.str.replace(r'^["\']|["\']$', "", regex=True).str.strip()
    if remove_suffix:
        s = s.str.replace(r"-1$", "", regex=True)
    if uppercase:
        s = s.str.upper()
    return s


def audit_attempt(name, left, right):
    left = pd.Series(left, dtype="object")
    right = pd.Series(right, dtype="object")
    lset, rset = set(left), set(right)
    overlap = lset & rset
    return {
        "attempt": name,
        "metadata_rows": int(len(left)),
        "expression_rows": int(len(right)),
        "metadata_unique_ids": int(left.nunique()),
        "expression_unique_ids": int(right.nunique()),
        "metadata_duplicate_ids": int(left.duplicated().sum()),
        "expression_duplicate_ids": int(right.duplicated().sum()),
        "matched_unique_ids": int(len(overlap)),
        "metadata_match_fraction": float(len(overlap) / max(left.nunique(), 1)),
        "expression_match_fraction": float(len(overlap) / max(right.nunique(), 1)),
    }


def zscore_frame(frame):
    means = frame.mean(axis=0)
    sds = frame.std(axis=0, ddof=0).replace(0, np.nan)
    return (frame - means) / sds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=DEFAULT_ROOT)
    ap.add_argument("--metadata", default=DEFAULT_META)
    ap.add_argument("--expression", default=DEFAULT_EXPR)
    args = ap.parse_args()

    out = os.path.join(args.data_root, "diagnostics_v3_remaining_datasets", "GSE154763")
    if os.path.exists(out):
        raise RuntimeError("Output already exists; refusing to overwrite: " + out)
    os.makedirs(out)

    metadata = pd.read_csv(args.metadata, compression="infer", dtype={"index": str, "barcode": str})
    expr_ids = pd.read_csv(args.expression, compression="infer", usecols=[0], dtype=str).iloc[:, 0]
    meta_id_col = "index" if "index" in metadata.columns else metadata.columns[0]
    author_col = "MajorCluster"
    meta_ids = metadata[meta_id_col].astype(str)

    attempts = []
    candidates = []
    variants = [
        ("exact", meta_ids, expr_ids),
        ("strip_quotes_whitespace", clean_id(meta_ids), clean_id(expr_ids)),
        ("remove_terminal_dash1", clean_id(meta_ids, True), clean_id(expr_ids, True)),
        ("uppercase", clean_id(meta_ids, True, True), clean_id(expr_ids, True, True)),
        ("dot_hyphen_normalization", clean_id(meta_ids, True, True).str.replace(".", "-", regex=False),
         clean_id(expr_ids, True, True).str.replace(".", "-", regex=False)),
    ]
    for name, left, right in variants:
        row = audit_attempt(name, left, right)
        attempts.append(row)
        candidates.append((name, left, right, row))

    for prefix_name, columns in [
        ("sample_plus_barcode", ["library_id", "barcode"]),
        ("patient_plus_barcode", ["patient", "barcode"]),
    ]:
        if all(c in metadata.columns for c in columns):
            left = clean_id(metadata[columns[0]], True, True) + "__" + clean_id(metadata[columns[1]], True, True)
            right = clean_id(expr_ids, True, True)
            row = audit_attempt(prefix_name, left, right)
            attempts.append(row)
            candidates.append((prefix_name, left, right, row))
        else:
            attempts.append({"attempt": prefix_name, "status": "NOT_AVAILABLE_MISSING_METADATA_COLUMN"})

    valid = [x for x in candidates if x[3]["metadata_duplicate_ids"] == 0 and
             x[3]["expression_duplicate_ids"] == 0 and
             x[3]["metadata_match_fraction"] >= 0.99 and x[3]["expression_match_fraction"] >= 0.99]
    selected = valid[0] if valid else None
    for row in attempts:
        row["selected"] = bool(selected and row.get("attempt") == selected[0])
        if "status" not in row:
            row["status"] = "PASS_ONE_TO_ONE" if row["selected"] else "AUDITED_NOT_SELECTED"
    audit = pd.DataFrame(attempts)
    audit.to_csv(os.path.join(out, "id_matching_audit.csv"), index=False)

    if selected is None:
        summary = [
            "# GSE154763 normalized-expression reference", "",
            "- execution_status: BLOCKED_ID_MATCH",
            "- analysis_role: myeloid_reference_only",
            "- reason: no unique near-100% metadata-expression ID match was found.",
            "- normalized expression was not treated as raw counts; no count-QC, doublet detection, decontamination, raw marker testing or CNV analysis was run.",
        ]
        with open(os.path.join(out, "analysis_summary.md"), "w", encoding="utf-8") as handle:
            handle.write("\n".join(summary) + "\n")
        metadata[author_col if author_col in metadata.columns else metadata.columns[0]].value_counts().rename_axis("author_cell_type").reset_index(name="n_cells").to_csv(
            os.path.join(out, "author_metadata_celltype_summary.csv"), index=False)
        shutil.copyfile(os.path.join(out, "analysis_summary.md"), os.path.join(out, "blocked_analysis_summary.md"))
        with open(os.path.join(out, "BLOCKED_ID_MATCH"), "w", encoding="utf-8") as handle:
            handle.write("See id_match_audit.csv\n")
        return 2

    selected_name, meta_key, expr_key, selected_row = selected
    if len(meta_key) != len(metadata) or len(expr_key) != len(expr_ids):
        raise RuntimeError("ID vector length mismatch")
    metadata = metadata.copy()
    metadata["cell_id_harmonized"] = list(meta_key)
    expr_map = pd.DataFrame({"cell_id_harmonized": list(expr_key), "expr_row": np.arange(len(expr_key))})
    ordered = metadata[["cell_id_harmonized"]].merge(expr_map, on="cell_id_harmonized", how="left", validate="one_to_one")
    if ordered["expr_row"].isna().any():
        raise RuntimeError("Selected ID mapping did not cover every metadata row")

    header = pd.read_csv(args.expression, compression="infer", nrows=0).columns.tolist()
    gene_columns = header[1:]
    gene_lookup = {str(g).upper(): g for g in gene_columns}
    modules = {
        "SPP1_program": ["SPP1", "APOC1", "GPNMB", "TREM2", "LPL", "CTSD"],
        "C1QC_program": ["C1QA", "C1QB", "C1QC", "APOE", "MRC1", "SELENOP"],
        "FOLR2_program": ["FOLR2", "MRC1", "SELENOP", "C1QC", "LYVE1", "CD163"],
        "lipid_associated_program": ["TREM2", "APOC1", "GPNMB", "LPL", "CTSD", "LGALS3"],
        "GPNMB_hypoxia_program": ["GPNMB", "SPP1", "VEGFA", "BNIP3", "NDRG1", "CTSL"],
        "inflammatory_monocyte_program": ["S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "IL1B"],
    }
    required_genes = sorted(set(g for vals in modules.values() for g in vals) |
                            {"SPP1", "CD44", "ITGB1", "C1QA", "C1QB", "C1QC", "FOLR2", "TREM2"})
    present = [gene_lookup[g] for g in required_genes if g in gene_lookup]
    raw = pd.read_csv(args.expression, compression="infer", usecols=[header[0]] + present)
    raw["cell_id_harmonized"] = list(expr_key)
    raw = raw.set_index("cell_id_harmonized", verify_integrity=True)
    expr = raw.loc[metadata["cell_id_harmonized"], present].copy()
    expr.columns = [str(c).upper() for c in expr.columns]
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if expr.isna().values.all():
        raise RuntimeError("Normalized expression values could not be parsed")
    zexpr = zscore_frame(expr)

    score = pd.DataFrame(index=expr.index)
    module_gene_counts = {}
    for name, genes in modules.items():
        use = [g for g in genes if g in zexpr.columns]
        module_gene_counts[name] = len(use)
        score[name] = zexpr[use].mean(axis=1) if use else np.nan
    for g in ["SPP1", "CD44", "ITGB1", "C1QA", "C1QB", "C1QC", "FOLR2", "TREM2"]:
        if g in expr.columns:
            score[g + "_expression"] = expr[g].values
            score[g + "_positive"] = expr[g].values > 0

    def broad_type(x):
        u = str(x).upper()
        if "MACRO" in u: return "Macrophage"
        if "MONO" in u: return "Monocyte"
        if "CDC1" in u: return "cDC1"
        if "CDC2" in u: return "cDC2"
        if "CDC3" in u or "LAMP3" in u: return "cDC3_LAMP3"
        if "PDC" in u: return "pDC"
        if "MAST" in u: return "Mast"
        return "Author_myeloid_other"

    if author_col not in metadata.columns:
        raise RuntimeError("MajorCluster author annotation is missing")
    annotation = metadata.copy()
    annotation["dataset_id"] = "GSE154763"
    annotation["analysis_role"] = "myeloid_reference_only"
    annotation["cell_type_original"] = annotation[author_col].astype(str)
    annotation["cell_type_major"] = annotation[author_col].map(broad_type)
    annotation["cell_type_subtype"] = annotation[author_col].astype(str).str.replace(r"^M\d+_", "", regex=True)
    annotation["final_cell_type"] = annotation["cell_type_major"]
    annotation["cell_subtype"] = annotation["cell_type_subtype"]
    state_cols = list(modules)
    best_state = score[state_cols].idxmax(axis=1).str.replace("_program", "", regex=False)
    annotation["cell_state"] = list(best_state)
    annotation["annotation_status"] = "AUTHOR_ANNOTATION_NORMALIZED_EXPRESSION_SCORE"
    annotation["id_match_method"] = selected_name
    annotation["id_match_status"] = "PASS_ONE_TO_ONE"
    annotation["exclusion_note"] = "Not used for raw-count QC, doublet/decontamination, raw FindAllMarkers, inferCNV or tumor-wide integration"

    score_out = score.reset_index().rename(columns={"index": "cell_id_harmonized"})
    score_out.insert(0, "dataset_id", "GSE154763")
    score_out["cell_type_original"] = annotation["cell_type_original"].values
    score_out["final_cell_type"] = annotation["final_cell_type"].values
    score_out["cell_subtype"] = annotation["cell_subtype"].values
    score_out["cell_state"] = annotation["cell_state"].values
    score_out["sample_id"] = annotation["library_id"].astype(str).values if "library_id" in annotation else "NA"
    score_out["patient_id"] = annotation["patient"].astype(str).values if "patient" in annotation else "NA"

    annotation.to_csv(os.path.join(out, "author_annotation_harmonized.csv.gz"), index=False, compression="gzip")
    score_out.to_csv(os.path.join(out, "myeloid_state_scores.csv.gz"), index=False, compression="gzip")

    group_cols = ["cell_type_original", "final_cell_type", "cell_subtype"]
    joined = score_out.copy()
    summaries = []
    for keys, group in joined.groupby(group_cols, dropna=False):
        row = dict(zip(group_cols, keys))
        row["n_cells"] = len(group)
        row["n_samples"] = group["sample_id"].nunique()
        row["n_patients"] = group["patient_id"].nunique()
        for c in state_cols:
            row[c + "_mean"] = group[c].mean()
        summaries.append(row)
    subtype = pd.DataFrame(summaries).sort_values("n_cells", ascending=False)
    subtype.to_csv(os.path.join(out, "myeloid_state_by_subtype.csv"), index=False)

    spp1_rows = []
    for keys, group in joined.groupby(["cell_type_original", "sample_id"], dropna=False):
        row = {"cell_type_original": keys[0], "sample_id": keys[1], "n_cells": len(group),
               "SPP1_program_mean": group["SPP1_program"].mean()}
        if "SPP1_expression" in group:
            row["SPP1_average_expression"] = group["SPP1_expression"].mean()
            row["SPP1_positive_fraction"] = group["SPP1_positive"].mean()
        spp1_rows.append(row)
    spp1 = pd.DataFrame(spp1_rows)
    spp1.to_csv(os.path.join(out, "SPP1_reference_summary.csv"), index=False)

    counts = annotation["final_cell_type"].value_counts()
    summary = [
        "# GSE154763 normalized-expression myeloid reference", "",
        "- execution_status: SUCCESS_NORMALIZED_EXPRESSION_REFERENCE",
        "- analysis_role: myeloid_reference_only",
        "- id_match_status: PASS_ONE_TO_ONE",
        "- id_match_method: " + selected_name,
        "- matched_cells: " + str(len(annotation)),
        "- normalized_expression_handling: author-annotation-driven module scoring only",
        "- prohibited_operations_not_run: scDblFinder; SoupX; DecontX; count-based QC; raw FindAllMarkers; inferCNV",
        "- module_genes_detected: " + json.dumps(module_gene_counts, sort_keys=True), "",
        "## Author-driven broad myeloid counts", "",
    ] + ["- {}: {}".format(k, int(v)) for k, v in counts.items()]
    with open(os.path.join(out, "analysis_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
