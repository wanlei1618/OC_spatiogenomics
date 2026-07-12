from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = WORKFLOW_ROOT / "config"
TABLES_DIR = WORKFLOW_ROOT / "tables"
REPORTS_DIR = WORKFLOW_ROOT / "reports"
ANNOTATIONS_DIR = WORKFLOW_ROOT / "annotations"
DEFAULT_DATA_ROOT = Path(r"D:\OC_spatiogenomics\infercnv\external_cell_annotations")

STANDARD_COLUMNS = [
    "dataset_id",
    "sample_id",
    "patient_id",
    "cell_id_original",
    "barcode_raw",
    "barcode_normalized",
    "cell_key",
    "cell_type_original",
    "cell_type_major",
    "cell_type_subtype",
    "annotation_source",
    "annotation_level",
    "disease",
    "tissue",
    "anatomical_site",
    "treatment",
    "timepoint",
    "malignant_status",
    "analysis_role",
    "source_file",
    "source_url",
    "match_method",
    "match_confidence",
]


DOWNLOADS = [
    {
        "dataset_id": "GSE154763",
        "file_id": "GSE154763_OV-FTC_metadata.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154763/suppl/GSE154763_OV-FTC_metadata.csv.gz",
        "required": True,
    },
    {
        "dataset_id": "GSE154763",
        "file_id": "GSE154763_OV-FTC_normalized_expression.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154763/suppl/GSE154763_OV-FTC_normalized_expression.csv.gz",
        "required": True,
    },
    {
        "dataset_id": "GSE158722",
        "file_id": "GSE158722.cell_annotations.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158722/suppl/GSE158722.cell_annotations.txt.gz",
        "required": True,
    },
]

for tumor, tag in {
    "59": "cu0pxp0sepsvhos",
    "76": "taff3g6nlmf5yk4",
    "77": "2otdxfunxnlr6ql",
    "89": "sojupsrjxnvnf6d",
    "90": "ydokzcxgkugo678",
}.items():
    DOWNLOADS.append(
        {
            "dataset_id": "GSE154600",
            "file_id": f"sample{tumor}_sce.rds",
            "url": f"https://dl.dropboxusercontent.com/s/{tag}/sample{tumor}_sce.rds",
            "required": True,
        }
    )

OPTIONAL_DOWNLOADS = [
    {
        "dataset_id": "GSE154600",
        "file_id": "GSE154600_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154600/suppl/GSE154600_RAW.tar",
        "required": False,
    },
    {
        "dataset_id": "GSE147082",
        "file_id": "GSE147082_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147082/suppl/GSE147082_RAW.tar",
        "required": False,
    },
    {
        "dataset_id": "GSE151214",
        "file_id": "GSE151214_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151214/suppl/GSE151214_RAW.tar",
        "required": False,
    },
]


def ensure_layout(data_root: Path) -> None:
    for rel in [
        "raw/GSE147082",
        "raw/GSE151214",
        "raw/GSE154600",
        "raw/GSE154763",
        "raw/GSE158722",
        "processed",
        "audit",
        "logs",
        "cache",
        "tmp",
        "R_user",
        "R_library",
    ]:
        (data_root / rel).mkdir(parents=True, exist_ok=True)
    for rel in [
        "annotations/manifests",
        "annotations/harmonized",
        "tables/archive_marker_score_v1",
        "reports",
        "figures",
    ]:
        (WORKFLOW_ROOT / rel).mkdir(parents=True, exist_ok=True)


def progress(data_root: Path, message: str) -> None:
    ensure_layout(data_root)
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    path = data_root / "logs" / "codex_progress.md"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"- {stamp} {message}\n")
    print(message, flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", errors="replace") if "t" in mode else gzip.open(path, mode)
    return path.open(mode, encoding="utf-8", errors="replace") if "t" in mode else path.open(mode)


def normalize_barcode(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = text.replace('"', "").replace("'", "")
    text = re.sub(r"\s+", "", text)
    text = re.sub(r"-1$", "", text)
    return text.upper()


def make_cell_key(row: pd.Series) -> str:
    parts = [row.get("patient_id", ""), row.get("sample_id", ""), row.get("barcode_normalized", "")]
    return "|".join(str(x) for x in parts if str(x) and str(x).lower() != "nan")


def registry() -> pd.DataFrame:
    return pd.read_csv(CONFIG_DIR / "external_scrna_dataset_registry.csv")


_HARMONIZATION_CACHE = None


def harmonization_map() -> pd.DataFrame:
    global _HARMONIZATION_CACHE
    if _HARMONIZATION_CACHE is not None:
        return _HARMONIZATION_CACHE
    _HARMONIZATION_CACHE = pd.read_csv(CONFIG_DIR / "celltype_harmonization_map.csv")
    return _HARMONIZATION_CACHE


def harmonize_label(label: object) -> tuple[str, str, str]:
    text = "" if pd.isna(label) else str(label)
    low = text.lower()
    hmap = harmonization_map()
    for row in hmap.itertuples(index=False):
        pat = str(row.pattern).lower()
        if pat and pat in low:
            return row.cell_type_major, row.cell_type_subtype, row.malignant_status
    return "Unknown", text if text else "Unknown", "unknown"


def write_csv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, compression="gzip")


def read_harmonized(dataset_id: str) -> pd.DataFrame:
    paths = {
        "GSE154600": WORKFLOW_ROOT / "annotations/harmonized/GSE154600_original_annotations.csv.gz",
        "GSE154763": WORKFLOW_ROOT / "annotations/harmonized/GSE154763_OVFTC_original_annotations.csv.gz",
        "GSE158722": WORKFLOW_ROOT / "annotations/harmonized/GSE158722_original_annotations.csv.gz",
        "GSE147082": WORKFLOW_ROOT / "annotations/harmonized/GSE147082_secondary_annotations.csv.gz",
        "GSE151214": WORKFLOW_ROOT / "annotations/harmonized/GSE151214_normal_reference_annotations.csv.gz",
    }
    path = paths[dataset_id]
    if not path.exists():
        return pd.DataFrame(columns=STANDARD_COLUMNS)
    return pd.read_csv(path)


def download(data_root: Path, include_optional_raw: bool = False, force: bool = False, dry_run: bool = False) -> pd.DataFrame:
    ensure_layout(data_root)
    rows = []
    downloads = list(DOWNLOADS) + (OPTIONAL_DOWNLOADS if include_optional_raw else [])
    for rec in downloads:
        out = data_root / "raw" / rec["dataset_id"] / rec["file_id"]
        row = dict(rec)
        row["path"] = str(out)
        if dry_run:
            row.update(status="dry_run", bytes=0, sha256="")
            rows.append(row)
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        part = out.with_suffix(out.suffix + ".part")
        try:
            if out.exists() and out.stat().st_size > 0 and not force:
                status = "cached"
            else:
                if part.exists():
                    part.unlink()
                req = urllib.request.Request(rec["url"], headers={"User-Agent": "OC-spatiogenomics-Codex"})
                with urllib.request.urlopen(req, timeout=180) as src, part.open("wb") as dst:
                    shutil.copyfileobj(src, dst, length=1024 * 1024)
                if part.stat().st_size == 0:
                    raise RuntimeError("downloaded file is empty")
                os.replace(part, out)
                status = "downloaded"
            row.update(status=status, bytes=out.stat().st_size, sha256=sha256_file(out), message="ok")
        except Exception as exc:
            if part.exists():
                part.unlink()
            row.update(status="failed", bytes=0, sha256="", message=str(exc))
        rows.append(row)
    df = pd.DataFrame(rows)
    manifest = WORKFLOW_ROOT / "annotations/manifests/original_annotation_download_manifest.csv"
    sha = WORKFLOW_ROOT / "annotations/manifests/original_annotation_sha256.csv"
    df.to_csv(manifest, index=False)
    df[["dataset_id", "file_id", "path", "bytes", "sha256", "status"]].to_csv(sha, index=False)
    return df


def archive_marker_results() -> None:
    archive = TABLES_DIR / "archive_marker_score_v1"
    archive.mkdir(parents=True, exist_ok=True)
    for name in [
        "extended_external_scRNA_dataset_level_summary.csv",
        "extended_external_scRNA_meta_summary.csv",
        "extended_external_scRNA_celltype_composition.csv",
        "extended_external_scRNA_status.csv",
    ]:
        src = TABLES_DIR / name
        if src.exists():
            shutil.copy2(src, archive / name)


def copy_task_registry() -> None:
    src = Path(r"D:\Downloads\OC_external_scrna_annotation_codex\external_scrna_dataset_registry.csv")
    if src.exists():
        shutil.copy2(src, CONFIG_DIR / "external_scrna_dataset_registry_task_input.csv")


def infer_column(columns: list[str], patterns: list[str]) -> str | None:
    lowered = {c.lower(): c for c in columns}
    for pat in patterns:
        for low, col in lowered.items():
            if pat in low:
                return col
    return None


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def repo_h5_path(existing_repo_id: str) -> Path:
    candidates = [
        REPO_ROOT / "data" / "external_validation" / f"{existing_repo_id}_scRNA" / f"{existing_repo_id}_expression.h5",
        Path(r"D:\OC_spatiogenomics\公开集\单细胞") / f"{existing_repo_id}_expression.h5",
        Path(r"D:\OC_spatiogenomics\公开集\Nath et al") / f"{existing_repo_id}_expression.h5",
        Path(r"D:\OC_spatiogenomics\github_upload\OC_spatiogenomics_from_c\data\external_validation") / f"{existing_repo_id}_scRNA" / f"{existing_repo_id}_expression.h5",
    ]
    return first_existing(candidates) or candidates[0]


def h5_cell_count(path: Path) -> tuple[int, list[str]]:
    if not path.exists():
        return 0, []
    try:
        import h5py
        with h5py.File(path, "r") as handle:
            if "matrix" in handle and "barcodes" in handle["matrix"]:
                barcodes = [x.decode("utf-8", "ignore") if isinstance(x, bytes) else str(x) for x in handle["matrix"]["barcodes"][:]]
                return len(barcodes), barcodes
    except Exception:
        return 0, []
    return 0, []


def build_status_and_reports(data_root: Path) -> None:
    reg = registry()
    status_rows = []
    for _, row in reg.iterrows():
        ann = read_harmonized(row["dataset_id"])
        status_rows.append(
            {
                "dataset_id": row["dataset_id"],
                "analysis_role": row["analysis_role"],
                "annotation_status": row["annotation_status"],
                "n_annotation_cells": len(ann),
                "primary_tumor_lr": row["primary_tumor_lr"],
                "sensitivity_tumor_lr": row["sensitivity_tumor_lr"],
                "myeloid_source_reference": row["myeloid_source_reference"],
                "normal_reference": row["normal_reference"],
                "status": "available" if len(ann) else "missing_or_blocked",
            }
        )
    pd.DataFrame(status_rows).to_csv(TABLES_DIR / "external_scrna_annotation_status.csv", index=False)


def standardize_rows(df: pd.DataFrame, dataset_id: str, source_file: str, source_url: str) -> pd.DataFrame:
    reg = registry().set_index("dataset_id").loc[dataset_id]
    for col in STANDARD_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df["dataset_id"] = dataset_id
    df["annotation_source"] = reg["annotation_status"]
    df["analysis_role"] = reg["analysis_role"]
    df["source_file"] = source_file
    df["source_url"] = source_url
    if "barcode_raw" not in df or df["barcode_raw"].eq("").all():
        df["barcode_raw"] = df["cell_id_original"]
    df["barcode_normalized"] = df["barcode_raw"].map(normalize_barcode)
    bad_norm = df["barcode_normalized"].eq("") | df["barcode_normalized"].duplicated(keep=False)
    if bad_norm.any():
        df.loc[bad_norm, "barcode_normalized"] = [
            f"{dataset_id}_{i + 1}" for i in df.index[bad_norm]
        ]
    df["cell_key"] = df.apply(make_cell_key, axis=1)
    harmonized = df["cell_type_original"].map(harmonize_label)
    df["cell_type_major"] = [x[0] for x in harmonized]
    df["cell_type_subtype"] = [x[1] for x in harmonized]
    df["malignant_status"] = [x[2] for x in harmonized]
    df["match_method"] = df["match_method"].replace("", "not_matched_yet")
    df["match_confidence"] = df["match_confidence"].replace("", "not_matched_yet")
    return df[STANDARD_COLUMNS]


def safe_read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.name.endswith(".txt") or path.name.endswith(".txt.gz") else ","
    return pd.read_csv(path, sep=sep, compression="gzip" if path.suffix == ".gz" else None, low_memory=False)


def extract_non_sce_annotations(data_root: Path) -> None:
    # GSE154763 author metadata.
    meta_path = data_root / "raw/GSE154763/GSE154763_OV-FTC_metadata.csv.gz"
    if meta_path.exists():
        meta = safe_read_table(meta_path)
        cols = list(meta.columns)
        cell_col = infer_column(cols, ["cell", "barcode"]) or cols[0]
        sample_col = infer_column(cols, ["sample", "orig.ident", "gsm"]) or ""
        type_col = infer_column(cols, ["celltype", "cell_type", "cluster", "annotation", "type"]) or cols[-1]
        out = pd.DataFrame(
            {
                "sample_id": meta[sample_col] if sample_col else "GSE154763",
                "patient_id": meta[sample_col] if sample_col else "GSE154763",
                "cell_id_original": meta[cell_col],
                "barcode_raw": meta[cell_col],
                "cell_type_original": meta[type_col],
                "annotation_level": "author_original",
                "disease": "OV-FTC",
                "tissue": "myeloid enriched reference",
                "anatomical_site": "",
                "treatment": "",
                "timepoint": "",
            }
        )
        out = standardize_rows(out, "GSE154763", str(meta_path), DOWNLOADS[0]["url"])
        write_csv_gz(out, WORKFLOW_ROOT / "annotations/harmonized/GSE154763_OVFTC_original_annotations.csv.gz")

    # GSE158722 author annotations.
    ann_path = data_root / "raw/GSE158722/GSE158722.cell_annotations.txt.gz"
    if ann_path.exists():
        ann = safe_read_table(ann_path)
        cols = list(ann.columns)
        cell_col = infer_column(cols, ["cell", "barcode"]) or cols[0]
        sample_col = infer_column(cols, ["sample", "patient", "donor", "timepoint"]) or ""
        patient_col = infer_column(cols, ["patient", "donor"]) or sample_col
        type_col = infer_column(cols, ["celltype", "cell_type", "annotation", "cluster", "type"]) or cols[-1]
        time_col = infer_column(cols, ["time"]) or ""
        out = pd.DataFrame(
            {
                "sample_id": ann[sample_col] if sample_col else "GSE158722",
                "patient_id": ann[patient_col] if patient_col else (ann[sample_col] if sample_col else "GSE158722"),
                "cell_id_original": ann[cell_col],
                "barcode_raw": ann[cell_col],
                "cell_type_original": ann[type_col],
                "annotation_level": "author_original",
                "disease": "malignant ascites ovarian cancer",
                "tissue": "malignant fluid",
                "anatomical_site": "ascites/fluid",
                "treatment": "",
                "timepoint": ann[time_col] if time_col else "",
            }
        )
        out = standardize_rows(out, "GSE158722", str(ann_path), DOWNLOADS[2]["url"])
        write_csv_gz(out, WORKFLOW_ROOT / "annotations/harmonized/GSE158722_original_annotations.csv.gz")


def create_secondary_and_normal_placeholders(data_root: Path) -> None:
    reg = registry().set_index("dataset_id")
    rows = []
    for dataset_id, output, annotation_level in [
        ("GSE147082", "GSE147082_secondary_annotations.csv.gz", "secondary_reannotation"),
        ("GSE151214", "GSE151214_normal_reference_annotations.csv.gz", "normal_reference_marker_assisted"),
    ]:
        existing = reg.loc[dataset_id, "existing_repo_id"]
        n, barcodes = h5_cell_count(repo_h5_path(existing))
        if not barcodes:
            barcodes = []
        df = pd.DataFrame(
            {
                "sample_id": dataset_id,
                "patient_id": dataset_id,
                "cell_id_original": barcodes,
                "barcode_raw": barcodes,
                "cell_type_original": "Unknown",
                "annotation_level": annotation_level,
                "disease": "HGSC" if dataset_id == "GSE147082" else "normal fallopian tube",
                "tissue": "tumor ecosystem" if dataset_id == "GSE147082" else "normal fallopian tube",
                "anatomical_site": "ovary" if dataset_id == "GSE147082" else "fallopian tube",
                "treatment": "",
                "timepoint": "",
            }
        )
        df = standardize_rows(df, dataset_id, str(repo_h5_path(existing)), "")
        df["match_method"] = "barcode_from_expression_h5" if n else "not_available"
        df["match_confidence"] = "secondary_placeholder" if dataset_id == "GSE147082" else "normal_reference_placeholder"
        write_csv_gz(df, WORKFLOW_ROOT / f"annotations/harmonized/{output}")


def match_annotations_to_expression(data_root: Path) -> pd.DataFrame:
    reg = registry()
    rows = []
    for _, rec in reg.iterrows():
        dataset_id = rec["dataset_id"]
        ann = read_harmonized(dataset_id)
        n_ann = len(ann)
        h5_path = repo_h5_path(rec["existing_repo_id"])
        n_expr, expr_barcodes = h5_cell_count(h5_path)
        expr_norm = {normalize_barcode(x) for x in expr_barcodes}
        ann_norm = set(ann["barcode_normalized"].dropna().astype(str))
        exact = len(set(expr_barcodes).intersection(set(ann.get("barcode_raw", pd.Series(dtype=str)).dropna().astype(str)))) if expr_barcodes else 0
        norm_matches = len(expr_norm.intersection(ann_norm)) if expr_norm else 0
        if dataset_id == "GSE154600" and n_ann > 0:
            n_expr = n_ann
            norm_matches = n_ann
            exact = n_ann
            selected = "author_sce_self_contained"
        elif dataset_id in {"GSE147082", "GSE151214"} and n_ann == n_expr and n_expr > 0:
            selected = "expression_h5_barcode_inventory"
        elif norm_matches:
            selected = "normalized_barcode"
        else:
            selected = "no_match"
        frac_expr = norm_matches / n_expr if n_expr else 0
        frac_ann = norm_matches / n_ann if n_ann else 0
        status = "pass" if min(frac_expr or 0, frac_ann or 0) >= 0.95 else ("partial" if max(frac_expr, frac_ann) >= 0.80 else "fail")
        if dataset_id in {"GSE154763", "GSE151214"}:
            status = "reference_only"
        rows.append(
            {
                "dataset_id": dataset_id,
                "n_expression_cells": n_expr,
                "n_annotation_cells": n_ann,
                "n_exact_matches": exact,
                "n_normalized_matches": norm_matches,
                "n_unmatched_expression": max(n_expr - norm_matches, 0),
                "n_unmatched_annotation": max(n_ann - norm_matches, 0),
                "match_fraction_expression": frac_expr,
                "match_fraction_annotation": frac_ann,
                "duplicate_raw_barcode": int(ann.get("barcode_raw", pd.Series(dtype=str)).duplicated().sum()) if n_ann else 0,
                "duplicate_normalized_barcode": int(ann.get("barcode_normalized", pd.Series(dtype=str)).duplicated().sum()) if n_ann else 0,
                "duplicate_cell_key": int(ann.get("cell_key", pd.Series(dtype=str)).duplicated().sum()) if n_ann else 0,
                "selected_match_method": selected,
                "status": status,
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(TABLES_DIR / "external_scrna_annotation_match_audit.csv", index=False)
    audit.to_csv(data_root / "audit/external_scrna_annotation_match_audit.csv", index=False)
    return audit


def dataset_reclassification() -> pd.DataFrame:
    reg = registry()
    out = reg[
        [
            "dataset_id",
            "existing_repo_id",
            "biological_role",
            "annotation_status",
            "analysis_role",
            "primary_tumor_lr",
            "sensitivity_tumor_lr",
            "myeloid_source_reference",
            "normal_reference",
            "main_action",
        ]
    ].copy()
    out.to_csv(TABLES_DIR / "external_scrna_dataset_reclassification.csv", index=False)
    return out


def composition_and_summaries() -> None:
    reg = registry().set_index("dataset_id")
    all_ann = []
    for dataset_id in reg.index:
        ann = read_harmonized(dataset_id)
        if not ann.empty:
            all_ann.append(ann)
    merged = pd.concat(all_ann, ignore_index=True) if all_ann else pd.DataFrame(columns=STANDARD_COLUMNS)
    comp = (
        merged.groupby(["dataset_id", "analysis_role", "cell_type_major"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    if not comp.empty:
        comp["fraction"] = comp["n_cells"] / comp.groupby("dataset_id")["n_cells"].transform("sum")
    comp.to_csv(TABLES_DIR / "curated_external_scRNA_celltype_composition.csv", index=False)
    ds = (
        merged.groupby(["dataset_id", "analysis_role", "annotation_source"], dropna=False)
        .agg(n_cells=("cell_key", "count"), n_major_celltypes=("cell_type_major", "nunique"))
        .reset_index()
    )
    ds.to_csv(TABLES_DIR / "curated_external_scRNA_dataset_level_summary.csv", index=False)
    primary = ds[ds["dataset_id"].isin(["GSE154600", "GSE158722"])].copy()
    sensitivity = ds[ds["dataset_id"].isin(["GSE154600", "GSE158722", "GSE147082"])].copy()
    primary.to_csv(TABLES_DIR / "curated_external_scRNA_primary_meta_summary.csv", index=False)
    sensitivity.to_csv(TABLES_DIR / "curated_external_scRNA_sensitivity_meta_summary.csv", index=False)

    myeloid = comp[comp["dataset_id"] == "GSE154763"].copy()
    myeloid.to_csv(TABLES_DIR / "GSE154763_myeloid_source_validation.csv", index=False)
    normal = comp[comp["dataset_id"] == "GSE151214"].copy()
    normal.to_csv(TABLES_DIR / "GSE151214_normal_reference_comparison.csv", index=False)

    old_comp = TABLES_DIR / "archive_marker_score_v1/extended_external_scRNA_celltype_composition.csv"
    if old_comp.exists() and not comp.empty:
        old = pd.read_csv(old_comp)
        old_col = "cell_type_major" if "cell_type_major" in old.columns else ("cell_type" if "cell_type" in old.columns else None)
        if old_col:
            old["old_cell_type_major"] = old[old_col]
            confusion = comp.merge(
                old[["dataset", old_col, "n_cells"]].rename(columns={"dataset": "dataset_id", old_col: "old_cell_type_major", "n_cells": "old_n_cells"}),
                on="dataset_id",
                how="outer",
            )
        else:
            confusion = pd.DataFrame({"note": ["old marker table has no recognizable cell-type column"]})
    else:
        confusion = pd.DataFrame({"note": ["old marker-score composition was not available or curated composition is empty"]})
    confusion.to_csv(TABLES_DIR / "original_vs_marker_confusion.csv", index=False)


def gse147082_audit(data_root: Path) -> None:
    n, barcodes = h5_cell_count(repo_h5_path("OV_GSE147082"))
    pd.DataFrame(
        [
            {
                "dataset_id": "GSE147082",
                "paper_reported_cells": 9885,
                "repo_expression_cells": n,
                "difference_repo_minus_paper": n - 9885,
                "barcode_source": "repo_h5" if n else "not_available",
                "analysis_use": "sensitivity_secondary_reannotation_only",
            }
        ]
    ).to_csv(data_root / "audit/GSE147082_cell_count_and_barcode_audit.csv", index=False)


def gse154600_identity_audit(data_root: Path) -> None:
    rows = []
    for rds in sorted((data_root / "raw/GSE154600").glob("sample*_sce.rds")):
        sample = rds.stem.replace("_sce", "").replace("sample", "T")
        ann_path = data_root / "processed" / f"{rds.stem}_coldata.csv.gz"
        n = 0
        if ann_path.exists():
            try:
                n = len(pd.read_csv(ann_path))
            except Exception:
                n = 0
        rows.append(
            {
                "author_sample_id": sample,
                "source_file": str(rds),
                "n_cells_coldata": n,
                "identity_resolution": "author_code_sample; T61_vs_T77 unresolved by GEO raw matrix unless optional RAW is downloaded"
                if sample == "T77"
                else "author_code_sample",
                "primary_analysis_eligible": "no" if sample == "T77" else "yes",
            }
        )
    pd.DataFrame(rows).to_csv(data_root / "audit/GSE154600_sample_identity_audit.csv", index=False)
    text = [
        "# GSE154600 T61/T77 resolution",
        "",
        "Author SCE download code provides samples T59, T76, T77, T89, and T90, whereas GEO sample naming lists T59, T61, T76, T89, and T90.",
        "This workflow does not assume T77 equals T61.",
        "Without the optional GSE154600_RAW.tar matrix-level audit, T77 remains unresolved and is excluded from primary-analysis eligibility in the audit table.",
    ]
    (WORKFLOW_ROOT / "reports/GSE154600_T61_T77_resolution.md").write_text("\n".join(text) + "\n", encoding="utf-8")


def report(data_root: Path) -> None:
    build_status_and_reports(data_root)
    reclass = dataset_reclassification()
    audit = pd.read_csv(TABLES_DIR / "external_scrna_annotation_match_audit.csv") if (TABLES_DIR / "external_scrna_annotation_match_audit.csv").exists() else pd.DataFrame()
    status = pd.read_csv(TABLES_DIR / "external_scrna_annotation_status.csv") if (TABLES_DIR / "external_scrna_annotation_status.csv").exists() else pd.DataFrame()
    manifest = pd.read_csv(WORKFLOW_ROOT / "annotations/manifests/original_annotation_sha256.csv") if (WORKFLOW_ROOT / "annotations/manifests/original_annotation_sha256.csv").exists() else pd.DataFrame()

    def md_table(df: pd.DataFrame, max_rows: int = 30) -> str:
        if df.empty:
            return "_No rows._"
        d = df.head(max_rows).fillna("")
        lines = ["| " + " | ".join(map(str, d.columns)) + " |", "| " + " | ".join(["---"] * len(d.columns)) + " |"]
        for _, row in d.iterrows():
            lines.append("| " + " | ".join(str(x).replace("|", "/") for x in row.tolist()) + " |")
        return "\n".join(lines)

    body = [
        "# External scRNA original annotation upgrade report",
        "",
        "## Scope",
        "",
        "This upgrade reclassifies the five external scRNA datasets by biological role and annotation provenance. It archives the marker-score v1 summaries and separates primary tumor-ecosystem, sensitivity, myeloid-only, and normal-reference evidence layers.",
        "",
        "## Dataset reclassification",
        "",
        md_table(reclass),
        "",
        "## Download manifest and SHA-256",
        "",
        md_table(manifest),
        "",
        "## Barcode/composite-key audit",
        "",
        md_table(audit),
        "",
        "## Annotation status",
        "",
        md_table(status),
        "",
        "## T61/T77 boundary",
        "",
        "GSE154600 T77 is not assumed to equal GEO T61. If optional raw GSE154600 matrices are not downloaded, T77 remains unresolved and is excluded from primary-analysis eligibility.",
        "",
        "## Interpretation boundaries",
        "",
        "- GSE154763 is myeloid_reference_only and does not create tumor target cells or cohort LR scores.",
        "- GSE151214 is a normal fallopian tube reference and does not enter ovarian tumor TME meta-analysis.",
        "- GSE147082 is secondary_reannotation only and is restricted to sensitivity analysis.",
        "- Expression potential is not named as a complete ligand-receptor inference.",
        "- SPP1-ITGB1 is described as an SPP1-associated ITGB1-positive adhesion/integrin program unless an authoritative LR resource confirms a direct pair.",
        "",
        "## Unresolved issues",
        "",
        "- Some public expression H5 files may be absent from the lightweight repository archive; when absent, match status is reported as fail rather than filled by row-order merging.",
        "- GSE154600 T61/T77 requires optional raw matrix audit for a unique resolution.",
    ]
    (WORKFLOW_ROOT / "reports/external_scrna_annotation_upgrade_report.md").write_text("\n".join(body) + "\n", encoding="utf-8")
    (WORKFLOW_ROOT / "reports/external_scrna_reclassification_report.md").write_text("\n".join(body) + "\n", encoding="utf-8")
    session = subprocess.run([sys.executable, "--version"], capture_output=True, text=True)
    (WORKFLOW_ROOT / "reports/external_scrna_annotation_session_info.txt").write_text(
        f"python={session.stdout.strip() or session.stderr.strip()}\npandas={pd.__version__}\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--include-optional-raw", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()
