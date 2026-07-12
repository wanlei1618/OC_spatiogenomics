import csv
import os
import re
import shutil
from pathlib import Path


BASE = Path("data") / "ovarian_spatial_geo"


def link_or_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def choose(files, patterns):
    for pat in patterns:
        rx = re.compile(pat, re.I)
        hits = [p for p in files if rx.search(p.name)]
        if hits:
            return sorted(hits, key=lambda p: len(str(p)))[0]
    return None


def build_for_sample(accession: str, gsm: str):
    suppl = BASE / accession / "suppl"
    out = BASE / accession / "seurat_visium" / gsm
    spatial = out / "spatial"
    files = [p for p in suppl.rglob("*") if p.is_file() and gsm in str(p)]
    if not files:
        return "", "no_gsm_files_found"

    status = []

    h5 = choose(files, [r"filtered_feature_bc_matrix\.h5$"])
    if h5:
        link_or_copy(h5, out / "filtered_feature_bc_matrix.h5")
        status.append("h5")
    else:
        matrix = choose(files, [r"matrix.*\.mtx\.gz$", r"matrix\.mtx$"])
        features = choose(files, [r"features.*\.tsv\.gz$", r"features\.tsv$", r"genes.*\.tsv\.gz$", r"genes\.tsv$"])
        barcodes = choose(files, [r"barcodes.*\.tsv\.gz$", r"barcodes\.tsv$"])
        if matrix and features and barcodes:
            fdir = out / "filtered_feature_bc_matrix"
            link_or_copy(matrix, fdir / ("matrix.mtx.gz" if matrix.name.endswith(".gz") else "matrix.mtx"))
            link_or_copy(features, fdir / ("features.tsv.gz" if features.name.endswith(".gz") else "features.tsv"))
            link_or_copy(barcodes, fdir / ("barcodes.tsv.gz" if barcodes.name.endswith(".gz") else "barcodes.tsv"))
            status.append("mtx_triplet")

    spatial_map = {
        "tissue_hires_image.png": [r"tissue_hires_image\.png$"],
        "tissue_lowres_image.png": [r"tissue_lowres_image\.png$"],
        "scalefactors_json.json": [r"scalefactors.*\.json$"],
        "tissue_positions_list.csv": [r"tissue_positions_list\.csv$", r"tissue_positions\.csv$"],
        "aligned_fiducials.jpg": [r"aligned_fiducials\.jpg$"],
        "detected_tissue_image.jpg": [r"detected_tissue_image\.jpg$"],
    }
    for std_name, pats in spatial_map.items():
        src = choose(files, pats)
        if src:
            link_or_copy(src, spatial / std_name)
            status.append(std_name)
            if std_name == "tissue_positions_list.csv" and src.name.lower() == "tissue_positions.csv":
                link_or_copy(src, spatial / "tissue_positions.csv")

    matrix_ready = (out / "filtered_feature_bc_matrix.h5").exists() or (out / "filtered_feature_bc_matrix" / "matrix.mtx.gz").exists() or (out / "filtered_feature_bc_matrix" / "matrix.mtx").exists()
    spatial_ready = (spatial / "scalefactors_json.json").exists() and (spatial / "tissue_positions_list.csv").exists() and ((spatial / "tissue_hires_image.png").exists() or (spatial / "tissue_lowres_image.png").exists())
    ready = bool(matrix_ready and spatial_ready)
    return str(out).replace("\\", "/"), "ready" if ready else "incomplete:" + ";".join(status)


def main():
    manifest_path = BASE / "ovarian_spatial_geo_sample_manifest.csv"
    with open(manifest_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    for row in rows:
        if row.get("source") == "GEO" and row.get("visium_compatible_Load10X_Spatial") == "True":
            out, status = build_for_sample(row["accession"], row["gsm"])
            row["seurat_load10x_spatial_dir"] = out
            row["seurat_load10x_spatial_dir_status"] = status
        else:
            row["seurat_load10x_spatial_dir"] = ""
            row["seurat_load10x_spatial_dir_status"] = ""

    fields = list(rows[0].keys())
    with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    with open(BASE / "seurat_visium_dirs.csv", "w", newline="", encoding="utf-8") as handle:
        fields2 = ["accession", "gsm", "sample_name", "seurat_load10x_spatial_dir", "seurat_load10x_spatial_dir_status"]
        writer = csv.DictWriter(handle, fieldnames=fields2)
        writer.writeheader()
        for row in rows:
            if row.get("seurat_load10x_spatial_dir"):
                writer.writerow({k: row.get(k, "") for k in fields2})

    print(f"Updated {manifest_path}")
    print(f"Wrote {BASE / 'seurat_visium_dirs.csv'}")


if __name__ == "__main__":
    main()
