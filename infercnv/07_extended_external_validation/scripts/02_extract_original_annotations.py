#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

from external_annotation_common import (
    DEFAULT_DATA_ROOT,
    WORKFLOW_ROOT,
    extract_non_sce_annotations,
    gse154600_identity_audit,
    progress,
    standardize_rows,
    write_csv_gz,
)


def extract_gse154600(data_root: Path) -> None:
    rows = []
    for path in sorted((data_root / "processed").glob("sample*_sce_coldata.csv.gz")):
        cd = pd.read_csv(path)
        sample = path.name.replace("_sce_coldata.csv.gz", "").replace("sample", "T")
        cols = list(cd.columns)
        type_col = next((c for c in cols if c.lower() in {"celltype", "cell_type", "subtype", "cluster"}), None)
        if type_col is None:
            type_col = next((c for c in cols if "cell" in c.lower() and "type" in c.lower()), None)
        if type_col is None:
            type_col = next((c for c in cols if "cluster" in c.lower()), None)
        if type_col is None:
            type_col = cols[0]
        cell_col = "cell_id_original" if "cell_id_original" in cd.columns else cols[0]
        cell_ids = cd[cell_col].astype(str)
        if cell_ids.duplicated().any() or cell_ids.str.lower().isin(["", "nan", "none"]).any():
            cell_ids = [f"{sample}_cell_{i + 1}" for i in range(len(cd))]
        out = pd.DataFrame(
            {
                "sample_id": sample,
                "patient_id": sample,
                "cell_id_original": cell_ids,
                "barcode_raw": cell_ids,
                "cell_type_original": cd[type_col],
                "annotation_level": "author_original",
                "disease": "ovarian cancer",
                "tissue": "tumor ecosystem",
                "anatomical_site": "ovary",
                "treatment": "",
                "timepoint": "",
            }
        )
        rows.append(out)
    if rows:
        merged = pd.concat(rows, ignore_index=True)
        merged = standardize_rows(merged, "GSE154600", str(data_root / "raw/GSE154600"), "https://github.com/waldronlab/subtypeHeterogeneity")
        write_csv_gz(merged, WORKFLOW_ROOT / "annotations/harmonized/GSE154600_original_annotations.csv.gz")


def main() -> int:
    data_root = DEFAULT_DATA_ROOT
    progress(data_root, "02 extract original annotations started")
    extract_gse154600(data_root)
    extract_non_sce_annotations(data_root)
    gse154600_identity_audit(data_root)
    progress(data_root, "02 extract original annotations finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
