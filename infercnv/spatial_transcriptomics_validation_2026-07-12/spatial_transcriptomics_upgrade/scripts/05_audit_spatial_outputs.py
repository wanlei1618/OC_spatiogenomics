#!/usr/bin/env python3
r"""Audit curated spatial-transcriptomics metadata and output tables.

The checks prevent three recurrent interpretation errors:
1. excluded/non-ovarian samples entering ovarian summaries;
2. expression-only datasets being used for coordinate-aware claims;
3. missing or duplicated curated samples.

Usage:
    python 05_audit_spatial_outputs.py
    python 05_audit_spatial_outputs.py --results D:\OC_spatiogenomics\spatial_data\results\spatial_curated
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_bool(value: str) -> bool:
    return value.strip().upper() in {"TRUE", "T", "1", "YES", "Y"}


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=here.parent / "metadata" / "spatial_sample_manifest.csv",
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=Path(r"D:\OC_spatiogenomics\spatial_data\results\spatial_curated"),
    )
    parser.add_argument("--strict-results", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = load_csv(args.manifest)
    errors: list[str] = []
    warnings: list[str] = []

    sample_ids = [row["sample_id"] for row in samples]
    duplicates = sorted({sample for sample in sample_ids if sample_ids.count(sample) > 1})
    if duplicates:
        errors.append(f"Duplicate sample IDs in manifest: {duplicates}")

    included = {
        row["sample_id"]: row
        for row in samples
        if parse_bool(row["include_in_ovarian_analysis"])
    }
    excluded = {
        row["sample_id"]: row
        for row in samples
        if not parse_bool(row["include_in_ovarian_analysis"])
    }

    expected_gse203612 = {"GSM6177614", "GSM6177617"}
    observed_gse203612 = {
        row["sample_id"]
        for row in included.values()
        if row["dataset"] == "GSE203612"
    }
    if observed_gse203612 != expected_gse203612:
        errors.append(
            f"GSE203612 ovarian set must be {sorted(expected_gse203612)}, "
            f"observed {sorted(observed_gse203612)}"
        )

    if "GSM6177618" not in excluded:
        errors.append("GSM6177618 must be explicitly excluded as PDAC")

    expected_gse189843 = {f"GSM{gsm}" for gsm in range(5708485, 5708497)}
    observed_gse189843 = {
        row["sample_id"]
        for row in included.values()
        if row["dataset"] == "GSE189843"
    }
    if observed_gse189843 != expected_gse189843:
        errors.append(
            "GSE189843 must contain all 12 HGSC samples "
            f"({sorted(expected_gse189843 - observed_gse189843)} missing)"
        )

    for row in included.values():
        if row["dataset"] == "GSE189843" and row["analysis_level"] != "expression_only":
            errors.append(
                f"{row['sample_id']} must remain expression_only until coordinates are obtained"
            )

    result_files = {
        "correlation": args.results / "spatial_correlation_curated.csv",
        "neighborhood": args.results / "spatial_neighborhood_enrichment_curated.csv",
    }
    for label, path in result_files.items():
        if not path.exists():
            message = f"Result file not found: {path}"
            (errors if args.strict_results else warnings).append(message)
            continue

        rows = load_csv(path)
        seen = {row.get("sample_id", "") for row in rows}
        leaked = sorted(seen.intersection(excluded))
        if leaked:
            errors.append(f"Excluded samples present in {label} results: {leaked}")

        if label == "neighborhood":
            invalid = [
                row.get("sample_id", "")
                for row in rows
                if included.get(row.get("sample_id", ""), {}).get("coordinate_status")
                != "available"
            ]
            if invalid:
                errors.append(
                    f"Coordinate-unavailable samples present in neighborhood results: {invalid}"
                )

    report = {
        "manifest": str(args.manifest),
        "results": str(args.results),
        "included_samples": sorted(included),
        "excluded_samples": sorted(excluded),
        "errors": errors,
        "warnings": warnings,
        "status": "failed" if errors else "passed",
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
