#!/usr/bin/env python3
from external_annotation_common import DEFAULT_DATA_ROOT, gse147082_audit, match_annotations_to_expression, progress


def main() -> int:
    progress(DEFAULT_DATA_ROOT, "04 barcode/composite-key matching started")
    audit = match_annotations_to_expression(DEFAULT_DATA_ROOT)
    gse147082_audit(DEFAULT_DATA_ROOT)
    progress(DEFAULT_DATA_ROOT, f"04 matching finished for {len(audit)} datasets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
