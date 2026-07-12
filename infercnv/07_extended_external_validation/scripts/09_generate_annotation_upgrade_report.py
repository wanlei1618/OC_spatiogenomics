#!/usr/bin/env python3
from external_annotation_common import DEFAULT_DATA_ROOT, report, progress


def main() -> int:
    progress(DEFAULT_DATA_ROOT, "09 generate annotation upgrade report started")
    report(DEFAULT_DATA_ROOT)
    progress(DEFAULT_DATA_ROOT, "09 generate annotation upgrade report finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
