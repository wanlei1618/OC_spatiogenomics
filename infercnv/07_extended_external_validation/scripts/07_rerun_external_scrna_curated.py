#!/usr/bin/env python3
from external_annotation_common import DEFAULT_DATA_ROOT, archive_marker_results, composition_and_summaries, dataset_reclassification, progress


def main() -> int:
    progress(DEFAULT_DATA_ROOT, "07 rerun curated external scRNA summaries started")
    archive_marker_results()
    dataset_reclassification()
    composition_and_summaries()
    progress(DEFAULT_DATA_ROOT, "07 rerun curated external scRNA summaries finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
