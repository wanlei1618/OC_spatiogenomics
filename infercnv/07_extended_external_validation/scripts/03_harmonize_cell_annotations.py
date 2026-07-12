#!/usr/bin/env python3
from external_annotation_common import DEFAULT_DATA_ROOT, copy_task_registry, create_secondary_and_normal_placeholders, progress


def main() -> int:
    progress(DEFAULT_DATA_ROOT, "03 harmonize cell annotations started")
    copy_task_registry()
    create_secondary_and_normal_placeholders(DEFAULT_DATA_ROOT)
    progress(DEFAULT_DATA_ROOT, "03 harmonize cell annotations finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
