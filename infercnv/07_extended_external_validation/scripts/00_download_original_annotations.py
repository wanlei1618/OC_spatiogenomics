#!/usr/bin/env python3
from external_annotation_common import download, parse_args, progress


def main() -> int:
    args = parse_args()
    progress(args.data_root, "00 download original annotations started")
    df = download(args.data_root, args.include_optional_raw, args.force, args.dry_run)
    failed = df[(df["status"] == "failed") & (df["required"] == True)]
    progress(args.data_root, f"00 download original annotations finished: {len(df) - len(failed)}/{len(df)} available")
    return 1 if len(failed) else 0


if __name__ == "__main__":
    raise SystemExit(main())
