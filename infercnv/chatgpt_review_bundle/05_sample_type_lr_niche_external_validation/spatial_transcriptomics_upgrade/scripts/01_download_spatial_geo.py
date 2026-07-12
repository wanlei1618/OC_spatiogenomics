#!/usr/bin/env python3
r"""Download and curate ovarian spatial-transcriptomics inputs from GEO.

Default output root:
    D:\OC_spatiogenomics\spatial_data

The script is intentionally dependency-light and uses only the Python standard
library. It reads metadata/download_manifest.csv, downloads only curated targets,
decompresses GSE203612 Visium sidecar files into SpaceRanger-compatible folders,
and safely extracts the GSE189843 series archive.

Examples
--------
python 01_download_spatial_geo.py --dry-run
python 01_download_spatial_geo.py --root D:\OC_spatiogenomics\spatial_data
python 01_download_spatial_geo.py --root /data/OC_spatiogenomics/spatial_data
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tarfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable


DEFAULT_ROOT = Path(r"D:\OC_spatiogenomics\spatial_data")
USER_AGENT = "OC-spatiogenomics-spatial-curation/1.0"


@dataclass
class DownloadRecord:
    dataset: str
    target_id: str
    url: str
    output_path: str
    status: str
    bytes: int
    sha256: str
    message: str


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    default_manifest = here.parent / "metadata" / "download_manifest.csv"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=default_manifest)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        chunk = handle.read(chunk_size)
        while chunk:
            digest.update(chunk)
            chunk = handle.read(chunk_size)
    return digest.hexdigest()


def download_file(
    url: str,
    destination: Path,
    retries: int,
    timeout: int,
    force: bool,
) -> tuple[str, int, str]:
    """Download to a temporary file and atomically rename on success."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return "cached", destination.stat().st_size, sha256_file(destination)

    part = destination.with_suffix(destination.suffix + ".part")
    if part.exists():
        part.unlink()

    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response, part.open("wb") as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            if part.stat().st_size == 0:
                raise IOError("Downloaded file is empty")
            os.replace(part, destination)
            return "downloaded", destination.stat().st_size, sha256_file(destination)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            last_error = exc
            if part.exists():
                part.unlink()
            if attempt < retries:
                time.sleep(min(2 ** attempt, 30))

    raise RuntimeError(f"Failed to download {url}: {last_error}")


def gunzip_to(source: Path, destination: Path, force: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    temp = destination.with_suffix(destination.suffix + ".part")
    with gzip.open(source, "rb") as src, temp.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    os.replace(temp, destination)


def safe_extract_tar(archive: Path, destination: Path, force: bool) -> None:
    """Extract a TAR archive while blocking path traversal."""
    destination.mkdir(parents=True, exist_ok=True)
    marker = destination / ".extracted.ok"
    if marker.exists() and not force:
        return

    root = destination.resolve()
    with tarfile.open(archive, "r:*") as tar:
        for member in tar.getmembers():
            target = (destination / member.name).resolve()
            if target != root and root not in target.parents:
                raise RuntimeError(f"Unsafe path in TAR: {member.name}")
        tar.extractall(destination)

    for gz_path in destination.rglob("*.gz"):
        output = gz_path.with_suffix("")
        gunzip_to(gz_path, output, force=force)

    marker.write_text("ok\n", encoding="utf-8")


def read_manifest(path: Path) -> Iterable[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"dataset", "target_id", "include", "url", "relative_path", "compression"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
        for row in reader:
            if row["include"].strip().upper() == "TRUE":
                yield row


def main() -> int:
    args = parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    records: list[DownloadRecord] = []

    for row in read_manifest(args.manifest):
        output = args.root / Path(row["relative_path"])
        compressed_download = output
        if row["compression"] == "gzip":
            compressed_download = output.with_suffix(output.suffix + ".gz")

        if args.dry_run:
            print(f"[DRY RUN] {row['target_id']}: {row['url']} -> {output}")
            continue

        try:
            status, size, checksum = download_file(
                row["url"],
                compressed_download,
                retries=args.retries,
                timeout=args.timeout,
                force=args.force,
            )
            if row["compression"] == "gzip":
                gunzip_to(compressed_download, output, force=args.force)
                final_size = output.stat().st_size
                final_sha = sha256_file(output)
            elif row["compression"] == "tar":
                extract_dir = args.root / "raw" / row["dataset"] / "extracted"
                safe_extract_tar(output, extract_dir, force=args.force)
                final_size = size
                final_sha = checksum
            else:
                final_size = size
                final_sha = checksum

            records.append(
                DownloadRecord(
                    dataset=row["dataset"],
                    target_id=row["target_id"],
                    url=row["url"],
                    output_path=str(output),
                    status=status,
                    bytes=final_size,
                    sha256=final_sha,
                    message="ok",
                )
            )
            print(f"[OK] {row['target_id']} -> {output}")
        except Exception as exc:
            records.append(
                DownloadRecord(
                    dataset=row["dataset"],
                    target_id=row["target_id"],
                    url=row["url"],
                    output_path=str(output),
                    status="failed",
                    bytes=0,
                    sha256="",
                    message=str(exc),
                )
            )
            print(f"[FAILED] {row['target_id']}: {exc}")

    if not args.dry_run:
        log_dir = args.root / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / "download_audit.json"
        log_path.write_text(
            json.dumps([asdict(record) for record in records], indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        failed = [record for record in records if record.status == "failed"]
        print(f"Audit log: {log_path}")
        return 1 if failed else 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
