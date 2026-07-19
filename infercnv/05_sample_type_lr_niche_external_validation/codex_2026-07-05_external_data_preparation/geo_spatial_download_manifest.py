import csv
import gzip
import html.parser
import os
import re
import shutil
import tarfile
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path


ACCESSIONS = ["GSE211956", "GSE189843", "GSE203612", "GSE227019"]
BASE = Path("data") / "ovarian_spatial_geo"
TENX_OFFICIAL = [
    "https://www.10xgenomics.com/datasets/human-ovarian-cancer-1-standard",
    "https://www.10xgenomics.com/datasets/human-ovarian-cancer-11-mm-capture-area-ffpe-2-standard",
    "https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer",
]


class LinkParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        attrs = dict(attrs)
        href = attrs.get("href")
        if href and href not in ("../", "/"):
            self.links.append(href)


def gse_prefix(accession: str) -> str:
    return accession[:6] + "nnn"


def geo_url(accession: str, kind: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_prefix(accession)}/{accession}/{kind}/"


def list_url(url: str):
    with urllib.request.urlopen(url, timeout=60) as response:
        html = response.read().decode("utf-8", errors="replace")
    parser = LinkParser()
    parser.feed(html)
    files = []
    for href in parser.links:
        if href.endswith("/"):
            continue
        parsed = urllib.parse.urlparse(href)
        if parsed.scheme and "ftp.ncbi.nlm.nih.gov" not in parsed.netloc:
            continue
        name = urllib.parse.unquote(href)
        files.append((name, urllib.parse.urljoin(url, href)))
    return files


def remote_size(url: str):
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as response:
            val = response.headers.get("Content-Length")
        return int(val) if val else None
    except Exception:
        return None


def download(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return "existing"
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as response, open(tmp, "wb") as out:
                shutil.copyfileobj(response, out)
            tmp.replace(dest)
            return "downloaded"
        except Exception:
            if tmp.exists():
                tmp.unlink()
            if attempt == 2:
                raise
            time.sleep(3 * (attempt + 1))
    return "failed"


def safe_extract_tar(path: Path, dest: Path):
    dest = dest.resolve()
    with tarfile.open(path, "r:*") as tar:
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if not str(target).startswith(str(dest)):
                raise RuntimeError(f"Unsafe tar path: {member.name}")
        tar.extractall(dest)


def safe_extract_zip(path: Path, dest: Path):
    dest = dest.resolve()
    with zipfile.ZipFile(path) as zf:
        for member in zf.namelist():
            target = (dest / member).resolve()
            if not str(target).startswith(str(dest)):
                raise RuntimeError(f"Unsafe zip path: {member}")
        zf.extractall(dest)


def decompress_file(path: Path):
    name = path.name.lower()
    try:
        if name.endswith((".tar.gz", ".tgz", ".tar")):
            marker = path.with_name(path.name + ".extracted")
            if not marker.exists():
                safe_extract_tar(path, path.parent)
                marker.write_text("extracted\n", encoding="utf-8")
            return "tar_extracted"
        if name.endswith(".zip"):
            marker = path.with_name(path.name + ".extracted")
            if not marker.exists():
                safe_extract_zip(path, path.parent)
                marker.write_text("extracted\n", encoding="utf-8")
            return "zip_extracted"
        if name.endswith(".gz"):
            if name.endswith(("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz", "genes.tsv.gz")):
                return "kept_10x_gz"
            out = path.with_suffix("")
            if not out.exists():
                with gzip.open(path, "rb") as src, open(out, "wb") as dst:
                    shutil.copyfileobj(src, dst)
            return "gz_decompressed"
    except Exception as exc:
        return f"decompress_failed:{exc}"
    return "not_archive"


def postprocess_nested_archives(root: Path):
    statuses = []
    for path in list(root.rglob("*")):
        if not path.is_file():
            continue
        name = path.name.lower()
        try:
            if name.endswith(".zip"):
                out_dir = path.with_suffix("")
                marker = out_dir / ".extracted"
                if not marker.exists():
                    out_dir.mkdir(parents=True, exist_ok=True)
                    safe_extract_zip(path, out_dir)
                    marker.write_text("extracted\n", encoding="utf-8")
                statuses.append((str(path.relative_to(root)), "zip_extracted_nested"))
            elif name.endswith(".gz") and not name.endswith((".tar.gz", ".tgz")):
                if name.endswith(("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz", "genes.tsv.gz")):
                    statuses.append((str(path.relative_to(root)), "kept_10x_gz"))
                    continue
                out = path.with_suffix("")
                if not out.exists():
                    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                statuses.append((str(path.relative_to(root)), "gz_decompressed_nested"))
        except Exception as exc:
            statuses.append((str(path.relative_to(root)), f"nested_decompress_failed:{exc}"))
    return statuses


def parse_soft(soft_dir: Path):
    soft_files = list(soft_dir.glob("*.soft")) + list(soft_dir.glob("*.soft.gz"))
    samples = {}
    current = None
    for soft in soft_files:
        opener = gzip.open if soft.suffix == ".gz" else open
        mode = "rt"
        with opener(soft, mode, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith("^SAMPLE = "):
                    gsm = line.split("=", 1)[1].strip()
                    current = samples.setdefault(gsm, {"gsm": gsm, "title": "", "platform": "", "supplementary": []})
                elif current is not None and line.startswith("!Sample_title = "):
                    current["title"] = line.split("=", 1)[1].strip()
                elif current is not None and line.startswith("!Sample_platform_id = "):
                    current["platform"] = line.split("=", 1)[1].strip()
                elif current is not None and line.startswith("!Sample_supplementary_file"):
                    current["supplementary"].append(line.split("=", 1)[1].strip())
    return samples


def classify_files(root: Path):
    matrix_patterns = re.compile(r"(filtered_feature_bc_matrix|raw_feature_bc_matrix|matrix\.mtx|features\.tsv|genes\.tsv|barcodes\.tsv|\.h5$)", re.I)
    image_patterns = re.compile(r"(spatial|tissue_hires|tissue_lowres|\.tif$|\.tiff$|\.jpg$|\.jpeg$|\.png$|scalefactors|tissue_positions)", re.I)
    files = [p for p in root.rglob("*") if p.is_file() and not p.name.endswith(".part")]
    matrix = [str(p.relative_to(root)) for p in files if matrix_patterns.search(p.name) or matrix_patterns.search(str(p.parent))]
    images = [str(p.relative_to(root)) for p in files if image_patterns.search(p.name) or image_patterns.search(str(p.parent))]
    return matrix, images


def is_visium_compatible(matrix_files, image_files):
    joined_m = "\n".join(matrix_files).lower()
    joined_i = "\n".join(image_files).lower()
    has_matrix = ("filtered_feature_bc_matrix.h5" in joined_m) or (
        "matrix.mtx" in joined_m and "barcodes.tsv" in joined_m and ("features.tsv" in joined_m or "genes.tsv" in joined_m)
    )
    has_spatial = ("spatial" in joined_i and "scalefactors" in joined_i and "tissue_positions" in joined_i) or (
        "tissue_hires_image" in joined_i and "scalefactors" in joined_i
    )
    return bool(has_matrix and has_spatial)


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    download_rows = []
    for accession in ACCESSIONS:
        for kind in ("suppl", "soft"):
            out_dir = BASE / accession / kind
            out_dir.mkdir(parents=True, exist_ok=True)
            url = geo_url(accession, kind)
            print(f"Listing {url}")
            try:
                files = list_url(url)
            except Exception as exc:
                download_rows.append({"accession": accession, "kind": kind, "file": "", "url": url, "status": f"list_failed:{exc}", "bytes": ""})
                continue
            for name, file_url in files:
                size = remote_size(file_url)
                dest = out_dir / name
                print(f"Downloading {accession}/{kind}/{name}")
                try:
                    status = download(file_url, dest)
                    decomp = decompress_file(dest)
                except Exception as exc:
                    status = f"download_failed:{exc}"
                    decomp = ""
                download_rows.append({
                    "accession": accession,
                    "kind": kind,
                    "file": name,
                    "url": file_url,
                    "status": status,
                    "decompress_status": decomp,
                    "bytes": size if size is not None else "",
                })

    with open(BASE / "download_manifest.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "kind", "file", "url", "status", "decompress_status", "bytes"])
        writer.writeheader()
        writer.writerows(download_rows)

    nested_rows = []
    for accession in ACCESSIONS:
        for rel, status in postprocess_nested_archives(BASE / accession / "suppl"):
            nested_rows.append({"accession": accession, "file": rel, "status": status})
        for rel, status in postprocess_nested_archives(BASE / accession / "soft"):
            nested_rows.append({"accession": accession, "file": rel, "status": status})
    with open(BASE / "decompression_manifest.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "file", "status"])
        writer.writeheader()
        writer.writerows(nested_rows)

    sample_rows = []
    for accession in ACCESSIONS:
        acc_dir = BASE / accession
        samples = parse_soft(acc_dir / "soft")
        matrix_files, image_files = classify_files(acc_dir / "suppl")
        visium = is_visium_compatible(matrix_files, image_files)
        if samples:
            for gsm, sample in samples.items():
                supp_hint = " ".join(sample.get("supplementary", []))
                sample_matrix = [f for f in matrix_files if gsm in f or sample["title"] in f or Path(f).name in supp_hint]
                sample_images = [f for f in image_files if gsm in f or sample["title"] in f or Path(f).name in supp_hint]
                if not sample_matrix:
                    sample_matrix = matrix_files
                if not sample_images:
                    sample_images = image_files
                sample_rows.append({
                    "source": "GEO",
                    "accession": accession,
                    "sample_name": sample["title"] or gsm,
                    "gsm": gsm,
                    "platform": sample["platform"],
                    "available_matrix_files": ";".join(sample_matrix),
                    "spatial_image_files": ";".join(sample_images),
                    "visium_compatible_Load10X_Spatial": is_visium_compatible(sample_matrix, sample_images),
                    "notes": "",
                })
        else:
            sample_rows.append({
                "source": "GEO",
                "accession": accession,
                "sample_name": "",
                "gsm": "",
                "platform": "",
                "available_matrix_files": ";".join(matrix_files),
                "spatial_image_files": ";".join(image_files),
                "visium_compatible_Load10X_Spatial": visium,
                "notes": "SOFT sample records not parsed",
            })

    for url in TENX_OFFICIAL:
        sample_rows.append({
            "source": "10x official",
            "accession": "",
            "sample_name": url.rstrip("/").split("/")[-1],
            "gsm": "",
            "platform": "10x Genomics public dataset page",
            "available_matrix_files": "",
            "spatial_image_files": "",
            "visium_compatible_Load10X_Spatial": "manual_check_required",
            "notes": url,
        })

    with open(BASE / "ovarian_spatial_geo_sample_manifest.csv", "w", newline="", encoding="utf-8") as handle:
        fields = ["source", "accession", "sample_name", "gsm", "platform", "available_matrix_files",
                  "spatial_image_files", "visium_compatible_Load10X_Spatial", "notes"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sample_rows)

    print(f"Wrote {BASE / 'download_manifest.csv'}")
    print(f"Wrote {BASE / 'ovarian_spatial_geo_sample_manifest.csv'}")


if __name__ == "__main__":
    main()
