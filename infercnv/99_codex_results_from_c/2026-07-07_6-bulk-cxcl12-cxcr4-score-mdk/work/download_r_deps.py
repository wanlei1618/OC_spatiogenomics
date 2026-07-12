import os
import re
import shutil
import urllib.request
from pathlib import Path


ROOT = Path("work/r_packages/localrepo/bin/windows/contrib/4.0")
ROOT.mkdir(parents=True, exist_ok=True)

CRAN_PACKAGES = Path("work/r_packages/PACKAGES.txt")
BIOC_ROOTS = [
    "https://bioconductor.org/packages/3.12/bioc/bin/windows/contrib/4.0",
    "https://bioconductor.org/packages/3.12/data/annotation/bin/windows/contrib/4.0",
    "https://bioconductor.org/packages/3.12/data/experiment/bin/windows/contrib/4.0",
]

CELLCHAT_DEPS = [
    "future", "future.apply", "pbapply", "irlba", "NMF", "ggalluvial",
    "stringr", "svglite", "expm", "Rtsne", "ggrepel", "circlize",
    "RColorBrewer", "cowplot", "ComplexHeatmap", "RSpectra", "Rcpp",
    "RcppEigen", "reticulate", "scales", "sna", "forcats", "reshape2",
    "FNN", "shape", "BiocGenerics", "magrittr", "patchwork", "colorspace",
    "plyr", "ggpubr", "ggnetwork", "BiocNeighbors",
]

BASE = set("""
base compiler datasets graphics grDevices grid methods parallel splines stats
stats4 tcltk tools utils Matrix MASS lattice nlme mgcv survival boot class
cluster codetools foreign KernSmooth nnet rpart spatial
""".split())


def parse_packages_text(text):
    records = {}
    for block in re.split(r"\n\s*\n", text.strip()):
        rec = {}
        current = None
        for line in block.splitlines():
            if not line:
                continue
            if line[0].isspace() and current:
                rec[current] += " " + line.strip()
            elif ":" in line:
                k, v = line.split(":", 1)
                current = k
                rec[k] = v.strip()
        if "Package" in rec:
            records[rec["Package"]] = rec
    return records


def read_url(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def dep_names(value):
    if not value:
        return []
    value = re.sub(r"\([^)]*\)", "", value)
    out = []
    for part in value.replace("\n", " ").split(","):
        name = part.strip()
        if name:
            out.append(name)
    return out


def collect(seed, repos):
    need = set(seed)
    done = set()
    changed = True
    while changed:
        changed = False
        for pkg in list(need):
            if pkg in done or pkg in BASE:
                continue
            rec = None
            for repo in repos:
                if pkg in repo:
                    rec = repo[pkg]
                    break
            done.add(pkg)
            if not rec:
                print("MISSING", pkg)
                continue
            for field in ["Depends", "Imports", "LinkingTo"]:
                for dep in dep_names(rec.get(field, "")):
                    if dep not in need and dep not in BASE:
                        need.add(dep)
                        changed = True
    return sorted(p for p in need if p not in BASE)


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return
    print("download", url)
    with urllib.request.urlopen(url, timeout=120) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)


def main():
    cran = parse_packages_text(CRAN_PACKAGES.read_text(encoding="utf-8", errors="replace"))
    bioc_repos = []
    bioc_urls = []
    for root in BIOC_ROOTS:
        try:
            text = read_url(root + "/PACKAGES")
            repo = parse_packages_text(text)
            bioc_repos.append(repo)
            bioc_urls.append(root)
            print("bioc repo", root, len(repo))
        except Exception as exc:
            print("bioc repo failed", root, exc)
    repos = [cran] + bioc_repos
    pkgs = collect(CELLCHAT_DEPS, repos)
    print("packages", len(pkgs))
    missing = []
    for pkg in pkgs:
        if pkg in cran:
            rec = cran[pkg]
            fname = f"{pkg}_{rec['Version']}.zip"
            url = f"https://cran.r-project.org/bin/windows/contrib/4.0/{fname}"
        else:
            found = False
            for repo, root in zip(bioc_repos, bioc_urls):
                if pkg in repo:
                    rec = repo[pkg]
                    fname = f"{pkg}_{rec['Version']}.zip"
                    url = f"{root}/{fname}"
                    found = True
                    break
            if not found:
                missing.append(pkg)
                continue
        try:
            download(url, ROOT / fname)
        except Exception as exc:
            missing.append(pkg)
            print("DOWNLOAD_FAILED", pkg, exc)
    if missing:
        print("missing/download failed:", sorted(set(missing)))


if __name__ == "__main__":
    main()
