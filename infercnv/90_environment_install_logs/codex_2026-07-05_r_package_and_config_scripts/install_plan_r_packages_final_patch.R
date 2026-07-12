options(download.file.method = "libcurl")
options(timeout = 1200)
Sys.setenv(TMP = "D:/TEMP", TEMP = "D:/TEMP", TMPDIR = "D:/TEMP")
dir.create("D:/TEMP", recursive = TRUE, showWarnings = FALSE)

lib <- "D:/Documents/R/win-library/4.0"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

repos <- c(
  BioCsoft = "https://bioconductor.org/packages/3.12/bioc",
  BioCann = "https://bioconductor.org/packages/3.12/data/annotation",
  BioCexp = "https://bioconductor.org/packages/3.12/data/experiment",
  CRAN = "https://cran.rstudio.com"
)
options(repos = repos)

log_file <- "D:/OC_spatiogenomics/infercnv/install_plan_tools_R_final_patch_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

install_archive <- function(pkg, version) {
  if (pkg %in% rownames(installed.packages())) {
    cat(pkg, "already installed", as.character(packageVersion(pkg)), "\n")
    return(invisible(TRUE))
  }
  cat("\nArchive:", pkg, version, "\n")
  tryCatch(remotes::install_version(pkg, version = version, lib = lib,
                                    repos = "https://cran.rstudio.com",
                                    upgrade = "never", dependencies = TRUE),
           error = function(e) cat("FAILED archive", pkg, ":", conditionMessage(e), "\n"))
}

install_current <- function(pkgs) {
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (!length(miss)) return(invisible(TRUE))
  cat("\nInstall current:", paste(miss, collapse = ", "), "\n")
  tryCatch(install.packages(miss, lib = lib, dependencies = TRUE, type = "source", Ncpus = 2),
           error = function(e) cat("FAILED current:", conditionMessage(e), "\n"))
}

install_github <- function(pkg, repo) {
  if (pkg %in% rownames(installed.packages())) {
    cat(pkg, "already installed", as.character(packageVersion(pkg)), "\n")
    return(invisible(TRUE))
  }
  cat("\nGitHub:", pkg, repo, "\n")
  tryCatch(remotes::install_github(repo, lib = lib, upgrade = "never",
                                  dependencies = TRUE, build_vignettes = FALSE),
           error = function(e) cat("FAILED github", pkg, ":", conditionMessage(e), "\n"))
}

cat("Patch install for plan packages\n")
cat("R:", R.version.string, "\n")
cat("lib:", lib, "\n")

install_archive("dplyr", "1.1.4")
install_archive("transport", "0.14-7")
install_archive("msigdbr", "7.5.1")
install_archive("scTenifoldNet", "1.2.4")
install_archive("scTenifoldKnk", "1.0.1")

install_current(c("SingleCellExperiment", "edgeR", "UCell", "scDblFinder"))

install_github("copykat", "navinlabcode/copykat")
install_github("liana", "saezlab/liana")

cat("\nFinal status\n")
pkgs <- c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "pheatmap", "remotes",
          "SingleCellExperiment", "scDblFinder", "ComplexHeatmap", "msigdbr", "UCell",
          "limma", "edgeR", "fgsea", "clusterProfiler", "org.Hs.eg.db",
          "liana", "CellChat", "nichenetr", "scTenifoldKnk", "copykat", "CaSpER")
ip <- installed.packages()
status <- data.frame(
  package = pkgs,
  installed = pkgs %in% rownames(ip),
  version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
  stringsAsFactors = FALSE
)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_plan_tools_R_status_final.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
