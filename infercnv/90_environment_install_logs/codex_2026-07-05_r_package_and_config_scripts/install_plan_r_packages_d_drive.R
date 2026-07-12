options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "libcurl")
options(timeout = 1200)

Sys.setenv(TMP = "D:/TEMP", TEMP = "D:/TEMP", TMPDIR = "D:/TEMP")
dir.create("D:/TEMP", recursive = TRUE, showWarnings = FALSE)

lib <- "D:/Documents/R/win-library/4.0"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_plan_tools_R_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("R version:", R.version.string, "\n")
cat("R library:", lib, "\n")
cat("Temp:", tempdir(), "\n")

install_cran_current <- function(pkgs) {
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (length(miss) == 0) return(invisible(TRUE))
  cat("\nCRAN current:", paste(miss, collapse = ", "), "\n")
  tryCatch(install.packages(miss, lib = lib, dependencies = TRUE, type = "source", Ncpus = 2),
           error = function(e) cat("FAILED current:", conditionMessage(e), "\n"))
}

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

install_bioc <- function(pkgs) {
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (length(miss) == 0) return(invisible(TRUE))
  cat("\nBioconductor:", paste(miss, collapse = ", "), "\n")
  tryCatch(BiocManager::install(miss, lib = lib, ask = FALSE, update = FALSE),
           error = function(e) cat("FAILED bioc:", conditionMessage(e), "\n"))
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

cat("\nBase CRAN/Bioc packages from plan\n")
install_cran_current(c("dplyr", "tidyr", "ggplot2", "patchwork", "pheatmap", "remotes",
                       "msigdbr", "scTenifoldKnk"))
install_bioc(c("SingleCellExperiment", "scDblFinder", "ComplexHeatmap", "UCell",
               "limma", "edgeR", "fgsea", "clusterProfiler", "org.Hs.eg.db"))

cat("\nCompatibility dependencies\n")
install_archive("transport", "0.15-2")
install_archive("httr2", "1.0.7")
install_archive("randomForest", "4.7-1.1")
install_archive("DiagrammeR", "1.0.11")
install_archive("ggnewscale", "0.4.9")
install_archive("ggridges", "0.5.4")
install_archive("bslib", "0.6.2")
install_archive("shiny", "1.8.1.1")
install_archive("mclust", "6.0.1")
install_archive("signal", "0.7-7")
install_cran_current(c("leiden", "miniUI", "plotly", "spatstat.core", "sass", "cli"))

cat("\nSpecial packages from plan\n")
install_github("copykat", "navinlabcode/copykat")
install_github("CaSpER", "akdess/CaSpER")
install_github("liana", "saezlab/liana")
install_github("nichenetr", "saeyslab/nichenetr")
install_github("CellChat", "sqjin/CellChat")

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
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_plan_tools_R_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
