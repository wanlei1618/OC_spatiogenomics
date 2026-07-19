options(download.file.method = "libcurl")
options(timeout = 1200)
Sys.setenv(TMP = "D:/TEMP", TEMP = "D:/TEMP", TMPDIR = "D:/TEMP")
lib <- "D:/Documents/R/win-library/4.0"
.libPaths(c(lib, .libPaths()))

repos <- c(
  BioCsoft = "https://bioconductor.org/packages/3.12/bioc",
  BioCann = "https://bioconductor.org/packages/3.12/data/annotation",
  BioCexp = "https://bioconductor.org/packages/3.12/data/experiment",
  CRAN = "https://cran.rstudio.com"
)

log_file <- "D:/OC_spatiogenomics/infercnv/install_remaining_r_plan_packages_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("Install remaining R packages\n")
cat("R:", R.version.string, "\n")

install_from_repos <- function(pkgs) {
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (!length(miss)) return(invisible(TRUE))
  cat("\nInstalling:", paste(miss, collapse = ", "), "\n")
  tryCatch(install.packages(miss, lib = lib, repos = repos, type = "source",
                            dependencies = c("Depends", "Imports", "LinkingTo"), Ncpus = 2),
           error = function(e) cat("FAILED:", conditionMessage(e), "\n"))
}

install_archive <- function(pkg, version) {
  cat("\nArchive:", pkg, version, "\n")
  tryCatch(remotes::install_version(pkg, version = version, lib = lib,
                                    repos = "https://cran.rstudio.com",
                                    upgrade = "never",
                                    dependencies = c("Depends", "Imports", "LinkingTo")),
           error = function(e) cat("FAILED archive", pkg, ":", conditionMessage(e), "\n"))
}

install_github <- function(pkg, repo) {
  if (pkg %in% rownames(installed.packages())) {
    cat(pkg, "already installed", as.character(packageVersion(pkg)), "\n")
    return(invisible(TRUE))
  }
  cat("\nGitHub:", pkg, repo, "\n")
  tryCatch(remotes::install_github(repo, lib = lib, upgrade = "never",
                                  dependencies = c("Depends", "Imports", "LinkingTo"),
                                  build_vignettes = FALSE),
           error = function(e) cat("FAILED github", pkg, ":", conditionMessage(e), "\n"))
}

install_archive("dplyr", "1.1.4")
install_from_repos(c("SummarizedExperiment", "SingleCellExperiment", "edgeR", "scDblFinder"))
install_github("UCell", "carmonalab/UCell")
install_github("liana", "saezlab/liana")

cat("\nFinal status\n")
pkgs <- c("dplyr", "SingleCellExperiment", "scDblFinder", "UCell", "edgeR", "liana")
ip <- installed.packages()
status <- data.frame(
  package = pkgs,
  installed = pkgs %in% rownames(ip),
  version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
  stringsAsFactors = FALSE
)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_remaining_r_plan_packages_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
