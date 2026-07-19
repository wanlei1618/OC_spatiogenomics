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
options(repos = repos)

log_file <- "D:/OC_spatiogenomics/infercnv/install_remaining_r_plan_packages_round2_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("Remaining package install round2\n")

install_current <- function(pkgs) {
  cat("\nInstalling current:", paste(pkgs, collapse = ", "), "\n")
  tryCatch(install.packages(pkgs, lib = lib, repos = repos, type = "source",
                            dependencies = c("Depends", "Imports", "LinkingTo"), Ncpus = 2),
           error = function(e) cat("FAILED current:", conditionMessage(e), "\n"))
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
  cat("\nGitHub:", pkg, repo, "\n")
  tryCatch(remotes::install_github(repo, lib = lib, upgrade = "never",
                                  dependencies = c("Depends", "Imports", "LinkingTo"),
                                  build_vignettes = FALSE),
           error = function(e) cat("FAILED github", pkg, ":", conditionMessage(e), "\n"))
}

install_current(c("rlang", "cli", "lifecycle", "vctrs", "pillar", "tibble"))
install_archive("dplyr", "1.1.4")
install_current(c("edgeR", "scDblFinder"))
install_github("UCell", "carmonalab/UCell@v1.0")
install_github("liana", "saezlab/liana")

cat("\nFinal status\n")
pkgs <- c("rlang", "vctrs", "pillar", "tibble", "dplyr",
          "SingleCellExperiment", "scDblFinder", "UCell", "edgeR", "liana")
ip <- installed.packages()
status <- data.frame(
  package = pkgs,
  installed = pkgs %in% rownames(ip),
  version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
  stringsAsFactors = FALSE
)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_remaining_r_plan_packages_round2_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
