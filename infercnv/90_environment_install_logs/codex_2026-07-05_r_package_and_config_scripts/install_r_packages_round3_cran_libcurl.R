options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "libcurl")
options(timeout = 1200)
lib <- "D:/Documents/R/win-library/4.0"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_CopyKAT_CaSpER_LIANA_NicheNet_round3_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("R version:", R.version.string, "\n")
cat("Library paths:\n")
print(.libPaths())
cat("CRAN connectivity rows:\n")
ap <- available.packages()
cat(nrow(ap), "\n")

install_cran <- function(pkgs) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing) == 0) {
    cat("All already installed:", paste(pkgs, collapse = ", "), "\n")
    return(invisible(TRUE))
  }
  cat("\nInstalling CRAN source packages:", paste(missing, collapse = ", "), "\n")
  tryCatch({
    install.packages(missing, lib = lib, dependencies = TRUE, type = "source", Ncpus = 2)
  }, error = function(e) {
    cat("CRAN install failed:", conditionMessage(e), "\n")
  })
}

copykat_deps <- c("parallelDist", "dlm", "mixtools", "MCMCpack", "transport")
liana_deps <- c("httr2", "later", "logger")
nichenetr_deps <- c("fdrtool", "ROCR", "randomForest", "DiagrammeR",
                    "mlrMBO", "parallelMap", "emoa", "DiceKriging",
                    "ggnewscale")
seurat_deps <- c("fitdistrplus", "ggridges", "ica", "leiden", "miniUI", "plotly",
                 "RANN", "RcppAnnoy", "scattermore", "sctransform", "shiny",
                 "spatstat.core", "spatstat.geom", "uwot", "RcppProgress")
casper_deps <- c("signal", "pheatmap", "reshape", "mclust", "intergraph", "ape")

install_cran(unique(c(copykat_deps, liana_deps, nichenetr_deps, seurat_deps, casper_deps)))

cat("\nInstalling Bioconductor deps for CaSpER and NicheNet if available...\n")
tryCatch({
  BiocManager::install(c("GOstats", "graph", "RBGL", "GenomicRanges", "GenomeInfoDb",
                         "IRanges", "S4Vectors", "AnnotationDbi"),
                       lib = lib, ask = FALSE, update = FALSE)
}, error = function(e) cat("Bioc deps failed:", conditionMessage(e), "\n"))

cat("\nRetry Seurat 4.0.6...\n")
if (!("Seurat" %in% rownames(installed.packages()))) {
  tryCatch({
    remotes::install_version("Seurat", version = "4.0.6", lib = lib,
                             repos = "https://cran.rstudio.com", upgrade = "never",
                             dependencies = TRUE)
  }, error = function(e) cat("Seurat 4.0.6 failed:", conditionMessage(e), "\n"))
}

install_github_pkg <- function(pkg, repo) {
  cat("\n==== GitHub install:", pkg, repo, "====\n")
  tryCatch({
    remotes::install_github(repo, lib = lib, upgrade = "never",
                            dependencies = TRUE, build_vignettes = FALSE)
  }, error = function(e) {
    cat("FAILED:", conditionMessage(e), "\n")
  })
  cat(pkg, "installed?", pkg %in% rownames(installed.packages()), "\n")
  if (pkg %in% rownames(installed.packages())) {
    cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
  }
}

install_github_pkg("copykat", "navinlabcode/copykat")
install_github_pkg("nichenetr", "saeyslab/nichenetr")
install_github_pkg("liana", "saezlab/liana")

cat("\nCaSpER retry...\n")
tryCatch({
  remotes::install_github("akdess/CaSpER", lib = lib, upgrade = "never",
                          dependencies = TRUE, build_vignettes = FALSE)
}, error = function(e) cat("GitHub CaSpER failed:", conditionMessage(e), "\n"))

cat("\n==== Final status ====\n")
pkgs <- c("copykat", "CaSpER", "liana", "nichenetr", "Seurat", "SeuratObject",
          copykat_deps, liana_deps, nichenetr_deps, seurat_deps, casper_deps, "GOstats")
ip <- installed.packages()
status <- data.frame(
  package = pkgs,
  installed = pkgs %in% rownames(ip),
  version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
  stringsAsFactors = FALSE
)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_CopyKAT_CaSpER_LIANA_NicheNet_status_round3.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
