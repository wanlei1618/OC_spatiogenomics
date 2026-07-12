options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "libcurl")
options(timeout = 1200)
lib <- "D:/Documents/R/win-library/4.0"
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_CopyKAT_CaSpER_LIANA_NicheNet_round4_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("R version:", R.version.string, "\n")
cat("Library paths:\n")
print(.libPaths())

pkgs <- c("transport", "httr2", "randomForest", "DiagrammeR", "ggnewscale",
          "ggridges", "leiden", "miniUI", "plotly", "shiny", "spatstat.core",
          "signal", "mclust", "bslib", "sass")

missing <- setdiff(pkgs, rownames(installed.packages()))
cat("Missing targeted deps:", paste(missing, collapse = ", "), "\n")
if (length(missing)) {
  tryCatch({
    install.packages(missing, lib = lib, dependencies = TRUE, type = "binary")
  }, error = function(e) cat("binary install failed:", conditionMessage(e), "\n"))
}

install_github_pkg <- function(pkg, repo) {
  cat("\n==== GitHub install:", pkg, repo, "====\n")
  tryCatch({
    remotes::install_github(repo, lib = lib, upgrade = "never",
                            dependencies = TRUE, build_vignettes = FALSE)
  }, error = function(e) cat("FAILED:", conditionMessage(e), "\n"))
  cat(pkg, "installed?", pkg %in% rownames(installed.packages()), "\n")
  if (pkg %in% rownames(installed.packages())) {
    cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
  }
}

install_github_pkg("copykat", "navinlabcode/copykat")
install_github_pkg("nichenetr", "saeyslab/nichenetr")
install_github_pkg("liana", "saezlab/liana")
install_github_pkg("CaSpER", "akdess/CaSpER")

cat("\n==== Load tests ====\n")
target <- c("copykat", "CaSpER", "liana", "nichenetr")
load_status <- lapply(target, function(pkg) {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  }, error = function(e) {
    cat(pkg, "load failed:", conditionMessage(e), "\n")
    FALSE
  })
  data.frame(package = pkg, load_ok = ok, stringsAsFactors = FALSE)
})
load_status <- do.call(rbind, load_status)
print(load_status)
write.csv(load_status, "D:/OC_spatiogenomics/infercnv/install_CopyKAT_CaSpER_LIANA_NicheNet_load_status.csv", row.names = FALSE)

cat("\n==== Final status ====\n")
all_pkgs <- unique(c(target, "Seurat", "SeuratObject", pkgs))
ip <- installed.packages()
status <- data.frame(
  package = all_pkgs,
  installed = all_pkgs %in% rownames(ip),
  version = ifelse(all_pkgs %in% rownames(ip), ip[match(all_pkgs, rownames(ip)), "Version"], NA),
  stringsAsFactors = FALSE
)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_CopyKAT_CaSpER_LIANA_NicheNet_status_round4.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
