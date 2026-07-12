options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 1200)

lib <- "D:/Documents/R/win-library/4.0"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_CopyKAT_CaSpER_LIANA_NicheNet_log.txt"
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("R version:", R.version.string, "\n")
cat("Library paths:\n")
print(.libPaths())

install_if_missing <- function(pkg, installer) {
  cat("\n====", pkg, "====\n")
  if (pkg %in% rownames(installed.packages())) {
    cat(pkg, "already installed:", as.character(packageVersion(pkg)), "\n")
    return(TRUE)
  }
  ok <- tryCatch({
    installer()
    TRUE
  }, error = function(e) {
    cat("FAILED:", conditionMessage(e), "\n")
    FALSE
  })
  if (pkg %in% rownames(installed.packages())) {
    cat(pkg, "installed:", as.character(packageVersion(pkg)), "\n")
    TRUE
  } else {
    cat(pkg, "not installed after attempt.\n")
    ok && FALSE
  }
}

install_if_missing("BiocManager", function() {
  install.packages("BiocManager", lib = lib, dependencies = TRUE)
})

install_if_missing("remotes", function() {
  install.packages("remotes", lib = lib, dependencies = TRUE)
})

cat("\nBioconductor version info:\n")
print(tryCatch(BiocManager::version(), error = function(e) conditionMessage(e)))

results <- list()

results$copykat <- install_if_missing("copykat", function() {
  remotes::install_github("navinlabcode/copykat", lib = lib, upgrade = "never",
                          dependencies = TRUE, build_vignettes = FALSE)
})

results$CaSpER <- install_if_missing("CaSpER", function() {
  BiocManager::install("CaSpER", lib = lib, ask = FALSE, update = FALSE)
})

results$liana <- install_if_missing("liana", function() {
  install.packages("liana", lib = lib, dependencies = TRUE)
  if (!("liana" %in% rownames(installed.packages()))) {
    remotes::install_github("saezlab/liana", lib = lib, upgrade = "never",
                            dependencies = TRUE, build_vignettes = FALSE)
  }
})

results$nichenetr <- install_if_missing("nichenetr", function() {
  remotes::install_github("saeyslab/nichenetr", lib = lib, upgrade = "never",
                          dependencies = TRUE, build_vignettes = FALSE)
})

cat("\n==== Final status ====\n")
pkgs <- c("copykat", "CaSpER", "liana", "nichenetr")
ip <- installed.packages()
status <- data.frame(
  package = pkgs,
  installed = pkgs %in% rownames(ip),
  version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
  stringsAsFactors = FALSE
)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_CopyKAT_CaSpER_LIANA_NicheNet_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
