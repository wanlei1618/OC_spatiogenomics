options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "libcurl")
options(timeout = 1200)
Sys.setenv(TMP = "D:/TEMP", TEMP = "D:/TEMP", TMPDIR = "D:/TEMP")
lib <- "D:/Documents/R/win-library/4.0"
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_liana_no_multiarch_retry_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

install_url_nm <- function(pkg, url) {
  cat("\nInstall URL no-multiarch:", pkg, "\n", url, "\n")
  tryCatch({
    remotes::install_url(url, lib = lib, upgrade = "never", dependencies = TRUE,
                         build_vignettes = FALSE, INSTALL_opts = "--no-multiarch")
  }, error = function(e) cat(pkg, "failed:", conditionMessage(e), "\n"))
  cat(pkg, "installed?", pkg %in% rownames(installed.packages()), "\n")
  if (pkg %in% rownames(installed.packages())) cat(pkg, as.character(packageVersion(pkg)), "\n")
}

install_url_nm("OmnipathR", "https://github.com/saezlab/OmnipathR/archive/refs/heads/master.zip")
install_url_nm("liana", "https://github.com/saezlab/liana/archive/refs/heads/master.zip")

pkgs <- c("OmnipathR", "liana", "copykat", "CaSpER", "nichenetr", "dplyr", "sessioninfo", "purrr")
load_status <- do.call(rbind, lapply(pkgs, function(pkg) {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  }, error = function(e) {
    cat(pkg, "load failed:", conditionMessage(e), "\n")
    FALSE
  })
  data.frame(package = pkg, load_ok = ok, stringsAsFactors = FALSE)
}))
print(load_status)
write.csv(load_status, "D:/OC_spatiogenomics/infercnv/install_liana_no_multiarch_retry_load_status.csv", row.names = FALSE)

ip <- installed.packages()
status <- data.frame(package = pkgs,
                     installed = pkgs %in% rownames(ip),
                     version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
                     stringsAsFactors = FALSE)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_liana_no_multiarch_retry_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
