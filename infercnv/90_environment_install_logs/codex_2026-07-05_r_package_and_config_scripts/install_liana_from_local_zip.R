options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "libcurl")
options(timeout = 1200)
Sys.setenv(TMP = "D:/TEMP", TEMP = "D:/TEMP", TMPDIR = "D:/TEMP")
lib <- "D:/Documents/R/win-library/4.0"
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_liana_from_local_zip_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

install_local <- function(pkg, path) {
  cat("\nLocal install:", pkg, path, "\n")
  tryCatch({
    remotes::install_local(path, lib = lib, upgrade = "never",
                           dependencies = TRUE, build_vignettes = FALSE,
                           INSTALL_opts = "--no-multiarch")
  }, error = function(e) cat(pkg, "failed:", conditionMessage(e), "\n"))
  cat(pkg, "installed?", pkg %in% rownames(installed.packages()), "\n")
  if (pkg %in% rownames(installed.packages())) cat(pkg, as.character(packageVersion(pkg)), "\n")
}

install_local("OmnipathR", "D:/TEMP/OmnipathR-master.zip")
install_local("liana", "D:/TEMP/liana-master.zip")

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
write.csv(load_status, "D:/OC_spatiogenomics/infercnv/install_liana_from_local_zip_load_status.csv", row.names = FALSE)

ip <- installed.packages()
status <- data.frame(package = pkgs,
                     installed = pkgs %in% rownames(ip),
                     version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
                     stringsAsFactors = FALSE)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_liana_from_local_zip_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
