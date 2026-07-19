options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "libcurl")
options(timeout = 1200)
Sys.setenv(TMP = "D:/TEMP", TEMP = "D:/TEMP", TMPDIR = "D:/TEMP")
lib <- "D:/Documents/R/win-library/4.0"
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_liana_final_retry_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("R:", R.version.string, "\n")
cat("Before dplyr:", if ("dplyr" %in% rownames(installed.packages())) as.character(packageVersion("dplyr")) else "not installed", "\n")

cat("\nForce installing dplyr 1.1.4...\n")
tryCatch({
  remotes::install_version("dplyr", version = "1.1.4", lib = lib,
                           repos = "https://cran.rstudio.com",
                           upgrade = "never", dependencies = TRUE)
}, error = function(e) cat("dplyr install failed:", conditionMessage(e), "\n"))

cat("After dplyr:", if ("dplyr" %in% rownames(installed.packages())) as.character(packageVersion("dplyr")) else "not installed", "\n")

cat("\nInstalling OmnipathR...\n")
tryCatch({
  remotes::install_github("saezlab/OmnipathR", lib = lib, upgrade = "never",
                          dependencies = TRUE, build_vignettes = FALSE)
}, error = function(e) cat("OmnipathR failed:", conditionMessage(e), "\n"))

cat("\nInstalling liana...\n")
tryCatch({
  remotes::install_github("saezlab/liana", lib = lib, upgrade = "never",
                          dependencies = TRUE, build_vignettes = FALSE)
}, error = function(e) cat("liana failed:", conditionMessage(e), "\n"))

cat("\nLoad tests:\n")
pkgs <- c("dplyr", "OmnipathR", "liana", "copykat", "CaSpER", "nichenetr")
load_status <- lapply(pkgs, function(pkg) {
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
write.csv(load_status, "D:/OC_spatiogenomics/infercnv/install_liana_final_retry_load_status.csv", row.names = FALSE)

ip <- installed.packages()
status <- data.frame(
  package = pkgs,
  installed = pkgs %in% rownames(ip),
  version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
  stringsAsFactors = FALSE
)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_liana_final_retry_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
