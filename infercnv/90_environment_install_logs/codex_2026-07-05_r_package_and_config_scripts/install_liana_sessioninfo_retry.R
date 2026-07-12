options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "libcurl")
options(timeout = 1200)
Sys.setenv(TMP = "D:/TEMP", TEMP = "D:/TEMP", TMPDIR = "D:/TEMP")
lib <- "D:/Documents/R/win-library/4.0"
.libPaths(c(lib, .libPaths()))

log_file <- "D:/OC_spatiogenomics/infercnv/install_liana_sessioninfo_retry_log.txt"
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("Before sessioninfo:", as.character(packageVersion("sessioninfo")), "\n")
tryCatch({
  remotes::install_version("sessioninfo", version = "1.2.2", lib = lib,
                           repos = "https://cran.rstudio.com",
                           upgrade = "never", dependencies = TRUE)
}, error = function(e) cat("sessioninfo failed:", conditionMessage(e), "\n"))
cat("After sessioninfo:", as.character(packageVersion("sessioninfo")), "\n")
print(args(sessioninfo::session_info))

tryCatch({
  remotes::install_github("saezlab/OmnipathR", lib = lib, upgrade = "never",
                          dependencies = TRUE, build_vignettes = FALSE)
}, error = function(e) cat("OmnipathR failed:", conditionMessage(e), "\n"))

tryCatch({
  remotes::install_github("saezlab/liana", lib = lib, upgrade = "never",
                          dependencies = TRUE, build_vignettes = FALSE)
}, error = function(e) cat("liana failed:", conditionMessage(e), "\n"))

pkgs <- c("sessioninfo", "OmnipathR", "liana", "copykat", "CaSpER", "nichenetr")
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
write.csv(load_status, "D:/OC_spatiogenomics/infercnv/install_liana_sessioninfo_retry_load_status.csv", row.names = FALSE)

ip <- installed.packages()
status <- data.frame(package = pkgs,
                     installed = pkgs %in% rownames(ip),
                     version = ifelse(pkgs %in% rownames(ip), ip[match(pkgs, rownames(ip)), "Version"], NA),
                     stringsAsFactors = FALSE)
print(status)
write.csv(status, "D:/OC_spatiogenomics/infercnv/install_liana_sessioninfo_retry_status.csv", row.names = FALSE)

sink(type = "message")
sink()
close(log_con)
