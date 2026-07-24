options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

local_root <- "D:/OC_spatiogenomics/spatial_data/spatial_analysis_correction_v2"
dir.create(local_root, recursive = TRUE, showWarnings = FALSE)
cnv_root <- paste0(
  "D:/OC_spatiogenomics/infercnv/external_seurat_preannotation/"
)
calibrated_file <- file.path(
  cnv_root, "research_spatial_transition",
  "GSE154600_calibrated_cnv_by_cell.csv.gz"
)
counts <- readRDS(file.path(
  cnv_root, "research_validation_independent_cnv",
  "GSE154600_complete_final_epithelial_counts.rds"
))

if (file.exists(calibrated_file)) {
  cnv <- fread(calibrated_file)
  selected <- cnv[
    integrated_calibrated_cnv_evidence == "CALIBRATED_DUAL_METHOD_SUPPORT"
  ]
  evidence_tier <- "CALIBRATED_DUAL_METHOD_SUPPORT"
} else {
  consensus <- fread(file.path(
    cnv_root, "research_validation_independent_cnv",
    "GSE154600_copykat_infercnv_consensus_by_cell.csv.gz"
  ))
  status_col <- intersect(
    c("integrated_cnv_evidence", "consensus_class"), names(consensus)
  )[[1L]]
  selected <- consensus[
    get(status_col) %chin% c("DUAL_METHOD_MALIGNANT_SUPPORT",
                             "DUAL_METHOD_SUPPORT")
  ]
  cnv <- consensus
  evidence_tier <- "PROVISIONAL_DUAL_METHOD_SUPPORT"
}
stopifnot(nrow(selected) > 0L)

patient_col <- intersect(c("patient_id", "sample_id"), names(cnv))[[1L]]
lfc <- rbindlist(lapply(unique(cnv[[patient_col]]), function(pid) {
  all_ids <- intersect(
    cnv[get(patient_col) == pid, cell_id], colnames(counts)
  )
  sel_ids <- intersect(
    selected[get(patient_col) == pid, cell_id], colnames(counts)
  )
  other_ids <- setdiff(all_ids, sel_ids)
  if (length(sel_ids) < 20L || length(other_ids) < 20L) return(NULL)
  s1 <- Matrix::rowSums(counts[, sel_ids, drop = FALSE])
  s0 <- Matrix::rowSums(counts[, other_ids, drop = FALSE])
  data.table(
    patient_id = pid, gene = rownames(counts),
    log2_cpm_difference =
      log2(s1 / pmax(sum(s1), 1) * 1e6 + 1) -
      log2(s0 / pmax(sum(s0), 1) * 1e6 + 1),
    n_dual = length(sel_ids), n_other = length(other_ids)
  )
}))

forbidden_exact <- c(
  "SPP1", "CD44", "ITGB1", "ITGA4", "ITGA5", "ITGAV", "ITGA8", "ITGA9",
  "MKI67", "TOP2A", "UBE2C", "CENPF", "TYMS", "PCNA", "STMN1",
  "TUBA1B", "HMGB2", "CDK1", "CCNB1", "CCNB2"
)
audit <- lfc[, .(
  n_patients_positive = sum(log2_cpm_difference > 0),
  n_patients_evaluated = .N,
  mean_log2_cpm_difference = mean(log2_cpm_difference),
  min_log2_cpm_difference = min(log2_cpm_difference)
), by = gene]
audit[, excluded := gene %chin% forbidden_exact |
  grepl("^(MT-|RPL|RPS|HLA-|IG[HKL]|PTPRC$|CD74$|LST1$|TYROBP$|FCER1G$|C1Q)",
        gene)]
audit[, selection_rule_pass := !excluded &
        n_patients_positive >= 3L & mean_log2_cpm_difference > .25]
setorder(audit, -selection_rule_pass, -n_patients_positive,
         -mean_log2_cpm_difference)
malignant <- head(audit[selection_rule_pass == TRUE, gene], 30L)
if (length(malignant) < 5L) {
  stop("Updated CNV-supported receiver signature has fewer than five genes")
}

signature_root <- paste0(
  "D:/OC_spatiogenomics/infercnv/",
  "05_sample_type_lr_niche_external_validation/",
  "sample_type_LR_niche_analysis/tables"
)
read_signature <- function(name) {
  x <- fread(file.path(signature_root, name))
  unique(as.character(x[[1L]]))
}
sub02 <- read_signature("signature_Subclone02_like.csv")
sub04 <- read_signature("signature_Subclone04_like.csv")
common <- read_signature("signature_Subclone02_04_common.csv")
kras_hypoxia <- c(
  "ATF3", "DDIT4", "VEGFA", "CA9", "EGLN3", "BNIP3", "NDRG1",
  "ADM", "LDHA", "JUN", "FOS"
)
fixed <- list(
  macrophage_identity =
    c("LST1", "TYROBP", "FCER1G", "CTSD", "C1QA", "C1QB", "C1QC"),
  SPP1_macrophage_program =
    c("SPP1", "APOC1", "GPNMB", "TREM2", "LPL", "CTSD"),
  C1QC_macrophage_control =
    c("C1QA", "C1QB", "C1QC", "APOE", "MRC1", "SELENOP"),
  epithelial_identity =
    c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "PAX8", "MSLN", "MUC16"),
  CNV_supported_malignant_epithelial_signature = malignant,
  Subclone02_like = sub02,
  Subclone04_like = sub04,
  Subclone02_04_common = common,
  Subclone02_04_KRAS_hypoxia =
    unique(c(common, kras_hypoxia))
)
signature <- rbindlist(lapply(names(fixed), function(nm) {
  data.table(
    signature_name = nm, gene = unique(fixed[[nm]]),
    evidence_tier = if (nm == "CNV_supported_malignant_epithelial_signature") {
      evidence_tier
    } else {
      "PREDEFINED_OR_FROZEN_SIGNATURE"
    }
  )
}))
fwrite(signature, file.path(local_root, "updated_spatial_signatures.csv"),
       na = "NA")
audit[, selected_in_top30 := gene %chin% malignant]
audit[, malignant_signature_evidence_tier := evidence_tier]
fwrite(audit, file.path(local_root, "updated_spatial_signature_audit.csv"),
       na = "NA")
message("Updated sender/receiver signatures built from patient pseudobulk")
