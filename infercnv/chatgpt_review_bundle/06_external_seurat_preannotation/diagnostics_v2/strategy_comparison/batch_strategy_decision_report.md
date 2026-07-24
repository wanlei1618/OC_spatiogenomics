# Batch strategy decision report

Strategies were compared within each provisional broad lineage. Lower sample dominance and sample silhouette were balanced against cluster/lineage preservation; mixing alone was never the selection rule.

Epithelial-like cells always retain the uncorrected result as the primary patient-specific tumor-state view. Corrected epithelial results are sensitivity views only.

- iLISI: SKIPPED_PACKAGE_UNAVAILABLE
- kBET: SKIPPED_PACKAGE_UNAVAILABLE
- RNA marker coherence is finalized in step 07; this step uses embedding coherence and broad-marker lineage purity as an explicitly labeled proxy.

## Recommendations

- GSE154600 / B_Plasma_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.8558.
- GSE154600 / Cycling_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.7617.
- GSE154600 / Endothelial_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.9615.
- GSE154600 / Epithelial_like: `A_uncorrected` - Preserve the uncorrected patient-specific epithelial structure as the primary result. Harmony/RPCA remain parallel shared-state sensitivity outputs when completed.
- GSE154600 / Fibroblast_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.196.
- GSE154600 / Myeloid_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.3755.
- GSE154600 / T_NK_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.4825.
- GSE158722 / B_Plasma_like: `B_harmony_patient` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.5241.
- GSE158722 / Cycling_like: `B_harmony_sample` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.3964.
- GSE158722 / Epithelial_like: `A_uncorrected` - Preserve the uncorrected patient-specific epithelial structure as the primary result. Harmony/RPCA remain parallel shared-state sensitivity outputs when completed.
- GSE158722 / Fibroblast_like: `A_uncorrected` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.8188.
- GSE158722 / Myeloid_like: `A_uncorrected` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=-0.008601.
- GSE158722 / T_NK_like: `B_harmony_patient` - Highest balanced score among completed strategies; score combines sample mixing, cluster preservation, broad-lineage purity, and QC association. Score=0.1891.
