library(survival)
bulk_meta <- read.csv('bulk_meta_with_scores.csv')
fit <- coxph(Surv(OS_time, OS_event) ~ SPP1_TAM_score * ITGB1_CD44_tumor_score +
               KRAS_Hypoxia_score + macrophage_fraction + tumor_purity +
               stage + grade + residual_disease, data = bulk_meta)
summary(fit)
