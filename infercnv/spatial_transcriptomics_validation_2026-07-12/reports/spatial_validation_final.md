# Spatial transcriptomics validation report

Generated: 2026-07-12 11:36:15

## Data sources and sample audit

- GSE203612 coordinate-aware ovarian Visium analysis is restricted to GSM6177614 and GSM6177617.
- GSM6177618 is excluded from ovarian results because it is PDAC.
- GSE189843 contributes 12 pretreatment HGSC expression-level samples only; no coordinate-neighborhood claim is made for this series.

## Scoring definitions

- Scores are read from `config/spatial_config.yml`, including SPP1_myeloid, target_subclone_02_04, KRAS_hypoxia, target_core_without_receptors, and SPP1_myeloid_without_SPP1.
- SPP1-CD44 is treated as a candidate ligand-receptor axis. SPP1-ITGB1 is reported as an SPP1-associated ITGB1-positive adhesion/integrin program, not as a proven direct ligand-receptor interaction.

## Statistical boundaries

- Spot-level spatial tests describe within-sample patterns only.
- Patient/sample-level summaries are the unit for cross-sample statements.
- Two coordinate-aware ovarian samples are insufficient to prove a universal mechanism.
- Low-confidence reference transfer spots are retained as uncertain rather than forced into a unique CNV label.

## Evidence layers

                      evidence_layer n_samples
1: coordinate-aware spatial evidence         2
2:      expression-level replication        12
                                      scope
1:    GSE203612 ovarian Visium samples only
2: GSE189843 pretreatment HGSC samples only
                                                            conclusion_boundary
1: Descriptive because only two coordinate-aware ovarian samples are available.
2:     Patient-level expression replication; no coordinate neighborhood claims.
   evidence_grade
1:        limited
2:       moderate

## QC sensitivity

    sample_id n_parameter_sets median_retained_spots positive_direction_rate
1: GSM6177614               81                  1760                       0
2: GSM6177617               81                  1204                       0
   unstable_parameter_rate
1:                       0
2:                       0

## GSE189843 response analysis

                                       metric excellent_median   poor_median
 1:                        SPP1_myeloid_score     7.725050e-18  1.897210e-18
 2:                                      CD44     5.848566e-01  2.572328e-01
 3:                                     ITGB1     9.968805e-01  9.785534e-01
 4:                    SPP1_CD44_expr_product     5.263687e-01  3.795269e-01
 5:                   SPP1_ITGB1_expr_product     8.480665e-01  1.352910e+00
 6:               target_subclone_02_04_score     2.152125e-17  3.806918e-18
 7:                        KRAS_hypoxia_score    -2.662857e-18  1.617324e-18
 8:          SPP1_myeloid_score_high_fraction     2.727273e-01  2.727273e-01
 9: target_subclone_02_04_score_high_fraction     2.727273e-01  2.727273e-01
10:          KRAS_hypoxia_score_high_fraction     2.727273e-01  2.727273e-01
11:                      score_score_spearman    -2.181818e-01 -9.090909e-03
    difference_poor_minus_excellent wilcoxon_p cliffs_delta_poor_vs_excellent
 1:                   -5.827840e-18  0.8101812                     -0.1111111
 2:                   -3.276239e-01  0.1282053                     -0.5555556
 3:                   -1.832712e-02  1.0000000                      0.0000000
 4:                   -1.468419e-01  0.5751735                     -0.2222222
 5:                    5.048435e-01  0.2979531                      0.3888889
 6:                   -1.771433e-17  0.2979531                     -0.3888889
 7:                    4.280181e-18  0.8101812                      0.1111111
 8:                    0.000000e+00         NA                      0.0000000
 9:                    0.000000e+00         NA                      0.0000000
10:                    0.000000e+00         NA                             NA
11:                    2.090909e-01  0.8101812                      0.1111111
    bootstrap_ci_low bootstrap_ci_high n_excellent n_poor       fdr
 1:    -3.209002e-17      1.401330e-17           6      6 0.9259214
 2:    -6.687919e-01     -7.629347e-02           6      6 0.7945415
 3:    -3.747131e-01      3.489147e-01           6      6 1.0000000
 4:    -1.032583e+00      1.705057e-01           6      6 0.9259214
 5:    -2.860990e-01      1.069558e+00           6      6 0.7945415
 6:    -5.653862e-17      1.411501e-17           6      6 0.7945415
 7:    -1.056303e-17      2.998626e-17           6      6 0.9259214
 8:     0.000000e+00      0.000000e+00           6      6        NA
 9:     0.000000e+00      0.000000e+00           6      6        NA
10:     0.000000e+00      0.000000e+00           6      6        NA
11:    -4.106818e-01      4.804167e-01           6      6 0.9259214

## Output manifest

See `reports/spatial_validation_result_manifest.csv` for generated file sizes and SHA-256 hashes when the digest package is available.
