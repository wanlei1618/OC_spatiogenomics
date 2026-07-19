# sample_type association for SPP1_TAM and ITGB1_CD44 tumor scores

Input object: `D:/OC_spatiogenomics/infercnv/integrated_oc.RData`

Metadata: `D:/OC_spatiogenomics/infercnv/integrated_oc_plan_analysis/tables/integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv`

## Score definitions

- `SPP1_TAM_score`: mean z-scored RNA log-normalized expression of present TAM genes in myeloid/macrophage cells.
- `ITGB1_CD44_tumor_score`: mean z-scored RNA log-normalized expression of ITGB1 and CD44 in CNV-subclone tumor cells.

## Gene coverage

                    score requested_gene present
1          SPP1_TAM_score           SPP1    TRUE
2          SPP1_TAM_score           CD68    TRUE
3          SPP1_TAM_score          CD163    TRUE
4          SPP1_TAM_score           CD14    TRUE
5          SPP1_TAM_score           LST1    TRUE
6          SPP1_TAM_score         TYROBP    TRUE
7          SPP1_TAM_score           C1QA    TRUE
8          SPP1_TAM_score           C1QB    TRUE
9          SPP1_TAM_score           C1QC    TRUE
10         SPP1_TAM_score           APOE    TRUE
11         SPP1_TAM_score           MRC1    TRUE
12         SPP1_TAM_score           MSR1    TRUE
13         SPP1_TAM_score         FCGR3A    TRUE
14         SPP1_TAM_score          ITGAM    TRUE
15         SPP1_TAM_score          CSF1R    TRUE
16 ITGB1_CD44_tumor_score          ITGB1    TRUE
17 ITGB1_CD44_tumor_score           CD44    TRUE

## sample_type summary

                     analysis      sample_type    n        mean      median
1:     myeloid_SPP1_TAM_score          ascites 4296 -0.02860035  0.03568061
2:     myeloid_SPP1_TAM_score            blood 1276 -0.63649807 -0.60142063
3:     myeloid_SPP1_TAM_score pleural_effusion 4769  0.06538000  0.11162136
4:     myeloid_SPP1_TAM_score            tumor 5651  0.11028869  0.16453233
5: CNV_tumor_ITGB1_CD44_score          ascites  562  0.43412027  0.33613677
6: CNV_tumor_ITGB1_CD44_score            blood   52  0.68767089  0.71450537
7: CNV_tumor_ITGB1_CD44_score pleural_effusion  124  0.39648985 -0.50171919
8: CNV_tumor_ITGB1_CD44_score            tumor 5769 -0.05701148 -0.50171919
          sd        q25        q75
1: 0.5926329 -0.4751790  0.4655655
2: 0.2587302 -0.7750689 -0.4573658
3: 0.4869919 -0.2596936  0.4347144
4: 0.4533840 -0.1265926  0.4252224
5: 0.9529368 -0.5017192  1.1237893
6: 1.2543650 -0.5017192  1.6631904
7: 1.0766472 -0.5017192  1.2559292
8: 0.6721138 -0.5017192  0.3145727

## Kruskal-Wallis tests

                    analysis statistic      p_value n_groups n_cells
1     myeloid_SPP1_TAM_score 2093.5036 0.000000e+00        4   15992
2 CNV_tumor_ITGB1_CD44_score  193.7938 9.248065e-42        4    6507
      p_adj_BH
1 0.000000e+00
2 9.248065e-42

## Pairwise Wilcoxon tests

                      analysis           group1           group2   n1   n2
 1:     myeloid_SPP1_TAM_score          ascites            blood 4296 1276
 2:     myeloid_SPP1_TAM_score          ascites pleural_effusion 4296 4769
 3:     myeloid_SPP1_TAM_score          ascites            tumor 4296 5651
 4:     myeloid_SPP1_TAM_score            blood pleural_effusion 1276 4769
 5:     myeloid_SPP1_TAM_score            blood            tumor 1276 5651
 6:     myeloid_SPP1_TAM_score pleural_effusion            tumor 4769 5651
 7: CNV_tumor_ITGB1_CD44_score          ascites            blood  562   52
 8: CNV_tumor_ITGB1_CD44_score          ascites pleural_effusion  562  124
 9: CNV_tumor_ITGB1_CD44_score          ascites            tumor  562 5769
10: CNV_tumor_ITGB1_CD44_score            blood pleural_effusion   52  124
11: CNV_tumor_ITGB1_CD44_score            blood            tumor   52 5769
12: CNV_tumor_ITGB1_CD44_score pleural_effusion            tumor  124 5769
        median1    median2 delta_median_group2_minus_group1       p_value
 1:  0.03568061 -0.6014206                      -0.63710124 8.965251e-245
 2:  0.03568061  0.1116214                       0.07594076  3.150708e-10
 3:  0.03568061  0.1645323                       0.12885172  5.304794e-24
 4: -0.60142063  0.1116214                       0.71304199  0.000000e+00
 5: -0.60142063  0.1645323                       0.76595296  0.000000e+00
 6:  0.11162136  0.1645323                       0.05291096  7.072923e-06
 7:  0.33613677  0.7145054                       0.37836860  3.394967e-01
 8:  0.33613677 -0.5017192                      -0.83785596  2.820185e-01
 9:  0.33613677 -0.5017192                      -0.83785596  1.161642e-38
10:  0.71450537 -0.5017192                      -1.21622456  1.747786e-01
11:  0.71450537 -0.5017192                      -1.21622456  2.208961e-05
12: -0.50171919 -0.5017192                       0.00000000  7.482266e-05
         p_adj_BH
 1: 1.793050e-244
 2:  3.780849e-10
 3:  7.957190e-24
 4:  0.000000e+00
 5:  0.000000e+00
 6:  7.072923e-06
 7:  3.394967e-01
 8:  3.384222e-01
 9:  6.969851e-38
10:  2.621680e-01
11:  6.626883e-05
12:  1.496453e-04

## Output files

- `tables/cell_level_SPP1_TAM_score_myeloid_by_sample_type.csv`
- `tables/cell_level_ITGB1_CD44_tumor_score_CNV_cells_by_sample_type.csv`
- `tables/sample_type_score_summary.csv`
- `tables/sample_type_score_kruskal_wallis_tests.csv`
- `tables/sample_type_score_pairwise_wilcoxon_tests.csv`
- `figures/SPP1_TAM_score_by_sample_type_myeloid.png`
- `figures/ITGB1_CD44_tumor_score_by_sample_type_CNV_tumor.png`
