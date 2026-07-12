# Reproduction status

Generated: 2026-07-12

## Environment

- Git: 2.40.0.windows.1
- Python: 3.7.4
- Rscript: D:/R/R-4.0.3/bin/Rscript.exe, R 4.0.3
- R libraries used from D:/OC_spatiogenomics/spatial_data/R_library, D:/Documents/R/win-library/4.0, and D:/R/R-4.0.3/library

## Step status

| Step | Status | Notes |
|---|---|---|
| 01 download dry-run | success | Listed only GSM6177614, GSM6177617, and GSE189843 archive. |
| 01 download | success | GEO files downloaded to D:/OC_spatiogenomics/spatial_data/raw. |
| 02 build objects | partial success | Built GSM6177614 and GSM6177617. GSE189843 Seurat object creation failed because the installed Matrix/Seurat binary combination raises `...names is not a BUILTIN function`; expression-only analysis was recovered in step 10 by direct MTX parsing without coordinates. |
| 03 score/statistics | success | Coordinate-aware results generated only for GSM6177614 and GSM6177617. |
| 04 reference mapping | success with fallback | Reference file was an RDS object despite .RData extension. Seurat anchor transfer failed under the local Matrix/Seurat binary combination, so predictions are marked `fallback_score_based` and must not be interpreted as direct CNV measurements. |
| 05 strict audit | passed | No excluded samples leaked; no GSE189843 coordinate-neighborhood results. |
| 06 QC sensitivity | success | Wrote grid and summary tables plus heatmap. |
| 07 autocorrelation/multiscale | success | Moran's I, Geary's C, and kNN enrichment for coordinate-aware samples. |
| 08 directional niche statistics | success | Source-target and negative-control tables generated. |
| 09 mapping stability | success | Low-confidence predictions retained as Uncertain. |
| 10 GSE189843 response | success | 12 samples analyzed at sample level from raw MTX; no coordinate analysis. |
| 11 patient-level meta | success | Spatial and expression evidence layers kept separate. |
| 12 figures | success | PDF and SVG figures generated. |
| 13 report | success | Final report, manifest, and session info generated. |

## Sample audit

- Included coordinate-aware GSE203612 ovarian samples: GSM6177614, GSM6177617.
- Explicitly excluded: GSM6177618 (PDAC).
- Included expression-only GSE189843 samples: GSM5708485-GSM5708496, with GSM5708485-GSM5708490 Excellent and GSM5708491-GSM5708496 Poor.
- Built Seurat objects: 2 coordinate-aware samples.
- GSE203612 raw spots: GSM6177614 = 1762, GSM6177617 = 1661.
- GSE203612 QC-retained spots: GSM6177614 = 1760, GSM6177617 = 1204.

## Key SHA-256

| File | SHA-256 |
|---|---|
| logs/download_audit.json | 6D653DDEA42678CB33A4229207BF4B7B02CF1A603C4A70018C647CA631C4A6CA |
| processed/spatial_objects_curated.rds | 4D57DFE1414333C86FD18E0F9202921CDFA57D5FD752918D1B64CE93BD11F44F |
| results/spatial_curated/spatial_spot_scores_curated.csv.gz | F8A53C1397954EB3A967DE30610304FE9DBFCC91A1B8FB125647FF2CDC11C9E8 |
| reports/spatial_validation_final.md | C5B19EFD03CB92F1EF3C64410F71F4F47376FDFACE79A472F32D92E18930DD63 |

## Boundaries

- SPP1-CD44 is a candidate ligand-receptor axis.
- SPP1-ITGB1 is reported as an SPP1-associated ITGB1-positive adhesion/integrin program, not a proven direct ligand-receptor interaction.
- Spot-level P values describe within-sample patterns only.
- Two coordinate-aware ovarian samples are not sufficient to claim a universal mechanism.
