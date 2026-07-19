# Spatial transcriptomics data upgrade

This directory corrects and extends the ovarian-cancer spatial-transcriptomics layer used by the `SPP1-CD44/ITGB1` niche analysis.

## What was corrected

### GSE203612

Only two Visium samples in GSE203612 are ovarian carcinoma:

- `GSM6177614` — `NYU_OVCA1_Vis`
- `GSM6177617` — `NYU_OVCA3_Vis`

`GSM6177618` is titled `NYU_PDAC1_Vis` and its GEO source is *primary pancreatic ductal adenocarcinoma*. GEO contains internally inconsistent characteristic fields that say “ovarian carcinoma/ovary”; the sample title and source are more specific and it is excluded from every ovarian analysis.

Official records:

- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203612>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6177614>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6177617>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6177618>

### GSE189843

GSE189843 contains 12 pretreatment high-grade serous ovarian carcinoma samples:

- `GSM5708485`–`GSM5708490`: excellent response to neoadjuvant chemotherapy
- `GSM5708491`–`GSM5708496`: poor response to neoadjuvant chemotherapy

The GEO supplementary archive contains count matrices, barcodes, features and tissue images, but no `tissue_positions*` or `scalefactors_json.json` files. Therefore:

- the samples can be used for expression-level scoring and patient-level comparisons;
- they cannot be used for coordinate-aware neighborhood, distance or spatial ligand–receptor claims unless coordinates are obtained from the authors or reconstructed from an authoritative source.

Official record:

- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189843>

## Directory contents

```text
spatial_transcriptomics_upgrade/
├── README.md
├── config/
│   └── spatial_config.yml
├── metadata/
│   ├── download_manifest.csv
│   └── spatial_sample_manifest.csv
└── scripts/
    ├── 01_download_spatial_geo.py
    ├── 02_build_spatial_objects.R
    ├── 03_score_and_spatial_statistics.R
    ├── 04_reference_mapping_to_cnv_niches.R
    └── 05_audit_spatial_outputs.py
```

## Recommended run order

Run from Windows PowerShell or the Codex terminal. Data and intermediate objects are placed on the D drive by default.

### 1. Download and curate GEO files

```powershell
python scripts\01_download_spatial_geo.py `
  --root D:\OC_spatiogenomics\spatial_data
```

This step:

- downloads the two valid GSE203612 ovarian Visium samples directly from GEO;
- creates SpaceRanger-compatible folders;
- downloads and safely extracts the GSE189843 series archive;
- records file sizes and SHA-256 hashes in `logs/download_audit.json`;
- never downloads `GSM6177618`.

A dry run is available:

```powershell
python scripts\01_download_spatial_geo.py --dry-run
```

### 2. Build curated Seurat objects

Required R packages:

```r
install.packages(c("Seurat", "Matrix", "data.table", "yaml", "RANN"))
```

Then run:

```powershell
Rscript scripts\02_build_spatial_objects.R D:\OC_spatiogenomics\spatial_data
```

Outputs:

- `processed/spatial_objects_curated.rds`
- `processed/spatial_qc_raw_summary.csv`
- `processed/spatial_object_build_log.csv`

GSE203612 objects retain Visium coordinates. GSE189843 objects are explicitly marked `expression_only`.

### 3. Perform QC, scoring and spatial statistics

Review `config/spatial_config.yml`, then run:

```powershell
Rscript scripts\03_score_and_spatial_statistics.R `
  config\spatial_config.yml
```

Outputs:

- `results/spatial_curated/spatial_correlation_curated.csv`
- `results/spatial_curated/spatial_neighborhood_enrichment_curated.csv`
- `results/spatial_curated/spatial_spot_scores_curated.csv.gz`
- `results/spatial_curated/spatial_qc_filtered_summary.csv`

The script keeps expression-level correlation and coordinate-aware neighborhood analysis separate. Neighborhood analysis is generated only for samples with released coordinates.

The current neighborhood test is an exploratory k-nearest-neighbor label-permutation analysis. It does not fully model spatial autocorrelation and should not be presented as definitive physical interaction evidence.

### 4. Map single-cell and CNV reference states

```powershell
Rscript scripts\04_reference_mapping_to_cnv_niches.R `
  D:\OC_spatiogenomics\spatial_data `
  D:\OC_spatiogenomics\infercnv\integrated_oc.RData `
  D:\OC_spatiogenomics\infercnv\integrated_oc_plan_analysis\tables\integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv
```

This transfers broad myeloid states and CNV subclone labels to Visium spots. Prediction scores represent mixed-spot composition or state evidence; they do not identify one physical cell per spot.

### 5. Audit the outputs

```powershell
python scripts\05_audit_spatial_outputs.py `
  --results D:\OC_spatiogenomics\spatial_data\results\spatial_curated `
  --strict-results
```

The audit fails if:

- `GSM6177618` enters an ovarian result;
- a coordinate-unavailable sample enters a neighborhood result;
- curated samples are missing or duplicated.

## Interpretation rules

Use the following language consistently:

- **GSE203612:** coordinate-aware but limited to two ovarian samples with discordant directions; evidence is exploratory and heterogeneous.
- **GSE189843:** expression-level validation only until authoritative coordinates are available.
- **SPP1-CD44:** a candidate communication axis supported by expression and LR databases.
- **SPP1-ITGB1:** an SPP1-associated ITGB1-positive adhesion/integrin program; not proof of direct SPP1–ITGB1 binding.
- **Spatial expression product:** potential co-occurrence, not physical ligand–receptor interaction.
- **Reference mapping:** probabilistic spot composition/state inference, not direct CNV measurement in tissue.

## Remaining data gaps

For manuscript-level spatial validation, obtain or generate:

1. original GSE189843 spot coordinates and scalefactors from the authors;
2. pathology annotations for tumor, stroma, necrosis and immune-rich regions;
3. additional ovarian Visium/CosMx/Xenium cohorts with patient-level replication;
4. a formal deconvolution model such as cell2location or RCTD using `integrated_oc`;
5. patient-level meta-analysis rather than pooling spots as independent replicates;
6. orthogonal validation by multiplex immunofluorescence or RNAscope.
