[CmdletBinding()]
param(
    [string]$RepoRoot = 'D:\OC_spatiogenomics\sync_compare\clean_clone',
    [string]$DataRoot = 'D:\OC_spatiogenomics\infercnv\external_seurat_preannotation',
    [string]$Datasets = 'GSE154600,GSE158722',
    [string]$Rscript = 'D:\R\R-4.6.1\bin\Rscript.exe',
    [int]$KnnK = 30,
    [double]$PredictionScore = 0.70,
    [double]$PredictionMargin = 0.20,
    [double]$ClusterFraction = 0.80,
    [int]$ClusterSupport = 30,
    [double]$GlobalResolution = 0.8,
    [double]$Resolution = 0.6,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$Script = Join-Path $RepoRoot `
    'infercnv\chatgpt_review_bundle\06_external_seurat_preannotation\workflow\scripts\10_refine_annotation_ready_clusters.R'
$Config = Join-Path $RepoRoot `
    'infercnv\chatgpt_review_bundle\06_external_seurat_preannotation\workflow\config\diagnostics_v2.yaml'
$LogRoot = Join-Path $DataRoot 'diagnostics_v2_marker_ready_refined\logs'

if (-not (Test-Path -LiteralPath $Rscript)) {
    throw "Rscript not found: $Rscript"
}
if (-not (Test-Path -LiteralPath $Script)) {
    throw "R script not found: $Script"
}
if (-not (Test-Path -LiteralPath $Config)) {
    throw "Config not found: $Config"
}
if (-not (Test-Path -LiteralPath $DataRoot)) {
    throw "DataRoot not found: $DataRoot"
}

New-Item -ItemType Directory -Force -Path $LogRoot | Out-Null
$stdout = Join-Path $LogRoot '10_refinement.stdout.log'
$stderr = Join-Path $LogRoot '10_refinement.stderr.log'

$arguments = @(
    $Script,
    '--config', $Config,
    '--data-root', $DataRoot,
    '--repo-root', $RepoRoot,
    '--datasets', $Datasets,
    '--knn-k', $KnnK,
    '--prediction-score', $PredictionScore,
    '--prediction-margin', $PredictionMargin,
    '--cluster-fraction', $ClusterFraction,
    '--cluster-support', $ClusterSupport,
    '--global-resolution', $GlobalResolution,
    '--resolution', $Resolution
)
if ($Force) {
    $arguments += '--force'
}

Write-Host "Running refined marker-ready workflow for $Datasets"
$process = Start-Process `
    -FilePath $Rscript `
    -ArgumentList $arguments `
    -WorkingDirectory $RepoRoot `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr `
    -Wait `
    -PassThru

if ($process.ExitCode -ne 0) {
    throw "Refined marker-ready workflow failed with exit code $($process.ExitCode). See $stderr"
}

$datasetList = @(
    $Datasets -split ',' |
    ForEach-Object { $_.Trim() } |
    Where-Object { $_ }
)

foreach ($dataset in $datasetList) {
    $template = Join-Path $DataRoot `
        "diagnostics_v2_marker_ready_refined\$dataset\annotation_ready_cluster_template_refined.csv"
    $assignments = Join-Path $DataRoot `
        "diagnostics_v2_marker_ready_refined\$dataset\annotation_ready_full_cell_assignments.csv.gz"
    if (-not (Test-Path -LiteralPath $template)) {
        throw "Expected refined template missing: $template"
    }
    if (-not (Test-Path -LiteralPath $assignments)) {
        throw "Expected full assignment missing: $assignments"
    }
    Write-Host "Ready: $template"
}

Write-Host 'Refined marker-ready workflow completed.'
