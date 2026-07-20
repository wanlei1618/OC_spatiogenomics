[CmdletBinding()]
param(
    [string]$RepoRoot = 'D:\OC_spatiogenomics\sync_compare\clean_clone',
    [string]$DataRoot = 'D:\OC_spatiogenomics\infercnv\external_seurat_preannotation',
    [string]$Datasets = 'GSE154600,GSE158722',
    [string]$Rscript = 'D:\R\R-4.6.1\bin\Rscript.exe',
    [int]$MaxGlobalCells = 30000,
    [double]$Resolution = 0.6,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$Script = Join-Path $RepoRoot `
    'infercnv\chatgpt_review_bundle\06_external_seurat_preannotation\workflow\scripts\09_build_marker_ready_annotation.R'
$Config = Join-Path $RepoRoot `
    'infercnv\chatgpt_review_bundle\06_external_seurat_preannotation\workflow\config\diagnostics_v2.yaml'
$LogRoot = Join-Path $DataRoot 'diagnostics_v2_marker_ready\logs'

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
$stdout = Join-Path $LogRoot '09_marker_ready.stdout.log'
$stderr = Join-Path $LogRoot '09_marker_ready.stderr.log'

$args = @(
    $Script,
    '--config', $Config,
    '--data-root', $DataRoot,
    '--repo-root', $RepoRoot,
    '--datasets', $Datasets,
    '--max-global-cells', $MaxGlobalCells,
    '--resolution', $Resolution
)
if ($Force) {
    $args += '--force'
}

Write-Host "Running marker-ready annotation workflow for $Datasets"
$process = Start-Process `
    -FilePath $Rscript `
    -ArgumentList $args `
    -WorkingDirectory $RepoRoot `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr `
    -Wait `
    -PassThru

if ($process.ExitCode -ne 0) {
    throw "Marker-ready workflow failed with exit code $($process.ExitCode). See $stderr"
}

$datasetList = @(
    $Datasets -split ',' |
    ForEach-Object { $_.Trim() } |
    Where-Object { $_ }
)

foreach ($dataset in $datasetList) {
    $template = Join-Path $DataRoot `
        "diagnostics_v2_marker_ready\$dataset\annotation_ready_cluster_template.csv"
    if (-not (Test-Path -LiteralPath $template)) {
        throw "Expected annotation template missing: $template"
    }
    Write-Host "Ready: $template"
}

Write-Host "Marker-ready annotation workflow completed."
