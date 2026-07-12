param(
    [string]$RepoRoot = "D:\OC_spatiogenomics\repo\OC_spatiogenomics",
    [string]$DataRoot = "D:\OC_spatiogenomics\infercnv\external_cell_annotations",
    [switch]$SkipDownload,
    [switch]$IncludeOptionalRaw,
    [string]$OnlyDataset = "",
    [switch]$Resume,
    [switch]$Force
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-Step([string]$Message) {
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    Write-Host "[$stamp] $Message"
}

function Invoke-Logged {
    param([string]$Id, [scriptblock]$Command)
    $log = Join-Path $DataRoot "logs\$Id.log"
    $marker = Join-Path $DataRoot "logs\$Id.done"
    if ($Resume -and (Test-Path $marker) -and -not $Force) {
        Write-Step "SKIP $Id"
        return
    }
    Write-Step "START $Id"
    & $Command *>&1 | Tee-Object -FilePath $log
    if ($LASTEXITCODE -ne 0) { throw "$Id failed with exit code $LASTEXITCODE" }
    Set-Content -Path $marker -Value (Get-Date).ToString("o") -Encoding UTF8
    Write-Step "DONE $Id"
}

$WorkflowRoot = Split-Path -Parent $PSScriptRoot
$env:TEMP = Join-Path $DataRoot "tmp"
$env:TMP = Join-Path $DataRoot "tmp"
$env:PIP_CACHE_DIR = Join-Path $DataRoot "cache\pip"
$env:R_USER = Join-Path $DataRoot "R_user"
$env:R_LIBS_USER = Join-Path $DataRoot "R_library"
New-Item -ItemType Directory -Force -Path $env:TEMP, $env:PIP_CACHE_DIR, $env:R_USER, $env:R_LIBS_USER, (Join-Path $DataRoot "logs") | Out-Null

$Rscript = "D:\R\R-4.0.3\bin\Rscript.exe"
if (-not (Test-Path $Rscript)) { $Rscript = "Rscript" }

Push-Location $WorkflowRoot
try {
    if (-not $SkipDownload) {
        $downloadArgs = @("scripts\00_download_original_annotations.py", "--data-root", $DataRoot)
        if ($IncludeOptionalRaw) { $downloadArgs += "--include-optional-raw" }
        if ($Force) { $downloadArgs += "--force" }
        Invoke-Logged "00_download_original_annotations" { python @downloadArgs }
    }
    Invoke-Logged "01_download_and_extract_gse154600_sce" { & $Rscript scripts\01_download_and_extract_gse154600_sce.R $DataRoot }
    Invoke-Logged "02_extract_original_annotations" { python scripts\02_extract_original_annotations.py }
    Invoke-Logged "03_harmonize_cell_annotations" { python scripts\03_harmonize_cell_annotations.py }
    Invoke-Logged "04_match_annotations_to_expression" { python scripts\04_match_annotations_to_expression.py }
    Invoke-Logged "05_secondary_annotation_gse147082" { & $Rscript scripts\05_secondary_annotation_gse147082.R $DataRoot }
    Invoke-Logged "06_prepare_normal_reference_gse151214" { & $Rscript scripts\06_prepare_normal_reference_gse151214.R $DataRoot }
    Invoke-Logged "07_rerun_external_scrna_curated" { python scripts\07_rerun_external_scrna_curated.py }
    Invoke-Logged "08_compare_original_vs_marker_annotations" { & $Rscript scripts\08_compare_original_vs_marker_annotations.R }
    Invoke-Logged "09_generate_annotation_upgrade_report" { python scripts\09_generate_annotation_upgrade_report.py }
} finally {
    Pop-Location
}
