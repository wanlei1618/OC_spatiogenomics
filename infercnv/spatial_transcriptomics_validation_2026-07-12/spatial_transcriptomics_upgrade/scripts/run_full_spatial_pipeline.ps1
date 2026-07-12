param(
    [string]$RepoRoot = "D:\OC_spatiogenomics\repo\OC_spatiogenomics",
    [string]$DataRoot = "D:\OC_spatiogenomics\spatial_data",
    [string]$InfercnvRoot = "D:\OC_spatiogenomics\infercnv",
    [switch]$SkipDownload,
    [switch]$SkipInstall,
    [switch]$SkipReferenceMapping,
    [string]$OnlyStep = "",
    [switch]$Resume,
    [switch]$Force
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-Step([string]$Message) {
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    Write-Host "[$stamp] $Message"
}

function Assert-Command([string]$Name) {
    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        throw "Required command not found: $Name"
    }
}

function Invoke-Step {
    param(
        [string]$Id,
        [string]$Name,
        [scriptblock]$Command
    )
    if ($OnlyStep -and $OnlyStep -ne $Id) { return }
    $marker = Join-Path $StateRoot "$Id.done"
    $log = Join-Path $LogRoot "$Id.log"
    if ($Resume -and (Test-Path $marker) -and -not $Force) {
        Write-Step "SKIP $Name"
        $script:Steps += [pscustomobject]@{ step = $Id; status = "skipped_resume"; log = $log }
        return
    }
    Write-Step "START $Name"
    try {
        & $Command *>&1 | Tee-Object -FilePath $log
        if ($LASTEXITCODE -ne 0) { throw "$Name failed with exit code $LASTEXITCODE" }
        Set-Content -Path $marker -Value (Get-Date).ToString("o") -Encoding UTF8
        $script:Steps += [pscustomobject]@{ step = $Id; status = "success"; log = $log }
        Write-Step "DONE $Name"
    } catch {
        $script:Steps += [pscustomobject]@{ step = $Id; status = "failed"; log = $log; message = $_.Exception.Message }
        throw
    }
}

$RepoRoot = [System.IO.Path]::GetFullPath($RepoRoot)
$DataRoot = [System.IO.Path]::GetFullPath($DataRoot)
$InfercnvRoot = [System.IO.Path]::GetFullPath($InfercnvRoot)
$UpgradeRoot = Split-Path -Parent $PSScriptRoot

$dirs = @(
    $DataRoot, "$DataRoot\raw", "$DataRoot\extracted", "$DataRoot\processed",
    "$DataRoot\reference_mapping", "$DataRoot\results", "$DataRoot\figures",
    "$DataRoot\reports", "$DataRoot\logs", "$DataRoot\tmp", "$DataRoot\cache\pip",
    "$DataRoot\R_user", "$DataRoot\R_library"
)
foreach ($dir in $dirs) { New-Item -ItemType Directory -Force -Path $dir | Out-Null }

$env:TEMP = "$DataRoot\tmp"
$env:TMP = "$DataRoot\tmp"
$env:PIP_CACHE_DIR = "$DataRoot\cache\pip"
$env:R_USER = "$DataRoot\R_user"
$env:R_LIBS_USER = "$DataRoot\R_library"

$RunId = Get-Date -Format "yyyyMMdd_HHmmss"
$LogRoot = Join-Path $DataRoot "logs\spatial_pipeline_$RunId"
$StateRoot = Join-Path $DataRoot "logs\spatial_pipeline_state"
New-Item -ItemType Directory -Force -Path $LogRoot, $StateRoot | Out-Null
$Steps = @()
$Status = "success"

try {
    Assert-Command git
    Assert-Command python
    Assert-Command Rscript
    $freeGB = [math]::Round((Get-PSDrive D).Free / 1GB, 2)
    if ($freeGB -lt 20) { throw "D drive free space is below 20 GB: $freeGB GB" }

    Push-Location $UpgradeRoot
    try {
        if (-not $SkipInstall) {
            $installer = "D:\Downloads\OC_spatial_codex_package\install_spatial_dependencies.R"
            if (Test-Path $installer) {
                Invoke-Step "00_install" "Install R dependencies" { Rscript $installer $DataRoot }
            }
        }
        Invoke-Step "01_dry_run" "Download dry run" { python scripts\01_download_spatial_geo.py --dry-run --root $DataRoot }
        if (-not $SkipDownload) {
            Invoke-Step "01_download" "Download GEO data" { python scripts\01_download_spatial_geo.py --root $DataRoot }
        }
        Invoke-Step "02_build" "Build spatial objects" { Rscript scripts\02_build_spatial_objects.R $DataRoot }
        Invoke-Step "03_score" "Score and spatial statistics" { Rscript scripts\03_score_and_spatial_statistics.R config\spatial_config.yml }
        if (-not $SkipReferenceMapping) {
            $integrated = Join-Path $InfercnvRoot "00_raw_objects_and_infercnv\integrated_oc.RData"
            if (-not (Test-Path $integrated)) { $integrated = Join-Path $InfercnvRoot "integrated_oc.RData" }
            $metadata = Join-Path $InfercnvRoot "integrated_oc_plan_analysis\tables\integrated_oc_metadata_with_cnv_TNK_myeloid_B_subtypes.csv"
            Invoke-Step "04_reference" "Reference mapping" { Rscript scripts\04_reference_mapping_to_cnv_niches.R $DataRoot $integrated $metadata }
        }
        Invoke-Step "05_audit" "Strict audit" { python scripts\05_audit_spatial_outputs.py --results "$DataRoot\results\spatial_curated" --strict-results }
        Invoke-Step "06_qc_sensitivity" "QC sensitivity" { Rscript scripts\06_spatial_qc_sensitivity.R config\spatial_config.yml }
        Invoke-Step "07_autocorrelation" "Autocorrelation and multiscale neighborhoods" { Rscript scripts\07_spatial_autocorrelation_multiscale.R config\spatial_config.yml }
        Invoke-Step "08_directional" "Directional niche statistics" { Rscript scripts\08_directional_niche_statistics.R config\spatial_config.yml }
        Invoke-Step "09_mapping_stability" "Reference mapping stability" { Rscript scripts\09_reference_mapping_stability.R config\spatial_config.yml }
        Invoke-Step "10_response" "GSE189843 response analysis" { Rscript scripts\10_gse189843_response_analysis.R config\spatial_config.yml }
        Invoke-Step "11_meta" "Patient-level meta analysis" { Rscript scripts\11_patient_level_meta_analysis.R config\spatial_config.yml }
        Invoke-Step "12_figures" "Spatial figures" { Rscript scripts\12_spatial_figures.R config\spatial_config.yml }
        Invoke-Step "13_report" "Final report" { Rscript scripts\13_generate_spatial_report.R config\spatial_config.yml }
        Invoke-Step "99_final_audit" "Final strict audit" { python scripts\05_audit_spatial_outputs.py --results "$DataRoot\results\spatial_curated" --strict-results }
    } finally {
        Pop-Location
    }
} catch {
    $Status = "failed"
    $Steps += [pscustomobject]@{ step = "pipeline"; status = "failed"; message = $_.Exception.Message }
    throw
} finally {
    $summary = [pscustomobject]@{
        run_id = $RunId
        status = $Status
        repo_root = $RepoRoot
        data_root = $DataRoot
        infercnv_root = $InfercnvRoot
        upgrade_root = $UpgradeRoot
        steps = $Steps
        completed_at = (Get-Date).ToString("o")
    }
    $summaryPath = Join-Path $DataRoot "reports\pipeline_run_summary.json"
    $summary | ConvertTo-Json -Depth 8 | Set-Content -Path $summaryPath -Encoding UTF8
    Write-Step "Run summary: $summaryPath"
}
