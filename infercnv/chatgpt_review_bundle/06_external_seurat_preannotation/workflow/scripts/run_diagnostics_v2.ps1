[CmdletBinding()]
param(
    [string]$RepoRoot = 'D:\OC_spatiogenomics\sync_compare\clean_clone',
    [string]$DataRoot = 'D:\OC_spatiogenomics\infercnv\external_seurat_preannotation',
    [string]$OnlyDataset = '',
    [switch]$AuditOnly,
    [switch]$Resume,
    [switch]$Force,
    [switch]$SkipRPCA,
    [switch]$SkipHarmony
)

$ErrorActionPreference = 'Stop'
$Rscript = 'D:\R\R-4.6.1\bin\Rscript.exe'
$WorkflowRoot = Join-Path $RepoRoot 'infercnv\chatgpt_review_bundle\06_external_seurat_preannotation\workflow'
$ScriptRoot = Join-Path $WorkflowRoot 'scripts'
$Config = Join-Path $WorkflowRoot 'config\diagnostics_v2.yaml'
$DiagnosticsRoot = Join-Path $DataRoot 'diagnostics_v2'
$LogRoot = Join-Path $DiagnosticsRoot 'logs'

if (-not (Test-Path -LiteralPath $Rscript)) { throw "Rscript not found: $Rscript" }
if (-not (Test-Path -LiteralPath $Config)) { throw "Config not found: $Config" }
if (-not (Test-Path -LiteralPath $RepoRoot)) { throw "RepoRoot not found: $RepoRoot" }
if (-not (Test-Path -LiteralPath $DataRoot)) { throw "DataRoot not found: $DataRoot" }

$Allowed = @('GSE147082', 'GSE154600', 'GSE158722')
$Datasets = if ([string]::IsNullOrWhiteSpace($OnlyDataset)) {
    $Allowed
} else {
    @($OnlyDataset -split ',' | ForEach-Object { $_.Trim() } | Where-Object { $_ })
}
$Invalid = @($Datasets | Where-Object { $_ -notin $Allowed })
if ($Invalid.Count -gt 0) { throw "Forbidden or unknown dataset(s): $($Invalid -join ', '). GSE154763 is never permitted." }

if ($Force -and (Test-Path -LiteralPath $DiagnosticsRoot)) {
    $resolvedData = [IO.Path]::GetFullPath($DataRoot).TrimEnd('\')
    $resolvedDiagnostics = [IO.Path]::GetFullPath($DiagnosticsRoot)
    if (-not $resolvedDiagnostics.StartsWith($resolvedData + '\', [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing archive move outside DataRoot: $resolvedDiagnostics"
    }
    $archive = Join-Path $DataRoot ("diagnostics_v2_archive_" + (Get-Date -Format 'yyyyMMdd_HHmmss'))
    Move-Item -LiteralPath $resolvedDiagnostics -Destination $archive
    Write-Host "Archived previous diagnostics_v2 to $archive"
}
New-Item -ItemType Directory -Force -Path $LogRoot | Out-Null

function Get-StepFingerprint {
    param([string]$ScriptPath, [string[]]$StepArguments)
    $scriptHash = (Get-FileHash -LiteralPath $ScriptPath -Algorithm SHA256).Hash
    $configHash = (Get-FileHash -LiteralPath $Config -Algorithm SHA256).Hash
    $text = "$scriptHash|$configHash|$RepoRoot|$DataRoot|$($Datasets -join ',')|$($StepArguments -join '|')|$SkipRPCA|$SkipHarmony"
    $bytes = [Text.Encoding]::UTF8.GetBytes($text)
    $sha = [Security.Cryptography.SHA256]::Create()
    try { ([BitConverter]::ToString($sha.ComputeHash($bytes))).Replace('-', '') } finally { $sha.Dispose() }
}

function Run-Step {
    param(
        [string]$Name,
        [string]$ScriptName,
        [string[]]$Arguments,
        [string]$MarkerPath
    )
    $scriptPath = Join-Path $ScriptRoot $ScriptName
    if (-not (Test-Path -LiteralPath $scriptPath)) { throw "Missing step script: $scriptPath" }
    $statePath = Join-Path $LogRoot "$Name.state.json"
    $fingerprint = Get-StepFingerprint -ScriptPath $scriptPath -StepArguments $Arguments
    if ($Resume -and (Test-Path -LiteralPath $statePath) -and (Test-Path -LiteralPath $MarkerPath)) {
        $state = Get-Content -LiteralPath $statePath -Raw | ConvertFrom-Json
        $markerHash = (Get-FileHash -LiteralPath $MarkerPath -Algorithm SHA256).Hash
        if ($state.status -eq 'COMPLETE' -and $state.fingerprint -eq $fingerprint -and $state.marker_sha256 -eq $markerHash) {
            Write-Host "[$Name] resume check passed; verified fingerprint and marker SHA-256."
            return
        }
        Write-Host "[$Name] stale/incomplete state detected; rerunning."
    }
    $log = Join-Path $LogRoot "$Name.log"
    $stdoutLog = Join-Path $LogRoot "$Name.stdout.log"
    $stderrLog = Join-Path $LogRoot "$Name.stderr.log"
    $started = Get-Date
    Write-Host "[$Name] starting"
    Remove-Item -LiteralPath $stdoutLog, $stderrLog, $log -ErrorAction SilentlyContinue
    $processArguments = @($scriptPath, '--config', $Config, '--data-root', $DataRoot,
                          '--repo-root', $RepoRoot) + $Arguments
    $process = Start-Process -FilePath $Rscript -ArgumentList $processArguments `
        -WorkingDirectory $RepoRoot -WindowStyle Hidden `
        -RedirectStandardOutput $stdoutLog -RedirectStandardError $stderrLog `
        -Wait -PassThru
    $exitCode = $process.ExitCode
    Get-Content -LiteralPath $stdoutLog, $stderrLog -ErrorAction SilentlyContinue |
        Set-Content -LiteralPath $log -Encoding UTF8
    if ($exitCode -ne 0) { throw "[$Name] failed with exit code $exitCode; see $log" }
    if (-not (Test-Path -LiteralPath $MarkerPath)) { throw "[$Name] marker output missing: $MarkerPath" }
    $state = [ordered]@{
        step = $Name
        status = 'COMPLETE'
        started_at = $started.ToString('o')
        finished_at = (Get-Date).ToString('o')
        fingerprint = $fingerprint
        marker_path = $MarkerPath
        marker_sha256 = (Get-FileHash -LiteralPath $MarkerPath -Algorithm SHA256).Hash
        script_sha256 = (Get-FileHash -LiteralPath $scriptPath -Algorithm SHA256).Hash
        config_sha256 = (Get-FileHash -LiteralPath $Config -Algorithm SHA256).Hash
    }
    $state | ConvertTo-Json | Set-Content -LiteralPath $statePath -Encoding UTF8
    Write-Host "[$Name] complete"
}

$datasetArg = $Datasets -join ','
Run-Step -Name '01_forensic_audit' -ScriptName '01_forensic_audit_current_results.R' `
    -Arguments @('--datasets', $datasetArg) `
    -MarkerPath (Join-Path $DiagnosticsRoot '00_forensic\current_result_snapshot.md')

if ($AuditOnly) {
    Write-Host 'AuditOnly requested; stopped after the frozen forensic audit.'
    exit 0
}

$MtDatasets = @($Datasets | Where-Object { $_ -in @('GSE147082', 'GSE158722') })
if ($MtDatasets.Count -gt 0) {
    $args02 = @('--datasets', ($MtDatasets -join ','))
    if ($Resume) { $args02 += '--resume' }
    $mtMarkerDataset = $MtDatasets[-1]
    Run-Step -Name '02_mitochondrial_qc' -ScriptName '02_diagnose_and_fix_mito_features.R' `
        -Arguments $args02 `
        -MarkerPath (Join-Path $DiagnosticsRoot "$mtMarkerDataset\01_mt_audit\mt_qc_decision.md")
}

$DominanceDatasets = @($Datasets | Where-Object { $_ -in @('GSE154600', 'GSE158722') })
if ($DominanceDatasets.Count -gt 0) {
    $domArg = $DominanceDatasets -join ','
    $markerDataset = $DominanceDatasets[-1]
    Run-Step -Name '03_sample_dominance' -ScriptName '03_diagnose_sample_dominance.R' `
        -Arguments @('--datasets', $domArg) `
        -MarkerPath (Join-Path $DiagnosticsRoot "$markerDataset\02_dominance\cluster_dominance_diagnostic_table.csv")
    foreach ($dataset in $DominanceDatasets) {
        Run-Step -Name "04_broad_lineages_$dataset" -ScriptName '04_build_provisional_broad_lineages.R' `
            -Arguments @('--datasets', $dataset) `
            -MarkerPath (Join-Path $DiagnosticsRoot "$dataset\03_broad_lineages\lineage_strategy_input_audit.csv")
        $args05 = @('--datasets', $dataset)
        if ($Resume) { $args05 += '--resume' }
        if ($SkipRPCA) { $args05 += '--skip-rpca' }
        if ($SkipHarmony) { $args05 += '--skip-harmony' }
        Run-Step -Name "05_batch_strategies_$dataset" -ScriptName '05_compare_batch_strategies.R' `
            -Arguments $args05 `
            -MarkerPath (Join-Path $DiagnosticsRoot "strategy_comparison\strategy_run_status_$dataset.csv")
    }
    Run-Step -Name '06_strategy_evaluation' -ScriptName '06_evaluate_batch_strategies.R' `
        -Arguments @('--datasets', $domArg) `
        -MarkerPath (Join-Path $DiagnosticsRoot 'strategy_comparison\recommended_strategy_by_dataset_and_lineage.csv')
    foreach ($dataset in $DominanceDatasets) {
        Run-Step -Name "07_rna_markers_$dataset" -ScriptName '07_rerun_markers_and_export.R' `
            -Arguments @('--datasets', $dataset) `
            -MarkerPath (Join-Path $DiagnosticsRoot "$dataset\05_markers\manual_annotation_template.csv")
    }
}

Run-Step -Name '08_review_package' -ScriptName '08_generate_review_package.R' `
    -Arguments @() `
    -MarkerPath (Join-Path $DiagnosticsRoot 'run_summary.json')

Write-Host "diagnostics_v2 workflow complete: $DiagnosticsRoot"
