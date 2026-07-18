param(
    [string]$Config = "config\five_external_datasets.yaml",
    [string]$Datasets = "",
    [switch]$AuditOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$SkillRoot = Split-Path -Parent $PSScriptRoot
Push-Location $SkillRoot
try {
    $Candidates = @(
        "D:\R\R-4.6.1\bin\Rscript.exe",
        "D:\R\R-4.3.1\bin\Rscript.exe",
        "D:\R\R-4.3.0\bin\Rscript.exe",
        "D:\R\R-4.0.3\bin\Rscript.exe"
    )

    $Rscript = $null
    foreach ($candidate in $Candidates) {
        if (Test-Path $candidate) {
            $Rscript = $candidate
            break
        }
    }

    if (-not $Rscript) {
        $cmd = Get-Command Rscript -ErrorAction SilentlyContinue
        if ($cmd) { $Rscript = $cmd.Source }
    }

    if (-not $Rscript) {
        throw "Rscript was not found. Edit scripts\run_skill.ps1 with the correct R installation path."
    }

    $RunArgs = @(
        "scripts\run_seurat_preannotation.R",
        "--config", $Config
    )
    if ($AuditOnly) { $RunArgs += "--audit-only" }
    if ($Datasets) { $RunArgs += @("--datasets", $Datasets) }

    & $Rscript @RunArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Workflow stopped with exit code $LASTEXITCODE. Check each dataset's logs\run_status.json."
    }
}
finally {
    Pop-Location
}
