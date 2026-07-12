$ErrorActionPreference = "Stop"

$destRoot = "D:\Codex_migrated_from_C"
New-Item -ItemType Directory -Force -Path $destRoot | Out-Null

$targets = @(
    @{ Source = "$env:USERPROFILE\.codex"; Dest = "$destRoot\user_profile\.codex" },
    @{ Source = "$env:APPDATA\Codex"; Dest = "$destRoot\AppData_Roaming\Codex" },
    @{ Source = "$env:LOCALAPPDATA\OpenAI"; Dest = "$destRoot\AppData_Local\OpenAI" },
    @{ Source = "$env:USERPROFILE\Documents\Codex"; Dest = "$destRoot\Documents\Codex" }
)

function Get-FolderStats {
    param([Parameter(Mandatory=$true)][string]$Path)

    $files = Get-ChildItem -LiteralPath $Path -Recurse -Force -File -ErrorAction SilentlyContinue
    $sum = ($files | Measure-Object Length -Sum).Sum
    [PSCustomObject]@{
        Count = @($files).Count
        Bytes = [int64]$sum
    }
}

foreach ($target in $targets) {
    $src = $target.Source
    $dst = $target.Dest

    if (-not (Test-Path -LiteralPath $src)) {
        Write-Host "Skip missing: $src"
        continue
    }

    $srcItem = Get-Item -LiteralPath $src -Force
    if (($srcItem.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {
        Write-Host "Already linked: $src"
        continue
    }

    Write-Host "Copying: $src -> $dst"
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $dst) | Out-Null
    robocopy $src $dst /MIR /R:1 /W:1 /XJ /NFL /NDL /NP | Out-Null
    if ($LASTEXITCODE -ge 8) {
        throw "Copy failed for $src with robocopy exit code $LASTEXITCODE"
    }

    $srcStats = Get-FolderStats -Path $src
    $dstStats = Get-FolderStats -Path $dst
    if ($srcStats.Count -ne $dstStats.Count -or $srcStats.Bytes -ne $dstStats.Bytes) {
        throw "Verification failed for $src"
    }

    $backup = "$src.C_migration_backup_$(Get-Date -Format yyyyMMddHHmmss)"
    Write-Host "Linking original path to D: $src"
    Rename-Item -LiteralPath $src -NewName (Split-Path -Leaf $backup)
    New-Item -ItemType Junction -Path $src -Target $dst | Out-Null
    Remove-Item -LiteralPath $backup -Recurse -Force
}

Write-Host ""
Write-Host "Migration complete. Current drive space:"
Get-PSDrive -PSProvider FileSystem |
    Select-Object Name,
        @{n="UsedGB";e={[math]::Round($_.Used / 1GB, 2)}},
        @{n="FreeGB";e={[math]::Round($_.Free / 1GB, 2)}},
        @{n="TotalGB";e={[math]::Round(($_.Used + $_.Free) / 1GB, 2)}} |
    Format-Table -AutoSize
