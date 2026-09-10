# Re-run all PDX analysis scripts from the repository root.

$ErrorActionPreference = "Stop"
$RepoRoot = $PSScriptRoot
$Src = Join-Path $RepoRoot "src"
Set-Location $RepoRoot

$RunStamp = Get-Date -Format "yyyy_MM_dd__HH_mm_ss"
$LogDir = Join-Path $Src "logs"
New-Item -ItemType Directory -Path $LogDir -Force | Out-Null

function Invoke-Step {
    param(
        [string]$Name,
        [scriptblock]$Command
    )
    Write-Host "==> $Name"
    $log = Join-Path $LogDir "$Name.log"
    $prevErrorAction = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try {
        $output = & $Command 2>&1
        $exit = $LASTEXITCODE
        $output | Tee-Object -FilePath $log
    } finally {
        $ErrorActionPreference = $prevErrorAction
    }
    if ($exit -ne 0) {
        throw "$Name failed with exit code $exit (see $log)"
    }
}

$sharedArgs = @("--project-root", $Src, "--run-stamp", $RunStamp)

Invoke-Step "scores_dist_plots.py" {
    python (Join-Path $Src "scores_dist_plots.py") @sharedArgs
}
Invoke-Step "Waterfall_DDA_score_final.py" {
    python (Join-Path $Src "Waterfall_DDA_score_final.py") @sharedArgs
}
Invoke-Step "Case_level_distribution.R" {
    Rscript (Join-Path $Src "Case_level_distribution.R") @sharedArgs
}
Invoke-Step "KM_PDX.py" {
    python (Join-Path $Src "KM_PDX.py") --output-dir $Src
}
Invoke-Step "clinical_utility_ranking.py" {
    python (Join-Path $Src "clinical_utility_ranking.py") --output-dir $Src
}
Invoke-Step "HR_PDX.py" {
    python (Join-Path $Src "HR_PDX.py") --output-dir $Src
}
Invoke-Step "data_modelling.R" {
    Rscript (Join-Path $Src "data_modelling.R")
}
Invoke-Step "DCR_ORR_DB.R" {
    Rscript (Join-Path $Src "DCR_ORR_DB.R") --output-dir $Src
}
Invoke-Step "pie_chart_final.py" {
    python (Join-Path $Src "pie_chart_final.py")
}
Invoke-Step "box_plot_final.py" {
    python (Join-Path $Src "box_plot_final.py")
}
Invoke-Step "moa_pie_chart_final.py" {
    python (Join-Path $Src "moa_pie_chart_final.py")
}
Invoke-Step "flowchart_final.py" {
    python (Join-Path $Src "flowchart_final.py")
}
Invoke-Step "build_figure_panel_pack.py" {
    python (Join-Path $Src "build_figure_panel_pack.py") --output-dir $Src
}

Write-Host "All steps finished. Logs: $LogDir"
