# install_R_packages.R
#
# Installs all R packages required by data_modelling.R.
# Run once before executing the analysis:
#   Rscript install_R_packages.R
#
# Package versions tested with R 4.4.x on Windows and macOS.

required_packages <- c(
  "ggplot2",       # >= 3.5.0  -- publication figures
  "dplyr",         # >= 1.1.4  -- data manipulation
  "tidyr",         # >= 1.3.1  -- data reshaping
  "readxl",        # >= 1.4.3  -- Excel input
  "writexl",       # >= 1.5.0  -- Excel output (data_modelling.R)
  "openxlsx",      # >= 4.2.6  -- Excel output with formatting
  "stringr",       # >= 1.5.1  -- string helpers (data_modelling.R)
  "tibble",        # >= 3.2.1  -- tidy data frames (data_modelling.R)
  "gridExtra",     # >= 2.3   -- multi-panel figures (data_modelling.R)
  "lme4",          # >= 1.1-35 -- mixed-effects models
  "lmerTest",      # >= 3.1-3  -- Satterthwaite p-values for lmer
  "scam",          # >= 1.2-16 -- shape-constrained additive models (SCAM)
  "mgcv",          # >= 1.9-1  -- unconstrained GAM (GAM reference)
  "survival",      # >= 3.7-0  -- Kaplan-Meier and Cox PH
  "survminer",     # >= 0.4.9  -- KM plot utilities
  "pROC",          # >= 1.18.5 -- ROC / AUC performance metrics
  "ggradar",       # >= 0.2    -- radar / spider charts
  "ggrepel",       # >= 0.9.5  -- non-overlapping text labels
  "patchwork",     # >= 1.2.0  -- figure composition
  "scales",        # >= 1.3.0  -- axis scale helpers
  "RColorBrewer",  # >= 1.1-3  -- colour palettes
  "tools",         # base R   -- file_path_sans_ext (no install needed)
  "utils"          # base R   -- write.csv (no install needed)
)

# Packages available only from GitHub
github_packages <- list(
  # ggradar from GitHub (CRAN version may be outdated)
  # list(repo = "ricardo-bion/ggradar", pkg = "ggradar")
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) == 0L) {
    message("All CRAN packages already installed.")
    return(invisible(NULL))
  }
  message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org", quiet = FALSE)
}

install_if_missing(required_packages)

# Uncomment to install GitHub-only packages:
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# for (gh in github_packages) {
#   if (!requireNamespace(gh$pkg, quietly = TRUE)) remotes::install_github(gh$repo)
# }

message("R package installation complete.")
message("Verify with: sessionInfo()")
