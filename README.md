# DDA-PDX Prediction of targeted drug efficacy in cancer avatars with computational reasoning

Reproducibility repository for the manuscript:

> **Prediction of targeted drug efficacy in cancer avatars with computational reasoning**

The DOI will be added upon publication.

---

## Overview

This repository contains the analysis scripts and input data that were used to generate
the figures for the manuscript. All analyses were performed on patient-derived
xenograft (PDX) pharmacology data. The DDA (Digital Drug Assignment) scoring
framework assigned each compound–patient pair a numerical score, referred to as LEVEL.

> **Note:** The DDA algorithm is proprietary and cannot be shared. This
> repository contains only the precomputed LEVEL scores that were used as input to the
> statistical analyses; no part of the DDA scoring engine or its underlying
> methodology was included.

**Primary input file (all analyses):**  
`data/pdx_curve_metrics_single_treatments_dda_scores.xlsx`

Additional input files used by some scripts are listed in `data/DATA_MANIFEST.md`.

---

## Repository Structure

```
PDX_MS/
├── data/
│   ├── DATA_MANIFEST.md           # Column definitions and file usage
│   ├── pdx_curve_metrics_single_treatments_dda_scores.xlsx
│   ├── chemo_data_2024_12_06.csv
│   └── non_assigned_merged.xlsx
├── src/
│   ├── scores_dist_plots.py       # Generated DDA score distributions + waterfall plots
│   ├── KM_PDX.py                  # Performed Kaplan-Meier survival analysis
│   ├── HR_PDX.py                  # Generated hazard-ratio forest plots
│   ├── DCR_ORR_DB_PDX.py          # Performed DCR, ORR, and durable benefit (DB) analysis
│   ├── pie_chart_final.py         # Generated primary tumour site distribution pie chart
│   └── data_modelling.R             # Linear mixed models, SCAM/GAM regression, Cochran-Armitage trend test
├── requirements.txt
├── install_R_packages.R
├── software_references.md
└── README.md
```

---

## Requirements

Minimum install versions are listed below. Exact versions used in the manuscript are
recorded in [`software_references.md`](software_references.md).

### Python (≥ 3.10)

```bash
pip install -r requirements.txt
```

| Package | Version | Purpose |
|---------|---------|---------|
| pandas | ≥ 2.2 | Data manipulation |
| numpy | ≥ 1.26 | Numerical operations |
| lifelines | ≥ 0.29 | KM and Cox PH models |
| matplotlib | ≥ 3.9 | Figure rendering |
| seaborn | ≥ 0.13 | Statistical graphics |
| scipy | ≥ 1.13 | Statistical tests |
| statsmodels | ≥ 0.14 | Benjamini–Hochberg FDR correction |
| openpyxl | ≥ 3.1 | Excel I/O |

### R (≥ 4.4)

```r
Rscript install_R_packages.R
```

The following R packages were used: `ggplot2`, `lme4`, `lmerTest`, `scam`, `mgcv`, `survival`,
`survminer`, `pROC`, `patchwork`, `ggradar`, `readxl`, `openxlsx`.

---

## Running the Analyses

In the manuscript, all scripts were run from the repository root. To reproduce the figures, run:

```bash
# 1. DDA score distribution + waterfall plots
python src/scores_dist_plots.py

# 2. Kaplan-Meier survival analysis
python src/KM_PDX.py

# 3. Hazard-ratio forest plots
python src/HR_PDX.py

# 4. DCR / ORR / durable benefit analysis
python src/DCR_ORR_DB_PDX.py

# 5. Primary tumour site distribution pie chart
python src/pie_chart_final.py

# 6. Data modelling (linear mixed models, SCAM/GAM, Cochran-Armitage trend test)
Rscript src/data_modelling.R
```

---

## Figure Reproduction Table

| Figure | Script | Key output |
|--------|--------|------------|
| DDA score compound distribution | `scores_dist_plots.py` | `scores_dist_plots_{timestamp}/figures/PDX_compound_plot_{date}.jpg` |
| DDA score patient distribution | `scores_dist_plots.py` | `scores_dist_plots_{timestamp}/figures/PDX_score_dist_compound_plot_by_case_{date}.jpg` |
| Waterfall — BestResponse (Fig. S2A) | `scores_dist_plots.py` | `scores_dist_plots_{timestamp}/figures/PDX_waterfall_Figure_S2A_{date}_BestResponse.png` |
| Waterfall — BestAvgResponse (Fig. S2B) | `scores_dist_plots.py` | `scores_dist_plots_{timestamp}/figures/PDX_waterfall_Figure_S2B_{date}_BestAvgResponse.png` |
| KM — SHIVA tiers | `KM_PDX.py` | `KM_SHIVA/` |
| KM — on-label utility | `KM_PDX.py` | `KM_utility/` |
| KM — compound ranking | `KM_PDX.py` | `KM_ranking/` |
| HR forest — on-label | `HR_PDX.py` | `HR_utility/forest_plot.jpg` |
| HR forest — ranking | `HR_PDX.py` | `HR_Ranking/forest_plot.jpg` |
| HR forest — SHIVA | `HR_PDX.py` | `HR_SHIVA/forest_plot.jpg` |
| DCR / ORR / DB — on-label, ranking, SHIVA | `DCR_ORR_DB_PDX.py` | `src/DCR_ORR_DB/` |
| Primary tumour site pie chart | `pie_chart_final.py` | `Primary_Tumor_Site_distribution/distribution_Primary_Tumor_Site_{date}.jpeg` |
| Data modelling — SCAM residuals, effect forest, SCAM vs GAM, mRECIST, durable benefit | `data_modelling.R` | `src/data_modeling_outputs/` |

`{date}` denotes the run date stamp (`YYYY_MM_DD`) appended to output filenames.
`{timestamp}` is the run folder suffix created by `scores_dist_plots.py`.

**Score distribution and waterfall plots:** `scores_dist_plots.py` generates all four
figures (compound distribution, patient distribution, and both waterfall panels).

**Data modelling outputs:** `data_modelling.R` writes figures under `src/data_modeling_outputs/`,
including `SCAM_residual_diagnostics.png`, `Effect_forest_enhanced.png`, `SCAM_vs_GAM_all_ranks.png`,
mRECIST response panels, durable-benefit Cochran-Armitage panels, and leave-one-TARGET-out
SCAM vs GAM panels per treatment target (section 3). The script applies linear mixed models,
SCAM/GAM regression models, and Cochran-Armitage trend tests.

All scripts use `data/pdx_curve_metrics_single_treatments_dda_scores.xlsx` as an input.

---

## Reproducibility Notes

- Random seed was set to **42** in all Python scripts and R modules.
- P-values were reported with raw values and Benjamini–Hochberg false discovery rate
  (BH-FDR) corrected values. BH-FDR is the multiple-testing correction applied via
  `statsmodels.stats.multitest.multipletests` (`method="fdr_bh"`).

---

## Software references

See [`software_references.md`](software_references.md) for package versions, urls, and
rrid identifiers used in the manuscript analyses.

---

## License

This repository was made available for academic reproducibility purposes.
Please cite the manuscript when using this code or data.
