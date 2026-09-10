# DDA-PDX

Reproducibility repository for the manuscript:

> **Prediction of targeted drug efficacy in cancer avatars with computational reasoning**

The DOI will be added upon publication.

This package contains the analysis scripts and input tables used to generate the manuscript figures from patient-derived xenograft (PDX) pharmacology data. The Digital Drug Assignment (DDA) framework assigns each compound–patient pair a numerical score (`LEVEL`).

> **Note:** The DDA algorithm is proprietary and is not included. This repository provides the precomputed DDA scores that were used as input to the statistical analyses. No part of the scoring engine or its methodology is shared.

Generated figures, tables, logs, and assembled figure panels are **not** stored in this repository. Running the scripts from a clean clone recreates them locally.

---

## Overview of the analysis

The primary cohort is the PCT drug-response table in Supplementary File 1. Scripts score, stratify, and model treatment outcomes by DDA rank, then compare those ranks with standard-of-care chemotherapy and with a non-assigned PDX reference set.

**Primary input (all analyses):**  
`data/Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx` (sheet `PCT_DRUG_RESPONSE`)

Additional inputs are listed in [`data/DATA_MANIFEST.md`](data/DATA_MANIFEST.md).

---

## Repository structure

```
DDA-PDX/
├── data/                              # Input data and source-data builders
│   ├── DATA_MANIFEST.md
│   ├── INVENTORY.md                   # How to build Nature-style source-data workbooks
│   ├── Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx
│   ├── chemo_data_2024_12_06.csv
│   ├── non_assigned_merged.xlsx
│   ├── Fig1a_PDX_Graphical_abstract.jpg
│   ├── build_source_data.py
│   └── build_negative_scores_supplementary_table.py
├── src/                               # Analysis scripts (outputs are written here when run)
│   ├── scores_dist_plots.py           # Fig. 2d — DDA score by compound
│   ├── Case_level_distribution.R      # Case-level score ranking panel
│   ├── Waterfall_DDA_score_final.py   # Fig. 2e–f, Fig. 3d–e waterfalls
│   ├── KM_PDX.py                      # Fig. 3f, Fig. 4a, 4b, 4g Kaplan–Meier
│   ├── HR_PDX.py                      # Fig. 3g hazard-ratio forest
│   ├── DCR_ORR_DB.R                   # Fig. 3a–c DCR, ORR, durable benefit
│   ├── data_modelling.R               # Fig. 3h–i and Extended Fig. 2
│   ├── clinical_utility_ranking.py    # Fig. 4c–f
│   ├── pie_chart_final.py             # Fig. 2a tumour-site pie
│   ├── box_plot_final.py              # Fig. 2b alteration-count plots
│   ├── moa_pie_chart_final.py         # Fig. 2c mechanism-of-action pie
│   ├── flowchart_final.py             # Fig. 1b data-processing flowchart
│   ├── build_figure_panel_pack.py     # Assembles numbered manuscript panels
│   ├── plot_typography.py / .R        # Shared figure typography
│   └── cochran_armitage_trend.R       # Shared Cochran–Armitage helper
├── requirements.txt
├── install_R_packages.R
├── run_all_analyses.ps1
├── software_references.md
└── README.md
```

---

## Requirements

Minimum install versions are listed below. Exact versions used for the manuscript analyses are recorded in [`software_references.md`](software_references.md).

### Python (≥ 3.10)

```bash
pip install -r requirements.txt
```

| Package     | Version | Purpose                           |
| ----------- | ------- | --------------------------------- |
| pandas      | ≥ 2.2   | Data manipulation                 |
| numpy       | ≥ 1.26  | Numerical operations              |
| lifelines   | ≥ 0.29  | Kaplan–Meier and Cox PH models    |
| matplotlib  | ≥ 3.9   | Figure rendering                  |
| seaborn     | ≥ 0.13  | Statistical graphics              |
| scipy       | ≥ 1.13  | Statistical tests                 |
| statsmodels | ≥ 0.14  | Benjamini–Hochberg FDR correction |
| openpyxl    | ≥ 3.1   | Excel I/O                         |
| Pillow      | ≥ 10.0  | Raster figure assembly            |

### R (≥ 4.4)

```bash
Rscript install_R_packages.R
```

The following R packages are used: `ggplot2`, `dplyr`, `tidyr`, `readxl`, `writexl`, `openxlsx`, `lme4`, `lmerTest`, `scam`, `mgcv`, `survival`, `survminer`, `pROC`, `patchwork`, `ggrepel`, `gridExtra`, `scales`.

---

## Running the analyses

Clone the repository and run all scripts **from the repository root**. On Windows:

```powershell
.\run_all_analyses.ps1
```

Or run the steps individually:

```bash
python src/scores_dist_plots.py --project-root src
python src/Waterfall_DDA_score_final.py --project-root src
Rscript src/Case_level_distribution.R --project-root src
python src/KM_PDX.py --output-dir src
python src/clinical_utility_ranking.py --output-dir src
python src/HR_PDX.py --output-dir src
Rscript src/data_modelling.R
Rscript src/DCR_ORR_DB.R --output-dir src
python src/pie_chart_final.py
python src/box_plot_final.py
python src/moa_pie_chart_final.py
python src/flowchart_final.py
python src/build_figure_panel_pack.py --output-dir src
```

After the analysis scripts have written their figures, Nature-style source-data workbooks can be built with:

```bash
python data/build_source_data.py
```

Those workbooks are written to `data/` (`Source_Data_*.xlsx`) and are not part of the committed tree.

---

## Generated output (created locally; not in git)

| Folder | Contents |
|--------|----------|
| `src/score_distributions/` | Fig. 2d compound plot, case-level ranking panel, Fig. 2e–f and Fig. 3d–e waterfalls |
| `src/kaplan_meier/ranking/` | Fig. 3f and Fig. 4a, 4b, 4g Kaplan–Meier curves |
| `src/hazard_ratios/` | Fig. 3g Cox forest plot and pairwise HR tables |
| `src/dcr_orr_db/` | Fig. 3a–c DCR, ORR, and durable-benefit figures |
| `src/data_modeling/` | Fig. 3h–i SCAM/GAM/LMM figures and Extended Fig. 2 leave-one-out panels |
| `src/clinical_utility_ranking/` | Fig. 4c–f clinical-utility panels |
| `src/tumor_site/` | Fig. 2a primary tumour site pie chart |
| `src/variant_counts/` | Fig. 2b alteration-count plots |
| `src/moa_distribution/` | Fig. 2c mechanism-of-action pie chart |
| `src/flowchart/` | Fig. 1b data-processing flowchart |
| `src/figure_panels/` | Assembled manuscript figure panels |
| `data/Source_Data_*.xlsx` | Nature-style source-data workbooks from `build_source_data.py` |
| `src/logs/` | Per-step logs from `run_all_analyses.ps1` |

Fig. 1a is source artwork (`data/Fig1a_PDX_Graphical_abstract.jpg`), not a scripted plot.

---

## Figure reproduction table

| Figure | Script | Local output after a run |
|--------|--------|--------------------------|
| Fig. 1a | — (source artwork) | `data/Fig1a_PDX_Graphical_abstract.jpg` |
| Fig. 1b | `flowchart_final.py` | `src/flowchart/` |
| Fig. 2a | `pie_chart_final.py` | `src/tumor_site/` |
| Fig. 2b | `box_plot_final.py` | `src/variant_counts/` |
| Fig. 2c | `moa_pie_chart_final.py` | `src/moa_distribution/` |
| Fig. 2d | `scores_dist_plots.py` | `src/score_distributions/` |
| Case-level ranking panel | `Case_level_distribution.R` | `src/score_distributions/` |
| Fig. 2e, 2f, Fig. 3d, 3e | `Waterfall_DDA_score_final.py` | `src/score_distributions/` |
| Fig. 3a–c | `DCR_ORR_DB.R` | `src/dcr_orr_db/figures/` |
| Fig. 3f, Fig. 4a, 4b, 4g | `KM_PDX.py` | `src/kaplan_meier/ranking/` |
| Fig. 3g | `HR_PDX.py` | `src/hazard_ratios/` |
| Fig. 3h–i, Extended Fig. 2 | `data_modelling.R` | `src/data_modeling/` |
| Fig. 4c–f | `clinical_utility_ranking.py` | `src/clinical_utility_ranking/` |
| Assembled panels | `build_figure_panel_pack.py` | `src/figure_panels/` |
| Source-data workbooks | `data/build_source_data.py` | `data/Source_Data_*.xlsx` |

---

## Reproducibility notes

- Random seed is set to **42** in the Python scripts and R modules.
- *P*-values are reported as raw values and as Benjamini–Hochberg false-discovery-rate (BH-FDR) corrected values (`statsmodels.stats.multitest.multipletests`, `method="fdr_bh"`).
- Input checksums are recorded in `run_metadata.json` inside each analysis output folder after a run.
- This clone is intentionally free of analysis output so reviewers can confirm that the committed scripts regenerate the figures.

---

## License

This repository is made available for academic reproducibility and peer review. Please cite the manuscript when using this code or data.

---

## Contact

- **Róbert Dóczi** — robert.doczi@genomate.health
- **István Peták** — istvan.petak@genomate.health
