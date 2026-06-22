# Data Manifest — PDX_MS Repository

All input data files required to reproduce the manuscript analyses are
located in this folder (`data/`).

---

## Files

### 1. `pdx_curve_metrics_single_treatments_dda_scores.xlsx`
**Primary source file — used by all analysis scripts.**

| Column | Type | Description |
|--------|------|-------------|
| `Model` | string | Patient-derived xenograft (PDX) model identifier |
| `COMPOUND` | string | Treatment compound name |
| `LEVEL` | numeric | DDA (Digital Drug Assignment) score |
| `TimeToDouble` | numeric | Tumour doubling time (days); primary survival endpoint |
| `Day_Last` | numeric | Last observation day |
| `APPROVED` | string | Regulatory approval status (`IGAZ` = approved, `HAMIS` = not approved) |
| `On_label` | string | On-label use for this indication (`IGAZ` = on-label, `HAMIS` = off-label) |
| `TARGET` / `Treatment target` | string | Molecular target class of the compound |
| `BestResponse` | numeric | Best percentage change from baseline in tumour volume |
| `BestAvgResponse` | numeric | Average percentage change from baseline across measurement time points |
| `ResponseCategory` | string | Clinical response category |

**Used by:** all scripts in `src/`

---

### 2. `chemo_data_2024_12_06.csv`
**Standard-of-care chemotherapy reference cohort (CSV format).**

| Column | Type | Description |
|--------|------|-------------|
| `Model` | string | PDX model identifier |
| `COMPOUND` | string | Chemotherapy agent name |
| `Treatment target` | string | `chemotherapy` or `Tubulin` |
| `TimeToDouble` | numeric | Tumour doubling time (days) |

**Used by:** `HR_PDX.py`, `KM_PDX.py`, `DCR_ORR_DB_PDX.py`, `pie_chart_final.py`

---

### 3. `non_assigned_merged.xlsx`
**PDX records without molecular-profiling assignment (Non_assigned group).**

| Column | Type | Description |
|--------|------|-------------|
| `Model` | string | PDX model identifier |
| `COMPOUND` | string | Treatment compound name |
| `TimeToDouble` | numeric | Tumour doubling time (days) |

**Used by:** `HR_PDX.py`, `KM_PDX.py`, `DCR_ORR_DB_PDX.py`, `pie_chart_final.py`

---

## Notes

- PDX models are identified by alphanumeric codes only; no personally identifiable information is included.
