# Breast Cancer Survival Analysis — TCGA BRCA

**Multivariate survival analysis of 1,035 breast cancer patients using clinical and receptor biomarker data from TCGA.**

---

## Overview

This project applies Kaplan-Meier survival analysis and Cox proportional hazards modelling to the TCGA Breast Cancer (BRCA) cohort. The goal is to identify clinical and molecular features that independently predict patient survival outcomes.

All data is sourced from the public TCGA repository via the `TCGAbiolinks` R package — no proprietary data is used.

---

## Key Findings

| Variable | Hazard Ratio | 95% CI | p-value |
|---|---|---|---|
| Age (per SD increase) | 1.84 | 1.38 – 2.44 | **< 0.001** |
| ER positive vs negative | 0.46 | 0.25 – 0.84 | **0.011** |
| HER2 positive vs negative | 1.82 | 0.97 – 3.43 | 0.064 |
| Stage IV vs Stage I | 3.46 | 0.95 – 12.54 | 0.059 |

**Model concordance: 0.856** (the model correctly ranks survival outcomes in 85.6% of patient pairs)

### Interpretation

- **Age** is the strongest independent predictor of mortality. Each standard deviation increase in age at diagnosis raises the hazard of death by 84% (HR=1.84, p<0.001), after controlling for stage and receptor status.
- **ER-positive status** is strongly protective. ER-positive patients have 54% lower hazard of death compared to ER-negative patients (HR=0.46, p=0.011). This reflects the availability of hormone-targeted therapies (e.g. tamoxifen) for ER-positive tumours.
- **HER2-positive status** shows a trend toward worse survival (HR=1.82, p=0.064) but does not reach conventional significance in this cohort, likely because TCGA data predates widespread Herceptin (trastuzumab) use — the survival disadvantage of HER2+ disease is more muted in modern treated cohorts.
- **Stage IV** disease shows a 3.46× hazard vs Stage I (p=0.059), with wide confidence intervals reflecting the small Stage IV sample (n=20). The log-rank test across all stages is highly significant (p=2×10⁻⁶).

---

## Figures

| Figure | Description |
|---|---|
| `01_age_by_stage.png` | Age at diagnosis distribution by cancer stage |
| `02_receptor_status.png` | ER, PR, HER2 receptor status proportions |
| `03_followup_distribution.png` | Patient follow-up time distribution |
| `04_km_by_stage.png` | Kaplan-Meier curves by cancer stage (log-rank p=2×10⁻⁶) |
| `05_km_by_er_status.png` | Kaplan-Meier curves by ER status (log-rank p=0.02) |
| `06_km_by_her2_status.png` | Kaplan-Meier curves by HER2 status (log-rank p=0.005) |
| `07_cox_forest_plot.png` | Multivariate Cox model forest plot |

---

## Dataset

- **Source:** The Cancer Genome Atlas (TCGA) — Breast Cancer (BRCA) cohort
- **Access:** Public, via `TCGAbiolinks` R/Bioconductor package
- **Patients:** 1,035 after quality filtering (66 removed: missing survival time or OS_time ≤ 0)
- **Events (deaths):** 103 (9.9% mortality in follow-up period)
- **Follow-up:** Median 463 days, range 1–7,067 days
- **Multivariate cohort:** 666 patients with complete data across all covariates

---

## Methods

### Data acquisition
Clinical data downloaded from TCGA-BRCA using `TCGAbiolinks::GDCquery()`. Two metadata header rows injected by TCGA were removed before analysis.

### Variables used
- `OS_time`: Days from diagnosis to death (deceased) or last contact (censored)
- `OS_event`: 1 = deceased, 0 = censored
- `age_at_diagnosis`: Continuous, standardised (z-scored) for Cox model
- `ajcc_pathologic_tumor_stage`: Collapsed to Stage I–IV
- `er_status_by_ihc`, `her2_status_by_ihc`: Positive / Negative

### Survival analysis
- **Kaplan-Meier curves** with log-rank tests (`survival`, `survminer`)
- **Cox proportional hazards model** — multivariate, all covariates simultaneously (`coxph`)
- **Forest plot** of hazard ratios with 95% CI (`ggplot2`)

---

## Repository Structure

```
tcga-brca-survival/
├── analysis.R          # Full analysis script (data pull → Cox model)
├── figures/            # All output plots
│   ├── 01_age_by_stage.png
│   ├── 02_receptor_status.png
│   ├── 03_followup_distribution.png
│   ├── 04_km_by_stage.png
│   ├── 05_km_by_er_status.png
│   ├── 06_km_by_her2_status.png
│   └── 07_cox_forest_plot.png
└── README.md
```

---

## Requirements

```r
# Bioconductor
BiocManager::install("TCGAbiolinks")

# CRAN
install.packages(c("survival", "survminer", "dplyr", "tidyr", "ggplot2"))
```

R version 4.3+ recommended.

---

## How to reproduce

```r
# 1. Clone the repo
# 2. Open analysis.R in RStudio
# 3. Run top to bottom — data downloads automatically from TCGA (~200MB)
# 4. All figures save to figures/
```

---

## Author

**Anshika Gupta**  
M.Sc. Bioinformatics and Genomics, Gautam Buddha University  
[LinkedIn](https://www.linkedin.com/in/anshika-gupta-0162886280) · anshikaaggupta@gmail.com

---

## License

Data: TCGA data is publicly available under TCGA data use policy.  
Code: MIT License — free to use and adapt with attribution.
