# COMISA Phenotyping Study — Analysis Pipeline

**Muroi K, et al.** "Phenotypic heterogeneity in comorbid insomnia and sleep apnea identified by polysomnographic and psychological variables"  
Registry: [jRCT1030250339](https://jrct.niph.go.jp/en-latest-detail/jRCT1030250339)

---

## Overview

This repository contains the complete R analysis pipeline for a retrospective multi-center cohort study identifying phenotypes of comorbid insomnia and sleep apnea (COMISA) using unsupervised machine learning (UMAP + hierarchical clustering) applied to real-world sleep clinic data (N = 599).

**Key finding:** 84.8% (95% CI 72.6–92.5%) of COMISA patients (AHI ≥ 15 + insomnia diagnosis) clustered with an insomnia-predominant phenotype characterised by elevated dysfunctional beliefs, sleep reactivity, and hyperarousal.

---

## Repository structure

```
comisa_pipeline/
├── 00_run_all.R          Master script — sources all steps in order
├── 01_setup.R            Package loading and global configuration
├── 02_preprocess.R       Data loading, age adjustment, MICE imputation
├── 03_clustering.R       UMAP grid search, hierarchical clustering, DBSCAN reference
├── 04_validation.R       Expert diagnosis validation, ROC analysis (Figure 4)
├── 05_stability.R        Consensus clustering, ARI, Jaccard stability
├── 06_figures.R          Publication figures 1–3, Supplementary Figures S1 and S2
├── 07_tables.R           Tables 1, 2, and Supplementary Table S2
└── README.md
```

---

## Requirements

- **R** ≥ 4.2.0
- Packages (auto-installed by `01_setup.R`):

| Category | Packages |
|---|---|
| Data manipulation | dplyr, tidyr, tidyselect |
| Visualisation | ggplot2, viridis, RColorBrewer, gridExtra, cowplot |
| Missing data | mice |
| UMAP | uwot |
| Clustering | cluster, dbscan, factoextra, mclust |
| Cluster evaluation | fpc, clusterSim |
| ROC analysis | pROC |
| Effect sizes | effsize |

---

## Usage

### Quick start

1. Clone the repository
2. Open `01_setup.R` and set the two paths:
   ```r
   DATA_PATH  <- "/path/to/comisa.csv"   # raw data (UTF-8 BOM)
   OUTPUT_DIR <- "/path/to/output"       # all outputs saved here
   ```
3. Run the master script:
   ```r
   source("00_run_all.R")
   ```

### Step-by-step (RStudio)

Each script can be run independently in sequence:
```r
source("01_setup.R")       # ~1 sec
source("02_preprocess.R")  # ~2 min (MICE imputation)
source("03_clustering.R")  # ~30–60 min (UMAP grid search × 5 datasets)
source("04_validation.R")  # ~1 min
source("05_stability.R")   # ~20–40 min (200 bootstrap iterations)
source("06_figures.R")     # ~2 min
source("07_tables.R")      # ~1 min
```

**Important:** Each script depends on objects created by prior scripts. Run in order within the same R session, or use `00_run_all.R`.

---

## Data

The raw data file (`comisa.csv`) contains de-identified records from Mates Sleep Clinic (three sites, Gifu, Japan), collected April 2010–December 2024. Data are not publicly shared due to patient privacy; requests may be directed to the corresponding author.

The file is UTF-8 BOM encoded and includes the following key columns:

| Column | Description |
|---|---|
| `Diagnosis` | Expert clinical diagnosis (Japanese text; "不眠症" = insomnia) |
| `AHI` | Apnoea–hypopnoea index (events/h) |
| `ISI` | Insomnia Severity Index (0–28) |
| `DBAS` | Dysfunctional Beliefs and Attitudes about Sleep (0–160) |
| `FIRST` | Ford Insomnia Response to Stress Test (9–36) |
| `HAS` | Hyperarousal Scale (0–72) |
| `SDS` | Self-rating Depression Scale |
| `JESS` | Japanese version Epworth Sleepiness Scale (0–24) |
| `TST` | Total sleep time (min) |
| `SE` | Sleep efficiency (proportion, 0–1) |
| `WASO` | Wake after sleep onset (min) |
| `SL` | Sleep latency (min) |
| `Arousal_Index` | Arousal index (events/h) |
| `ODI3` | Oxygen desaturation index ≥ 3% (events/h) |
| `Age` | Age (years) |
| `BMI` | Body mass index (kg/m²) |

---

## Methods summary

### Preprocessing
- Age adjustment of PSG variables (TST, SE, WASO, SL) using linear regression residuals when R² > 0.01 and p < 0.05 (following Seol et al., 2025, *npj Digital Medicine*)
- Missing data handled by multiple imputation (MICE, predictive mean matching, m = 5, maxit = 10, seed = 42)

### Clustering
- UMAP dimensionality reduction (grid search over n_neighbors, min_dist, metric)
- Hierarchical clustering on UMAP embedding (Euclidean distance, Ward.D2 linkage)
- Optimal k determined by convergent evidence across three cluster quality indices (silhouette coefficient, Calinski-Harabasz index, Davies-Bouldin index); k=3 showed best values on all three indices
- Best imputed dataset selected by maximum silhouette coefficient
- DBSCAN included as reference comparison (not used for primary results)

### Validation
- Criterion validity: cluster assignment vs expert clinical diagnosis (ROC/AUC)
- Exact binomial 95% CI for COMISA cluster proportions

### Stability
- Bootstrap consensus clustering (n = 200 iterations, 80% subsampling)
- Subsampling ARI (n = 100 iterations)
- Jaccard index by cluster

---

## Output

All outputs are saved to `OUTPUT_DIR` (flat structure):

**Figures** (PNG preview + TIFF 300 DPI for submission):
- `Figure1_COMISA_Phenotypes` — UMAP projections (4-panel)
- `Figure2_ISI_vs_AHI_Scatter` — ISI vs AHI scatter
- `Figure3_Clinical_Variables` — Clinical boxplots (9-panel)
- `Figure4_ROC_Expert_Diagnosis` — ROC curves
- `FigureS1_ROC_AHI_Threshold_Comparison` — Supplementary ROC: AHI ≥15 vs AHI ≥5
- `FigureS2_Clinical_Heatmap` — Supplementary heatmap

**Tables** (CSV):
- `Table1_Sample_Characteristics.csv`
- `Table2_Cluster_Characteristics.csv`
- `TableS2_AHI_Sensitivity_Analysis.csv`
- `FigureS1_AUC_Comparison.csv` — AUC summary with DeLong test

---

## Ethics and registration

- IRB approval: R06-258 (University of Tsukuba Hospital)
- Trial registry: jRCT1030250339 (retrospectively registered)
- Study design: retrospective multi-center cohort, opt-out consent

---

## Citation

> Muroi K, et al. Phenotypic heterogeneity in comorbid insomnia and sleep apnea identified by polysomnographic and psychological variables. (manuscript under review).

---

## License

MIT License — see `LICENSE` for details.
