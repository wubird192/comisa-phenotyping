# COMISA Phenotyping Study — Analysis Pipeline

**Muroi K, Noma H, Maeda S, Ishii A, Tanaka H, Yamamoto H, Seol J, Yanagisawa M, Nakajima S.** "Phenotypic heterogeneity in comorbid insomnia and sleep apnea identified by polysomnographic and psychological variables." *PLOS Digital Health* (manuscript PDIG-D-26-01140, under review; revised version submitted November 2026).
Registry: [jRCT1030250339](https://jrct.niph.go.jp/en-latest-detail/jRCT1030250339)

---

## Overview

This repository contains the complete R analysis pipeline for a retrospective cross-sectional study of 599 consecutive adult patients who underwent diagnostic polysomnography at a Japanese sleep clinic (2010–2024). Fifteen polysomnographic, psychological and demographic variables were embedded with UMAP and clustered with Ward's hierarchical clustering, and the resulting groups were compared with the clinical classification of comorbid insomnia and sleep apnea (COMISA) under a clinical definition (insomnia diagnosis + AHI ≥ 5) and a symptom-based definition (ISI ≥ 15 + AHI ≥ 5).

**Main result (revised manuscript):** the prespecified selection rule chose a two-cluster solution: Cluster 1, low AHI with high insomnia-related psychological burden (n = 247, 41.2%), and Cluster 2, high AHI (n = 352, 58.8%). 89.4% of patients with clinically defined COMISA and 62.6% of those with symptom-defined COMISA were in Cluster 1; 29.6% of patients who presented with suspected sleep apnea had the Cluster 1 profile without an insomnia diagnosis. A group with the Cluster 1 profile was reproduced across resampling, imputations, alternative algorithms and the sleep-apnea stream alone, but internal validity indices did not distinguish the solution from one obtained in simulated data without cluster structure, so the two clusters are reported as the ends of a continuum rather than as discrete types.

> **Note on the original submission.** The originally submitted analysis (June 2026) reported a three-cluster solution. During revision, two data-entry errors in two patients' questionnaire scores were corrected after verification with the clinic; re-applying the unchanged selection rule to the corrected data selected two clusters. The scripts in `revision/` reproduce the corrected analysis and all robustness analyses added at revision; the original pipeline (`00`–`07`) is unchanged except for the corrected input file and reproduces the submitted solution when run on the uncorrected data.

---

## Repository structure

```
comisa-phenotyping/
├── 00_run_all.R                 Master script — sources 01–07 in order
├── 01_setup.R                   Package loading, paths, global configuration (seeds)
├── 02_preprocess.R              Data loading, age adjustment, MICE imputation (m = 5, maxit = 10, seed 42)
├── 03_clustering.R              UMAP grid search (40 settings × 5 imputations), Ward.D2 clustering, selection rule
├── 04_validation.R              Agreement with the clinical classification (contingency tables, classification metrics)
├── 05_stability.R               Consensus clustering (200 subsamples of 80%), subsampling ARI, Jaccard
├── 06_figures.R                 Original figures (superseded by revision/36_main_figures_v4.R)
├── 07_tables.R                  Original tables (superseded by revision/ scripts)
├── revision/                    Analyses added at revision (run after 00–07; see "Order of execution")
│   ├── RUN_ALL_REVISION.R       Sources 01, 00 and 08–18 in order (see "Order of execution")
│   ├── 08_revision_setup.R … 39_save_consensus.R   See the table below (no 25, 27, 30–33 as separate files; 30/31/33 are combined in 30_31_33_rewrite_v1.R)
│   ├── build_app_v3.py          Builds S1 Appendix (Markdown) from the output CSVs
│   ├── assemble_manuscript.py   Assembles the manuscript Markdown from ms_parts/ and the reference list
│   └── make_docx.py             Converts Markdown to PLOS-style .docx
├── data/
│   └── synthetic_comisa.csv     Synthetic dataset (synthpop CART); see "Data"
├── output_synthetic/            Outputs of the full pipeline run on the synthetic dataset
├── LICENSE                      MIT
└── README.md
```

### Scripts added at revision (`revision/`)

Each script reads the output of the preceding steps from `OUTPUT_DIR` (and `OUTPUT_DIR/revision/`) and writes CSV tables and figures there. The purpose, inputs and outputs are also given in the header comment of each script. "Manuscript location" refers to the revised manuscript; reviewer comments addressed by each script (R1-M1 etc.) are listed in the script headers.

| Script | Purpose | Main outputs (`OUTPUT_DIR/revision/` unless stated) | Manuscript location |
|---|---|---|---|
| `08_revision_setup.R` | Loads the workspace saved by the original pipeline (00–07) and defines shared helpers; prerequisite for all later scripts | — | — |
| `09_cluster_validity.R` | Cluster-validity indices (silhouette, Calinski–Harabasz, Davies–Bouldin) for k = 2–6, in the embedding and in the original space, for the final configuration and each imputed dataset | `TableS2_ClusterValidity_by_k_final.csv`, `TableS2_ClusterValidity_by_k_imputations.csv`, `FigS2_ClusterValidity_by_k.*` | Results, "Cluster solution"; Table 4; S1 Appendix C1 |
| `10_embedding_sensitivity.R` | Clustering without UMAP (Ward's in the original space, PCA + Ward's, k-means, GMM) and ARI across the 40-setting UMAP grid | `TableS4_AltMethods_ARI.csv`, `TableS4_AltMethods_COMISA_in_C1.csv`, `TableS3_UMAP_grid_ARI.csv`, `FigS3_UMAP_grid_ARI.*` | Table 4; S1 Appendix C2, C5 |
| `11_imputation_consistency.R` | Agreement of cluster solutions across the five imputed datasets (own best setting and final setting); observed vs imputed distributions | `TableS6_Imputation_ARI_summary.csv`, `TableS6_Imputation_ARI_pairwise.csv`, `FigS6_Imputation_densities.*` | Table 4; S1 Appendix A6–A7 |
| `12_split_validation.R` | Repeated split-half replication (projection and independent re-clustering of the held-out half) and prediction strength for k = 2–6 | `TableS7_split_half_summary.csv`, `TableS7_prediction_strength.csv` | Table 4; S1 Appendix C4 |
| `13_site_era.R` | Cluster membership by site, calendar period, PSG system and AASM manual version; era-stratified and leave-one-era-out analyses | `TableS9_cluster_by_era3.csv`, `TableS9_cluster_by_device.csv`, `TableS9_cluster_by_aasm.csv`, `TableS9_era_tests.csv`, `TableS9_n_by_year.csv`, `TableS9_missingness_by_year.csv` | Results, "Calendar period, site, and insomnia duration"; S1 Appendix A2–A3, A8 |
| `14_age_adjustment_sensitivity.R` | Pipeline repeated with no age adjustment and with uniform adjustment of all seven PSG variables | `TableS10_age_adjustment_sensitivity.csv` | Table 4; S1 Appendix B1 |
| `15_collinearity.R` | Correlation matrix and variance inflation factors of the 15 clustering variables | `TableS11_correlation_matrix.csv`, `TableS11_VIF.csv`, `FigS11_correlation_heatmap.*` | Methods, "Selection of clustering variables"; S1 Appendix A5 |
| `16_pairwise_stats.R` | Between-cluster comparisons: Kruskal–Wallis/Wilcoxon with BH-FDR, Cohen's d with 95% CI, η² | `TableS12_compact_d_by_pair.csv`, `TableS12_omnibus_by_variable.csv`, `TableS12_pairwise_Tukey_CohenD.csv`, `TableS12_AHI5_vs_AHI15_effect_sizes.csv` | Table 2; S1 Appendix D3, E |
| `17_relabel_ahi5.R` | Rebuilds the clinical classification at AHI ≥ 5 (primary, registered) and AHI ≥ 15 (sensitivity); cross-tabulations, classification metrics and balanced accuracy | `Table2_crosstab_by_threshold.csv`, `Table3_metrics_by_threshold.csv`, `COMISA_distribution_by_threshold.csv` | Table 3; Results; S1 Appendix D1 |
| `18_synthetic_data.R` | Generates the public synthetic dataset with synthpop (CART) and checks it against the real data (marginals, correlations, missingness, no identical records) | `data/synthetic_comisa.csv`, `synthetic_utility_marginals.csv` | Data availability |
| `19_stream_adjusted.R` | Entry stream × cluster; logistic regression of Cluster 1 membership on stream adjusted for ISI, AHI and DBAS (z-scored covariates); cluster proportion by year within the sleep-apnea stream; sleep-apnea stream re-clustered alone | `TableS9_cluster_by_stream.csv`, `TableS9b_*.csv` | Results, "Patients in the sleep apnea stream assigned to Cluster 1"; Table 4; S1 Appendix C3, D4 |
| `20_symptom_based_comisa.R` | Symptom-based definition of COMISA (ISI ≥ 15 with AHI ≥ 5 or ≥ 15), independent of entry stream | `Table_symptomCOMISA_crosstab_AHI5.csv`, `Table_symptomCOMISA_crosstab_AHI15.csv` | Results, "Distribution of COMISA across the clusters"; S1 Appendix D2 |
| `21_recovery_simulation.R` | Cluster-recovery simulation: multivariate-normal data with the observed covariance, no structure or 2–4 planted clusters at 1–6 SD separation; forced-k silhouette distribution under the null | `TableS13_recovery_simulation_summary.csv`, `TableS13_recovery_simulation_raw.csv`, `TableS13_real_centroid_separation.csv`, `FigS13_recovery_simulation.*` | Results, "Robustness and stability"; Table 4; S1 Appendix C7 |
| `22_pipeline_on_synthetic.R` | Runs the published pipeline end to end on the synthetic dataset (reduced UMAP grid, about 5 min) and stores the outputs separately | `output_synthetic/` | Data availability |
| `23_complete_case.R` | Complete-case clustering (n = 446) with the final configuration | `TableS6_complete_case.csv` | Table 4; S1 Appendix A7 |
| `24_appendix_figures.R` | Appendix figures A2 (observed vs imputed densities) and C1 (consensus matrix) | `FigA2_imputation_densities.png`, `FigC1_consensus.png` | S1 Appendix A6, C3 |
| `26_plausibility_check.R` | Range check of every clustering variable against the instrument and PSG ranges; produces the list of values queried with the clinic. (The cell-by-cell comparison of the analysis file with the source spreadsheet was done locally against the clinic's records and is not distributed.) | `clinic_query_list.csv` | Methods, "Data quality" |
| `28_compare_old_new.R` | Cluster assignments and summary CSVs before vs after the data correction (ARI, cross-tabulation, best mapping, cell-level diffs); aggregate outputs only. **Cannot be run on the public data**: it needs the outputs of the pipeline run on the uncorrected analysis file (`output_old/`), which is not distributed. It is included because it produced the correspondence tables (R1, R2) in the response to reviewers | `28_cluster_crosstab_old_new.csv`, `28_cluster_mapping_old_to_new.csv`, `28_compare_summary.csv`, `28_corrected_cases_cluster.csv`, `28_summary_csv_diffs.csv` | Response letter, Additional corrections 22 |
| `29_forced_k3_sensitivity.R` | k = 3 forced on the corrected data with the final configuration; nesting within the two-cluster solution | `29_*.csv` | S1 Appendix C6 |
| `30_31_33_rewrite_v2.R` | Three descriptive analyses in one script, separated by headers (`30_sex_stream`: sex by cluster, entry stream by cluster, sex within the sleep-apnea stream; `31_delong_definitive`: unpaired DeLong test of balanced accuracy, AHI ≥ 5 vs ≥ 15, definitive-diagnosis sets; `33_C1_composition`: diagnostic composition of each cluster). Ends with an automatic check against the values reported in the manuscript | `TableS1_sex_by_cluster.csv`, `TableS1_stream_by_cluster.csv`, `TableS1_sex_within_SAS.csv`, `TableS12_DeLong_def5_vs_def15_definitive.csv`, `Table3A_rowcounts.csv`, `Table3A_rowpct.csv` | Table 2; Table 3A footnote and 3B; Results; S1 Appendix D4c |
| `34_stream_logistic_CI.R` | 95% confidence intervals for the logistic models of script 19 | `TableS9b_stream_logistic_CI.csv` | S1 Appendix D4b |
| `35_supp_figures_v3.R` | Supplementary figures redrawn in journal style: imputation densities, UMAP-grid ARI, consensus matrix, recovery simulation | `figs/FigS6_Imputation_densities_v2.*`, `figs/FigS3_UMAP_grid_ARI_v2.*`, `figs/FigureS_Consensus_Matrix_v2.*`, `figs/FigS13_recovery_simulation_v2.*` | S1 Appendix A6, C2, C3, C7 |
| `36_main_figures_v4.R` | Figures 1–3 (PNG, PDF and 600-dpi TIFF) | `figs/Fig1_UMAP_v3.*`, `figs/Fig2_ISI_vs_AHI_v3.*`, `figs/Fig3_Variables_by_cluster_v3.*` | Figs 1–3 |
| `38_save_imputation_embedding.R` | Saves the mice object and the final UMAP coordinates (individual-level; not distributed) | `OUTPUT_DIR/mice_imputation.rds`, `OUTPUT_DIR/UMAP_Embedding_Final.csv` | — |
| `39_save_consensus.R` | Re-runs the consensus clustering of `05_stability.R` in the revision workspace and saves the 599 × 599 consensus matrix (individual-level; not distributed) | `OUTPUT_DIR/Consensus_Matrix.rds`, `Cluster_Stability_Summary.csv` (unchanged values) | S1 Appendix C3 |

`37` is a three-line snippet that is pasted into the scripts that create the imputation, embedding and consensus objects; it is not a standalone script and is reproduced here for completeness:

```r
saveRDS(imp, "output/mice_imputation.rds")                                # after mice()
write.csv(data.frame(UMAP1 = emb$layout[, 1], UMAP2 = emb$layout[, 2]),
          "output/UMAP_Embedding_Final.csv", row.names = FALSE)           # after the final umap()
saveRDS(cons, "output/Consensus_Matrix.rds")                             # after the consensus matrix is built
```

---

## Order of execution

1. `source("00_run_all.R")` — original pipeline (01–07). Produces the reference solution (`Cluster_Assignments.csv`, `UMAP_Final_Parameters.csv`), the stability outputs and the original tables. The scripts save the workspace that `08_revision_setup.R` reloads.
2. `source("revision/RUN_ALL_REVISION.R")` — sources, in this order, `01_setup.R`, `00_run_all.R`, `08_revision_setup.R`, `09`, `10`, `11`, `13`, `15`, `16`, `17`, `12`, `14`, `18`. `RUN_ALL_REVISION.R` expects the original scripts (00–07) and the revision scripts in the same working directory (flat layout), so run it from a working copy that contains both, or adjust the paths at the top of the file.
3. `19` → `20` → `21` → `22` → `23` → `24` → `26` → `28` → `29` → `30_31_33` → `34` (each depends on earlier outputs; 28 additionally needs the outputs of the pipeline run on the uncorrected data in `output_old/`).
4. `38` and `39` (save intermediate objects), then `35_supp_figures_v3.R` and `36_main_figures_v4.R` (figures).
5. Optional: `python3 -I revision/build_app_v3.py <OUTPUT_DIR/revision> S1_Appendix.md <figs_dir>` to rebuild the S1 Appendix tables from the CSVs.

All random seeds are fixed (`set.seed(42)` for imputation; the seed is reset to 42 for each UMAP configuration; resampling and simulation seeds are set in each script), so rerunning the pipeline reproduces every value reported in the manuscript. Approximate run times: `03_clustering.R` 30–60 min (UMAP grid search × 5 imputations) and `05_stability.R` 20–40 min (200 consensus subsamples); the other original scripts take a few minutes at most. Of the revision scripts, `12` (split-half, 100 splits), `21` (recovery simulation) and `39` (consensus matrix, 20–40 min) take tens of minutes; the rest run in seconds to a few minutes (`22` about 5 min on the reduced grid; `30_31_33` about 4 min).

## Requirements

- **R 4.5.2** (the original analysis used R 4.4.2 with the same package versions; the revision analyses reproduced the submitted solution exactly under 4.5.2 before the data correction)
- Packages (installed automatically by `01_setup.R` if missing), with the versions used:

| Category | Packages |
|---|---|
| Data manipulation | dplyr 1.2.1, tidyr 1.3.2 |
| Visualisation | ggplot2 4.0.3, viridis 0.6.5, cowplot 1.2.0, gridExtra 2.3, RColorBrewer 1.1.3, patchwork 1.3.2 |
| Missing data | mice 3.19.0 |
| UMAP | uwot 0.2.4 |
| Clustering | cluster 2.1.8.2, dbscan 1.2.4, factoextra 1.0.7, mclust 6.1.2 |
| Cluster evaluation | fpc 2.2.14, clusterSim 0.51.6, clue 0.3.66 |
| Agreement metrics | pROC 1.19.0.1, effsize 0.8.1 |
| Synthetic data | synthpop 1.9.3 |

- Python 3 with `python-docx` and pandoc (only for `build_app_v3.py`, `assemble_manuscript.py`, `make_docx.py`)

## Usage

1. Clone the repository.
2. Open `01_setup.R` and set the two paths:
   ```r
   DATA_PATH  <- "/path/to/comisa.csv"        # raw data (UTF-8 BOM); or "data/synthetic_comisa.csv"
   OUTPUT_DIR <- "/path/to/output"            # all outputs are written here (output/ and output/revision/)
   ```
3. Run the steps in the order given above.

To check that the code runs end to end without the patient data, point `DATA_PATH` at `data/synthetic_comisa.csv`. The expected outputs of that run are in `output_synthetic/`.

---

## Data

The individual-level data are patient records obtained under a waiver of consent. The Institutional Review Board of the Tsukuba Clinical Research & Development Organization (T-CReDO; approval R06-258) restricts their public release to protect patient privacy, and the data are owned by Gifu Mates Sleep Clinic. De-identified individual-level data can be made available to qualified researchers on reasonable request, subject to approval by the data custodian and a data-use agreement. Requests should be sent to the Digital Behavioral Science Group, University of Tsukuba (**cbt-bsm@behavior-science.jp**; https://behavior-science.jp), which forwards them to the custodian.

The cohort comprises 597 patients from Gifu Mates Sleep Clinic and 2 from Iwata Mates Sleep Clinic (Gifu and Shizuoka, Japan; April 2010–December 2024), and is therefore effectively single-centre.

### Synthetic dataset

`data/synthetic_comisa.csv` was generated with **synthpop** (CART method) from the corrected analysis file. It reproduces the marginal distributions, the missingness rates and the correlation structure of the original data (maximum absolute difference in pairwise Spearman correlations 0.13, mean 0.04) and contains no record identical to any real patient.

**The synthetic data are a statistical imitation, not the data analysed in the paper.** When the pipeline is run on them, the selection rule chooses a *three*-cluster solution (Euclidean distance, n_neighbors = 5, min_dist = 0.05, imputed dataset 5; silhouette 0.509), and the proportions of COMISA patients in the low-AHI cluster (87.0% under AHI ≥ 5; 80.4% under AHI ≥ 15) differ from those in the manuscript. The synthetic dataset is provided to verify that the code runs end to end, not to reproduce the reported numbers.

### Columns of the analysis file

| Column | Description |
|---|---|
| — | The entry stream is not a separate column: a recorded insomnia diagnosis (`Diagnosis` containing "不眠症") identifies the insomnia stream (n = 113); all other patients entered through the sleep-apnea stream (n = 486) |
| `Diagnosis` | Clinical diagnosis recorded at the first visit (Japanese text; "不眠症" = insomnia disorder) |
| `AHI` | Apnea–hypopnea index (events/h) |
| `ISI` | Insomnia Severity Index (0–28) |
| `DBAS` | Dysfunctional Beliefs and Attitudes about Sleep, DBAS-16 (0–160) |
| `FIRST` | Ford Insomnia Response to Stress Test (9–36) |
| `HAS` | Hyperarousal Scale, 26-item (0–78) |
| `SDS` | Self-rating Depression Scale (20–80) |
| `JESS` | Japanese version of the Epworth Sleepiness Scale (0–24) |
| `TST` | Total sleep time (min) |
| `SE` | Sleep efficiency (proportion, 0–1) |
| `WASO` | Wake after sleep onset (min) |
| `SL` | Sleep latency (min) |
| `Arousal_Index` | Arousal index (events/h) |
| `ODI3` | 3% oxygen desaturation index (events/h) |
| `Age` | Age (years) |
| `BMI` | Body mass index (kg/m²) |
| `Sex` | 女 (female) / 男 (male); descriptive analyses only, not a clustering variable |
| `PSG_Year` | Year of polysomnography; descriptive analyses only |

The patient identifier in `Cluster_Assignments.csv` is the row number of the analysis file.

---

## Methods summary

- **Preprocessing.** TST, SE and WASO were replaced by age-regression residuals (criterion R² > 0.01 and p < 0.05; sleep latency did not meet it). Missing values were imputed with MICE (predictive mean matching, m = 5, maxit = 10, seed 42). Within each imputed dataset, all 15 variables were z-standardised.
- **Clustering.** UMAP (uwot, 500 epochs) over a grid of 40 settings (metric ∈ {Manhattan, Euclidean} × n_neighbors ∈ {5, 10, 15, 20, 30} × min_dist ∈ {0.001, 0.01, 0.05, 0.1}), Ward.D2 clustering of the embedding for k = 2–6, k chosen by maximum silhouette within each setting, settings ranked by 0.6 × silhouette + 0.3 × Calinski–Harabasz/1000 + 0.1 × 1/Davies–Bouldin (the latter two in the original space), and the imputed dataset whose best setting has the highest silhouette taken as reference. On the corrected data this rule selected imputed dataset 4, Manhattan, n_neighbors = 20, min_dist = 0.01, k = 2.
- **Agreement with the clinical classification.** Contingency tables, sensitivity, specificity, PPV, NPV (exact binomial CIs) and balanced accuracy; cluster membership is a binary predictor, so no ROC curves are drawn.
- **Robustness.** Consensus clustering (200 subsamples of 80%), subsampling ARI and Jaccard, agreement across imputations, hyperparameter-grid sensitivity, split-half replication, prediction strength, alternative algorithms without UMAP, forced k = 3, sensitivity to age adjustment, complete-case analysis, stream-restricted clustering, and a recovery simulation on multivariate-normal data with the observed covariance (null and planted clusters).

Full details are in the Methods and S1 Appendix of the manuscript.

---

## Output

All outputs are written to `OUTPUT_DIR` (original pipeline) and `OUTPUT_DIR/revision/` (revision analyses). Key files:

- `Cluster_Assignments.csv`, `UMAP_Final_Parameters.csv`, `UMAP_Embedding_Final.csv`
- `Table1_*`, `Table2_*`, `Table3_*` (manuscript tables; regenerated from the CSVs)
- `TableS1_*` … `TableS13_*` (S1 Appendix tables), `Cluster_Stability_Summary.csv`, `Consensus_Matrix.rds`
- `revision/figs/`: `Fig1_UMAP_v3.*`, `Fig2_ISI_vs_AHI_v3.*`, `Fig3_Variables_by_cluster_v3.*` (PNG, PDF and 600-dpi TIFF), supplementary figures

## What is not in this repository

No patient-level data of any kind are distributed. Specifically, the following are **never** committed (they are listed in `.gitignore`):

- the analysis file `comisa.csv` and the clinic's source spreadsheet;
- `output/`, `output_old/` and any per-patient file produced by the scripts: `Cluster_Assignments.csv`, `UMAP_Embedding_Final.csv`, `mice_imputation.rds`, `Consensus_Matrix.rds`, `28_*_cases_*.csv`, `clinic_query_list.csv`;
- the uncorrected analysis file and its outputs (`output_old/`), which `28_compare_old_new.R` needs.

What is distributed: the code, `data/synthetic_comisa.csv` (synthetic; no real record), and `output_synthetic/` (outputs of the pipeline on the synthetic data). The aggregate tables reported in the manuscript are in the manuscript and its S1 Appendix. Scripts `35`, `36`, `38` and `39` read or write per-patient intermediate files on the authors' machine; they run on the synthetic data in the same way, and their real-data outputs are not distributed.

---

## Ethics and registration

- IRB approval: R06-258 (Tsukuba Clinical Research & Development Organization, University of Tsukuba); individual consent waived; opt-out notice posted on the clinic website since February 2025
- Registry: jRCT1030250339 (retrospectively registered)
- Design: retrospective cross-sectional study of routine clinical records

---

## Citation

> Muroi K, Noma H, Maeda S, Ishii A, Tanaka H, Yamamoto H, Seol J, Yanagisawa M, Nakajima S. Phenotypic heterogeneity in comorbid insomnia and sleep apnea identified by polysomnographic and psychological variables. PLOS Digital Health (under review).

## License

MIT License — see `LICENSE`.
