# ============================================================
# COMISA Phenotyping Study
# 00_run_all.R — Master pipeline
#
# Muroi K, et al. "Phenotypic heterogeneity in comorbid insomnia and sleep apnea
# identified by polysomnographic and psychological variables"
# Target: npj Digital Medicine
# Registry: jRCT1030250339
#
# ── Pipeline steps ────────────────────────────────────────────
#   01_setup.R      Packages + global configuration
#   02_preprocess.R Data, age adjustment (Seol 2025), MICE
#   03_clustering.R UMAP grid search + hierarchical clustering
#                   + DBSCAN reference analysis
#   04_validation.R Expert diagnosis validation + ROC (Figure 4)
#   05_stability.R  Consensus, ARI, Jaccard stability
#   06_figures.R    Figures 1–3 and Supplementary Figure S2
#   07_tables.R     Tables 1, 2, and Supplementary Table S2
#
# ── Before running ───────────────────────────────────────────
#   1. Open 01_setup.R and set DATA_PATH and OUTPUT_DIR
#   2. Run this script: source("00_run_all.R")
#      or in terminal:  Rscript 00_run_all.R
#
# ── Estimated run time ────────────────────────────────────────
#   UMAP grid search × 5 imputed datasets: ~30–60 min
#   Stability analysis (200 bootstrap):    ~20–40 min
#   Total: approximately 1–2 hours
#
# ── Output files ──────────────────────────────────────────────
#   Figures (PNG + TIFF, 300 DPI):
#     Figure1_COMISA_Phenotypes_LDH
#     Figure2_ISI_vs_AHI_Scatter_LDH
#     Figure3_Clinical_Variables_LDH
#     Figure4_ROC_Expert_Diagnosis_FINAL
#     FigureS2_Clinical_Heatmap_LDH
#     FigureS_Consensus_Matrix
#     FigureS_ARI_Distribution
#     FigureS_Jaccard_by_Cluster
#   Tables (CSV):
#     Table1_Sample_Characteristics.csv
#     Table2_Cluster_Characteristics.csv
#     TableS2_AHI_Sensitivity_Analysis.csv
#   Intermediate data:
#     Cluster_Assignments.csv
#     UMAP_Final_Parameters.csv
#     Imputation_Comparison.csv
#     Age_Adjustment_Summary.csv
#     Missing_Data_Summary.csv
#     ROC_Results_Summary.csv
#     Cluster_Stability_Summary.csv
#     COMISA_Cluster_Distribution.csv
#     DBSCAN_Reference_Results.csv
#     ANOVA_Results_Clinical_Variables.csv
# ============================================================

pipeline_start <- Sys.time()

# Detect script directory robustly:
#   - When run via source("00_run_all.R") in RStudio → uses sys.frames
#   - When run via Rscript 00_run_all.R in terminal  → uses commandArgs
#   - Fallback: use current working directory
get_script_dir <- function() {
  # Try sys.frames (source() in RStudio)
  for (i in sys.nframe():1) {
    f <- sys.frame(i)$ofile
    if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f)))
  }
  # Try commandArgs (Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  m    <- regmatches(args, regexpr("(?<=--file=).+", args, perl = TRUE))
  if (length(m) > 0 && nzchar(m[1])) return(dirname(normalizePath(m[1])))
  # Fallback: working directory
  message("⚠ Could not detect script directory; using getwd().\n",
          "  If scripts are not in the working directory, set\n",
          "  script_dir manually below.")
  getwd()
}

script_dir <- get_script_dir()
cat(sprintf("Script directory: %s\n\n", script_dir))

# ── Override here if auto-detection fails ──────────────────
# script_dir <- "/path/to/comisa_pipeline"

source_step <- function(file) {
  path <- file.path(script_dir, file)
  cat(sprintf("\n%s\n  Sourcing: %s\n%s\n",
              strrep("=", 60), file, strrep("=", 60)))
  source(path, echo = FALSE)
}

source_step("01_setup.R")
source_step("02_preprocess.R")
source_step("03_clustering.R")
source_step("04_validation.R")
source_step("05_stability.R")
source_step("06_figures.R")
source_step("07_tables.R")

elapsed <- difftime(Sys.time(), pipeline_start, units = "mins")

cat("\n")
cat(strrep("=", 60), "\n")
cat("  PIPELINE COMPLETE\n")
cat(sprintf("  Total run time: %.1f minutes\n", as.numeric(elapsed)))
cat(sprintf("  Output directory: %s\n", OUTPUT_DIR))
cat(strrep("=", 60), "\n\n")

# Session info for reproducibility
cat("Session info:\n")
print(sessionInfo())
