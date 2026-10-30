# ============================================================
# RUN_ALL_REVISION.R — the original scripts (00–07) and the revision scripts must be in THIS folder (flat layout)
# Run from RStudio:  source("RUN_ALL_REVISION.R")
# ============================================================
# Make this folder the working directory regardless of where RStudio started
.this <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (!is.null(.this)) setwd(.this)
cat("Working directory:", getwd(), "\n")

## ---- Step A (run ONCE, 1–2 h): original pipeline + save workspace ----
source("01_setup.R")   # defines DATA_PATH / OUTPUT_DIR (absolute paths on this Mac)
if (!file.exists(file.path(OUTPUT_DIR, "pipeline_workspace.RData")) &&
    !exists("final_result")) {
  source("00_run_all.R")
  save.image(file.path(OUTPUT_DIR, "pipeline_workspace.RData"))
  cat("\n>>> Workspace saved. Restart R (Session > Restart R) and run this file again for Step B.\n")
} else {
## ---- Step B: revision scripts ----
  source("08_revision_setup.R")
  source("09_cluster_validity.R")
  source("10_embedding_sensitivity.R")
  source("11_imputation_consistency.R")
  source("13_site_era.R")
  source("15_collinearity.R")
  source("16_pairwise_stats.R")
  source("17_relabel_ahi5.R")
  source("12_split_validation.R")            # ~15–30 min
  source("14_age_adjustment_sensitivity.R")  # ~8 min
  source("18_synthetic_data.R")
  cat("\nAll revision scripts finished. Share output/revision/, data/synthetic_comisa.csv and output/FigureS_Consensus_Matrix.png\n")
}
