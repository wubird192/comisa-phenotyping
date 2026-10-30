# ============================================================
# 22_pipeline_on_synthetic.R — Run the published pipeline end-to-end on the
#  synthetic dataset and store outputs separately (for GitHub; R1-m4, R2-m22, JR-4)
#  Uses a reduced UMAP grid to keep runtime ~5 min; full grid is optional.
# ============================================================
cat("\n=== 22: Published pipeline on synthetic data ===\n")
stopifnot(file.exists("data/synthetic_comisa.csv"))
.DATA_PATH_REAL <- DATA_PATH; .OUTPUT_DIR_REAL <- OUTPUT_DIR
DATA_PATH  <- "data/synthetic_comisa.csv"
OUTPUT_DIR <- "output_synthetic"; if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)
# Reduced grid for speed (set FULL_GRID <- TRUE to reproduce the manuscript grid)
FULL_GRID <- FALSE
if (!FULL_GRID) { UMAP_N_NEIGHBORS <- c(5, 15); UMAP_MIN_DIST <- c(0.001, 0.05); BOOT_N_ITER <- 50; ARI_N_SUBSAMPLE <- 30 }
# Site column in synthetic data is already named per 08; nothing else to adapt
source("02_preprocess.R"); source("03_clustering.R"); source("04_validation.R"); source("05_stability.R")
source("06_figures.R"); source("07_tables.R")
cat("\nSynthetic run complete → output_synthetic/. Restoring real-data workspace...\n")
DATA_PATH <- .DATA_PATH_REAL; OUTPUT_DIR <- .OUTPUT_DIR_REAL
.ws_real <- file.path(.OUTPUT_DIR_REAL, "pipeline_workspace.RData")
if (file.exists(.ws_real)) {
  load(.ws_real, envir = .GlobalEnv)          # 実データのオブジェクト（comisa, final_result, ...）を読み戻す
  DATA_PATH <- .DATA_PATH_REAL; OUTPUT_DIR <- .OUTPUT_DIR_REAL
  cat("✓ Real-data workspace reloaded from", .ws_real, "\n")
} else {
  cat("WARNING: real workspace not found at", .ws_real, "— objects now refer to SYNTHETIC data. Re-run 00_run_all.R before real-data work.\n")
}
