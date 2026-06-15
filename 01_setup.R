# ============================================================
# COMISA Phenotyping Study
# 01_setup.R — Package installation, loading, and global config
#
# Muroi K, et al. "Phenotypic heterogeneity in comorbid insomnia and sleep apnea
# identified by polysomnographic and psychological variables"
# Manuscript: see citation above
# ============================================================

cat("\n============================================================\n")
cat("  COMISA PHENOTYPING PIPELINE\n")
cat("  01: Setup — packages & configuration\n")
cat("============================================================\n\n")

# ============================================================
# USER CONFIGURATION — edit these paths before running
# ============================================================

# Path to raw data CSV (UTF-8 BOM encoded)
DATA_PATH <- "data/comisa.csv"

# Working directory — all output files will be saved here
OUTPUT_DIR <- "output"

# Random seeds
SEED_MICE  <- 42
SEED_UMAP  <- 42
SEED_BOOT  <- 123

# MICE parameters
MICE_M     <- 5    # number of imputed datasets
MICE_MAXIT <- 10   # maximum iterations

# UMAP grid search ranges
UMAP_N_NEIGHBORS <- c(5, 10, 15, 20, 30)
UMAP_MIN_DIST    <- c(0.001, 0.01, 0.05, 0.1)
UMAP_METRIC      <- c("euclidean", "manhattan")
UMAP_MAX_K       <- 6   # maximum k evaluated in grid search

# Consensus / stability
BOOT_N_ITER        <- 200  # bootstrap iterations
BOOT_SUBSAMPLE     <- 0.8  # subsample proportion
ARI_N_SUBSAMPLE    <- 100  # iterations for ARI stability

# ============================================================
# PACKAGE INSTALLATION & LOADING
# ============================================================

required_packages <- c(
  # Data manipulation
  "dplyr", "tidyr", "tidyselect",
  # Visualisation
  "ggplot2", "viridis", "RColorBrewer",
  "gridExtra", "cowplot", "pheatmap",
  # Missing data
  "mice",
  # UMAP
  "uwot",
  # Clustering
  "cluster", "dbscan", "factoextra", "mclust",
  # Cluster evaluation
  "fpc",         # calinhara()
  "clusterSim",  # index.DB()
  # ROC analysis
  "pROC",
  # Effect sizes (sensitivity analysis)
  "effsize"
)

cat("Loading packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
cat(sprintf("✓ %d packages loaded\n\n", length(required_packages)))

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("✓ Created output directory: %s\n\n", OUTPUT_DIR))
} else {
  cat(sprintf("✓ Output directory: %s\n\n", OUTPUT_DIR))
}

# Colour palette (colour-blind safe)
CLUSTER_COLORS <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73",
                    "4" = "#F0E442", "5" = "#0072B2", "6" = "#D55E00")

cat("============================================================\n")
cat("  Setup complete. Proceed to 02_preprocess.R\n")
cat("============================================================\n\n")
