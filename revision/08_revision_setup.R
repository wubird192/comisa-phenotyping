# ============================================================
# COMISA Phenotyping — REVISION (PDIG-D-26-01140)
# 08_revision_setup.R — Load pipeline workspace + shared helpers
#
# PREREQUISITE (run once in the original pipeline session):
#   source("00_run_all.R")
#   save.image(file.path(OUTPUT_DIR, "pipeline_workspace.RData"))
#
# Then, in a fresh session:
#   source("01_setup.R")            # packages + config
#   source("revision/08_revision_setup.R")
#   source("revision/09_cluster_validity.R")   # etc.
# ============================================================

cat("\n=== 08: Revision setup ===\n")

WS_PATH <- file.path(OUTPUT_DIR, "pipeline_workspace.RData")
if (!exists("final_result")) {
  stopifnot(file.exists(WS_PATH))
  load(WS_PATH)
  cat("✓ Workspace loaded from", WS_PATH, "\n")
}

REV_DIR <- file.path(OUTPUT_DIR, "revision")
if (!dir.exists(REV_DIR)) dir.create(REV_DIR, recursive = TRUE)

# ---- USER CONFIG: column names for site / date (Script 13) ----
SITE_COL     <- "Site"             # G = Gifu, I = Iwata  (rename from 施設コード)
YEAR_COL     <- "PSG_Year"         # integer calendar year (rename from PSG実施年)
SEX_COL      <- "Sex"              # 1 = male per 07_tables.R
ID_COL       <- "ID"               # "IS-1", "SA-1", ... → intake stream derived from the prefix
STREAM_COL   <- "Stream"           # created below from ID_COL; no need to add it to the CSV
DUR_COL      <- "Duration_Months"  # insomnia duration, IS stream only (rename from 罹病期間（月）)
DUR_CERT_COL <- "Duration_Certain" # 明確 / 不明確 (rename from 罹病期間 明確／不明確)

# ---- Derive intake stream from ID prefix (IS = insomnia-based, SA = OSA-based) ----
if (ID_COL %in% names(comisa_raw) && !(STREAM_COL %in% names(comisa_raw))) {
  comisa_raw[[STREAM_COL]] <- ifelse(grepl("^IS", comisa_raw[[ID_COL]]), "IS",
                               ifelse(grepl("^SA", comisa_raw[[ID_COL]]), "SAS", NA))
  cat("Stream derived from ID prefix:\n"); print(table(comisa_raw[[STREAM_COL]], useNA = "ifany"))
}

# ---- Reference objects ----
REF_CL      <- final_result$clusters            # final 3-cluster assignment
REF_SCALED  <- final_result$scaled_data         # z-standardized 15 vars (best imputation)
REF_PARAMS  <- final_result$best_params         # metric, n_neighbors, min_dist, k
K_FINAL     <- length(unique(REF_CL))
N           <- nrow(REF_SCALED)
ari         <- mclust::adjustedRandIndex

cat(sprintf("Reference: n=%d, k=%d, UMAP(metric=%s, nn=%d, min_dist=%.3f), imputation=%d\n",
            N, K_FINAL, REF_PARAMS$metric, REF_PARAMS$n_neighbors,
            REF_PARAMS$min_dist, best_imp))

# ---- Helpers ----
fit_umap <- function(X, params = REF_PARAMS, seed = SEED_UMAP,
                     ret_model = FALSE, n_epochs = 500) {
  set.seed(seed)
  uwot::umap(X, n_components = 2, n_neighbors = params$n_neighbors,
             min_dist = params$min_dist, metric = params$metric,
             n_epochs = n_epochs, ret_model = ret_model, verbose = FALSE)
}

hc_cut <- function(emb, k = K_FINAL) {
  cutree(hclust(dist(emb, method = "euclidean"), method = "ward.D2"), k = k)
}

# Cluster validity indices. Silhouette/DB on the space used for clustering
# (`space`), CH on the same space. `orig` = original standardized space, so
# indices can be reported in both spaces (embedding vs original features).
cvi <- function(space, cl) {
  dm <- dist(space)
  c(silhouette = mean(cluster::silhouette(cl, dm)[, 3]),
    CH         = fpc::calinhara(space, cl),
    DB         = clusterSim::index.DB(space, cl)$DB)
}

# Relabel `cl` to best match `ref` (Hungarian) — for cross-tabs only; ARI is label-invariant
match_labels <- function(ref, cl) {
  tab  <- table(ref, cl)
  cost <- max(tab) - tab
  perm <- clue::solve_LSAP(matrix(cost, nrow(tab), ncol(tab)))
  map  <- setNames(seq_len(nrow(tab)), colnames(tab)[perm])
  as.integer(map[as.character(cl)])
}

save_csv <- function(df, name) {
  write.csv(df, file.path(REV_DIR, name), row.names = FALSE)
  cat("✓", name, "saved\n")
}

theme_rev <- ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", size = 10))

if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")
cat("✓ Revision helpers ready. Outputs →", REV_DIR, "\n\n")
