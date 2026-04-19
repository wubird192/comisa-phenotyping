# ============================================================
# COMISA Phenotyping Study
# 05_stability.R — Cluster stability evaluation:
#                   (1) Consensus clustering (bootstrap, n=200)
#                   (2) Subsampling ARI
#                   (3) Jaccard index by cluster
# ============================================================

cat("\n============================================================\n")
cat("  05: Cluster Stability Analysis\n")
cat("============================================================\n\n")

best_scaled  <- final_result$scaled_data
best_cl      <- final_result$clusters
n_samp       <- nrow(best_scaled)
n_cl         <- length(unique(best_cl))
best_params  <- final_result$best_params

cat(sprintf("n = %d  k = %d\n", n_samp, n_cl))
cat(sprintf("UMAP params: metric=%s  n_neighbors=%d  min_dist=%.3f\n\n",
            best_params$metric, best_params$n_neighbors, best_params$min_dist))

# ============================================================
# 1. CONSENSUS CLUSTERING (bootstrap, n = 200, subsample = 80%)
# ============================================================

cat("=== 1. Bootstrap consensus clustering ===\n")
cat(sprintf("  Iterations: %d  Subsample: %.0f%%\n\n",
            BOOT_N_ITER, 100 * BOOT_SUBSAMPLE))

n_sub              <- floor(n_samp * BOOT_SUBSAMPLE)
co_cluster_mat     <- matrix(0L, n_samp, n_samp)
pair_count_mat     <- matrix(0L, n_samp, n_samp)

set.seed(SEED_BOOT)
pb <- txtProgressBar(min = 0, max = BOOT_N_ITER, style = 3)

for (iter in seq_len(BOOT_N_ITER)) {
  idx    <- sample(seq_len(n_samp), n_sub, replace = FALSE)
  Xsub   <- best_scaled[idx, ]

  set.seed(SEED_UMAP + iter)
  usub   <- tryCatch(
    umap(Xsub, n_components = 2,
         n_neighbors = best_params$n_neighbors,
         min_dist    = best_params$min_dist,
         metric      = best_params$metric,
         n_epochs    = 200, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(usub)) { setTxtProgressBar(pb, iter); next }

  dsub  <- dist(usub, method = "euclidean")
  hcsub <- hclust(dsub, method = "ward.D2")
  clsub <- cutree(hcsub, k = n_cl)

  # Vectorised co-occurrence update
  for (i in seq_len(length(idx) - 1)) {
    for (j in (i + 1):length(idx)) {
      ii <- idx[i]; jj <- idx[j]
      pair_count_mat[ii, jj] <- pair_count_mat[ii, jj] + 1L
      pair_count_mat[jj, ii] <- pair_count_mat[jj, ii] + 1L
      if (clsub[i] == clsub[j]) {
        co_cluster_mat[ii, jj] <- co_cluster_mat[ii, jj] + 1L
        co_cluster_mat[jj, ii] <- co_cluster_mat[jj, ii] + 1L
      }
    }
  }
  setTxtProgressBar(pb, iter)
}
close(pb)

# Consensus matrix (proportion of co-clustering)
consensus_mat <- ifelse(pair_count_mat > 0,
                        co_cluster_mat / pair_count_mat, 0)

# Per-sample consensus score = mean of within-cluster values
consensus_scores <- numeric(n_cl)
for (k in seq_len(n_cl)) {
  idx_k <- which(best_cl == k)
  if (length(idx_k) < 2) next
  sub_mat <- consensus_mat[idx_k, idx_k]
  consensus_scores[k] <- mean(sub_mat[upper.tri(sub_mat)])
}

cat(sprintf("\nConsensus scores (within-cluster co-occurrence):\n"))
for (k in seq_len(n_cl)) {
  cat(sprintf("  Cluster %d: %.3f\n", k, consensus_scores[k]))
}
cat(sprintf("  Overall:   %.3f\n\n", mean(consensus_scores)))

# Visualise consensus matrix using ggplot2 (avoids pheatmap annotation_colors
# bug which causes "subscript out of bounds" on some R/pheatmap versions).
#
# Rows/columns are sorted by cluster assignment so that the block structure
# is visible. A colour sidebar strip is added via a separate panel.

ord      <- order(best_cl)
mat_ord  <- consensus_mat[ord, ord]
cl_ord   <- factor(best_cl[ord], levels = seq_len(n_cl))
n_obs    <- length(cl_ord)

# Long-format data frame for geom_tile
cm_long <- data.frame(
  x     = rep(seq_len(n_obs), times = n_obs),
  y     = rep(seq_len(n_obs), each  = n_obs),
  value = as.vector(mat_ord)
)

# Cluster strip (left sidebar)
strip_df <- data.frame(
  x       = seq_len(n_obs),
  Cluster = cl_ord
)

cl_pal <- CLUSTER_COLORS[seq_len(n_cl)]
names(cl_pal) <- as.character(seq_len(n_cl))

p_heat <- ggplot(cm_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("white", "#2166AC"))(100),
    limits  = c(0, 1),
    name    = "Co-cluster\nproportion"
  ) +
  labs(title = "Consensus Matrix (200 bootstrap iterations)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 10)
  )

p_strip <- ggplot(strip_df, aes(x = x, y = 1, fill = Cluster)) +
  geom_tile() +
  scale_fill_manual(values = cl_pal, name = "Cluster") +
  theme_void() +
  theme(legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.text     = element_text(size = 8))

# Stack: cluster strip (thin) on top of heatmap
p_consensus <- plot_grid(
  p_strip, p_heat,
  ncol        = 1,
  rel_heights = c(0.06, 0.94),
  align       = "v",
  axis        = "lr"
)

ggsave(file.path(OUTPUT_DIR, "FigureS_Consensus_Matrix.png"),
       p_consensus,
       width = 160, height = 155, units = "mm", dpi = 300)
cat("\u2713 FigureS_Consensus_Matrix.png saved\n\n")

# ============================================================
# 2. SUBSAMPLING ARI STABILITY
# ============================================================

cat("=== 2. Subsampling ARI ===\n")
cat(sprintf("  Iterations: %d  Subsample: %.0f%%\n\n",
            ARI_N_SUBSAMPLE, 100 * BOOT_SUBSAMPLE))

ari_values <- numeric(ARI_N_SUBSAMPLE)

set.seed(SEED_BOOT + 1)
pb2 <- txtProgressBar(min = 0, max = ARI_N_SUBSAMPLE, style = 3)

for (iter in seq_len(ARI_N_SUBSAMPLE)) {
  idx   <- sample(seq_len(n_samp), n_sub, replace = FALSE)
  Xsub  <- best_scaled[idx, ]

  set.seed(SEED_UMAP + 1000 + iter)
  usub  <- tryCatch(
    umap(Xsub, n_components = 2,
         n_neighbors = best_params$n_neighbors,
         min_dist    = best_params$min_dist,
         metric      = best_params$metric,
         n_epochs    = 200, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(usub)) { setTxtProgressBar(pb2, iter); next }

  dsub   <- dist(usub, method = "euclidean")
  hcsub  <- hclust(dsub, method = "ward.D2")
  clsub  <- cutree(hcsub, k = n_cl)
  ari_values[iter] <- adjustedRandIndex(best_cl[idx], clsub)
  setTxtProgressBar(pb2, iter)
}
close(pb2)

mean_ari <- mean(ari_values, na.rm = TRUE)
sd_ari   <- sd(ari_values,   na.rm = TRUE)
ari_interp <- if (mean_ari >= 0.80) "Excellent" else
              if (mean_ari >= 0.60) "Good"       else "Moderate"

cat(sprintf("\nARI: mean=%.3f  SD=%.3f  [%s]\n\n",
            mean_ari, sd_ari, ari_interp))

# ARI distribution plot
ari_df <- data.frame(ARI = ari_values[!is.na(ari_values)])
p_ari  <- ggplot(ari_df, aes(x = ARI)) +
  geom_histogram(binwidth = 0.05, fill = "#56B4E9", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean_ari, linetype = "dashed",
             color = "#D55E00", linewidth = 1) +
  annotate("text", x = mean_ari + 0.02, y = Inf,
           label = sprintf("mean = %.3f", mean_ari),
           vjust = 1.5, hjust = 0, size = 3.5, color = "#D55E00") +
  labs(title = "Subsampling ARI Distribution",
       x = "Adjusted Rand Index", y = "Count") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "FigureS_ARI_Distribution.png"), p_ari,
       width = 120, height = 100, units = "mm", dpi = 300)
cat("✓ FigureS_ARI_Distribution.png saved\n\n")

# ============================================================
# 3. JACCARD INDEX BY CLUSTER
# ============================================================

cat("=== 3. Jaccard index by cluster ===\n")

jaccard_by_cl <- matrix(NA_real_, nrow = ARI_N_SUBSAMPLE, ncol = n_cl)
set.seed(SEED_BOOT + 2)
pb3 <- txtProgressBar(min = 0, max = ARI_N_SUBSAMPLE, style = 3)

for (iter in seq_len(ARI_N_SUBSAMPLE)) {
  idx   <- sample(seq_len(n_samp), n_sub, replace = FALSE)
  Xsub  <- best_scaled[idx, ]
  ref   <- best_cl[idx]

  set.seed(SEED_UMAP + 2000 + iter)
  usub  <- tryCatch(
    umap(Xsub, n_components = 2,
         n_neighbors = best_params$n_neighbors,
         min_dist    = best_params$min_dist,
         metric      = best_params$metric,
         n_epochs    = 200, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(usub)) { setTxtProgressBar(pb3, iter); next }

  dsub   <- dist(usub, method = "euclidean")
  hcsub  <- hclust(dsub, method = "ward.D2")
  clsub  <- cutree(hcsub, k = n_cl)

  # Jaccard per cluster (Hungarian matching via clue::solve_LSAP)
  # solve_LSAP *minimises* cost and requires nonneg entries.
  # We convert the co-occurrence table to a cost matrix: cost = max - count.
  co_tab  <- table(ref, clsub)
  cost_mat <- matrix(max(co_tab) - co_tab,
                     nrow = nrow(co_tab), ncol = ncol(co_tab))
  perm    <- clue::solve_LSAP(cost_mat)
  matched <- as.integer(perm[seq_len(n_cl)])

  for (k in seq_len(n_cl)) {
    a    <- as.integer(ref == k)
    b    <- as.integer(clsub == matched[k])
    tp   <- sum(a == 1 & b == 1)
    fp   <- sum(a == 0 & b == 1)
    fn   <- sum(a == 1 & b == 0)
    jaccard_by_cl[iter, k] <- tp / (tp + fp + fn)
  }
  setTxtProgressBar(pb3, iter)
}
close(pb3)

# clue not always available — fall back to no matching if needed
if (any(is.na(jaccard_by_cl))) {
  cat("  Note: clue package unavailable for Hungarian matching;\n")
  cat("  using direct cluster index (valid when cluster labels are stable)\n")
  set.seed(SEED_BOOT + 3)
  for (iter in seq_len(ARI_N_SUBSAMPLE)) {
    idx   <- sample(seq_len(n_samp), n_sub, replace = FALSE)
    Xsub  <- best_scaled[idx, ]
    ref   <- best_cl[idx]

    set.seed(SEED_UMAP + 3000 + iter)
    usub  <- tryCatch(
      umap(Xsub, n_components = 2,
           n_neighbors = best_params$n_neighbors,
           min_dist    = best_params$min_dist,
           metric      = best_params$metric,
           n_epochs    = 200, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(usub)) next
    dsub  <- dist(usub, method = "euclidean")
    hcsub <- hclust(dsub, method = "ward.D2")
    clsub <- cutree(hcsub, k = n_cl)
    for (k in seq_len(n_cl)) {
      a  <- as.integer(ref == k)
      b  <- as.integer(clsub == k)
      tp <- sum(a == 1 & b == 1)
      fp <- sum(a == 0 & b == 1)
      fn <- sum(a == 1 & b == 0)
      jaccard_by_cl[iter, k] <- tp / (tp + fp + fn)
    }
  }
}

mean_jacc <- colMeans(jaccard_by_cl, na.rm = TRUE)
sd_jacc   <- apply(jaccard_by_cl, 2, sd, na.rm = TRUE)
cat(sprintf("\nJaccard index (mean ± SD):\n"))
for (k in seq_len(n_cl)) {
  interp <- if (mean_jacc[k] >= 0.75) "Stable" else "Moderate"
  cat(sprintf("  Cluster %d: %.3f ± %.3f [%s]\n",
              k, mean_jacc[k], sd_jacc[k], interp))
}

# Jaccard box-plot
jacc_long <- data.frame(
  Jaccard = as.vector(jaccard_by_cl),
  Cluster = factor(rep(seq_len(n_cl), each = ARI_N_SUBSAMPLE))
)
p_jacc <- ggplot(jacc_long, aes(x = Cluster, y = Jaccard, fill = Cluster)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = CLUSTER_COLORS[seq_len(n_cl)]) +
  geom_hline(yintercept = 0.75, linetype = "dashed",
             color = "gray40", linewidth = 0.8) +
  labs(title = "Jaccard Index by Cluster", x = "Cluster", y = "Jaccard Index") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "FigureS_Jaccard_by_Cluster.png"), p_jacc,
       width = 120, height = 100, units = "mm", dpi = 300)
cat("✓ FigureS_Jaccard_by_Cluster.png saved\n\n")

# ============================================================
# 4. STABILITY SUMMARY
# ============================================================

stability_summary <- data.frame(
  Metric = c("Consensus — overall",
             paste0("Consensus — Cluster ", seq_len(n_cl)),
             "Subsampling ARI — mean",
             paste0("Jaccard — Cluster ", seq_len(n_cl))),
  Mean   = round(c(mean(consensus_scores), consensus_scores,
                   mean_ari, mean_jacc), 3),
  SD     = round(c(sd(consensus_scores), rep(NA, n_cl),
                   sd_ari, sd_jacc), 3),
  Interpretation = c(
    if (mean(consensus_scores) > 0.90) "Excellent" else "Good",
    ifelse(consensus_scores > 0.90, "Excellent", "Good"),
    ari_interp,
    ifelse(mean_jacc >= 0.75, "Stable", "Moderate")
  )
)

print(stability_summary)
write.csv(stability_summary,
          file.path(OUTPUT_DIR, "Cluster_Stability_Summary.csv"),
          row.names = FALSE)
cat("✓ Cluster_Stability_Summary.csv saved\n\n")

cat("============================================================\n")
cat("  Stability analysis complete. Proceed to 06_figures.R\n")
cat("============================================================\n\n")
