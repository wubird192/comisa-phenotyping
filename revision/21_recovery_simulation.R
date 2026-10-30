# ============================================================
# 21_recovery_simulation.R  (v2) — Cluster-recovery simulation (R2-2/3/4)
#  Does UMAP+HC recover known structure, and does it fabricate clusters when
#  none exist?  n=599, p=15, covariance = that of the real standardized data.
#  Clusters are offset along the real data's leading principal directions
#  (where between-cluster variance actually lives), delta in SD units of the
#  first PC; the real C1–C3 centroid separation is reported on the same scale.
#  No re-standardization after generation (so delta is not shrunk).
# ============================================================
cat("\n=== 21 (v2): Cluster-recovery simulation ===\n")
set.seed(SEED_BOOT + 21)
p <- ncol(REF_SCALED); n <- nrow(REF_SCALED)
Sigma <- cov(REF_SCALED); L <- chol(Sigma)
pc <- prcomp(REF_SCALED); V <- pc$rotation[, 1:3]; sd1 <- pc$sdev[1]
# Real-data reference: centroid separation of C1..C_K projected on PC1–PC2, in units of sd1
cent_real <- t(sapply(1:K_FINAL, function(k) colMeans(pc$x[REF_CL == k, 1:2])))
d_real <- as.vector(dist(cent_real)) / sd1
cat(sprintf("Real data: pairwise centroid separations on PC1–2 = %s (units of PC1 SD = %.2f)\n",
            paste(round(d_real, 2), collapse = ", "), sd1))

sil_k <- function(hc, dm) { s <- sapply(2:6, function(k) mean(cluster::silhouette(cutree(hc, k), dm)[, 3])); (2:6)[which.max(s)] }
simulate_once <- function(k_true, delta) {
  X <- matrix(rnorm(n * p), n, p) %*% L
  if (k_true > 1) {
    lab <- rep(seq_len(k_true), length.out = n)[sample(n)]
    # place centroids on a regular simplex in the PC1–PC2 (and PC3 for k=4) plane, edge length = delta*sd1
    if (k_true == 2) C <- rbind(c(-.5, 0, 0), c(.5, 0, 0))
    if (k_true == 3) C <- rbind(c(-.5, -sqrt(3)/6, 0), c(.5, -sqrt(3)/6, 0), c(0, sqrt(3)/3, 0))
    if (k_true == 4) C <- rbind(c(-.5, -sqrt(3)/6, -sqrt(6)/12), c(.5, -sqrt(3)/6, -sqrt(6)/12), c(0, sqrt(3)/3, -sqrt(6)/12), c(0, 0, sqrt(6)/4))
    shift <- (C * delta * sd1) %*% t(V)       # k x p
    X <- X + shift[lab, ]
  } else lab <- rep(1L, n)
  out <- list()
  for (m in c("UMAP_HC", "HC_orig", "PCA2_HC")) {
    emb <- switch(m, UMAP_HC = fit_umap(X, seed = sample.int(1e6, 1), n_epochs = 200), HC_orig = X, PCA2_HC = prcomp(X)$x[, 1:2])
    dm <- dist(emb); hc <- hclust(dm, method = "ward.D2")
    out[[m]] <- data.frame(Method = m, k_selected = sil_k(hc, dm),
                           ARI_at_true_k = if (k_true > 1) ari(lab, cutree(hc, k_true)) else NA,
                           Sil_k2 = mean(cluster::silhouette(cutree(hc, 2), dm)[, 3]),
                           Sil_k3 = mean(cluster::silhouette(cutree(hc, 3), dm)[, 3]))
  }
  cbind(k_true = k_true, delta = delta, do.call(rbind, out))
}
DELTAS <- c(1, 2, 3, 4, 6)
grid <- expand.grid(k_true = 1:4, delta = DELTAS, rep = 1:20); grid <- grid[!(grid$k_true == 1 & grid$delta > 1), ]
res <- list(); t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  res[[i]] <- cbind(rep = grid$rep[i], simulate_once(grid$k_true[i], grid$delta[i]))
  if (i %% 40 == 0) cat(sprintf("  %d/%d (%.1f min)\n", i, nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, res); rownames(res) <- NULL
save_csv(res, "TableS13_recovery_simulation_raw.csv")
summ <- aggregate(cbind(ARI_at_true_k, Sil_k2, Sil_k3) ~ k_true + delta + Method, res, mean, na.action = na.pass)
summ$Pct_k_correct  <- mapply(function(kt, d, m) 100 * mean(res$k_selected[res$k_true == kt & res$delta == d & res$Method == m] == kt), summ$k_true, summ$delta, summ$Method)
summ$Pct_k3_selected <- mapply(function(kt, d, m) 100 * mean(res$k_selected[res$k_true == kt & res$delta == d & res$Method == m] == 3), summ$k_true, summ$delta, summ$Method)
summ$Pct_k2_selected <- mapply(function(kt, d, m) 100 * mean(res$k_selected[res$k_true == kt & res$delta == d & res$Method == m] == 2), summ$k_true, summ$delta, summ$Method)
summ <- summ[order(summ$Method, summ$k_true, summ$delta), ]
print(summ, digits = 3); save_csv(summ, "TableS13_recovery_simulation_summary.csv")
pairs <- combn(K_FINAL, 2); save_csv(data.frame(Pair = paste0("C", pairs[1, ], "-C", pairs[2, ]), Separation_PC1SD = d_real), "TableS13_real_centroid_separation.csv")
cat("\nNull data (k_true=1): % runs selecting k=2 / k=3, mean silhouette of forced k=2 / k=3:\n"); print(summ[summ$k_true == 1, c("Method", "Pct_k2_selected", "Pct_k3_selected", "Sil_k2", "Sil_k3")])
null_sil2 <- res$Sil_k2[res$k_true == 1 & res$Method == "UMAP_HC"]
cat(sprintf("Null UMAP_HC silhouette at k=2: mean %.3f, 95th percentile %.3f, max %.3f (real data: %.3f)\n", mean(null_sil2), quantile(null_sil2, 0.95), max(null_sil2), final_result$silhouette))
cat(sprintf("Real-data silhouette at k=%d (embedding) = %.3f\n", K_FINAL, final_result$silhouette))

library(ggplot2)
pd <- summ[summ$k_true > 1, ]
vl <- data.frame(x = d_real)
p1 <- ggplot(pd, aes(delta, ARI_at_true_k, colour = Method)) + geom_vline(data = vl, aes(xintercept = x), linetype = "dotted", colour = "grey55") +
  geom_line() + geom_point() + facet_wrap(~k_true, labeller = labeller(k_true = function(x) paste0("true k = ", x))) +
  labs(x = "Between-centroid separation (units of PC1 SD); dotted = real between-cluster separation", y = "ARI vs. true labels") + theme_rev + theme(legend.position = "bottom")
p2 <- ggplot(pd, aes(delta, Pct_k_correct, colour = Method)) + geom_vline(data = vl, aes(xintercept = x), linetype = "dotted", colour = "grey55") +
  geom_line() + geom_point() + facet_wrap(~k_true, labeller = labeller(k_true = function(x) paste0("true k = ", x))) +
  labs(x = "Between-centroid separation (units of PC1 SD)", y = "Silhouette selects true k (%)") + theme_rev + theme(legend.position = "bottom")
g <- cowplot::plot_grid(p1, p2, ncol = 1, labels = c("A", "B"))
ggsave(file.path(REV_DIR, "FigS13_recovery_simulation.png"), g, width = 180, height = 170, units = "mm", dpi = 300)
ggsave(file.path(REV_DIR, "FigS13_recovery_simulation.tiff"), g, width = 180, height = 170, units = "mm", dpi = 300, compression = "lzw")
cat("✓ FigS13 saved\n")
