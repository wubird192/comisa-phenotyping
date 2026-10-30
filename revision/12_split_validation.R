# ============================================================
# 12_split_validation.R — Script D
#  Internal cross-validation of the 3-cluster structure (R1-M3, R2-3, R2-9)
#  (1) Repeated split-half: fit UMAP(ret_model)+HC on half A; project half B
#      with umap_transform; assign B by nearest centroid; ARI vs reference on B.
#      Also: cluster B independently; ARI vs reference on B.
#  (2) Prediction strength (Tibshirani & Walther) on original z-space, k=2..6
# Output: S7 Table
# ============================================================
cat("\n=== 12: Split-sample validation ===\n")
N_SPLIT <- 100
res <- data.frame(iter = 1:N_SPLIT, ARI_projected = NA, ARI_independent = NA, k_indep_by_sil = NA)
set.seed(SEED_BOOT + 10)
for (it in 1:N_SPLIT) {
  a <- sample(N, floor(N / 2)); b <- setdiff(1:N, a)
  mA <- fit_umap(REF_SCALED[a, ], ret_model = TRUE, seed = SEED_UMAP + it)
  clA <- hc_cut(mA$embedding, K_FINAL)
  cent <- sapply(1:K_FINAL, function(k) colMeans(mA$embedding[clA == k, , drop = FALSE]))
  embB <- uwot::umap_transform(REF_SCALED[b, ], mA, verbose = FALSE)
  clB_proj <- apply(embB, 1, function(p) which.min(colSums((cent - p)^2)))
  res$ARI_projected[it] <- ari(REF_CL[b], clB_proj)
  embB2 <- fit_umap(REF_SCALED[b, ], seed = SEED_UMAP + 500 + it)
  hcB <- hclust(dist(embB2), method = "ward.D2")
  res$ARI_independent[it] <- ari(REF_CL[b], cutree(hcB, K_FINAL))
  sils <- sapply(2:6, function(k) mean(cluster::silhouette(cutree(hcB, k), dist(embB2))[, 3]))
  res$k_indep_by_sil[it] <- (2:6)[which.max(sils)]
  if (it %% 10 == 0) cat(sprintf("  %d/%d  ARI proj=%.3f indep=%.3f\n", it, N_SPLIT,
                                  mean(res$ARI_projected[1:it]), mean(res$ARI_independent[1:it])))
}
summ <- data.frame(
  Metric = c("ARI, held-out half projected via UMAP model + nearest centroid",
             "ARI, held-out half clustered independently",
             "Held-out half: silhouette-optimal k == 3 (%)"),
  Mean = c(mean(res$ARI_projected), mean(res$ARI_independent), 100 * mean(res$k_indep_by_sil == 3)),
  SD   = c(sd(res$ARI_projected), sd(res$ARI_independent), NA),
  Q025 = c(quantile(res$ARI_projected, .025), quantile(res$ARI_independent, .025), NA),
  Q975 = c(quantile(res$ARI_projected, .975), quantile(res$ARI_independent, .975), NA))
print(summ, digits = 3)
save_csv(res, "TableS7_split_half_iterations.csv"); save_csv(summ, "TableS7_split_half_summary.csv")
cat("k chosen on held-out halves:\n"); print(table(res$k_indep_by_sil))

# (2) Prediction strength on original z-space (method-agnostic check of k)
ps <- fpc::prediction.strength(REF_SCALED, Gmin = 2, Gmax = 6, M = 50,
                               clustermethod = fpc::hclustCBI, method = "ward.D2")
ps_df <- data.frame(k = 2:6, prediction_strength = ps$mean.pred[2:6])
print(ps_df); save_csv(ps_df, "TableS7_prediction_strength.csv")
