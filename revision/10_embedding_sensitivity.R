# ============================================================
# 10_embedding_sensitivity.R — Script B
#  (1) Clustering WITHOUT UMAP: HC on original z-space, PCA+HC, k-means, GMM
#      → ARI vs reference, cross-tabs                     (R2-2, R2-13)
#  (2) UMAP hyperparameter grid → ARI of k=3 solution vs reference,
#      using embeddings already stored in final_result$all_results (R2-3)
# Output: S3 Fig (grid heatmap), S4 Table (ARI across methods)
# ============================================================
cat("\n=== 10: Embedding sensitivity ===\n")
X <- REF_SCALED

# ---- (1) Alternative clustering pipelines ----
alt <- list()
alt$HC_original     <- hc_cut(X, K_FINAL)
pca <- prcomp(X, center = FALSE, scale. = FALSE)
cum <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
n80 <- which(cum >= 0.80)[1]
alt$PCA2_HC         <- hc_cut(pca$x[, 1:2], K_FINAL)
alt[[sprintf("PCA%d_HC_80pct", n80)]] <- hc_cut(pca$x[, 1:n80], K_FINAL)
set.seed(SEED_UMAP); alt$kmeans_original <- kmeans(X, K_FINAL, nstart = 50)$cluster
set.seed(SEED_UMAP); alt$GMM_original    <- mclust::Mclust(X, G = K_FINAL, verbose = FALSE)$classification
alt$UMAP_HC_reference <- REF_CL

method_df <- do.call(rbind, lapply(names(alt), function(nm) {
  cl <- alt[[nm]]; e <- cvi(X, cl)
  data.frame(Method = nm, ARI_vs_reference = ari(REF_CL, cl),
             Sil_original = e["silhouette"], CH_original = e["CH"], DB_original = e["DB"],
             Cluster_sizes = paste(sort(table(cl), decreasing = TRUE), collapse = "/"))
}))
rownames(method_df) <- NULL; print(method_df, digits = 3)
save_csv(method_df, "TableS4_AltMethods_ARI.csv")

# Cross-tabs (labels matched to reference) — for the response letter
for (nm in setdiff(names(alt), "UMAP_HC_reference")) {
  cat("\n", nm, "(rows = reference UMAP+HC, cols = matched alt):\n")
  print(table(Reference = REF_CL, Alt = match_labels(REF_CL, alt[[nm]])))
}
# COMISA distribution under each alternative (does the key finding survive?)
comisa_idx <- which(!is.na(vd$Expert_Label) & vd$Expert_Label == "COMISA")
key <- do.call(rbind, lapply(names(alt), function(nm) {
  m <- match_labels(REF_CL, alt[[nm]]); tb <- table(factor(m[comisa_idx], levels = 1:K_FINAL))
  data.frame(Method = nm, N_COMISA = length(comisa_idx),
             Pct_in_matched_C1 = round(100 * tb[1] / length(comisa_idx), 1))
}))
save_csv(key, "TableS4_AltMethods_COMISA_in_C1.csv")

# ---- (2) UMAP grid: ARI of k=3 solution vs reference ----
grid <- do.call(rbind, lapply(final_result$all_results, function(r) {
  if (is.null(r)) return(NULL)
  cl3 <- cutree(r$hclust_object, k = K_FINAL)
  data.frame(metric = r$metric, n_neighbors = r$n_neighbors, min_dist = r$min_dist,
             best_k_by_silhouette = r$best_k, silhouette_at_best_k = r$silhouette,
             composite = r$composite_score,
             ARI_k3_vs_reference = ari(REF_CL, cl3),
             Sil_k3_embedding = mean(cluster::silhouette(cl3, dist(r$umap_result))[, 3]))
}))
save_csv(grid, "TableS3_UMAP_grid_ARI.csv")
cat(sprintf("\nUMAP grid (%d configs): ARI vs reference — median %.3f, min %.3f; best_k==3 in %d/%d\n",
            nrow(grid), median(grid$ARI_k3_vs_reference), min(grid$ARI_k3_vs_reference),
            sum(grid$best_k_by_silhouette == 3), nrow(grid)))

library(ggplot2)
p <- ggplot(grid, aes(factor(n_neighbors), factor(min_dist), fill = ARI_k3_vs_reference)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI_k3_vs_reference)), size = 2.6) +
  facet_wrap(~metric) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), name = "ARI") +
  labs(x = "n_neighbors", y = "min_dist",
       title = "Adjusted Rand index of the 3-cluster solution vs. final configuration") + theme_rev
ggsave(file.path(REV_DIR, "FigS3_UMAP_grid_ARI.png"), p, width = 180, height = 80, units = "mm", dpi = 300)
ggsave(file.path(REV_DIR, "FigS3_UMAP_grid_ARI.tiff"), p, width = 180, height = 80, units = "mm", dpi = 300, compression = "lzw")
cat("✓ FigS3 saved\n")
