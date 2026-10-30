# ============================================================
# 09_cluster_validity.R — Script A
# Cluster validity indices for k = 2..6  (R2-4, R1-M3)
#   (a) Final UMAP configuration, best imputation — embedding & original space
#   (b) All 5 imputed datasets at their own best configuration (embedding)
# Output: S2 Table (csv) + S2 Fig (png/tiff)
# ============================================================
cat("\n=== 09: Cluster validity across k ===\n")
K_RANGE <- 2:6

# (a) Final configuration
hc_ref <- hclust(dist(final_result$umap_result), method = "ward.D2")
rows <- lapply(K_RANGE, function(k) {
  cl <- cutree(hc_ref, k = k)
  e  <- cvi(final_result$umap_result, cl)
  o  <- cvi(REF_SCALED, cl)
  data.frame(k = k, n_min_cluster = min(table(cl)),
             Sil_embedding = e["silhouette"], CH_embedding = e["CH"], DB_embedding = e["DB"],
             Sil_original  = o["silhouette"], CH_original  = o["CH"], DB_original  = o["DB"])
})
validity_final <- do.call(rbind, rows); rownames(validity_final) <- NULL
print(round(validity_final, 3))
save_csv(validity_final, "TableS2_ClusterValidity_by_k_final.csv")

# (b) Per imputation (embedding space, each at its own best UMAP params)
rows <- list()
for (imp in seq_len(MICE_M)) {
  cr <- clustering_results[[imp]]
  hc <- hclust(dist(cr$umap_result), method = "ward.D2")
  for (k in K_RANGE) {
    cl <- cutree(hc, k = k); e <- cvi(cr$umap_result, cl)
    rows[[length(rows) + 1]] <- data.frame(Imputation = imp, k = k,
      Silhouette = e["silhouette"], CH = e["CH"], DB = e["DB"])
  }
}
validity_imp <- do.call(rbind, rows); rownames(validity_imp) <- NULL
save_csv(validity_imp, "TableS2_ClusterValidity_by_k_imputations.csv")

# Figure S2: three panels (embedding space, final config), k=3 highlighted
library(ggplot2)
long <- tidyr::pivot_longer(validity_final[, c("k", "Sil_embedding", "CH_embedding", "DB_embedding")],
                            -k, names_to = "Index", values_to = "Value")
long$Index <- factor(long$Index, levels = c("Sil_embedding", "CH_embedding", "DB_embedding"),
                     labels = c("Silhouette (higher = better)",
                                "Calinski\u2013Harabasz (higher = better)",
                                "Davies\u2013Bouldin (lower = better)"))
p <- ggplot(long, aes(k, Value)) +
  geom_vline(xintercept = K_FINAL, linetype = "dotted", colour = "grey55") +
  geom_line(colour = "grey30") + geom_point(size = 2, colour = "grey20") +
  facet_wrap(~Index, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "Number of clusters (k)", y = NULL) + theme_rev
ggsave(file.path(REV_DIR, "FigS2_ClusterValidity_by_k.png"), p, width = 180, height = 60, units = "mm", dpi = 300)
ggsave(file.path(REV_DIR, "FigS2_ClusterValidity_by_k.tiff"), p, width = 180, height = 60, units = "mm", dpi = 300, compression = "lzw")
cat("✓ FigS2 saved\n")
