# ============================================================
# 11_imputation_consistency.R — Script C
#  (1) Cluster agreement across the 5 imputed datasets (R2-8)
#      a. each dataset at its own best UMAP config (as in pipeline)
#      b. each dataset at the FINAL config (isolates imputation effect)
#  (2) Observed vs imputed distributions per variable (R1-M1)
# Output: S6 Table (pairwise ARI), S6 Fig (density plots)
# ============================================================
cat("\n=== 11: Imputation consistency ===\n")
vars_imputed <- names(which(colSums(is.na(imputed$data)) > 0))

# (1a) own best config, k = 3
cl_own <- lapply(seq_len(MICE_M), function(i) {
  cr <- clustering_results[[i]]
  hc_cut(cr$umap_result, K_FINAL)
})
# (1b) final config on every completed dataset
cl_fin <- lapply(seq_len(MICE_M), function(i) {
  X <- scale(mice::complete(imputed, i)); hc_cut(fit_umap(X), K_FINAL)
})
pair_mat <- function(cls) {
  m <- matrix(NA, MICE_M, MICE_M, dimnames = list(1:MICE_M, 1:MICE_M))
  for (i in 1:MICE_M) for (j in 1:MICE_M) m[i, j] <- ari(cls[[i]], cls[[j]])
  m
}
A_own <- pair_mat(cl_own); A_fin <- pair_mat(cl_fin)
cat("Pairwise ARI — own best config:\n"); print(round(A_own, 3))
cat("Pairwise ARI — final config:\n");    print(round(A_fin, 3))
off <- function(m) m[upper.tri(m)]
summ <- data.frame(
  Setting = c("Own best config per imputation", "Final config applied to all"),
  Mean_pairwise_ARI = c(mean(off(A_own)), mean(off(A_fin))),
  Min_pairwise_ARI  = c(min(off(A_own)),  min(off(A_fin))),
  ARI_vs_reference  = c(paste(round(sapply(cl_own, ari, REF_CL), 3), collapse = "/"),
                        paste(round(sapply(cl_fin, ari, REF_CL), 3), collapse = "/")))
print(summ); save_csv(summ, "TableS6_Imputation_ARI_summary.csv")
save_csv(as.data.frame(as.table(A_fin)), "TableS6_Imputation_ARI_pairwise_finalconfig.csv")

# Majority-vote membership across imputations (final config): % identical to reference
matched <- sapply(cl_fin, match_labels, ref = REF_CL)
vote    <- apply(matched, 1, function(r) as.integer(names(which.max(table(r)))))
agree   <- rowMeans(matched == vote)
cat(sprintf("Majority-vote cluster == reference for %.1f%% of participants; %.1f%% unanimous across 5 imputations\n",
            100 * mean(vote == REF_CL), 100 * mean(agree == 1)))
save_csv(data.frame(ID = seq_len(N), Reference = REF_CL, MajorityVote = vote, Agreement = agree),
         "Imputation_membership_by_participant.csv")

# (2) Observed vs imputed densities
library(ggplot2)
long <- do.call(rbind, lapply(vars_imputed, function(v) {
  obs <- imputed$data[[v]]; miss <- is.na(obs)
  rbind(data.frame(Variable = v, Source = "Observed", Value = obs[!miss]),
        do.call(rbind, lapply(1:MICE_M, function(i)
          data.frame(Variable = v, Source = paste0("Imputed ", i),
                     Value = mice::complete(imputed, i)[[v]][miss]))))
}))
p <- ggplot(long, aes(Value, colour = Source, linetype = Source)) +
  geom_density(linewidth = 0.5) +
  scale_colour_manual(values = c(Observed = "black", setNames(rep("#D55E00", MICE_M), paste0("Imputed ", 1:MICE_M)))) +
  scale_linetype_manual(values = c(Observed = "solid", setNames(rep("dashed", MICE_M), paste0("Imputed ", 1:MICE_M)))) +
  facet_wrap(~Variable, scales = "free") + labs(x = NULL, y = "Density") + theme_rev +
  theme(legend.position = "bottom")
ggsave(file.path(REV_DIR, "FigS6_Imputation_densities.png"), p, width = 180, height = 120, units = "mm", dpi = 300)
ggsave(file.path(REV_DIR, "FigS6_Imputation_densities.tiff"), p, width = 180, height = 120, units = "mm", dpi = 300, compression = "lzw")
cat("✓ FigS6 saved\n")
