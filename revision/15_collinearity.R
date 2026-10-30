# ============================================================
# 15_collinearity.R — Script G
#  Collinearity diagnostics among the 15 clustering variables (R1-M1)
#  Correlation matrix (complete cases + best imputation), VIF
# Output: S11 Fig (heatmap), S11 Table (VIF)
# ============================================================
cat("\n=== 15: Collinearity ===\n")
X  <- as.data.frame(REF_SCALED)
R  <- cor(X, method = "spearman")
vif <- sapply(names(X), function(v) 1 / (1 - summary(lm(reformulate(setdiff(names(X), v), v), data = X))$r.squared))
vif_df <- data.frame(Variable = names(vif), VIF = round(vif, 2),
                     Max_abs_r = round(sapply(names(X), function(v) max(abs(R[v, setdiff(names(X), v)]))), 2))
print(vif_df); save_csv(vif_df, "TableS11_VIF.csv")
save_csv(cbind(Variable = rownames(R), round(as.data.frame(R), 2)), "TableS11_correlation_matrix.csv")
pairs_hi <- which(abs(R) > 0.7 & upper.tri(R), arr.ind = TRUE)
cat("|r| > 0.7 pairs:\n"); print(data.frame(A = rownames(R)[pairs_hi[,1]], B = colnames(R)[pairs_hi[,2]], r = round(R[pairs_hi], 2)))

library(ggplot2)
long <- as.data.frame(as.table(R)); names(long) <- c("A", "B", "r")
long$A <- factor(long$A, levels = rownames(R)); long$B <- factor(long$B, levels = rev(rownames(R)))
p <- ggplot(long, aes(A, B, fill = r)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", limits = c(-1, 1), name = "Spearman r") +
  labs(x = NULL, y = NULL) + theme_rev + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(REV_DIR, "FigS11_correlation_heatmap.png"), p, width = 160, height = 150, units = "mm", dpi = 300)
ggsave(file.path(REV_DIR, "FigS11_correlation_heatmap.tiff"), p, width = 160, height = 150, units = "mm", dpi = 300, compression = "lzw")
cat("✓ FigS11 saved\n")
