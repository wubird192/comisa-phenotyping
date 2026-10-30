# ============================================================
# 24_appendix_figures.R — S1 Appendix Figure A2 / C1 を誌面用に再描画
# Run after: source("01_setup.R"); source("08_revision_setup.R")
# 出力：output/revision/FigA2_imputation_densities.png / FigC1_consensus.png
# （集計済みの図のみ。個人単位データは出力しない）
# ============================================================
library(ggplot2)
# ---- Figure C1: consensus matrix（タイトルなし・凡例付き） ----
stopifnot(exists("consensus_mat"))
ord <- order(REF_CL); m <- consensus_mat[ord, ord]
sz  <- as.vector(table(REF_CL)); K <- length(sz); b <- cumsum(sz)[seq_len(K - 1)] + 0.5  # k-generic (2026-10-18)
png(file.path(REV_DIR, "FigC1_consensus.png"), width = 160, height = 140, units = "mm", res = 300)
layout(matrix(1:2, 1), widths = c(5, 1))
par(mar = c(2, 2, 1, 0.5), family = "sans")
image(seq_len(N), seq_len(N), m, col = colorRampPalette(c("white", "#374E55"))(100),
      zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
abline(v = b, h = b, col = "grey40", lwd = 0.6)
mid <- cumsum(sz) - sz / 2
axis(1, at = mid, labels = paste("Cluster", seq_len(K)), tick = FALSE, line = -0.8, cex.axis = 0.8)
axis(2, at = mid, labels = paste("Cluster", seq_len(K)), tick = FALSE, line = -0.8, cex.axis = 0.8)
par(mar = c(6, 1, 6, 3))
image(1, seq(0, 1, length.out = 100), matrix(seq(0, 1, length.out = 100), 1),
      col = colorRampPalette(c("white", "#374E55"))(100), axes = FALSE, xlab = "", ylab = "")
axis(4, las = 1, cex.axis = 0.7); mtext("Proportion co-clustered", side = 4, line = 2, cex = 0.7)
dev.off()
# ---- Figure A2: observed vs imputed densities ----
vars_imputed <- names(which(colSums(is.na(imputed$data)) > 0))
long <- do.call(rbind, lapply(vars_imputed, function(v) {
  obs <- imputed$data[[v]]; miss <- is.na(obs)
  rbind(data.frame(Variable = v, Source = "Observed", Set = "obs", Value = obs[!miss]),
        do.call(rbind, lapply(1:MICE_M, function(i)
          data.frame(Variable = v, Source = "Imputed", Set = paste0("imp", i),
                     Value = mice::complete(imputed, i)[[v]][miss]))))
}))
p <- ggplot(long, aes(Value, group = Set, colour = Source, linewidth = Source)) +
  geom_density() +
  scale_colour_manual(values = c(Observed = "#374E55", Imputed = "#DF8F44")) +
  scale_linewidth_manual(values = c(Observed = 0.8, Imputed = 0.35)) +
  facet_wrap(~Variable, scales = "free", nrow = 2) +
  labs(x = NULL, y = "Density", colour = NULL, linewidth = NULL) +
  theme_classic(base_size = 9) + theme(legend.position = "top", strip.background = element_blank())
ggsave(file.path(REV_DIR, "FigA2_imputation_densities.png"), p, width = 180, height = 90, units = "mm", dpi = 300)
cat("✓ FigA2 / FigC1 saved\n")
