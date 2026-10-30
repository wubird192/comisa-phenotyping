# ============================================================
# 16_pairwise_stats.R — Script H
#  (1) Between-cluster comparisons: omnibus (ANOVA/Kruskal), Tukey HSD
#      pairwise p, Cohen's d per pair, eta², BH-FDR   (R1-m2, R2-18, R2-19, R2-m17)
#  (2) COMISA defined by AHI>=5 vs AHI>=15: standardized effect sizes (R1-m2)
#  (3) DeLong test C2 vs C3 for OSA (paired, same reference)   (R1-m2)
# Output: S12 Tables; Table 1 effect-size columns
# ============================================================
cat("\n=== 16: Pairwise statistics ===\n")
cl  <- factor(REF_CL)
dat <- comisa  # raw (non-imputed) values, as in Table 1
VARS <- c("Age", "BMI", "AHI", "ODI3", "Arousal_Index", "TST", "SE", "WASO", "SL",
          "ISI", "DBAS", "FIRST", "HAS", "JESS", "SDS")
VARS <- VARS[VARS %in% names(dat)]
pairs <- combn(levels(cl), 2, simplify = FALSE)

rows <- list(); omni <- list()
for (v in VARS) {
  y <- dat[[v]]; ok <- !is.na(y); yy <- y[ok]; cc <- droplevels(cl[ok])
  ao <- aov(yy ~ cc); tuk <- TukeyHSD(ao)$cc
  eta2 <- summary(ao)[[1]][["Sum Sq"]][1] / sum(summary(ao)[[1]][["Sum Sq"]])
  omni[[v]] <- data.frame(Variable = v, n = length(yy),
    Shapiro_p = tryCatch(shapiro.test(if (length(yy) > 5000) sample(yy, 5000) else yy)$p.value, error = function(e) NA),
    ANOVA_p = summary(ao)[[1]][["Pr(>F)"]][1], Kruskal_p = kruskal.test(yy ~ cc)$p.value, eta2 = eta2)
  for (p in pairs) {
    a <- yy[cc == p[1]]; b <- yy[cc == p[2]]
    d <- effsize::cohen.d(a, b)
    rows[[length(rows) + 1]] <- data.frame(Variable = v, Cluster_A = p[1], Cluster_B = p[2],
      Mean_A = mean(a), SD_A = sd(a), Median_A = median(a), Mean_B = mean(b), SD_B = sd(b), Median_B = median(b),
      Diff = mean(a) - mean(b), Tukey_p = tuk[paste(p[2], p[1], sep = "-"), "p adj"],
      Wilcoxon_p = wilcox.test(a, b)$p.value,
      Cohen_d = d$estimate, d_CI_low = d$conf.int[1], d_CI_high = d$conf.int[2])
  }
}
pw <- do.call(rbind, rows); om <- do.call(rbind, omni); rownames(pw) <- rownames(om) <- NULL
pw$Tukey_p_FDR <- p.adjust(pw$Tukey_p, "BH"); pw$Wilcoxon_p_FDR <- p.adjust(pw$Wilcoxon_p, "BH")
om$ANOVA_p_FDR <- p.adjust(om$ANOVA_p, "BH"); om$Kruskal_p_FDR <- p.adjust(om$Kruskal_p, "BH")
om$Normal <- om$Shapiro_p > 0.05
print(om, digits = 3); save_csv(om, "TableS12_omnibus_by_variable.csv")
save_csv(pw, "TableS12_pairwise_Tukey_CohenD.csv")
# Compact version for the supplement: d (95% CI) per pair, p_FDR
pw$Cell <- sprintf("%.2f (%.2f to %.2f)%s", pw$Cohen_d, pw$d_CI_low, pw$d_CI_high,
                   ifelse(pw$Tukey_p_FDR < 0.001, "***", ifelse(pw$Tukey_p_FDR < 0.01, "**", ifelse(pw$Tukey_p_FDR < 0.05, "*", ""))))
save_csv(pw[, c("Variable", "Cluster_A", "Cluster_B", "Cell")], "TableS12_compact_d_by_pair.csv")

# (2) AHI>=5 vs AHI>=15 COMISA definitions: standardized differences on psych scales
ins <- dat$Diagnosis == "\u4e0d\u7720\u75c7" & !is.na(dat$AHI)
c5  <- dat[ins & dat$AHI >= 5, ]; c15 <- dat[ins & dat$AHI >= 15, ]
c5_15 <- dat[ins & dat$AHI >= 5 & dat$AHI < 15, ]   # non-overlapping band
es <- do.call(rbind, lapply(c("ISI", "DBAS", "FIRST", "HAS", "SDS", "JESS"), function(v) {
  d1 <- effsize::cohen.d(c5[[v]], c15[[v]], na.rm = TRUE)          # overlapping cohorts (descriptive)
  d2 <- effsize::cohen.d(c5_15[[v]], c15[[v]], na.rm = TRUE)       # 5<=AHI<15 vs AHI>=15 (independent)
  data.frame(Variable = v, n_AHI5 = sum(!is.na(c5[[v]])), n_AHI15 = sum(!is.na(c15[[v]])),
             Mean_AHI5 = mean(c5[[v]], na.rm = TRUE), Mean_AHI15 = mean(c15[[v]], na.rm = TRUE),
             d_AHI5_vs_AHI15 = d1$estimate,
             d_band5to15_vs_AHI15 = d2$estimate, d_band_CI_low = d2$conf.int[1], d_band_CI_high = d2$conf.int[2],
             p_band = wilcox.test(c5_15[[v]], c15[[v]])$p.value)
}))
print(es, digits = 3); save_csv(es, "TableS12_AHI5_vs_AHI15_effect_sizes.csv")

# (3) DeLong test: C2 vs C3 as predictors of OSA (paired; same reference set)
valid <- !is.na(vd$Expert_Label)
actual <- as.integer(vd$Expert_Label[valid] == "OSA")
r2 <- pROC::roc(actual, as.integer(vd$Cluster[valid] == 2), quiet = TRUE)
r3 <- pROC::roc(actual, as.integer(vd$Cluster[valid] == 3), quiet = TRUE)
dl <- pROC::roc.test(r2, r3, method = "delong", paired = TRUE)
cat(sprintf("DeLong C2 vs C3 (OSA): AUC %.3f vs %.3f, p = %.3g\n", pROC::auc(r2), pROC::auc(r3), dl$p.value))
save_csv(data.frame(Comparison = "C2 vs C3 as predictors of OSA (paired DeLong)",
                    AUC_C2 = as.numeric(pROC::auc(r2)), AUC_C3 = as.numeric(pROC::auc(r3)), p = dl$p.value),
         "TableS12_DeLong_C2_vs_C3.csv")
