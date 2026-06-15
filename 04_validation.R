# ============================================================
# COMISA Phenotyping Study
# 04_validation.R — Criterion validity: cluster assignment vs
#                   expert clinical diagnosis (ROC/AUC analysis)
#
# Figure 4 is generated here.
# Diagonal line fix: uses plot.new() + lines() to avoid the
# double diagonal caused by pROC::plot.roc() internal abline().
# ============================================================

cat("\n============================================================\n")
cat("  04: Expert Diagnosis Validation\n")
cat("============================================================\n\n")

# ============================================================
# 1. BUILD EXPERT DIAGNOSTIC LABELS
#
#    Insomnia      : Diagnosis == "不眠症"  (all)
#    OSA           : AHI ≥ 15             (all)
#    COMISA        : Diagnosis == "不眠症" AND AHI ≥ 15  (overwrites above)
#    Insomnia_Only : Diagnosis == "不眠症" AND (AHI < 15 or NA)
# ============================================================

cat("=== 1. Expert diagnostic labels ===\n\n")

vd <- comisa
vd$Expert_Label  <- NA_character_
vd$Insomnia_Only <- NA_character_

vd$Expert_Label[vd$Diagnosis == "\u4e0d\u7720\u75c7"]                                       <- "Insomnia"
vd$Expert_Label[vd$AHI >= 15]                                               <- "OSA"
vd$Expert_Label[vd$Diagnosis == "\u4e0d\u7720\u75c7" & vd$AHI >= 15]                        <- "COMISA"
vd$Insomnia_Only[vd$Diagnosis == "\u4e0d\u7720\u75c7" &
                 (is.na(vd$AHI) | vd$AHI < 15)]                            <- "Insomnia_Only"

cat("Expert label distribution:\n")
print(table(vd$Expert_Label, useNA = "always"))
cat(sprintf("\n  Insomnia (all):           n = %d\n",
            sum(vd$Expert_Label == "Insomnia",  na.rm = TRUE)))
cat(sprintf("    of which Insomnia-only: n = %d\n",
            sum(!is.na(vd$Insomnia_Only))))
cat(sprintf("    of which COMISA:        n = %d\n",
            sum(vd$Expert_Label == "COMISA",    na.rm = TRUE)))
cat(sprintf("  OSA (AHI ≥ 15):           n = %d\n\n",
            sum(vd$Expert_Label == "OSA",       na.rm = TRUE)))

# Cross-tabulation
cat("Cluster × Expert Label cross-tabulation:\n")
print(table(Cluster = vd$Cluster, Expert = vd$Expert_Label, useNA = "always"))

# Chi-square test
valid <- !is.na(vd$Expert_Label) & !is.na(vd$Cluster)
chi   <- chisq.test(vd$Cluster[valid], vd$Expert_Label[valid])
cat(sprintf("\nChi-square test: χ² = %.2f, df = %d, p = %.3e\n\n",
            chi$statistic, chi$parameter, chi$p.value))

# ============================================================
# 2. ROC ANALYSIS FUNCTION
# ============================================================

perform_roc <- function(data, cluster_val, expert_lbl, label) {

  valid  <- !is.na(data$Cluster)
  if (expert_lbl == "Insomnia_Only") {
    actual <- as.integer(!is.na(data$Insomnia_Only[valid]))
  } else {
    valid  <- valid & !is.na(data$Expert_Label)
    actual <- as.integer(data$Expert_Label[valid] == expert_lbl)
  }
  pred <- as.integer(data$Cluster[valid] == cluster_val)

  n_pos <- sum(actual == 1)
  n_neg <- sum(actual == 0)
  cat(sprintf("\n--- %s ---\n  Positive: n=%d  Negative: n=%d\n",
              label, n_pos, n_neg))

  if (n_pos < 5 || n_neg < 5) {
    cat("  ⚠ Insufficient sample — skipped\n")
    return(NULL)
  }

  robj  <- roc(actual, pred, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_v <- auc(robj)
  ci_v  <- ci.auc(robj)

  cm    <- table(Predicted = pred, Actual = actual)
  cat("  Confusion matrix:\n"); print(cm)

  sens <- spec <- ppv <- npv <- NA_real_
  if (nrow(cm) == 2 && ncol(cm) == 2) {
    sens <- cm[2, 2] / sum(cm[, 2])
    spec <- cm[1, 1] / sum(cm[, 1])
    ppv  <- cm[2, 2] / sum(cm[2, ])
    npv  <- cm[1, 1] / sum(cm[1, ])
    cat(sprintf("  Sensitivity: %.3f  Specificity: %.3f\n", sens, spec))
    cat(sprintf("  PPV: %.3f  NPV: %.3f\n", ppv, npv))
  }
  cat(sprintf("  AUC: %.3f (95%% CI: %.3f–%.3f)\n", auc_v, ci_v[1], ci_v[3]))

  list(roc_obj     = robj,
       auc         = as.numeric(auc_v),
       ci          = as.numeric(ci_v),
       confusion   = cm,
       sensitivity = sens,
       specificity = spec,
       ppv         = ppv,
       npv         = npv,
       label       = label)
}

# ============================================================
# 3. RUN ROC ANALYSES (3-cluster solution)
# ============================================================

cat("=== 2. ROC analyses ===\n")

n_cl      <- length(unique(vd$Cluster[!is.na(vd$Cluster)]))
roc_results <- list()

if (n_cl == 3) {

  # Cluster 1 (insomnia-predominant) — two comparisons
  roc_results$Cluster1_InsomniaOnly <- perform_roc(
    vd, 1, "Insomnia_Only", "Cluster 1 vs Insomnia only (AHI < 15)")
  roc_results$Cluster1_COMISA       <- perform_roc(
    vd, 1, "COMISA",        "Cluster 1 vs COMISA (AHI ≥ 15)")

  # Cluster 2 (moderate OSA)
  roc_results$Cluster2_OSA          <- perform_roc(
    vd, 2, "OSA",           "Cluster 2 vs OSA (AHI ≥ 15)")

  # Cluster 3 (severe OSA)
  roc_results$Cluster3_OSA          <- perform_roc(
    vd, 3, "OSA",           "Cluster 3 vs OSA (AHI ≥ 15)")
}

# ============================================================
# 4. FIGURE 4 — ROC curves (single diagonal reference line)
#
#    Implementation note:
#    pROC::plot.roc() calls abline(0,1) internally. To guarantee
#    exactly one reference diagonal, we use plot.new() +
#    plot.window() and draw all curves via lines(). This
#    completely bypasses pROC's rendering path.
# ============================================================

cat("\n=== 3. Figure 4 — ROC curves ===\n")

if (length(roc_results) >= 4) {

  r1a   <- roc_results$Cluster1_InsomniaOnly$roc_obj
  r1b   <- roc_results$Cluster1_COMISA$roc_obj
  r2    <- roc_results$Cluster2_OSA$roc_obj
  r3    <- roc_results$Cluster3_OSA$roc_obj
  auc1a <- roc_results$Cluster1_InsomniaOnly$auc
  auc1b <- roc_results$Cluster1_COMISA$auc
  auc2  <- roc_results$Cluster2_OSA$auc
  auc3  <- roc_results$Cluster3_OSA$auc

  draw_fig4 <- function(cex_leg = 0.70) {
    # Empty canvas — no pROC rendering
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
    axis(1); axis(2); box()
    title(main = "ROC Curves: Cluster vs Expert Diagnosis",
          xlab = "1 \u2212 Specificity", ylab = "Sensitivity")

    # Single diagonal reference line
    segments(0, 0, 1, 1, lty = 2, col = "gray50", lwd = 1.5)

    # ROC curves via lines() — specificities are descending, so rev()
    lines(1 - rev(r1a$specificities), rev(r1a$sensitivities),
          col = CLUSTER_COLORS["1"], lwd = 2.5, lty = 1)
    lines(1 - rev(r1b$specificities), rev(r1b$sensitivities),
          col = CLUSTER_COLORS["1"], lwd = 2.5, lty = 2)
    lines(1 - rev(r2$specificities),  rev(r2$sensitivities),
          col = CLUSTER_COLORS["2"], lwd = 2.5, lty = 1)
    lines(1 - rev(r3$specificities),  rev(r3$sensitivities),
          col = CLUSTER_COLORS["3"], lwd = 2.5, lty = 1)

    legend("bottomright",
           legend = c(sprintf("C1 vs Insomnia only (AUC = %.3f)", auc1a),
                      sprintf("C1 vs COMISA (AUC = %.3f)",         auc1b),
                      sprintf("C2 vs OSA (AUC = %.3f)",             auc2),
                      sprintf("C3 vs OSA (AUC = %.3f)",             auc3)),
           col  = c(CLUSTER_COLORS["1"], CLUSTER_COLORS["1"],
                    CLUSTER_COLORS["2"], CLUSTER_COLORS["3"]),
           lty  = c(1, 2, 1, 1), lwd = 2.5, cex = cex_leg, bty = "n")
  }

  # PNG (preview)
  png(file.path(OUTPUT_DIR, "Figure4_ROC_Expert_Diagnosis.png"),
      width = 180, height = 140, units = "mm", res = 300)
  par(pty = "s", mar = c(5, 4, 4, 2))
  draw_fig4(0.70)
  dev.off()

  # TIFF (journal submission, 300 DPI)
  # macOS quartz device does not support LZW compression; use cairo on
  # Linux/Windows for LZW. On macOS the TIFF is saved uncompressed.
  tiff_path <- file.path(OUTPUT_DIR, "Figure4_ROC_Expert_Diagnosis.tiff")
  on_mac    <- Sys.info()[["sysname"]] == "Darwin"
  if (on_mac) {
    tiff(tiff_path, width = 180, height = 140, units = "mm", res = 300)
  } else {
    tiff(tiff_path, width = 180, height = 140, units = "mm", res = 300,
         compression = "lzw", type = "cairo")
  }
  par(pty = "s", mar = c(5, 4, 4, 2), family = "Arial", ps = 9)
  draw_fig4(0.65)
  dev.off()
  if (on_mac)
    cat("  Note: macOS — TIFF saved without LZW. Convert for submission if needed.\n")
  cat("\u2713 Figure4_ROC_Expert_Diagnosis.png / .tiff saved\n\n")
}

# ============================================================
# 5. SAVE ROC SUMMARY
# ============================================================

if (length(roc_results) >= 4) {
  roc_df <- data.frame(
    Comparison  = c("C1 vs Insomnia only", "C1 vs COMISA",
                    "C2 vs OSA",            "C3 vs OSA"),
    AUC         = sapply(roc_results, function(x) x$auc),
    CI_Lower    = sapply(roc_results, function(x) x$ci[1]),
    CI_Upper    = sapply(roc_results, function(x) x$ci[3]),
    Sensitivity = sapply(roc_results, function(x) x$sensitivity),
    Specificity = sapply(roc_results, function(x) x$specificity),
    PPV         = sapply(roc_results, function(x) x$ppv),
    NPV         = sapply(roc_results, function(x) x$npv)
  )
  write.csv(roc_df,
            file.path(OUTPUT_DIR, "ROC_Results_Summary.csv"),
            row.names = FALSE)
  cat("✓ ROC_Results_Summary.csv saved\n\n")
  print(roc_df)
}

# ============================================================
# 6. COMISA CLUSTER DISTRIBUTION ANALYSIS
#
#    Key finding: most COMISA patients (AHI ≥ 15 + insomnia
#    diagnosis) cluster with the insomnia-predominant phenotype
#    (Cluster 1), supporting the hypothesis that insomnia
#    psychopathology dominates the COMISA clinical presentation.
# ============================================================

cat("\n=== 4. COMISA cluster distribution ===\n\n")

# Primary threshold: AHI ≥ 15
comisa_ahi15 <- vd[!is.na(vd$Expert_Label) &
                   vd$Expert_Label == "COMISA" &
                   !is.na(vd$Cluster), ]
n_comisa_15  <- nrow(comisa_ahi15)
dist_ahi15   <- table(comisa_ahi15$Cluster)

cat(sprintf("COMISA patients (AHI ≥ 15): n = %d\n", n_comisa_15))
cat("Cluster distribution:\n")
print(dist_ahi15)
cat("\nCluster proportions (%):\n")
print(round(100 * prop.table(dist_ahi15), 1))

n_c1_comisa <- as.integer(dist_ahi15["1"])
pct_c1      <- 100 * n_c1_comisa / n_comisa_15

cat(sprintf("\nCOMISA patients in Cluster 1: n=%d (%.1f%%)\n",
            n_c1_comisa, pct_c1))

# Exact binomial 95% CI
binom_ci <- binom.test(n_c1_comisa, n_comisa_15)
cat(sprintf("95%% CI (exact binomial): %.1f%%–%.1f%%\n\n",
            100 * binom_ci$conf.int[1],
            100 * binom_ci$conf.int[2]))

cat("Interpretation:\n")
cat(sprintf("  %.1f%% of COMISA patients (AHI ≥ 15) cluster with the\n", pct_c1))
cat("  insomnia-predominant phenotype (Cluster 1).\n")
cat("  This suggests that insomnia psychopathology—not OSA severity—\n")
cat("  drives the clustering of most COMISA patients.\n\n")

# Sensitivity threshold: AHI ≥ 5
comisa_ahi5 <- vd[vd$Diagnosis == "\u4e0d\u7720\u75c7" &
                  !is.na(vd$AHI) & vd$AHI >= 5 &
                  !is.na(vd$Cluster), ]
n_comisa_5  <- nrow(comisa_ahi5)
n_c1_ahi5   <- sum(comisa_ahi5$Cluster == 1)
pct_c1_ahi5 <- 100 * n_c1_ahi5 / n_comisa_5

cat(sprintf("Sensitivity check (AHI ≥ 5): n=%d total,\n", n_comisa_5))
cat(sprintf("  %d in C1 (%.1f%%)\n\n", n_c1_ahi5, pct_c1_ahi5))

# Save COMISA distribution summary
comisa_dist_df <- rbind(
  data.frame(Threshold = "AHI >= 15", N_COMISA = n_comisa_15,
             N_in_C1 = n_c1_comisa, Pct_in_C1 = round(pct_c1, 1),
             CI_Lower = round(100 * binom_ci$conf.int[1], 1),
             CI_Upper = round(100 * binom_ci$conf.int[2], 1)),
  data.frame(Threshold = "AHI >= 5",  N_COMISA = n_comisa_5,
             N_in_C1 = n_c1_ahi5, Pct_in_C1 = round(pct_c1_ahi5, 1),
             CI_Lower = NA, CI_Upper = NA)
)
write.csv(comisa_dist_df,
          file.path(OUTPUT_DIR, "COMISA_Cluster_Distribution.csv"),
          row.names = FALSE)
cat("✓ COMISA_Cluster_Distribution.csv saved\n\n")

cat("============================================================\n")
cat("  Validation complete. Proceed to 05_stability.R\n")
cat("============================================================\n\n")
