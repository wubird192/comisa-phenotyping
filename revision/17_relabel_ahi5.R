# ============================================================
# 17_relabel_ahi5.R — Script J  (Plan A: primary definition = registered AHI>=5)
#  Rebuild expert labels at a given AHI threshold and regenerate
#  Table 2 (cross-tab), Table 3 (contingency metrics + balanced accuracy/AUC),
#  and the COMISA cluster distribution with exact binomial CI.
#  Runs for BOTH thresholds (5 = primary, 15 = sensitivity) so the two
#  can be reported side by side (R2-7).
#  Cluster assignments are unchanged; only reference labels move.
#  2026-10-18: made k-generic (works for any number of clusters; C1 is the
#  low-AHI/high-burden cluster, all other clusters are treated as OSA-type).
# ============================================================
cat("\n=== 17: Relabel at AHI thresholds ===\n")
INS <- comisa$Diagnosis == "不眠症"
K_LEV <- sort(unique(REF_CL)); K <- length(K_LEV)
cat(sprintf("Clusters present: %s (k = %d)\n", paste(K_LEV, collapse = ", "), K))

build_labels <- function(thr) {
  osa <- !is.na(comisa$AHI) & comisa$AHI >= thr
  lab <- rep("Neither", N)
  lab[INS & !osa] <- "Insomnia_only"; lab[!INS & osa] <- "OSA_only"; lab[INS & osa] <- "COMISA"
  factor(lab, levels = c("COMISA", "Insomnia_only", "Neither", "OSA_only"))
}

# Binary-predictor metrics with exact-binomial CIs; AUC = balanced accuracy,
# CI via Hanley–McNeil-equivalent variance used in the original pipeline (matches pROC/DeLong)
bci <- function(x, k) {
  if (is.na(k) || k == 0) return("NA")
  ci <- binom.test(x, k)$conf.int
  sprintf("%.1f (%.1f–%.1f)", 100 * x / k, 100 * ci[1], 100 * ci[2])
}
metrics <- function(pred, actual) {
  tp <- sum(pred & actual); fp <- sum(pred & !actual); fn <- sum(!pred & actual); tn <- sum(!pred & !actual)
  m <- tp + fn; n <- fp + tn
  sens <- tp / m; spec <- tn / n
  auc <- (sens + spec) / 2
  v   <- sens * (1 - sens) / (4 * (m - 1)) + spec * (1 - spec) / (4 * (n - 1))
  data.frame(TP = tp, FP = fp, FN = fn, TN = tn, N = m + n,
             AUC = sprintf("%.2f (%.2f–%.2f)", auc, auc - 1.959964 * sqrt(v), auc + 1.959964 * sqrt(v)),
             Sensitivity = bci(tp, m), Specificity = bci(tn, n),
             PPV = bci(tp, tp + fp), NPV = bci(tn, tn + fn), Balanced_accuracy = round(auc, 3))
}

out_t3 <- list(); out_t2 <- list(); out_dist <- list()
for (thr in c(5, 15)) {
  lab <- build_labels(thr); cl <- REF_CL
  cat(sprintf("\n----- AHI >= %d -----\n", thr))
  t2 <- table(Cluster = cl, Expert = lab); print(t2); print(addmargins(t2))
  chi <- suppressWarnings(chisq.test(t2))
  cat(sprintf("Chi2 = %.1f, df = %d, p = %.3g\n", chi$statistic, chi$parameter, chi$p.value))
  t2df <- as.data.frame.matrix(t2); t2df <- cbind(Threshold = thr, Cluster = rownames(t2df), t2df)
  out_t2[[as.character(thr)]] <- t2df

  # Table 3 — definitive-diagnosis reference set (exclude Neither), as in Methods
  def <- lab != "Neither"
  osa_any <- lab %in% c("OSA_only", "COMISA")
  ins_any <- lab %in% c("Insomnia_only", "COMISA")
  rows <- rbind(
    cbind(Comparison = "Cluster 1 vs COMISA",                 metrics((cl == 1)[def], (lab == "COMISA")[def])),
    cbind(Comparison = "Cluster 1 vs insomnia disorder only", metrics((cl == 1)[def], (lab == "Insomnia_only")[def])),
    cbind(Comparison = "Cluster 1 vs any insomnia disorder",  metrics((cl == 1)[def], ins_any[def])))
  for (k in setdiff(K_LEV, 1))
    rows <- rbind(rows, cbind(Comparison = sprintf("Cluster %s vs OSA (any)", k), metrics((cl == k)[def], osa_any[def])))
  if (K > 2)
    rows <- rbind(rows, cbind(Comparison = sprintf("Cluster %s vs OSA (any)", paste(setdiff(K_LEV, 1), collapse = " or ")),
                              metrics((cl != 1)[def], osa_any[def])))
  # Full-cohort versions (n = 599) for the supplement
  rows_full <- rbind(
    cbind(Comparison = "Cluster 1 vs COMISA",                 metrics(cl == 1, lab == "COMISA")),
    cbind(Comparison = "Cluster 1 vs insomnia disorder only", metrics(cl == 1, lab == "Insomnia_only")))
  for (k in setdiff(K_LEV, 1))
    rows_full <- rbind(rows_full, cbind(Comparison = sprintf("Cluster %s vs OSA (any)", k), metrics(cl == k, osa_any)))
  rows$Threshold <- thr; rows$Reference_set <- "definitive diagnosis"
  rows_full$Threshold <- thr; rows_full$Reference_set <- "full cohort"
  out_t3[[as.character(thr)]] <- rbind(rows, rows_full)
  print(rows[, c("Comparison", "N", "AUC", "Sensitivity", "Specificity", "PPV", "NPV")])

  # COMISA distribution across clusters (columns C1..CK generated dynamically)
  n_com <- sum(lab == "COMISA")
  ci <- binom.test(sum(cl[lab == "COMISA"] == 1), n_com)
  d <- data.frame(Threshold = thr, N_COMISA = n_com)
  for (k in K_LEV) d[[paste0("C", k)]] <- sum(cl[lab == "COMISA"] == k)
  d$Pct_C1 <- round(100 * ci$estimate, 1); d$CI_low <- round(100 * ci$conf.int[1], 1); d$CI_high <- round(100 * ci$conf.int[2], 1)
  d$N_OSA_only <- sum(lab == "OSA_only"); d$OSA_only_in_C1 <- sum(cl[lab == "OSA_only"] == 1)
  d$Pct_OSA_only_in_C1 <- round(100 * mean(cl[lab == "OSA_only"] == 1), 1)
  out_dist[[as.character(thr)]] <- d
  print(d)
}
save_csv(do.call(rbind, out_t2),   "Table2_crosstab_by_threshold.csv")
save_csv(do.call(rbind, out_t3),   "Table3_metrics_by_threshold.csv")
save_csv(do.call(rbind, out_dist), "COMISA_distribution_by_threshold.csv")
