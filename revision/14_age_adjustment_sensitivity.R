# ============================================================
# 14_age_adjustment_sensitivity.R — Script F
#  Age-adjustment strategy sensitivity (R2-14)
#   (i)   none            — raw PSG variables
#   (ii)  uniform         — residualize ALL 7 PSG variables on age
#   (iii) threshold-based — current pipeline (R2>0.01 & p<0.05)
#  Each: MICE (same seed) → take the same imputation index as the main
#  analysis → scale → UMAP (final params) → HC k=3 → ARI vs reference
# Output: S10 Table
# ============================================================
cat("\n=== 14: Age-adjustment sensitivity ===\n")
PSG_ALL <- c("TST", "SE", "WASO", "SL", "AHI", "Arousal_Index", "ODI3")
BASE_VARS <- c("TST", "SE", "WASO", "SL", "AHI", "Arousal_Index", "ODI3",
               "ISI", "DBAS", "FIRST", "HAS", "JESS", "SDS", "Age", "BMI")

residualize <- function(data, vars) {
  for (v in vars) {
    ok <- complete.cases(data[, c(v, "Age")])
    fit <- lm(reformulate("Age", v), data = data[ok, ])
    r <- rep(NA_real_, nrow(data)); r[ok] <- data[[v]][ok] - predict(fit)
    data[[v]] <- r + mean(data[[v]], na.rm = TRUE)
  }
  data
}
build <- function(variant) {
  d <- comisa_raw[, BASE_VARS]
  if (variant == "uniform") d <- residualize(d, PSG_ALL)
  if (variant == "threshold") {
    adj <- adjustment_log$Variable[adjustment_log$Adjusted]
    d <- residualize(d, adj)
  }
  set.seed(SEED_MICE)
  imp <- mice::mice(d, m = MICE_M, method = "pmm", maxit = MICE_MAXIT, seed = SEED_MICE, printFlag = FALSE)
  X <- scale(mice::complete(imp, best_imp))
  hc_cut(fit_umap(X), K_FINAL)
}
variants <- c("none", "uniform", "threshold")
cls <- lapply(variants, build); names(cls) <- variants
out <- data.frame(Variant = variants,
                  ARI_vs_reference = sapply(cls, ari, REF_CL),
                  Cluster_sizes = sapply(cls, function(c) paste(sort(table(c), decreasing = TRUE), collapse = "/")))
comisa_idx <- which(!is.na(vd$Expert_Label) & vd$Expert_Label == "COMISA")
out$Pct_COMISA_in_matched_C1 <- sapply(cls, function(c) {
  m <- match_labels(REF_CL, c); round(100 * mean(m[comisa_idx] == 1), 1) })
print(out, digits = 3); save_csv(out, "TableS10_age_adjustment_sensitivity.csv")
cat("Note: 'threshold' should reproduce the reference (ARI ~1) — serves as a check of the rebuild.\n")
