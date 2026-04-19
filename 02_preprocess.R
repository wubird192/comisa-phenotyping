# ============================================================
# COMISA Phenotyping Study
# 02_preprocess.R — Data loading, age adjustment (Seol 2025),
#                   missing data, and MICE imputation
# ============================================================

cat("\n============================================================\n")
cat("  02: Preprocessing\n")
cat("============================================================\n\n")

# ============================================================
# 1. DATA LOADING
# ============================================================

cat("=== 1. Data loading ===\n")
comisa_raw <- read.csv(DATA_PATH, fileEncoding = "UTF-8-BOM")
cat(sprintf("✓ N=%d, variables=%d\n\n", nrow(comisa_raw), ncol(comisa_raw)))

# ============================================================
# 2. CLUSTERING VARIABLE DEFINITION (15 variables)
# ============================================================

cat("=== 2. Clustering variables ===\n")

# PSG variables subject to age adjustment
PSG_VARS_TO_ADJUST <- c("TST", "SE", "WASO", "SL")

# All 15 clustering variables (will be updated after age adjustment)
clustering_vars <- c(
  # PSG — age-adjusted candidates
  "TST", "SE", "WASO", "SL",
  # PSG — not age-adjusted
  "AHI", "Arousal_Index", "ODI3",
  # Psychological scales
  "ISI", "DBAS", "FIRST", "HAS",
  # Additional clinical
  "JESS", "SDS",
  # Demographics
  "Age", "BMI"
)

cat(sprintf("Clustering variables (n=%d):\n", length(clustering_vars)))
cat(paste(" -", clustering_vars, collapse = "\n"), "\n\n")

# ============================================================
# 3. AGE ADJUSTMENT OF PSG VARIABLES (Seol et al., 2025)
#
#    Criterion: linear regression on Age; if R² > 0.01 AND p < 0.05
#               → replace with residuals + grand mean
# ============================================================

cat("=== 3. Age adjustment of PSG variables (Seol et al., 2025) ===\n")

age_adjust_psg <- function(data,
                           psg_vars    = c("TST", "SE", "WASO", "SL"),
                           min_rsq     = 0.01,
                           max_p       = 0.05) {
  adjusted_data    <- data
  adjustment_log   <- data.frame()

  for (var in psg_vars) {
    if (!var %in% names(data)) {
      cat(sprintf("  ⚠ %s not found — skipped\n", var))
      next
    }
    ok <- complete.cases(data[, c(var, "Age")])
    if (mean(ok) < 0.8) {
      cat(sprintf("  ⚠ %s: >20%% missing — skipped\n", var))
      next
    }
    tryCatch({
      fit   <- lm(reformulate("Age", var), data = data[ok, ])
      s     <- summary(fit)
      r2    <- s$r.squared
      pval  <- s$coefficients[2, 4]
      coef  <- s$coefficients[2, 1]

      if (r2 > min_rsq && pval < max_p) {
        resid_full        <- rep(NA_real_, nrow(data))
        resid_full[ok]    <- data[[var]][ok] - predict(fit)
        new_var           <- paste0(var, "_AgeAdj")
        adjusted_data[[new_var]] <- resid_full + mean(data[[var]], na.rm = TRUE)
        row_new <- data.frame(Variable = var, R2 = r2, P = pval,
                              Coef = coef, Adjusted = TRUE,
                              New_Variable = new_var, stringsAsFactors = FALSE)
        cat(sprintf("  ✓ %s: R²=%.3f, p=%.3e → %s\n", var, r2, pval, new_var))
      } else {
        row_new <- data.frame(Variable = var, R2 = r2, P = pval,
                              Coef = coef, Adjusted = FALSE,
                              New_Variable = var, stringsAsFactors = FALSE)
        cat(sprintf("  - %s: R²=%.3f, p=%.3e → no adjustment needed\n", var, r2, pval))
      }
      adjustment_log <- rbind(adjustment_log, row_new)
    }, error = function(e) {
      cat(sprintf("  ✗ %s: error — %s\n", var, e$message))
    })
  }

  cat(sprintf("\nAdjusted: %d/%d PSG variables\n\n",
              sum(adjustment_log$Adjusted), nrow(adjustment_log)))
  list(data = adjusted_data, log = adjustment_log)
}

adj_result     <- age_adjust_psg(comisa_raw, PSG_VARS_TO_ADJUST)
comisa         <- adj_result$data
adjustment_log <- adj_result$log

# Update clustering variable names to age-adjusted versions
for (i in seq_len(nrow(adjustment_log))) {
  if (adjustment_log$Adjusted[i]) {
    clustering_vars[clustering_vars == adjustment_log$Variable[i]] <-
      adjustment_log$New_Variable[i]
  }
}

cat("Updated clustering variables:\n")
cat(paste(" -", clustering_vars, collapse = "\n"), "\n\n")

write.csv(adjustment_log,
          file.path(OUTPUT_DIR, "Age_Adjustment_Summary.csv"),
          row.names = FALSE)
cat("✓ Age_Adjustment_Summary.csv saved\n\n")

# ============================================================
# 4. MISSING DATA SUMMARY
# ============================================================

cat("=== 4. Missing data ===\n")

missing_df <- data.frame(
  Variable        = clustering_vars,
  N_Missing       = sapply(clustering_vars, function(v) sum(is.na(comisa[[v]]))),
  Percent_Missing = sapply(clustering_vars, function(v)
    round(100 * mean(is.na(comisa[[v]])), 1))
)
missing_df <- missing_df[order(-missing_df$N_Missing), ]
rownames(missing_df) <- NULL
print(missing_df)

cat(sprintf(
  "\nComplete cases: %d/%d (%.1f%%)\n\n",
  sum(complete.cases(comisa[, clustering_vars])),
  nrow(comisa),
  100 * mean(complete.cases(comisa[, clustering_vars]))
))

write.csv(missing_df,
          file.path(OUTPUT_DIR, "Missing_Data_Summary.csv"),
          row.names = FALSE)
cat("✓ Missing_Data_Summary.csv saved\n\n")

# ============================================================
# 5. MICE MULTIPLE IMPUTATION (PMM, m=5, maxit=10, seed=42)
# ============================================================

cat("=== 5. MICE multiple imputation ===\n")
cat(sprintf("  m=%d, maxit=%d, method=PMM, seed=%d\n\n",
            MICE_M, MICE_MAXIT, SEED_MICE))

set.seed(SEED_MICE)
imputed <- mice(
  comisa[, clustering_vars, drop = FALSE],
  m        = MICE_M,
  method   = "pmm",
  maxit    = MICE_MAXIT,
  seed     = SEED_MICE,
  printFlag = TRUE
)

cat("\n✓ MICE imputation complete\n\n")

cat("============================================================\n")
cat("  Preprocessing complete. Proceed to 03_clustering.R\n")
cat("============================================================\n\n")
