# ============================================================
# COMISA Phenotyping Study
# 07_tables.R — Generate Table 1 (demographics), Table 2
#               (cluster characteristics by phenotype), and
#               Supplementary Table S2 (AHI threshold sensitivity)
# ============================================================

cat("\n============================================================\n")
cat("  07: Tables\n")
cat("============================================================\n\n")

# ============================================================
# HELPER FUNCTIONS
# ============================================================

fmt_mean_sd <- function(x)
  sprintf("%.1f \u00b1 %.1f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

fmt_median_iqr <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
  sprintf("%.1f (%.1f\u2013%.1f)", q[2], q[1], q[3])
}

fmt_n_pct <- function(x, total)
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE),
          100 * sum(x, na.rm = TRUE) / total)

is_normal <- function(x, threshold = 0.05) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(FALSE)
  tryCatch(shapiro.test(x)$p.value > threshold, error = function(e) FALSE)
}

anova_or_kruskal <- function(data, var, cl_var = "Cluster") {
  df  <- data.frame(y = data[[var]], cl = factor(data[[cl_var]]))
  df  <- df[complete.cases(df), ]
  if (is_normal(df$y)) {
    p <- summary(aov(y ~ cl, data = df))[[1]]$`Pr(>F)`[1]
  } else {
    p <- kruskal.test(y ~ cl, data = df)$p.value
  }
  p
}

fmt_p <- function(p) {
  if (is.na(p)) return("—")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# ============================================================
# TABLE 1 — Overall sample characteristics
# ============================================================

cat("=== Table 1: Overall sample characteristics ===\n\n")

N <- nrow(comisa)

# Categorical: sex
male_n <- sum(comisa$Sex == 1, na.rm = TRUE)  # adjust column name if needed

# Continuous variables
t1_vars <- list(
  list(var = "Age",          label = "Age (years)"),
  list(var = "BMI",          label = "BMI (kg/m²)"),
  list(var = "AHI",          label = "AHI (events/h)"),
  list(var = "ISI",          label = "ISI score"),
  list(var = "DBAS",         label = "DBAS score"),
  list(var = "FIRST",        label = "FIRST score"),
  list(var = "HAS",          label = "HAS score"),
  list(var = "SDS",          label = "SDS score"),
  list(var = "JESS",         label = "ESS score"),
  list(var = "TST",          label = "TST (min)"),
  list(var = "SE",           label = "Sleep efficiency (proportion)"),
  list(var = "WASO",         label = "WASO (min)"),
  list(var = "Arousal_Index",label = "Arousal index (events/h)"),
  list(var = "ODI3",         label = "ODI3 (events/h)")
)
t1_vars <- t1_vars[sapply(t1_vars, function(x) x$var %in% names(comisa))]

table1_rows <- lapply(t1_vars, function(x) {
  v   <- comisa[[x$var]]
  nrm <- is_normal(v)
  data.frame(
    Variable    = x$label,
    N_nonmissing = sum(!is.na(v)),
    Value       = if (nrm) fmt_mean_sd(v) else fmt_median_iqr(v),
    Distribution = if (nrm) "mean ± SD" else "median (IQR)",
    stringsAsFactors = FALSE
  )
})
table1 <- do.call(rbind, table1_rows)

cat(sprintf("N = %d\n", N))
print(table1)

write.csv(table1, file.path(OUTPUT_DIR, "Table1_Sample_Characteristics.csv"),
          row.names = FALSE)
cat("✓ Table1_Sample_Characteristics.csv saved\n\n")

# ============================================================
# TABLE 2 — Cluster characteristics
# ============================================================

cat("=== Table 2: Cluster characteristics ===\n\n")

t2_vars <- list(
  list(var = "Age",           label = "Age (years)"),
  list(var = "BMI",           label = "BMI (kg/m²)"),
  list(var = "AHI",           label = "AHI (events/h)"),
  list(var = "ODI3",          label = "ODI3 (events/h)"),
  list(var = "Arousal_Index", label = "Arousal index (events/h)"),
  list(var = "ISI",           label = "ISI score"),
  list(var = "DBAS",          label = "DBAS score"),
  list(var = "FIRST",         label = "FIRST score"),
  list(var = "HAS",           label = "HAS score"),
  list(var = "SDS",           label = "SDS score"),
  list(var = "JESS",          label = "ESS score"),
  list(var = "TST",           label = "TST (min)"),
  list(var = "SE",            label = "Sleep efficiency (proportion)"),
  list(var = "WASO",          label = "WASO (min)")
)
t2_vars <- t2_vars[sapply(t2_vars, function(x) x$var %in% names(comisa))]

cl_sizes <- table(comisa$Cluster)

table2_rows <- lapply(t2_vars, function(x) {
  v   <- comisa[[x$var]]
  cls <- comisa$Cluster
  nrm <- is_normal(v)

  row <- data.frame(
    Variable     = x$label,
    Distribution = if (nrm) "mean ± SD" else "median (IQR)",
    stringsAsFactors = FALSE
  )
  for (k in seq_len(n_clusters)) {
    vk <- v[cls == k & !is.na(cls)]
    row[[paste0("C", k, sprintf("_n%d", cl_sizes[k]))]] <-
      if (nrm) fmt_mean_sd(vk) else fmt_median_iqr(vk)
  }
  row$P_value <- fmt_p(anova_or_kruskal(
    data.frame(setNames(list(v), x$var), Cluster = cls), x$var))
  row
})

table2 <- do.call(rbind, table2_rows)
print(table2)

write.csv(table2, file.path(OUTPUT_DIR, "Table2_Cluster_Characteristics.csv"),
          row.names = FALSE)
cat("✓ Table2_Cluster_Characteristics.csv saved\n\n")

# ============================================================
# SUPPLEMENTARY TABLE S2 — AHI threshold sensitivity analysis
#   Compare COMISA defined by AHI >=15 vs AHI >=5
#   Structure matches Supplementary_Materials_v8 exactly:
#   - mean (SD) for continuous; n (%) for categorical
#   - t-test p-values; Cohen's d effect sizes
#   - Sections: Sample/Cluster, Demographics, Insomnia Severity,
#     Psychological Vulnerability, Respiratory Parameters,
#     Criterion Validity (ROC, N=599)
# ============================================================

cat("=== Supplementary Table S2: AHI sensitivity analysis ===\n\n")

library(effsize)

# --- Define COMISA groups ---
comisa_ahi15 <- comisa$Diagnosis == "\u4e0d\u7720\u75c7" &
                !is.na(comisa$AHI) & comisa$AHI >= 15
comisa_ahi5  <- comisa$Diagnosis == "\u4e0d\u7720\u75c7" &
                !is.na(comisa$AHI) & comisa$AHI >=  5

d15 <- comisa[comisa_ahi15 & !is.na(comisa$Cluster), ]
d5  <- comisa[comisa_ahi5  & !is.na(comisa$Cluster), ]
n15 <- nrow(d15)
n5  <- nrow(d5)

cat(sprintf("COMISA (AHI >=15): n=%d (%.1f%%)\n", n15, 100*n15/nrow(comisa)))
cat(sprintf("COMISA (AHI >=5):  n=%d (%.1f%%)\n\n", n5,  100*n5/nrow(comisa)))

# --- Helper: mean(SD) ---
fmt_mean_sd_s2 <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

# --- Helper: continuous row with t-test and Cohen's d ---
s2_cont_row <- function(label, v15, v5) {
  v15 <- v15[!is.na(v15)]; v5 <- v5[!is.na(v5)]
  p   <- tryCatch(t.test(v15, v5)$p.value, error = function(e) NA)
  d   <- tryCatch(cohen.d(v15, v5)$estimate, error = function(e) NA)
  data.frame(
    Characteristic = label,
    COMISA_AHI15   = fmt_mean_sd_s2(v15),
    COMISA_AHI5    = fmt_mean_sd_s2(v5),
    P_value        = fmt_p(p),
    Cohens_d       = if (is.na(d)) "---" else sprintf("%.2f", d),
    stringsAsFactors = FALSE
  )
}

# --- Helper: section header row ---
s2_header <- function(label) {
  data.frame(Characteristic = label,
             COMISA_AHI15 = "", COMISA_AHI5 = "",
             P_value = "", Cohens_d = "",
             stringsAsFactors = FALSE)
}

# --- Cluster distribution ---
n_c1_15 <- sum(d15$Cluster == 1, na.rm = TRUE)
n_c2_15 <- sum(d15$Cluster == 2, na.rm = TRUE)
n_c3_15 <- sum(d15$Cluster == 3, na.rm = TRUE)
n_c1_5  <- sum(d5$Cluster  == 1, na.rm = TRUE)
n_c2_5  <- sum(d5$Cluster  == 2, na.rm = TRUE)
n_c3_5  <- sum(d5$Cluster  == 3, na.rm = TRUE)

p_c1 <- prop.test(c(n_c1_15, n_c1_5), c(n15, n5))$p.value

# --- ROC with N=599 reference population ---
all_valid     <- comisa[!is.na(comisa$Cluster), ]
lbl15_599     <- as.integer(all_valid$Diagnosis == "\u4e0d\u7720\u75c7" &
                              !is.na(all_valid$AHI) & all_valid$AHI >= 15)
lbl5_599      <- as.integer(all_valid$Diagnosis == "\u4e0d\u7720\u75c7" &
                              !is.na(all_valid$AHI) & all_valid$AHI >= 5)
pred_c1_599   <- as.integer(all_valid$Cluster == 1)

roc15_599 <- roc(lbl15_599, pred_c1_599, levels=c(0,1), direction="<", quiet=TRUE)
roc5_599  <- roc(lbl5_599,  pred_c1_599, levels=c(0,1), direction="<", quiet=TRUE)
auc15_599 <- as.numeric(auc(roc15_599)); ci15_599 <- as.numeric(ci.auc(roc15_599))
auc5_599  <- as.numeric(auc(roc5_599));  ci5_599  <- as.numeric(ci.auc(roc5_599))
delong_p  <- roc.test(roc15_599, roc5_599, method="delong")$p.value

# Sensitivity / specificity / PPV / NPV for each threshold
get_perf <- function(lbl, pred) {
  cm <- table(pred, lbl)
  if (!all(dim(cm) == c(2,2))) return(rep(NA,4))
  c(sens = cm[2,2]/sum(cm[,2]),
    spec = cm[1,1]/sum(cm[,1]),
    ppv  = cm[2,2]/sum(cm[2,]),
    npv  = cm[1,1]/sum(cm[1,]))
}
perf15 <- get_perf(lbl15_599, pred_c1_599)
perf5  <- get_perf(lbl5_599,  pred_c1_599)

# --- Assemble table ---
table_s2 <- rbind(

  # Section 1: Sample size and cluster distribution
  s2_header("Sample Size and Cluster Distribution"),
  data.frame(
    Characteristic = "Total COMISA patients, n (%)",
    COMISA_AHI15   = sprintf("%d (%.1f%%)", n15, 100*n15/nrow(comisa)),
    COMISA_AHI5    = sprintf("%d (%.1f%%)", n5,  100*n5/nrow(comisa)),
    P_value = "< 0.001", Cohens_d = "---", stringsAsFactors = FALSE),
  data.frame(
    Characteristic = "Cluster 1 (Insomnia-predominant), n (%)",
    COMISA_AHI15   = sprintf("%d (%.1f%%)", n_c1_15, 100*n_c1_15/n15),
    COMISA_AHI5    = sprintf("%d (%.1f%%)", n_c1_5,  100*n_c1_5/n5),
    P_value = fmt_p(p_c1), Cohens_d = "---", stringsAsFactors = FALSE),
  data.frame(
    Characteristic = "Cluster 2 (Moderate OSA), n (%)",
    COMISA_AHI15   = sprintf("%d (%.1f%%)", n_c2_15, 100*n_c2_15/n15),
    COMISA_AHI5    = sprintf("%d (%.1f%%)", n_c2_5,  100*n_c2_5/n5),
    P_value = "---", Cohens_d = "---", stringsAsFactors = FALSE),
  data.frame(
    Characteristic = "Cluster 3 (Severe OSA), n (%)",
    COMISA_AHI15   = sprintf("%d (%.1f%%)", n_c3_15, 100*n_c3_15/n15),
    COMISA_AHI5    = sprintf("%d (%.1f%%)", n_c3_5,  100*n_c3_5/n5),
    P_value = "---", Cohens_d = "---", stringsAsFactors = FALSE),

  # Section 2: Demographics
  s2_header("Demographic Characteristics"),
  s2_cont_row("Age, years",            d15$Age, d5$Age),
  s2_cont_row("Body mass index, kg/m2", d15$BMI, d5$BMI),

  # Section 3: Insomnia Symptom Severity
  s2_header("Insomnia Symptom Severity"),
  s2_cont_row("Insomnia Severity Index",   d15$ISI,  d5$ISI),
  s2_cont_row("Sleep latency, min",
    if ("SL" %in% names(d15)) d15$SL else rep(NA, n15),
    if ("SL" %in% names(d5))  d5$SL  else rep(NA, n5)),
  s2_cont_row("Wake after sleep onset, min", d15$WASO, d5$WASO),
  s2_cont_row("Total sleep time, min",       d15$TST,  d5$TST),
  s2_cont_row("Sleep efficiency",
    if ("SE_AgeAdj" %in% names(d15)) d15$SE_AgeAdj else d15$SE,
    if ("SE_AgeAdj" %in% names(d5))  d5$SE_AgeAdj  else d5$SE),

  # Section 4: Psychological Vulnerability
  s2_header("Psychological Vulnerability"),
  s2_cont_row("DBAS-16 score",   d15$DBAS,  d5$DBAS),
  s2_cont_row("Hyperarousal Scale score",   d15$HAS,   d5$HAS),
  s2_cont_row("Ford Insomnia Response to Stress Test", d15$FIRST, d5$FIRST),
  s2_cont_row("Self-rating Depression Scale", d15$SDS, d5$SDS),
  s2_cont_row("Japanese Epworth Sleepiness Scale",
    if ("JESS" %in% names(d15)) d15$JESS else rep(NA, n15),
    if ("JESS" %in% names(d5))  d5$JESS  else rep(NA, n5)),

  # Section 5: Respiratory Parameters
  s2_header("Respiratory Parameters"),
  s2_cont_row("Apnoea-hypopnoea index, events/h", d15$AHI,  d5$AHI),
  s2_cont_row("Oxygen desaturation index (3%), events/h",
    if ("ODI3" %in% names(d15)) d15$ODI3 else rep(NA, n15),
    if ("ODI3" %in% names(d5))  d5$ODI3  else rep(NA, n5)),
  s2_cont_row("Arousal index, events/h",
    if ("Arousal_Index" %in% names(d15)) d15$Arousal_Index else rep(NA, n15),
    if ("Arousal_Index" %in% names(d5))  d5$Arousal_Index  else rep(NA, n5)),

  # Section 6: Criterion Validity (N=599)
  s2_header("Criterion Validity (Cluster 1 vs COMISA, N=599)"),
  data.frame(
    Characteristic = "Area under ROC curve (95% CI)",
    COMISA_AHI15   = sprintf("%.2f (%.2f-%.2f)", auc15_599, ci15_599[1], ci15_599[3]),
    COMISA_AHI5    = sprintf("%.2f (%.2f-%.2f)", auc5_599,  ci5_599[1],  ci5_599[3]),
    P_value = if (delong_p >= 0.05) "> 0.05" else fmt_p(delong_p),
    Cohens_d = "---", stringsAsFactors = FALSE),
  data.frame(
    Characteristic = "Sensitivity, %",
    COMISA_AHI15   = sprintf("%.1f", 100*perf15["sens"]),
    COMISA_AHI5    = sprintf("%.1f", 100*perf5["sens"]),
    P_value = "---", Cohens_d = "---", stringsAsFactors = FALSE),
  data.frame(
    Characteristic = "Specificity, %",
    COMISA_AHI15   = sprintf("%.1f", 100*perf15["spec"]),
    COMISA_AHI5    = sprintf("%.1f", 100*perf5["spec"]),
    P_value = "---", Cohens_d = "---", stringsAsFactors = FALSE),
  data.frame(
    Characteristic = "Positive predictive value, %",
    COMISA_AHI15   = sprintf("%.1f", 100*perf15["ppv"]),
    COMISA_AHI5    = sprintf("%.1f", 100*perf5["ppv"]),
    P_value = "---", Cohens_d = "---", stringsAsFactors = FALSE),
  data.frame(
    Characteristic = "Negative predictive value, %",
    COMISA_AHI15   = sprintf("%.1f", 100*perf15["npv"]),
    COMISA_AHI5    = sprintf("%.1f", 100*perf5["npv"]),
    P_value = "---", Cohens_d = "---", stringsAsFactors = FALSE)
)

print(table_s2)
write.csv(table_s2, file.path(OUTPUT_DIR, "TableS2_AHI_Sensitivity_Analysis.csv"),
          row.names = FALSE)
cat("\u2713 TableS2_AHI_Sensitivity_Analysis.csv saved\n\n")
cat("  Note: Continuous data presented as mean (SD); t-test p-values;\n")
cat("  AUC/ROC calculated using all 599 participants as reference.\n\n")

cat("============================================================\n")
cat("  Tables complete. All analyses done --", OUTPUT_DIR, "\n")
cat("============================================================\n\n")
