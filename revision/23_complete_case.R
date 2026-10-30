# ============================================================
# 23_complete_case.R — Complete-case sensitivity analysis (R2-8)
#  Re-runs the final pipeline (final UMAP params, Ward k=3) on
#  participants with no missing clustering variables, and compares
#  with the primary (imputed) solution.
# Run after: source("01_setup.R"); source("08_revision_setup.R")
# ============================================================
cat("\n=== 23: Complete-case analysis ===\n")
cc  <- complete.cases(comisa[, clustering_vars])
Xcc <- scale(comisa[cc, clustering_vars])
cl_cc <- hc_cut(fit_umap(Xcc), K_FINAL)
m_cc  <- match_labels(REF_CL[cc], cl_cc)
com5  <- comisa$Diagnosis == "不眠症" & comisa$AHI >= 5
com15 <- comisa$Diagnosis == "不眠症" & comisa$AHI >= 15
out <- data.frame(
  N_complete = sum(cc),
  ARI_vs_reference = round(ari(REF_CL[cc], cl_cc), 3),
  Agreement_pct = round(100 * mean(m_cc == REF_CL[cc]), 1),
  Cluster_sizes = paste(table(m_cc), collapse = "/"),
  N_COMISA5_cc = sum(com5[cc]),
  Pct_COMISA5_in_C1 = round(100 * mean(m_cc[com5[cc]] == 1), 1),
  N_COMISA15_cc = sum(com15[cc]),
  Pct_COMISA15_in_C1 = round(100 * mean(m_cc[com15[cc]] == 1), 1))
print(out)
save_csv(out, "TableS6_complete_case.csv")
