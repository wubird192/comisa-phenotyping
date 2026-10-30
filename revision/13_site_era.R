# ============================================================
# 13_site_era.R — Script E (v2, after clinic reply 2026-09-29)
#  Facts: sites = Gifu (G, n~597) + Iwata (I, n=2) → effectively single-centre;
#  Alice 5 until 2017/4, Alice 6 from 2017/5; AASM 2007→2.0 (2013/7)→2.1 (2014/10)
#  →2.5 (2019/5)→2.6 (2022/8); hypopnea rule (3%/arousal) constant;
#  FIRST introduced 2014/7, DBAS 2010/9; two intake streams (IS / SAS).
#  Only PSG_Year is available → eras are year-based approximations.
#
#  (1) Site: descriptive only (R2-12, m16)
#  (2) Era × cluster: 3 calendar eras, device era, AASM-version era (R2-12)
#  (3) Missingness of FIRST/DBAS by year → MAR by introduction date (R1-M1, R2-8)
#  (4) Intake stream × cluster (new Methods detail)
#  (5) Insomnia duration (IS stream only) by cluster (R1-M2)
#  (6) Sex × cluster (R1-M1)
# ============================================================
cat("\n=== 13: Site / era / missingness / duration ===\n")
stopifnot(SITE_COL %in% names(comisa_raw), YEAR_COL %in% names(comisa_raw))
cl   <- factor(REF_CL)
site <- factor(comisa_raw[[SITE_COL]])
yr   <- as.integer(comisa_raw[[YEAR_COL]])
era3   <- cut(yr, c(2009, 2014, 2019, 2024), labels = c("2010\u20132014", "2015\u20132019", "2020\u20132024"))
device <- factor(ifelse(yr <= 2016, "Alice 5 (\u22642016)", ifelse(yr >= 2018, "Alice 6 (\u22652018)", "2017 (transition)")))
aasm   <- cut(yr, c(2009, 2012, 2014, 2018, 2021, 2024),
              labels = c("AASM 2007 (\u22642012)", "v2.0\u20132.1 (2013\u20132014)", "v2.1 (2015\u20132018)",
                         "v2.5 (2019\u20132021)", "v2.6 (2022\u20132024)"))
cv <- function(tab) { ch <- suppressWarnings(chisq.test(tab)); c(chi2 = unname(ch$statistic), df = unname(ch$parameter),
       p = ch$p.value, V = unname(sqrt(ch$statistic / (sum(tab) * (min(dim(tab)) - 1))))) }
fisher_p <- function(tab) tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 1e4)$p.value, error = function(e) NA)

# (1) Site — descriptive
ts <- table(Site = site, Cluster = cl); print(ts)
save_csv(as.data.frame(ts), "TableS9_cluster_by_site.csv")
cat(sprintf("Iwata n = %d → site-stratified analysis not feasible; report as single-centre limitation.\n", sum(site == "I")))

# (2) Era × cluster
tests <- list()
for (nm in c("era3", "device", "aasm")) {
  f  <- get(nm); tab <- table(Era = f, Cluster = cl)
  cat("\n", nm, ":\n"); print(tab); print(round(100 * prop.table(tab, 1), 1))
  t  <- cv(tab); tests[[nm]] <- data.frame(Stratifier = nm, chi2 = t["chi2"], df = t["df"], p = t["p"],
                                           Cramers_V = t["V"], Fisher_p = fisher_p(tab))
  save_csv(cbind(Stratifier = nm, as.data.frame(tab)), sprintf("TableS9_cluster_by_%s.csv", nm))
}
tests <- do.call(rbind, tests); rownames(tests) <- NULL; print(tests, digits = 3)
save_csv(tests, "TableS9_era_tests.csv")
save_csv(as.data.frame(table(Year = yr)), "TableS9_n_by_year.csv")

# Era-stratified re-clustering (calendar eras only, n permitting): does k=3 structure replicate within era?
strat <- do.call(rbind, lapply(levels(era3), function(e) {
  idx <- which(era3 == e); if (length(idx) < 80) return(data.frame(Era = e, n = length(idx), ARI = NA))
  data.frame(Era = e, n = length(idx), ARI = ari(REF_CL[idx], hc_cut(fit_umap(REF_SCALED[idx, ]), K_FINAL)))
}))
print(strat); save_csv(strat, "TableS9_era_stratified_ARI.csv")

# Leave-one-era-out projection
loeo <- do.call(rbind, lapply(levels(era3), function(e) {
  te <- which(era3 == e); tr <- setdiff(1:N, te)
  m  <- fit_umap(REF_SCALED[tr, ], ret_model = TRUE); clT <- hc_cut(m$embedding, K_FINAL)
  cent <- sapply(1:K_FINAL, function(k) colMeans(m$embedding[clT == k, , drop = FALSE]))
  pr <- apply(uwot::umap_transform(REF_SCALED[te, ], m, verbose = FALSE), 1, function(p) which.min(colSums((cent - p)^2)))
  data.frame(Held_out_era = e, n = length(te), ARI_projected = ari(REF_CL[te], pr))
}))
print(loeo); save_csv(loeo, "TableS9_leave_one_era_out.csv")

# (3) Missingness by year (introduction dates: DBAS 2010/9, FIRST 2014/7)
miss <- do.call(rbind, lapply(c("FIRST", "DBAS", "HAS", "SDS", "ISI", "JESS"), function(v)
  data.frame(Variable = v, Year = sort(unique(yr)),
             N = as.vector(table(yr)), Pct_missing = round(100 * tapply(is.na(comisa_raw[[v]]), yr, mean), 1))))
save_csv(miss, "TableS9_missingness_by_year.csv")
first_pre <- mean(is.na(comisa_raw$FIRST[yr <= 2014])); first_post <- mean(is.na(comisa_raw$FIRST[yr >= 2015]))
cat(sprintf("FIRST missing: %.1f%% in \u22642014 vs %.1f%% in \u22652015 (introduced 2014/7)\n", 100 * first_pre, 100 * first_post))
# Are pre-introduction patients different on observed variables? (MAR vs MNAR argument)
miss_flag <- is.na(comisa_raw$FIRST)
mar <- do.call(rbind, lapply(c("ISI", "DBAS", "AHI", "Age", "BMI", "SDS"), function(v) {
  d <- effsize::cohen.d(comisa_raw[[v]][miss_flag], comisa_raw[[v]][!miss_flag], na.rm = TRUE)
  data.frame(Variable = v, d_missingFIRST_vs_observed = d$estimate, CI_low = d$conf.int[1], CI_high = d$conf.int[2]) }))
print(mar, digits = 3); save_csv(mar, "TableS9_FIRST_missing_vs_observed_d.csv")
tm <- table(FIRST_missing = miss_flag, Cluster = cl); print(tm); print(cv(tm))

# (4) Intake stream × cluster
if (STREAM_COL %in% names(comisa_raw)) {
  st <- table(Stream = comisa_raw[[STREAM_COL]], Cluster = cl); print(st); print(round(100 * prop.table(st, 1), 1)); print(cv(st))
  save_csv(as.data.frame(st), "TableS9_cluster_by_stream.csv")
} else cat("Stream column absent — add 'Stream' (IS/SAS) to the CSV.\n")

# (5) Insomnia duration (IS stream only) — R1-M2 exploratory
if (DUR_COL %in% names(comisa_raw)) {
  dur <- suppressWarnings(as.numeric(comisa_raw[[DUR_COL]]))
  cert <- if (DUR_CERT_COL %in% names(comisa_raw)) comisa_raw[[DUR_CERT_COL]] == "\u660e\u78ba" else NA
  ok <- !is.na(dur)
  cat(sprintf("Duration available: n = %d (certain: %d)\n", sum(ok), sum(cert & ok, na.rm = TRUE)))
  dd <- data.frame(Cluster = cl[ok], Months = dur[ok], Certain = cert[ok])
  summ <- do.call(rbind, lapply(split(dd, dd$Cluster), function(x) data.frame(Cluster = x$Cluster[1], n = nrow(x),
    Median = median(x$Months), Q1 = quantile(x$Months, .25), Q3 = quantile(x$Months, .75), n_certain = sum(x$Certain, na.rm = TRUE))))
  rownames(summ) <- NULL; print(summ)
  cat(sprintf("Kruskal p = %.3g (all); %.3g (certain only)\n", kruskal.test(Months ~ Cluster, dd)$p.value,
              tryCatch(kruskal.test(Months ~ Cluster, dd[dd$Certain %in% TRUE, ])$p.value, error = function(e) NA)))
  # Within-cluster: does duration relate to psychological load? (Spearman)
  rho <- sapply(c("ISI", "DBAS", "FIRST", "HAS"), function(v) cor(dur[ok], comisa_raw[[v]][ok], method = "spearman", use = "complete.obs"))
  print(round(rho, 3))
  save_csv(summ, "TableS9_insomnia_duration_by_cluster.csv")
  save_csv(data.frame(Variable = names(rho), Spearman_rho_with_duration = rho), "TableS9_duration_correlations.csv")
} else cat("Duration column absent (IS stream only) — skipped.\n")

# (6) Sex × cluster
if (SEX_COL %in% names(comisa_raw)) {
  sx <- table(Cluster = cl, Male = comisa_raw[[SEX_COL]] == 1); print(sx); print(cv(sx))
  save_csv(as.data.frame(sx), "TableS9_sex_by_cluster.csv")
}
