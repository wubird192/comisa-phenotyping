# ============================================================
# 18_synthetic_data.R — Public synthetic dataset for pipeline verification
#  (Journal Requirement 4; Reviewer 1 minor 4; Reviewer 2 data-availability "No")
#
#  Generates a synthetic cohort (n = 599) with the same variables, marginal
#  distributions, missingness pattern and multivariate structure as the real
#  data, using sequential CART synthesis (synthpop). No real record is copied.
#  Then re-runs the clustering pipeline on the synthetic data to confirm it
#  reproduces a 3-cluster solution of similar shape (not identical values).
#
#  Output: data/synthetic_comisa.csv (for GitHub), output/revision/synthetic_*
# ============================================================
cat("\n=== 18: Synthetic dataset ===\n")
if (!requireNamespace("synthpop", quietly = TRUE)) install.packages("synthpop")
library(synthpop)

# ---- 1. Variables to release ----
# Clinical variables + labels + the revision covariates (site, year, stream, duration).
# Exclude free text / IDs. Duration kept only as a coarse variable.
REL <- c("Diagnosis", SEX_COL, "Age", "BMI", "AHI", "ODI3", "Arousal_Index", "TST", "SE",
         "WASO", "SL", "ISI", "DBAS", "FIRST", "HAS", "JESS", "SDS",
         SITE_COL, YEAR_COL, STREAM_COL, DUR_COL, DUR_CERT_COL)
REL <- REL[REL %in% names(comisa_raw)]
real <- comisa_raw[, REL]
for (v in c("Diagnosis", SEX_COL, SITE_COL, STREAM_COL, DUR_CERT_COL)) if (v %in% REL) real[[v]] <- factor(real[[v]])
if (DUR_COL %in% REL) real[[DUR_COL]] <- suppressWarnings(as.numeric(real[[DUR_COL]]))

# ---- 2. Synthesis ----
# Visit sequence: structural variables first so clinical values are conditioned on them.
seq_order <- intersect(c(YEAR_COL, STREAM_COL, SITE_COL, "Diagnosis", SEX_COL, "Age", "BMI",
                         "AHI", "ODI3", "Arousal_Index", "TST", "SE", "WASO", "SL",
                         "ISI", "DBAS", "FIRST", "HAS", "JESS", "SDS", DUR_CERT_COL, DUR_COL), REL)
set.seed(SEED_MICE + 99)
syn_obj <- syn(real, method = "cart", visit.sequence = seq_order, m = 1,
               minnumlevels = 5, print.flag = FALSE)
syn_dat <- syn_obj$syn

# Enforce logical constraints the CART may violate
syn_dat$SE   <- pmin(pmax(syn_dat$SE, 0), 100)
syn_dat$AHI  <- pmax(syn_dat$AHI, 0); syn_dat$ODI3 <- pmax(syn_dat$ODI3, 0)
syn_dat$WASO <- pmax(syn_dat$WASO, 0); syn_dat$SL <- pmax(syn_dat$SL, 0)
syn_dat$ISI  <- round(pmin(pmax(syn_dat$ISI, 0), 28))
syn_dat$DBAS <- pmin(pmax(syn_dat$DBAS, 0), 160)   # DBAS-16 total, 0–160
syn_dat$HAS  <- round(pmin(pmax(syn_dat$HAS, 0), 72)); syn_dat$JESS <- round(pmin(pmax(syn_dat$JESS, 0), 24))
syn_dat$SDS  <- round(pmin(pmax(syn_dat$SDS, 20), 80)); syn_dat$FIRST <- round(pmin(pmax(syn_dat$FIRST, 9), 36))
# FIRST / DBAS must be missing before their introduction (keeps the MAR story consistent)
if (YEAR_COL %in% REL) {
  syn_dat$FIRST[syn_dat[[YEAR_COL]] <= 2014] <- NA
  syn_dat$DBAS [syn_dat[[YEAR_COL]] <= 2010] <- NA
}
syn_dat <- cbind(ID = sprintf("SYN-%03d", seq_len(nrow(syn_dat))), syn_dat)

# ---- 3. Disclosure check: no synthetic record identical to a real one on the 15 clinical vars ----
cv15 <- c("Age", "BMI", "AHI", "ODI3", "Arousal_Index", "TST", "SE", "WASO", "SL", "ISI", "DBAS", "FIRST", "HAS", "JESS", "SDS")
key_real <- do.call(paste, c(round(real[, cv15], 1), sep = "|"))
key_syn  <- do.call(paste, c(round(syn_dat[, cv15], 1), sep = "|"))
cat(sprintf("Synthetic rows identical to a real row on all 15 clinical variables: %d (must be 0)\n", sum(key_syn %in% key_real)))
rep_tab <- replicated.uniques(syn_obj, real)
cat(sprintf("synthpop replicated uniques: %d\n", rep_tab$no.replications))

# ---- 4. Utility: marginals and correlations ----
comp <- compare(syn_obj, real, vars = cv15, print.flag = FALSE, plot = FALSE)
util <- data.frame(Variable = cv15,
  Real_mean = sapply(cv15, function(v) mean(real[[v]], na.rm = TRUE)),
  Syn_mean  = sapply(cv15, function(v) mean(syn_dat[[v]], na.rm = TRUE)),
  Real_SD   = sapply(cv15, function(v) sd(real[[v]], na.rm = TRUE)),
  Syn_SD    = sapply(cv15, function(v) sd(syn_dat[[v]], na.rm = TRUE)),
  Real_pct_missing = sapply(cv15, function(v) 100 * mean(is.na(real[[v]]))),
  Syn_pct_missing  = sapply(cv15, function(v) 100 * mean(is.na(syn_dat[[v]]))))
print(util, digits = 3); save_csv(util, "synthetic_utility_marginals.csv")
cor_diff <- cor(real[, cv15], use = "pairwise", method = "spearman") - cor(syn_dat[, cv15], use = "pairwise", method = "spearman")
cat(sprintf("Max |Δ Spearman r| real vs synthetic: %.3f; mean |Δ|: %.3f\n", max(abs(cor_diff)), mean(abs(cor_diff[upper.tri(cor_diff)]))))
cat("Diagnosis:\n"); print(rbind(Real = table(real$Diagnosis), Syn = table(syn_dat$Diagnosis)))

# ---- 5. Does the pipeline reproduce a 3-cluster structure on synthetic data? ----
set.seed(SEED_MICE)
imp_s <- mice::mice(syn_dat[, cv15], m = 1, method = "pmm", maxit = MICE_MAXIT, seed = SEED_MICE, printFlag = FALSE)
Xs <- scale(mice::complete(imp_s, 1))
emb_s <- fit_umap(Xs); hc_s <- hclust(dist(emb_s), method = "ward.D2")
sil_s <- sapply(2:6, function(k) mean(cluster::silhouette(cutree(hc_s, k), dist(emb_s))[, 3]))
cat("Synthetic data — silhouette by k (2..6):", round(sil_s, 3), "\n")
cl_s <- cutree(hc_s, 3)
ins_s <- syn_dat$Diagnosis == "\u4e0d\u7720\u75c7"; osa_s <- !is.na(syn_dat$AHI) & syn_dat$AHI >= 5
lab_s <- ifelse(ins_s & osa_s, "COMISA", ifelse(ins_s, "Insomnia_only", ifelse(osa_s, "OSA_only", "Neither")))
cat("Synthetic cluster x label:\n"); print(table(Cluster = cl_s, Label = lab_s))
cat("Mean AHI / ISI by synthetic cluster:\n")
print(aggregate(cbind(AHI, ISI, DBAS) ~ cl_s, data = cbind(syn_dat, cl_s), FUN = function(x) round(mean(x, na.rm = TRUE), 1)))

# ---- 6. Write ----
if (!dir.exists("data")) dir.create("data")
write.csv(syn_dat, "data/synthetic_comisa.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ data/synthetic_comisa.csv written (n =", nrow(syn_dat), ") — review the checks above before committing.\n")
