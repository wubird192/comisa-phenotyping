# ============================================================
# 19_stream_adjusted.R — Disentangle era vs intake stream (follow-up to 13)
#  Finding from 13: 2010–2019 (n=84) is almost entirely C1, and the IS stream
#  (n=113) is 91% C1. Era and stream are confounded. This script:
#  (1) Stream × year table
#  (2) Era × cluster WITHIN the SAS stream only (2020–2024)
#  (3) Re-cluster the SAS stream alone (n=486): does the k=K_FINAL structure replicate?
#  (4) Cluster × stream adjusted for AHI/ISI (logistic: P(C1) ~ stream + z-ISI + z-AHI)
#  (5) COMISA-in-C1 (AHI>=5) by stream
# ============================================================
cat("\n=== 19: Stream-adjusted analyses ===\n")
cl  <- factor(REF_CL); yr <- as.integer(comisa_raw[[YEAR_COL]])
st  <- factor(comisa_raw[[STREAM_COL]])
era3 <- cut(yr, c(2009, 2014, 2019, 2024), labels = c("2010-2014", "2015-2019", "2020-2024"))

t1 <- table(Stream = st, Year = yr); print(t1); save_csv(as.data.frame(t1), "TableS9b_stream_by_year.csv")

sas <- which(st == "SAS")
t2 <- table(Era = droplevels(era3[sas]), Cluster = cl[sas]); print(t2); print(round(100*prop.table(t2,1),1))
if (nrow(t2) > 1) { ch <- suppressWarnings(chisq.test(t2)); cat(sprintf("Within SAS: chi2=%.1f df=%d p=%.3g V=%.3f\n", ch$statistic, ch$parameter, ch$p.value, sqrt(ch$statistic/(sum(t2)*(min(dim(t2))-1))))) }
t2y <- table(Year = yr[sas], Cluster = cl[sas]); print(round(100*prop.table(t2y,1),1))
save_csv(as.data.frame(t2y), "TableS9b_cluster_by_year_within_SAS.csv")

# (3) Re-cluster SAS stream alone
emb <- fit_umap(REF_SCALED[sas, ]); hc <- hclust(dist(emb), method = "ward.D2")
sil <- sapply(2:6, function(k) mean(cluster::silhouette(cutree(hc, k), dist(emb))[, 3]))
cat("SAS-only silhouette k=2..6:", round(sil, 3), "\n")
cl_s <- cutree(hc, K_FINAL)
cat(sprintf("SAS-only ARI vs reference (k=%d): %.3f\n", K_FINAL, ari(REF_CL[sas], cl_s)))
m <- match_labels(REF_CL[sas], cl_s); print(table(Reference = REF_CL[sas], SAS_only = m))
ins <- comisa$Diagnosis == "不眠症"; com5 <- ins & !is.na(comisa$AHI) & comisa$AHI >= 5
cat(sprintf("SAS-only: COMISA(AHI>=5) in matched C1 = %.1f%% (n=%d)\n", 100*mean(m[com5[sas]] == 1), sum(com5[sas])))
save_csv(data.frame(k = 2:6, silhouette = sil), "TableS9b_SAS_only_silhouette.csv")

# (4) Is stream associated with C1 beyond ISI/AHI?
z <- function(x) as.numeric(scale(x))
d <- data.frame(C1 = as.integer(cl == 1), Stream = st, ISI = z(comisa$ISI), AHI = z(comisa$AHI), DBAS = z(comisa$DBAS))
f0 <- glm(C1 ~ Stream, d, family = binomial); f1 <- glm(C1 ~ Stream + ISI + AHI, d, family = binomial)
f2 <- glm(C1 ~ Stream + ISI + AHI + DBAS, d, family = binomial)
or <- function(f) { s <- summary(f)$coefficients; data.frame(Term = rownames(s), OR = round(exp(s[,1]),2), p = signif(s[,4],3)) }
cat("Crude:\n"); print(or(f0)); cat("Adjusted ISI+AHI:\n"); print(or(f1)); cat("Adjusted +DBAS:\n"); print(or(f2))
save_csv(rbind(cbind(Model="crude", or(f0)), cbind(Model="ISI+AHI", or(f1)), cbind(Model="ISI+AHI+DBAS", or(f2))), "TableS9b_stream_logistic.csv")

# (5) COMISA in C1 by stream
for (s in levels(st)) { i <- com5 & st == s; cat(sprintf("%s: COMISA n=%d, in C1 %.1f%%\n", s, sum(i), 100*mean(cl[i]==1))) }
