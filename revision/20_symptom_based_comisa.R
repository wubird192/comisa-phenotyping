# ============================================================
# 20_symptom_based_comisa.R — Stream-independent (symptom-based) reference
#  Finding from 19: Diagnosis=="不眠症" coincides exactly with the IS intake
#  stream; no SAS-stream patient carries an insomnia diagnosis. The clinical
#  label is therefore presentation-based. Here COMISA is defined by SYMPTOMS:
#     insomnia_sx = ISI >= 15 ;  OSA = AHI >= 5 (and >= 15 as sensitivity)
#  (1) Symptom-based cross-tabs and C1 capture, overall and by stream
#  (2) Profile of SAS-stream patients in C1 vs other SAS patients
#  (3) Table-3-style metrics for symptom-based labels
# ============================================================
cat("\n=== 20: Symptom-based COMISA ===\n")
cl <- REF_CL; st <- comisa_raw[[STREAM_COL]]
isi_hi <- !is.na(comisa$ISI) & comisa$ISI >= 15
for (thr in c(5, 15)) {
  osa <- !is.na(comisa$AHI) & comisa$AHI >= thr
  lab <- ifelse(isi_hi & osa, "COMISA_sx", ifelse(isi_hi, "Insomnia_sx_only", ifelse(osa, "OSA_only", "Neither")))
  cat(sprintf("\n--- ISI>=15 & AHI>=%d ---\n", thr)); t <- table(Cluster = cl, Label = lab); print(addmargins(t))
  ci <- binom.test(sum(cl[lab == "COMISA_sx"] == 1), sum(lab == "COMISA_sx"))
  cat(sprintf("COMISA_sx n=%d: in C1 %.1f%% (%.1f-%.1f)\n", sum(lab=="COMISA_sx"), 100*ci$estimate, 100*ci$conf.int[1], 100*ci$conf.int[2]))
  for (s in c("IS","SAS")) { i <- lab=="COMISA_sx" & st==s; cat(sprintf("  %s stream: n=%d, in C1 %.1f%%\n", s, sum(i), 100*mean(cl[i]==1))) }
  m <- metrics((cl == 1), (lab == "COMISA_sx"))   # from 17 (full cohort)
  print(m[, c("N","AUC","Sensitivity","Specificity","PPV","NPV")])
  save_csv(cbind(Threshold = thr, as.data.frame.matrix(t)), sprintf("Table_symptomCOMISA_crosstab_AHI%d.csv", thr))
  save_csv(cbind(Threshold = thr, m), sprintf("Table_symptomCOMISA_metrics_AHI%d.csv", thr))
}
# (2) SAS-stream patients in C1: who are they?
sasC1 <- st=="SAS" & cl==1; sasOther <- st=="SAS" & cl!=1
prof <- do.call(rbind, lapply(c("ISI","DBAS","FIRST","HAS","SDS","AHI","ODI3","Arousal_Index","Age","BMI"), function(v) {
  d <- effsize::cohen.d(comisa[[v]][sasC1], comisa[[v]][sasOther], na.rm=TRUE)
  data.frame(Variable=v, SAS_C1_mean=mean(comisa[[v]][sasC1],na.rm=TRUE), SAS_other_mean=mean(comisa[[v]][sasOther],na.rm=TRUE), d=d$estimate) }))
print(prof, digits=3); save_csv(prof, "Table_SAS_C1_profile.csv")
cat(sprintf("SAS-stream C1 (n=%d): ISI>=15 in %.1f%%; ISI>=10 in %.1f%%; AHI>=5 in %.1f%%\n", sum(sasC1),
  100*mean(comisa$ISI[sasC1]>=15,na.rm=TRUE), 100*mean(comisa$ISI[sasC1]>=10,na.rm=TRUE), 100*mean(comisa$AHI[sasC1]>=5,na.rm=TRUE)))
cat(sprintf("SAS-stream other (n=%d): ISI>=15 in %.1f%%\n", sum(sasOther), 100*mean(comisa$ISI[sasOther]>=15,na.rm=TRUE)))
