# ============================================================
# 26_plausibility_check.R — 解析変数の妥当性チェック（クリニック照会リスト作成）
#  対象：comisa_raw（comisa.csv をそのまま読み込んだもの）
#  出力：
#   (1) 画面：変数ごとの該当件数のみ（集計値。Claude に貼ってよいのはこれだけ）
#   (2) output/revision/clinic_query_list.csv：ID・変数・値・許容範囲
#       （個人単位。メイツクリニックへの照会用。Claude には送らない）
# Run after: source("01_setup.R"); source("08_revision_setup.R")
# ============================================================
cat("\n=== 26: Plausibility check ===\n")
d <- comisa_raw

# 許容範囲（質問紙は採点上の範囲、PSG・身体計測は生理的にあり得る範囲）
RULES <- list(
  ISI   = c(0, 28),   DBAS = c(0, 160), FIRST = c(9, 36), HAS = c(0, 60),
  JESS  = c(0, 24),   SDS  = c(20, 80),
  Age   = c(18, 80),  BMI  = c(13, 60),
  AHI   = c(0, 150),  ODI3 = c(0, 150), Arousal_Index = c(0, 150),
  TST   = c(0, 720),  SE   = c(0, 1),   WASO = c(0, 720), SL = c(0, 480)
)
QUESTIONNAIRES <- c("ISI", "DBAS", "FIRST", "HAS", "JESS", "SDS")

rows <- list(); summ <- list()
for (v in names(RULES)) {
  x <- suppressWarnings(as.numeric(d[[v]])); r <- RULES[[v]]
  out_rng <- !is.na(x) & (x < r[1] | x > r[2])
  non_int <- if (v %in% QUESTIONNAIRES) !is.na(x) & abs(x - round(x)) > 1e-8 else rep(FALSE, length(x))
  flag <- out_rng | non_int
  summ[[v]] <- data.frame(Variable = v, Range = paste0(r[1], "–", r[2]),
                          N_observed = sum(!is.na(x)), Out_of_range = sum(out_rng),
                          Non_integer = sum(non_int))
  w <- which(flag)
  if (length(w)) rows[[v]] <- data.frame(ID = d$ID[w], Variable = v, Value = x[w],
                                         Allowed = paste0(r[1], "–", r[2]),
                                         Reason = ifelse(out_rng[w], "out of range", "non-integer total"))
}

# 論理チェック：受診経路（ID 接頭辞）と診断ラベルの整合
st  <- ifelse(grepl("^IS", d$ID), "IS", ifelse(grepl("^SA", d$ID), "SAS", NA))
dx  <- d$Diagnosis == "不眠症"
inc <- which((st == "IS" & !dx) | (st == "SAS" & dx))
summ[["Stream_vs_Diagnosis"]] <- data.frame(Variable = "Stream vs Diagnosis", Range = "IS = 不眠症",
                                            N_observed = sum(!is.na(st)), Out_of_range = length(inc), Non_integer = 0)
if (length(inc)) rows[["inc"]] <- data.frame(ID = d$ID[inc], Variable = "Diagnosis", Value = NA,
                                             Allowed = "IS→不眠症 / SA→不眠症以外", Reason = "stream/diagnosis mismatch")

S <- do.call(rbind, summ); rownames(S) <- NULL
cat("\n=== 該当件数（この表だけ共有してください） ===\n"); print(S)
if (length(rows)) {
  Q <- do.call(rbind, rows); rownames(Q) <- NULL
  cat(sprintf("\n照会対象：%d 件（患者 %d 人）\n", nrow(Q), length(unique(Q$ID))))
  f <- file.path(REV_DIR, "clinic_query_list.csv")
  write.csv(Q, f, row.names = FALSE, fileEncoding = "UTF-8")
  cat("✓ 照会リスト（ID 付き・共有不可）:", f, "\n")
} else cat("\n該当なし\n")
