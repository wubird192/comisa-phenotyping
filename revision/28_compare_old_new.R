# 28_compare_old_new.R
# 目的：データ訂正（FIRST 59→23、SDS 86→46）の前後でクラスタ所属・指標がどう変わったかを集計する。
#       出力は集計値のみ（個人単位の行は書き出さない）。回答書「Additional corrections」に引用する。
# 使い方：先生の Mac で、再実行前に旧出力を output_old/ に退避してから、再実行後に本スクリプトを実行。
#   1) 旧 Cluster_Assignments.csv を output_old/ にコピー（再実行前に！）
#   2) 00–07、RUN_ALL_REVISION.R、23・24 を再実行
#   3) source("28_compare_old_new.R")
# ---- 設定 ---------------------------------------------------------------------------
if (!exists("DATA_PATH") || !exists("OUTPUT_DIR")) stop("Run source(\"01_setup.R\") first (defines DATA_PATH and OUTPUT_DIR).")
BASE <- dirname(OUTPUT_DIR)   # folder that contains output/ (corrected run) and output_old/ (run on the uncorrected data)
OLD_FILE <- file.path(BASE, "output_old/Cluster_Assignments.csv")   # 再実行前（退避コピー）
NEW_FILE <- file.path(BASE, file.path(OUTPUT_DIR, "Cluster_Assignments.csv"))       # 再実行後
OUT_DIR  <- file.path(BASE, file.path(OUTPUT_DIR, "revision"))
ID_COL   <- NULL   # NULL なら自動検出（"id" を含む列名）。違えば "patient_id" のように指定
CL_COL   <- NULL   # NULL なら自動検出（"cluster" を含む列名）
# ------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(mclust) })

old <- read.csv(OLD_FILE, check.names = FALSE)
new <- read.csv(NEW_FILE, check.names = FALSE)
cat("columns:", paste(names(new), collapse = ", "), "\n")

pick <- function(cols, pattern, given) {
  if (!is.null(given)) return(given)
  hit <- grep(pattern, cols, ignore.case = TRUE, value = TRUE)
  if (length(hit) != 1) stop(sprintf("列を特定できません（候補: %s）。設定で指定してください。",
                                     paste(hit, collapse = ", ")))
  hit
}
ID_COL <- pick(names(new), "^id$|_id$|^id_|patient", ID_COL)
CL_COL <- pick(names(new), "cluster", CL_COL)
cat(sprintf("ID_COL = %s, CL_COL = %s\n", ID_COL, CL_COL))
stopifnot(ID_COL %in% names(old), CL_COL %in% names(old))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

m <- merge(old[, c(ID_COL, CL_COL)], new[, c(ID_COL, CL_COL)],
           by = ID_COL, suffixes = c("_old", "_new"))
cat(sprintf("n old = %d, n new = %d, matched = %d\n", nrow(old), nrow(new), nrow(m)))

co <- factor(m[[paste0(CL_COL, "_old")]]); cn <- factor(m[[paste0(CL_COL, "_new")]])

# 1. 新旧クロス表（ラベルの入れ替わりがあっても対角で読めるよう、そのまま出力）
tab <- table(old = co, new = cn)
write.csv(as.data.frame.matrix(tab), file.path(OUT_DIR, "28_cluster_crosstab_old_new.csv"))

# 2. 要約指標
ari   <- mclust::adjustedRandIndex(co, cn)
# ラベル対応を最適化したうえでの一致率（ハンガリー法：clue）
# k が新旧で違ってもよいように、旧クラスタごとに「最も多く流れ込んだ新クラスタ」を対応先とみなす
# （旧 C2 と旧 C3 が同じ新クラスタに対応することを許す）。一致率 = 対応先に入った人の割合
mt <- as.matrix(tab)
best_new <- apply(mt, 1, which.max)
agree <- sum(mt[cbind(seq_len(nrow(mt)), best_new)]) / sum(mt)
moved <- sum(mt) - sum(mt[cbind(seq_len(nrow(mt)), best_new)])
mapping <- data.frame(old_cluster = rownames(mt), maps_to_new = colnames(mt)[best_new],
                      n_old = rowSums(mt), n_in_mapped = mt[cbind(seq_len(nrow(mt)), best_new)])
mapping$pct_retained <- round(100 * mapping$n_in_mapped / mapping$n_old, 1)
write.csv(mapping, file.path(OUT_DIR, "28_cluster_mapping_old_to_new.csv"), row.names = FALSE)
print(mapping, row.names = FALSE)

summ <- data.frame(
  metric = c("n_matched", "ARI_old_vs_new", "agreement_under_best_mapping", "n_members_outside_mapping",
             paste0("n_old_C", levels(co)), paste0("n_new_C", levels(cn))),
  value  = c(nrow(m), round(ari, 3), round(agree, 3), moved,
             as.integer(table(co)), as.integer(table(cn)))
)
write.csv(summ, file.path(OUT_DIR, "28_compare_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

# 3. 訂正した 2 名の新旧所属（ID は出力しない。行番号も出さず、所属の組だけ）
#    訂正対象 ID をここに入れる（パスワード付きファイルから。コミットしない）
CORRECTED_IDS <- c()   # 例: c("A123", "B456")
if (length(CORRECTED_IDS) > 0) {
  mm <- m[m[[ID_COL]] %in% CORRECTED_IDS, ]
  mv <- data.frame(case = seq_len(nrow(mm)),
                   cluster_old = mm[[paste0(CL_COL, "_old")]],
                   cluster_new = mm[[paste0(CL_COL, "_new")]])
  write.csv(mv, file.path(OUT_DIR, "28_corrected_cases_cluster.csv"), row.names = FALSE)
  print(mv, row.names = FALSE)
}

# 4. 集計 CSV の新旧比較：output_old/ と output/（revision/ を含む）で同名の CSV を総当たりで突き合わせ、
#    数値列の差分を 1 ファイルにまとめる（個人単位ファイルは行数 > 100 で除外）。
old_dir <- file.path(BASE, "output_old"); new_dir <- file.path(BASE, OUTPUT_DIR)
files <- intersect(list.files(old_dir, "\\.csv$", recursive = TRUE),
                   list.files(new_dir, "\\.csv$", recursive = TRUE))
diffs <- list()
for (f in files) {
  a <- tryCatch(read.csv(file.path(old_dir, f), check.names = FALSE), error = function(e) NULL)
  b <- tryCatch(read.csv(file.path(new_dir, f), check.names = FALSE), error = function(e) NULL)
  if (is.null(a) || is.null(b) || nrow(a) > 100 || !identical(dim(a), dim(b))) next
  num <- names(b)[sapply(b, is.numeric) & sapply(a, is.numeric)]
  for (v in num) {
    d <- which(abs(a[[v]] - b[[v]]) > 1e-9 | xor(is.na(a[[v]]), is.na(b[[v]])))
    if (length(d) == 0) next
    lab <- if (ncol(a) > 0 && !is.numeric(a[[1]])) as.character(a[[1]][d]) else as.character(d)
    diffs[[length(diffs) + 1]] <- data.frame(file = f, column = v, row = lab,
                                             old = a[[v]][d], new = b[[v]][d])
  }
}
if (length(diffs)) {
  dd <- do.call(rbind, diffs)
  write.csv(dd, file.path(OUT_DIR, "28_summary_csv_diffs.csv"), row.names = FALSE)
  cat(sprintf("summary CSVs with changed values: %d files, %d cells\n",
              length(unique(dd$file)), nrow(dd)))
} else cat("no differences found in summary CSVs\n")

cat("done: ", file.path(OUT_DIR, "28_compare_summary.csv"), "\n")
