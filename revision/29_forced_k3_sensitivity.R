# ============================================================
# 29_forced_k3_sensitivity.R — 訂正後データで k=3 を強制した感度分析
#  目的：事前規定の規則では k=2 が選ばれる（訂正後）。投稿版の k=3 構造が
#        訂正後データでどう再現されるかを示す。
#  出力（集計のみ、output/revision/）：
#    29_forced_k3_sizes_indices.csv  強制 k=3 のサイズ・Silhouette/CH/DB
#    29_forced_k3_vs_new_k2.csv      強制 k=3 × 主解析 k=2 のクロス表
#    29_forced_k3_vs_old_k3.csv      強制 k=3 × 投稿版 k=3（output_old）のクロス表＋ARI
#    29_forced_k3_COMISA_distribution.csv  COMISA（AHI≥5／≥15）のクラスタ分布
# Run after: source("01_setup.R"); source("08_revision_setup.R")
# ============================================================
cat("\n=== 29: Forced k=3 sensitivity on corrected data ===\n")
suppressPackageStartupMessages({ library(cluster); library(mclust) })

# ---- 設定 ----------------------------------------------------------------------
EMB_OBJ  <- "umap_mat"   # 最終 UMAP 埋め込み（N×2）。08/09 が使う行列。NULL なら自動検出（emb_s は別物なので不可）
OLD_FILE <- file.path(dirname(OUTPUT_DIR), "output_old", "Cluster_Assignments.csv")
NEW_FILE <- file.path(OUTPUT_DIR, "Cluster_Assignments.csv")
HC_METHOD <- "ward.D2"
# ------------------------------------------------------------------------------

if (is.null(EMB_OBJ)) {
  cand <- Filter(function(o) {
    x <- get(o); (is.matrix(x) || is.data.frame(x)) && nrow(x) == N && ncol(x) == 2 && all(sapply(as.data.frame(x), is.numeric))
  }, ls(envir = .GlobalEnv))
  cat("N x 2 numeric objects found:", paste(cand, collapse = ", "), "\n")
  if (length(cand) == 0) stop("埋め込みが見つかりません。EMB_OBJ に名前を指定してください（09 で使っている埋め込み）。")
  EMB_OBJ <- cand[grepl("umap|emb", cand, ignore.case = TRUE)][1]
  if (is.na(EMB_OBJ)) EMB_OBJ <- cand[1]
}
emb <- as.matrix(get(EMB_OBJ)); cat("Using embedding object:", EMB_OBJ, "\n")

# 1. 強制 k=3（主解析と同じ HC：ward.D2、埋め込み上のユークリッド距離）
d   <- dist(emb)
hc  <- hclust(d, method = HC_METHOD)
cl3 <- cutree(hc, k = 3)
# クラスタ番号を AHI 中央値の昇順に振り直す（C1 = 低 AHI）
ord <- order(tapply(comisa$AHI, cl3, median, na.rm = TRUE)); cl3 <- match(cl3, ord)
sil <- mean(silhouette(cl3, d)[, 3])
idx <- data.frame(k = 3, C1 = sum(cl3 == 1), C2 = sum(cl3 == 2), C3 = sum(cl3 == 3),
                  Silhouette_embedding = round(sil, 3),
                  CH_embedding = round(clusterSim::index.G1(emb, cl3), 1),
                  DB_embedding = round(clusterSim::index.DB(emb, cl3)$DB, 3),
                  ARI_vs_primary_k2 = round(adjustedRandIndex(cl3, REF_CL), 3))
print(idx); save_csv(idx, "29_forced_k3_sizes_indices.csv")

# 2. 強制 k=3 × 主解析 k=2
t_new <- table(Forced_k3 = cl3, Primary_k2 = REF_CL); print(t_new)
save_csv(cbind(Forced_k3 = rownames(t_new), as.data.frame.matrix(t_new)), "29_forced_k3_vs_new_k2.csv")

# 3. 強制 k=3 × 投稿版 k=3
if (file.exists(OLD_FILE)) {
  # Cluster_Assignments.csv の ID は連番（comisa$ID とは別体系）。行順で対応させ、
  # 新ファイルの Cluster が REF_CL と一致することで行順の同一性を確認する
  newf <- read.csv(NEW_FILE); oldf <- read.csv(OLD_FILE)
  stopifnot(nrow(newf) == N, nrow(oldf) == N, all(newf$Cluster == REF_CL))
  old <- oldf[match(newf$ID, oldf$ID), ]
  stopifnot(!any(is.na(old$Cluster)))
  t_old <- table(Forced_k3 = cl3, Submitted_k3 = old$Cluster); print(t_old)
  ari_old <- adjustedRandIndex(cl3, old$Cluster)
  cat(sprintf("ARI forced k=3 vs submitted k=3: %.3f\n", ari_old))
  tdf <- cbind(Forced_k3 = rownames(t_old), as.data.frame.matrix(t_old)); tdf$ARI <- round(ari_old, 3)
  save_csv(tdf, "29_forced_k3_vs_old_k3.csv")
} else cat("output_old が見つからないため投稿版との比較はスキップ\n")

# 4. COMISA の分布
INS <- comisa$Diagnosis == "不眠症"
res <- do.call(rbind, lapply(c(5, 15), function(thr) {
  com <- INS & !is.na(comisa$AHI) & comisa$AHI >= thr
  ci  <- binom.test(sum(cl3[com] == 1), sum(com))
  data.frame(Threshold = thr, N_COMISA = sum(com), C1 = sum(cl3[com] == 1), C2 = sum(cl3[com] == 2), C3 = sum(cl3[com] == 3),
             Pct_C1 = round(100 * ci$estimate, 1), CI_low = round(100 * ci$conf.int[1], 1), CI_high = round(100 * ci$conf.int[2], 1))
}))
print(res); save_csv(res, "29_forced_k3_COMISA_distribution.csv")
cat("✓ 29 done\n")
