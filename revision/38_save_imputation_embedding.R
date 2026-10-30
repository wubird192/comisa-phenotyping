# ==============================================================================
# 38_save_imputation_embedding.R — mice オブジェクトと最終設定の UMAP 埋め込みを保存（検証付き）
#   前提:  の 01_setup.R, 02_preprocess.R, 08_revision_setup.R
#   出力: output/mice_imputation.rds, output/UMAP_Embedding_Final.csv（個人単位。共有・GitHub 不可）
#   検証: 埋め込みを Ward k=2 で切ると Cluster_Assignments.csv と ARI = 1 になること
# ==============================================================================
if (!exists("DATA_PATH") || !exists("OUTPUT_DIR")) stop("Run source(\"01_setup.R\") first (defines DATA_PATH and OUTPUT_DIR).")
REV <- "."   # the revision scripts are in this folder
rm(list = intersect(c("imp", "emb", "imputed"), ls()))

# ---- 1. mice オブジェクト（01 → 02 を再実行。seed 固定なので同じ補完） ----
source(file.path(REV, "01_setup.R"))
source(file.path(REV, "02_preprocess.R"))
stopifnot(exists("imputed"), inherits(imputed, "mids"), imputed$m == 5, nrow(imputed$data) == 599)
saveRDS(imputed, file.path(OUTPUT_DIR, "mice_imputation.rds"))
cat("saved output/mice_imputation.rds  (m =", imputed$m, ", n =", nrow(imputed$data), ")\n")

# ---- 2. 参照埋め込み（08 が REF_SCALED / REF_CL / fit_umap / hc_cut を用意する） ----
source(file.path(REV, "08_revision_setup.R"))
stopifnot(exists("REF_SCALED"), exists("REF_CL"), exists("fit_umap"), exists("hc_cut"), exists("K_FINAL"))
emb_ref <- fit_umap(REF_SCALED)                     # 既定引数 = 最終設定（manhattan, 20, 0.01, seed）
if (!is.matrix(emb_ref)) emb_ref <- emb_ref$embedding
stopifnot(nrow(emb_ref) == 599, ncol(emb_ref) >= 2)

cl_chk <- hc_cut(emb_ref, K_FINAL)
cl_csv <- read.csv(file.path(OUTPUT_DIR, "Cluster_Assignments.csv"))
clc    <- intersect(c("Cluster", "cluster", "Cluster_k2"), names(cl_csv))[1]
ari1 <- mclust::adjustedRandIndex(cl_chk, REF_CL)
ari2 <- mclust::adjustedRandIndex(cl_chk, cl_csv[[clc]])
print(table(embedding_k2 = cl_chk, Cluster_Assignments = cl_csv[[clc]]))
cat(sprintf("ARI vs REF_CL = %.3f ; ARI vs Cluster_Assignments.csv = %.3f\n", ari1, ari2))
if (ari1 < 1 || ari2 < 1) stop("埋め込みが主解析を再現していません。保存しません。上の表を Claude に貼ってください。")

write.csv(data.frame(UMAP1 = emb_ref[, 1], UMAP2 = emb_ref[, 2]),
          file.path(OUTPUT_DIR, "UMAP_Embedding_Final.csv"), row.names = FALSE)
cat("saved output/UMAP_Embedding_Final.csv\n")
cat("次: source('36_main_figures_v3.R') で Fig 1、35_v3 で FigS6 が出ます。\n")
