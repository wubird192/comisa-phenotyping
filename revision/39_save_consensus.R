# ==============================================================================
# 39_save_consensus.R — 05_stability.R を 08 の環境で再実行し、コンセンサス行列を保存
#   所要: 200 回の UMAP ≈ 20–40 分。seed 固定なので Cluster_Stability_Summary.csv は同じ値で上書きされる
#   出力: output/Consensus_Matrix.rds（599×599、個人単位。共有・GitHub 不可）
#   検証: 全体平均 0.815、C1 0.848、C2 0.783（SUMMARY v2 §8）と一致すること
# ==============================================================================
if (!exists("DATA_PATH") || !exists("OUTPUT_DIR")) stop("Run source(\"01_setup.R\") first (defines DATA_PATH and OUTPUT_DIR).")
REV <- "."   # the revision scripts are in this folder

source(file.path(REV, "01_setup.R"))
source(file.path(REV, "08_revision_setup.R"))          # final_result, REF_CL, SEED_* などが揃う
stopifnot(exists("final_result"), exists("BOOT_N_ITER"), exists("SEED_BOOT"))

res <- tryCatch(source(file.path(REV, "05_stability.R")), error = function(e) { cat("05 でエラー:", conditionMessage(e), "\n"); NULL })
stopifnot(exists("consensus_mat"))                      # 05 の 74 行目で作られる
stopifnot(nrow(consensus_mat) == 599)
saveRDS(consensus_mat, file.path(OUTPUT_DIR, "Consensus_Matrix.rds"))
cat("saved output/Consensus_Matrix.rds\n")

# ---- 検証 ----
sc <- sapply(sort(unique(REF_CL)), function(k) { i <- which(REF_CL == k); m <- consensus_mat[i, i]; mean(m[upper.tri(m)]) })
cat(sprintf("consensus C1 %.3f, C2 %.3f, overall %.3f  (期待 0.848 / 0.783 / 0.815)\n", sc[1], sc[2], mean(sc)))
cat("次: source('35_supp_figures_v3.R') でコンセンサス図が出ます。\n")
