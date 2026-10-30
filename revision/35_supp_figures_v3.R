# ==============================================================================
# 35_supp_figures_v3.R — 補足図 4 点をラベル修正のうえ journal-figures 規範で再描画
#   FigS6_Imputation_densities_v2   : DBAS・FIRST・HAS の 3 尺度のみ（観測 vs 補完 5 組）
#   FigS3_UMAP_grid_ARI_v2          : タイトルを「best-k solution vs reference」に
#   FigureS_Consensus_Matrix_v2     : 「bootstrap」→「200 subsamples of 80%」
#   FigS13_recovery_simulation_v2   : 凡例をコード名から手法名に
#   出力: output/revision/figs/*.png (600 dpi) と *.pdf（ベクタ）
# ------------------------------------------------------------------------------
# 規範（journal-figures）: 白背景、Arial/Helvetica、theme_classic、グリッドなし、
#   上・右の軸線なし、色は系列数ぶんだけ（JAMA palette：slate が主、非焦点は grey）、
#   凡例は上に枠なし or 直接ラベル、装飾なし。
# ⟨要確認⟩ と書いた入力パス・列名は先生の環境に合わせて直してください。
#   走らせる前に  list.files(file.path(OUTPUT_DIR, "revision"), pattern = "TableS13|Consensus")  で確認。
# ==============================================================================

if (!exists("DATA_PATH") || !exists("OUTPUT_DIR")) stop("Run source(\"01_setup.R\") first (defines DATA_PATH and OUTPUT_DIR).")
set.seed(20261020)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

FIGDIR <- file.path(OUTPUT_DIR, "revision/figs"); dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ---- JAMA palette (ggsci 実測値) --------------------------------------------
JAMA <- c(slate = "#374E55", gold = "#DF8F44", cerulean = "#00A1D5", red = "#B24745",
          sage = "#79AF97", violet = "#6A6599", grey = "#80796B")

theme_journal <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          axis.text = element_text(size = base_size - 1, colour = "black"),
          axis.title = element_text(size = base_size),
          plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0), plot.tag = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = base_size - 1, hjust = 0, colour = "grey30"),
          strip.background = element_blank(), strip.text = element_text(size = base_size, face = "bold", hjust = 0),
          legend.position = "top", legend.title = element_blank(), legend.key = element_blank(),
          legend.text = element_text(size = base_size - 1), legend.key.height = unit(3, "mm"),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
}
save_fig <- function(p, name, w, h) {               # w, h in mm（印刷実寸）
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(FIGDIR, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat("wrote", file.path(FIGDIR, name), ".pdf/.png\n")
}
# 列名を候補から拾う小道具
pick_col <- function(df, cands) { x <- intersect(cands, names(df)); if (!length(x)) stop("列が見つかりません: ", paste(cands, collapse = "/"), "\n  ある列: ", paste(names(df), collapse = ", ")); x[1] }
pick_file <- function(dir, pattern) { f <- list.files(dir, pattern = pattern, full.names = TRUE); if (!length(f)) stop("ファイルなし: ", pattern, " in ", dir); f[1] }

# ==============================================================================
# 1. FigS6 — 観測値 vs 補完値の密度（DBAS・FIRST・HAS のみ）
# ==============================================================================
IMP_RDS <- file.path(OUTPUT_DIR, "mice_imputation.rds")   # 02 の mice() 直後に saveRDS(imp, file.path(OUTPUT_DIR, "mice_imputation.rds")) を足して再実行
if (file.exists(IMP_RDS)) {
  imp  <- readRDS(IMP_RDS)
  long <- mice::complete(imp, "long", include = TRUE)        # .imp = 0 が元データ
  vars <- c("DBAS", "FIRST", "HAS")
  miss <- is.na(long[long$.imp == 0, vars])                  # 元データの欠測位置
  dens <- bind_rows(lapply(vars, function(v) {
    obs <- data.frame(Variable = v, Set = "Observed", Value = long[long$.imp == 0, v][!miss[, v]])
    imps <- bind_rows(lapply(1:imp$m, function(m) {
      data.frame(Variable = v, Set = paste0("Imputed ", m), Value = long[long$.imp == m, v][miss[, v]])
    }))
    bind_rows(obs, imps)
  }))
  dens$Variable <- factor(dens$Variable, levels = vars,
                          labels = c("DBAS (0–160)", "FIRST (9–36)", "HAS (0–78)"))
  dens$Set <- factor(dens$Set, levels = c("Observed", paste0("Imputed ", 1:imp$m)))
  n_lab <- dens %>% group_by(Variable, Set) %>% summarise(n = n(), .groups = "drop") %>%
    filter(Set %in% c("Observed", "Imputed 1")) %>%
    mutate(lab = ifelse(Set == "Observed", paste0("observed n = ", n), paste0("imputed n = ", n, " (× ", imp$m, ")"))) %>%
    group_by(Variable) %>% summarise(lab = paste(lab, collapse = "\n"), .groups = "drop")
  p6 <- ggplot(dens, aes(Value, colour = Set, linewidth = Set, group = Set)) +
    geom_density(adjust = 1.2, key_glyph = "path") +
    scale_colour_manual(values = c(Observed = JAMA[["slate"]], setNames(rep("#B9B5AE", imp$m), paste0("Imputed ", 1:imp$m))),
                        breaks = c("Observed", "Imputed 1"), labels = c("Observed values", "Imputed values (5 datasets)")) +
    scale_linewidth_manual(values = c(Observed = 0.8, setNames(rep(0.35, imp$m), paste0("Imputed ", 1:imp$m))), guide = "none") +
    geom_text(data = n_lab, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.3, size = 2.4, colour = "grey30", lineheight = 0.9) +
    facet_wrap(~ Variable, scales = "free", nrow = 1) +
    labs(x = "Score", y = "Density") + theme_journal()
  save_fig(p6, "FigS6_Imputation_densities_v2", 180, 65)
} else cat("[skip] FigS6: ", IMP_RDS, " が見つかりません。IMP_RDS を直してください。\n")

# ==============================================================================
# 2. FigS3 — UMAP グリッド 40 設定：best-k 解 vs 参照解の ARI
# ==============================================================================
prm <- read.csv(file.path(OUTPUT_DIR, "UMAP_Final_Parameters.csv")); getp <- function(k) prm$Value[prm$Parameter == k]
REF_METRIC <- tolower(getp("Metric")); REF_NN <- as.numeric(getp("N_neighbors")); REF_MD <- as.numeric(getp("Min_dist")); OBS_SIL <- as.numeric(getp("Silhouette"))
g <- read.csv(file.path(OUTPUT_DIR, "revision/TableS3_UMAP_grid_ARI.csv"))
ari_col <- pick_col(g, c("ARI_bestk_vs_reference", "ARI_k3_vs_reference"))   # 17 で改名予定
g2 <- data.frame(metric = tools::toTitleCase(g$metric), n_neighbors = factor(g$n_neighbors),
                 min_dist = factor(g$min_dist), best_k = g$best_k_by_silhouette, ARI = g[[ari_col]],
                 ref = tolower(g$metric) == REF_METRIC & g$n_neighbors == REF_NN & abs(g$min_dist - REF_MD) < 1e-9)
stopifnot(nrow(g2) == 40, sum(g2$ref) == 1)
p3 <- ggplot(g2, aes(min_dist, n_neighbors)) +
  geom_tile(aes(fill = ARI), colour = "white", linewidth = 0.6) +
  geom_tile(data = filter(g2, ref), fill = NA, colour = JAMA[["red"]], linewidth = 0.8) +
  geom_text(aes(label = best_k, colour = ARI > 0.55), size = 2.8, family = "Helvetica") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradient(low = "#F2F2F0", high = JAMA[["slate"]], limits = c(0, 1), breaks = c(0, 0.5, 1),
                      name = "ARI vs reference") +
  facet_wrap(~ metric, nrow = 1) +
  labs(title = "ARI of each configuration's best-k solution against the reference solution",
       subtitle = sprintf("Cell label = silhouette-optimal k; red outline = reference configuration (%s, %g, %g)", tools::toTitleCase(REF_METRIC), REF_NN, REF_MD),
       x = "min_dist", y = "n_neighbors") +
  theme_journal() + theme(legend.position = "right", legend.title = element_text(size = 8),
                          legend.key.height = unit(6, "mm"), legend.key.width = unit(3, "mm"),
                          axis.line = element_blank(), axis.ticks = element_blank())
save_fig(p3, "FigS3_UMAP_grid_ARI_v2", 170, 75)

# ==============================================================================
# 3. コンセンサス行列（200 subsamples of 80%）
# ==============================================================================
CONS_FILE <- file.path(OUTPUT_DIR, "Consensus_Matrix.rds")   # 安定性スクリプトでコンセンサス行列を saveRDS(<行列>, file.path(OUTPUT_DIR, "Consensus_Matrix.rds"))
if (file.exists(CONS_FILE)) {
  M <- if (grepl("\\.rds$", CONS_FILE)) readRDS(CONS_FILE) else as.matrix(read.csv(CONS_FILE, row.names = 1))
  stopifnot(nrow(M) == 599, ncol(M) == 599)
  cl  <- read.csv(file.path(OUTPUT_DIR, "Cluster_Assignments.csv")); clcol <- pick_col(cl, c("Cluster", "cluster", "Cluster_k2"))
  ord <- order(cl[[clcol]], -rowMeans(M))                   # クラスタ順 → クラスタ内は平均一致率の高い順
  Mo  <- M[ord, ord]; sizes <- table(cl[[clcol]]); cuts <- head(cumsum(sizes), -1); mids <- (c(0, cuts) + cumsum(sizes)) / 2
  long <- as.data.frame(as.table(Mo)); names(long) <- c("i", "j", "Consensus")
  long$i <- as.integer(long$i); long$j <- as.integer(long$j)
  pc <- ggplot(long, aes(j, i, fill = Consensus)) + geom_raster() +
    scale_fill_gradient(low = "white", high = JAMA[["slate"]], limits = c(0, 1), breaks = c(0, 0.5, 1),
                        name = "Consensus") +
    geom_hline(yintercept = cuts + 0.5, colour = JAMA[["red"]], linewidth = 0.4) +
    geom_vline(xintercept = cuts + 0.5, colour = JAMA[["red"]], linewidth = 0.4) +
    annotate("text", x = mids, y = -12, label = paste("Cluster", names(sizes)),
             size = 2.8, family = "Helvetica") +
    scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + coord_fixed(clip = "off") +
    labs(title = "Consensus matrix: 200 subsamples of 80% of patients",
         x = "Patients, ordered by cluster", y = "Patients, ordered by cluster") +
    theme_journal() + theme(legend.position = "right", legend.title = element_text(size = 8),
                            legend.key.height = unit(6, "mm"), legend.key.width = unit(3, "mm"),
                            axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
                            plot.margin = margin(4, 4, 4, 4, "mm"))
  save_fig(pc, "FigureS_Consensus_Matrix_v2", 120, 110)
} else cat("[skip] Consensus: ", CONS_FILE, " が見つかりません。CONS_FILE を直してください。\n")

# ==============================================================================
# 4. FigS13 — 回収シミュレーション（3 パネル；凡例は手法名）
#    入力: TableS13_recovery_simulation_raw.csv (rep, k_true, delta, Method, k_selected, ARI_at_true_k, Sil_k2, Sil_k3)
#          TableS13_recovery_simulation_summary.csv (k_true, delta, Method, ARI_at_true_k, ..., Pct_k_correct)
# ==============================================================================
METHOD_LAB <- function(x) {                        # コード名 → 手法名
  x <- as.character(x)
  dplyr::case_when(grepl("umap", x, ignore.case = TRUE) ~ "UMAP + Ward's",
                   grepl("pca",  x, ignore.case = TRUE) ~ "Ward's on PC1–2",
                   TRUE                                  ~ "Ward's, original space")
}
METHOD_LEVELS <- c("UMAP + Ward's", "Ward's, original space", "Ward's on PC1–2")
METHOD_COLS   <- setNames(c(JAMA[["slate"]], JAMA[["grey"]], JAMA[["gold"]]), METHOD_LEVELS)

raw13 <- read.csv(file.path(OUTPUT_DIR, "revision/TableS13_recovery_simulation_raw.csv"))
sum13 <- read.csv(file.path(OUTPUT_DIR, "revision/TableS13_recovery_simulation_summary.csv"))
cat("Method のコード名:", paste(unique(raw13$Method), collapse = ", "), "\n")   # 対応を確認
raw13$method <- factor(METHOD_LAB(raw13$Method), levels = METHOD_LEVELS)
sum13$method <- factor(METHOD_LAB(sum13$Method), levels = METHOD_LEVELS)
stopifnot(length(unique(raw13$method)) == length(unique(raw13$Method)))     # 3 手法が別ラベルに落ちたか

# 4a. 無構造データ（k_true = 1）で選ばれた k（手法別、20 反復の %）
nk <- raw13 %>% dplyr::filter(k_true == 1) %>% dplyr::mutate(k = k_selected) %>% dplyr::count(method, k) %>%
  dplyr::group_by(method) %>% dplyr::mutate(prop = 100 * n / sum(n)) %>% dplyr::ungroup() %>%
  tidyr::complete(method, k = 2:6, fill = list(n = 0, prop = 0))
pA <- ggplot(nk, aes(factor(k), prop, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = METHOD_COLS) + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "Null data: k selected", x = "Selected k", y = "Datasets, %", tag = "A") +
  theme_journal() + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1),
                          legend.direction = "vertical", legend.key.size = unit(3, "mm"))

# 4b. k = 2 を強制したときの null silhouette 分布（UMAP）vs 実データ 0.457
sil <- raw13 %>% dplyr::filter(k_true == 1, method == "UMAP + Ward's")
if (!exists("OBS_SIL")) OBS_SIL <- 0.457
cat(sprintf("null Sil_k2 (UMAP): mean %.3f, 95th pct %.3f, max %.3f (n = %d)\n",
            mean(sil$Sil_k2), quantile(sil$Sil_k2, .95), max(sil$Sil_k2), nrow(sil)))   # SUMMARY v2 §13 と照合
pB <- ggplot(sil, aes(Sil_k2)) +
  geom_histogram(bins = 10, fill = JAMA[["grey"]], colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = OBS_SIL, colour = JAMA[["red"]], linewidth = 0.6) +
  annotate("text", x = OBS_SIL, y = Inf, label = sprintf("Observed %.3f", OBS_SIL), hjust = 1.05, vjust = 1.5,
           size = 2.6, colour = JAMA[["red"]], family = "Helvetica") +
  labs(title = "Null data: silhouette with k = 2 forced", x = "Silhouette coefficient (UMAP + Ward's)", y = "Datasets, n", tag = "B") +
  theme_journal()

# 4c. 真の構造の回収：ARI vs 重心間距離（UMAP + Ward's；真 k = 2, 3, 4）
SEP <- tryCatch(read.csv(file.path(OUTPUT_DIR, "revision/TableS13_real_centroid_separation.csv"))$Separation_PC1SD[1], error = function(e) 1.46)
rcU <- sum13 %>% dplyr::filter(method == "UMAP + Ward's", k_true %in% 2:4) %>% dplyr::mutate(k_true = factor(k_true))
pC <- ggplot(rcU, aes(delta, ARI_at_true_k, colour = k_true, group = k_true)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
  geom_vline(xintercept = SEP, linetype = "22", colour = "grey50", linewidth = 0.4) +
  annotate("text", x = SEP, y = 0.02, label = sprintf("Observed\n%.2f SD", SEP), hjust = -0.1, vjust = 0, size = 2.4,
           colour = "grey30", family = "Helvetica", lineheight = 0.9) +
  scale_colour_manual(values = c(`2` = JAMA[["slate"]], `3` = JAMA[["cerulean"]], `4` = JAMA[["gold"]]),
                      labels = c("True k = 2", "True k = 3", "True k = 4")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 6)) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Recovery of known structure", x = "Between-centroid separation (SD of PC1)",
       y = "ARI vs true labels at the true k (UMAP + Ward's)", tag = "C") +
  theme_journal() + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.05), legend.justification = c(1, 0),
                          legend.direction = "vertical", legend.key.size = unit(3, "mm"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p13 <- (pA | pB | pC) + plot_layout(widths = c(1, 1, 1))
  save_fig(p13, "FigS13_recovery_simulation_v2", 240, 80)
} else {
  save_fig(pA, "FigS13A_null_k_v2", 85, 70); save_fig(pB, "FigS13B_null_silhouette_v2", 85, 70); save_fig(pC, "FigS13C_recovery_v2", 85, 70)
  cat("patchwork がないため 3 枚に分割保存しました（install.packages('patchwork') で 1 枚に）。\n")
}
cat("\n完了。output/revision/figs を確認し、PNG を添付してください。\n")
