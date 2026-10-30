# ==============================================================================
# 36_main_figures_v4.R（dplyr:: 明示版）— Fig 1–3 の再作成（訂正後データ・k=2、journal-figures 規範、R1-m1）
#   v4（2026-10-24）: v3 からの変更は「色」だけ。構成・注記・閾値・▲ の処理は v3 と同一。
#     * クラスタ 2 色: C1 = JAMA slate (#374E55), C2 = JAMA gold (#DF8F44)（v3 の warm-grey は対比不足）
#     * Fig 1 B–D の連続色: viridis "mako" を反転（暗い = 高値）。凡例 "darker points indicate higher values" と整合、色覚バリアフリー
#     * 点をやや大きく・濃く（Fig 1: size 1.3, alpha 0.85; Fig 2: size 1.4, alpha 0.7）
#     * Fig 3 の箱の塗りも 2 色の淡色（slate 淡 / gold 淡）に変更
#   Fig 1: UMAP 埋め込み（A クラスタ、B AHI、C ISI、D DBAS）
#   Fig 2: ISI vs AHI 散布図（クラスタ別）、閾値線と COMISA 領域
#   Fig 3: 9 変数のクラスタ別箱ひげ（AHI, ISI, BMI, DBAS, FIRST, SDS, TST, SE, WASO）
#   出力: output/revision/figs/Fig1_UMAP_v3 / Fig2_ISI_vs_AHI_v3 / Fig3_Variables_by_cluster_v3 (.pdf, 600dpi .png, .tiff)
# ------------------------------------------------------------------------------
# 入力（個人単位。Mac 上でのみ使用、共有しない）
#   EMB_CSV : 06 の出力。列 UMAP1, UMAP2。行順 = comisa.csv
#   comisa.csv（訂正済み）, output/Cluster_Assignments.csv（ID = 行番号）
#   output/revision/TableS12_pairwise_Tukey_CohenD.csv（Fig 3 の d 注記）
# ==============================================================================

if (!exists("DATA_PATH") || !exists("OUTPUT_DIR")) stop("Run source(\"01_setup.R\") first (defines DATA_PATH and OUTPUT_DIR).")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })
FIGDIR <- file.path(OUTPUT_DIR, "revision/figs"); dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
EMB_CSV <- file.path(OUTPUT_DIR, "UMAP_Embedding_Final.csv")      # 無ければ Fig 1 をスキップ
SUFFIX  <- "_v3"                                   # 出力ファイル名の版（v2 は残す）

JAMA <- c(slate = "#374E55", gold = "#DF8F44", cerulean = "#00A1D5", red = "#B24745",
          sage = "#79AF97", violet = "#6A6599", grey = "#80796B")
CLCOL  <- c(`Cluster 1` = JAMA[["slate"]], `Cluster 2` = JAMA[["gold"]])      # v4: gold
CLFILL <- c(`Cluster 1` = "#AAB8BE",        `Cluster 2` = "#F3D1AE")           # 淡色（Fig 3 の箱）
theme_journal <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          axis.text = element_text(size = base_size - 1, colour = "black"), axis.title = element_text(size = base_size),
          plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
          plot.tag = element_text(size = 12, face = "bold"),
          strip.background = element_blank(), strip.text = element_text(size = base_size, hjust = 0),
          legend.position = "top", legend.title = element_blank(), legend.key = element_blank(),
          legend.text = element_text(size = base_size - 1),
          plot.background = element_rect(fill = "white", colour = NA))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")),  p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(FIGDIR, paste0(name, ".png")),  p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(FIGDIR, paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white", compression = "lzw")
  cat("wrote", name, "\n")
}

# ---- データ ------------------------------------------------------------------
raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
cl  <- read.csv(file.path(OUTPUT_DIR, "Cluster_Assignments.csv")); clcol <- intersect(c("Cluster","cluster","Cluster_k2"), names(cl))[1]
HAVE_EMB <- file.exists(EMB_CSV)
stopifnot(nrow(raw) == 599, nrow(cl) == 599)
d <- raw %>% dplyr::mutate(Cluster = factor(cl[[clcol]], levels = 1:2, labels = c("Cluster 1", "Cluster 2")))
if (HAVE_EMB) { emb <- read.csv(EMB_CSV); stopifnot(nrow(emb) == 599); d$UMAP1 <- emb$UMAP1; d$UMAP2 <- emb$UMAP2 }
stopifnot(table(d$Cluster)[["Cluster 1"]] == 247)
n_lab <- as.data.frame(table(Cluster = d$Cluster), responseName = "n"); n_lab$lab <- paste0(n_lab$Cluster, " (n = ", n_lab$n, ")")
CLLAB <- setNames(n_lab$lab, n_lab$Cluster)

# ==============================================================================
# Fig 1 — UMAP embedding
# ==============================================================================
if (HAVE_EMB) {
base1 <- function(p) p + labs(x = "UMAP 1", y = "UMAP 2") + theme_journal() +
  theme(legend.position = "right", legend.title = element_text(size = 8), legend.key.height = unit(5, "mm"),
        legend.key.width = unit(3, "mm"), axis.text = element_blank(), axis.ticks = element_blank())
pt <- list(size = 1.3, alpha = 0.85, stroke = 0)
p1A <- base1(ggplot(d, aes(UMAP1, UMAP2, colour = Cluster)) + geom_point(size = pt$size, alpha = pt$alpha, stroke = 0) +
               scale_colour_manual(values = CLCOL, labels = CLLAB) +
               guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))) + labs(tag = "A") +
  theme(legend.position = "top", legend.title = element_blank(), legend.justification = "left")
# v4: 連続色 = viridis "mako" 反転（淡いミント = 低値 → 濃紺 = 高値）。凡例 "darker points indicate higher values" と整合
grad <- function(p, name) p + scale_colour_viridis_c(option = "mako", direction = -1, begin = 0.08, end = 0.92, name = name, na.value = "white")
p1B <- base1(grad(ggplot(d, aes(UMAP1, UMAP2, colour = AHI)) + geom_point(size = pt$size, alpha = pt$alpha, stroke = 0), "AHI\n(events/h)")) + labs(tag = "B")
p1C <- base1(grad(ggplot(d, aes(UMAP1, UMAP2, colour = ISI)) + geom_point(size = pt$size, alpha = pt$alpha, stroke = 0), "ISI\nscore")) + labs(tag = "C")
p1D <- base1(grad(ggplot(d %>% dplyr::arrange(!is.na(DBAS)), aes(UMAP1, UMAP2, colour = DBAS)) +
                    geom_point(data = dplyr::filter(d, is.na(DBAS)), colour = "grey60", shape = 1, size = pt$size, stroke = 0.3) +
                    geom_point(data = dplyr::filter(d, !is.na(DBAS)), size = pt$size, alpha = pt$alpha, stroke = 0), "DBAS\nscore")) + labs(tag = "D")
fig1 <- (p1A | p1B) / (p1C | p1D)
save_fig(fig1, paste0("Fig1_UMAP", SUFFIX), 180, 150)
} else cat("[skip] Fig 1: ", EMB_CSV, " がありません（Fig 2・3 は続行）\n")

# ==============================================================================
# Fig 2 — ISI vs AHI（v3 と同一構成。色と点のみ変更）
# ==============================================================================
XMAX <- 120
d2 <- d %>% dplyr::mutate(over = AHI > XMAX, AHI_plot = pmin(AHI, XMAX))
n_over <- sum(d2$over)
p2 <- ggplot(d2, aes(AHI_plot, ISI, colour = Cluster)) +
  annotate("rect", xmin = 5, xmax = Inf, ymin = 15, ymax = Inf, fill = "grey95", colour = NA) +
  geom_vline(xintercept = c(5, 15), linetype = c("22", "solid"), colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 15, linetype = "22", colour = "grey40", linewidth = 0.4) +
  geom_point(data = dplyr::filter(d2, !over), size = 1.4, alpha = 0.7, stroke = 0, position = position_jitter(width = 0, height = 0.15, seed = 1)) +
  geom_point(data = dplyr::filter(d2, over), shape = 17, size = 1.8, show.legend = FALSE) +
  annotate("text", x = 5,  y = 29.5, label = "AHI 5",  hjust = -0.1, size = 2.6, colour = "grey30", family = "Helvetica") +
  annotate("text", x = 15, y = 29.5, label = "AHI 15", hjust = -0.1, size = 2.6, colour = "grey30", family = "Helvetica") +
  annotate("text", x = XMAX, y = 15, label = "ISI 15", hjust = 1, vjust = -0.4, size = 2.6, colour = "grey30", family = "Helvetica") +
  annotate("text", x = XMAX, y = 28.5, hjust = 1, size = 2.6, colour = "grey30", family = "Helvetica",
           label = "Symptom-based COMISA\n(ISI ≥ 15 and AHI ≥ 5)") +
  scale_colour_manual(values = CLCOL, labels = CLLAB) +
  scale_x_continuous(limits = c(0, XMAX), breaks = seq(0, XMAX, 20), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(limits = c(-0.5, 30), breaks = seq(0, 28, 7)) +
  labs(x = "Apnea–hypopnea index (events/h)", y = "Insomnia Severity Index",
       caption = if (n_over > 0) paste0("Filled triangle: ", n_over, ifelse(n_over == 1, " patient", " patients"), " with AHI > ", XMAX, " shown at ", XMAX) else NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  theme_journal() + theme(plot.caption = element_text(size = 7, hjust = 0, colour = "grey30"))
save_fig(p2, paste0("Fig2_ISI_vs_AHI", SUFFIX), 120, 95)

# ==============================================================================
# Fig 3 — 9 variables by cluster（◇ = mean、d を注記、p は載せない。塗りのみ変更）
# ==============================================================================
VARS <- c(AHI = "AHI (events/h)", ISI = "Insomnia Severity Index",
          BMI = "Body mass index (kg/m²)", DBAS = "DBAS score", FIRST = "FIRST score", SDS = "SDS score",
          TST = "Total sleep time (min)", SE = "Sleep efficiency (proportion)", WASO = "Wake after sleep onset (min)")
VARS_LAB <- setNames(paste0(LETTERS[seq_along(VARS)], "   ", VARS), names(VARS))
es <- read.csv(file.path(OUTPUT_DIR, "revision/TableS12_pairwise_Tukey_CohenD.csv")) %>%
  dplyr::transmute(Variable, lab = sprintf("d = %s%.2f", ifelse(Cohen_d < 0, "−", ""), abs(Cohen_d)))
stopifnot(all(names(VARS) %in% names(d)))
long3 <- d %>% dplyr::select(Cluster, dplyr::all_of(names(VARS))) %>% dplyr::mutate(dplyr::across(-Cluster, ~ as.numeric(.x))) %>%
  tidyr::pivot_longer(-Cluster, names_to = "Variable") %>%
  dplyr::filter(!is.na(value)) %>% dplyr::mutate(Variable = factor(Variable, levels = names(VARS)))
es <- es %>% dplyr::filter(Variable %in% names(VARS)) %>% dplyr::mutate(Variable = factor(Variable, levels = names(VARS)))
p3 <- ggplot(long3, aes(Cluster, value, fill = Cluster)) +
  geom_boxplot(width = 0.55, outlier.size = 0.5, outlier.alpha = 0.4, outlier.colour = "grey40", linewidth = 0.35, colour = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white", colour = "black", stroke = 0.4) +
  geom_text(data = es, aes(x = 2.5, y = Inf, label = lab), inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 2.5, family = "Helvetica") +
  scale_fill_manual(values = CLFILL, guide = "none") +
  scale_x_discrete(labels = c("C1", "C2")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3, labeller = as_labeller(VARS_LAB)) +
  labs(x = NULL, y = NULL,
       caption = paste0("C1, Cluster 1 (n = 247; low AHI, high psychological burden); C2, Cluster 2 (n = 352; high AHI). AHI, apnea–hypopnea index.\n",
                        "◇ mean; box, median and IQR; whiskers, 1.5 × IQR; points, values beyond the whiskers. d, Cohen's d (C1 − C2). Observed values; n varies with missingness.")) +
  theme_journal() + theme(plot.caption = element_text(size = 7, hjust = 0, colour = "grey30"),
                          strip.text = element_text(size = 8, hjust = 0, face = "bold"), panel.spacing = unit(4, "mm"))
save_fig(p3, paste0("Fig3_Variables_by_cluster", SUFFIX), 180, 180)
cat("完了。output/revision/figs の *_v3.png を確認してください（v2 は上書きしません）。\n")
