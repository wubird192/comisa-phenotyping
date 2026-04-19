# ============================================================
# COMISA Phenotyping Study
# 06_figures.R — Publication-ready figures for Lancet Digital Health
#
#   Figure 1: UMAP phenotype projections (4-panel)
#   Figure 2: ISI vs AHI scatter plot (with AHI = 15 primary line)
#   Figure 3: Clinical characteristics boxplots (9-panel)
#             Panel H y-axis: "Sleep efficiency (proportion)"
#   Figure S2: Clinical heatmap (supplementary)
#
# Note: Figure 4 (ROC curves) is generated in 04_validation.R
# All figures saved as PNG (300 DPI) and TIFF (300 DPI, LZW)
# Dimensions follow Lancet Digital Health artwork guidelines:
#   single column: 89 mm; full width: 180 mm
# ============================================================

cat("\n============================================================\n")
cat("  06: Publication Figures\n")
cat("============================================================\n\n")

n_clusters <- length(unique(comisa$Cluster[!is.na(comisa$Cluster)]))

# Shared ggplot theme
lancet_theme <- theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray50", fill = NA, linewidth = 0.5)
  )

# save_fig: saves PNG (preview) + TIFF (Lancet submission, 300 DPI).
# macOS quartz TIFF does not support LZW; cairo is used on Linux/Windows.
save_fig <- function(plot_obj, name, w, h) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), plot_obj,
         width = w, height = h, units = "mm", dpi = 300)

  tiff_path <- file.path(OUTPUT_DIR, paste0(name, ".tiff"))
  on_mac    <- Sys.info()[["sysname"]] == "Darwin"
  if (on_mac) {
    # quartz: no compression argument
    suppressWarnings(
      ggsave(tiff_path, plot_obj,
             width = w, height = h, units = "mm", dpi = 300,
             device = "tiff")
    )
  } else {
    ggsave(tiff_path, plot_obj,
           width = w, height = h, units = "mm", dpi = 300,
           compression = "lzw", device = "tiff")
  }
  cat(sprintf("  \u2713 %s.png / .tiff saved\n", name))
}

# ============================================================
# FIGURE 1 — UMAP phenotype projections (4-panel)
#   A: cluster labels
#   B: AHI continuous colour
#   C: ISI continuous colour
#   D: DBAS continuous colour
# ============================================================

cat("=== Figure 1: UMAP projections ===\n")

umap_df <- data.frame(
  UMAP1   = final_result$umap_result[, 1],
  UMAP2   = final_result$umap_result[, 2],
  Cluster = factor(comisa$Cluster),
  AHI     = comisa$AHI,
  ISI     = comisa$ISI,
  DBAS    = comisa$DBAS
)

common_umap_theme <- lancet_theme +
  theme(legend.position = "right")

# Panel A — cluster labels
pA <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Cluster)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = CLUSTER_COLORS[seq_len(n_clusters)],
                     name = "Cluster",
                     labels = paste0("C", seq_len(n_clusters))) +
  labs(title = "A", x = "UMAP 1", y = "UMAP 2") +
  common_umap_theme

# Panel B — AHI
pB <- ggplot(umap_df, aes(UMAP1, UMAP2, color = AHI)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_viridis_c(name = "AHI\n(events/h)", option = "plasma") +
  labs(title = "B", x = "UMAP 1", y = "UMAP 2") +
  common_umap_theme

# Panel C — ISI
pC <- ggplot(umap_df, aes(UMAP1, UMAP2, color = ISI)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_viridis_c(name = "ISI\nscore", option = "viridis") +
  labs(title = "C", x = "UMAP 1", y = "UMAP 2") +
  common_umap_theme

# Panel D — DBAS
pD <- ggplot(umap_df, aes(UMAP1, UMAP2, color = DBAS)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_viridis_c(name = "DBAS\nscore", option = "inferno") +
  labs(title = "D", x = "UMAP 1", y = "UMAP 2") +
  common_umap_theme

fig1 <- plot_grid(pA, pB, pC, pD, ncol = 2, nrow = 2, align = "hv")
save_fig(fig1, "Figure1_COMISA_Phenotypes_LDH", 180, 160)
cat("\n")

# ============================================================
# FIGURE 2 — ISI vs AHI scatter plot
#   Reference lines: ISI = 15 (dashed red), AHI = 5 (dashed blue),
#                    AHI = 15 (solid dark blue — primary threshold)
#   COMISA zone shading: ISI ≥ 15 & AHI ≥ 15
# ============================================================

cat("=== Figure 2: ISI vs AHI scatter ===\n")

scatter_df <- data.frame(
  AHI     = comisa$AHI,
  ISI     = comisa$ISI,
  Cluster = factor(comisa$Cluster)
)
scatter_df    <- scatter_df[!is.na(scatter_df$Cluster) &
                             !is.na(scatter_df$AHI) &
                             !is.na(scatter_df$ISI), ]
cl_sizes      <- table(scatter_df$Cluster)
cl_labels     <- paste0(seq_len(n_clusters),
                        " (n=", cl_sizes[seq_len(n_clusters)], ")")

p2 <- ggplot(scatter_df, aes(AHI, ISI, color = Cluster)) +

  # COMISA zone shading
  annotate("rect", xmin = 15, xmax = Inf, ymin = 15, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.3) +
  annotate("text", x = 85, y = 27,
           label = "COMISA zone\n(ISI \u2265 15 & AHI \u2265 15)",
           color = "gray40", size = 3, fontface = "italic") +

  # Data points
  geom_point(size = 2.5, alpha = 0.6) +

  # Threshold lines
  geom_hline(yintercept = 15, linetype = "dashed",
             color = "#D32F2F", linewidth = 0.8) +
  geom_vline(xintercept = 5,  linetype = "dashed",
             color = "#1976D2", linewidth = 0.6) +
  geom_vline(xintercept = 15, linetype = "solid",
             color = "#0D47A1", linewidth = 0.9) +

  # Threshold annotations
  annotate("text", x = 1.5,  y = 16.2,
           label = "ISI = 15", color = "#D32F2F",
           size = 3.2, hjust = 0, fontface = "bold") +
  annotate("text", x = 0.5,  y = 2.5,
           label = "AHI = 5\n(sensitivity)",
           color = "#1976D2", size = 2.8, hjust = 0) +
  annotate("text", x = 16.5, y = 2.5,
           label = "AHI = 15\n(primary)",
           color = "#0D47A1", size = 2.8, hjust = 0, fontface = "bold") +

  scale_color_manual(values = CLUSTER_COLORS[seq_len(n_clusters)],
                     name = "Cluster", labels = cl_labels) +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 150, 25)) +
  scale_y_continuous(limits = c(0, 28), breaks = seq(0, 28, 7)) +
  labs(title = "ISI vs AHI by Cluster",
       x = "AHI (events/h)", y = "ISI score") +
  lancet_theme +
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.title = element_text(face = "bold", size = 10))

save_fig(p2, "Figure2_ISI_vs_AHI_Scatter_LDH", 140, 140)
cat("\n")

# ============================================================
# FIGURE 3 — Clinical characteristics boxplots (9-panel)
#   Variables: AHI, ISI, BMI, DBAS, FIRST, SDS, TST, SE, WASO
#   Panel H y-axis: "Sleep efficiency (proportion)"  ← FIXED
# ============================================================

cat("=== Figure 3: Clinical characteristics (9-panel) ===\n")

key_vars   <- c("AHI", "ISI", "BMI", "DBAS", "FIRST", "SDS",
                "TST", "SE",  "WASO")
var_labels <- c(
  AHI   = "AHI (events/h)",
  ISI   = "ISI score",
  BMI   = "BMI (kg/m\u00b2)",
  DBAS  = "DBAS score",
  FIRST = "FIRST score",
  SDS   = "SDS score",
  TST   = "TST (min)",
  SE    = "Sleep efficiency (proportion)",   # NOTE: values are 0–1
  WASO  = "WASO (min)"
)

# ANOVA + Tukey HSD
anova_res <- lapply(setNames(key_vars, key_vars), function(var) {
  df  <- data.frame(y = comisa[[var]], cl = factor(comisa$Cluster))
  df  <- df[complete.cases(df), ]
  ao  <- aov(y ~ cl, data = df)
  tuk <- TukeyHSD(ao)
  list(
    anova_p = summary(ao)[[1]]$`Pr(>F)`[1],
    f_stat  = summary(ao)[[1]]$`F value`[1],
    tukey   = tuk$cl,
    means   = tapply(df$y, df$cl, mean,  na.rm = TRUE),
    sds     = tapply(df$y, df$cl, sd,    na.rm = TRUE)
  )
})

# Save ANOVA summary
anova_df <- do.call(rbind, lapply(names(anova_res), function(v) {
  r <- anova_res[[v]]
  data.frame(
    Variable    = v,
    F_statistic = round(r$f_stat, 2),
    P_value     = signif(r$anova_p, 3),
    C1_mean_sd  = sprintf("%.1f ± %.1f", r$means["1"], r$sds["1"]),
    C2_mean_sd  = sprintf("%.1f ± %.1f", r$means["2"], r$sds["2"]),
    C3_mean_sd  = sprintf("%.1f ± %.1f", r$means["3"], r$sds["3"])
  )
}))
write.csv(anova_df, file.path(OUTPUT_DIR, "ANOVA_Results_Clinical_Variables.csv"),
          row.names = FALSE)

# Build plots
panel_list <- lapply(seq_along(key_vars), function(i) {
  var  <- key_vars[i]
  pd   <- data.frame(value   = comisa[[var]],
                     cluster = factor(comisa$Cluster))
  pd   <- pd[complete.cases(pd), ]
  mn   <- data.frame(cluster = factor(seq_len(n_clusters)),
                     mean    = anova_res[[var]]$means[seq_len(n_clusters)])
  plab <- if (anova_res[[var]]$anova_p < 0.001) "p<0·001"
          else sprintf("p=%.3f", anova_res[[var]]$anova_p)

  ggplot(pd, aes(cluster, value, fill = cluster)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.alpha = 0.5) +
    geom_point(data = mn, aes(cluster, mean),
               shape = 23, size = 3, fill = "white", color = "black") +
    scale_fill_manual(values = CLUSTER_COLORS[seq_len(n_clusters)]) +
    annotate("text", x = Inf, y = Inf, label = plab,
             hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") +
    labs(title = LETTERS[i], x = "Cluster", y = var_labels[var]) +
    lancet_theme +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9))
})

fig3 <- plot_grid(plotlist = panel_list, ncol = 3, nrow = 3, align = "hv")
save_fig(fig3, "Figure3_Clinical_Variables_LDH", 180, 240)
cat("\n")

# ============================================================
# SUPPLEMENTARY FIGURE S2 — Clinical heatmap
# ============================================================

cat("=== Supplementary Figure S2: Clinical heatmap ===\n")

# pheatmap annotation_colors causes "subscript out of bounds" on some
# R/pheatmap versions; replaced with ggplot2 geom_tile implementation.

heatmap_vars <- c("AHI", "ISI", "DBAS", "FIRST", "HAS",
                  "SDS", "JESS", "TST", "SE", "WASO",
                  "Arousal_Index", "ODI3", "BMI")
heatmap_vars <- heatmap_vars[heatmap_vars %in% names(comisa)]

valid_idx <- !is.na(comisa$Cluster) &
             complete.cases(comisa[, heatmap_vars])
hm_raw    <- comisa[valid_idx, heatmap_vars]
hm_cl     <- factor(comisa$Cluster[valid_idx], levels = seq_len(n_clusters))

hm_scaled <- as.data.frame(scale(hm_raw))
ord        <- order(hm_cl)
hm_sorted  <- hm_scaled[ord, ]
cl_sorted  <- hm_cl[ord]
n_pat      <- nrow(hm_sorted)
n_var      <- ncol(hm_sorted)

var_labels_hm <- c(
  AHI = "AHI", ISI = "ISI", DBAS = "DBAS", FIRST = "FIRST",
  HAS = "HAS", SDS = "SDS", JESS = "ESS", TST = "TST",
  SE = "Sleep efficiency", WASO = "WASO",
  Arousal_Index = "Arousal index", ODI3 = "ODI3", BMI = "BMI"
)
var_names <- names(hm_sorted)
var_disp  <- ifelse(var_names %in% names(var_labels_hm),
                    var_labels_hm[var_names], var_names)

s2_long <- data.frame(
  patient  = rep(seq_len(n_pat), times = n_var),
  variable = rep(factor(var_disp, levels = rev(var_disp)), each = n_pat),
  z_score  = unlist(hm_sorted, use.names = FALSE)
)
s2_long$z_score <- pmax(pmin(s2_long$z_score, 3), -3)

cl_pal_s2 <- CLUSTER_COLORS[seq_len(n_clusters)]
names(cl_pal_s2) <- as.character(seq_len(n_clusters))

strip_s2 <- data.frame(patient = seq_len(n_pat), Cluster = cl_sorted)

p_strip_s2 <- ggplot(strip_s2, aes(x = patient, y = 1, fill = Cluster)) +
  geom_tile() +
  scale_fill_manual(values = cl_pal_s2, name = "Cluster") +
  theme_void() +
  theme(legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.text     = element_text(size = 8),
        legend.title    = element_text(size = 8, face = "bold"))

p_hm_s2 <- ggplot(s2_long, aes(x = patient, y = variable, fill = z_score)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits  = c(-3, 3), name = "z-score",
    breaks  = c(-3, -1.5, 0, 1.5, 3),
    labels  = c("-3", "-1.5", "0", "1.5", "3")
  ) +
  labs(title = "Clinical Characteristics by Cluster", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 8),
        panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold", size = 10))

fig_s2 <- plot_grid(p_strip_s2, p_hm_s2,
                    ncol = 1, rel_heights = c(0.04, 0.96),
                    align = "v", axis = "lr")

ggsave(file.path(OUTPUT_DIR, "FigureS2_Clinical_Heatmap_LDH.png"),
       fig_s2, width = 180, height = 220, units = "mm", dpi = 300)
cat("  \u2713 FigureS2_Clinical_Heatmap_LDH.png saved\n\n")

cat("============================================================\n")
cat("  Figures complete. Proceed to 07_tables.R\n")
cat("============================================================\n\n")

# ============================================================
# SUPPLEMENTARY FIGURE S1 — ROC Curve Comparison: AHI >=15 vs AHI >=5
#
#   Reference population: ALL 599 participants (not limited to
#   those with expert diagnosis). Allows direct threshold comparison.
#   No diagonal reference line (per manuscript style).
#
#   NOTE (v4, April 2026): In the submitted manuscript (v4+), this figure
#   has been converted to Table S2 (numeric AUC/sensitivity/specificity
#   table) following co-author feedback. This code block still generates
#   the ROC figure and FigureS1_AUC_Comparison.csv for internal reference,
#   but the PNG/TIFF outputs are NOT used in the submitted paper.
#   The underlying data (AUC values, DeLong p) feed directly into
#   TableS2_AHI_Sensitivity_Analysis.csv generated by 07_tables.R.
# ============================================================

cat("=== Supplementary Figure S1: ROC curves AHI >=15 vs AHI >=5 ===\n")

comisa_s1    <- comisa[!is.na(comisa$Cluster), ]
label_ahi15  <- as.integer(comisa_s1$Diagnosis == "\u4e0d\u7720\u75c7" &
                            !is.na(comisa_s1$AHI) & comisa_s1$AHI >= 15)
label_ahi5   <- as.integer(comisa_s1$Diagnosis == "\u4e0d\u7720\u75c7" &
                            !is.na(comisa_s1$AHI) & comisa_s1$AHI >= 5)
pred_c1_s1   <- as.integer(comisa_s1$Cluster == 1)

roc_s1_15 <- roc(label_ahi15, pred_c1_s1, levels = c(0,1), direction = "<", quiet = TRUE)
roc_s1_5  <- roc(label_ahi5,  pred_c1_s1, levels = c(0,1), direction = "<", quiet = TRUE)
auc_15    <- as.numeric(auc(roc_s1_15))
auc_5     <- as.numeric(auc(roc_s1_5))
ci_15     <- as.numeric(ci.auc(roc_s1_15))
ci_5      <- as.numeric(ci.auc(roc_s1_5))
delong    <- roc.test(roc_s1_15, roc_s1_5, method = "delong")

cat(sprintf("  AHI>=15: AUC=%.3f (95%% CI %.3f-%.3f)\n", auc_15, ci_15[1], ci_15[3]))
cat(sprintf("  AHI>=5:  AUC=%.3f (95%% CI %.3f-%.3f)\n", auc_5,  ci_5[1],  ci_5[3]))
cat(sprintf("  DeLong: p=%.3f\n\n", delong$p.value))

draw_fig_s1 <- function(cex_leg = 0.72) {
  plot.new()
  plot.window(xlim = c(0,1), ylim = c(0,1), asp = 1)
  axis(1, at = seq(0,1,0.2)); axis(2, at = seq(0,1,0.2), las = 1); box()
  title(main = "Supplementary Figure S1",
        xlab = "1 - Specificity", ylab = "Sensitivity",
        cex.main = 1.1, font.main = 2)
  # NO diagonal line
  lines(1 - rev(roc_s1_15$specificities), rev(roc_s1_15$sensitivities),
        col = "#E69F00", lwd = 2.5)
  lines(1 - rev(roc_s1_5$specificities),  rev(roc_s1_5$sensitivities),
        col = "#0072B2", lwd = 2.5)
  legend("bottomright",
    legend = c(
      sprintf("AHI >=15 events/h  AUC=%.2f (95%% CI %.2f-%.2f)", auc_15, ci_15[1], ci_15[3]),
      sprintf("AHI >=5  events/h  AUC=%.2f (95%% CI %.2f-%.2f)", auc_5,  ci_5[1],  ci_5[3]),
      sprintf("DeLong test: p=%.3f", delong$p.value)),
    col = c("#E69F00","#0072B2", NA), lty = c(1,1,NA), lwd = c(2.5,2.5,NA),
    cex = cex_leg, bty = "n", seg.len = 1.5)
}

png(file.path(OUTPUT_DIR, "FigureS1_ROC_AHI_Threshold_Comparison.png"),
    width = 140, height = 140, units = "mm", res = 300)
par(pty = "s", mar = c(5,4,4,2), family = "Arial", ps = 9)
draw_fig_s1(0.72); dev.off()

on_mac_s1    <- Sys.info()[["sysname"]] == "Darwin"
tiff_path_s1 <- file.path(OUTPUT_DIR, "FigureS1_ROC_AHI_Threshold_Comparison.tiff")
if (on_mac_s1) {
  tiff(tiff_path_s1, width = 140, height = 140, units = "mm", res = 300)
} else {
  tiff(tiff_path_s1, width = 140, height = 140, units = "mm", res = 300,
       compression = "lzw", type = "cairo")
}
par(pty = "s", mar = c(5,4,4,2), family = "Arial", ps = 9)
draw_fig_s1(0.68); dev.off()

cat("  \u2713 FigureS1_ROC_AHI_Threshold_Comparison.png / .tiff saved\n\n")

write.csv(
  data.frame(
    Threshold = c("AHI >= 15","AHI >= 5"),
    N_COMISA  = c(sum(label_ahi15), sum(label_ahi5)),
    AUC       = round(c(auc_15, auc_5), 3),
    CI_Lower  = round(c(ci_15[1], ci_5[1]), 3),
    CI_Upper  = round(c(ci_15[3], ci_5[3]), 3),
    DeLong_p  = c(round(delong$p.value, 3), NA)),
  file.path(OUTPUT_DIR, "FigureS1_AUC_Comparison.csv"), row.names = FALSE)
cat("  \u2713 FigureS1_AUC_Comparison.csv saved\n\n")
