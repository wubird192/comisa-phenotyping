# ==============================================================================
# 34_stream_logistic_CI.R
# Cluster 1 membership のロジスティック回帰（Appendix D4b）に 95% CI を付ける
#   出力: output/revision/TableS9b_stream_logistic_CI.csv
#         列 = Model, Term, n, OR, CI_low, CI_high, p   （OR・CI は小数 3 桁）
# ------------------------------------------------------------------------------
# 方針
#  * 既存スクリプト（TableS9b_stream_logistic.csv を出したもの）と同じ 3 モデル
#      crude        : C1 ~ Stream
#      ISI+AHI      : C1 ~ Stream + ISI + AHI
#      ISI+AHI+DBAS : C1 ~ Stream + ISI + AHI + DBAS
#    ISI / AHI / DBAS は 1 SD 単位に標準化（scale）。Stream の参照は IS。
#  * 95% CI は profile likelihood（confint.glm）。Wald CI ではない。
#  * ⟨要確認⟩ 既存スクリプトが「観測値（DBAS 欠測 53 を listwise 除外）」を使ったか
#    「補完データセット 4」を使ったかで n が変わる。既定は観測値。下の USE_IMPUTED で切替。
#    実行後に表示される n と OR が既存出力（StreamSAS: 0.04 → 0.22 (p=0.009) → 0.28 (p=0.062)）
#    と一致することを確認してから CSV を採用する。
# ==============================================================================

if (!exists("DATA_PATH") || !exists("OUTPUT_DIR")) stop("Run source(\"01_setup.R\") first (defines DATA_PATH and OUTPUT_DIR).")
set.seed(20261020)
suppressPackageStartupMessages({ library(dplyr) })

USE_IMPUTED <- FALSE                      # TRUE にすると補完 4 の ISI/AHI/DBAS を使う
IMP_RDS     <- file.path(OUTPUT_DIR, "mice_imputation.rds")  # ⟨要確認⟩ 02 で保存した mice オブジェクトのパス

# ---- 1. データ ---------------------------------------------------------------
raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
cl  <- read.csv(file.path(OUTPUT_DIR, "Cluster_Assignments.csv"))
stopifnot(nrow(raw) == 599, nrow(cl) == 599)
clcol <- intersect(c("Cluster", "cluster", "Cluster_k2"), names(cl))[1]
stopifnot(!is.na(clcol))

dat <- raw %>%
  mutate(row_id = row_number(),
         Stream = factor(ifelse(Diagnosis == "不眠症", "IS", "SAS"), levels = c("IS", "SAS")),
         C1     = as.integer(cl[[clcol]][row_id] == 1))

if (USE_IMPUTED) {
  imp4 <- mice::complete(readRDS(IMP_RDS), 4)        # 行順は comisa.csv と同じ前提
  dat$ISI  <- imp4$ISI;  dat$AHI <- imp4$AHI;  dat$DBAS <- imp4$DBAS
}
if (table(dat$Stream)[["IS"]] != 113 || sum(dat$C1) != 247) message("Note: stream/cluster counts differ from the manuscript (expected on synthetic data).")

dat <- dat %>% mutate(ISI_z = as.numeric(scale(ISI)),
                      AHI_z = as.numeric(scale(AHI)),
                      DBAS_z = as.numeric(scale(DBAS)))
# ⟨要確認⟩ 既存スクリプトが「モデルごとの完全ケースで scale」していた場合は
#           scale() を各 glm の直前（na.omit 後）に移す。

# ---- 2. モデルと CI ----------------------------------------------------------
fit_one <- function(formula, label) {
  d   <- na.omit(model.frame(formula, dat))
  fit <- glm(formula, data = d, family = binomial)
  ci  <- suppressMessages(exp(confint(fit)))           # profile likelihood
  co  <- summary(fit)$coefficients
  out <- data.frame(Model   = label,
                    Term    = rownames(co),
                    n       = nrow(d),
                    OR      = exp(co[, "Estimate"]),
                    CI_low  = ci[, 1],
                    CI_high = ci[, 2],
                    p       = co[, "Pr(>|z|)"],
                    row.names = NULL)
  out[out$Term != "(Intercept)", ]
}

res <- bind_rows(
  fit_one(C1 ~ Stream,                               "crude"),
  fit_one(C1 ~ Stream + ISI_z + AHI_z,               "ISI+AHI"),
  fit_one(C1 ~ Stream + ISI_z + AHI_z + DBAS_z,      "ISI+AHI+DBAS")
)
res$Term <- sub("_z$", "", res$Term)                  # ISI_z → ISI（既存 CSV と同じ Term 名）
res <- res %>% mutate(across(c(OR, CI_low, CI_high), ~ round(.x, 3)),
                      p = signif(p, 3))

print(res, digits = 3)
dir.create(file.path(OUTPUT_DIR, "revision"), showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(OUTPUT_DIR, "revision/TableS9b_stream_logistic_CI.csv"), row.names = FALSE)
cat("\nwrote output/revision/TableS9b_stream_logistic_CI.csv\n")
cat("確認: StreamSAS の OR が 0.04 / 0.22 (p=0.009) / 0.28 (p=0.062) に一致するか。\n")
