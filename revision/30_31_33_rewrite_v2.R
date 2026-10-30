# ==============================================================================
# 30_31_33_rewrite_v2.R — Sex by cluster / entry stream by cluster (30), unpaired DeLong test
#   AHI >= 5 vs >= 15 on the definitive-diagnosis sets (31), and diagnostic composition of each
#   cluster, Table 3A (33). One script, three sections. GitHub version: no local paths — uses
#   DATA_PATH and OUTPUT_DIR from 01_setup.R, so run source("01_setup.R") first.
#   Inputs : DATA_PATH (corrected analysis file), OUTPUT_DIR/Cluster_Assignments.csv (row order = analysis file, k = 2)
#   Outputs: OUTPUT_DIR/revision/TableS1_sex_by_cluster.csv, TableS1_stream_by_cluster.csv, TableS1_sex_within_SAS.csv,
#            TableS12_DeLong_def5_vs_def15_definitive.csv, Table3A_rowcounts.csv, Table3A_rowpct.csv
#   The final block checks the results against the values reported in the manuscript (real data only;
#   the checks are expected to report MISMATCH on the synthetic dataset).
# ==============================================================================
if (!exists("DATA_PATH") || !exists("OUTPUT_DIR")) stop("Run source(\"01_setup.R\") first (defines DATA_PATH and OUTPUT_DIR).")
suppressPackageStartupMessages({ library(dplyr); library(pROC) })
OUT <- file.path(OUTPUT_DIR, "revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
cl  <- read.csv(file.path(OUTPUT_DIR, "Cluster_Assignments.csv"))
clcol <- intersect(c("Cluster", "cluster", "Cluster_k2"), names(cl))[1]
stopifnot(nrow(raw) == 599, nrow(cl) == 599)
pick_col <- function(df, cands) { x <- intersect(cands, names(df)); if (!length(x)) stop("列が見つかりません: ", paste(cands, collapse = "/"), "\n  ある列: ", paste(names(df), collapse = ", ")); x[1] }
SEX  <- pick_col(raw, c("Sex", "sex", "Gender", "gender", "性別"))
DIAG <- pick_col(raw, c("Diagnosis", "diagnosis", "診断"))
FEMALE <- c("女", "女性", "F", "f", "Female", "female")       # the analysis file codes sex as 女 / 男
message("Levels of Sex: ", paste(names(table(raw[[SEX]])), table(raw[[SEX]]), sep = "=", collapse = ", "))

d <- raw %>% dplyr::mutate(
  Cluster  = factor(cl[[clcol]], levels = 1:2, labels = c("Cluster 1", "Cluster 2")),
  Female   = as.character(.data[[SEX]]) %in% FEMALE,
  InsDx    = grepl("不眠症", .data[[DIAG]]),                      # recorded insomnia diagnosis = insomnia stream (n = 113)
  Stream   = factor(ifelse(InsDx, "Insomnia stream", "Sleep apnea stream"), levels = c("Insomnia stream", "Sleep apnea stream")),
  COMISA5  = InsDx & AHI >= 5,  COMISA15 = InsDx & AHI >= 15,
  Group5   = dplyr::case_when(InsDx & AHI >= 5 ~ "COMISA", InsDx & AHI < 5 ~ "Insomnia disorder only",
                              !InsDx & AHI >= 5 ~ "OSA only", TRUE ~ "Neither"),
  Group15  = dplyr::case_when(InsDx & AHI >= 15 ~ "COMISA", InsDx & AHI < 15 ~ "Insomnia disorder only",
                              !InsDx & AHI >= 15 ~ "OSA only", TRUE ~ "Neither"))
if (sum(d$Cluster == "Cluster 1") != 247 || sum(d$InsDx) != 113) message("Note: cluster sizes differ from the manuscript (expected on synthetic data).")
cramer_v <- function(tab) { x <- suppressWarnings(chisq.test(tab)); c(chisq = unname(x$statistic), df = unname(x$parameter), p = x$p.value, V = sqrt(unname(x$statistic) / sum(tab) / (min(dim(tab)) - 1))) }
pct <- function(n, N) round(100 * n / N, 1)

# ---------------------------------------------------------------- 30_sex_stream
tab_sex <- table(d$Cluster, d$Female)
sex_by_cluster <- d %>% dplyr::group_by(Cluster) %>% dplyr::summarise(n = dplyr::n(), women = sum(Female), women_pct = pct(sum(Female), dplyr::n()), men = sum(!Female), .groups = "drop") %>%
  dplyr::mutate(chisq = cramer_v(tab_sex)[["chisq"]], df = cramer_v(tab_sex)[["df"]], p = cramer_v(tab_sex)[["p"]], cramer_V = cramer_v(tab_sex)[["V"]])
write.csv(sex_by_cluster, file.path(OUT, "TableS1_sex_by_cluster.csv"), row.names = FALSE)
tab_str <- table(d$Stream, d$Cluster)
stream_by_cluster <- d %>% dplyr::group_by(Stream) %>% dplyr::summarise(n = dplyr::n(), C1 = sum(Cluster == "Cluster 1"), C1_pct = pct(sum(Cluster == "Cluster 1"), dplyr::n()), C2 = sum(Cluster == "Cluster 2"), .groups = "drop") %>%
  dplyr::mutate(chisq = cramer_v(tab_str)[["chisq"]], p = cramer_v(tab_str)[["p"]], cramer_V = cramer_v(tab_str)[["V"]])
write.csv(stream_by_cluster, file.path(OUT, "TableS1_stream_by_cluster.csv"), row.names = FALSE)
sas <- d %>% dplyr::filter(Stream == "Sleep apnea stream")
sex_within_sas <- sas %>% dplyr::group_by(Cluster) %>% dplyr::summarise(n = dplyr::n(), women = sum(Female), women_pct = pct(sum(Female), dplyr::n()), .groups = "drop")
write.csv(sex_within_sas, file.path(OUT, "TableS1_sex_within_SAS.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 31_delong_definitive
# definitive set = patients classified as "Neither" excluded. AUC of the binary C1 indicator = balanced accuracy; the two sets differ in composition, hence unpaired
def5  <- d %>% dplyr::filter(Group5  != "Neither"); def15 <- d %>% dplyr::filter(Group15 != "Neither")
r5  <- pROC::roc(def5$COMISA5,   as.integer(def5$Cluster  == "Cluster 1"), quiet = TRUE, direction = "<")
r15 <- pROC::roc(def15$COMISA15, as.integer(def15$Cluster == "Cluster 1"), quiet = TRUE, direction = "<")
dl  <- pROC::roc.test(r5, r15, method = "delong", paired = FALSE)
ci5 <- pROC::ci.auc(r5, method = "delong"); ci15 <- pROC::ci.auc(r15, method = "delong")
delong <- data.frame(definition = c("AHI>=5 (primary)", "AHI>=15 (sensitivity)"), n = c(nrow(def5), nrow(def15)),
                     COMISA_n = c(sum(def5$COMISA5), sum(def15$COMISA15)),
                     AUC = c(as.numeric(pROC::auc(r5)), as.numeric(pROC::auc(r15))), AUC_lo = c(ci5[1], ci15[1]), AUC_hi = c(ci5[3], ci15[3]),
                     DeLong_D = unname(dl$statistic), DeLong_p = dl$p.value, test = "unpaired DeLong")
write.csv(delong, file.path(OUT, "TableS12_DeLong_def5_vs_def15_definitive.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 33_C1_composition
lev <- c("COMISA", "Insomnia disorder only", "Neither", "OSA only")
ct  <- table(d$Cluster, factor(d$Group5, levels = lev))
rowcounts <- as.data.frame.matrix(ct); rowcounts$Total <- rowSums(ct); rowcounts <- cbind(Cluster = rownames(rowcounts), rowcounts)
rowpct <- as.data.frame.matrix(round(100 * prop.table(ct, 1), 1)); rowpct <- cbind(Cluster = rownames(rowpct), rowpct)
write.csv(rowcounts, file.path(OUT, "Table3A_rowcounts.csv"), row.names = FALSE)
write.csv(rowpct,    file.path(OUT, "Table3A_rowpct.csv"),    row.names = FALSE)

# ---------------------------------------------------------------- Check against the values reported in the manuscript (real data)
chk <- function(label, got, want) message(sprintf("%-38s %s  (expected %s) %s", label, paste(got, collapse = "/"), paste(want, collapse = "/"), ifelse(isTRUE(all.equal(got, want, tolerance = 0.01)), "OK", "MISMATCH")))
chk("women in C1 / C2",           sex_by_cluster$women,        c(112, 46))
chk("sex Cramer V",               round(sex_by_cluster$cramer_V[1], 2), 0.36)
chk("C1 by stream (IS / SAS)",    stream_by_cluster$C1,         c(103, 144))
chk("stream Cramer V",            round(stream_by_cluster$cramer_V[1], 2), 0.48)
chk("women within SAS (C1 / C2)", sex_within_sas$women,         c(54, 42))
chk("AUC def5 / def15",           round(delong$AUC, 2),         c(0.80, 0.77))
chk("DeLong p",                   round(delong$DeLong_p[1], 2), 0.37)
chk("Table 3A row C1",            as.numeric(rowcounts[1, lev]), c(84, 19, 21, 123))
chk("Table 3A row C2",            as.numeric(rowcounts[2, lev]), c(10, 0, 1, 341))
cat("Done. Six files written to", OUT, "\n")
