# ============================================================
# TCGA-PAAD analysis: GDF15, cytotoxic immune genes, CYT, DDIT3, and survival
#
# Maps to:
#   Fig. 5a: CD8A, GZMA, GZMB, PRF1 by GDF15-high/low
#   Fig. 5b: CYT by GDF15-high/low
#   Fig. 5c: GDF15 vs CD8A
#   Fig. 5d: GDF15 vs DDIT3/CHOP
#   Supplementary Fig. S7: overall survival by GDF15-high/low
#
# Source:
#   UCSC Xena GDC Hub, accessed 2026-04-28
#
# Input files:
#   data/TCGA/TCGA-PAAD.star_fpkm-uq.tsv.gz
#   data/TCGA/TCGA-PAAD.survival.tsv.gz
#   data/TCGA/gencode.v36.annotation.gtf.gene.probemap
#
# Expression values:
#   log2(FPKM-UQ + 1)
#
# Expected manuscript cohort:
#   n = 177 primary tumors, 93 OS events
#   GDF15-low n = 88, GDF15-high n = 89
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

# ---- Paths ----
data_dir      <- "data/TCGA"
expr_file     <- file.path(data_dir, "TCGA-PAAD.star_fpkm-uq.tsv.gz")
survival_file <- file.path(data_dir, "TCGA-PAAD.survival.tsv.gz")
probe_map     <- file.path(data_dir, "gencode.v36.annotation.gtf.gene.probemap")

stopifnot(file.exists(expr_file))
stopifnot(file.exists(survival_file))
stopifnot(file.exists(probe_map))

# ---- Load and clean expression matrix ----
expr_raw <- read.table(
  gzfile(expr_file),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

id_col <- colnames(expr_raw)[1]
expr_raw$id_clean <- sub("\\..*", "", expr_raw[[id_col]])
expr_raw <- expr_raw[!duplicated(expr_raw$id_clean), ]
rownames(expr_raw) <- expr_raw$id_clean

expr_mat <- expr_raw[, !colnames(expr_raw) %in% c(id_col, "id_clean")]

# ---- Probe map: Ensembl ID to gene symbol ----
probe <- read.table(
  probe_map,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

probe$id_clean <- sub("\\..*", "", probe$id)

# ---- Helper function: extract log2(FPKM-UQ + 1) expression ----
get_log2 <- function(gene, samples) {
  ids <- probe$id_clean[probe$gene == gene]

  if (length(ids) == 0) {
    stop(paste("Gene not found in probe map:", gene))
  }

  eid <- ids[1]

  if (!(eid %in% rownames(expr_mat))) {
    stop(paste("Ensembl ID not found in expression matrix:", eid, "for gene:", gene))
  }

  vals <- log2(as.numeric(expr_mat[eid, samples]) + 1)
  names(vals) <- samples
  vals
}

# ---- Primary tumor samples ----
tumor_expr <- colnames(expr_mat)[grepl("-01A$", colnames(expr_mat))]

# ---- Load survival data ----
surv <- read.table(
  gzfile(survival_file),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

surv <- surv[grepl("-01A$", surv$sample), ]
surv$patient <- substr(surv$sample, 1, 12)

# ---- Build expression data frame ----
expr_df <- data.frame(
  sample  = tumor_expr,
  patient = substr(tumor_expr, 1, 12),
  GDF15   = get_log2("GDF15", tumor_expr),
  CD8A    = get_log2("CD8A",  tumor_expr),
  GZMA    = get_log2("GZMA",  tumor_expr),
  GZMB    = get_log2("GZMB",  tumor_expr),
  PRF1    = get_log2("PRF1",  tumor_expr),
  DDIT3   = get_log2("DDIT3", tumor_expr),
  stringsAsFactors = FALSE
)

stopifnot(!any(duplicated(expr_df$patient)))

# ---- Merge expression and survival data ----
df <- merge(
  expr_df,
  surv[, c("patient", "OS.time", "OS")],
  by = "patient"
)

df <- df[!is.na(df$OS.time) & !is.na(df$OS), ]

stopifnot(!any(duplicated(df$patient)))

cat("Analytic cohort: n =", nrow(df),
    "| events =", sum(df$OS), "\n")

# ---- GDF15 high/low group based on cohort median ----
df$GDF15_group <- factor(
  ifelse(df$GDF15 >= median(df$GDF15), "High", "Low"),
  levels = c("Low", "High")
)

cat("GDF15 group counts:\n")
print(table(df$GDF15_group))

# ---- CYT score ----
# CYT was calculated as the arithmetic mean of log2-transformed
# GZMA and PRF1 expression values, adapted from the cytolytic
# activity score described by Rooney et al.
df$CYT <- (df$GZMA + df$PRF1) / 2

# ---- Overall survival time in months ----
df$OS_months <- df$OS.time / 30.44

# ============================================================
# Fig. 5a: immune effector genes by GDF15 group
# Two-sided Mann-Whitney U test / Wilcoxon rank-sum test
# R wilcox.test uses normal approximation for this sample size.
# ============================================================

effector_genes <- c("CD8A", "GZMA", "GZMB", "PRF1")

fig5a_results <- do.call(rbind, lapply(effector_genes, function(g) {
  wt <- wilcox.test(
    df[[g]] ~ df$GDF15_group,
    exact = FALSE,
    correct = TRUE
  )

  data.frame(
    Analysis = "Fig5a",
    Gene = g,
    Low_median = median(df[[g]][df$GDF15_group == "Low"], na.rm = TRUE),
    High_median = median(df[[g]][df$GDF15_group == "High"], na.rm = TRUE),
    P_value = wt$p.value
  )
}))

cat("\nFig. 5a results:\n")
print(fig5a_results)

# Expected approximate manuscript values:
# CD8A p = 0.00245
# GZMA p = 0.0068
# GZMB p = 0.00861
# PRF1 p = 0.000414

# ============================================================
# Fig. 5b: CYT by GDF15 group
# ============================================================

wt_cyt <- wilcox.test(
  CYT ~ GDF15_group,
  data = df,
  exact = FALSE,
  correct = TRUE
)

fig5b_result <- data.frame(
  Analysis = "Fig5b",
  Variable = "CYT",
  Low_median = median(df$CYT[df$GDF15_group == "Low"], na.rm = TRUE),
  High_median = median(df$CYT[df$GDF15_group == "High"], na.rm = TRUE),
  P_value = wt_cyt$p.value
)

cat("\nFig. 5b result:\n")
print(fig5b_result)

# ============================================================
# Fig. 5c: GDF15 vs CD8A Spearman correlation
# ============================================================

rho_cd8 <- cor.test(
  df$GDF15,
  df$CD8A,
  method = "spearman",
  exact = FALSE
)

fig5c_result <- data.frame(
  Analysis = "Fig5c",
  Comparison = "GDF15_vs_CD8A",
  Spearman_rho = unname(rho_cd8$estimate),
  P_value = rho_cd8$p.value
)

cat("\nFig. 5c result:\n")
print(fig5c_result)

# Expected approximate manuscript value:
# rho ≈ -0.28, p ≈ 0.009

# ============================================================
# Fig. 5d: GDF15 vs DDIT3/CHOP Spearman correlation
# ============================================================

rho_ddit3 <- cor.test(
  df$GDF15,
  df$DDIT3,
  method = "spearman",
  exact = FALSE
)

fig5d_result <- data.frame(
  Analysis = "Fig5d",
  Comparison = "GDF15_vs_DDIT3_CHOP",
  Spearman_rho = unname(rho_ddit3$estimate),
  P_value = rho_ddit3$p.value
)

cat("\nFig. 5d result:\n")
print(fig5d_result)

# Expected approximate manuscript value:
# rho ≈ 0.47, p < 0.001

# ============================================================
# Supplementary Fig. S7: overall survival by GDF15 group
# ============================================================

lr <- survdiff(Surv(OS_months, OS) ~ GDF15_group, data = df)
logrank_p <- 1 - pchisq(lr$chisq, df = 1)

s7_result <- data.frame(
  Analysis = "Supplementary_Fig_S7",
  Test = "Log-rank",
  Chi_square = unname(lr$chisq),
  P_value = logrank_p
)

cat("\nSupplementary Fig. S7 result:\n")
print(s7_result)

# Expected manuscript value:
# log-rank p ≈ 0.154

# ---- Optional: save results ----
dir.create("results", showWarnings = FALSE)

write.csv(fig5a_results, "results/TCGA_Fig5a_results.csv", row.names = FALSE)
write.csv(fig5b_result,  "results/TCGA_Fig5b_CYT_result.csv", row.names = FALSE)
write.csv(fig5c_result,  "results/TCGA_Fig5c_GDF15_CD8A_result.csv", row.names = FALSE)
write.csv(fig5d_result,  "results/TCGA_Fig5d_GDF15_DDIT3_result.csv", row.names = FALSE)
write.csv(s7_result,     "results/TCGA_SuppFigS7_survival_result.csv", row.names = FALSE)

# ---- Record R/package versions ----
sink("results/sessionInfo_TCGA_PAAD_GDF15_analysis.txt")
sessionInfo()
sink()
