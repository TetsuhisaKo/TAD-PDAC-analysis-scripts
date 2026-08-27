#!/usr/bin/env Rscript
# ============================================================
# 06_tcga_paad_analysis.R
#
# TCGA-PAAD GDF15 analyses.
#
# Computes the statistics underlying:
#   Fig. 4c                    GDF15 high/low vs CD8A/GZMA/GZMB/PRF1
#   Fig. 4d                    continuous GDF15 vs CD8A
#   Fig. 4e                    continuous GDF15 vs DDIT3
#   Supplementary Table S10A-D
#   Supplementary Fig. S5      GDF15 high/low overall survival
#
# This script outputs analysis results; final figure assembly is performed separately.
#
# Expression cohort:
#   n=163 primary tumors
#   143 infiltrating duct carcinoma, NOS
#    20 adenocarcinoma, NOS
#
# GDF15 split:
#   expression >= median -> High
#   n=163: Low 81 / High 82
#   OS n=162: Low 81 / High 81
#
# ABSOLUTE purity:
#   n=143
#   duct-only n=143; purity n=123
#
# BH families:
#   six median-split comparisons
#   six purity-adjusted partial correlations
#   six duct-only purity-adjusted partial correlations
#
# Expression-scale handling:
# The UCSC Xena matrix used here already contains the log2-transformed
# FPKM-UQ values used in the manuscript; values are therefore used directly
# without an additional log2(x+1) transformation.
#
# Sources:
#   UCSC Xena GDC Hub, accessed 2026-04-28
#   PanCanAtlas ABSOLUTE, accessed 2026-07-31
#   GDC UUID 4f277128-f793-4354-a13d-30cc7fe9f6b5
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

data_dir <- "data/TCGA"
out_dir <- "results/TCGA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

expr_file <- file.path(data_dir, "TCGA-PAAD.star_fpkm-uq.tsv.gz")
clin_file <- file.path(data_dir, "TCGA-PAAD.clinical.tsv.gz")
surv_file <- file.path(data_dir, "TCGA-PAAD.survival.tsv.gz")
probe_file <- file.path(data_dir, "gencode.v36.annotation.gtf.gene.probemap")
purity_file <- file.path(data_dir, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")

input_files <- c(
  expression = expr_file,
  clinical = clin_file,
  survival = surv_file,
  probe_map = probe_file,
  ABSOLUTE_purity = purity_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}

# ------------------------------------------------------------
# Input-file MD5 manifest
# ------------------------------------------------------------

md5_manifest <- data.frame(
  File = basename(input_files),
  MD5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)

write.csv(
  md5_manifest,
  file.path(out_dir, "TCGA_input_file_md5.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Expression matrix
# ------------------------------------------------------------

expr_raw <- read.delim(
  gzfile(expr_file),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

id_col <- names(expr_raw)[1]

# Remove PAR_Y records BEFORE stripping version suffixes.
expr_raw <- expr_raw[
  !grepl("_PAR_Y", expr_raw[[id_col]], fixed = TRUE),
  ,
  drop = FALSE
]

expr_raw$Ensembl_ID_clean <- sub("\\..*$", "", expr_raw[[id_col]])

dup_n <- sum(duplicated(expr_raw$Ensembl_ID_clean))
cat("Duplicated cleaned Ensembl IDs removed:", dup_n, "\n")

expr_raw <- expr_raw[
  !duplicated(expr_raw$Ensembl_ID_clean),
  ,
  drop = FALSE
]

all_sample_cols <- setdiff(
  names(expr_raw),
  c(id_col, "Ensembl_ID_clean")
)

# Explicit 01A restriction avoids unintentionally including 01B/01C.
sample_cols <- all_sample_cols[grepl("-01A$", all_sample_cols)]

cat(
  "Expression columns:",
  length(all_sample_cols), "total;",
  length(sample_cols), "with -01A suffix\n"
)

# ------------------------------------------------------------
# GENCODE v36 probe map
# ------------------------------------------------------------

probe <- read.delim(
  probe_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

stopifnot(all(c("id", "gene") %in% names(probe)))

probe <- probe[
  !grepl("_PAR_Y", probe$id, fixed = TRUE),
  ,
  drop = FALSE
]

probe$Ensembl_ID_clean <- sub("\\..*$", "", probe$id)

genes_needed <- c("GDF15", "CD8A", "GZMA", "GZMB", "PRF1", "DDIT3")

resolve_gene <- function(symbol) {
  ids <- unique(probe$Ensembl_ID_clean[probe$gene == symbol])
  ids <- ids[ids %in% expr_raw$Ensembl_ID_clean]

  if (length(ids) == 0) {
    stop("No expression Ensembl ID found for ", symbol)
  }

  if (length(ids) > 1) {
    warning(
      "Multiple Ensembl IDs available for ", symbol,
      "; using ", ids[1]
    )
  }

  ids[1]
}

gene_map <- data.frame(
  Gene = genes_needed,
  Ensembl_ID = vapply(genes_needed, resolve_gene, character(1)),
  stringsAsFactors = FALSE
)

write.csv(
  gene_map,
  file.path(out_dir, "TCGA_gene_mapping_used.csv"),
  row.names = FALSE
)

cat("\nEnsembl IDs used:\n")
print(gene_map)

extract_gene <- function(symbol) {
  eid <- gene_map$Ensembl_ID[gene_map$Gene == symbol]

  x <- expr_raw[
    expr_raw$Ensembl_ID_clean == eid,
    sample_cols,
    drop = FALSE
  ]

  if (nrow(x) != 1) {
    stop("Expected exactly one expression row for ", symbol)
  }

  data.frame(
    sample = sample_cols,
    value = as.numeric(x[1, ]),
    stringsAsFactors = FALSE
  )
}

expr_gene <- lapply(genes_needed, extract_gene)
names(expr_gene) <- genes_needed

expr_df <- data.frame(sample = sample_cols, stringsAsFactors = FALSE)

for (g in genes_needed) {
  tmp <- expr_gene[[g]]
  names(tmp)[2] <- g
  expr_df <- inner_join(expr_df, tmp, by = "sample")
}

# ------------------------------------------------------------
# Clinical histology restriction
# ------------------------------------------------------------

clin <- read.delim(
  gzfile(clin_file),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

stopifnot(
  all(c("sample", "primary_diagnosis.diagnoses") %in% names(clin))
)

pdac_histology <- c(
  "Infiltrating duct carcinoma, NOS",
  "Adenocarcinoma, NOS"
)

clin_pdac <- clin %>%
  filter(
    grepl("-01A$", sample),
    primary_diagnosis.diagnoses %in% pdac_histology
  ) %>%
  select(sample, primary_diagnosis.diagnoses) %>%
  distinct(sample, .keep_all = TRUE)

tcga <- inner_join(expr_df, clin_pdac, by = "sample") %>%
  mutate(
    CYT = (GZMA + PRF1) / 2
  )

stopifnot(nrow(tcga) == 163)

stopifnot(
  sum(
    tcga$primary_diagnosis.diagnoses ==
      "Infiltrating duct carcinoma, NOS"
  ) == 143
)

stopifnot(
  sum(
    tcga$primary_diagnosis.diagnoses ==
      "Adenocarcinoma, NOS"
  ) == 20
)

cat("\nExpression cohort histology:\n")
print(table(tcga$primary_diagnosis.diagnoses))

cat("\nCD8A range:", paste(range(tcga$CD8A), collapse = " to "), "\n")
cat("GDF15 range:", paste(range(tcga$GDF15), collapse = " to "), "\n")

# ------------------------------------------------------------
# GDF15 median split in n=163
# Median itself belongs to High.
# ------------------------------------------------------------

gdf15_median <- median(tcga$GDF15, na.rm = TRUE)

tcga <- tcga %>%
  mutate(
    GDF15_group = factor(
      ifelse(GDF15 >= gdf15_median, "High", "Low"),
      levels = c("Low", "High")
    )
  )

stopifnot(sum(tcga$GDF15_group == "Low") == 81)
stopifnot(sum(tcga$GDF15_group == "High") == 82)

cat("\nGDF15 median:", gdf15_median, "\n")
print(table(tcga$GDF15_group))

outcomes <- c("CD8A", "GZMA", "GZMB", "PRF1", "CYT", "DDIT3")

# ------------------------------------------------------------
# S10A: six median-split comparisons; BH across six
# ------------------------------------------------------------

median_split_results <- bind_rows(
  lapply(outcomes, function(v) {

    x_low <- tcga[[v]][tcga$GDF15_group == "Low"]
    x_high <- tcga[[v]][tcga$GDF15_group == "High"]

    wt <- wilcox.test(
      x_low,
      x_high,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE
    )

    data.frame(
      Transcript_or_score = v,
      Low_n = sum(!is.na(x_low)),
      High_n = sum(!is.na(x_high)),
      Low_median = median(x_low, na.rm = TRUE),
      High_median = median(x_high, na.rm = TRUE),
      P_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
) %>%
  mutate(q_value = p.adjust(P_value, method = "BH"))

write.csv(
  median_split_results,
  file.path(out_dir, "Supplementary_Table_S10A_median_split.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# S10B: continuous Spearman correlations, n=163
# ------------------------------------------------------------

spearman_results <- bind_rows(
  lapply(outcomes, function(v) {

    ct <- cor.test(
      tcga$GDF15,
      tcga[[v]],
      method = "spearman",
      exact = FALSE
    )

    data.frame(
      Transcript_or_score = v,
      n = sum(complete.cases(tcga[, c("GDF15", v)])),
      Spearman_rho = unname(ct$estimate),
      P_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  spearman_results,
  file.path(out_dir, "Supplementary_Table_S10B_continuous_Spearman.csv"),
  row.names = FALSE
)

# Reference estimates reported in the manuscript:
# CD8A  rho -0.371, P 1.12e-6
# GZMA  rho -0.363, P 1.94e-6
# GZMB  rho -0.347, P 5.53e-6
# PRF1  rho -0.408, P 6.45e-8
# CYT   rho -0.404, P 8.81e-8
# DDIT3 rho  0.529, P 4.14e-13

# ------------------------------------------------------------
# ABSOLUTE tumor purity
# ------------------------------------------------------------

abs_dat <- read.delim(
  purity_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

required_abs <- c("array", "call status", "purity")
stopifnot(all(required_abs %in% names(abs_dat)))

purity <- abs_dat %>%
  filter(
    `call status` == "called",
    !is.na(purity),
    purity != ""
  ) %>%
  transmute(
    sample_key = as.character(array),
    ABSOLUTE_purity = as.numeric(purity)
  ) %>%
  filter(!is.na(ABSOLUTE_purity)) %>%
  distinct(sample_key, .keep_all = TRUE)

tcga <- tcga %>%
  mutate(sample_key = substr(sample, 1, 15)) %>%
  left_join(purity, by = "sample_key")

stopifnot(sum(!is.na(tcga$ABSOLUTE_purity)) == 143)

purity_keep <- !is.na(tcga$ABSOLUTE_purity)

purity_cor <- cor.test(
  tcga$GDF15[purity_keep],
  tcga$ABSOLUTE_purity[purity_keep],
  method = "spearman",
  exact = FALSE
)

purity_gdf15 <- data.frame(
  n = sum(purity_keep),
  Spearman_rho = unname(purity_cor$estimate),
  P_value = purity_cor$p.value
)

write.csv(
  purity_gdf15,
  file.path(out_dir, "GDF15_vs_ABSOLUTE_purity.csv"),
  row.names = FALSE
)

cat("\nGDF15 vs ABSOLUTE purity:\n")
print(purity_gdf15)
cat("Reference estimate: rho ~0.18; P~0.031.\n")

# ------------------------------------------------------------
# Partial Spearman by rank residualization.
#
# With k adjustment covariates:
#   t-test df = n - k - 2
#   Fisher-z CI SE = 1 / sqrt(n - k - 3)
#
# Here k=1 (ABSOLUTE purity):
#   df = n - 3
#   Fisher-z SE = 1 / sqrt(n - 4)
# ------------------------------------------------------------

partial_spearman <- function(x, y, z) {

  keep <- complete.cases(x, y, z)

  xr <- rank(x[keep], ties.method = "average")
  yr <- rank(y[keep], ties.method = "average")
  zr <- rank(z[keep], ties.method = "average")

  rx <- resid(lm(xr ~ zr))
  ry <- resid(lm(yr ~ zr))

  rho <- cor(rx, ry, method = "pearson")
  n <- length(rx)
  k <- 1L

  if (n <= k + 3) {
    stop("Insufficient observations for partial-correlation CI.")
  }

  dfree <- n - k - 2
  t_value <- rho * sqrt(dfree / (1 - rho^2))
  p_value <- 2 * pt(-abs(t_value), df = dfree)

  bounded_rho <- max(min(rho, 0.999999), -0.999999)
  z_rho <- atanh(bounded_rho)

  se_z <- 1 / sqrt(n - k - 3)
  zcrit <- qnorm(0.975)

  ci_low <- tanh(z_rho - zcrit * se_z)
  ci_high <- tanh(z_rho + zcrit * se_z)

  data.frame(
    n = n,
    k_covariates = k,
    Partial_rho = rho,
    CI_low = ci_low,
    CI_high = ci_high,
    P_value = p_value
  )
}

# ------------------------------------------------------------
# S10C: purity-adjusted partial correlations; BH across six
# ------------------------------------------------------------

purity_subset <- tcga %>%
  filter(!is.na(ABSOLUTE_purity))

stopifnot(nrow(purity_subset) == 143)

partial_results <- bind_rows(
  lapply(outcomes, function(v) {

    x <- partial_spearman(
      purity_subset$GDF15,
      purity_subset[[v]],
      purity_subset$ABSOLUTE_purity
    )

    cbind(
      data.frame(Transcript_or_score = v),
      x
    )
  })
) %>%
  mutate(q_value = p.adjust(P_value, method = "BH"))

write.csv(
  partial_results,
  file.path(out_dir, "Supplementary_Table_S10C_purity_adjusted.csv"),
  row.names = FALSE
)

# Reference point estimates and P values reported in the manuscript:
# CD8A  -0.335, P 4.64e-5
# GZMA  -0.332, P 5.48e-5
# GZMB  -0.341, P 3.30e-5
# PRF1  -0.360, P 1.11e-5
# CYT   -0.363, P 8.92e-6
# DDIT3  0.506, P 1.36e-10
# all q<0.001.
#
# Fisher-z confidence intervals for the manuscript-highlighted values:
# CD8A approximately -0.47 to -0.18
# CYT  approximately -0.50 to -0.21
# DDIT3 approximately 0.37 to 0.62

# ------------------------------------------------------------
# S10D: infiltrating duct carcinoma, NOS only
# ------------------------------------------------------------

duct <- tcga %>%
  filter(
    primary_diagnosis.diagnoses ==
      "Infiltrating duct carcinoma, NOS"
  )

stopifnot(nrow(duct) == 143)
stopifnot(sum(!is.na(duct$ABSOLUTE_purity)) == 123)

duct_results <- bind_rows(
  lapply(outcomes, function(v) {

    ct <- cor.test(
      duct$GDF15,
      duct[[v]],
      method = "spearman",
      exact = FALSE
    )

    pa <- partial_spearman(
      duct$GDF15,
      duct[[v]],
      duct$ABSOLUTE_purity
    )

    data.frame(
      Transcript_or_score = v,
      Spearman_n = sum(complete.cases(duct[, c("GDF15", v)])),
      Spearman_rho = unname(ct$estimate),
      Spearman_P = ct$p.value,
      Partial_n = pa$n,
      Partial_rho = pa$Partial_rho,
      CI_low = pa$CI_low,
      CI_high = pa$CI_high,
      Partial_P = pa$P_value,
      stringsAsFactors = FALSE
    )
  })
) %>%
  mutate(Partial_q = p.adjust(Partial_P, method = "BH"))

write.csv(
  duct_results,
  file.path(out_dir, "Supplementary_Table_S10D_duct_only_sensitivity.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Fig. S5: overall survival
# Median split remains defined in n=163 before survival filtering.
# ------------------------------------------------------------

surv <- read.delim(
  gzfile(surv_file),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

stopifnot(all(c("sample", "OS.time", "OS") %in% names(surv)))

surv_primary <- surv %>%
  filter(grepl("-01A$", sample)) %>%
  select(sample, OS.time, OS) %>%
  distinct(sample, .keep_all = TRUE)

surv_analysis <- tcga %>%
  left_join(surv_primary, by = "sample") %>%
  filter(!is.na(OS.time), !is.na(OS)) %>%
  mutate(
    OS = as.integer(OS),
    OS_months = as.numeric(OS.time) / 30.44
  )

stopifnot(nrow(surv_analysis) == 162)
stopifnot(sum(surv_analysis$GDF15_group == "Low") == 81)
stopifnot(sum(surv_analysis$GDF15_group == "High") == 81)

missing_os <- tcga %>%
  anti_join(
    surv_analysis %>% select(sample),
    by = "sample"
  )

stopifnot(nrow(missing_os) == 1)
stopifnot(as.character(missing_os$GDF15_group[1]) == "High")

fit_lr <- survdiff(
  Surv(OS_months, OS) ~ GDF15_group,
  data = surv_analysis
)

lr_df <- length(fit_lr$n) - 1
logrank_p <- 1 - pchisq(fit_lr$chisq, df = lr_df)

s5 <- data.frame(
  n = nrow(surv_analysis),
  Low_n = sum(surv_analysis$GDF15_group == "Low"),
  High_n = sum(surv_analysis$GDF15_group == "High"),
  Events = sum(surv_analysis$OS),
  Logrank_chisq = unname(fit_lr$chisq),
  Logrank_df = lr_df,
  Logrank_P = logrank_p
)

write.csv(
  s5,
  file.path(out_dir, "Supplementary_Figure_S5_OS_statistics.csv"),
  row.names = FALSE
)

write.csv(
  tcga,
  file.path(out_dir, "TCGA_PAAD_histology_restricted_n163.csv"),
  row.names = FALSE
)

write.csv(
  surv_analysis,
  file.path(out_dir, "TCGA_PAAD_OS_evaluable_n162.csv"),
  row.names = FALSE
)

cat("\n============================================\n")
cat("TCGA ANALYSIS SUMMARY\n")
cat("============================================\n")
cat("Expression n =", nrow(tcga), "\n")
print(table(tcga$primary_diagnosis.diagnoses))
cat("GDF15 median =", gdf15_median, "\n")
print(table(tcga$GDF15_group))
cat("\nS10B continuous correlations:\n")
print(spearman_results)
cat("\nS10C purity-adjusted:\n")
print(partial_results)
cat("\nS10A median split:\n")
print(median_split_results)
cat("\nS10D duct-only:\n")
print(duct_results)
cat("\nSupplementary Fig. S5 OS:\n")
print(s5)
cat("Reference log-rank P approximately 0.75.\n")
cat("============================================\n")

sink(file.path(out_dir, "sessionInfo_TCGA.txt"))
sessionInfo()
sink()
