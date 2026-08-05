# ============================================================
# 05_tcga_gdf15_analysis.R
#
# Current TCGA-PAAD analyses for the revised manuscript.
#
# Maps to:
#   Fig. 4c:
#     GDF15-high vs low: CD8A, GZMA, GZMB, PRF1
#   Fig. 4d:
#     continuous GDF15 vs CD8A
#   Fig. 4e:
#     continuous GDF15 vs DDIT3/CHOP
#   Supplementary Table S8:
#     A. unadjusted continuous Spearman correlations
#     B. ABSOLUTE-purity-adjusted partial Spearman correlations
#     C. six GDF15 median-split comparisons
#     D. infiltrating duct carcinoma, NOS sensitivity analysis
#   Supplementary Fig. S5:
#     GDF15-high vs low overall survival
#
# Cohorts:
#   expression n=163:
#     143 Infiltrating duct carcinoma, NOS
#      20 Adenocarcinoma, NOS
#   OS-evaluable n=162:
#      81 GDF15-low / 81 GDF15-high
#   ABSOLUTE purity-evaluable n=143
#   duct-carcinoma-only n=143; purity-evaluable n=123
#
# Multiple testing:
#   BH across SIX median-split comparisons:
#     CD8A, GZMA, GZMB, PRF1, CYT, DDIT3
#   BH separately across SIX purity-adjusted partial correlations.
#   Duct-only purity-adjusted analyses form another separate family of six.
#
# IMPORTANT EXPRESSION-SCALE NOTE:
# The UCSC Xena matrix used here contains the log2-transformed
# FPKM-UQ values used in the manuscript. Values are used directly.
# DO NOT apply another log2(x+1) transformation.
#
# Public sources:
#   UCSC Xena GDC Hub, accessed 2026-04-28
#   PanCanAtlas ABSOLUTE file, accessed 2026-07-31
#   GDC UUID:
#   4f277128-f793-4354-a13d-30cc7fe9f6b5
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
# Input-file checksums
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

# If cleaning creates duplicates, keep the first and record the fact.
dup_n <- sum(duplicated(expr_raw$Ensembl_ID_clean))
cat("Duplicated cleaned Ensembl IDs removed:", dup_n, "\n")

expr_raw <- expr_raw[
  !duplicated(expr_raw$Ensembl_ID_clean),
  ,
  drop = FALSE
]

sample_cols <- setdiff(
  names(expr_raw),
  c(id_col, "Ensembl_ID_clean")
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

extract_gene <- function(symbol) {
  eid <- gene_map$Ensembl_ID[gene_map$Gene == symbol]
  x <- expr_raw[
    expr_raw$Ensembl_ID_clean == eid,
    sample_cols,
    drop = FALSE
  ]
  if (nrow(x) != 1) stop("Expected one expression row for ", symbol)

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
    substr(sample, 14, 15) == "01",
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
  sum(tcga$primary_diagnosis.diagnoses ==
        "Infiltrating duct carcinoma, NOS") == 143
)
stopifnot(
  sum(tcga$primary_diagnosis.diagnoses ==
        "Adenocarcinoma, NOS") == 20
)

cat("\nExpression cohort histology:\n")
print(table(tcga$primary_diagnosis.diagnoses))

# Expression-scale sanity check against the current analysis snapshot.
cat("\nCD8A range:", paste(range(tcga$CD8A), collapse = " to "), "\n")
cat("GDF15 range:", paste(range(tcga$GDF15), collapse = " to "), "\n")
# Current CD8A range is approximately 0.2028878 to 4.6617610.

# ------------------------------------------------------------
# GDF15 median split in n=163
# ------------------------------------------------------------

gdf15_median <- median(tcga$GDF15, na.rm = TRUE)

tcga <- tcga %>%
  mutate(
    GDF15_group = factor(
      ifelse(GDF15 >= gdf15_median, "High", "Low"),
      levels = c("Low", "High")
    )
  )

cat("\nGDF15 median:", gdf15_median, "\n")
print(table(tcga$GDF15_group))

stopifnot(sum(tcga$GDF15_group == "Low") == 81)
stopifnot(sum(tcga$GDF15_group == "High") == 82)

# Current dataset snapshot: median approximately 4.19693.

outcomes <- c("CD8A", "GZMA", "GZMB", "PRF1", "CYT", "DDIT3")

# ------------------------------------------------------------
# S8C: median-split comparisons; BH across six
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
      P_value = wt$p.value
    )
  })
) %>%
  mutate(
    q_value = p.adjust(P_value, method = "BH")
  )

write.csv(
  median_split_results,
  file.path(out_dir, "Supplementary_Table_S8C_median_split.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# S8A: continuous Spearman correlations in n=163
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
      P_value = ct$p.value
    )
  })
)

write.csv(
  spearman_results,
  file.path(out_dir, "Supplementary_Table_S8A_continuous_Spearman.csv"),
  row.names = FALSE
)

# Current expected values:
# CD8A rho -0.371, P 1.12e-6
# GZMA rho -0.363, P 1.94e-6
# GZMB rho -0.347, P 5.53e-6
# PRF1 rho -0.408, P 6.45e-8
# CYT  rho -0.404, P 8.81e-8
# DDIT3 rho 0.529, P 4.14e-13

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

purity_cor <- cor.test(
  tcga$GDF15[!is.na(tcga$ABSOLUTE_purity)],
  tcga$ABSOLUTE_purity[!is.na(tcga$ABSOLUTE_purity)],
  method = "spearman",
  exact = FALSE
)

purity_gdf15 <- data.frame(
  n = sum(!is.na(tcga$ABSOLUTE_purity)),
  Spearman_rho = unname(purity_cor$estimate),
  P_value = purity_cor$p.value
)

write.csv(
  purity_gdf15,
  file.path(out_dir, "GDF15_vs_ABSOLUTE_purity.csv"),
  row.names = FALSE
)

# Current expected: rho ~0.18, P~0.031.

# ------------------------------------------------------------
# Partial Spearman by rank residualization
# One covariate (purity): t-test uses df=n-3.
# Fisher-z CI uses SE=1/sqrt(n-3).
# ------------------------------------------------------------

partial_spearman <- function(x, y, z) {
  keep <- complete.cases(x, y, z)
  x <- rank(x[keep], ties.method = "average")
  y <- rank(y[keep], ties.method = "average")
  z <- rank(z[keep], ties.method = "average")

  rx <- resid(lm(x ~ z))
  ry <- resid(lm(y ~ z))

  rho <- cor(rx, ry, method = "pearson")
  n <- length(rx)

  dfree <- n - 3
  t_value <- rho * sqrt(dfree / (1 - rho^2))
  p_value <- 2 * pt(-abs(t_value), df = dfree)

  z_rho <- atanh(max(min(rho, 0.999999), -0.999999))
  se <- 1 / sqrt(n - 3)
  zcrit <- qnorm(0.975)

  ci_low <- tanh(z_rho - zcrit * se)
  ci_high <- tanh(z_rho + zcrit * se)

  data.frame(
    n = n,
    Partial_rho = rho,
    CI_low = ci_low,
    CI_high = ci_high,
    P_value = p_value
  )
}

# ------------------------------------------------------------
# S8B: purity-adjusted correlations in n=143; BH across six
# ------------------------------------------------------------

purity_subset <- tcga %>% filter(!is.na(ABSOLUTE_purity))
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
  mutate(
    q_value = p.adjust(P_value, method = "BH")
  )

write.csv(
  partial_results,
  file.path(out_dir, "Supplementary_Table_S8B_purity_adjusted.csv"),
  row.names = FALSE
)

# Current expected:
# CD8A partial rho -0.335, P 4.64e-5
# GZMA partial rho -0.332, P 5.48e-5
# GZMB partial rho -0.341, P 3.30e-5
# PRF1 partial rho -0.360, P 1.11e-5
# CYT  partial rho -0.363, P 8.92e-6
# DDIT3 partial rho  0.506, P 1.36e-10
# all BH q<0.001.

# ------------------------------------------------------------
# S8D: infiltrating duct carcinoma, NOS only
# n=143; purity n=123
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
      Partial_P = pa$P_value
    )
  })
) %>%
  mutate(
    Partial_q = p.adjust(Partial_P, method = "BH")
  )

write.csv(
  duct_results,
  file.path(out_dir, "Supplementary_Table_S8D_duct_only_sensitivity.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Fig. S5: OS
# Median split remains the split defined above in n=163.
# ------------------------------------------------------------

surv <- read.delim(
  gzfile(surv_file),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

stopifnot(all(c("sample", "OS.time", "OS") %in% names(surv)))

surv_analysis <- tcga %>%
  left_join(
    surv %>%
      select(sample, OS.time, OS) %>%
      distinct(sample, .keep_all = TRUE),
    by = "sample"
  ) %>%
  filter(!is.na(OS.time), !is.na(OS)) %>%
  mutate(OS_months = as.numeric(OS.time) / 30.44)

stopifnot(nrow(surv_analysis) == 162)
stopifnot(sum(surv_analysis$GDF15_group == "Low") == 81)
stopifnot(sum(surv_analysis$GDF15_group == "High") == 81)

fit_lr <- survdiff(
  Surv(OS_months, OS) ~ GDF15_group,
  data = surv_analysis
)

logrank_p <- 1 - pchisq(fit_lr$chisq, df = 1)

s5 <- data.frame(
  n = nrow(surv_analysis),
  Low_n = sum(surv_analysis$GDF15_group == "Low"),
  High_n = sum(surv_analysis$GDF15_group == "High"),
  Events = sum(surv_analysis$OS),
  Logrank_chisq = unname(fit_lr$chisq),
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

sink(file.path(out_dir, "sessionInfo_TCGA.txt"))
sessionInfo()
sink()

cat("\n============================================\n")
cat("FINAL TCGA MANUSCRIPT CHECK\n")
cat("============================================\n")
cat("Expression n =", nrow(tcga), "\n")
print(table(tcga$primary_diagnosis.diagnoses))
cat("GDF15 median =", gdf15_median, "\n")
print(table(tcga$GDF15_group))
cat("\nS8A continuous correlations:\n")
print(spearman_results)
cat("\nS8B purity-adjusted:\n")
print(partial_results)
cat("\nS8C median split:\n")
print(median_split_results)
cat("\nS8D duct-only:\n")
print(duct_results)
cat("\nSupplementary Fig. S5 OS:\n")
print(s5)
cat("Expected log-rank P approximately 0.75.\n")
cat("============================================\n")
