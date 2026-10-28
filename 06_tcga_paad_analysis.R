
#!/usr/bin/env Rscript
# ============================================================
# 06_tcga_paad_analysis.R
#
# Outputs (results/TCGA/):
#   Fig. 4c-e          GDF15 and immune-related expression
#                      (values in Supplementary Tables S11A/S11B)
#   Suppl. Table S11A  Median-split analysis
#   Suppl. Table S11B  Spearman correlation
#   Suppl. Table S11C  Purity-adjusted correlation
#   Suppl. Table S11D  Duct-only sensitivity
#   Suppl. Fig. S6     Overall survival (log-rank)
#   Also saved: cohort summary, gene mapping, input MD5,
#   analysis datasets, sessionInfo
#
# Expression: UCSC Xena GDC hub STAR FPKM-UQ matrix, already
#   log2(FPKM-UQ + 1); not transformed again
# GDF15 high: >= median of the 163-tumour cohort
# CYT = mean(GZMA, PRF1) on the supplied log2 scale
# Partial Spearman: rank residualization on ABSOLUTE purity;
#   Fisher-z CI with k = 1 covariate
# BH adjustment across the six transcripts/scores within
#   S11A, S11C and the partial correlations of S11D
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

data_dir <- "data/TCGA"
out_dir <- "results/TCGA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(
  data_dir, "TCGA-PAAD.star_fpkm-uq.tsv.gz"
)
clin_file <- file.path(
  data_dir, "TCGA-PAAD.clinical.tsv.gz"
)
surv_file <- file.path(
  data_dir, "TCGA-PAAD.survival.tsv.gz"
)
probe_file <- file.path(
  data_dir, "gencode.v36.annotation.gtf.gene.probemap"
)
purity_file <- file.path(
  data_dir, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
)

input_files <- c(
  expression = expr_file,
  clinical = clin_file,
  survival = surv_file,
  probe_map = probe_file,
  ABSOLUTE_purity = purity_file
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing inputs: ", paste(missing_files, collapse = ", "))
}

write.csv(
  data.frame(
    File = basename(input_files),
    MD5 = unname(tools::md5sum(input_files))
  ),
  file.path(out_dir, "TCGA_input_file_md5.csv"),
  row.names = FALSE
)

read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  read.delim(
    con,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

# ------------------------------------------------------------
# Expression
# ------------------------------------------------------------
expr_raw <- read_tsv(expr_file)
id_col <- names(expr_raw)[1]

expr_raw <- expr_raw[
  !grepl("_PAR_Y", expr_raw[[id_col]], fixed = TRUE),
  ,
  drop = FALSE
]
expr_raw$Ensembl_ID_clean <- sub(
  "\\..*$", "", expr_raw[[id_col]]
)
expr_raw <- expr_raw[
  !duplicated(expr_raw$Ensembl_ID_clean),
  ,
  drop = FALSE
]

sample_cols <- setdiff(
  names(expr_raw),
  c(id_col, "Ensembl_ID_clean")
)
sample_cols <- sample_cols[grepl("-01A$", sample_cols)]

# GENCODE v36 gene mapping
probe <- read_tsv(probe_file)
stopifnot(all(c("id", "gene") %in% names(probe)))

probe <- probe[
  !grepl("_PAR_Y", probe$id, fixed = TRUE),
  ,
  drop = FALSE
]
probe$Ensembl_ID_clean <- sub("\\..*$", "", probe$id)

genes <- c(
  "GDF15", "CD8A", "GZMA",
  "GZMB", "PRF1", "DDIT3"
)

resolve_gene <- function(symbol) {
  ids <- unique(
    probe$Ensembl_ID_clean[probe$gene == symbol]
  )
  ids <- ids[ids %in% expr_raw$Ensembl_ID_clean]

  if (!length(ids)) {
    stop("Expression gene missing: ", symbol)
  }
  if (length(ids) > 1) {
    warning("Multiple IDs for ", symbol, "; using first.")
  }
  c(ids[1], length(ids))
}

gene_ids <- vapply(genes, resolve_gene, character(2))
gene_map <- data.frame(
  Gene = genes,
  Ensembl_ID = gene_ids[1, ],
  Candidate_IDs = as.integer(gene_ids[2, ])
)

write.csv(
  gene_map,
  file.path(out_dir, "TCGA_gene_mapping_used.csv"),
  row.names = FALSE
)

expr_df <- data.frame(
  sample = sample_cols,
  stringsAsFactors = FALSE
)

for (g in genes) {
  eid <- gene_map$Ensembl_ID[gene_map$Gene == g]

  x <- expr_raw[
    expr_raw$Ensembl_ID_clean == eid,
    sample_cols,
    drop = FALSE
  ]

  if (nrow(x) != 1) {
    stop("Expected one gene row: ", g)
  }

  expr_df[[g]] <- as.numeric(unlist(x[1, ], use.names = FALSE))
}

# Scale guard: values must be on the log2(FPKM-UQ + 1) scale.
# CYT averages log2 values, so a matrix on the linear scale
# would change the CYT ranking.
expr_values <- as.matrix(expr_df[genes])
stopifnot(
  !anyNA(expr_values),
  min(expr_values) >= 0,
  max(expr_values) < 30
)

# ------------------------------------------------------------
# Histology restriction
# ------------------------------------------------------------
clin <- read_tsv(clin_file)

stopifnot(all(
  c("sample", "primary_diagnosis.diagnoses") %in%
    names(clin)
))

allowed_histology <- c(
  "Infiltrating duct carcinoma, NOS",
  "Adenocarcinoma, NOS"
)

clin_pdac <- clin %>%
  filter(
    grepl("-01A$", sample),
    primary_diagnosis.diagnoses %in% allowed_histology
  ) %>%
  select(sample, primary_diagnosis.diagnoses) %>%
  distinct(sample, .keep_all = TRUE)

tcga <- inner_join(
  expr_df, clin_pdac, by = "sample"
) %>%
  mutate(CYT = (GZMA + PRF1) / 2)

stopifnot(
  nrow(tcga) == 163,
  sum(tcga$primary_diagnosis.diagnoses ==
        "Infiltrating duct carcinoma, NOS") == 143,
  sum(tcga$primary_diagnosis.diagnoses ==
        "Adenocarcinoma, NOS") == 20
)

# GDF15 median dichotomization
# (reference: median approximately 4.197)
gdf15_median <- median(tcga$GDF15)

tcga <- tcga %>%
  mutate(
    GDF15_group = factor(
      ifelse(GDF15 >= gdf15_median, "High", "Low"),
      levels = c("Low", "High")
    )
  )

stopifnot(
  sum(tcga$GDF15_group == "Low") == 81,
  sum(tcga$GDF15_group == "High") == 82
)

outcomes <- c(
  "CD8A", "GZMA", "GZMB",
  "PRF1", "CYT", "DDIT3"
)

# ------------------------------------------------------------
# S11A: median-split comparison
# ------------------------------------------------------------
median_split_results <- bind_rows(
  lapply(outcomes, function(v) {
    low <- tcga[[v]][tcga$GDF15_group == "Low"]
    high <- tcga[[v]][tcga$GDF15_group == "High"]

    wt <- wilcox.test(
      low, high,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE
    )

    data.frame(
      Transcript_or_score = v,
      Low_n = sum(!is.na(low)),
      High_n = sum(!is.na(high)),
      Low_median = median(low, na.rm = TRUE),
      High_median = median(high, na.rm = TRUE),
      P_value = wt$p.value
    )
  })
) %>%
  mutate(q_value = p.adjust(P_value, method = "BH"))

write.csv(
  median_split_results,
  file.path(out_dir, "Supplementary_Table_S11A_median_split.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# S11B: continuous Spearman correlations
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
  file.path(
    out_dir,
    "Supplementary_Table_S11B_continuous_Spearman.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# ABSOLUTE purity
# ------------------------------------------------------------
abs_dat <- read_tsv(purity_file)
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

# GDF15 versus purity
keep <- !is.na(tcga$ABSOLUTE_purity)
ct_purity <- cor.test(
  tcga$GDF15[keep],
  tcga$ABSOLUTE_purity[keep],
  method = "spearman",
  exact = FALSE
)

write.csv(
  data.frame(
    n = sum(keep),
    Spearman_rho = unname(ct_purity$estimate),
    P_value = ct_purity$p.value
  ),
  file.path(out_dir, "GDF15_vs_ABSOLUTE_purity.csv"),
  row.names = FALSE
)

# Partial Spearman using ranked residuals
partial_spearman <- function(x, y, z) {
  keep <- complete.cases(x, y, z)

  xr <- rank(x[keep], ties.method = "average")
  yr <- rank(y[keep], ties.method = "average")
  zr <- rank(z[keep], ties.method = "average")

  rx <- resid(lm(xr ~ zr))
  ry <- resid(lm(yr ~ zr))

  rho <- cor(rx, ry)
  n <- length(rx)
  k <- 1L

  if (n <= k + 3) {
    stop("Insufficient observations for partial Spearman.")
  }

  dfree <- n - k - 2
  tval <- rho * sqrt(dfree / (1 - rho^2))
  pval <- 2 * pt(-abs(tval), df = dfree)

  bounded <- min(max(rho, -0.999999), 0.999999)
  zrho <- atanh(bounded)
  se_z <- 1 / sqrt(n - k - 3)

  lower <- tanh(zrho - qnorm(0.975) * se_z)
  upper <- tanh(zrho + qnorm(0.975) * se_z)

  data.frame(
    n = n,
    k_covariates = k,
    Partial_rho = rho,
    CI_low = lower,
    CI_high = upper,
    P_value = pval
  )
}

# ------------------------------------------------------------
# S11C: purity-adjusted correlations
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
  file.path(
    out_dir,
    "Supplementary_Table_S11C_purity_adjusted.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# S11D: histology-restricted sensitivity
# ------------------------------------------------------------
duct <- tcga %>%
  filter(
    primary_diagnosis.diagnoses ==
      "Infiltrating duct carcinoma, NOS"
  )

stopifnot(
  nrow(duct) == 143,
  sum(!is.na(duct$ABSOLUTE_purity)) == 123
)

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
      Spearman_n = sum(
        complete.cases(duct[, c("GDF15", v)])
      ),
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
  mutate(Partial_q = p.adjust(Partial_P, method = "BH"))

write.csv(
  duct_results,
  file.path(
    out_dir,
    "Supplementary_Table_S11D_duct_only_sensitivity.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Fig. S6: overall survival
# ------------------------------------------------------------
surv <- read_tsv(surv_file)
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

stopifnot(
  nrow(surv_analysis) == 162,
  sum(surv_analysis$GDF15_group == "Low") == 81,
  sum(surv_analysis$GDF15_group == "High") == 81
)

fit_lr <- survdiff(
  Surv(OS_months, OS) ~ GDF15_group,
  data = surv_analysis
)

lr_df <- length(fit_lr$n) - 1

os_stats <- data.frame(
  n = nrow(surv_analysis),
  Low_n = sum(surv_analysis$GDF15_group == "Low"),
  High_n = sum(surv_analysis$GDF15_group == "High"),
  Events = sum(surv_analysis$OS),
  Logrank_chisq = unname(fit_lr$chisq),
  Logrank_df = lr_df,
  Logrank_P = pchisq(
    fit_lr$chisq, df = lr_df, lower.tail = FALSE
  )
)

write.csv(
  os_stats,
  file.path(
    out_dir,
    "Supplementary_Figure_S6_OS_statistics.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Cohort summary and analysis datasets
# ------------------------------------------------------------
cohort_summary <- data.frame(
  Item = c(
    "Histology-restricted primary tumours",
    "Infiltrating duct carcinoma, NOS",
    "Adenocarcinoma, NOS",
    "GDF15 median (log2 FPKM-UQ + 1)",
    "GDF15 low / high",
    "ABSOLUTE purity available (all / duct only)",
    "OS evaluable (low / high)",
    "OS events"
  ),
  Value = c(
    nrow(tcga),
    nrow(duct),
    sum(tcga$primary_diagnosis.diagnoses == "Adenocarcinoma, NOS"),
    format(gdf15_median, digits = 6),
    paste(sum(tcga$GDF15_group == "Low"),
          sum(tcga$GDF15_group == "High"), sep = " / "),
    paste(nrow(purity_subset),
          sum(!is.na(duct$ABSOLUTE_purity)), sep = " / "),
    paste(os_stats$Low_n, os_stats$High_n, sep = " / "),
    os_stats$Events
  )
)

write.csv(
  cohort_summary,
  file.path(out_dir, "TCGA_cohort_summary.csv"),
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

cat("\nTCGA-PAAD results:\n")
print(cohort_summary)
print(median_split_results)
print(spearman_results)
print(partial_results)
print(duct_results)
print(os_stats)

sink(file.path(out_dir, "sessionInfo_TCGA.txt"))
print(sessionInfo())
sink()

cat("\nTCGA analysis completed.\n")
