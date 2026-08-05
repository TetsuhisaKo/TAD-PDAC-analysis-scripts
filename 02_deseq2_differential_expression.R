# ============================================================
# 02_deseq2_differential_expression.R
#
# Patient-derived PDAC organoids:
# TAD (n=6) vs non-TAD (n=30)
#
# Current manuscript:
#   Primary DESeq2 design: ~ DM_group
#   Reference: non-TAD
#   Row-summed raw count <=10 excluded
#   apeglm-shrunken log2FC
#   Differential expression: BH-adjusted P<0.05
#   No fold-change threshold
#   PCA: variance-stabilized counts
#
# Sensitivity:
#   ~ platform + DM_group
#   KYK015 and KYK019 = MiSeq
#   all remaining lines = HiSeq 2500
#
# Maps to:
#   Fig. 3a
#   Supplementary Fig. S4
#   Supplementary Table S7A
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(readxl)
  library(dplyr)
  library(tibble)
  library(openxlsx)
})

count_file <- "data/count_matrix.xlsx"
out_dir <- "results/organoid_RNAseq"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(count_file))

counts_df <- read_excel(count_file)
gene_col <- names(counts_df)[1]

counts <- counts_df %>%
  column_to_rownames(var = gene_col)

counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

if (any(is.na(counts))) stop("Count matrix contains NA values.")
if (any(counts < 0)) stop("Count matrix contains negative values.")

# DESeq2 requires integer counts.
if (any(abs(counts - round(counts)) > 1e-8)) {
  warning("Non-integer values detected; rounding to nearest integer for DESeq2.")
}
counts <- round(counts)

tad_cases <- c(
  "KYK019", "KYK020", "KYK067",
  "KYK084", "KYK090", "KYK093"
)

sample_id <- sub("_.*", "", colnames(counts))

meta <- data.frame(
  Sample = colnames(counts),
  case = sample_id,
  DM_group = ifelse(sample_id %in% tad_cases, "TAD", "non-TAD"),
  platform = ifelse(sample_id %in% c("KYK015", "KYK019"), "MiSeq", "HiSeq2500"),
  stringsAsFactors = FALSE
)

meta$DM_group <- factor(meta$DM_group, levels = c("non-TAD", "TAD"))
meta$platform <- factor(meta$platform, levels = c("HiSeq2500", "MiSeq"))
rownames(meta) <- meta$Sample

stopifnot(identical(colnames(counts), rownames(meta)))
stopifnot(ncol(counts) == 36)
stopifnot(sum(meta$DM_group == "TAD") == 6)
stopifnot(sum(meta$DM_group == "non-TAD") == 30)

write.csv(
  meta,
  file.path(out_dir, "organoid_sample_metadata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Helper: preserve unshrunken P/padj but replace log2FC with
# apeglm-shrunken log2FC.
# ------------------------------------------------------------

build_result_table <- function(dds, coef_name) {
  raw_res <- results(dds, name = coef_name, independentFiltering = TRUE)
  shr_res <- lfcShrink(dds, coef = coef_name, type = "apeglm")

  raw_df <- as.data.frame(raw_res) %>%
    rownames_to_column("SYMBOL") %>%
    transmute(
      SYMBOL,
      baseMean,
      raw_log2FoldChange = log2FoldChange,
      lfcSE_raw = lfcSE,
      stat,
      pvalue,
      padj
    )

  shr_df <- as.data.frame(shr_res) %>%
    rownames_to_column("SYMBOL") %>%
    transmute(
      SYMBOL,
      log2FoldChange = log2FoldChange,
      lfcSE_shrunken = lfcSE
    )

  left_join(raw_df, shr_df, by = "SYMBOL") %>%
    arrange(padj, pvalue)
}

# ------------------------------------------------------------
# Primary model
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ DM_group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

coef_primary <- "DM_group_TAD_vs_non.TAD"
if (!(coef_primary %in% resultsNames(dds))) {
  # DESeq2 can sanitize the hyphen differently across versions.
  coef_primary <- grep(
    "^DM_group_.*TAD.*vs.*non",
    resultsNames(dds),
    value = TRUE
  )[1]
}
if (is.na(coef_primary) || length(coef_primary) == 0) {
  stop("Could not identify the TAD vs non-TAD coefficient.")
}

res_primary <- build_result_table(dds, coef_primary)

res_plottable <- res_primary %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

sig_primary <- res_plottable %>%
  filter(padj < 0.05)

cat("Genes after row-sum >10 filter:", nrow(res_primary), "\n")
cat("Genes with non-missing padj:", nrow(res_plottable), "\n")
cat("Significant DEGs (BH-adjusted P<0.05):", nrow(sig_primary), "\n")

# Current dataset snapshot:
#   genes with non-missing padj: approximately 16,272
#   significant DEGs: 307
# If package/input versions differ, inspect before forcing these counts.

write.xlsx(
  res_primary,
  file.path(out_dir, "DESeq2_primary_all_results.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  res_plottable,
  file.path(out_dir, "DESeq2_primary_plottable_results.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  sig_primary,
  file.path(out_dir, "DESeq2_primary_significant_padj_lt_0.05.xlsx"),
  rowNames = FALSE
)

# ------------------------------------------------------------
# Variance-stabilized PCA
# Supplementary Fig. S4
# ------------------------------------------------------------

vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "DM_group", returnData = TRUE)
percent_var <- round(100 * attr(pca, "percentVar"), 1)

pca_out <- pca %>%
  rownames_to_column("Sample") %>%
  mutate(
    PC1_percent = percent_var[1],
    PC2_percent = percent_var[2]
  )

write.csv(
  pca_out,
  file.path(out_dir, "Supplementary_Figure_S4_PCA_coordinates.csv"),
  row.names = FALSE
)

cat(
  "PCA variance: PC1 =", percent_var[1],
  "%; PC2 =", percent_var[2], "%\n"
)
# Current figure: PC1 11.2%, PC2 10.5%.

# ------------------------------------------------------------
# Platform-adjusted sensitivity model
# ------------------------------------------------------------

dds_platform <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ platform + DM_group
)

dds_platform <- dds_platform[rowSums(counts(dds_platform)) > 10, ]
dds_platform <- DESeq(dds_platform)

coef_platform <- "DM_group_TAD_vs_non.TAD"
if (!(coef_platform %in% resultsNames(dds_platform))) {
  coef_platform <- grep(
    "^DM_group_.*TAD.*vs.*non",
    resultsNames(dds_platform),
    value = TRUE
  )[1]
}
if (is.na(coef_platform) || length(coef_platform) == 0) {
  stop("Could not identify the platform-adjusted TAD coefficient.")
}

res_platform <- build_result_table(dds_platform, coef_platform)

write.xlsx(
  res_platform,
  file.path(out_dir, "DESeq2_platform_adjusted_all_results.xlsx"),
  rowNames = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S7A
# Current selected genes
# ------------------------------------------------------------

selected_s7a <- c("GDF15", "DDIT3", "MARS1", "XBP1")

s7a <- res_primary %>%
  select(
    SYMBOL,
    primary_log2FC = log2FoldChange,
    primary_adjusted_P = padj
  ) %>%
  filter(SYMBOL %in% selected_s7a) %>%
  left_join(
    res_platform %>%
      select(
        SYMBOL,
        platform_adjusted_log2FC = log2FoldChange,
        platform_adjusted_P = padj
      ),
    by = "SYMBOL"
  ) %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = selected_s7a)
  ) %>%
  arrange(SYMBOL) %>%
  mutate(SYMBOL = as.character(SYMBOL))

write.csv(
  s7a,
  file.path(out_dir, "Supplementary_Table_S7A_platform_sensitivity.csv"),
  row.names = FALSE
)

cat("\nSupplementary Table S7A selected genes:\n")
print(s7a)

# Expected current values:
# GDF15 primary log2FC 2.30, adjusted P 0.044;
#       platform-adjusted log2FC 2.46, adjusted P 0.025
# DDIT3 primary 3.61, 0.0002; platform 3.57, 0.0003
# MARS1 primary 2.16, 0.0002; platform 2.15, 0.0002
# XBP1  primary 0.95, 0.044;  platform 0.84, 0.062

sink(file.path(out_dir, "sessionInfo_DESeq2.txt"))
sessionInfo()
sink()
