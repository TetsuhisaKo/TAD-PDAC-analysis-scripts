
#!/usr/bin/env Rscript
# ============================================================
# 02_deseq2_differential_expression.R
#
# Outputs (results/organoid_RNAseq/):
#   Fig. 3a           Full and plottable DESeq2 results
#   Fig. 3c           ER stress/UPR/Golgi-related genes (log2FC)
#   Fig. 3d           Secreted and secretion-associated genes
#                     (log2FC; TPM sensitivity in
#                     04_secreted_factor_tpm.py, Table S9B)
#   Suppl. Fig. S5    PCA coordinates (VST, top 500 genes)
#   Suppl. Table S9A  Primary vs platform-adjusted DESeq2
#
# Primary DESeq2 design: ~ DM_group
# Sensitivity design:    ~ platform + DM_group
# The stage IV-adjusted model underlying Supplementary
# Table S8C is not implemented: it requires individual-level
# clinical stage annotations that are not publicly distributed.
# Effect sizes: apeglm-shrunken log2FC.
# P values and BH-adjusted P values: unshrunken Wald tests,
# with DESeq2 default independent filtering (alpha = 0.1).
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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(count_file)) {
  stop("Missing controlled-access count matrix: ", count_file)
}

counts_df <- read_excel(count_file)
gene_col <- names(counts_df)[1]
gene_ids <- as.character(counts_df[[gene_col]])

if (anyNA(gene_ids) || any(gene_ids == "")) {
  stop("Missing gene identifiers in count matrix.")
}
if (anyDuplicated(gene_ids)) {
  stop("Duplicate gene identifiers in count matrix.")
}
# Flag identifiers that may have been altered by spreadsheet software.
suspect_ids <- gene_ids[grepl(
  "^[0-9]+(\\.[0-9]+)?$|^[0-9]{1,2}-[A-Za-z]{3}$|^[A-Za-z]{3}-[0-9]{1,2}$",
  gene_ids
)]
if (length(suspect_ids)) {
  warning(
    length(suspect_ids),
    " gene identifiers look numeric or date-like: ",
    paste(head(suspect_ids, 10), collapse = ", ")
  )
}

counts <- counts_df %>%
  column_to_rownames(var = gene_col) %>%
  as.matrix()

storage.mode(counts) <- "numeric"

if (anyNA(counts) || any(counts < 0)) {
  stop("Invalid or missing read counts.")
}
if (any(abs(counts - round(counts)) > 1e-8)) {
  warning("Non-integer counts rounded for DESeq2.")
}
counts <- round(counts)

# ------------------------------------------------------------
# Sample metadata
# ------------------------------------------------------------
tad_cases <- c(
  "KYK019", "KYK020", "KYK067",
  "KYK084", "KYK090", "KYK093"
)
miseq_cases <- c("KYK015", "KYK019")

sample_id <- sub("_.*", "", colnames(counts))

meta <- data.frame(
  Sample = colnames(counts),
  case = sample_id,
  DM_group = ifelse(
    sample_id %in% tad_cases, "TAD", "non-TAD"
  ),
  platform = ifelse(
    sample_id %in% miseq_cases, "MiSeq", "HiSeq2500"
  ),
  stringsAsFactors = FALSE
)

meta$DM_group <- factor(
  meta$DM_group,
  levels = c("non-TAD", "TAD")
)
meta$platform <- factor(
  meta$platform,
  levels = c("HiSeq2500", "MiSeq")
)
rownames(meta) <- meta$Sample

stopifnot(
  ncol(counts) == 36,
  !anyDuplicated(sample_id),
  all(tad_cases %in% sample_id),
  all(miseq_cases %in% sample_id),
  sum(meta$DM_group == "TAD") == 6,
  sum(meta$DM_group == "non-TAD") == 30,
  sum(meta$platform == "MiSeq") == 2,
  identical(colnames(counts), rownames(meta))
)

write.csv(
  meta,
  file.path(out_dir, "organoid_sample_metadata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
# Nominal statistics from unshrunken results;
# effect sizes from apeglm-shrunken log2FC.
build_result_table <- function(dds, coef_name) {
  raw_res <- results(
    dds, name = coef_name,
    independentFiltering = TRUE
  )
  shr_res <- lfcShrink(
    dds, coef = coef_name, type = "apeglm"
  )

  raw_df <- as.data.frame(raw_res) %>%
    rownames_to_column("SYMBOL") %>%
    transmute(
      SYMBOL, baseMean,
      raw_log2FoldChange = log2FoldChange,
      lfcSE_raw = lfcSE,
      stat, pvalue, padj
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

identify_coef <- function(dds) {
  hits <- grep(
    "^DM_group_.*TAD.*vs.*non",
    resultsNames(dds),
    value = TRUE
  )
  if (length(hits) != 1) {
    stop("Could not uniquely identify TAD coefficient.")
  }
  hits[1]
}

run_deseq <- function(design) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = design
  )
  dds <- dds[rowSums(counts(dds)) > 10, ]
  DESeq(dds)
}

# ------------------------------------------------------------
# Primary DESeq2 (Fig. 3a)
# ------------------------------------------------------------
dds <- run_deseq(~ DM_group)

res_primary <- build_result_table(
  dds, identify_coef(dds)
)

res_plottable <- res_primary %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

sig_primary <- res_plottable %>%
  filter(padj < 0.05)

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
  file.path(
    out_dir,
    "DESeq2_primary_significant_padj_lt_0.05.xlsx"
  ),
  rowNames = FALSE
)

cat(
  "Genes with non-missing padj:", nrow(res_plottable), "\n",
  "Significant genes (padj < 0.05):", nrow(sig_primary),
  "(up", sum(sig_primary$raw_log2FoldChange > 0),
  "/ down", sum(sig_primary$raw_log2FoldChange < 0), ")\n"
)

# ------------------------------------------------------------
# Supplementary Fig. S5: PCA
# ------------------------------------------------------------
PCA_NTOP <- 500L
vsd <- vst(dds, blind = TRUE)

pca <- plotPCA(
  vsd,
  intgroup = "DM_group",
  ntop = PCA_NTOP,
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca, "percentVar"), 1
)

pca_out <- pca %>%
  rownames_to_column("Sample") %>%
  mutate(
    PCA_ntop = PCA_NTOP,
    PC1_percent = percent_var[1],
    PC2_percent = percent_var[2]
  )

write.csv(
  pca_out,
  file.path(
    out_dir,
    "Supplementary_Figure_S5_PCA_coordinates.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Platform-adjusted sensitivity
# ------------------------------------------------------------
dds_platform <- run_deseq(~ platform + DM_group)

res_platform <- build_result_table(
  dds_platform, identify_coef(dds_platform)
)

write.xlsx(
  res_platform,
  file.path(
    out_dir, "DESeq2_platform_adjusted_all_results.xlsx"
  ),
  rowNames = FALSE
)

# ------------------------------------------------------------
# Selected genes
# ------------------------------------------------------------
gene_estimates <- function(res, genes) {
  out <- res %>%
    filter(SYMBOL %in% genes) %>%
    mutate(SYMBOL = factor(SYMBOL, levels = genes)) %>%
    arrange(SYMBOL) %>%
    mutate(SYMBOL = as.character(SYMBOL))
  missing <- setdiff(genes, out$SYMBOL)
  if (length(missing)) {
    stop("Genes absent after filtering: ",
         paste(missing, collapse = ", "))
  }
  out
}

panel_table <- function(genes) {
  gene_estimates(res_primary, genes) %>%
    select(
      SYMBOL, baseMean,
      log2FC = log2FoldChange,
      lfcSE = lfcSE_shrunken,
      pvalue, padj
    )
}

# Fig. 3c: ER stress / UPR / Golgi-related genes
fig3c_genes <- c(
  "DDIT3", "OS9", "MARS1", "GOLT1B",
  "COPA", "COPB2", "XBP1"
)

# Fig. 3d: secreted and secretion-associated genes
# (same genes as 04_secreted_factor_tpm.py)
fig3d_genes <- c(
  "GDF15", "PRSS2", "CPZ", "NXPH4",
  "LYZ", "P3H3", "HHIPL1", "SCG2",
  "IL33", "COPA", "BMP4", "FGF19"
)

fig3c <- panel_table(fig3c_genes)
fig3d <- panel_table(fig3d_genes)

write.csv(
  fig3c,
  file.path(out_dir, "Figure_3c_ER_UPR_Golgi_genes.csv"),
  row.names = FALSE
)
write.csv(
  fig3d,
  file.path(out_dir, "Figure_3d_secreted_genes.csv"),
  row.names = FALSE
)

# Supplementary Table S9A: primary vs platform-adjusted model.
# padj columns are BH-adjusted P values within each model.
s9_genes <- c("GDF15", "DDIT3", "MARS1", "XBP1")

s9 <- gene_estimates(res_primary, s9_genes) %>%
  select(
    SYMBOL,
    primary_log2FC = log2FoldChange,
    primary_padj = padj
  ) %>%
  left_join(
    gene_estimates(res_platform, s9_genes) %>%
      select(
        SYMBOL,
        platform_model_log2FC = log2FoldChange,
        platform_model_padj = padj
      ),
    by = "SYMBOL"
  )
stopifnot(!anyNA(s9))

write.csv(
  s9,
  file.path(
    out_dir,
    "Supplementary_Table_S9A_platform_sensitivity.csv"
  ),
  row.names = FALSE
)

print(fig3c)
print(fig3d)
print(s9)

sink(file.path(out_dir, "sessionInfo_DESeq2.txt"))
print(sessionInfo())
sink()

cat("DESeq2 analyses completed.\n")
