# ============================================================
# Script 1: DESeq2 differential expression (TAD vs non-TAD)
#   Primary design: ~ DM_group   (reference = non-TAD / "Normal")
#   Sensitivity:    ~ platform + DM_group
#
# Maps to: Fig. 3a (volcano); platform-adjusted results -> Suppl. Table S6
# Group labels: "DM" = tumor-associated diabetes (TAD); "Normal" = non-TAD.
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(readxl)
  library(dplyr)
  library(tibble)
  library(openxlsx)
})

# ---- Data import ----
# count_matrix.xlsx: genes in the first column, samples in remaining columns.
counts <- read_excel("count_matrix.xlsx") %>%
  column_to_rownames(var = names(.)[1])

# ---- Sample metadata ----
dm_cases <- c("KYK019","KYK020","KYK067","KYK084","KYK090","KYK093")
meta <- data.frame(Sample = colnames(counts)) %>%
  mutate(
    DM_group = ifelse(sub("_.*", "", Sample) %in% dm_cases, "DM", "Normal"),
    # Sequencing platform: KYK015 and KYK019 on MiSeq; all others on HiSeq 2500
    platform = ifelse(sub("_.*", "", Sample) %in% c("KYK015","KYK019"),
                      "MiSeq", "HiSeq")
  )
meta$DM_group <- factor(meta$DM_group, levels = c("Normal", "DM"))
meta$platform <- factor(meta$platform)
rownames(meta) <- meta$Sample

# Ensure count-matrix columns and metadata rows are in the same sample order.
stopifnot(all(colnames(counts) == rownames(meta)))

# ---- Primary model: ~ DM_group ----
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = meta,
  design    = ~ DM_group
)
dds <- dds[rowSums(counts(dds)) > 10, ]      # low-count filter (row-sum > 10)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "DM_group_DM_vs_Normal", type = "apeglm")

res_df <- as.data.frame(res) %>%
  rownames_to_column("SYMBOL") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)
sig_df <- res_df %>% filter(padj < 0.05)
cat("Primary model - significant DEGs (padj<0.05):", nrow(sig_df), "\n")

write.xlsx(res_df, "DESeq2_all_genes.xlsx",        rowNames = FALSE)
write.xlsx(sig_df, "DESeq2_significant_genes.xlsx", rowNames = FALSE)

# ---- Sensitivity model: ~ platform + DM_group ----
dds_p <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = meta,
  design    = ~ platform + DM_group
)
dds_p <- dds_p[rowSums(counts(dds_p)) > 10, ]
dds_p <- DESeq(dds_p)
res_p <- lfcShrink(dds_p, coef = "DM_group_DM_vs_Normal", type = "apeglm")

res_p_df <- as.data.frame(res_p) %>%
  rownames_to_column("SYMBOL") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)
write.xlsx(res_p_df, "DESeq2_platform_adjusted.xlsx", rowNames = FALSE)

# ---- Example genes reported in the text, under each model ----
for (gene in c("GDF15", "DDIT3", "MARS1")) {
  cat(paste0(gene, " (primary):\n"))
  print(res_df[res_df$SYMBOL == gene, ])
  cat(paste0(gene, " (platform-adj):\n"))
  print(res_p_df[res_p_df$SYMBOL == gene, ])
}

# ---- Record session/package versions for reproducibility ----
sink("sessionInfo_DESeq2.txt"); print(sessionInfo()); sink() 
