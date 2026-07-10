#!/usr/bin/env Rscript
# ============================================================
# RNAseq-5: ER stress/UPR leave-one-out overlap sensitivity analysis
#
# Purpose:
#   Assess whether over-representation of a curated ER stress/UPR
#   gene set among genes upregulated in TAD-associated organoids is
#   preserved after excluding individual TAD-derived organoid lines.
#   For each iteration one TAD-derived line is removed, DESeq2 is
#   refitted (apeglm shrinkage), and hypergeometric over-representation
#   of the ER/UPR gene set among upregulated genes is tested.
#
# Maps to:
#   Supplementary Table S7A: nominally upregulated genes
#                            (log2FoldChange > 0 and nominal P < 0.05)
#   Supplementary Table S7B: adjusted-P-based upregulated genes
#                            (log2FoldChange > 0 and adjusted P < 0.05)
#
# Manuscript interpretation:
#   ER stress/UPR-related enrichment was directionally preserved when
#   nominally upregulated genes were used (S7A), whereas enrichment was
#   not retained under the stricter adjusted-P-based criterion after
#   exclusion of individual TAD-derived lines (S7B). This supports a
#   directionally preserved but exploratory ER stress/UPR-related signal.
#
# Expected full-cohort value (S7A, nominal criterion):
#   ~1041 upregulated genes, 17 observed / 4.44 expected,
#   fold over expected 3.83, hypergeometric P = 1.01e-06.
#
# Gene universe:
# Restrict the count matrix to genes included in the primary DESeq2 analysis.
# In each leave-one-out iteration, the gene universe comprises all genes
# tested after applying the same low-count filter.
#
# ER/UPR gene set:
#   A fixed, manually curated ER stress / UPR / ERAD / ER-Golgi-transport gene set is defined in this script
#   It is embedded here so the manuscript values reproduce exactly and
#   are not dependent on an external database version.
#
# Group labels:
#   DM_group == 1  -> "TAD" (tumor-associated diabetes)
#   DM_group == 0  -> "nonTAD"
#
# Input data (NOT distributed here; controlled access per the manuscript
# Data Availability statement). Place in data/ :
#   data/count_analysis.csv              (gene_id + one count column per line)
#   data/organoid_DM_status.xlsx         (columns: case, DM_group)
#   data/DESeq2_all_results.xlsx         (main DESeq2 universe + annotation:
#                                         gene_id, SYMBOL, ENTREZID, GENENAME)
#
# Usage:
#   Rscript RNAseq_Script5_ER_UPR_leave_one_out_overlap.R
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(readxl)
})

# ---- Paths (relative; edit to match local layout) ----
base_dir   <- "data"
count_file <- file.path(base_dir, "count_analysis.csv")
meta_file  <- file.path(base_dir, "organoid_DM_status.xlsx")
main_file  <- file.path(base_dir, "DESeq2_all_results.xlsx")
out_dir    <- "results/leave_one_out_ER_UPR"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_files <- c(count_file, meta_file, main_file)
missing_files  <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

# ---- Genes of interest (reported per-gene in the leave-one-out output) ----
genes_of_interest <- c(
  "GDF15", "DDIT3", "XBP1", "MARS1", "EIF2AK3", "ATF4", "ATF3",
  "HSPA5", "HSP90B1", "HERPUD1", "DNAJB9", "DNAJC3", "EDEM3",
  "SLC7A11", "IL33", "AREG", "EREG", "SERPINE1", "CXCL8", "CXCL2"
)

# ---- Fixed curated ER stress / UPR / ERAD / ER-Golgi gene set ----
# This exact list defines the ER/UPR target set for the over-representation
# test and reproduces the manuscript S7A/S7B values.
er_genes <- c(
  "ATF3", "ATF4", "ATF6", "ATF6B", "DDIT3", "PPP1R15A", "HSPA5",
  "HSP90B1", "HYOU1", "XBP1", "ERN1", "EIF2AK3", "DNAJB9", "DNAJC3",
  "HERPUD1", "EDEM1", "EDEM2", "EDEM3", "SEL1L", "DERL1", "DERL2",
  "SYVN1", "PDIA3", "PDIA4", "PDIA5", "PDIA6", "P4HB", "CALR", "CANX",
  "MANF", "CREB3L1", "CREB3L2", "MBTPS1", "MBTPS2", "SEC61A1",
  "SEC61B", "SEC61G", "SEC62", "SEC63", "SSR1", "SSR2", "SSR3",
  "SSR4", "TMED2", "TMED10", "COPA", "COPB1", "COPB2", "COPG1",
  "COPG2", "ARCN1", "SAR1A", "SAR1B", "MARS1", "WARS1", "ASNS",
  "TRIB3", "SLC7A5", "SLC7A11", "CEBPB", "GADD45A", "DNAJA4",
  "DDIT4", "EIF4EBP1"
)

sample_id_from_col <- function(x) sub("^(KYK[0-9]+).*", "\\1", x)
rank_metric <- function(log2fc, pvalue) sign(log2fc) * -log10(pvalue)

# ---- Load main DESeq2 result table (defines the shared gene universe) ----
main_results <- as.data.frame(read_excel(main_file))
main_results <- main_results[!duplicated(main_results$gene_id), ]
main_gene_ids <- main_results$gene_id
annotation <- unique(main_results[, c("gene_id", "SYMBOL", "ENTREZID", "GENENAME")])

# ---- Load counts and restrict to the main gene universe ----
counts_raw <- read.csv(count_file, check.names = FALSE)
stopifnot("gene_id" %in% colnames(counts_raw))
counts <- counts_raw[counts_raw$gene_id %in% main_gene_ids, , drop = FALSE]
counts <- counts[match(main_gene_ids, counts$gene_id), , drop = FALSE]
if (any(is.na(counts$gene_id))) {
  stop("Count matrix is missing genes from the main DESeq2 universe.")
}
rownames(counts) <- counts$gene_id

count_cols <- setdiff(colnames(counts), "gene_id")
sample_ids <- sample_id_from_col(count_cols)
colnames(counts)[match(count_cols, colnames(counts))] <- sample_ids
count_mat <- as.matrix(counts[, sample_ids, drop = FALSE])
storage.mode(count_mat) <- "integer"

# ---- Metadata ----
metadata <- as.data.frame(read_excel(meta_file))
metadata$case <- as.character(metadata$case)
metadata <- metadata[match(sample_ids, metadata$case), , drop = FALSE]
if (any(is.na(metadata$case))) {
  stop("Some count-matrix sample names were not found in metadata.")
}
metadata$DM_group <- factor(
  ifelse(metadata$DM_group == 1, "TAD", "nonTAD"),
  levels = c("nonTAD", "TAD")
)
rownames(metadata) <- sample_ids
message("nonTAD n: ", sum(metadata$DM_group == "nonTAD"),
        "; TAD n: ", sum(metadata$DM_group == "TAD"))

# ============================================================
# DESeq2 + ER/UPR over-representation for one (sub)cohort
# ============================================================
run_deseq <- function(count_matrix, coldata, label) {
  message("Running ", label)
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData   = coldata,
    design    = ~ DM_group
  )
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)

  res_unshrunken <- results(dds, contrast = c("DM_group", "TAD", "nonTAD"))
  coef_name <- "DM_group_TAD_vs_nonTAD"
  stopifnot(coef_name %in% resultsNames(dds))
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

  res_df <- as.data.frame(res_shrunk)
  un_df  <- as.data.frame(res_unshrunken)
  res_df$stat   <- un_df$stat
  res_df$pvalue <- un_df$pvalue
  res_df$padj   <- un_df$padj
  res_df$gene_id <- rownames(res_df)
  res_df$rank_metric <- rank_metric(res_df$log2FoldChange, res_df$pvalue)
  res_df <- merge(res_df, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
  res_df <- res_df[order(res_df$pvalue), ]
  res_df$rank_pvalue <- seq_len(nrow(res_df))
  write.csv(res_df, file.path(out_dir, paste0("DESeq2_", label, ".csv")), row.names = FALSE)

  universe <- res_df[!is.na(res_df$pvalue) & !is.na(res_df$SYMBOL), ]
  universe$is_er <- universe$SYMBOL %in% er_genes

  summarize_er <- function(selected, name) {
    n_universe   <- nrow(universe)
    k_er         <- sum(universe$is_er)
    n_selected   <- sum(selected, na.rm = TRUE)
    a_er_selected <- sum(universe$is_er & selected, na.rm = TRUE)
    expected     <- n_selected * k_er / n_universe
    p_hyper <- phyper(a_er_selected - 1, k_er, n_universe - k_er,
                      n_selected, lower.tail = FALSE)
    data.frame(
      analysis = label, set = name,
      universe_genes = n_universe, er_genes_in_universe = k_er,
      selected_genes = n_selected, er_genes_selected = a_er_selected,
      expected_er_genes = expected,
      fold_over_expected = ifelse(expected > 0, a_er_selected / expected, NA_real_),
      hypergeom_p = p_hyper
    )
  }

  er_summary <- rbind(
    summarize_er(universe$pvalue < 0.05 & universe$log2FoldChange > 0, "S7A_p_lt_0.05_up"),
    summarize_er(universe$padj   < 0.05 & universe$log2FoldChange > 0, "S7B_padj_lt_0.05_up")
  )

  goi <- res_df[res_df$SYMBOL %in% genes_of_interest,
                c("gene_id","SYMBOL","GENENAME","baseMean","log2FoldChange",
                  "lfcSE","stat","pvalue","padj","rank_pvalue","rank_metric")]
  goi$analysis <- label
  list(goi = goi, er_summary = er_summary)
}

# ============================================================
# Full cohort + leave-one-out over each TAD-derived line
# ============================================================
tad_samples <- rownames(metadata)[metadata$DM_group == "TAD"]
all_goi <- list(); all_er <- list()

full <- run_deseq(count_mat, metadata, "full")
all_goi[["full"]] <- full$goi
all_er[["full"]]  <- full$er_summary

for (s in tad_samples) {
  label <- paste0("leave_out_", s)
  keep  <- setdiff(colnames(count_mat), s)
  tmp   <- run_deseq(count_mat[, keep, drop = FALSE], metadata[keep, , drop = FALSE], label)
  all_goi[[label]] <- tmp$goi
  all_er[[label]]  <- tmp$er_summary
}

goi_summary <- do.call(rbind, all_goi)
er_summary  <- do.call(rbind, all_er)

write.csv(goi_summary, file.path(out_dir, "summary_genes_of_interest.csv"), row.names = FALSE)
write.csv(er_summary,  file.path(out_dir, "summary_ER_UPR_overlap.csv"),     row.names = FALSE)

cat("\nER/UPR over-representation summary (S7A = nominal, S7B = adjustedP):\n")
print(er_summary)

# ---- Record R/package versions ----
sink(file.path(out_dir, "sessionInfo_RNAseq_ER_UPR_leave_one_out_overlap.txt"))
sessionInfo()
sink()
