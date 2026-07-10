# ============================================================
# Script 2: Pre-ranked GSEA (Hallmark and Reactome)
#   Ranking metric: sign(log2FC) * -log10(p)
#
# Maps to: Fig. 3b (Hallmark), Fig. 3c (Reactome)
# Note: msigdbr argument names (category/subcategory) are version-dependent;
#       the versions used are recorded via sessionInfo() at the end.
# ============================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(readxl)
  library(openxlsx)
})

deg_all <- read_excel("DESeq2_all_genes.xlsx") %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue))

# ---- Ranking score (signed -log10 p) ----
# Secondary sort on log2FoldChange breaks ties deterministically.
deg_all <- deg_all %>%
  mutate(rank_score = sign(log2FoldChange) * (-log10(pvalue + 1e-300))) %>%
  arrange(desc(rank_score), desc(log2FoldChange))
gene_list <- setNames(deg_all$rank_score, deg_all$SYMBOL)

# ---- Hallmark ----
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
gsea_hallmark <- GSEA(
  geneList = gene_list, TERM2GENE = hallmark_sets,
  pvalueCutoff = 1.0, minGSSize = 15, maxGSSize = 500, seed = 42)
write.xlsx(as.data.frame(gsea_hallmark),
           "GSEA_Hallmark_results_full.xlsx", rowNames = FALSE)

# ---- Reactome ----
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2",
                         subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)
gsea_reactome <- GSEA(
  geneList = gene_list, TERM2GENE = reactome_sets,
  pvalueCutoff = 1.0, minGSSize = 15, maxGSSize = 500, seed = 42)
write.xlsx(as.data.frame(gsea_reactome),
           "GSEA_Reactome_results_full.xlsx", rowNames = FALSE)

# ---- Record session/package versions for reproducibility ----
# Captures clusterProfiler and msigdbr versions (msigdbr API is version-sensitive).
sink("sessionInfo_GSEA.txt"); print(sessionInfo()); sink()           
