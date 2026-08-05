# ============================================================
# 03_preranked_gsea.R
#
# Pre-ranked GSEA for patient-derived PDAC organoids.
#
# Current manuscript:
#   rank = sign(apeglm-shrunken log2FC) * -log10(nominal P)
#   collections = MSigDB Hallmark + Reactome
#   gene-set size = 15-500
#   random seed = 42
#   prespecified GSEA threshold = q < 0.25
#
# Maps to:
#   Fig. 3b: selected Hallmark signals
#   Fig. 3c: selected Reactome signals
#   Supplementary Table S6:
#     A. Hallmark gene sets with nominal P <0.05
#     B. Reactome pathways with q <0.25
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(readxl)
  library(openxlsx)
})

in_file <- "results/organoid_RNAseq/DESeq2_primary_all_results.xlsx"
out_dir <- "results/organoid_GSEA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(in_file))

deg <- read_excel(in_file) %>%
  filter(
    !is.na(SYMBOL),
    !is.na(log2FoldChange),
    !is.na(pvalue)
  ) %>%
  mutate(
    rank_score = sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))
  ) %>%
  arrange(desc(rank_score), desc(log2FoldChange))

# clusterProfiler requires unique gene names.
# If duplicate symbols exist, retain the row with the largest
# absolute ranking score.
deg <- deg %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(rank_score), desc(log2FoldChange))

gene_list <- setNames(deg$rank_score, deg$SYMBOL)
gene_list <- sort(gene_list, decreasing = TRUE)

get_msig <- function(collection_name) {
  # msigdbr changed argument names across releases.
  # Try the current API first, then fall back to the older API.
  if (collection_name == "H") {
    x <- tryCatch(
      msigdbr(species = "Homo sapiens", collection = "H"),
      error = function(e) NULL
    )
    if (is.null(x)) {
      x <- msigdbr(species = "Homo sapiens", category = "H")
    }
  } else if (collection_name == "REACTOME") {
    x <- tryCatch(
      msigdbr(
        species = "Homo sapiens",
        collection = "C2",
        subcollection = "CP:REACTOME"
      ),
      error = function(e) NULL
    )
    if (is.null(x)) {
      x <- msigdbr(
        species = "Homo sapiens",
        category = "C2",
        subcategory = "CP:REACTOME"
      )
    }
  } else {
    stop("Unknown collection.")
  }

  x %>% select(gs_name, gene_symbol) %>% distinct()
}

hallmark_sets <- get_msig("H")
reactome_sets <- get_msig("REACTOME")

run_gsea <- function(term2gene) {
  set.seed(42)
  GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 15,
    maxGSSize = 500,
    verbose = FALSE,
    seed = TRUE
  )
}

gsea_hallmark <- run_gsea(hallmark_sets)
gsea_reactome <- run_gsea(reactome_sets)

hallmark_full <- as.data.frame(gsea_hallmark)
reactome_full <- as.data.frame(gsea_reactome)

write.xlsx(
  hallmark_full,
  file.path(out_dir, "GSEA_Hallmark_results_full.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  reactome_full,
  file.path(out_dir, "GSEA_Reactome_results_full.xlsx"),
  rowNames = FALSE
)

# Supplementary Table S6A:
# ALL Hallmark gene sets with nominal P<0.05.
s6a <- hallmark_full %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue)

# Supplementary Table S6B:
# ALL Reactome pathways meeting the prespecified q<0.25 threshold.
s6b <- reactome_full %>%
  filter(p.adjust < 0.25) %>%
  arrange(p.adjust, pvalue)

write.xlsx(
  s6a,
  file.path(out_dir, "Supplementary_Table_S6A_Hallmark_nominal_P_lt_0.05.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  s6b,
  file.path(out_dir, "Supplementary_Table_S6B_Reactome_q_lt_0.25.xlsx"),
  rowNames = FALSE
)

# Current Fig. 3b highlights these four Hallmark sets.
fig3b_ids <- c(
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_MYC_TARGETS_V2"
)

fig3b <- hallmark_full %>%
  filter(ID %in% fig3b_ids | Description %in% fig3b_ids)

write.csv(
  fig3b,
  file.path(out_dir, "Fig3b_selected_Hallmark_results.csv"),
  row.names = FALSE
)

cat("\nSupplementary Table S6A: Hallmark nominal P<0.05\n")
print(s6a[, intersect(
  c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust"),
  names(s6a)
)])

cat("\nSupplementary Table S6B: Reactome q<0.25\n")
print(s6b[, intersect(
  c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust"),
  names(s6b)
)])

# Current manuscript-table snapshot includes:
# Hallmark:
#   PI3K-AKT-mTOR signaling  NES ~1.52, P~0.010, q~0.24
#   Protein secretion       NES ~1.41, P~0.034, q~0.40 (nominal only)
#   Interferon-alpha        NES ~1.39, P~0.034, q~0.40 (nominal only)
#   MYC targets V2          NES ~-1.67, P~0.002, q~0.11
#
# Reactome S6B includes all q<0.25 pathways, including
# N-linked glycosylation, ER/Golgi transport and IRE1alpha
# chaperone-related signals.

sink(file.path(out_dir, "sessionInfo_GSEA.txt"))
sessionInfo()
sink()        
