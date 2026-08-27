#!/usr/bin/env Rscript
# ============================================================
# 03_preranked_gsea.R
#
# Pre-ranked GSEA for patient-derived PDAC organoids.
#
# Computes the statistics underlying:
#   Fig. 3b                    selected Reactome signals
#   Supplementary Table S7A    Hallmark results
#   Supplementary Table S7B    Reactome results
#
# This script outputs enrichment statistics; final figure assembly is performed separately.
#
# Ranking:
#   sign(apeglm-shrunken log2FC) * -log10(nominal P)
#
# Collections: MSigDB Hallmark and Reactome
# Gene-set size: 15-500
# Random seed: 42
# Prespecified significance threshold: FDR q < 0.25
#
# Ties in ranking score are resolved deterministically by
# log2FC and then gene symbol. No later sort() is applied.
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
    SYMBOL != "",
    !is.na(log2FoldChange),
    !is.na(pvalue)
  ) %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    rank_score = sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))
  )

# Unique gene symbols are required. If duplicates exist, retain
# the row with the greatest absolute rank score; tie-breaks are
# explicitly deterministic.
deg <- deg %>%
  arrange(
    SYMBOL,
    desc(abs(rank_score)),
    desc(log2FoldChange),
    pvalue
  ) %>%
  distinct(SYMBOL, .keep_all = TRUE)

ord <- order(
  -deg$rank_score,
  -deg$log2FoldChange,
  deg$SYMBOL,
  method = "radix"
)

deg <- deg[ord, , drop = FALSE]
gene_list <- setNames(deg$rank_score, deg$SYMBOL)

stopifnot(!anyDuplicated(names(gene_list)))
stopifnot(all(diff(gene_list) <= 0))

get_msig_raw <- function(collection_name) {
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
    stop("Unknown collection: ", collection_name)
  }

  x
}

hallmark_raw <- get_msig_raw("H")
reactome_raw <- get_msig_raw("REACTOME")

hallmark_sets <- hallmark_raw %>%
  select(gs_name, gene_symbol) %>%
  distinct()

reactome_sets <- reactome_raw %>%
  select(gs_name, gene_symbol) %>%
  distinct()

db_version_value <- function(x) {
  if ("db_version" %in% names(x)) {
    vals <- sort(unique(na.omit(as.character(x$db_version))))
    if (length(vals) == 0) return("not_available")
    return(paste(vals, collapse = ";"))
  }
  "not_exposed_by_installed_msigdbr"
}

version_manifest <- data.frame(
  Item = c(
    "clusterProfiler_version",
    "msigdbr_version",
    "Hallmark_MSigDB_db_version",
    "Reactome_MSigDB_db_version",
    "ranking_metric",
    "min_gene_set_size",
    "max_gene_set_size",
    "seed"
  ),
  Value = c(
    as.character(packageVersion("clusterProfiler")),
    as.character(packageVersion("msigdbr")),
    db_version_value(hallmark_raw),
    db_version_value(reactome_raw),
    "sign(shrunken log2FC) * -log10(nominal P)",
    "15",
    "500",
    "42"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  version_manifest,
  file.path(out_dir, "GSEA_version_and_database_manifest.csv"),
  row.names = FALSE
)

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

# Supplementary Table S7A: all Hallmark gene sets with nominal P<0.05.
s6a <- hallmark_full %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue, p.adjust)

# Supplementary Table S7B: all Reactome pathways with q<0.25.
s6b <- reactome_full %>%
  filter(p.adjust < 0.25) %>%
  arrange(p.adjust, pvalue)

write.xlsx(
  s6a,
  file.path(out_dir, "Supplementary_Table_S7A_Hallmark_nominal_P_lt_0.05.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  s6b,
  file.path(out_dir, "Supplementary_Table_S7B_Reactome_q_lt_0.25.xlsx"),
  rowNames = FALSE
)

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
  file.path(out_dir, "Fig3b_selected_Hallmark_statistics.csv"),
  row.names = FALSE
)

cat("\nGSEA version/database manifest:\n")
print(version_manifest)

cat("\nSupplementary Table S7A: Hallmark nominal P<0.05\n")
print(
  s6a[, intersect(
    c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust"),
    names(s6a)
  )]
)

cat("\nSupplementary Table S7B: Reactome q<0.25\n")
print(
  s6b[, intersect(
    c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust"),
    names(s6b)
  )]
)

# Reference enrichment results reported in the manuscript:
# PI3K-AKT-mTOR signaling: NES ~1.52, P~0.010, q~0.24
# Protein secretion:       NES ~1.41, P~0.034, q~0.40
# Interferon-alpha:        NES ~1.39, P~0.034, q~0.40
# MYC targets V2:          NES ~-1.67, P~0.002, q~0.11

sink(file.path(out_dir, "sessionInfo_GSEA.txt"))
sessionInfo()
sink()
