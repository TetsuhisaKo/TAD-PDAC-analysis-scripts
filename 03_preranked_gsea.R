
#!/usr/bin/env Rscript
# ============================================================
# 03_preranked_gsea.R
#
# Outputs (results/organoid_GSEA/):
#   Suppl. Table S8A  Hallmark, BH-adjusted P < 0.25 (8 sets)
#   Suppl. Table S8B  Reactome, BH-adjusted P < 0.05 (16 sets)
#   Fig. 3b           The 16 Reactome pathways of Table S8B
#   Not tabulated:    Reactome, BH-adjusted P < 0.25 (53 sets;
#                     exploratory, saved as a separate file)
#   Also saved: full results, version manifest, sessionInfo
#
# Supplementary Table S8C uses restricted donor-level stage data
# and is not included in this script.
#
# Input: DESeq2 results from 02_deseq2_differential_expression.R
# Ranking: sign(apeglm-shrunken log2FC) * -log10(nominal P)
# Gene-set size: 15-500; BH adjustment within each collection.
# Adjusted P values are from clusterProfiler's BH correction.
#
# Random seed: 42; MSigDB release: 2025.1.Hs.
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(readxl)
  library(openxlsx)
})

# Reporting thresholds and expected gene-set counts.
S8A_Q <- 0.25;          S8A_N <- 8L    # Hallmark
S8B_Q <- 0.05;          S8B_N <- 16L   # Reactome; Fig. 3b
REACTOME_EXPL_Q <- 0.25; REACTOME_EXPL_N <- 53L

# MSigDB release used in the analysis.
REQUIRED_MSIGDB <- "2025.1.Hs"

in_file <- "results/organoid_RNAseq/DESeq2_primary_all_results.xlsx"
out_dir <- "results/organoid_GSEA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop("Missing DESeq2 results: ", in_file)
}

# ------------------------------------------------------------
# Ranked gene list
# ------------------------------------------------------------
deg <- read_excel(in_file) %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != "",
    !is.na(log2FoldChange),
    !is.na(pvalue)
  ) %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    rank_score = sign(log2FoldChange) *
      (-log10(pmax(pvalue, 1e-300)))
  ) %>%
  arrange(
    SYMBOL,
    desc(abs(rank_score)),
    desc(log2FoldChange),
    pvalue
  ) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# Deterministic ordering for tied ranking scores
ord <- order(
  -deg$rank_score,
  -deg$log2FoldChange,
  deg$SYMBOL,
  method = "radix"
)

deg <- deg[ord, , drop = FALSE]
gene_list <- setNames(deg$rank_score, deg$SYMBOL)

stopifnot(
  length(gene_list) > 0,
  !anyDuplicated(names(gene_list)),
  all(diff(gene_list) <= 0)
)

cat("Ranked genes:", length(gene_list),
    "| tied scores:", sum(duplicated(gene_list)), "\n")

# ------------------------------------------------------------
# MSigDB gene sets (msigdbr >= 10 and older interfaces)
# ------------------------------------------------------------
get_msig <- function(type) {
  if (type == "H") {
    x <- tryCatch(
      msigdbr(
        species = "Homo sapiens",
        collection = "H"
      ),
      error = function(e) NULL
    )
    if (is.null(x)) {
      x <- msigdbr(
        species = "Homo sapiens",
        category = "H"
      )
    }
  } else {
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
  }
  x
}

hallmark_raw <- get_msig("H")
reactome_raw <- get_msig("REACTOME")

hallmark_sets <- hallmark_raw %>%
  select(gs_name, gene_symbol) %>%
  distinct()

reactome_sets <- reactome_raw %>%
  select(gs_name, gene_symbol) %>%
  distinct()

db_version <- function(x) {
  if (!"db_version" %in% names(x)) {
    return(paste0(
      "not exposed (msigdbr ",
      packageVersion("msigdbr"), ")"
    ))
  }
  paste(
    sort(unique(na.omit(as.character(x$db_version)))),
    collapse = ";"
  )
}

# Check database release.
for (nm in c("Hallmark", "Reactome")) {
  actual <- db_version(
    list(Hallmark = hallmark_raw, Reactome = reactome_raw)[[nm]]
  )
  if (!identical(actual, REQUIRED_MSIGDB)) {
    stop(
      "MSigDB version mismatch (", nm, "): expected ",
      REQUIRED_MSIGDB, ", obtained ", actual, ".",
      call. = FALSE
    )
  }
}

# ------------------------------------------------------------
# GSEA
# ------------------------------------------------------------
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

# Reported tables carry BH-adjusted P only
report_table <- function(x, q_cut) {
  x %>%
    filter(!is.na(p.adjust), p.adjust < q_cut) %>%
    select(-any_of(c("qvalue", "qvalues")))
}

manifest <- data.frame(
  Item = c(
    "R", "clusterProfiler", "DOSE", "fgsea", "msigdbr",
    "Hallmark_MSigDB", "Reactome_MSigDB",
    "ranking", "ranked_genes",
    "minimum_size", "maximum_size",
    "Hallmark_sets_tested", "Reactome_sets_tested",
    "multiple_testing", "seed_handling"
  ),
  Value = c(
    R.version.string,
    as.character(packageVersion("clusterProfiler")),
    as.character(packageVersion("DOSE")),
    as.character(packageVersion("fgsea")),
    as.character(packageVersion("msigdbr")),
    db_version(hallmark_raw),
    db_version(reactome_raw),
    "sign(shrunken log2FC)*-log10(nominal P)",
    as.character(length(gene_list)),
    "15", "500",
    as.character(nrow(hallmark_full)),
    as.character(nrow(reactome_full)),
    "BH within each collection",
    "set.seed(42); GSEA(seed = TRUE)"
  )
)

write.csv(
  manifest,
  file.path(out_dir, "GSEA_version_and_database_manifest.csv"),
  row.names = FALSE
)

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

# Warn if the number of enriched gene sets differs from the report.
check_count <- function(x, expected, label) {
  if (nrow(x) != expected) {
    warning(
      label, ": expected ", expected, " gene sets, found ",
      nrow(x), ". Confirm MSigDB and package versions.",
      call. = FALSE
    )
  }
}

# S8A: Hallmark, BH-adjusted P < 0.25
s8a <- report_table(hallmark_full, S8A_Q) %>%
  arrange(p.adjust, pvalue)

write.xlsx(
  s8a,
  file.path(
    out_dir,
    "Supplementary_Table_S8A_Hallmark_BH_lt_0.25.xlsx"
  ),
  rowNames = FALSE
)

# S8B: Reactome, BH-adjusted P < 0.05
s8b <- report_table(reactome_full, S8B_Q) %>%
  arrange(p.adjust, pvalue)

write.xlsx(
  s8b,
  file.path(
    out_dir,
    "Supplementary_Table_S8B_Reactome_BH_lt_0.05.xlsx"
  ),
  rowNames = FALSE
)

# Fig. 3b: the Table S8B pathways, ordered by NES
fig3b <- s8b %>%
  arrange(desc(NES))

write.csv(
  fig3b,
  file.path(
    out_dir,
    "Figure_3b_Reactome_BH_lt_0.05.csv"
  ),
  row.names = FALSE
)

# Exploratory Reactome set (BH-adjusted P < 0.25); not
# tabulated in the Supplementary Information
reactome_expl <- report_table(reactome_full, REACTOME_EXPL_Q) %>%
  mutate(In_Table_S8B = ID %in% s8b$ID) %>%
  arrange(p.adjust, pvalue)

write.xlsx(
  reactome_expl,
  file.path(
    out_dir,
    "Reactome_exploratory_BH_lt_0.25.xlsx"
  ),
  rowNames = FALSE
)

check_count(s8a, S8A_N, "Table S8A (Hallmark)")
check_count(s8b, S8B_N, "Table S8B (Reactome)")
check_count(
  reactome_expl, REACTOME_EXPL_N,
  "Reactome, BH-adjusted P < 0.25"
)

cat("GSEA environment:\n")
print(manifest)
cat("\nSupplementary Table S8A:\n")
print(s8a[, c("Description", "NES", "pvalue", "p.adjust")])
cat("\nSupplementary Table S8B / Fig. 3b:\n")
print(fig3b[, c("Description", "NES", "pvalue", "p.adjust")])
cat("\nReactome, BH-adjusted P < 0.25:", nrow(reactome_expl),
    "sets (", sum(!reactome_expl$In_Table_S8B),
    "not in Table S8B)\n")

sink(file.path(out_dir, "sessionInfo_GSEA.txt"))
print(sessionInfo())
sink()

cat(
  "\nSupplementary Table S8C (stage IV-adjusted) is not",
  "implemented: individual-level clinical stage annotations",
  "are not publicly distributed.\n"
)
