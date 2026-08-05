# ============================================================
# TCGA-PAAD analysis:
# GDF15, cytotoxic immune genes, CYT, DDIT3, and survival
#
# TAD-PDAC study
#
# Maps to:
#   Fig. 4c: CD8A, GZMA, GZMB, and PRF1 by GDF15-high/low
#   Fig. 4d: GDF15 vs CD8A
#   Fig. 4e: GDF15 vs DDIT3/CHOP
#   Supplementary Table S8: CYT and DDIT3 analyses
#   Supplementary Fig. S5: overall survival by GDF15-high/low
#
# Public data source:
#   UCSC Xena GDC Hub
#   Accessed: 2026-04-28
#
# Expected expression cohort:
#   n = 163
#   Infiltrating duct carcinoma, NOS = 143
#   Adenocarcinoma, NOS              = 20
#
# Expected GDF15 median split:
#   Low  = 81
#   High = 82
#
# Expected OS-evaluable cohort:
#   n = 162
#   Low  = 81
#   High = 81
#
# IMPORTANT:
# The TCGA-PAAD.star_fpkm-uq.tsv.gz matrix used in this study
# already contains the log2-transformed FPKM-UQ expression
# values used for the analyses.
#
# DO NOT apply an additional log2(x + 1) transformation.
#
# Multiple-testing correction:
# Benjamini-Hochberg correction is applied jointly across the
# six prespecified GDF15 median-split comparisons:
#   CD8A, GZMA, GZMB, PRF1, CYT, and DDIT3.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})


# ============================================================
# 1. Paths
# ============================================================

data_dir    <- "data/TCGA"
results_dir <- "results"

expr_file <- file.path(
  data_dir,
  "TCGA-PAAD.star_fpkm-uq.tsv.gz"
)

survival_file <- file.path(
  data_dir,
  "TCGA-PAAD.survival.tsv.gz"
)

probe_map <- file.path(
  data_dir,
  "gencode.v36.annotation.gtf.gene.probemap"
)

# Change this filename only if the locally downloaded
# UCSC Xena phenotype file has a different name.
clin_file <- file.path(
  data_dir,
  "TCGA-PAAD.GDC_phenotype.tsv.gz"
)

stopifnot(file.exists(expr_file))
stopifnot(file.exists(survival_file))
stopifnot(file.exists(probe_map))
stopifnot(file.exists(clin_file))

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. Record input-file checksums
#
# These checksums identify the exact input files used for the
# analyses and help distinguish future data-version changes
# from changes in the analysis code.
# ============================================================

input_files <- c(
  expression = expr_file,
  survival   = survival_file,
  phenotype  = clin_file,
  probe_map  = probe_map
)

md5_manifest <- data.frame(
  File = basename(input_files),
  MD5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)

write.csv(
  md5_manifest,
  file.path(
    results_dir,
    "TCGA_input_file_md5.csv"
  ),
  row.names = FALSE
)

cat("\nInput-file MD5 checksums:\n")
print(md5_manifest)


# ============================================================
# 3. Load expression matrix
# ============================================================

expr_raw <- read.table(
  gzfile(expr_file),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

id_col <- colnames(expr_raw)[1]

cat(
  "\nRaw expression matrix:",
  nrow(expr_raw),
  "rows x",
  ncol(expr_raw),
  "columns\n"
)


# ============================================================
# 4. Clean Ensembl IDs
#
# PAR_Y entries are removed explicitly before Ensembl version
# suffixes are stripped.
# ============================================================

expr_raw <- expr_raw[
  !grepl(
    "_PAR_Y",
    expr_raw[[id_col]],
    fixed = TRUE
  ),
  ,
  drop = FALSE
]

expr_raw$id_clean <- sub(
  "\\..*$",
  "",
  expr_raw[[id_col]]
)

# Retain one row per cleaned Ensembl ID
expr_raw <- expr_raw[
  !duplicated(expr_raw$id_clean),
  ,
  drop = FALSE
]

rownames(expr_raw) <- expr_raw$id_clean

expr_mat <- expr_raw[
  ,
  !colnames(expr_raw) %in%
    c(id_col, "id_clean"),
  drop = FALSE
]

cat(
  "Cleaned expression matrix:",
  nrow(expr_mat),
  "genes x",
  ncol(expr_mat),
  "samples\n"
)


# ============================================================
# 5. Probe map: Ensembl ID -> gene symbol
# ============================================================

probe <- read.table(
  probe_map,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

stopifnot(
  all(
    c("id", "gene") %in%
      colnames(probe)
  )
)

probe <- probe[
  !grepl(
    "_PAR_Y",
    probe$id,
    fixed = TRUE
  ),
  ,
  drop = FALSE
]

probe$id_clean <- sub(
  "\\..*$",
  "",
  probe$id
)


# ============================================================
# 6. Resolve gene symbols to Ensembl IDs
#
# The exact Ensembl IDs used in the analysis are recorded and
# written to a CSV file.
# ============================================================

genes_needed <- c(
  "GDF15",
  "CD8A",
  "GZMA",
  "GZMB",
  "PRF1",
  "DDIT3"
)

resolve_gene_id <- function(gene) {

  ids <- unique(
    probe$id_clean[
      probe$gene == gene
    ]
  )

  if (length(ids) == 0) {
    stop(
      paste(
        "Gene not found in probe map:",
        gene
      )
    )
  }

  ids_available <- ids[
    ids %in% rownames(expr_mat)
  ]

  if (length(ids_available) == 0) {
    stop(
      paste(
        "No expression-matrix Ensembl ID found for gene:",
        gene
      )
    )
  }

  if (length(ids_available) > 1) {
    warning(
      paste(
        "Multiple available Ensembl IDs found for",
        gene,
        "- using",
        ids_available[1]
      )
    )
  }

  ids_available[1]
}

gene_mapping_used <- data.frame(
  Gene = genes_needed,
  Ensembl_ID = vapply(
    genes_needed,
    resolve_gene_id,
    character(1)
  ),
  stringsAsFactors = FALSE
)

cat("\nGene-to-Ensembl mapping used:\n")
print(gene_mapping_used)

write.csv(
  gene_mapping_used,
  file.path(
    results_dir,
    "TCGA_gene_mapping_used.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. Helper function to extract expression
#
# IMPORTANT:
# The values in the downloaded Xena matrix are used directly.
# No additional log transformation is applied here.
# ============================================================

get_expr <- function(gene, samples) {

  eid <- gene_mapping_used$Ensembl_ID[
    gene_mapping_used$Gene == gene
  ]

  if (length(eid) != 1) {
    stop(
      paste(
        "Unable to resolve a unique Ensembl ID for:",
        gene
      )
    )
  }

  missing_samples <- setdiff(
    samples,
    colnames(expr_mat)
  )

  if (length(missing_samples) > 0) {
    stop(
      paste(
        "Expression samples not found:",
        paste(
          missing_samples,
          collapse = ", "
        )
      )
    )
  }

  vals <- as.numeric(
    expr_mat[
      eid,
      samples,
      drop = TRUE
    ]
  )

  names(vals) <- samples

  vals
}


# ============================================================
# 8. Primary tumor expression samples
# ============================================================

tumor_expr <- colnames(expr_mat)[
  grepl(
    "-01A$",
    colnames(expr_mat)
  )
]

cat(
  "\nPrimary-tumor expression samples available:",
  length(tumor_expr),
  "\n"
)


# ============================================================
# 9. Load clinical phenotype data
# ============================================================

clin <- read.table(
  gzfile(clin_file),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

required_clinical_columns <- c(
  "sample",
  "primary_diagnosis.diagnoses"
)

stopifnot(
  all(
    required_clinical_columns %in%
      colnames(clin)
  )
)


# ============================================================
# 10. Histology-restricted TCGA-PAAD cohort
# ============================================================

pdac_histology <- c(
  "Infiltrating duct carcinoma, NOS",
  "Adenocarcinoma, NOS"
)

clin_pdac <- clin %>%
  filter(
    sample %in% tumor_expr,
    primary_diagnosis.diagnoses %in%
      pdac_histology
  ) %>%
  select(
    sample,
    primary_diagnosis.diagnoses
  ) %>%
  distinct(
    sample,
    .keep_all = TRUE
  )

cat(
  "\nHistology-restricted TCGA-PAAD cohort:\n"
)

print(
  table(
    clin_pdac$
      primary_diagnosis.diagnoses
  )
)

cat(
  "Total n =",
  nrow(clin_pdac),
  "\n"
)

# Expected dataset snapshot
stopifnot(
  nrow(clin_pdac) == 163
)

stopifnot(
  sum(
    clin_pdac$
      primary_diagnosis.diagnoses ==
      "Infiltrating duct carcinoma, NOS"
  ) == 143
)

stopifnot(
  sum(
    clin_pdac$
      primary_diagnosis.diagnoses ==
      "Adenocarcinoma, NOS"
  ) == 20
)

stopifnot(
  !anyDuplicated(
    clin_pdac$sample
  )
)


# ============================================================
# 11. Build expression data frame for n = 163
# ============================================================

samples163 <- clin_pdac$sample

tcga <- data.frame(
  GZMB  = get_expr(
    "GZMB",
    samples163
  ),

  GDF15 = get_expr(
    "GDF15",
    samples163
  ),

  GZMA = get_expr(
    "GZMA",
    samples163
  ),

  CD8A = get_expr(
    "CD8A",
    samples163
  ),

  DDIT3 = get_expr(
    "DDIT3",
    samples163
  ),

  PRF1 = get_expr(
    "PRF1",
    samples163
  ),

  sample = samples163,

  stringsAsFactors = FALSE,
  row.names = NULL
)

tcga <- tcga %>%
  left_join(
    clin_pdac,
    by = "sample"
  )

stopifnot(
  nrow(tcga) == 163
)

stopifnot(
  !anyDuplicated(
    tcga$sample
  )
)


# ============================================================
# 12. Expression-scale sanity checks
#
# These checks are included specifically to detect accidental
# double log transformation.
# ============================================================

cat("\nExpression summaries:\n")

cat("\nGDF15:\n")
print(
  summary(tcga$GDF15)
)

cat(
  "Range:",
  paste(
    range(
      tcga$GDF15,
      na.rm = TRUE
    ),
    collapse = " to "
  ),
  "\n"
)

cat("\nCD8A:\n")
print(
  summary(tcga$CD8A)
)

cat(
  "Range:",
  paste(
    range(
      tcga$CD8A,
      na.rm = TRUE
    ),
    collapse = " to "
  ),
  "\n"
)

# Confirmed in the manuscript analytic dataset:
# CD8A range approximately:
# 0.2028878 to 4.6617610
#
# DO NOT use:
# log2(tcga$CD8A + 1)


# ============================================================
# 13. GDF15 median split
#
# GDF15 grouping is defined using all 163 histology-restricted
# tumors before restricting the cohort for survival analysis.
#
# The observation equal to the median is assigned to the
# GDF15-high group.
# ============================================================

gdf15_median <- median(
  tcga$GDF15,
  na.rm = TRUE
)

tcga$GDF15_group <- factor(
  ifelse(
    tcga$GDF15 >= gdf15_median,
    "High",
    "Low"
  ),
  levels = c(
    "Low",
    "High"
  )
)

cat(
  "\nGDF15 median =",
  gdf15_median,
  "\n"
)

cat(
  "\nGDF15 group counts:\n"
)

print(
  table(
    tcga$GDF15_group
  )
)

# Expected:
# median = 4.19693
# Low  = 81
# High = 82

stopifnot(
  sum(
    tcga$GDF15_group == "Low"
  ) == 81
)

stopifnot(
  sum(
    tcga$GDF15_group == "High"
  ) == 82
)


# ============================================================
# 14. CYT score
#
# CYT is calculated as the arithmetic mean of the
# log2-transformed GZMA and PRF1 expression values, adapted
# from the cytolytic activity score described by Rooney et al.
# ============================================================

tcga$CYT <- (
  tcga$GZMA +
    tcga$PRF1
) / 2


# ============================================================
# 15. Six GDF15 median-split comparisons
#
# Variables:
#   CD8A
#   GZMA
#   GZMB
#   PRF1
#   CYT
#   DDIT3
#
# Statistical test:
#   Two-sided Wilcoxon rank-sum / Mann-Whitney U test
#
# Multiplicity:
#   Benjamini-Hochberg correction across all six comparisons.
# ============================================================

median_split_vars <- c(
  "CD8A",
  "GZMA",
  "GZMB",
  "PRF1",
  "CYT",
  "DDIT3"
)

median_split_results <- do.call(
  rbind,
  lapply(
    median_split_vars,
    function(v) {

      wt <- wilcox.test(
        tcga[[v]] ~
          tcga$GDF15_group,
        exact = FALSE,
        correct = TRUE
      )

      data.frame(
        Variable = v,

        Low_n = sum(
          tcga$GDF15_group ==
            "Low" &
            !is.na(tcga[[v]])
        ),

        High_n = sum(
          tcga$GDF15_group ==
            "High" &
            !is.na(tcga[[v]])
        ),

        Low_median = median(
          tcga[[v]][
            tcga$GDF15_group ==
              "Low"
          ],
          na.rm = TRUE
        ),

        High_median = median(
          tcga[[v]][
            tcga$GDF15_group ==
              "High"
          ],
          na.rm = TRUE
        ),

        P_value = wt$p.value,

        stringsAsFactors = FALSE
      )
    }
  )
)

median_split_results$BH_q <- p.adjust(
  median_split_results$P_value,
  method = "BH"
)

cat(
  "\nSix GDF15 median-split comparisons:\n"
)

print(
  median_split_results,
  digits = 8
)


# ============================================================
# 16. Fig. 4c results:
# CD8A, GZMA, GZMB, and PRF1
# ============================================================

fig4c_results <- median_split_results %>%
  filter(
    Variable %in%
      c(
        "CD8A",
        "GZMA",
        "GZMB",
        "PRF1"
      )
  ) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "CD8A",
        "GZMA",
        "GZMB",
        "PRF1"
      )
    )
  ) %>%
  arrange(
    Variable
  )

cat(
  "\nFig. 4c results:\n"
)

print(
  fig4c_results,
  digits = 8
)


# ============================================================
# 17. CYT median-split result
# Supplementary Table S8
# ============================================================

cyt_group_result <- median_split_results %>%
  filter(
    Variable == "CYT"
  )

cat(
  "\nCYT by GDF15 group:\n"
)

print(
  cyt_group_result,
  digits = 8
)


# ============================================================
# 18. DDIT3 median-split result
# Supplementary Table S8
# ============================================================

ddit3_group_result <- median_split_results %>%
  filter(
    Variable == "DDIT3"
  )

cat(
  "\nDDIT3 by GDF15 group:\n"
)

print(
  ddit3_group_result,
  digits = 8
)


# ============================================================
# 19. GDF15 vs CYT
# Spearman correlation
# Supplementary Table S8
# ============================================================

rho_cyt <- cor.test(
  tcga$GDF15,
  tcga$CYT,
  method = "spearman",
  exact = FALSE
)

cyt_cor_result <- data.frame(
  Analysis = "GDF15_vs_CYT",
  N = sum(
    complete.cases(
      tcga[, c(
        "GDF15",
        "CYT"
      )]
    )
  ),
  Spearman_rho = unname(
    rho_cyt$estimate
  ),
  P_value = rho_cyt$p.value,
  stringsAsFactors = FALSE
)

cat(
  "\nGDF15 vs CYT:\n"
)

print(
  cyt_cor_result,
  digits = 8
)


# ============================================================
# 20. Fig. 4d:
# GDF15 vs CD8A
# Spearman correlation
# ============================================================

rho_cd8a <- cor.test(
  tcga$GDF15,
  tcga$CD8A,
  method = "spearman",
  exact = FALSE
)

fig4d_result <- data.frame(
  Analysis = "Fig4d",
  Comparison = "GDF15_vs_CD8A",
  N = sum(
    complete.cases(
      tcga[, c(
        "GDF15",
        "CD8A"
      )]
    )
  ),
  Spearman_rho = unname(
    rho_cd8a$estimate
  ),
  P_value = rho_cd8a$p.value,
  stringsAsFactors = FALSE
)

cat(
  "\nFig. 4d result:\n"
)

print(
  fig4d_result,
  digits = 8
)

# Confirmed using the current n = 163 dataset:
# rho = -0.3705219
# P   = 1.124e-06
#
# Reported:
# rho = -0.37, P < 0.001


# ============================================================
# 21. Fig. 4e:
# GDF15 vs DDIT3 / CHOP
# Spearman correlation
# ============================================================

rho_ddit3 <- cor.test(
  tcga$GDF15,
  tcga$DDIT3,
  method = "spearman",
  exact = FALSE
)

fig4e_result <- data.frame(
  Analysis = "Fig4e",
  Comparison = "GDF15_vs_DDIT3_CHOP",
  N = sum(
    complete.cases(
      tcga[, c(
        "GDF15",
        "DDIT3"
      )]
    )
  ),
  Spearman_rho = unname(
    rho_ddit3$estimate
  ),
  P_value = rho_ddit3$p.value,
  stringsAsFactors = FALSE
)

cat(
  "\nFig. 4e result:\n"
)

print(
  fig4e_result,
  digits = 8
)


# ============================================================
# 22. Load survival data
# ============================================================

surv <- read.table(
  gzfile(survival_file),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

required_survival_columns <- c(
  "sample",
  "OS.time",
  "OS"
)

stopifnot(
  all(
    required_survival_columns %in%
      colnames(surv)
  )
)

surv <- surv[
  grepl(
    "-01A$",
    surv$sample
  ),
  ,
  drop = FALSE
]

surv$patient <- substr(
  surv$sample,
  1,
  12
)

tcga$patient <- substr(
  tcga$sample,
  1,
  12
)

stopifnot(
  !anyDuplicated(
    tcga$patient
  )
)


# ============================================================
# 23. Build OS-evaluable cohort
#
# IMPORTANT:
# The GDF15 group defined in the full n = 163 expression cohort
# is retained.
#
# The GDF15 median is NOT recalculated in the survival subset.
# ============================================================

surv_analysis <- tcga %>%
  left_join(
    surv %>%
      select(
        patient,
        OS.time,
        OS
      ) %>%
      distinct(
        patient,
        .keep_all = TRUE
      ),
    by = "patient"
  ) %>%
  filter(
    !is.na(OS.time),
    !is.na(OS)
  )

cat(
  "\nOS-evaluable cohort: n =",
  nrow(surv_analysis),
  "| events =",
  sum(
    surv_analysis$OS,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nGDF15 groups in OS-evaluable cohort:\n"
)

print(
  table(
    surv_analysis$GDF15_group
  )
)

# Expected:
# n = 162
# Low  = 81
# High = 81

stopifnot(
  nrow(surv_analysis) == 162
)

stopifnot(
  sum(
    surv_analysis$GDF15_group ==
      "Low"
  ) == 81
)

stopifnot(
  sum(
    surv_analysis$GDF15_group ==
      "High"
  ) == 81
)


# ============================================================
# 24. Convert OS time to months
# ============================================================

surv_analysis$OS_months <-
  surv_analysis$OS.time / 30.44


# ============================================================
# 25. Supplementary Fig. S5:
# Overall survival according to GDF15 expression
# ============================================================

lr <- survdiff(
  Surv(
    OS_months,
    OS
  ) ~ GDF15_group,
  data = surv_analysis
)

logrank_p <- 1 -
  pchisq(
    lr$chisq,
    df = 1
  )

s5_result <- data.frame(
  Analysis = "Supplementary_Fig_S5",
  Test = "Log-rank",
  N = nrow(
    surv_analysis
  ),
  Low_n = sum(
    surv_analysis$GDF15_group ==
      "Low"
  ),
  High_n = sum(
    surv_analysis$GDF15_group ==
      "High"
  ),
  Events = sum(
    surv_analysis$OS,
    na.rm = TRUE
  ),
  Chi_square = unname(
    lr$chisq
  ),
  P_value = logrank_p,
  stringsAsFactors = FALSE
)

cat(
  "\nSupplementary Fig. S5 result:\n"
)

print(
  s5_result,
  digits = 8
)

# Expected current result:
# Low n  = 81
# High n = 81
# log-rank P approximately 0.75


# ============================================================
# 26. Save statistical results
# ============================================================

write.csv(
  median_split_results,
  file.path(
    results_dir,
    "TCGA_GDF15_six_median_split_comparisons.csv"
  ),
  row.names = FALSE
)

write.csv(
  fig4c_results,
  file.path(
    results_dir,
    "TCGA_Fig4c_effector_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  cyt_group_result,
  file.path(
    results_dir,
    "TCGA_CYT_group_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  ddit3_group_result,
  file.path(
    results_dir,
    "TCGA_DDIT3_group_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  cyt_cor_result,
  file.path(
    results_dir,
    "TCGA_GDF15_CYT_correlation.csv"
  ),
  row.names = FALSE
)

write.csv(
  fig4d_result,
  file.path(
    results_dir,
    "TCGA_Fig4d_GDF15_CD8A.csv"
  ),
  row.names = FALSE
)

write.csv(
  fig4e_result,
  file.path(
    results_dir,
    "TCGA_Fig4e_GDF15_DDIT3.csv"
  ),
  row.names = FALSE
)

write.csv(
  s5_result,
  file.path(
    results_dir,
    "TCGA_SuppFigS5_survival.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 27. Save exact analytic datasets
#
# These files contain only public TCGA-derived data.
# ============================================================

write.csv(
  tcga,
  file.path(
    results_dir,
    "TCGA_PAAD_histology_restricted_n163.csv"
  ),
  row.names = FALSE
)

write.csv(
  surv_analysis,
  file.path(
    results_dir,
    "TCGA_PAAD_survival_evaluable_n162.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 28. Save R and package versions
# ============================================================

sink(
  file.path(
    results_dir,
    "sessionInfo_TCGA_PAAD_GDF15_analysis.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 29. Final validation summary
# ============================================================

cat("\n")
cat("============================================\n")
cat("FINAL VALIDATION SUMMARY\n")
cat("============================================\n")

cat(
  "\nExpression cohort n =",
  nrow(tcga),
  "\n"
)

cat(
  "\nHistology:\n"
)

print(
  table(
    tcga$
      primary_diagnosis.diagnoses
  )
)

cat(
  "\nGDF15 median =",
  gdf15_median,
  "\n"
)

cat(
  "\nGDF15 groups:\n"
)

print(
  table(
    tcga$GDF15_group
  )
)

cat(
  "\nCD8A range =",
  paste(
    range(
      tcga$CD8A,
      na.rm = TRUE
    ),
    collapse = " to "
  ),
  "\n"
)

cat(
  "\nSix median-split comparisons:\n"
)

print(
  median_split_results,
  digits = 8
)

cat(
  "\nGDF15 vs CYT:\n"
)

print(
  cyt_cor_result,
  digits = 8
)

cat(
  "\nFig. 4d GDF15 vs CD8A:\n"
)

print(
  fig4d_result,
  digits = 8
)

cat(
  "\nFig. 4e GDF15 vs DDIT3:\n"
)

print(
  fig4e_result,
  digits = 8
)

cat(
  "\nOS-evaluable cohort n =",
  nrow(surv_analysis),
  "\n"
)

cat(
  "\nOS GDF15 groups:\n"
)

print(
  table(
    surv_analysis$GDF15_group
  )
)

cat(
  "\nSupplementary Fig. S5:\n"
)

print(
  s5_result,
  digits = 8
)

cat("\nInput-file MD5 checksums:\n")
print(md5_manifest)

cat("\nGene-to-Ensembl mapping used:\n")
print(gene_mapping_used)

cat("\n============================================\n")
cat("Analysis completed successfully.\n")
cat("============================================\n")
