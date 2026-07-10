# ============================================================
# Script 4: Driver gene alterations, TAD vs non-TAD (Fisher)
#   Per-gene alteration calls integrate mutation +
#   selected copy-number events (see Methods).
#
# Maps to: Supplementary Fig. S3
# Group labels: "DM" = tumor-associated diabetes (TAD);
#               "Normal" = non-TAD.
#
# Input data:
#   organoid_driver_alterations.xlsx
#     - one row per organoid line
#     - SampleID column
#     - KRAS_any, TP53_any, CDKN2A_any, SMAD4_any columns
#     - each *_any column is a 0/1 alteration call based on
#       the gene-specific definitions described in the Methods.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# ---- Data import ----
big4 <- read_excel("organoid_driver_alterations.xlsx")

# ---- Group assignment ----
dm_cases <- c("KYK019", "KYK020", "KYK067", "KYK084", "KYK090", "KYK093")

big4 <- big4 %>%
  mutate(
    DM_group = ifelse(SampleID %in% dm_cases, "DM", "Normal")
  )

# ---- Fisher's exact test for each driver gene ----
genes <- c("KRAS", "TP53", "CDKN2A", "SMAD4")

results <- do.call(rbind, lapply(genes, function(g) {
  col <- paste0(g, "_any")
  if (!(col %in% colnames(big4))) return(NULL)

  # Force a 2 x 2 table even if one alteration category is absent.
  tbl <- table(
    factor(big4$DM_group, levels = c("Normal", "DM")),
    factor(big4[[col]], levels = c(0, 1))
  )

  ft <- fisher.test(tbl)

  data.frame(
    Gene       = g,
    DM_pct     = mean(big4[[col]][big4$DM_group == "DM"], na.rm = TRUE) * 100,
    Normal_pct = mean(big4[[col]][big4$DM_group == "Normal"], na.rm = TRUE) * 100,
    p.value    = ft$p.value
  )
}))

print(results)

# ---- Record R/package versions ----
sessionInfo()
