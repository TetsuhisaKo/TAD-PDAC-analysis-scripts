# ============================================================
# 08_driver_gene_wgs_analysis.R
#
# Major PDAC driver-gene alterations in patient-derived
# organoids according to TAD status.
#
# Current manuscript:
#   WGS-derived organoid alteration calls are compared using
#   Fisher's exact test.
#
# Maps to Supplementary Fig. S3.
#
# Expected groups:
#   non-TAD n=30
#   TAD n=6
#
# Expected current figure:
#   KRAS    100% vs 100%, P=1.000
#   TP53     87% vs  67%, P=0.256
#   CDKN2A   77% vs  50%, P=0.317
#   SMAD4    70% vs  33%, P=0.161
#
# Input file contains PRE-DERIVED gene-level 0/1 alteration
# calls from WGS. This script does not process raw WGS data.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

in_file <- "data/organoid_driver_alterations.xlsx"
out_dir <- "results/organoid_driver_genes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(in_file))

big4 <- read_excel(in_file)

required <- c(
  "SampleID",
  "KRAS_any", "TP53_any", "CDKN2A_any", "SMAD4_any"
)

missing <- setdiff(required, names(big4))
if (length(missing) > 0) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

tad_cases <- c(
  "KYK019", "KYK020", "KYK067",
  "KYK084", "KYK090", "KYK093"
)

big4 <- big4 %>%
  mutate(
    DM_group = factor(
      ifelse(SampleID %in% tad_cases, "TAD", "non-TAD"),
      levels = c("non-TAD", "TAD")
    )
  )

stopifnot(nrow(big4) == 36)
stopifnot(sum(big4$DM_group == "TAD") == 6)
stopifnot(sum(big4$DM_group == "non-TAD") == 30)

genes <- c("KRAS", "TP53", "CDKN2A", "SMAD4")

results <- bind_rows(
  lapply(genes, function(g) {
    col <- paste0(g, "_any")

    if (!all(na.omit(big4[[col]]) %in% c(0, 1))) {
      stop(col, " must contain only 0/1 (or NA) calls.")
    }

    tab <- table(
      factor(big4$DM_group, levels = c("non-TAD", "TAD")),
      factor(big4[[col]], levels = c(0, 1))
    )

    ft <- fisher.test(tab)

    data.frame(
      Gene = g,
      Non_TAD_n = sum(big4$DM_group == "non-TAD" & !is.na(big4[[col]])),
      TAD_n = sum(big4$DM_group == "TAD" & !is.na(big4[[col]])),
      Non_TAD_altered_n = sum(
        big4$DM_group == "non-TAD" & big4[[col]] == 1,
        na.rm = TRUE
      ),
      TAD_altered_n = sum(
        big4$DM_group == "TAD" & big4[[col]] == 1,
        na.rm = TRUE
      ),
      Non_TAD_percent = mean(
        big4[[col]][big4$DM_group == "non-TAD"],
        na.rm = TRUE
      ) * 100,
      TAD_percent = mean(
        big4[[col]][big4$DM_group == "TAD"],
        na.rm = TRUE
      ) * 100,
      P_value = ft$p.value
    )
  })
)

print(results)

write.csv(
  results,
  file.path(out_dir, "Supplementary_Figure_S3_driver_gene_statistics.csv"),
  row.names = FALSE
)

sink(file.path(out_dir, "sessionInfo_driver_genes.txt"))
sessionInfo()
sink()
sessionInfo()
