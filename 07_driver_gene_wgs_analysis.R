
#!/usr/bin/env Rscript
# ============================================================
# 07_driver_gene_wgs_analysis.R
#
# Supplementary Fig. S4:
# KRAS / TP53 / CDKN2A / SMAD4 alterations
# Patient-derived organoids: TAD n=6, non-TAD n=30
#
# Fisher exact tests (two-sided) without multiplicity
# adjustment. Input consists of pre-derived WGS alteration
# calls (0/1 per gene, "_any" = any qualifying event as
# defined in the Supplementary Fig. S4 legend):
#   KRAS    mutation
#   TP53    mutation, heterozygous deletion or CN-LOH/UPD
#   CDKN2A  mutation or homozygous deletion
#   SMAD4   mutation, homozygous or heterozygous deletion
# This script does not process raw WGS reads.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

in_file <- "data/organoid_driver_alterations.xlsx"
meta_file <- "results/organoid_RNAseq/organoid_sample_metadata.csv"
out_dir <- "results/organoid_driver_genes"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(in_file)) {
  stop("Missing controlled-access WGS calls: ", in_file)
}

big4 <- read_excel(in_file)

genes <- c("KRAS", "TP53", "CDKN2A", "SMAD4")
call_cols <- paste0(genes, "_any")

required <- c("SampleID", call_cols)

missing <- setdiff(required, names(big4))
if (length(missing)) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# Validate values as binary 0/1 before converting to integer.
# This avoids silently truncating invalid values such as 0.5.
for (col in call_cols) {
  x <- big4[[col]]
  if (is.character(x)) {
    x <- trimws(toupper(x))
    x[x %in% c("TRUE", "YES")] <- "1"
    x[x %in% c("FALSE", "NO")] <- "0"
  }

  x_num <- suppressWarnings(as.numeric(x))

  if (
    anyNA(x_num) ||
    any(!is.finite(x_num)) ||
    !all(x_num %in% c(0, 1))
  ) {
    stop(
      "Alteration calls must be non-missing binary ",
      "0/1 values for ", col
    )
  }

  big4[[col]] <- as.integer(x_num)
}

tad_cases <- c(
  "KYK019", "KYK020", "KYK067",
  "KYK084", "KYK090", "KYK093"
)

big4 <- big4 %>%
  mutate(
    SampleID = trimws(as.character(SampleID)),
    DM_group = factor(
      ifelse(SampleID %in% tad_cases, "TAD", "non-TAD"),
      levels = c("non-TAD", "TAD")
    )
  )

stopifnot(
  nrow(big4) == 36,
  !anyDuplicated(big4$SampleID),
  all(tad_cases %in% big4$SampleID),
  sum(big4$DM_group == "TAD") == 6,
  sum(big4$DM_group == "non-TAD") == 30
)

# Cross-check organoids against the RNA-seq sample metadata
if (file.exists(meta_file)) {
  meta <- read.csv(meta_file, stringsAsFactors = FALSE)
  if (!setequal(meta$case, big4$SampleID)) {
    stop("Organoids differ from the RNA-seq sample metadata.")
  }
} else {
  warning(meta_file, " not found; sample cross-check skipped.")
}

results <- bind_rows(
  lapply(genes, function(g) {
    col <- paste0(g, "_any")

    tab <- table(
      big4$DM_group,
      factor(big4[[col]], levels = c(0, 1))
    )
    stopifnot(sum(tab) == sum(!is.na(big4[[col]])))

    ft <- fisher.test(tab)

    non <- big4[[col]][big4$DM_group == "non-TAD"]
    tad <- big4[[col]][big4$DM_group == "TAD"]

    data.frame(
      Gene = g,
      Non_TAD_n = sum(!is.na(non)),
      TAD_n = sum(!is.na(tad)),
      Non_TAD_altered_n = sum(non == 1, na.rm = TRUE),
      TAD_altered_n = sum(tad == 1, na.rm = TRUE),
      Non_TAD_percent = mean(non, na.rm = TRUE) * 100,
      TAD_percent = mean(tad, na.rm = TRUE) * 100,
      P_value = ft$p.value,
      Multiple_testing_adjustment = "None; exploratory",
      stringsAsFactors = FALSE
    )
  })
)

print(results)

write.csv(
  results,
  file.path(
    out_dir,
    "Supplementary_Figure_S4_driver_gene_statistics.csv"
  ),
  row.names = FALSE
)

sink(file.path(out_dir, "sessionInfo_driver_genes.txt"))
print(sessionInfo())
sink()

cat("\nDriver-gene WGS analysis completed.\n")

