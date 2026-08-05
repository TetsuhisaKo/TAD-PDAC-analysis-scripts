# ============================================================
# 05_tcell_adjusted_models.R
#
# Adjusted intratumoral T-cell-density analyses.
#
# Computes the statistics underlying:
#   Fig. 2c
#     adjusted geometric-mean ratios for CD4/CD8/FOXP3
#   Supplementary Table S5A
#     primary HC3 models + 2,000-resample bootstrap CIs
#   Supplementary Table S5B
#     additional adjustment for continuous HbA1c
#   Supplementary Table S5C
#     adjusted CD45RO supporting outcome
#   Supplementary Table S4B
#     four-diabetes-phenotype Kruskal-Wallis comparisons
#   Supplementary Table S4C
#     TAD-by-stage interaction for CD4/CD8/FOXP3
#
# Plotting code is intentionally not included.
# Unadjusted descriptive/IHC comparisons performed in SPSS are
# intentionally not reproduced here.
#
# Primary model:
#   log(marker density) ~
#     TAD + age + pathologic stage + log2(CA19-9)
#     + continuous NLR + neoadjuvant chemotherapy
#
# HC3 heteroscedasticity-robust covariance is used.
# Robust Wald z inference reproduces the manuscript analysis.
# BH correction is applied across the three primary outcomes:
# CD4, CD8, FOXP3.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(sandwich)
})

in_file <- "data/clinical_data.xlsx"
out_dir <- "results/tcell_models"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(in_file))

raw <- read_excel(in_file, sheet = 1)

required <- c(
  "TAD", "dm_group", "age", "stage",
  "CA199_raw", "NLR_raw", "neoadj", "HbA1c",
  "CD4", "CD8", "FOXP3", "CD45RO"
)

missing <- setdiff(required, names(raw))
if (length(missing) > 0) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

dat <- raw %>%
  transmute(
    TAD = as.numeric(TAD),

    dm_group = factor(
      as.integer(dm_group),
      levels = c(0, 1, 2, 3),
      labels = c(
        "No diabetes",
        "Long-standing stable",
        "New-onset",
        "Worsening"
      )
    ),

    age = as.numeric(age),
    stage = as.numeric(stage),
    CA199_raw = as.numeric(CA199_raw),
    NLR_raw = as.numeric(NLR_raw),
    neoadj = as.numeric(neoadj),
    HbA1c = as.numeric(HbA1c),

    CD4 = as.numeric(CD4),
    CD8 = as.numeric(CD8),
    FOXP3 = as.numeric(FOXP3),
    CD45RO = as.numeric(CD45RO)
  )

primary_model_cols <- c(
  "TAD", "age", "stage", "CA199_raw",
  "NLR_raw", "neoadj", "CD4", "CD8", "FOXP3"
)

stopifnot(nrow(dat) == 162)
stopifnot(all(complete.cases(dat[primary_model_cols])))
stopifnot(sum(dat$TAD == 0) == 98)
stopifnot(sum(dat$TAD == 1) == 64)
stopifnot(all(dat$CA199_raw > 0, na.rm = TRUE))
stopifnot(all(dat$CD4 > 0, na.rm = TRUE))
stopifnot(all(dat$CD8 > 0, na.rm = TRUE))
stopifnot(all(dat$FOXP3 > 0, na.rm = TRUE))

dat$log2_CA199 <- log2(dat$CA199_raw)

# ------------------------------------------------------------
# HC3 robust Wald helper
# ------------------------------------------------------------

hc3_wald <- function(fit) {
  V <- vcovHC(fit, type = "HC3")
  b <- coef(fit)
  se <- sqrt(diag(V))
  z <- b / se
  p <- 2 * pnorm(-abs(z))

  data.frame(
    term = names(b),
    estimate = unname(b),
    robust_se = unname(se),
    z = unname(z),
    p = unname(p),
    stringsAsFactors = FALSE
  )
}

fit_marker_hc3 <- function(marker, add_hba1c = FALSE, interaction = FALSE) {

  if (interaction) {
    rhs <- c(
      "TAD * stage",
      "age",
      "log2_CA199",
      "NLR_raw",
      "neoadj"
    )
  } else {
    rhs <- c(
      "TAD",
      "age",
      "stage",
      "log2_CA199",
      "NLR_raw",
      "neoadj"
    )
  }

  if (add_hba1c) {
    rhs <- c(rhs, "HbA1c")
  }

  f <- as.formula(
    paste0("log(", marker, ") ~ ", paste(rhs, collapse = " + "))
  )

  fit <- lm(f, data = dat)
  tab <- hc3_wald(fit)

  list(fit = fit, table = tab)
}

extract_tad <- function(marker, add_hba1c = FALSE) {
  z <- fit_marker_hc3(marker, add_hba1c = add_hba1c)
  row <- z$table[z$table$term == "TAD", , drop = FALSE]

  if (nrow(row) != 1) {
    stop("Could not identify TAD coefficient for ", marker)
  }

  zcrit <- qnorm(0.975)

  data.frame(
    Marker = marker,
    Ratio = exp(row$estimate),
    Lower95 = exp(row$estimate - zcrit * row$robust_se),
    Upper95 = exp(row$estimate + zcrit * row$robust_se),
    P_value = row$p,
    stringsAsFactors = FALSE
  )
}

primary_markers <- c("CD4", "CD8", "FOXP3")

primary <- bind_rows(
  lapply(primary_markers, extract_tad, add_hba1c = FALSE)
) %>%
  mutate(q_value = p.adjust(P_value, method = "BH"))

# ------------------------------------------------------------
# 2,000 patient-level bootstrap percentile confidence intervals
# Fixed seed ensures reproducibility.
# ------------------------------------------------------------

BOOT_B <- 2000L
BOOT_SEED <- 20260803L
set.seed(BOOT_SEED)

bootstrap_marker <- function(marker) {
  n <- nrow(dat)
  vals <- rep(NA_real_, BOOT_B)

  f <- as.formula(
    paste0(
      "log(", marker, ") ~ ",
      "TAD + age + stage + log2_CA199 + NLR_raw + neoadj"
    )
  )

  for (b in seq_len(BOOT_B)) {
    idx <- sample.int(n, size = n, replace = TRUE)

    fit_b <- try(
      lm(f, data = dat[idx, , drop = FALSE]),
      silent = TRUE
    )

    if (!inherits(fit_b, "try-error")) {
      cf <- coef(fit_b)

      if ("TAD" %in% names(cf) && is.finite(cf["TAD"])) {
        vals[b] <- exp(cf["TAD"])
      }
    }
  }

  vals <- vals[is.finite(vals)]

  if (length(vals) < 0.95 * BOOT_B) {
    warning("More than 5% of bootstrap fits failed for ", marker)
  }

  ci <- quantile(
    vals,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )

  data.frame(
    Marker = marker,
    Bootstrap_B_requested = BOOT_B,
    Bootstrap_B_successful = length(vals),
    Bootstrap_seed = BOOT_SEED,
    Bootstrap_Lower95 = ci[1],
    Bootstrap_Upper95 = ci[2],
    stringsAsFactors = FALSE
  )
}

bootstrap_results <- bind_rows(
  lapply(primary_markers, bootstrap_marker)
)

primary <- primary %>%
  left_join(bootstrap_results, by = "Marker")

write.csv(
  primary,
  file.path(out_dir, "Supplementary_Table_S5A_primary_HC3_and_bootstrap.csv"),
  row.names = FALSE
)

cat("\nPrimary adjusted immune models:\n")
print(primary)

# Current manuscript checks:
# CD4   0.74 (0.58-0.95), P=0.016, q=0.049
#        bootstrap 0.59-0.93
# CD8   0.79 (0.63-0.99), P=0.037, q=0.055
#        bootstrap 0.63-0.99
# FOXP3 0.80 (0.63-1.03), P=0.079, q=0.079
#        bootstrap 0.63-1.02

# ------------------------------------------------------------
# HbA1c-adjusted sensitivity
# ------------------------------------------------------------

if (!all(complete.cases(dat[, c(
  "TAD", "age", "stage", "log2_CA199",
  "NLR_raw", "neoadj", "HbA1c",
  "CD4", "CD8", "FOXP3"
)]))) {
  warning("HbA1c sensitivity models use complete cases.")
}

hba1c <- bind_rows(
  lapply(primary_markers, extract_tad, add_hba1c = TRUE)
) %>%
  mutate(q_value = p.adjust(P_value, method = "BH"))

write.csv(
  hba1c,
  file.path(out_dir, "Supplementary_Table_S5B_HbA1c_adjusted.csv"),
  row.names = FALSE
)

cat("\nHbA1c-adjusted immune models:\n")
print(hba1c)

# Current manuscript checks:
# CD4   0.73 (0.52-1.02), P=0.062, q=0.094
# CD8   0.69 (0.51-0.93), P=0.014, q=0.042
# FOXP3 0.78 (0.55-1.12), P=0.177, q=0.177

# ------------------------------------------------------------
# Additional CD45RO outcome: nominal P only
# ------------------------------------------------------------

if (all(complete.cases(dat[, c(
  "TAD", "age", "stage", "log2_CA199",
  "NLR_raw", "neoadj", "CD45RO"
)])) && all(dat$CD45RO > 0)) {

  cd45ro <- extract_tad("CD45RO", add_hba1c = FALSE)

  write.csv(
    cd45ro,
    file.path(out_dir, "Supplementary_Table_S5C_CD45RO_adjusted.csv"),
    row.names = FALSE
  )

  cat("\nCD45RO supporting outcome:\n")
  print(cd45ro)
  cat("Expected: ratio ~0.92 (0.73-1.16), P~0.494.\n")

} else {
  warning("CD45RO supporting outcome skipped because of missing/non-positive values.")
}

# ------------------------------------------------------------
# Four diabetes phenotypes: Kruskal-Wallis
# Supplementary Table S4B
# ------------------------------------------------------------

fourcat_results <- bind_rows(
  lapply(primary_markers, function(marker) {

    kt <- kruskal.test(dat[[marker]] ~ dat$dm_group)

    med <- tapply(dat[[marker]], dat$dm_group, median, na.rm = TRUE)

    q1 <- tapply(
      dat[[marker]], dat$dm_group,
      function(x) quantile(x, 0.25, na.rm = TRUE, names = FALSE)
    )

    q3 <- tapply(
      dat[[marker]], dat$dm_group,
      function(x) quantile(x, 0.75, na.rm = TRUE, names = FALSE)
    )

    data.frame(
      Marker = marker,

      No_DM_median = med["No diabetes"],
      No_DM_Q1 = q1["No diabetes"],
      No_DM_Q3 = q3["No diabetes"],

      Stable_DM_median = med["Long-standing stable"],
      Stable_DM_Q1 = q1["Long-standing stable"],
      Stable_DM_Q3 = q3["Long-standing stable"],

      New_onset_median = med["New-onset"],
      New_onset_Q1 = q1["New-onset"],
      New_onset_Q3 = q3["New-onset"],

      Worsening_median = med["Worsening"],
      Worsening_Q1 = q1["Worsening"],
      Worsening_Q3 = q3["Worsening"],

      Kruskal_Wallis_P = kt$p.value,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  fourcat_results,
  file.path(out_dir, "Supplementary_Table_S4B_four_diabetes_phenotypes.csv"),
  row.names = FALSE
)

cat("\nFour-category Kruskal-Wallis results:\n")
print(fourcat_results)
cat("Expected P: CD4 ~0.144; CD8 ~0.147; FOXP3 ~0.102.\n")

# ------------------------------------------------------------
# TAD-by-stage interaction in immune-density models:
# HC3 robust WALD test for the interaction coefficient
# Supplementary Table S4C
# ------------------------------------------------------------

interaction_results <- bind_rows(
  lapply(primary_markers, function(marker) {

    z <- fit_marker_hc3(marker, interaction = TRUE)

    int_rows <- z$table %>%
      filter(grepl("TAD:stage|stage:TAD", term))

    if (nrow(int_rows) != 1) {
      stop("Could not identify TAD-by-stage term for ", marker)
    }

    data.frame(
      Marker = marker,
      Test = "HC3 robust Wald test for TAD-by-stage interaction term",
      Interaction_term = int_rows$term,
      Wald_z = int_rows$z,
      P_value = int_rows$p,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  interaction_results,
  file.path(out_dir, "Supplementary_Table_S4C_immune_stage_interactions.csv"),
  row.names = FALSE
)

cat("\nImmune TAD-by-stage interaction (HC3 Wald):\n")
print(interaction_results)
cat("Expected P: CD4 ~0.168; CD8 ~0.384; FOXP3 ~0.687.\n")

cat("\n============================================\n")
cat("FINAL T-CELL MODEL CHECK\n")
cat("============================================\n")
cat(
  "n =", nrow(dat),
  "| non-TAD =", sum(dat$TAD == 0),
  "| TAD =", sum(dat$TAD == 1), "\n"
)
print(primary[, c(
  "Marker", "Ratio", "Lower95", "Upper95", "P_value", "q_value",
  "Bootstrap_Lower95", "Bootstrap_Upper95"
)])
cat("============================================\n")

sink(file.path(out_dir, "sessionInfo_tcell_models.txt"))
sessionInfo()
sink()
