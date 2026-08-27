# ============================================================
# 05_tcell_adjusted_models.R
#
# Adjusted intratumoral T-cell-density analyses.
#
# Computes the statistics underlying:
#   Fig. 2c
#     adjusted geometric-mean ratios for CD4/CD8/FOXP3
#   Supplementary Table S6B
#     primary HC3 models + 2,000-resample bootstrap CIs
#   Supplementary Table S6C
#     additional adjustment for continuous HbA1c
#   Supplementary Table S6D
#     additional adjustment for BMI, SMI, serum albumin and ECOG-PS
#   Supplementary Table S5C1
#     descriptive and adjusted four-diabetes-phenotype analyses
#   Supplementary Table S5B
#     TAD-by-stage interaction for CD4/CD8/FOXP3
#
# This script outputs model results; final figure assembly is performed separately.
# Unadjusted descriptive/IHC comparisons performed in SPSS are
# are not reproduced by this R script.
#
# Primary model:
#   log(marker density) ~
#     TAD + age + pathologic stage + log2(CA19-9)
#     + continuous NLR + neoadjuvant chemotherapy
#
# HC3 heteroscedasticity-robust covariance is used.
# Robust Wald z inference reproduces the manuscript analysis.
# BH correction is applied across the three T-cell markers:
# CD4, CD8, FOXP3.
# Supporting four-category analyses are not included in this
# Benjamini-Hochberg correction.
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
  "BMI", "SMI", "albumin", "ECOG_PS",
  "CD4", "CD8", "FOXP3"
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

    BMI = as.numeric(BMI),
    SMI = as.numeric(SMI),
    albumin = as.numeric(albumin),
    ECOG_PS = as.numeric(ECOG_PS),

    CD4 = as.numeric(CD4),
    CD8 = as.numeric(CD8),
    FOXP3 = as.numeric(FOXP3)
  )

primary_model_cols <- c(
  "TAD", "age", "stage", "CA199_raw",
  "NLR_raw", "neoadj", "CD4", "CD8", "FOXP3"
)

body_status_cols <- c(
  "BMI", "SMI", "albumin", "ECOG_PS"
)

stopifnot(nrow(dat) == 162)
stopifnot(all(complete.cases(dat[primary_model_cols])))
stopifnot(all(complete.cases(dat[body_status_cols])))
stopifnot(sum(dat$TAD == 0) == 98)
stopifnot(sum(dat$TAD == 1) == 64)
stopifnot(all(dat$CA199_raw > 0, na.rm = TRUE))
stopifnot(all(dat$CD4 > 0, na.rm = TRUE))
stopifnot(all(dat$CD8 > 0, na.rm = TRUE))
stopifnot(all(dat$FOXP3 > 0, na.rm = TRUE))

dat$log2_CA199 <- log2(dat$CA199_raw)
dat$ECOG_hi <- as.integer(dat$ECOG_PS >= 1)

# Validate ECOG-PS distribution used in the analysis.
stopifnot(
  sum(dat$TAD == 0 & dat$ECOG_hi == 0) == 82,
  sum(dat$TAD == 0 & dat$ECOG_hi == 1) == 16,
  sum(dat$TAD == 1 & dat$ECOG_hi == 0) == 53,
  sum(dat$TAD == 1 & dat$ECOG_hi == 1) == 11
)

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

fit_marker_hc3 <- function(
    marker,
    add_hba1c = FALSE,
    add_body_status = FALSE,
    interaction = FALSE
) {

  if (add_hba1c && add_body_status) {
    stop("HbA1c and body-status sensitivity flags should be run separately.")
  }

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

  if (add_body_status) {
    rhs <- c(rhs, "BMI", "SMI", "albumin", "ECOG_hi")
  }

  f <- as.formula(
    paste0("log(", marker, ") ~ ", paste(rhs, collapse = " + "))
  )

  fit <- lm(f, data = dat)
  tab <- hc3_wald(fit)

  list(fit = fit, table = tab)
}

extract_tad <- function(
    marker,
    add_hba1c = FALSE,
    add_body_status = FALSE
) {
  z <- fit_marker_hc3(
    marker,
    add_hba1c = add_hba1c,
    add_body_status = add_body_status
  )
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
  file.path(out_dir, "Supplementary_Table_S6B_primary_HC3_and_bootstrap.csv"),
  row.names = FALSE
)

cat("\nPrimary adjusted immune models:\n")
print(primary)

# Reference estimates reported in the manuscript:
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
  file.path(out_dir, "Supplementary_Table_S6C_HbA1c_adjusted.csv"),
  row.names = FALSE
)

cat("\nHbA1c-adjusted immune models:\n")
print(hba1c)

# Reference estimates reported in the manuscript:
# CD4   0.73 (0.52-1.02), P=0.062, q=0.094
# CD8   0.69 (0.51-0.93), P=0.014, q=0.042
# FOXP3 0.78 (0.55-1.12), P=0.177, q=0.177


# ------------------------------------------------------------
# Additional adjustment for BMI, SMI, serum albumin and ECOG-PS
# Supplementary Table S6D
# ------------------------------------------------------------

body_status <- bind_rows(
  lapply(
    primary_markers,
    extract_tad,
    add_body_status = TRUE
  )
) %>%
  mutate(q_value = p.adjust(P_value, method = "BH"))

# Validate the body-status sensitivity estimates reported in Supplementary Table S6D.
reference_body_status_ratio <- c(
  CD4 = 0.76640799,
  CD8 = 0.77958052,
  FOXP3 = 0.80411894
)

observed_body_status_ratio <- setNames(
  body_status$Ratio,
  body_status$Marker
)

stopifnot(
  all(
    abs(
      observed_body_status_ratio[names(reference_body_status_ratio)] -
        reference_body_status_ratio
    ) < 1e-6
  )
)

write.csv(
  body_status,
  file.path(
    out_dir,
    "Supplementary_Table_S6D_BMI_SMI_albumin_ECOG_adjusted.csv"
  ),
  row.names = FALSE
)

cat("\nBMI/SMI/albumin/ECOG-PS-adjusted immune models:\n")
print(body_status)

# Reference estimates reported in the manuscript:
# CD4   0.77 (0.60-0.98), P=0.035, q=0.068
# CD8   0.78 (0.61-1.00), P=0.045, q=0.068
# FOXP3 0.80 (0.63-1.04), P=0.090, q=0.090

# ------------------------------------------------------------
# Four diabetes phenotypes
# Supplementary Table S5C1
#
# B1: descriptive median/IQR + Kruskal-Wallis
# B2: adjusted four-category HC3 models
#
# No diabetes is the reference category.
#
# The adjusted four-category analyses are supporting analyses
# and are not included in the BH correction applied to the
# main marker-specific TAD/non-TAD analyses.
# ------------------------------------------------------------


# ============================================================
# B1. Descriptive T-cell densities
# ============================================================

fourcat_descriptive <- bind_rows(
  lapply(primary_markers, function(marker) {

    kt <- kruskal.test(dat[[marker]] ~ dat$dm_group)

    med <- tapply(
      dat[[marker]],
      dat$dm_group,
      median,
      na.rm = TRUE
    )

    q1 <- tapply(
      dat[[marker]],
      dat$dm_group,
      function(x) {
        quantile(
          x,
          0.25,
          na.rm = TRUE,
          names = FALSE
        )
      }
    )

    q3 <- tapply(
      dat[[marker]],
      dat$dm_group,
      function(x) {
        quantile(
          x,
          0.75,
          na.rm = TRUE,
          names = FALSE
        )
      }
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
  fourcat_descriptive,
  file.path(
    out_dir,
    "Supplementary_Table_S5C1_descriptive_four_diabetes_phenotypes.csv"
  ),
  row.names = FALSE
)

cat("\nFour-category descriptive analyses:\n")
print(fourcat_descriptive)

cat(
  paste0(
    "Reference values: ",
    "Kruskal-Wallis P: ",
    "CD4 ~0.144; ",
    "CD8 ~0.147; ",
    "FOXP3 ~0.102.\n"
  )
)


# ============================================================
# B2. Adjusted four-category analyses
#
# Same covariate adjustment as the main T-cell models:
#
# log(marker) ~
#   diabetes phenotype
#   + age
#   + pathologic stage
#   + log2(CA19-9)
#   + continuous NLR
#   + neoadjuvant chemotherapy
#
# HC3 heteroscedasticity-robust covariance.
# ============================================================

fit_fourcat_hc3 <- function(marker) {

  f <- as.formula(
    paste0(
      "log(", marker, ") ~ ",
      "dm_group + age + stage + log2_CA199 + NLR_raw + neoadj"
    )
  )

  fit <- lm(f, data = dat)

  V <- vcovHC(
    fit,
    type = "HC3"
  )

  b <- coef(fit)
  se <- sqrt(diag(V))

  phenotype_terms <- grep(
    "^dm_group",
    names(b),
    value = TRUE
  )

  if (length(phenotype_terms) != 3) {
    stop(
      "Expected three diabetes-phenotype coefficients for ",
      marker,
      "; found ",
      length(phenotype_terms)
    )
  }

  beta <- b[phenotype_terms]
  robust_se <- se[phenotype_terms]

  z <- beta / robust_se
  p <- 2 * pnorm(-abs(z))

  zcrit <- qnorm(0.975)

  contrast_results <- data.frame(
    Marker = marker,
    Contrast = sub(
      "^dm_group",
      "",
      phenotype_terms
    ),
    Adjusted_ratio = exp(beta),
    Lower95 = exp(
      beta - zcrit * robust_se
    ),
    Upper95 = exp(
      beta + zcrit * robust_se
    ),
    P_value = p,
    stringsAsFactors = FALSE
  )

  # HC3 robust omnibus Wald test:
  # H0: all three diabetes-phenotype coefficients = 0
  V_sub <- V[
    phenotype_terms,
    phenotype_terms,
    drop = FALSE
  ]

  W <- as.numeric(
    t(beta) %*%
      solve(V_sub) %*%
      beta
  )

  omnibus_p <- pchisq(
    W,
    df = length(beta),
    lower.tail = FALSE
  )

  list(
    fit = fit,
    covariance = V,
    contrasts = contrast_results,
    omnibus_P = omnibus_p
  )
}


fourcat_adjusted_list <- lapply(
  primary_markers,
  fit_fourcat_hc3
)

fourcat_adjusted <- bind_rows(
  lapply(
    fourcat_adjusted_list,
    function(x) x$contrasts
  )
)

fourcat_omnibus <- data.frame(
  Marker = primary_markers,
  Omnibus_P = sapply(
    fourcat_adjusted_list,
    function(x) x$omnibus_P
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Wide output matching Supplementary Table S5C2
# ------------------------------------------------------------

extract_fourcat_row <- function(marker) {

  z <- fourcat_adjusted %>%
    filter(Marker == marker)

  get_contrast <- function(name) {

    x <- z %>%
      filter(Contrast == name)

    if (nrow(x) != 1) {
      stop(
        "Could not identify contrast '",
        name,
        "' for ",
        marker
      )
    }

    x
  }

  stable <- get_contrast(
    "Long-standing stable"
  )

  new_onset <- get_contrast(
    "New-onset"
  )

  worsening <- get_contrast(
    "Worsening"
  )

  omnibus <- fourcat_omnibus %>%
    filter(Marker == marker)

  data.frame(
    Marker = marker,

    Stable_ratio = stable$Adjusted_ratio,
    Stable_Lower95 = stable$Lower95,
    Stable_Upper95 = stable$Upper95,
    Stable_P = stable$P_value,

    New_onset_ratio = new_onset$Adjusted_ratio,
    New_onset_Lower95 = new_onset$Lower95,
    New_onset_Upper95 = new_onset$Upper95,
    New_onset_P = new_onset$P_value,

    Worsening_ratio = worsening$Adjusted_ratio,
    Worsening_Lower95 = worsening$Lower95,
    Worsening_Upper95 = worsening$Upper95,
    Worsening_P = worsening$P_value,

    Omnibus_P = omnibus$Omnibus_P,

    stringsAsFactors = FALSE
  )
}


fourcat_adjusted_wide <- bind_rows(
  lapply(
    primary_markers,
    extract_fourcat_row
  )
)


write.csv(
  fourcat_adjusted,
  file.path(
    out_dir,
    "Supplementary_Table_S5C2_adjusted_four_diabetes_phenotypes_long.csv"
  ),
  row.names = FALSE
)

write.csv(
  fourcat_adjusted_wide,
  file.path(
    out_dir,
    "Supplementary_Table_S5C2_adjusted_four_diabetes_phenotypes.csv"
  ),
  row.names = FALSE
)


cat("\nAdjusted four-category T-cell analyses:\n")
print(fourcat_adjusted_wide)

cat(
  paste0(
    "\nReference values:\n",
    "CD4: stable 1.08, new-onset 0.76, worsening 0.75; ",
    "omnibus P ~0.127\n",
    "CD8: stable 0.97, new-onset 0.76, worsening 0.82; ",
    "omnibus P ~0.197\n",
    "FOXP3: stable 0.86, new-onset 0.80, worsening 0.75; ",
    "omnibus P ~0.234\n"
  )
)

cat(
  paste0(
    "\nSupporting four-category analyses are not included ",
    "in the Benjamini-Hochberg correction applied to the ",
    "main marker-specific TAD/non-TAD analyses.\n"
  )
)

# ------------------------------------------------------------
# TAD-by-stage interaction in immune-density models:
# HC3 robust WALD test for the interaction coefficient
# Supplementary Table S5B
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
  file.path(out_dir, "Supplementary_Table_S5B_immune_stage_interactions.csv"),
  row.names = FALSE
)

cat("\nImmune TAD-by-stage interaction (HC3 Wald):\n")
print(interaction_results)
cat("Reference interaction P values: CD4 ~0.168; CD8 ~0.384; FOXP3 ~0.687.\n")

cat("\n============================================\n")
cat("T-CELL MODEL SUMMARY\n")
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

cat("\nBody-status sensitivity:\n")
print(body_status[, c(
  "Marker", "Ratio", "Lower95", "Upper95", "P_value", "q_value"
)])

cat("============================================\n")

sink(file.path(out_dir, "sessionInfo_tcell_models.txt"))
sessionInfo()
sink()