# ============================================================
# 01_survival_analysis.R
#
# Survival analyses in the surgical cohort (n = 162)
#
# Outcomes:
#   Overall survival (OS)
#   Cancer-specific survival (CSS)
#
# Main / supplementary results:
#   Fig. 1b: TAD vs non-TAD OS
#   Fig. 1c: OS according to four diabetes phenotypes
#   Table 2: primary multivariable OS model
#   Supplementary Fig. S1: TAD vs non-TAD CSS
#   Supplementary Fig. S2: stage-stratified survival analyses
#   Supplementary Table S3: CSS sensitivity analyses
#   Supplementary Table S4: OS sensitivity analyses
#   Supplementary Table S5A: four diabetes phenotypes and
#                            diabetes-only survival analyses
#   Supplementary Table S5B: TAD-by-stage interaction tests
#
# Primary multivariable model:
#   TAD + pathologic stage + resection margin + NLR
#       + CA19-9 >37 U/mL
#
# Sensitivity analyses:
#   + continuous HbA1c
#   TAD replaced by continuous HbA1c
#   + adjuvant chemotherapy
#   + continuous age
#   + continuous surgery year
#   continuous log2 CA19-9 instead of dichotomised CA19-9
#   continuous NLR + continuous log2 CA19-9 (OS only)
#   continuous log2 CA19-9 + adjuvant chemotherapy
#   exclusion of neoadjuvant-chemotherapy patients
#   exclusion of treatment-only worsening diabetes
#   diabetes-only exclusion of treatment-only worsening diabetes
#   diabetes-only restriction of stable diabetes to patients with
#       available longitudinal glycaemic data
#   + BMI, SMI, serum albumin and ECOG-PS
#
# Additional analyses:
#   four diabetes phenotypes
#   TAD vs stable diabetes among patients with diabetes
#   same diabetes-only model + HbA1c
#   reduced diabetes-only model: stage + margin
#   stage-specific Cox models
#   TAD-by-stage interaction
#
# Cohort-count validation:
#   Total n = 162
#   non-TAD = 98; TAD = 64
#   no diabetes = 81
#   stable diabetes = 17
#   new-onset diabetes = 40
#   worsening diabetes = 24
#   worsening basis: HbA1c only / treatment only / both = 9 / 11 / 4
#   stable diabetes with / without longitudinal glycaemic data = 8 / 9
#   ECOG-PS 0 / >=1: non-TAD = 82 / 16; TAD = 53 / 11
#
# OS event-count validation:
#   total = 120
#   TAD = 52
#   stable diabetes = 12
#
# Stage-specific count validation:
#   <=IIa: n = 43; OS events = 24; CSS events = 20
#   >=IIb: n = 119; OS events = 96; CSS events = 90
# ============================================================


suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readxl)
  library(broom)
})


# ============================================================
# Helper for input-column names
# ============================================================

get_required_col <- function(data, candidates, label) {

  hit <- candidates[candidates %in% names(data)]

  if (length(hit) == 0) {
    stop(
      paste0(
        "Required column for ", label, " not found. Tried: ",
        paste(candidates, collapse = ", ")
      )
    )
  }

  data[[hit[1]]]
}


# ============================================================
# Data import
# ============================================================

raw <- read_excel("data/clinical_data.xlsx", sheet = 1)


# Existing repository variable names are listed first.
# Alternative capitalisations are allowed only to make the
# script robust to spreadsheet-header formatting.

df <- data.frame(

  os_time =
    as.numeric(
      get_required_col(
        raw,
        c("os_time"),
        "OS time"
      )
    ),

  os_event =
    as.integer(
      get_required_col(
        raw,
        c("os_event"),
        "OS event"
      )
    ),

  css_event =
    as.integer(
      get_required_col(
        raw,
        c("css_event"),
        "CSS event"
      )
    ),

  TAD_num =
    as.integer(
      get_required_col(
        raw,
        c("TAD"),
        "TAD"
      )
    ),

  stage_num =
    as.integer(
      get_required_col(
        raw,
        c("stage"),
        "pathologic stage"
      )
    ),

  margin_num =
    as.integer(
      get_required_col(
        raw,
        c("margin"),
        "resection margin"
      )
    ),

  NLR_num =
    as.integer(
      get_required_col(
        raw,
        c("NLR"),
        "NLR category"
      )
    ),

  NLR_raw =
    as.numeric(
      get_required_col(
        raw,
        c("NLR_raw"),
        "continuous NLR"
      )
    ),

  CA199_num =
    as.integer(
      get_required_col(
        raw,
        c("CA199"),
        "CA19-9 category"
      )
    ),

  CA199_raw =
    as.numeric(
      get_required_col(
        raw,
        c("CA199_raw"),
        "continuous CA19-9"
      )
    ),

  adj_num =
    as.integer(
      get_required_col(
        raw,
        c("adj_chemo"),
        "adjuvant chemotherapy"
      )
    ),

  neoadj_num =
    as.integer(
      get_required_col(
        raw,
        c("neoadj"),
        "neoadjuvant chemotherapy"
      )
    ),

  HbA1c =
    as.numeric(
      get_required_col(
        raw,
        c("HbA1c", "hba1c"),
        "HbA1c"
      )
    ),

  age =
    as.numeric(
      get_required_col(
        raw,
        c("age", "Age"),
        "age"
      )
    ),

  surgery_year =
    as.numeric(
      get_required_col(
        raw,
        c("surgery_year", "Surgery_year"),
        "surgery year"
      )
    ),

  dm_group_num =
    as.integer(
      get_required_col(
        raw,
        c("dm_group"),
        "four-category diabetes phenotype"
      )
    ),

  worsening_basis =
    as.integer(
      get_required_col(
        raw,
        c("worsening_basis"),
        "worsening-diabetes classification basis"
      )
    ),

  stable_longitudinal_data =
    as.integer(
      get_required_col(
        raw,
        c("stable_longitudinal_data"),
        "stable-diabetes longitudinal glycaemic data availability"
      )
    ),

  BMI =
    as.numeric(
      get_required_col(
        raw,
        c("BMI", "bmi"),
        "body mass index"
      )
    ),

  SMI =
    as.numeric(
      get_required_col(
        raw,
        c("SMI", "smi"),
        "skeletal muscle index"
      )
    ),

  albumin =
    as.numeric(
      get_required_col(
        raw,
        c("albumin", "Albumin"),
        "serum albumin"
      )
    ),

  ECOG_PS =
    as.integer(
      get_required_col(
        raw,
        c("ECOG_PS", "ECOG-PS", "ECOGPS"),
        "ECOG performance status"
      )
    )
)


# ============================================================
# Variable coding
# ============================================================

df <- df %>%
  mutate(

    TAD = factor(
      TAD_num,
      levels = c(0, 1),
      labels = c("non-TAD", "TAD")
    ),

    stage = factor(
      stage_num,
      levels = c(0, 1),
      labels = c("<=IIa", ">=IIb")
    ),

    margin = factor(
      margin_num,
      levels = c(0, 1),
      labels = c("R0", "R1")
    ),

    NLR = factor(
      NLR_num,
      levels = c(0, 1),
      labels = c("<2.1", ">=2.1")
    ),

    CA199 = factor(
      CA199_num,
      levels = c(0, 1),
      labels = c("<=37", ">37")
    ),

    adj_chemo = factor(
      adj_num,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),

    neoadj = factor(
      neoadj_num,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),

    dm4 = factor(
      dm_group_num,
      levels = c(0, 1, 2, 3),
      labels = c(
        "No diabetes",
        "Stable diabetes",
        "New-onset diabetes",
        "Worsening diabetes"
      )
    ),

    dm_compare = factor(
      case_when(
        dm_group_num == 1 ~ "Stable diabetes",
        dm_group_num %in% c(2, 3) ~ "TAD",
        TRUE ~ NA_character_
      ),
      levels = c("Stable diabetes", "TAD")
    ),

    ECOG_hi = as.integer(ECOG_PS >= 1),

    CA199_log2 = log2(CA199_raw)
  )


# ============================================================
# Data checks
# ============================================================

stopifnot(nrow(df) == 162)

stopifnot(sum(df$TAD == "TAD") == 64)
stopifnot(sum(df$TAD == "non-TAD") == 98)

stopifnot(sum(df$dm4 == "No diabetes") == 81)
stopifnot(sum(df$dm4 == "Stable diabetes") == 17)
stopifnot(sum(df$dm4 == "New-onset diabetes") == 40)
stopifnot(sum(df$dm4 == "Worsening diabetes") == 24)

# TAD must correspond exactly to new-onset + worsening diabetes.
stopifnot(
  all(
    df$TAD_num ==
      as.integer(df$dm_group_num %in% c(2, 3))
  )
)

# Worsening-diabetes classification basis:
#   1 = HbA1c increase only
#   2 = treatment intensification only
#   3 = both
stopifnot(
  all(is.na(df$worsening_basis[df$dm_group_num != 3])),
  all(df$worsening_basis[df$dm_group_num == 3] %in% c(1, 2, 3)),
  sum(df$worsening_basis == 1, na.rm = TRUE) == 9,
  sum(df$worsening_basis == 2, na.rm = TRUE) == 11,
  sum(df$worsening_basis == 3, na.rm = TRUE) == 4
)

# Long-standing stable diabetes:
#   1 = longitudinal glycaemic data available
#   0 = unavailable
stopifnot(
  all(is.na(df$stable_longitudinal_data[df$dm_group_num != 1])),
  all(df$stable_longitudinal_data[df$dm_group_num == 1] %in% c(0, 1)),
  sum(df$stable_longitudinal_data == 1, na.rm = TRUE) == 8,
  sum(df$stable_longitudinal_data == 0, na.rm = TRUE) == 9
)

# Body-status covariates are complete in the 162-patient cohort.
stopifnot(
  !anyNA(df$BMI),
  !anyNA(df$SMI),
  !anyNA(df$albumin),
  !anyNA(df$ECOG_PS),
  all(df$ECOG_PS >= 0)
)

# Correct ECOG-PS distribution.
stopifnot(
  sum(df$TAD_num == 0 & df$ECOG_hi == 0) == 82,
  sum(df$TAD_num == 0 & df$ECOG_hi == 1) == 16,
  sum(df$TAD_num == 1 & df$ECOG_hi == 0) == 53,
  sum(df$TAD_num == 1 & df$ECOG_hi == 1) == 11
)

stopifnot(sum(df$os_event) == 120)
stopifnot(sum(df$css_event) == 110)

stopifnot(all(df$CA199_raw > 0))
stopifnot(!anyNA(df$NLR_raw))
stopifnot(all(df$NLR_raw >= 0))

cat("\nCohort checks passed.\n")


# ============================================================
# Median follow-up: reverse Kaplan-Meier
# ============================================================

reverse_km <- survfit(
  Surv(os_time, 1 - os_event) ~ 1,
  data = df
)

cat("\n================ MEDIAN FOLLOW-UP ================\n")
print(reverse_km$table["median"])


# ============================================================
# Kaplan-Meier / log-rank
# ============================================================

cat("\n================ KM / LOG-RANK ================\n")

cat("\n-- OS: TAD vs non-TAD (Fig. 1b) --\n")
print(
  survdiff(
    Surv(os_time, os_event) ~ TAD,
    data = df
  )
)

cat("\n-- CSS: TAD vs non-TAD (Supplementary Fig. S1) --\n")
print(
  survdiff(
    Surv(os_time, css_event) ~ TAD,
    data = df
  )
)

cat("\n-- OS: four diabetes phenotypes (Fig. 1c) --\n")
print(
  survdiff(
    Surv(os_time, os_event) ~ dm4,
    data = df
  )
)


# ============================================================
# Helper functions
# ============================================================

tidy_cox <- function(fit) {

  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    select(
      term,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
}


extract_term <- function(fit, pattern) {

  tidy_cox(fit) %>%
    filter(grepl(pattern, term))
}


key_result <- function(
    fit,
    model,
    endpoint,
    pattern = "^TADTAD$"
) {

  x <- extract_term(fit, pattern)

  if (nrow(x) != 1) {
    stop(
      paste(
        "Expected one matching coefficient for",
        model,
        "but found",
        nrow(x)
      )
    )
  }

  data.frame(
    Endpoint = endpoint,
    Model = model,
    HR = x$estimate,
    CI_low = x$conf.low,
    CI_high = x$conf.high,
    P_value = x$p.value
  )
}


# ============================================================
# PRIMARY MULTIVARIABLE MODELS
#
# Table 2 / Supplementary Tables S2-S3
# ============================================================

cox_os_primary <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199,
  data = df
)

cox_css_primary <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR + CA199,
  data = df
)


cat("\n================ PRIMARY OS ================\n")
print(tidy_cox(cox_os_primary))

cat("\n================ PRIMARY CSS ================\n")
print(tidy_cox(cox_css_primary))


# Reference TAD estimates reported in the manuscript:
# OS:  HR 1.84 (1.25-2.72), P = 0.002
# CSS: HR 1.87 (1.25-2.81), P = 0.002


# ============================================================
# Proportional-hazards assumption
# ============================================================

cat("\n================ PH ASSUMPTION: OS ================\n")
print(cox.zph(cox_os_primary))

cat("\n================ PH ASSUMPTION: CSS ================\n")
print(cox.zph(cox_css_primary))


# ============================================================
# SENSITIVITY: + continuous HbA1c
# ============================================================

cox_os_hba1c <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199 + HbA1c,
  data = df
)

cox_css_hba1c <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR + CA199 + HbA1c,
  data = df
)

# Reference TAD estimates reported in the manuscript:
# OS:  HR 2.18 (1.31-3.62), P = 0.003
# CSS: HR 2.40 (1.41-4.07), P = 0.001


# ============================================================
# SENSITIVITY: HbA1c replacing TAD
#
# Same clinicopathological covariates as the primary model,
# with continuous HbA1c entered instead of TAD.
# ============================================================

cox_os_hba1c_replace <- coxph(
  Surv(os_time, os_event) ~
    HbA1c + stage + margin + NLR + CA199,
  data = df
)

cox_css_hba1c_replace <- coxph(
  Surv(os_time, css_event) ~
    HbA1c + stage + margin + NLR + CA199,
  data = df
)

# Reference HbA1c estimates reported in the manuscript (per 1% increase):
# OS:  HR 1.08 (0.96-1.21), P = 0.220
# CSS: HR 1.06 (0.94-1.20), P = 0.359


# ============================================================
# SENSITIVITY: + adjuvant chemotherapy
# ============================================================

cox_os_adj <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199 + adj_chemo,
  data = df
)

cox_css_adj <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR + CA199 + adj_chemo,
  data = df
)

# Reference TAD estimates reported in the manuscript:
# OS:  HR 1.64 (1.09-2.46), P = 0.018
# CSS: HR 1.68 (1.10-2.56), P = 0.016


# ============================================================
# SENSITIVITY: + continuous age
# ============================================================

cox_os_age <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199 + age,
  data = df
)

cox_css_age <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR + CA199 + age,
  data = df
)

# Reference TAD estimates reported in the manuscript:
# OS:  HR 1.81 (1.22-2.67), P = 0.003
# CSS: HR 1.85 (1.23-2.77), P = 0.003


# ============================================================
# SENSITIVITY: + continuous surgery year
# ============================================================

cox_os_year <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199 + surgery_year,
  data = df
)

cox_css_year <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR + CA199 + surgery_year,
  data = df
)

# Reference TAD estimates reported in the manuscript:
# OS:  HR 1.85 (1.26-2.74), P = 0.002
# CSS: HR 1.89 (1.26-2.83), P = 0.002
#
# Surgery-year coefficient itself:
# OS:  HR 0.97/year (0.91-1.05), P = 0.497
# CSS: HR 0.98/year (0.91-1.06), P = 0.624


# ============================================================
# SENSITIVITY: log2-transformed continuous CA19-9
# ============================================================

cox_os_logca <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199_log2,
  data = df
)

cox_css_logca <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR + CA199_log2,
  data = df
)

# Reference TAD estimates reported in the manuscript:
# OS:  HR 1.86 (1.26-2.74), P = 0.002
# CSS: HR 1.93 (1.29-2.89), P = 0.001


# ============================================================
# SENSITIVITY (OS ONLY):
# continuous NLR + log2-transformed continuous CA19-9
#
# This model corresponds to Supplementary Table S4 and is not
# applied to CSS / Supplementary Table S3.
# ============================================================

cox_os_continuous_nlr_logca <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR_raw + CA199_log2,
  data = df
)

stopifnot(nobs(cox_os_continuous_nlr_logca) == 162)

cat("\n================ OS: CONTINUOUS NLR + LOG2 CA19-9 ================\n")
print(tidy_cox(cox_os_continuous_nlr_logca))

# Reference TAD estimate reported in the manuscript:
# HR 1.93 (1.31-2.85), P = 0.001


# ============================================================
# SENSITIVITY:
# log2-transformed continuous CA19-9 + adjuvant chemotherapy
# ============================================================

cox_os_logca_adj <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR +
    CA199_log2 + adj_chemo,
  data = df
)

cox_css_logca_adj <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR +
    CA199_log2 + adj_chemo,
  data = df
)

# Reference TAD estimates reported in the manuscript:
# OS:  HR 1.63 (1.09-2.45), P = 0.018
# CSS: HR 1.71 (1.12-2.61), P = 0.012


# ============================================================
# SENSITIVITY: excluding neoadjuvant chemotherapy
# ============================================================

df_noNAC <- df %>%
  filter(neoadj == "No")

stopifnot(nrow(df_noNAC) == 144)

cox_os_noNAC <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199,
  data = df_noNAC
)

cox_css_noNAC <- coxph(
  Surv(os_time, css_event) ~
    TAD + stage + margin + NLR + CA199,
  data = df_noNAC
)

# Reference TAD estimates reported in the manuscript:
# OS:  HR 1.73 (1.13-2.62), P = 0.011
# CSS: HR 1.78 (1.16-2.73), P = 0.009


# ============================================================
# CLASSIFICATION SENSITIVITY:
# excluding treatment-only worsening diabetes
#
# worsening_basis:
#   1 = HbA1c increase only
#   2 = treatment intensification only
#   3 = both
# ============================================================

df_no_treatment_only <- df %>%
  filter(
    !(
      dm_group_num == 3 &
        coalesce(worsening_basis == 2L, FALSE)
    )
  )

stopifnot(nrow(df_no_treatment_only) == 151)
stopifnot(sum(df_no_treatment_only$os_event) == 110)

cox_os_no_treatment_only <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199,
  data = df_no_treatment_only
)

stopifnot(nobs(cox_os_no_treatment_only) == 151)

# Reference result reported in the manuscript:
# HR 1.70 (1.12-2.56), P = 0.012


# ============================================================
# CLASSIFICATION SENSITIVITY:
# diabetes only, excluding treatment-only worsening diabetes
# ============================================================

df_dm_no_treatment_only <- df %>%
  filter(
    !is.na(dm_compare),
    !(
      dm_group_num == 3 &
        coalesce(worsening_basis == 2L, FALSE)
    )
  )

stopifnot(nrow(df_dm_no_treatment_only) == 70)
stopifnot(sum(df_dm_no_treatment_only$os_event) == 54)
stopifnot(
  sum(df_dm_no_treatment_only$dm_compare == "Stable diabetes") == 17,
  sum(df_dm_no_treatment_only$dm_compare == "TAD") == 53
)

cox_dm_no_treatment_only <- coxph(
  Surv(os_time, os_event) ~
    dm_compare + stage + margin + NLR + CA199,
  data = df_dm_no_treatment_only
)

stopifnot(nobs(cox_dm_no_treatment_only) == 70)

# Reference result reported in the manuscript:
# HR 2.21 (1.12-4.36), P = 0.022


# ============================================================
# CLASSIFICATION SENSITIVITY:
# diabetes only, stable diabetes restricted to patients
# with available longitudinal glycaemic data
# ============================================================

df_dm_stable_longitudinal <- df %>%
  filter(
    dm_group_num %in% c(2, 3) |
      (
        dm_group_num == 1 &
          stable_longitudinal_data == 1L
      )
  )

stopifnot(nrow(df_dm_stable_longitudinal) == 72)
stopifnot(sum(df_dm_stable_longitudinal$os_event) == 59)
stopifnot(
  sum(df_dm_stable_longitudinal$dm_compare == "Stable diabetes") == 8,
  sum(df_dm_stable_longitudinal$dm_compare == "TAD") == 64
)

cox_dm_stable_longitudinal <- coxph(
  Surv(os_time, os_event) ~
    dm_compare + stage + margin + NLR + CA199,
  data = df_dm_stable_longitudinal
)

stopifnot(nobs(cox_dm_stable_longitudinal) == 72)

# Reference result reported in the manuscript:
# HR 2.13 (0.92-4.93), P = 0.076


# ============================================================
# SENSITIVITY:
# + BMI, SMI, serum albumin and ECOG-PS
# ============================================================

cox_os_body_status <- coxph(
  Surv(os_time, os_event) ~
    TAD + stage + margin + NLR + CA199 +
    BMI + SMI + albumin + ECOG_hi,
  data = df
)

stopifnot(nobs(cox_os_body_status) == 162)

# Reference result reported in the manuscript:
# HR 1.70 (1.14-2.51), P = 0.009


cat("\n================ CLASSIFICATION / BODY-STATUS SENSITIVITIES ================\n")

cat("\n-- Excluding treatment-only worsening --\n")
print(extract_term(cox_os_no_treatment_only, "^TADTAD$"))

cat("\n-- Diabetes only, excluding treatment-only worsening --\n")
print(extract_term(cox_dm_no_treatment_only, "^dm_compareTAD$"))

cat("\n-- Diabetes only, stable DM restricted to longitudinal data --\n")
print(extract_term(cox_dm_stable_longitudinal, "^dm_compareTAD$"))

cat("\n-- + BMI, SMI, serum albumin and ECOG-PS --\n")
print(extract_term(cox_os_body_status, "^TADTAD$"))


# ============================================================
# COLLINEARITY:
# primary model + HbA1c
# ============================================================

if (requireNamespace("car", quietly = TRUE)) {

  cat("\n================ VIF: OS + HbA1c ================\n")

  vif_values <- car::vif(cox_os_hba1c)
  print(vif_values)

  if (is.matrix(vif_values)) {
    # Generalised VIF output, if generated.
    adjusted_vif <- vif_values[, "GVIF"]^(
      1 / (2 * vif_values[, "Df"])
    )
    cat(
      "Maximum adjusted GVIF:",
      max(adjusted_vif),
      "\n"
    )
  } else {
    cat(
      "Maximum VIF:",
      max(vif_values),
      "\n"
    )
  }

} else {

  cat(
    "\nPackage 'car' is not installed; VIF calculation skipped.\n"
  )
}


# ============================================================
# SUPPLEMENTARY TABLE S3 / S4:
# collect sensitivity results
#
# S3 = cancer-specific survival (CSS)
# S4 = overall survival (OS)
#
# The continuous-NLR + log2-CA19-9 model is OS only and is
# therefore included in S4 but not S3.
# ============================================================

results_os <- bind_rows(

  key_result(
    cox_os_primary,
    "Primary multivariable model",
    "OS"
  ),

  key_result(
    cox_os_hba1c,
    "+ HbA1c",
    "OS"
  ),

  key_result(
    cox_os_adj,
    "+ Adjuvant chemotherapy",
    "OS"
  ),

  key_result(
    cox_os_age,
    "+ continuous age",
    "OS"
  ),

  key_result(
    cox_os_year,
    "+ continuous surgery year",
    "OS"
  ),

  key_result(
    cox_os_logca,
    "Using log2-transformed CA19-9",
    "OS"
  ),

  key_result(
    cox_os_continuous_nlr_logca,
    "Using continuous NLR and log2-transformed CA19-9",
    "OS"
  ),

  key_result(
    cox_os_logca_adj,
    "Using log2-transformed CA19-9 + Adjuvant chemotherapy",
    "OS"
  ),

  key_result(
    cox_os_noNAC,
    "Excluding neoadjuvant chemotherapy",
    "OS"
  ),

  key_result(
    cox_os_hba1c_replace,
    "Replacing TAD with HbA1c",
    "OS",
    "^HbA1c$"
  ),

  key_result(
    cox_os_no_treatment_only,
    "Excluding treatment-only worsening",
    "OS"
  ),

  key_result(
    cox_dm_no_treatment_only,
    "Diabetes only, excluding treatment-only worsening",
    "OS",
    "^dm_compareTAD$"
  ),

  key_result(
    cox_dm_stable_longitudinal,
    "Diabetes only, stable DM with longitudinal data",
    "OS",
    "^dm_compareTAD$"
  ),

  key_result(
    cox_os_body_status,
    "+ BMI, SMI, serum albumin and ECOG-PS",
    "OS"
  )
)


results_css <- bind_rows(

  key_result(
    cox_css_primary,
    "Primary multivariable model",
    "CSS"
  ),

  key_result(
    cox_css_hba1c,
    "+ HbA1c",
    "CSS"
  ),

  key_result(
    cox_css_adj,
    "+ Adjuvant chemotherapy",
    "CSS"
  ),

  key_result(
    cox_css_age,
    "+ continuous age",
    "CSS"
  ),

  key_result(
    cox_css_year,
    "+ continuous surgery year",
    "CSS"
  ),

  key_result(
    cox_css_logca,
    "Using log2-transformed CA19-9",
    "CSS"
  ),

  key_result(
    cox_css_logca_adj,
    "Using log2-transformed CA19-9 + Adjuvant chemotherapy",
    "CSS"
  ),

  key_result(
    cox_css_noNAC,
    "Excluding neoadjuvant chemotherapy",
    "CSS"
  ),

  key_result(
    cox_css_hba1c_replace,
    "Replacing TAD with HbA1c",
    "CSS",
    "^HbA1c$"
  )
)


cat("\n================ SUPPLEMENTARY TABLE S4: OS ================\n")
print(results_os)

cat("\n================ SUPPLEMENTARY TABLE S3: CSS ================\n")
print(results_css)

# Separate files are written deliberately to prevent accidental
# interchange of the S3 (CSS) and S4 (OS) results.
write.csv(
  results_css,
  "supplementary_table_S3_css_sensitivity.csv",
  row.names = FALSE
)

write.csv(
  results_os,
  "supplementary_table_S4_os_sensitivity.csv",
  row.names = FALSE
)

# Combined results file.
results <- bind_rows(results_os, results_css)

write.csv(
  results,
  "survival_sensitivity_results_all.csv",
  row.names = FALSE
)


# ============================================================
# FOUR DIABETES PHENOTYPES
#
# Fig. 1c / Supplementary Table S5A
#
# Reference: no diabetes
# ============================================================

cox_os_dm4 <- coxph(
  Surv(os_time, os_event) ~
    dm4 + stage + margin + NLR + CA199,
  data = df
)

cat(
  "\n================ FOUR DIABETES PHENOTYPES ================\n"
)

print(tidy_cox(cox_os_dm4))

# Reference HRs reported in the manuscript (no diabetes as reference):
#
# Stable diabetes:
#   HR 0.81 (0.43-1.52), P = 0.504
#
# New-onset diabetes:
#   HR 1.70 (1.08-2.69), P = 0.023
#
# Worsening diabetes:
#   HR 1.86 (1.08-3.22), P = 0.026


# ============================================================
# DIABETES-ONLY ANALYSES
#
# n = 81
# TAD = new-onset or worsening diabetes
# Reference = long-standing stable diabetes
#
# Supplementary Table S5A
# ============================================================

df_dm <- df %>%
  filter(!is.na(dm_compare))


stopifnot(nrow(df_dm) == 81)

stopifnot(
  sum(df_dm$dm_compare == "Stable diabetes") == 17
)

stopifnot(
  sum(df_dm$dm_compare == "TAD") == 64
)

stopifnot(
  sum(
    df_dm$os_event[
      df_dm$dm_compare == "Stable diabetes"
    ]
  ) == 12
)

stopifnot(
  sum(
    df_dm$os_event[
      df_dm$dm_compare == "TAD"
    ]
  ) == 52
)


# Primary diabetes-only model

cox_dm_primary <- coxph(
  Surv(os_time, os_event) ~
    dm_compare + stage + margin + NLR + CA199,
  data = df_dm
)

cat("\n-- Diabetes only: primary model --\n")
print(tidy_cox(cox_dm_primary))

# Reference result reported in the manuscript:
# HR 2.34 (1.20-4.54), P = 0.012


# + HbA1c

cox_dm_hba1c <- coxph(
  Surv(os_time, os_event) ~
    dm_compare + stage + margin + NLR +
    CA199 + HbA1c,
  data = df_dm
)

cat("\n-- Diabetes only: + HbA1c --\n")
print(tidy_cox(cox_dm_hba1c))

# Reference result reported in the manuscript:
# HR 2.51 (1.26-5.01), P = 0.009


# Reduced sensitivity model:
# stage + margin only

cox_dm_reduced <- coxph(
  Surv(os_time, os_event) ~
    dm_compare + stage + margin,
  data = df_dm
)

cat("\n-- Diabetes only: reduced model --\n")
print(tidy_cox(cox_dm_reduced))

# Reference result reported in the manuscript:
# HR 2.02 (1.07-3.82), P = 0.031


# ============================================================
# STAGE-SPECIFIC ANALYSES
#
# Supplementary Fig. S2
# ============================================================

df_low <- df %>%
  filter(stage == "<=IIa")

df_high <- df %>%
  filter(stage == ">=IIb")


# Verify sample sizes and event counts reported in S2.

stopifnot(nrow(df_low) == 43)
stopifnot(sum(df_low$os_event) == 24)
stopifnot(sum(df_low$css_event) == 20)

stopifnot(nrow(df_high) == 119)
stopifnot(sum(df_high$os_event) == 96)
stopifnot(sum(df_high$css_event) == 90)


# Stage <=IIa

cox_os_low <- coxph(
  Surv(os_time, os_event) ~
    TAD + margin + NLR + CA199,
  data = df_low
)

cox_css_low <- coxph(
  Surv(os_time, css_event) ~
    TAD + margin + NLR + CA199,
  data = df_low
)


# Stage >=IIb

cox_os_high <- coxph(
  Surv(os_time, os_event) ~
    TAD + margin + NLR + CA199,
  data = df_high
)

cox_css_high <- coxph(
  Surv(os_time, css_event) ~
    TAD + margin + NLR + CA199,
  data = df_high
)


cat("\n================ STAGE-SPECIFIC OS ================\n")

cat("\n-- <=IIa --\n")
print(extract_term(cox_os_low, "^TADTAD$"))

cat("\n-- >=IIb --\n")
print(extract_term(cox_os_high, "^TADTAD$"))


cat("\n================ STAGE-SPECIFIC CSS ================\n")

cat("\n-- <=IIa --\n")
print(extract_term(cox_css_low, "^TADTAD$"))

cat("\n-- >=IIb --\n")
print(extract_term(cox_css_high, "^TADTAD$"))


# ============================================================
# TAD-BY-STAGE INTERACTIONS
#
# Supplementary Table S5B
# Wald P value for the interaction coefficient
# ============================================================

cox_os_interaction <- coxph(
  Surv(os_time, os_event) ~
    TAD * stage + margin + NLR + CA199,
  data = df
)

cox_css_interaction <- coxph(
  Surv(os_time, css_event) ~
    TAD * stage + margin + NLR + CA199,
  data = df
)


os_int <- tidy_cox(cox_os_interaction) %>%
  filter(grepl("TAD.*:stage|stage.*:TAD", term))

css_int <- tidy_cox(cox_css_interaction) %>%
  filter(grepl("TAD.*:stage|stage.*:TAD", term))


cat("\n================ TAD x STAGE INTERACTION ================\n")

cat("\n-- OS --\n")
print(os_int)

cat("\n-- CSS --\n")
print(css_int)

# Reference interaction P values reported in the manuscript:
# OS:  0.200
# CSS: 0.386


# ============================================================
# Session information
# ============================================================

sink("sessionInfo_survival.txt")
print(sessionInfo())
sink()


cat("\nAll survival analyses completed successfully.\n")