#!/usr/bin/env Rscript
# ============================================================
# 01_survival_analysis.R
#
# Clinically defined tumor-associated diabetes, intratumoral
# T-cell densities, and survival in resected PDAC
#
# Current manuscript mapping:
#   Fig. 1a                  OS KM: TAD vs non-TAD
#   Fig. 1b                  OS KM: four diabetes phenotypes
#                            + adjusted HRs from Cox model
#   Supplementary Fig. S1    CSS KM: TAD vs non-TAD
#   Supplementary Fig. S2    Stage-stratified adjusted HRs
#                            (FOREST PLOT; not a KM figure)
#   Supplementary Table S2   OS sensitivity analyses
#   Supplementary Table S3   CSS sensitivity analyses
#   Supplementary Table S4A  Four-category OS + DM-only comparison
#   Supplementary Table S4C  TAD x stage interaction tests
#
# Primary multivariable model:
#   TAD + stage + margin + NLR + CA19-9 (>37 vs <=37 U/mL)
#
# Sensitivity models:
#   + continuous HbA1c
#   + adjuvant chemotherapy
#   CA19-9 as log2-continuous
#   log2-continuous CA19-9 + adjuvant chemotherapy
#   excluding neoadjuvant-chemotherapy patients
#
# Expected cohort:
#   n=162; non-TAD=98; TAD=64
#   no DM=81; stable DM=17; new-onset DM=40; worsening DM=24
#   stage >=IIb=119
#   excluding neoadjuvant chemotherapy=144
#   OS deaths=120; PDAC deaths=110
#
# NOTE: clinical_data.xlsx is not distributed publicly.
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readxl)
  library(broom)
})

clin_file <- "clinical_data.xlsx"
out_dir <- "results/clinical_survival"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(clin_file))

raw <- read_excel(clin_file, sheet = 1)

required_cols <- c(
  "os_time", "os_event", "css_event",
  "TAD", "dm_group",
  "stage", "margin", "NLR",
  "CA199", "CA199_raw",
  "HbA1c",
  "adj_chemo", "neoadj"
)

missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

df <- raw %>%
  transmute(
    os_time = as.numeric(os_time),
    os_event = as.integer(os_event),
    css_event = as.integer(css_event),

    TAD = factor(
      as.integer(TAD),
      levels = c(0, 1),
      labels = c("non-TAD", "TAD")
    ),

    dm_group = factor(
      as.integer(dm_group),
      levels = c(0, 1, 2, 3),
      labels = c(
        "No diabetes",
        "Long-standing stable diabetes",
        "New-onset diabetes",
        "Worsening diabetes"
      )
    ),

    stage = factor(
      as.integer(stage),
      levels = c(0, 1),
      labels = c("<=IIa", ">=IIb")
    ),

    margin = factor(
      as.integer(margin),
      levels = c(0, 1),
      labels = c("R0", "R1")
    ),

    NLR = factor(
      as.integer(NLR),
      levels = c(0, 1),
      labels = c("<2.1", ">=2.1")
    ),

    CA199 = factor(
      as.integer(CA199),
      levels = c(0, 1),
      labels = c("<=37", ">37")
    ),

    CA199_raw = as.numeric(CA199_raw),
    HbA1c = as.numeric(HbA1c),

    adj_chemo = factor(
      as.integer(adj_chemo),
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),

    neoadj = factor(
      as.integer(neoadj),
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

# ------------------------------------------------------------
# Cohort validation
# ------------------------------------------------------------

primary_cols <- c(
  "os_time", "os_event", "css_event",
  "TAD", "stage", "margin", "NLR", "CA199"
)

stopifnot(nrow(df) == 162)
stopifnot(all(complete.cases(df[, primary_cols])))
stopifnot(sum(df$TAD == "non-TAD") == 98)
stopifnot(sum(df$TAD == "TAD") == 64)
stopifnot(sum(df$os_event) == 120)
stopifnot(sum(df$css_event) == 110)

stopifnot(sum(df$dm_group == "No diabetes") == 81)
stopifnot(sum(df$dm_group == "Long-standing stable diabetes") == 17)
stopifnot(sum(df$dm_group == "New-onset diabetes") == 40)
stopifnot(sum(df$dm_group == "Worsening diabetes") == 24)

stopifnot(
  all(
    (df$dm_group %in% c("New-onset diabetes", "Worsening diabetes")) ==
      (df$TAD == "TAD")
  )
)

if (any(is.na(df$CA199_raw)) || any(df$CA199_raw <= 0)) {
  stop("CA199_raw must be non-missing and >0 for the log2-continuous analysis.")
}

df$CA199_log2 <- log2(df$CA199_raw)

df_noNAC <- df %>% filter(neoadj == "No")
stopifnot(nrow(df_noNAC) == 144)

df_dm <- df %>%
  filter(dm_group != "No diabetes") %>%
  mutate(
    dm_binary = factor(
      ifelse(
        dm_group == "Long-standing stable diabetes",
        "Stable diabetes",
        "TAD"
      ),
      levels = c("Stable diabetes", "TAD")
    )
  )
stopifnot(nrow(df_dm) == 81)

# ------------------------------------------------------------
# Reverse Kaplan-Meier median follow-up
# ------------------------------------------------------------

rev_km <- survfit(Surv(os_time, 1 - os_event) ~ 1, data = df)
followup_median <- unname(summary(rev_km)$table["median"])

cat("\nMedian follow-up (reverse KM):", followup_median, "months\n")
# Expected manuscript value: 71.7 months

# ------------------------------------------------------------
# KM / log-rank tests used in current figures
# ------------------------------------------------------------

logrank_p <- function(formula, data) {
  x <- survdiff(formula, data = data)
  1 - pchisq(x$chisq, df = length(x$n) - 1)
}

km_results <- data.frame(
  Analysis = c(
    "Fig1a_OS_TAD_vs_nonTAD",
    "Fig1b_OS_four_categories",
    "SupplFigS1_CSS_TAD_vs_nonTAD"
  ),
  P_value = c(
    logrank_p(Surv(os_time, os_event) ~ TAD, df),
    logrank_p(Surv(os_time, os_event) ~ dm_group, df),
    logrank_p(Surv(os_time, css_event) ~ TAD, df)
  )
)

print(km_results)
write.csv(
  km_results,
  file.path(out_dir, "KM_logrank_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Cox helpers
# ------------------------------------------------------------

fit_cox <- function(formula, data, endpoint, model_name, ph_check = TRUE) {
  fit <- coxph(formula, data = data, ties = "efron", x = TRUE)

  out <- tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      endpoint = endpoint,
      model = model_name,
      n = fit$n,
      events = fit$nevent,
      .before = 1
    )

  if (ph_check) {
    zph <- cox.zph(fit)
    write.csv(
      data.frame(
        term = rownames(zph$table),
        zph$table,
        check.names = FALSE
      ),
      file.path(
        out_dir,
        paste0("PH_", endpoint, "_", model_name, ".csv")
      ),
      row.names = FALSE
    )
  }

  list(fit = fit, table = out)
}

extract_term <- function(x, pattern) {
  x$table %>% filter(grepl(pattern, term))
}

# ------------------------------------------------------------
# Primary OS / CSS
# ------------------------------------------------------------

os_primary <- fit_cox(
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199,
  df, "OS", "primary"
)

css_primary <- fit_cox(
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199,
  df, "CSS", "primary"
)

# Expected TAD estimates:
# OS  HR 1.84 (1.25-2.72), P=0.002
# CSS HR 1.87 (1.25-2.81), P=0.002

# ------------------------------------------------------------
# Sensitivity: continuous HbA1c
# ------------------------------------------------------------

df_hba1c <- df %>%
  filter(complete.cases(HbA1c))

os_hba1c <- fit_cox(
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199 + HbA1c,
  df_hba1c, "OS", "plus_continuous_HbA1c"
)

css_hba1c <- fit_cox(
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199 + HbA1c,
  df_hba1c, "CSS", "plus_continuous_HbA1c"
)

# Expected:
# OS  TAD HR 2.18 (1.31-3.62), P=0.003
# CSS TAD HR 2.40 (1.41-4.07), P=0.001

# Collinearity check for the HbA1c-adjusted OS model.
# car::vif() is optional to avoid making the script fail if car is absent.
if (requireNamespace("car", quietly = TRUE)) {
  vif_out <- car::vif(os_hba1c$fit)
  write.csv(
    data.frame(term = names(vif_out), VIF = as.numeric(vif_out)),
    file.path(out_dir, "VIF_OS_plus_HbA1c.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# Sensitivity: + adjuvant chemotherapy
# ------------------------------------------------------------

os_adj <- fit_cox(
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199 + adj_chemo,
  df, "OS", "plus_adjuvant_chemotherapy"
)

css_adj <- fit_cox(
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199 + adj_chemo,
  df, "CSS", "plus_adjuvant_chemotherapy"
)

# Expected:
# OS  TAD HR 1.64 (1.09-2.46), P=0.018
# CSS TAD HR 1.68 (1.10-2.56), P=0.016

# ------------------------------------------------------------
# Sensitivity: log2-continuous CA19-9
# ------------------------------------------------------------

os_ca19log <- fit_cox(
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199_log2,
  df, "OS", "log2_continuous_CA199"
)

css_ca19log <- fit_cox(
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199_log2,
  df, "CSS", "log2_continuous_CA199"
)

# Expected:
# OS  TAD HR 1.86 (1.26-2.74), P=0.002
# CSS TAD HR 1.93 (1.29-2.89), P=0.001

# ------------------------------------------------------------
# Sensitivity: log2 CA19-9 + adjuvant chemotherapy
# ------------------------------------------------------------

os_ca19log_adj <- fit_cox(
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199_log2 + adj_chemo,
  df, "OS", "log2_CA199_plus_adjuvant"
)

css_ca19log_adj <- fit_cox(
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199_log2 + adj_chemo,
  df, "CSS", "log2_CA199_plus_adjuvant"
)

# Expected:
# OS  TAD HR 1.63 (1.09-2.45), P=0.018
# CSS TAD HR 1.71 (1.12-2.61), P=0.012

# ------------------------------------------------------------
# Sensitivity: exclude neoadjuvant chemotherapy
# ------------------------------------------------------------

os_noNAC <- fit_cox(
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199,
  df_noNAC, "OS", "excluding_neoadjuvant"
)

css_noNAC <- fit_cox(
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199,
  df_noNAC, "CSS", "excluding_neoadjuvant"
)

# Expected:
# OS  TAD HR 1.73 (1.13-2.62), P=0.011
# CSS TAD HR 1.78 (1.16-2.73), P=0.009

# ------------------------------------------------------------
# Four-category OS analysis: Supplementary Table S4A / Fig. 1b
# Reference = no diabetes
# ------------------------------------------------------------

fourcat_os <- fit_cox(
  Surv(os_time, os_event) ~ dm_group + stage + margin + NLR + CA199,
  df, "OS", "four_diabetes_categories"
)

# Expected adjusted HRs vs no diabetes:
# Long-standing stable diabetes 0.81 (0.43-1.52), P=0.504
# New-onset diabetes            1.70 (1.08-2.69), P=0.023
# Worsening diabetes            1.86 (1.08-3.22), P=0.026

# ------------------------------------------------------------
# Diabetes-only analysis: TAD vs long-standing stable diabetes
# ------------------------------------------------------------

dm_only <- fit_cox(
  Surv(os_time, os_event) ~ dm_binary + stage + margin + NLR + CA199,
  df_dm, "OS", "diabetes_only_TAD_vs_stable"
)

df_dm_hba1c <- df_dm %>% filter(complete.cases(HbA1c))

dm_only_hba1c <- fit_cox(
  Surv(os_time, os_event) ~ dm_binary + stage + margin + NLR + CA199 + HbA1c,
  df_dm_hba1c, "OS", "diabetes_only_TAD_vs_stable_plus_HbA1c"
)

# Expected:
# TAD vs stable DM HR 2.34 (1.20-4.54), P=0.012
# + HbA1c          HR 2.51 (1.26-5.01), P=0.009

# ------------------------------------------------------------
# Stage interaction and stage-stratified estimates
# Supplementary Fig. S2 / Supplementary Table S4C
# ------------------------------------------------------------

interaction_wald <- function(event, endpoint) {
  f <- as.formula(
    paste0(
      "Surv(os_time, ", event,
      ") ~ TAD * stage + margin + NLR + CA199"
    )
  )
  fit <- coxph(f, data = df, ties = "efron")
  sm <- summary(fit)$coefficients
  interaction_row <- grep("TAD.*:stage|stage.*:TAD", rownames(sm))
  if (length(interaction_row) != 1) {
    stop("Could not uniquely identify TAD-by-stage interaction term.")
  }
  data.frame(
    endpoint = endpoint,
    interaction = "TAD x pathologic stage",
    P_value = sm[interaction_row, "Pr(>|z|)"]
  )
}

stage_fit <- function(event, endpoint, stage_level) {
  d <- df %>% filter(stage == stage_level)
  fit <- coxph(
    as.formula(
      paste0(
        "Surv(os_time, ", event,
        ") ~ TAD + margin + NLR + CA199"
      )
    ),
    data = d,
    ties = "efron"
  )
  tab <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tad <- tab %>% filter(grepl("^TAD", term))
  data.frame(
    endpoint = endpoint,
    stage = stage_level,
    n = nrow(d),
    events = sum(d[[event]]),
    HR = tad$estimate,
    CI_low = tad$conf.low,
    CI_high = tad$conf.high,
    P_value = tad$p.value
  )
}

stage_results <- bind_rows(
  stage_fit("os_event",  "OS",  "<=IIa"),
  stage_fit("os_event",  "OS",  ">=IIb"),
  stage_fit("css_event", "CSS", "<=IIa"),
  stage_fit("css_event", "CSS", ">=IIb")
)

interaction_results <- bind_rows(
  interaction_wald("os_event", "OS"),
  interaction_wald("css_event", "CSS")
)

# Expected stage-stratified results:
# OS  <=IIa HR 2.06 (0.76-5.62), P=0.156
# OS  >=IIb HR 2.23 (1.44-3.46), P<0.001
# CSS <=IIa HR 1.98 (0.66-5.99), P=0.225
# CSS >=IIb HR 2.16 (1.37-3.39), P<0.001
# Interaction P: OS 0.200; CSS 0.386

write.csv(
  stage_results,
  file.path(out_dir, "Supplementary_Figure_S2_stage_stratified_HRs.csv"),
  row.names = FALSE
)

write.csv(
  interaction_results,
  file.path(out_dir, "Supplementary_Table_S4C_stage_interactions.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Save all Cox results
# ------------------------------------------------------------

all_models <- bind_rows(
  os_primary$table,
  css_primary$table,
  os_hba1c$table,
  css_hba1c$table,
  os_adj$table,
  css_adj$table,
  os_ca19log$table,
  css_ca19log$table,
  os_ca19log_adj$table,
  css_ca19log_adj$table,
  os_noNAC$table,
  css_noNAC$table,
  fourcat_os$table,
  dm_only$table,
  dm_only_hba1c$table
)

write.csv(
  all_models,
  file.path(out_dir, "clinical_cox_all_reported_models.csv"),
  row.names = FALSE
)

# Compact manuscript-check table: TAD / DM terms only
manuscript_check <- all_models %>%
  filter(
    grepl("^TAD", term) |
      grepl("^dm_group", term) |
      grepl("^dm_binary", term)
  )

write.csv(
  manuscript_check,
  file.path(out_dir, "clinical_manuscript_check.csv"),
  row.names = FALSE
)

sink(file.path(out_dir, "sessionInfo_survival.txt"))
sessionInfo()
sink()

cat("\n============================================\n")
cat("FINAL MANUSCRIPT CHECK\n")
cat("============================================\n")
cat("n =", nrow(df), "\n")
cat("OS events =", sum(df$os_event), "\n")
cat("CSS events =", sum(df$css_event), "\n")
cat("Median follow-up =", followup_median, "months\n\n")
print(manuscript_check)
cat("\nStage-stratified estimates:\n")
print(stage_results)
cat("\nInteraction tests:\n")
print(interaction_results)
cat("============================================\n")
