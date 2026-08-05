# ============================================================
# 01_survival_analysis.R
#
# Clinical survival analyses for the current TAD-PDAC manuscript.
#
# Computes the statistics underlying:
#   Fig. 1a                    OS, TAD vs non-TAD
#   Fig. 1b                    OS, four diabetes phenotypes
#   Supplementary Fig. S1      CSS, TAD vs non-TAD
#   Supplementary Fig. S2      stage-stratified adjusted HRs
#   Supplementary Tables S2-S4 survival/sensitivity analyses
#
# Plotting code is intentionally not included.
#
# Current cohort checks:
#   total n = 162
#   non-TAD = 98; TAD = 64
#   no DM = 81; stable DM = 17; new-onset = 40; worsening = 24
#   stage <=IIa = 43; stage >=IIb = 119
#   OS events = 120; CSS events = 110
#   stage <=IIa: OS events = 24; CSS events = 20
#   stage >=IIb: OS events = 96; CSS events = 90
#   neoadjuvant-excluded cohort = 144
#   diabetes-only subset = 81
#
# CSS uses the same follow-up time as OS; deaths from causes
# other than PDAC are censored at the date of death.
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readxl)
  library(broom)
})

clin_file <- "data/clinical_data.xlsx"
out_dir <- "results/clinical_survival"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(clin_file))

raw <- read_excel(clin_file, sheet = 1)

required_cols <- c(
  "os_time", "os_event", "css_event",
  "TAD", "dm_group", "stage", "margin",
  "NLR", "CA199", "CA199_raw", "HbA1c",
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
        "Long-standing stable",
        "New-onset",
        "Worsening"
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
stopifnot(all(complete.cases(df[primary_cols])))
stopifnot(all(df$CA199_raw > 0, na.rm = TRUE))

stopifnot(sum(df$TAD == "non-TAD") == 98)
stopifnot(sum(df$TAD == "TAD") == 64)

stopifnot(sum(df$dm_group == "No diabetes") == 81)
stopifnot(sum(df$dm_group == "Long-standing stable") == 17)
stopifnot(sum(df$dm_group == "New-onset") == 40)
stopifnot(sum(df$dm_group == "Worsening") == 24)

stopifnot(
  all(
    (df$dm_group %in% c("New-onset", "Worsening")) ==
      (df$TAD == "TAD")
  )
)

stopifnot(sum(df$os_event) == 120)
stopifnot(sum(df$css_event) == 110)

df$CA199_log2 <- log2(df$CA199_raw)

df_stage_low <- df %>% filter(stage == "<=IIa")
df_stage_high <- df %>% filter(stage == ">=IIb")
df_noNAC <- df %>% filter(neoadj == "No")

stopifnot(nrow(df_stage_low) == 43)
stopifnot(nrow(df_stage_high) == 119)
stopifnot(sum(df_stage_low$os_event) == 24)
stopifnot(sum(df_stage_low$css_event) == 20)
stopifnot(sum(df_stage_high$os_event) == 96)
stopifnot(sum(df_stage_high$css_event) == 90)
stopifnot(nrow(df_noNAC) == 144)

df_hba1c <- df %>% filter(!is.na(HbA1c))
cat("HbA1c-evaluable cohort: n =", nrow(df_hba1c), "\n")

df_dm <- df %>%
  filter(dm_group != "No diabetes") %>%
  mutate(
    dm_group = droplevels(dm_group),
    dm_binary = factor(
      ifelse(dm_group == "Long-standing stable", "Stable DM", "TAD"),
      levels = c("Stable DM", "TAD")
    )
  )

stopifnot(nrow(df_dm) == 81)

df_dm_hba1c <- df_dm %>% filter(!is.na(HbA1c))
cat("Diabetes-only HbA1c-evaluable cohort: n =", nrow(df_dm_hba1c), "\n")

# ------------------------------------------------------------
# Median follow-up: reverse Kaplan-Meier
# ------------------------------------------------------------

rev_km <- survfit(
  Surv(os_time, 1 - os_event) ~ 1,
  data = df
)

rev_tab <- summary(rev_km)$table

followup_summary <- data.frame(
  n = nrow(df),
  median_followup_months = unname(rev_tab["median"]),
  lower_95 = unname(rev_tab["0.95LCL"]),
  upper_95 = unname(rev_tab["0.95UCL"])
)

write.csv(
  followup_summary,
  file.path(out_dir, "reverse_KM_followup.csv"),
  row.names = FALSE
)

cat("\nReverse-KM follow-up:\n")
print(followup_summary)
cat("Current manuscript median follow-up: approximately 71.7 months.\n")

# ------------------------------------------------------------
# Kaplan-Meier / log-rank statistics
# ------------------------------------------------------------

logrank_result <- function(formula, data, label) {
  fit <- survdiff(formula, data = data)
  df_lr <- length(fit$n) - 1

  data.frame(
    Analysis = label,
    n = nrow(data),
    Chi_square = unname(fit$chisq),
    df = df_lr,
    P_value = 1 - pchisq(fit$chisq, df = df_lr),
    stringsAsFactors = FALSE
  )
}

logrank_results <- bind_rows(
  logrank_result(
    Surv(os_time, os_event) ~ TAD,
    df,
    "Fig1a_OS_TAD_vs_nonTAD"
  ),
  logrank_result(
    Surv(os_time, css_event) ~ TAD,
    df,
    "SuppFigS1_CSS_TAD_vs_nonTAD"
  ),
  logrank_result(
    Surv(os_time, os_event) ~ dm_group,
    df,
    "Fig1b_OS_four_diabetes_categories"
  )
)

write.csv(
  logrank_results,
  file.path(out_dir, "KM_logrank_statistics.csv"),
  row.names = FALSE
)

cat("\nKM / log-rank statistics:\n")
print(logrank_results)

# ------------------------------------------------------------
# Cox helper
# ------------------------------------------------------------

cox_tables <- list()
ph_tables <- list()

fit_cox <- function(formula, data, endpoint, model_name) {
  fit <- coxph(formula, data = data, x = TRUE)

  tab <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      Endpoint = endpoint,
      Model = model_name,
      n = fit$n,
      Events = fit$nevent,
      .before = 1
    )

  zph <- cox.zph(fit)
  global_p <- unname(zph$table["GLOBAL", "p"])

  ph <- data.frame(
    Endpoint = endpoint,
    Model = model_name,
    n = fit$n,
    Events = fit$nevent,
    Global_PH_P = global_p,
    stringsAsFactors = FALSE
  )

  cat(
    endpoint, "|", model_name,
    "| n =", fit$n,
    "| events =", fit$nevent,
    "| cox.zph GLOBAL P =", signif(global_p, 4), "\n"
  )

  list(fit = fit, table = tab, ph = ph)
}

add_cox <- function(name, formula, data, endpoint, model_name) {
  z <- fit_cox(formula, data, endpoint, model_name)
  cox_tables[[name]] <<- z$table
  ph_tables[[name]] <<- z$ph
  invisible(z$fit)
}

# ------------------------------------------------------------
# Primary multivariable models
# ------------------------------------------------------------

add_cox(
  "OS_primary",
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199,
  df, "OS", "Primary"
)

add_cox(
  "CSS_primary",
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199,
  df, "CSS", "Primary"
)

# ------------------------------------------------------------
# Sensitivity models
# ------------------------------------------------------------

add_cox(
  "OS_HbA1c",
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199 + HbA1c,
  df_hba1c, "OS", "Primary_plus_HbA1c"
)

add_cox(
  "CSS_HbA1c",
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199 + HbA1c,
  df_hba1c, "CSS", "Primary_plus_HbA1c"
)

add_cox(
  "OS_adjuvant",
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199 + adj_chemo,
  df, "OS", "Primary_plus_adjuvant"
)

add_cox(
  "CSS_adjuvant",
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199 + adj_chemo,
  df, "CSS", "Primary_plus_adjuvant"
)

add_cox(
  "OS_logCA199",
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199_log2,
  df, "OS", "log2_continuous_CA199"
)

add_cox(
  "CSS_logCA199",
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199_log2,
  df, "CSS", "log2_continuous_CA199"
)

add_cox(
  "OS_logCA199_adjuvant",
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR +
    CA199_log2 + adj_chemo,
  df, "OS", "log2_CA199_plus_adjuvant"
)

add_cox(
  "CSS_logCA199_adjuvant",
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR +
    CA199_log2 + adj_chemo,
  df, "CSS", "log2_CA199_plus_adjuvant"
)

add_cox(
  "OS_noNAC",
  Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199,
  df_noNAC, "OS", "Exclude_neoadjuvant"
)

add_cox(
  "CSS_noNAC",
  Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199,
  df_noNAC, "CSS", "Exclude_neoadjuvant"
)

# Four diabetes phenotypes; no diabetes is the reference.
add_cox(
  "OS_four_category",
  Surv(os_time, os_event) ~ dm_group + stage + margin + NLR + CA199,
  df, "OS", "Four_diabetes_categories"
)

# Diabetes-only: TAD vs long-standing stable diabetes.
add_cox(
  "OS_diabetes_only",
  Surv(os_time, os_event) ~ dm_binary + stage + margin + NLR + CA199,
  df_dm, "OS", "Diabetes_only"
)

add_cox(
  "OS_diabetes_only_HbA1c",
  Surv(os_time, os_event) ~ dm_binary + stage + margin + NLR + CA199 + HbA1c,
  df_dm_hba1c, "OS", "Diabetes_only_plus_HbA1c"
)

# ------------------------------------------------------------
# Supplementary Fig. S2: stage-stratified adjusted HRs
# Stage is constant within each subgroup and is therefore not
# included in the subgroup-specific Cox models.
# ------------------------------------------------------------

stage_stratified <- bind_rows(
  fit_cox(
    Surv(os_time, os_event) ~ TAD + margin + NLR + CA199,
    df_stage_low, "OS", "Stage_<=IIa"
  )$table %>% filter(grepl("^TAD", term)),

  fit_cox(
    Surv(os_time, os_event) ~ TAD + margin + NLR + CA199,
    df_stage_high, "OS", "Stage_>=IIb"
  )$table %>% filter(grepl("^TAD", term)),

  fit_cox(
    Surv(os_time, css_event) ~ TAD + margin + NLR + CA199,
    df_stage_low, "CSS", "Stage_<=IIa"
  )$table %>% filter(grepl("^TAD", term)),

  fit_cox(
    Surv(os_time, css_event) ~ TAD + margin + NLR + CA199,
    df_stage_high, "CSS", "Stage_>=IIb"
  )$table %>% filter(grepl("^TAD", term))
)

write.csv(
  stage_stratified,
  file.path(out_dir, "Supplementary_Figure_S2_stage_stratified_HRs.csv"),
  row.names = FALSE
)

# Current expected values:
# OS  <=IIa: HR 2.06 (0.76-5.62), P=0.156
# OS  >=IIb: HR 2.23 (1.44-3.46), P<0.001
# CSS <=IIa: HR 1.98 (0.66-5.99), P=0.225
# CSS >=IIb: HR 2.16 (1.37-3.39), P<0.001

# ------------------------------------------------------------
# TAD-by-stage interaction:
# Reported P is the WALD test for the interaction coefficient.
# ------------------------------------------------------------

wald_interaction <- function(event_var, endpoint) {
  f <- as.formula(
    paste0(
      "Surv(os_time, ", event_var, ") ~ ",
      "TAD * stage + margin + NLR + CA199"
    )
  )

  fit <- coxph(f, data = df, x = TRUE)
  sm <- summary(fit)$coefficients

  int_row <- grep(
    "TAD.*:stage|stage.*:TAD",
    rownames(sm),
    value = TRUE
  )

  if (length(int_row) != 1) {
    stop("Could not uniquely identify TAD-by-stage term for ", endpoint)
  }

  data.frame(
    Endpoint = endpoint,
    Test = "Wald test for TAD-by-stage interaction term",
    Interaction_term = int_row,
    Wald_z = sm[int_row, "z"],
    P_value = sm[int_row, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

interaction_results <- bind_rows(
  wald_interaction("os_event", "OS"),
  wald_interaction("css_event", "CSS")
)

write.csv(
  interaction_results,
  file.path(out_dir, "Supplementary_Table_S4C_survival_stage_interaction_Wald.csv"),
  row.names = FALSE
)

cat("\nStage interaction (Wald):\n")
print(interaction_results)
cat("Expected: OS P ~0.200; CSS P ~0.386.\n")

# ------------------------------------------------------------
# Save Cox and PH results
# ------------------------------------------------------------

cox_table <- bind_rows(cox_tables)
ph_table <- bind_rows(ph_tables)

write.csv(
  cox_table,
  file.path(out_dir, "clinical_cox_models.csv"),
  row.names = FALSE
)

write.csv(
  ph_table,
  file.path(out_dir, "cox_zph_global_tests.csv"),
  row.names = FALSE
)

cat("\nPrimary-model PH global tests:\n")
print(ph_table %>% filter(Model == "Primary"))

cat("\n============================================\n")
cat("FINAL CLINICAL SURVIVAL CHECK\n")
cat("============================================\n")
cat("Total n:", nrow(df), "\n")
cat("OS events:", sum(df$os_event), "\n")
cat("CSS events:", sum(df$css_event), "\n")
cat(
  "Stage <=IIa: n =", nrow(df_stage_low),
  "| OS events =", sum(df_stage_low$os_event),
  "| CSS events =", sum(df_stage_low$css_event), "\n"
)
cat(
  "Stage >=IIb: n =", nrow(df_stage_high),
  "| OS events =", sum(df_stage_high$os_event),
  "| CSS events =", sum(df_stage_high$css_event), "\n"
)
cat("HbA1c model n:", nrow(df_hba1c), "\n")
cat("Diabetes-only n:", nrow(df_dm), "\n")
cat("No-neoadjuvant n:", nrow(df_noNAC), "\n")
cat("============================================\n")

sink(file.path(out_dir, "sessionInfo_survival.txt"))
sessionInfo()
sink()
