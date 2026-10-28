
#!/usr/bin/env Rscript
# ============================================================
# 01_survival_analysis.R
#
# Outputs (results/survival/):
#   Table 2           All 13 univariable OS terms and primary multivariable OS
#   Fig. 1b-c         OS log-rank statistics (TAD; four categories)
#   Fig. 1c           Adjusted HRs for the four diabetes categories
#                     (Supplementary Table S6A)
#   Suppl. Fig. S2    CSS log-rank statistics
#   Suppl. Fig. S3    Stage-specific Cox models
#   Suppl. Table S4   CSS sensitivity models
#   Suppl. Table S5   OS sensitivity models (including CD8 adjustment)
#   Suppl. Table S6A  Four-category and diabetes-only models
#   Suppl. Table S6B  TAD-by-stage interaction, OS/CSS rows
#                     (Wald tests; T-cell rows: 05_tcell_adjusted_models.R)
#   Also saved: primary CSS model, proportional hazards tests,
#   median follow-up, session information
#
# Input: data/clinical_data.xlsx (see DATA_DICTIONARY.md)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(survival)
  library(broom)
})

# Two diabetes-only sensitivity models appear in both S5 and S6A.

out_dir <- "results/survival"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- "data/clinical_data.xlsx"
if (!file.exists(input_file)) {
  stop("Missing controlled-access input: ", input_file)
}
raw <- read_excel(input_file, sheet = 1)

# Additional Table 2 covariates are defined in DATA_DICTIONARY.md.
# Binary variables are coded 0 (reference) or 1.

get_col <- function(candidates) {
  hits <- candidates[candidates %in% names(raw)]
  if (!length(hits)) {
    stop("Missing input column: ", paste(candidates, collapse = "/"))
  }
  raw[[hits[1]]]
}

df <- data.frame(
  os_time = as.numeric(get_col("os_time")),
  os_event = as.integer(get_col("os_event")),
  css_event = as.integer(get_col("css_event")),
  TAD_num = as.integer(get_col("TAD")),
  stage_num = as.integer(get_col("stage")),
  margin_num = as.integer(get_col("margin")),
  NLR_num = as.integer(get_col("NLR")),
  NLR_raw = as.numeric(get_col("NLR_raw")),
  CA199_num = as.integer(get_col("CA199")),
  CA199_raw = as.numeric(get_col("CA199_raw")),
  adj_num = as.integer(get_col("adj_chemo")),
  neoadj_num = as.integer(get_col("neoadj")),
  HbA1c = as.numeric(get_col(c("HbA1c", "hba1c"))),
  age = as.numeric(get_col(c("age", "Age"))),
  surgery_year = as.numeric(get_col(
    c("surgery_year", "Surgery_year")
  )),
  dm_group_num = as.integer(get_col("dm_group")),
  # Defined for worsening diabetes only (dm_group = 3);
  # 2 = worsening identified by treatment intensification alone.
  # Full coding in DATA_DICTIONARY.md.
  worsening_basis = as.integer(get_col("worsening_basis")),
  # Defined for stable diabetes only (dm_group = 1);
  # 1 = longitudinal glycaemic data available.
  stable_longitudinal_data = as.integer(
    get_col("stable_longitudinal_data")
  ),
  BMI = as.numeric(get_col(c("BMI", "bmi"))),
  SMI = as.numeric(get_col(c("SMI", "smi"))),
  albumin = as.numeric(get_col(c("albumin", "Albumin"))),
  ECOG_PS = as.integer(get_col(
    c("ECOG_PS", "ECOG-PS", "ECOGPS")
  )),
  CD8 = as.numeric(get_col(c("CD8", "CD8_density")))
)

df <- df %>%
  mutate(
    TAD = factor(TAD_num, c(0, 1), c("non-TAD", "TAD")),
    stage = factor(stage_num, c(0, 1), c("<=IIa", ">=IIb")),
    margin = factor(margin_num, c(0, 1), c("R0", "R1")),
    NLR = factor(NLR_num, c(0, 1), c("<2.1", ">=2.1")),
    CA199 = factor(CA199_num, c(0, 1), c("<=37", ">37")),
    adj_chemo = factor(adj_num, c(0, 1), c("No", "Yes")),
    neoadj = factor(neoadj_num, c(0, 1), c("No", "Yes")),
    dm4 = factor(
      dm_group_num,
      levels = 0:3,
      labels = c(
        "No diabetes", "Stable diabetes",
        "New-onset diabetes", "Worsening diabetes"
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
    CA199_log2 = log2(CA199_raw),
    CD8_log2 = log2(CD8)
  )

# ------------------------------------------------------------
# Cohort validation
# ------------------------------------------------------------
stopifnot(
  nrow(df) == 162,
  !anyNA(df$dm_group_num),
  sum(df$TAD_num == 0) == 98,
  sum(df$TAD_num == 1) == 64,
  identical(
    as.integer(table(df$dm_group_num)),
    c(81L, 17L, 40L, 24L)
  ),
  sum(df$os_event) == 120,
  sum(df$css_event) == 110,
  all(df$CA199_raw > 0),
  all(df$NLR_raw >= 0),
  all(is.finite(df$CD8)),
  all(df$CD8 > 0),
  all(df$TAD_num ==
        as.integer(df$dm_group_num %in% c(2, 3)))
)

is_worsening <- df$dm_group_num == 3
is_stable <- df$dm_group_num == 1
wb <- df$worsening_basis
sld <- df$stable_longitudinal_data

stopifnot(
  # worsening_basis: coded 1-3 for all 24 worsening cases only
  all(wb[is_worsening] %in% 1:3),
  !any(wb[!is_worsening] %in% 1:3),
  sum(wb[is_worsening] == 1) == 9,
  sum(wb[is_worsening] == 2) == 11,
  sum(wb[is_worsening] == 3) == 4,
  # stable_longitudinal_data: 0/1 for all 17 stable cases;
  # value 1 occurs only within stable diabetes
  all(sld[is_stable] %in% 0:1),
  sum(sld[is_stable] == 1) == 8,
  sum(sld[is_stable] == 0) == 9,
  sum(sld == 1, na.rm = TRUE) == 8
)

model_cols <- c(
  "os_time", "os_event", "css_event", "TAD",
  "stage", "margin", "NLR", "CA199",
  "CA199_log2", "NLR_raw", "HbA1c",
  "age", "surgery_year", "adj_chemo",
  "BMI", "SMI", "albumin", "ECOG_hi", "CD8_log2"
)
stopifnot(all(complete.cases(df[model_cols])))

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
endpoint_label <- c(os_event = "OS", css_event = "CSS")

cox_fit <- function(endpoint, rhs, data = df) {
  f <- as.formula(paste0(
    "Surv(os_time, ", endpoint, ") ~ ", rhs
  ))
  coxph(f, data = data, ties = "efron")
}

tidy_cox <- function(fit) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
}

# One coefficient per model. Term identifies the estimand:
#   TADTAD          TAD vs non-TAD
#   HbA1c           per 1% increase in HbA1c
#   dm_compareTAD   TAD vs stable diabetes (diabetes-only models)
model_result <- function(fit, label, endpoint,
                         pattern = "^TADTAD$") {
  x <- tidy_cox(fit)
  x <- x[grepl(pattern, x$term), , drop = FALSE]
  if (nrow(x) != 1) {
    stop("Expected exactly one coefficient in: ", label)
  }
  data.frame(
    Endpoint = endpoint_label[[endpoint]],
    Model = label,
    Term = x$term,
    N = fit$n,
    Events = fit$nevent,
    HR = x$estimate,
    CI_low = x$conf.low,
    CI_high = x$conf.high,
    P_value = x$p.value
  )
}

group_counts <- function(data, label, endpoint = "os_event") {
  data %>%
    group_by(Group = dm_compare) %>%
    summarise(
      N = n(),
      Events = sum(.data[[endpoint]]),
      .groups = "drop"
    ) %>%
    mutate(
      Endpoint = endpoint_label[[endpoint]],
      Model = label,
      .before = 1
    )
}

ph_table <- function(fit, label) {
  z <- cox.zph(fit)$table
  data.frame(
    Model = label,
    Term = rownames(z),
    Chi_square = z[, "chisq"],
    df = z[, "df"],
    P_value = z[, "p"],
    row.names = NULL
  )
}

save_csv <- function(x, name) {
  write.csv(
    x, file.path(out_dir, name),
    row.names = FALSE
  )
}

base_rhs <- "TAD + stage + margin + NLR + CA199"
base_without_tad <- "stage + margin + NLR + CA199"

# ------------------------------------------------------------
# Primary OS/CSS models
# ------------------------------------------------------------
cox_os_primary <- cox_fit("os_event", base_rhs)
cox_css_primary <- cox_fit("css_event", base_rhs)

cat("Primary OS model:\n")
print(tidy_cox(cox_os_primary))
cat("Primary CSS model:\n")
print(tidy_cox(cox_css_primary))

# ------------------------------------------------------------
# Table 2: ALL 13 univariable covariates
# ------------------------------------------------------------
# Input columns for additional Table 2 terms are described in
# DATA_DICTIONARY.md. Missing or invalid values yield NA estimates.

optional_col <- function(candidates) {
  matches <- candidates[candidates %in% names(raw)]
  if (!length(matches)) return(NULL)
  raw[[matches[1]]]
}

strict_binary <- function(x, variable) {
  if (is.null(x)) return(NULL)
  values <- trimws(as.character(x))
  values[values == ""] <- NA_character_
  if (anyNA(values) || !all(values %in% c("0", "1"))) {
    warning(variable, ": expected complete 0/1 coding; not analysed.")
    return(NULL)
  }
  as.integer(values)
}

optional_category <- function(binary_candidates, numeric_candidates,
                              threshold, variable) {
  binary <- strict_binary(optional_col(binary_candidates), variable)
  if (!is.null(binary)) return(binary)
  numeric_values <- optional_col(numeric_candidates)
  if (is.null(numeric_values)) return(NULL)
  x <- suppressWarnings(as.numeric(numeric_values))
  if (length(x) != nrow(df) || any(!is.finite(x)) || any(x < 0)) {
    warning(variable, ": incomplete/invalid numeric input; not analysed.")
    return(NULL)
  }
  as.integer(x > threshold)
}

df$age_ge65 <- factor(as.integer(df$age >= 65), 0:1,
                       c("<65", ">=65"))

male_values <- strict_binary(optional_col(c("male", "male_sex")),
                             "Male sex")
if (is.null(male_values)) {
  sx <- optional_col(c("sex", "Sex"))
  if (!is.null(sx)) {
    sx <- trimws(tolower(as.character(sx)))
    if (!anyNA(sx) && all(sx %in% c("female", "male"))) {
      male_values <- as.integer(sx == "male")
    } else {
      warning("sex must contain Female/Male for this input format.")
    }
  }
}

tumor_values <- optional_category(
  c("tumor_gt20", "tumor_size_gt20"),
  c("tumor_size_mm", "tumor_size_mm_raw"),
  20, "Tumor size >20 mm"
)
cea_values <- optional_category(
  c("CEA_gt5", "CEA_over5"),
  c("CEA_raw", "CEA_ng_mL"),
  5, "CEA >5 ng/mL"
)

extra_data <- list(
  male = male_values,
  tumor_gt20 = tumor_values,
  KRAS = strict_binary(optional_col(
    c("KRAS_mutation", "KRAS_mut", "KRAS_any")), "KRAS mutation"),
  TP53 = strict_binary(optional_col(
    c("TP53_alteration", "TP53_alt", "TP53_any")), "TP53 alteration"),
  SMAD4 = strict_binary(optional_col(
    c("SMAD4_alteration", "SMAD4_alt", "SMAD4_any")), "SMAD4 alteration"),
  CDKN2A = strict_binary(optional_col(
    c("CDKN2A_alteration", "CDKN2A_alt", "CDKN2A_any")),
    "CDKN2A alteration"),
  CEA_gt5 = cea_values
)

for (nm in names(extra_data)) {
  x <- extra_data[[nm]]
  if (!is.null(x)) {
    if (length(x) != nrow(df) || length(unique(x)) < 2L) {
      warning(nm, ": invalid length or no category variation; not analysed.")
    } else {
      df[[nm]] <- factor(x, levels = 0:1, labels = c("No", "Yes"))
    }
  }
}

# Table 2 variables (manuscript order).
# In multivariable Cox regression, only five prespecified terms enter.
table2_specs <- data.frame(
  Variable = c(
    "Tumor-associated diabetes (TAD)", "Age >=65 years", "Male sex",
    "Pathologic stage >=IIb", "R1 resection margin",
    "Tumor size >20 mm", "KRAS mutation", "TP53 alteration",
    "SMAD4 alteration", "CDKN2A alteration", "NLR >=2.1",
    "CEA >5 ng/mL", "CA19-9 >37 U/mL"
  ),
  Column = c(
    "TAD", "age_ge65", "male", "stage", "margin",
    "tumor_gt20", "KRAS", "TP53", "SMAD4", "CDKN2A",
    "NLR", "CEA_gt5", "CA199"
  ),
  stringsAsFactors = FALSE
)

multi_os <- tidy_cox(cox_os_primary)

table2 <- bind_rows(lapply(seq_len(nrow(table2_specs)), function(i) {
  variable <- table2_specs$Variable[i]
  column <- table2_specs$Column[i]
  available <- column %in% names(df)
  if (available) {
    available <- all(complete.cases(df[, c("os_time", "os_event", column)])) &&
      length(unique(df[[column]])) == 2L
  }
  if (available) {
    fit <- cox_fit("os_event", column)
    u <- tidy_cox(fit)
    if (nrow(u) != 1L) {
      stop("Table 2 binary variable produced !=1 coefficient: ", column)
    }
    status <- "Calculated"
    u_hr <- u$estimate
    u_low <- u$conf.low
    u_high <- u$conf.high
    u_p <- u$p.value
  } else {
    status <- "INPUT_MISSING_OR_INVALID"
    u_hr <- u_low <- u_high <- u_p <- NA_real_
    warning("Table 2 univariable not calculated: ", variable,
            ". Supply a validated 0/1 column (or specified raw measurement).")
  }
  primary_terms <- c(TAD = "TADTAD", stage = "stage>=IIb",
                     margin = "marginR1", NLR = "NLR>=2.1",
                     CA199 = "CA199>37")
  target_term <- if (column %in% names(primary_terms)) {
    unname(primary_terms[column])
  } else {
    "_not_in_primary_model_"
  }
  m <- multi_os[multi_os$term == target_term, , drop = FALSE]
  if (column %in% c("TAD", "stage", "margin", "NLR", "CA199") &&
      nrow(m) != 1) stop("Missing primary-model coefficient: ", column)
  data.frame(
    Variable = variable, Input_column = column,
    Status = status,
    Uni_HR = u_hr, Uni_CI_low = u_low,
    Uni_CI_high = u_high, Uni_P = u_p,
    Multi_HR = if (nrow(m)) m$estimate else NA_real_,
    Multi_CI_low = if (nrow(m)) m$conf.low else NA_real_,
    Multi_CI_high = if (nrow(m)) m$conf.high else NA_real_,
    Multi_P = if (nrow(m)) m$p.value else NA_real_,
    stringsAsFactors = FALSE
  )
}))

save_csv(table2, "Table_2_OS_Cox_univariable_multivariable.csv")
cat("\nTable 2 analysis status:\n")
print(table2[, c("Variable", "Status")], row.names = FALSE)

save_csv(
  multi_os %>% mutate(Endpoint = "OS", .before = 1),
  "Primary_OS_Cox_model.csv"
)
save_csv(
  tidy_cox(cox_css_primary) %>% mutate(Endpoint = "CSS", .before = 1),
  "Primary_CSS_Cox_model.csv"
)

# Proportional hazards (scaled Schoenfeld residuals)
ph_results <- bind_rows(
  ph_table(cox_os_primary, "Primary OS model"),
  ph_table(cox_css_primary, "Primary CSS model")
)
print(ph_results)
save_csv(ph_results, "Proportional_hazards_tests.csv")

# Median follow-up (reverse Kaplan-Meier)
reverse_km <- survfit(
  Surv(os_time, 1 - os_event) ~ 1,
  data = df
)
fu <- summary(reverse_km)$table
follow_up <- data.frame(
  Method = "Reverse Kaplan-Meier",
  Median = unname(fu["median"]),
  CI_low = unname(fu["0.95LCL"]),
  CI_high = unname(fu["0.95UCL"])
)
cat("Median follow-up:\n")
print(follow_up)
save_csv(follow_up, "Median_follow_up_reverse_KM.csv")

# ------------------------------------------------------------
# Kaplan-Meier/log-rank statistics
# ------------------------------------------------------------
km_comparison <- function(endpoint, group, name) {
  f <- as.formula(
    paste0("Surv(os_time,", endpoint, ") ~ ", group)
  )
  x <- survdiff(f, data = df)
  data.frame(
    Figure = name,
    Endpoint = endpoint_label[[endpoint]],
    Groups = group,
    Chi_square = unname(x$chisq),
    df = length(x$n) - 1,
    P_value = pchisq(
      x$chisq, length(x$n) - 1,
      lower.tail = FALSE
    )
  )
}

km_results <- bind_rows(
  km_comparison("os_event", "TAD", "Fig. 1b"),
  km_comparison("os_event", "dm4", "Fig. 1c"),
  km_comparison("css_event", "TAD", "Supplementary Fig. S2")
)
save_csv(km_results, "Figure_1b-c_and_Supplementary_Figure_S2_logrank.csv")

# ------------------------------------------------------------
# Sensitivity models (Supplementary Tables S4 and S5)
# ------------------------------------------------------------
sensitivity_specs <- list(
  "Primary multivariable model" = base_rhs,
  "+ HbA1c" = paste(base_rhs, "+ HbA1c"),
  "+ Adjuvant chemotherapy" =
    paste(base_rhs, "+ adj_chemo"),
  "+ continuous age" = paste(base_rhs, "+ age"),
  "+ continuous surgery year" =
    paste(base_rhs, "+ surgery_year"),
  "Using log2-transformed CA19-9" =
    "TAD + stage + margin + NLR + CA199_log2",
  "Using log2-transformed CA19-9 + Adjuvant chemotherapy" =
    "TAD + stage + margin + NLR + CA199_log2 + adj_chemo",
  "Replacing TAD with HbA1c" =
    paste("HbA1c +", base_without_tad)
)

run_sensitivities <- function(endpoint) {
  bind_rows(lapply(names(sensitivity_specs), function(label) {
    fit <- cox_fit(endpoint, sensitivity_specs[[label]])
    term <- if (label == "Replacing TAD with HbA1c") {
      "^HbA1c$"
    } else {
      "^TADTAD$"
    }
    model_result(fit, label, endpoint, term)
  }))
}

results_os <- run_sensitivities("os_event")
results_css <- run_sensitivities("css_event")

# OS only: continuous NLR + continuous CA19-9
results_os <- bind_rows(
  results_os,
  model_result(
    cox_fit(
      "os_event",
      "TAD + stage + margin + NLR_raw + CA199_log2"
    ),
    "Continuous NLR + log2(CA19-9)",
    "os_event"
  )
)

# OS sensitivity model adjusted for log2-transformed CD8+ density.
# Manuscript estimate: HR 1.68 (95% CI 1.14–2.47), P = 0.009.
cox_os_cd8 <- cox_fit(
  "os_event", paste(base_rhs, "+ CD8_log2")
)
cd8_result <- model_result(
  cox_os_cd8,
  "+ log2(intratumoral CD8+ density)",
  "os_event"
)
results_os <- bind_rows(results_os, cd8_result)
save_csv(cd8_result, "Supplementary_Table_S5_CD8_adjusted_OS.csv")

if (round(cd8_result$HR, 2) != 1.68 ||
    round(cd8_result$CI_low, 2) != 1.14 ||
    round(cd8_result$CI_high, 2) != 2.47 ||
    round(cd8_result$P_value, 3) != 0.009) {
  warning("CD8-adjusted OS differs from reference estimates.")
}

# Excluding neoadjuvant-treated patients
df_noNAC <- df %>% filter(neoadj == "No")
stopifnot(nrow(df_noNAC) == 144)

results_os <- bind_rows(
  results_os,
  model_result(
    cox_fit("os_event", base_rhs, df_noNAC),
    "Excluding neoadjuvant chemotherapy",
    "os_event"
  )
)
results_css <- bind_rows(
  results_css,
  model_result(
    cox_fit("css_event", base_rhs, df_noNAC),
    "Excluding neoadjuvant chemotherapy",
    "css_event"
  )
)

# Classification sensitivity: exclude worsening diabetes
# identified by treatment intensification alone (n = 11)
df_no_treatment_only <- df %>%
  filter(!(
    dm_group_num == 3 &
      coalesce(worsening_basis == 2L, FALSE)
  ))
stopifnot(nrow(df_no_treatment_only) == 151)

results_os <- bind_rows(
  results_os,
  model_result(
    cox_fit("os_event", base_rhs, df_no_treatment_only),
    "Excluding treatment-only worsening",
    "os_event"
  )
)

# Body status sensitivity
results_os <- bind_rows(
  results_os,
  model_result(
    cox_fit(
      "os_event",
      paste(base_rhs, "+ BMI + SMI + albumin + ECOG_hi")
    ),
    "+ BMI, SMI, albumin and ECOG-PS",
    "os_event"
  )
)

# ------------------------------------------------------------
# Diabetes-only models (TAD vs stable diabetes)
# ------------------------------------------------------------
df_dm <- df %>% filter(!is.na(dm_compare))
stopifnot(
  nrow(df_dm) == 81,
  sum(df_dm$dm_compare == "Stable diabetes") == 17,
  sum(df_dm$dm_compare == "TAD") == 64
)

df_dm_no_treatment_only <- df_no_treatment_only %>%
  filter(!is.na(dm_compare))
stopifnot(nrow(df_dm_no_treatment_only) == 70)

df_dm_stable_longitudinal <- df_dm %>%
  filter(
    dm_compare == "TAD" |
      coalesce(stable_longitudinal_data == 1L, FALSE)
  )
stopifnot(nrow(df_dm_stable_longitudinal) == 72)

diabetes_rhs <- "dm_compare + stage + margin + NLR + CA199"

dm_specs <- list(
  "Diabetes only: primary" =
    list(rhs = diabetes_rhs, data = df_dm),
  "Diabetes only: + HbA1c" =
    list(rhs = paste(diabetes_rhs, "+ HbA1c"), data = df_dm),
  "Diabetes only: reduced model" =
    list(rhs = "dm_compare + stage + margin", data = df_dm),
  "Diabetes only: excluding treatment-only worsening" =
    list(rhs = diabetes_rhs, data = df_dm_no_treatment_only),
  "Diabetes only: stable DM with longitudinal data" =
    list(rhs = diabetes_rhs, data = df_dm_stable_longitudinal)
)

dm_results <- bind_rows(lapply(names(dm_specs), function(nm) {
  s <- dm_specs[[nm]]
  model_result(
    cox_fit("os_event", s$rhs, s$data), nm,
    "os_event", "^dm_compareTAD$"
  )
}))

# Group sizes and events, including the reference group
# (stable diabetes; n = 8 in the longitudinal-data model)
dm_counts <- bind_rows(lapply(names(dm_specs), function(nm) {
  group_counts(dm_specs[[nm]]$data, nm)
}))
print(dm_counts)

# Include only the two diabetes-only classification sensitivities
# listed in Supplementary Table S5; the other three belong to S6A.
s5_dm_models <- c(
  "Diabetes only: excluding treatment-only worsening",
  "Diabetes only: stable DM with longitudinal data"
)
s5_dm <- dm_results %>% filter(Model %in% s5_dm_models)
stopifnot(nrow(s5_dm) == 2L)
results_os <- bind_rows(results_os, s5_dm)
stopifnot(nrow(results_os) == 15L, nrow(results_css) == 9L)

# S4 = CSS, S5 = OS
save_csv(
  results_css,
  "Supplementary_Table_S4_CSS_sensitivity.csv"
)
save_csv(
  results_os,
  "Supplementary_Table_S5_OS_sensitivity.csv"
)

# ------------------------------------------------------------
# Supplementary Table S6A: four categories and diabetes only
# ------------------------------------------------------------
# Four-category model (reference: no diabetes); also used for
# the adjusted HRs annotated in Fig. 1c
cox_dm4 <- cox_fit(
  "os_event",
  "dm4 + stage + margin + NLR + CA199"
)
fourcat <- tidy_cox(cox_dm4)
fourcat <- fourcat[grepl("^dm4", fourcat$term), ]

# Direct contrast: worsening vs stable diabetes
b <- coef(cox_dm4)
V <- vcov(cox_dm4)
term_worse <- "dm4Worsening diabetes"
term_stable <- "dm4Stable diabetes"
delta <- b[term_worse] - b[term_stable]
se_delta <- sqrt(
  V[term_worse, term_worse] + V[term_stable, term_stable] -
    2 * V[term_worse, term_stable]
)
fourcat_contrast <- data.frame(
  Contrast = "Worsening vs stable diabetes",
  HR = unname(exp(delta)),
  CI_low = unname(exp(delta - qnorm(0.975) * se_delta)),
  CI_high = unname(exp(delta + qnorm(0.975) * se_delta)),
  P_value = unname(2 * pnorm(-abs(delta / se_delta)))
)

# Group counts and death counts shown in Table S6A
fourcat_counts <- df %>%
  group_by(Diabetes_phenotype = dm4) %>%
  summarise(
    N = n(),
    Deaths = sum(os_event),
    .groups = "drop"
  )
stopifnot(
  identical(fourcat_counts$N, c(81L, 17L, 40L, 24L)),
  identical(fourcat_counts$Deaths, c(56L, 12L, 33L, 19L))
)
save_csv(
  fourcat_counts,
  "Supplementary_Table_S6A_four_category_counts.csv"
)

save_csv(
  fourcat,
  "Supplementary_Table_S6A_four_category_model.csv"
)
save_csv(
  fourcat_contrast,
  "Supplementary_Table_S6A_worsening_vs_stable_contrast.csv"
)
save_csv(
  dm_results,
  "Supplementary_Table_S6A_diabetes_only_models.csv"
)
save_csv(
  dm_counts,
  "Supplementary_Table_S6A_diabetes_only_group_counts.csv"
)

# ------------------------------------------------------------
# Supplementary Fig. S3: stage-specific models
# ------------------------------------------------------------
df_low <- df %>% filter(stage == "<=IIa")
df_high <- df %>% filter(stage == ">=IIb")

stopifnot(
  nrow(df_low) == 43,
  sum(df_low$os_event) == 24,
  sum(df_low$css_event) == 20,
  nrow(df_high) == 119,
  sum(df_high$os_event) == 96,
  sum(df_high$css_event) == 90
)

stage_rhs <- "TAD + margin + NLR + CA199"
stage_results <- bind_rows(
  model_result(
    cox_fit("os_event", stage_rhs, df_low),
    "Stage <=IIa", "os_event"
  ),
  model_result(
    cox_fit("os_event", stage_rhs, df_high),
    "Stage >=IIb", "os_event"
  ),
  model_result(
    cox_fit("css_event", stage_rhs, df_low),
    "Stage <=IIa", "css_event"
  ),
  model_result(
    cox_fit("css_event", stage_rhs, df_high),
    "Stage >=IIb", "css_event"
  )
)
save_csv(
  stage_results,
  "Supplementary_Figure_S3_stage_specific_models.csv"
)

# ------------------------------------------------------------
# Supplementary Table S6B: TAD-by-stage interaction
# (Wald test of the interaction coefficient)
# ------------------------------------------------------------
interaction_rhs <- "TAD * stage + margin + NLR + CA199"
interaction_results <- bind_rows(lapply(
  c("os_event", "css_event"),
  function(ep) {
    fit <- cox_fit(ep, interaction_rhs)
    x <- tidy_cox(fit)
    x <- x[grepl("TAD.*:stage|stage.*:TAD", x$term), ]
    if (nrow(x) != 1) {
      stop("Missing TAD-by-stage interaction coefficient")
    }
    data.frame(
      Endpoint = endpoint_label[[ep]],
      Term = x$term,
      N = fit$n,
      Events = fit$nevent,
      HR_ratio = x$estimate,
      CI_low = x$conf.low,
      CI_high = x$conf.high,
      Wald_P = x$p.value
    )
  }
))
save_csv(
  interaction_results,
  "Supplementary_Table_S6B_TAD_by_stage_interaction_survival.csv"
)

# ------------------------------------------------------------
sink(file.path(out_dir, "sessionInfo_survival.txt"))
print(sessionInfo())
sink()

cat("\nSurvival analyses completed.\n")
print(stage_results)
print(interaction_results)
