# ============================================================
# Clinical Script 1: Survival analysis (KM, log-rank, Cox)
# OS and CSS, TAD vs non-TAD
# OS time measured from date of surgery.
# CSS: deaths from non-PDAC causes censored at date of death.
#
# Reproduces:
# Fig. 1a/1b (OS KM); Supplementary Fig. S1 (CSS KM)
# Supplementary Table S2 (OS: univariable + multivariable + sensitivity)
# Supplementary Table S3 (CSS: univariable + multivariable + sensitivity)
#
# Models reported in the manuscript (generated below, OS and CSS in parallel):
# (1) Univariable Cox for every covariate
# (2) Primary multivariable: TAD + stage + margin + NLR + CA19-9(>37 vs <=37)
# (3) Sensitivity: primary + adjuvant chemotherapy
# (4) Sensitivity: CA19-9 as log2-continuous instead of dichotomized
# (5) Sensitivity: log2-continuous CA19-9 + adjuvant chemotherapy
# (6) Sensitivity: excluding neoadjuvant-chemotherapy patients (primary covariates)
# Proportional-hazards assumption checked with scaled Schoenfeld residuals.
# ============================================================
suppressPackageStartupMessages({
library(survival)
library(dplyr)
library(readxl)
library(broom)
})
 
# ---- Data import ----
# clinical_data.xlsx contains one row per patient (n = 162).
# Binary covariates are coded 0/1; raw CA19-9 (U/mL) is retained for the
# log2-continuous sensitivity analysis. Adjust column references to match
# the accompanying public data dictionary.
raw <- read_excel("clinical_data.xlsx", sheet = 1)
 
df <- raw %>%
transmute(
os_time = os_time,
os_event = as.integer(os_event), # 1 = death (any cause), 0 = censored
css_event = as.integer(css_event), # 1 = PDAC death, 0 = censored (incl. other-cause death)
TAD = factor(TAD, levels = c(0, 1), labels = c("non-TAD", "TAD")),
stage = factor(stage, levels = c(0, 1), labels = c("<=IIa", ">=IIb")),
margin = factor(margin, levels = c(0, 1), labels = c("R0", "R1")),
NLR = factor(NLR, levels = c(0, 1), labels = c("<2.1", ">=2.1")),
CA199 = factor(CA199, levels = c(0, 1), labels = c("<=37", ">37")),
CA199_raw = as.numeric(CA199_raw), # raw CA19-9 in U/mL (>0 for all cases)
adj_chemo = factor(adj_chemo, levels = c(0, 1), labels = c("No", "Yes")),
neoadj = factor(neoadj, levels = c(0, 1), labels = c("No", "Yes"))
) %>%
filter(!is.na(os_time), !is.na(os_event))
 
# log2-transformed continuous CA19-9 (all raw values > 0; no offset needed)
df$CA199_log2 <- log2(df$CA199_raw)
 
df_adv <- df %>% filter(stage == ">=IIb") # stage >=IIb subgroup (n = 119)
df_noNAC <- df %>% filter(neoadj == "No") # excluding neoadjuvant chemo (n = 144)
 
# ============================================================
# (0) Kaplan-Meier + log-rank
# ============================================================
cat("\n==================== KM / log-rank ====================\n")
cat("\n-- OS, all cases (Fig. 1a) --\n"); print(survdiff(Surv(os_time, os_event) ~ TAD, data = df))
cat("\n-- OS, stage >=IIb (Fig. 1b) --\n"); print(survdiff(Surv(os_time, os_event) ~ TAD, data = df_adv))
cat("\n-- CSS, all cases (Suppl. Fig. S1a) --\n"); print(survdiff(Surv(os_time, css_event) ~ TAD, data = df))
cat("\n-- CSS, stage >=IIb (Suppl. Fig. S1b) --\n"); print(survdiff(Surv(os_time, css_event) ~ TAD, data = df_adv))
 
# ============================================================
# Helper: run univariable Cox for a set of covariates on a given endpoint
# ============================================================
run_univariable <- function(data, time, event, vars) {
do.call(rbind, lapply(vars, function(v) {
fit <- coxph(as.formula(sprintf("Surv(%s, %s) ~ %s", time, event, v)), data = data)
tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = event, model = "univariable")
}))
}
uni_vars <- c("TAD", "stage", "margin", "NLR", "CA199", "adj_chemo")
 
# ============================================================
# (1) Univariable Cox -- OS and CSS
# ============================================================
cat("\n==================== Univariable Cox ====================\n")
cat("\n-- OS univariable (Suppl. Table S2) --\n"); print(run_univariable(df, "os_time", "os_event", uni_vars))
cat("\n-- CSS univariable (Suppl. Table S3) --\n"); print(run_univariable(df, "os_time", "css_event", uni_vars))
 
# ============================================================
# (2) Primary multivariable -- OS and CSS
# Covariates: TAD + stage + margin + NLR + CA19-9(>37 vs <=37)
# (adjuvant chemotherapy excluded; see Methods)
# ============================================================
cat("\n==================== Primary multivariable ====================\n")
cox_os <- coxph(Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199, data = df)
cat("\n-- OS primary (Suppl. Table S2) --\n"); print(summary(cox_os))
cat("\n-- OS primary: PH assumption (scaled Schoenfeld) --\n"); print(cox.zph(cox_os))
 
cox_css <- coxph(Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199, data = df)
cat("\n-- CSS primary (Suppl. Table S3) --\n"); print(summary(cox_css))
cat("\n-- CSS primary: PH assumption (scaled Schoenfeld) --\n"); print(cox.zph(cox_css))
 
# ============================================================
# (3) Sensitivity: primary + adjuvant chemotherapy -- OS and CSS
# ============================================================
cat("\n==================== + adjuvant chemotherapy ====================\n")
cox_os_adj <- coxph(Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199 + adj_chemo, data = df)
cat("\n-- OS + adj chemo (Suppl. Table S2) --\n"); print(summary(cox_os_adj))
cox_css_adj <- coxph(Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199 + adj_chemo, data = df)
cat("\n-- CSS + adj chemo (Suppl. Table S3) --\n"); print(summary(cox_css_adj))
 
# ============================================================
# (4) Sensitivity: CA19-9 as log2-continuous (instead of >37 dichotomy)
# ============================================================
cat("\n==================== log2-continuous CA19-9 ====================\n")
cox_os_log <- coxph(Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199_log2, data = df)
cat("\n-- OS, log2 CA19-9 (Suppl. Table S2) --\n"); print(summary(cox_os_log))
cox_css_log <- coxph(Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199_log2, data = df)
cat("\n-- CSS, log2 CA19-9 (Suppl. Table S3) --\n"); print(summary(cox_css_log))
 
# ============================================================
# (5) Sensitivity: log2-continuous CA19-9 + adjuvant chemotherapy
# ============================================================
cat("\n==================== log2 CA19-9 + adjuvant chemotherapy ====================\n")
cox_os_log_adj <- coxph(Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199_log2 + adj_chemo, data = df)
cat("\n-- OS, log2 CA19-9 + adj chemo (Suppl. Table S2) --\n"); print(summary(cox_os_log_adj))
cox_css_log_adj <- coxph(Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199_log2 + adj_chemo, data = df)
cat("\n-- CSS, log2 CA19-9 + adj chemo (Suppl. Table S3) --\n"); print(summary(cox_css_log_adj))
 
# ============================================================
# (6) Sensitivity: excluding neoadjuvant-chemotherapy patients
# (same covariates as the primary model; n = 144)
# ============================================================
cat("\n==================== excluding neoadjuvant chemotherapy ====================\n")
cox_os_noNAC <- coxph(Surv(os_time, os_event) ~ TAD + stage + margin + NLR + CA199, data = df_noNAC)
cat("\n-- OS, excl. neoadjuvant --\n"); print(summary(cox_os_noNAC))
cox_css_noNAC <- coxph(Surv(os_time, css_event) ~ TAD + stage + margin + NLR + CA199, data = df_noNAC)
cat("\n-- CSS, excl. neoadjuvant --\n"); print(summary(cox_css_noNAC))
