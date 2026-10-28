
#!/usr/bin/env Rscript
# ============================================================
# 05_tcell_adjusted_models.R
#
# Outputs (results/tcell_models/):
#   Fig. 2c           Adjusted T-cell density ratios
#                     (values in Supplementary Table S7B file)
#   Suppl. Table S6B  TAD-by-stage interactions, T-cell rows
#                     (OS/CSS rows: 01_survival_analysis.R)
#   Suppl. Table S6C  Four diabetes phenotypes
#   Suppl. Table S7B  Primary HC3 models and bootstrap CIs
#   Suppl. Table S7C  HbA1c sensitivity
#   Suppl. Table S7D  Body-status sensitivity
#
# Selected manuscript estimates are checked with warnings.
#
# Statistical model:
# log(density) ~ TAD + age + stage +
#   log2(CA19-9) + continuous NLR + neoadjuvant chemotherapy
#
# HC3 robust Wald z inference
# BH adjustment across CD4, CD8 and FOXP3
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(sandwich)
})

in_file <- "data/clinical_data.xlsx"
out_dir <- "results/tcell_models"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop("Missing clinical data: ", in_file)
}

raw <- read_excel(in_file, sheet = 1)

required <- c(
  "TAD", "dm_group", "age", "stage",
  "CA199_raw", "NLR_raw", "neoadj", "HbA1c",
  "BMI", "SMI", "albumin", "ECOG_PS",
  "CD4", "CD8", "FOXP3"
)

missing <- setdiff(required, names(raw))
if (length(missing)) {
  stop("Missing clinical columns: ", paste(missing, collapse = ", "))
}

# Category labels match 01_survival_analysis.R
dm_labels <- c(
  "No diabetes", "Stable diabetes",
  "New-onset diabetes", "Worsening diabetes"
)

dat <- raw %>%
  transmute(
    TAD = as.numeric(TAD),
    dm_group_num = as.integer(dm_group),
    dm_group = factor(
      dm_group_num,
      levels = 0:3,
      labels = dm_labels
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

stopifnot(
  nrow(dat) == 162,
  all(complete.cases(dat)),
  sum(dat$TAD == 0) == 98,
  sum(dat$TAD == 1) == 64,
  identical(
    as.integer(table(dat$dm_group_num)),
    c(81L, 17L, 40L, 24L)
  ),
  all(dat$TAD == as.numeric(dat$dm_group_num %in% c(2, 3))),
  all(dat$stage %in% c(0, 1)),
  all(dat$neoadj %in% c(0, 1)),
  all(dat$CA199_raw > 0),
  all(dat$CD4 > 0),
  all(dat$CD8 > 0),
  all(dat$FOXP3 > 0)
)

dat$log2_CA199 <- log2(dat$CA199_raw)
dat$ECOG_hi <- as.integer(dat$ECOG_PS >= 1)

stopifnot(
  sum(dat$TAD == 0 & dat$ECOG_hi == 0) == 82,
  sum(dat$TAD == 0 & dat$ECOG_hi == 1) == 16,
  sum(dat$TAD == 1 & dat$ECOG_hi == 0) == 53,
  sum(dat$TAD == 1 & dat$ECOG_hi == 1) == 11
)

markers <- c("CD4", "CD8", "FOXP3")
zcrit <- qnorm(0.975)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
# HC3 robust inference
hc3_table <- function(fit) {
  V <- vcovHC(fit, type = "HC3")
  b <- coef(fit)
  se <- sqrt(diag(V))

  data.frame(
    term = names(b),
    estimate = unname(b),
    robust_se = unname(se),
    z = unname(b / se),
    p = unname(2 * pnorm(-abs(b / se))),
    stringsAsFactors = FALSE
  )
}

model_formula <- function(marker, type = "primary") {
  rhs <- switch(
    type,
    primary =
      "TAD + age + stage + log2_CA199 + NLR_raw + neoadj",
    hba1c =
      "TAD + age + stage + log2_CA199 + NLR_raw + neoadj + HbA1c",
    body =
      paste(
        "TAD + age + stage + log2_CA199 + NLR_raw + neoadj",
        "+ BMI + SMI + albumin + ECOG_hi"
      ),
    fourcat =
      "dm_group + age + stage + log2_CA199 + NLR_raw + neoadj",
    interaction =
      "TAD * stage + age + log2_CA199 + NLR_raw + neoadj",
    stop("Unknown model type")
  )

  as.formula(
    paste0("log(", marker, ") ~ ", rhs)
  )
}

fit_model <- function(marker, type = "primary", data = dat) {
  lm(model_formula(marker, type), data = data)
}

extract_tad <- function(marker, type) {
  fit <- fit_model(marker, type)
  tab <- hc3_table(fit)
  x <- tab[tab$term == "TAD", , drop = FALSE]

  if (nrow(x) != 1) {
    stop("TAD coefficient not found: ", marker)
  }

  data.frame(
    Marker = marker,
    Ratio = exp(x$estimate),
    Lower95 = exp(x$estimate - zcrit * x$robust_se),
    Upper95 = exp(x$estimate + zcrit * x$robust_se),
    P_value = x$p,
    stringsAsFactors = FALSE
  )
}

run_tad_models <- function(type) {
  bind_rows(
    lapply(markers, extract_tad, type = type)
  ) %>%
    mutate(q_value = p.adjust(P_value, method = "BH"))
}

# Compare ratios with manuscript values; warn on mismatch
check_reference <- function(results, reference, tol, label) {
  observed <- setNames(results$Ratio, results$Marker)
  diff <- abs(observed[names(reference)] - reference)
  if (any(!is.finite(diff)) || any(diff >= tol)) {
    warning(
      label, ": ratios differ from manuscript values (",
      paste(
        sprintf("%s %.4f vs %.4f", names(reference),
                observed[names(reference)], reference),
        collapse = "; "
      ),
      ")",
      call. = FALSE
    )
  }
}

# ------------------------------------------------------------
# Primary models (Fig. 2c, Supplementary Table S7B)
# ------------------------------------------------------------
primary <- run_tad_models("primary")

# Values reported to two decimals in Table S7B / Fig. 2c
check_reference(
  primary,
  c(CD4 = 0.74, CD8 = 0.79, FOXP3 = 0.80),
  tol = 0.005,
  label = "Table S7B primary model"
)

# 2000 patient-level percentile bootstrap resamples.
# The seed is set once; markers are resampled in the order
# CD4, CD8, FOXP3, so each marker's CI depends on that order.
BOOT_B <- 2000L
BOOT_SEED <- 20260803L
set.seed(BOOT_SEED)

bootstrap_marker <- function(marker) {
  vals <- rep(NA_real_, BOOT_B)
  f <- model_formula(marker, "primary")

  for (i in seq_len(BOOT_B)) {
    idx <- sample.int(
      nrow(dat), nrow(dat), replace = TRUE
    )

    fit <- try(
      lm(f, data = dat[idx, , drop = FALSE]),
      silent = TRUE
    )

    if (!inherits(fit, "try-error")) {
      b <- coef(fit)
      if ("TAD" %in% names(b) && is.finite(b["TAD"])) {
        vals[i] <- exp(b["TAD"])
      }
    }
  }

  vals <- vals[is.finite(vals)]
  if (length(vals) < 0.95 * BOOT_B) {
    warning("Bootstrap failures exceeded 5% for ", marker)
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
    Bootstrap_Upper95 = ci[2]
  )
}

bootstrap_results <- bind_rows(
  lapply(markers, bootstrap_marker)
)

primary <- left_join(
  primary, bootstrap_results,
  by = "Marker"
)

write.csv(
  primary,
  file.path(
    out_dir,
    "Supplementary_Table_S7B_primary_HC3_and_bootstrap.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S7C: HbA1c adjustment
# ------------------------------------------------------------
hba1c <- run_tad_models("hba1c")
write.csv(
  hba1c,
  file.path(
    out_dir,
    "Supplementary_Table_S7C_HbA1c_adjusted.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S7D: body-status sensitivity
# ------------------------------------------------------------
body_status <- run_tad_models("body")

# Reference values from the original analysis
check_reference(
  body_status,
  c(CD4 = 0.76640799, CD8 = 0.77958052, FOXP3 = 0.80411894),
  tol = 1e-6,
  label = "Table S7D body-status model"
)

write.csv(
  body_status,
  file.path(
    out_dir,
    "Supplementary_Table_S7D_BMI_SMI_albumin_ECOG_adjusted.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S6C: four diabetes phenotypes
# ------------------------------------------------------------
# Descriptive statistics and Kruskal-Wallis tests
fourcat_descriptive <- bind_rows(
  lapply(markers, function(marker) {
    values <- split(dat[[marker]], dat$dm_group)

    out <- data.frame(
      Marker = marker,
      Kruskal_Wallis_P = kruskal.test(
        dat[[marker]] ~ dat$dm_group
      )$p.value
    )

    for (nm in names(values)) {
      prefix <- make.names(nm)
      x <- values[[nm]]
      out[[paste0(prefix, "_n")]] <- length(x)
      out[[paste0(prefix, "_median")]] <- median(x)
      out[[paste0(prefix, "_Q1")]] <- unname(
        quantile(x, 0.25)
      )
      out[[paste0(prefix, "_Q3")]] <- unname(
        quantile(x, 0.75)
      )
    }
    out
  })
)

write.csv(
  fourcat_descriptive,
  file.path(
    out_dir,
    "Supplementary_Table_S6C_descriptive_four_diabetes_phenotypes.csv"
  ),
  row.names = FALSE
)

# Adjusted model: robust contrasts vs no diabetes and
# HC3 Wald omnibus test (3 df)
fit_fourcat <- function(marker) {
  fit <- fit_model(marker, "fourcat")
  V <- vcovHC(fit, type = "HC3")
  b <- coef(fit)

  terms <- grep("^dm_group", names(b), value = TRUE)
  if (length(terms) != 3) {
    stop("Expected 3 phenotype coefficients for ", marker)
  }

  beta <- b[terms]
  se <- sqrt(diag(V))[terms]
  z <- beta / se

  tab <- data.frame(
    Marker = marker,
    Category = sub("^dm_group", "", terms),
    Reference = levels(dat$dm_group)[1],
    Adjusted_ratio = unname(exp(beta)),
    Lower95 = unname(exp(beta - zcrit * se)),
    Upper95 = unname(exp(beta + zcrit * se)),
    P_value = unname(2 * pnorm(-abs(z))),
    stringsAsFactors = FALSE
  )

  Vs <- V[terms, terms, drop = FALSE]
  W <- as.numeric(t(beta) %*% solve(Vs) %*% beta)
  tab$Omnibus_P <- pchisq(
    W, df = length(beta),
    lower.tail = FALSE
  )
  tab
}

fourcat_adjusted <- bind_rows(
  lapply(markers, fit_fourcat)
)

write.csv(
  fourcat_adjusted,
  file.path(
    out_dir,
    "Supplementary_Table_S6C_adjusted_four_diabetes_phenotypes_long.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S6B (T-cell rows): TAD-by-stage
# interactions. Ratio = ratio of TAD effects
# (stage >=IIb vs <=IIa)
# ------------------------------------------------------------
interaction_results <- bind_rows(
  lapply(markers, function(marker) {
    fit <- fit_model(marker, "interaction")
    tab <- hc3_table(fit)

    x <- tab[
      grepl("TAD:stage|stage:TAD", tab$term),
      ,
      drop = FALSE
    ]

    if (nrow(x) != 1) {
      stop("Interaction term not found: ", marker)
    }

    data.frame(
      Marker = marker,
      Interaction_term = x$term,
      Ratio = exp(x$estimate),
      Lower95 = exp(x$estimate - zcrit * x$robust_se),
      Upper95 = exp(x$estimate + zcrit * x$robust_se),
      Wald_z = x$z,
      P_value = x$p,
      Test = "HC3 robust Wald"
    )
  })
)

write.csv(
  interaction_results,
  file.path(
    out_dir,
    "Supplementary_Table_S6B_TAD_by_stage_interaction_Tcell.csv"
  ),
  row.names = FALSE
)

cat("\nPrimary adjusted T-cell ratios:\n")
print(primary)
cat("\nHbA1c sensitivity:\n")
print(hba1c)
cat("\nBody-status sensitivity:\n")
print(body_status)
cat("\nFour-category analyses:\n")
print(fourcat_adjusted)
cat("\nStage-interaction tests:\n")
print(interaction_results)

sink(file.path(out_dir, "sessionInfo_tcell_models.txt"))
print(sessionInfo())
sink()

cat("\nT-cell analyses completed.\n")
