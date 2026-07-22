# =========================================================
# 05_timevarying_nlr_motor_severity.R
#
# Aim 2B / Reviewer 1 Comment 2:
# Does time-varying NLR track motor severity over time?
#
# PD-only longitudinal mixed-effects model
#
# Outcome:
#   ON-medication MDS-UPDRS Part III / UPDRS-III
#
# Main predictors:
#   NLR_between_z:
#       Between-person component.
#       Participant's mean NLR across available visits, standardized.
#
#   NLR_within_z:
#       Within-person component.
#       Visit-wise deviation from participant's own mean NLR, standardized.
#
# Interpretation:
#   NLR_between_z tests whether individuals with generally higher NLR
#   across follow-up have higher ON-medication UPDRS-III.
#
#   NLR_within_z tests whether visits with higher-than-usual NLR
#   are accompanied by higher ON-medication UPDRS-III.
#
# Model:
#   updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c
#                     + age + sex + bmi + duration_yrs + LEDD + DOMSIDE
#                     + (1 | SITE) + (1 | PATNO)
#
# Secondary interaction model:
#   updrs3_score_on ~ NLR_between_z + NLR_within_z * year_c
#                     + age + sex + bmi + duration_yrs + LEDD + DOMSIDE
#                     + (1 | SITE) + (1 | PATNO)
#
# Input:
#   Revised Excel file with NLR already calculated
#
# Expected input:
#   data/derived/PPMI_with_NLR_all_visits_updated.xlsx
#
# Default output directory:
#   outputs/05_timevarying_nlr_motor_severity
#
# Optional command-line usage:
#   Rscript R/05_timevarying_nlr_motor_severity.R \
#     path/to/PPMI_with_NLR_all_visits_updated.xlsx \
#     path/to/output_directory
# =========================================================

rm(list = ls())

# ---------------------------------------------------------
# Packages
# ---------------------------------------------------------
required_packages <- c(
  "dplyr", "readxl", "tidyr", "ggplot2",
  "lme4", "lmerTest", "car", "emmeans",
  "tibble", "performance"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script, for example:\n",
      "install.packages(c(",
      paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(tibble)
library(performance)

# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

file_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(
    "data",
    "derived",
    "PPMI_with_NLR_all_visits_updated.xlsx"
  )
}

out_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("outputs", "05_timevarying_nlr_motor_severity")
}

if (!file.exists(file_path)) {
  stop(
    paste0(
      "Input file not found:\n",
      normalizePath(file_path, winslash = "/", mustWork = FALSE),
      "\n\nRun 01_calculate_nlr_and_prepare_datasets.R first, ",
      "or provide the derived workbook as the first command-line argument."
    ),
    call. = FALSE
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Input file: ", normalizePath(file_path, winslash = "/"))
message(
  "Output directory: ",
  normalizePath(out_dir, winslash = "/", mustWork = FALSE)
)

# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------
find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

standardize_ci_names <- function(x) {
  if (!"lower.CL" %in% names(x) && "asymp.LCL" %in% names(x)) {
    x <- x %>% rename(lower.CL = asymp.LCL)
  }
  if (!"upper.CL" %in% names(x) && "asymp.UCL" %in% names(x)) {
    x <- x %>% rename(upper.CL = asymp.UCL)
  }
  if (!"lower.CL" %in% names(x) && "lowerCL" %in% names(x)) {
    x <- x %>% rename(lower.CL = lowerCL)
  }
  if (!"upper.CL" %in% names(x) && "upperCL" %in% names(x)) {
    x <- x %>% rename(upper.CL = upperCL)
  }
  x
}

make_fixef_table <- function(model) {
  out <- as.data.frame(coef(summary(model))) %>%
    tibble::rownames_to_column("term")
  
  out <- out %>%
    dplyr::rename(
      estimate = Estimate,
      SE = `Std. Error`,
      t_value = `t value`
    )
  
  if ("Pr(>|t|)" %in% names(out)) {
    out <- out %>% dplyr::rename(p_value = `Pr(>|t|)`)
  } else if ("Pr(>|z|)" %in% names(out)) {
    out <- out %>% dplyr::rename(p_value = `Pr(>|z|)`)
  } else {
    out$p_value <- NA_real_
  }
  
  if (!"df" %in% names(out)) {
    out$df <- NA_real_
  }
  
  out %>%
    dplyr::mutate(
      lower_95_CI = estimate - 1.96 * SE,
      upper_95_CI = estimate + 1.96 * SE,
      p_formatted = vapply(p_value, fmt_p, character(1))
    ) %>%
    dplyr::select(
      term,
      estimate,
      SE,
      lower_95_CI,
      upper_95_CI,
      df,
      t_value,
      p_value,
      p_formatted
    )
}

make_anova_table <- function(model) {
  out <- as.data.frame(car::Anova(model, type = 3, test.statistic = "Chisq")) %>%
    tibble::rownames_to_column("term")
  
  if ("Chisq" %in% names(out)) {
    out <- out %>% dplyr::rename(test_statistic = Chisq)
  } else if ("F" %in% names(out)) {
    out <- out %>% dplyr::rename(test_statistic = F)
  }
  
  if ("Pr(>Chisq)" %in% names(out)) {
    out <- out %>% dplyr::rename(p_value = `Pr(>Chisq)`)
  } else if ("Pr(>F)" %in% names(out)) {
    out <- out %>% dplyr::rename(p_value = `Pr(>F)`)
  }
  
  if ("Df" %in% names(out)) {
    out <- out %>% dplyr::rename(df = Df)
  }
  
  out %>%
    dplyr::mutate(
      p_formatted = vapply(p_value, fmt_p, character(1))
    )
}

# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------
sheet_name <- "All_data_with_NLR"
available_sheets <- readxl::excel_sheets(file_path)

if (!sheet_name %in% available_sheets) {
  stop(
    paste0(
      "Worksheet '", sheet_name, "' was not found. Available sheets: ",
      paste(available_sheets, collapse = ", ")
    ),
    call. = FALSE
  )
}

df <- readxl::read_excel(file_path, sheet = sheet_name)

cat("\nData loaded successfully.\n")
cat("Satır sayısı:", nrow(df), "\n")
cat("Sütun sayısı:", ncol(df), "\n\n")

# ---------------------------------------------------------
# Required NLR variables
# ---------------------------------------------------------
required_vars <- c(
  "PATNO",
  "EVENT_ID",
  "PRIMDIAG",
  "year_from_baseline",
  "NLR"
)

missing_vars <- setdiff(required_vars, names(df))

if (length(missing_vars) > 0) {
  stop(
    paste0(
      "The following required columns were not found:\n",
      paste(missing_vars, collapse = ", "),
      "\n\nRun 01_calculate_nlr_and_prepare_datasets.R first and check the derived workbook."
    )
  )
}

# ---------------------------------------------------------
# Detect key columns
# ---------------------------------------------------------
age_col <- find_col(
  df,
  c("age", "AGE", "Age", "age_bl", "AGE_AT_VISIT", "AGE_AT_BASELINE", "age_at_baseline")
)

sex_col <- find_col(
  df,
  c("sex", "SEX", "Sex", "gender", "GENDER")
)

bmi_col <- find_col(
  df,
  c("BMI", "bmi", "Bmi", "body_mass_index", "BodyMassIndex")
)

site_col <- find_col(
  df,
  c("SITE", "site", "Site", "siteid", "SITEID", "site_id")
)

duration_col <- find_col(
  df,
  c("duration_yrs", "disease_duration", "disease_duration_years", "DURATION_YRS")
)

updrs_col <- find_col(
  df,
  c("updrs3_score_on", "UPDRS3_score_on", "updrs3_on", "UPDRS3_ON")
)

ledd_col <- find_col(
  df,
  c(
    "LEDD", "ledd", "LEDDTOT", "LEDD_total", "ledd_total",
    "total_ledd", "LEDD_Total", "LEDD_TOT"
  )
)

domside_col <- find_col(
  df,
  c("DOMSIDE", "domside", "DomSide", "dominant_side")
)

cat("Detected columns:\n")
cat("age_col      =", age_col, "\n")
cat("sex_col      =", sex_col, "\n")
cat("bmi_col      =", bmi_col, "\n")
cat("site_col     =", site_col, "\n")
cat("duration_col =", duration_col, "\n")
cat("updrs_col    =", updrs_col, "\n")
cat("ledd_col     =", ledd_col, "\n")
cat("domside_col  =", domside_col, "\n\n")

if (any(is.na(c(age_col, sex_col, bmi_col, site_col, duration_col, updrs_col, ledd_col, domside_col)))) {
  stop(
    "Could not detect one or more required columns. Check names(df) and update candidate names."
  )
}

# ---------------------------------------------------------
# Prepare PD-only longitudinal dataset
# ---------------------------------------------------------
dat0 <- df %>%
  dplyr::filter(
    PRIMDIAG == 1,
    EVENT_ID %in% c("BL", "V04", "V06", "V08", "V10", "V12", "V13", "V14")
  ) %>%
  dplyr::rename(
    year = year_from_baseline,
    SITE = !!site_col,
    age = !!age_col,
    sex = !!sex_col,
    bmi = !!bmi_col,
    duration_yrs = !!duration_col,
    LEDD_raw = !!ledd_col,
    DOMSIDE = !!domside_col,
    updrs3_score_on = !!updrs_col
  ) %>%
  dplyr::mutate(
    PATNO = as.factor(PATNO),
    SITE = as.factor(SITE),
    sex = as.factor(sex),
    DOMSIDE = as.factor(DOMSIDE),
    EVENT_ID = factor(
      EVENT_ID,
      levels = c("BL", "V04", "V06", "V08", "V10", "V12", "V13", "V14")
    ),
    year_c = as.numeric(year),
    
    # Baseline PD participants are expected to be drug-naive.
    # If baseline LEDD is missing or not coded as zero, force baseline to zero.
    LEDD = ifelse(year_c == 0, 0, LEDD_raw),
    
    # Visit-wise NLR
    NLR = as.numeric(NLR),
    log_NLR = ifelse(!is.na(NLR) & NLR > 0, log(NLR), NA_real_)
  ) %>%
  dplyr::select(
    PATNO,
    SITE,
    EVENT_ID,
    year,
    year_c,
    updrs3_score_on,
    NLR,
    log_NLR,
    age,
    sex,
    bmi,
    duration_yrs,
    LEDD,
    LEDD_raw,
    DOMSIDE,
    dplyr::everything()
  )

dat0 <- droplevels(dat0)

# ---------------------------------------------------------
# Decompose visit-wise NLR into between-person and within-person components
# ---------------------------------------------------------
# NLR_between_raw:
#   Participant-specific mean NLR across available visits.
#
# NLR_within_raw:
#   Visit-wise deviation from participant-specific mean NLR.
#
# Standardization:
#   NLR_between_z is standardized across participants/rows.
#   NLR_within_z is standardized across observations.
#
# This is the key analysis requested by Reviewer 1:
#   Does time-varying NLR track motor severity?
# ---------------------------------------------------------

person_nlr <- dat0 %>%
  dplyr::group_by(PATNO) %>%
  dplyr::summarise(
    NLR_between_raw = mean(NLR, na.rm = TRUE),
    n_NLR_available = sum(!is.na(NLR)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    NLR_between_raw = ifelse(is.nan(NLR_between_raw), NA_real_, NLR_between_raw)
  )

dat0 <- dat0 %>%
  dplyr::left_join(person_nlr, by = "PATNO") %>%
  dplyr::mutate(
    NLR_within_raw = NLR - NLR_between_raw
  )

between_mean <- mean(person_nlr$NLR_between_raw, na.rm = TRUE)
between_sd <- sd(person_nlr$NLR_between_raw, na.rm = TRUE)

within_mean <- mean(dat0$NLR_within_raw, na.rm = TRUE)
within_sd <- sd(dat0$NLR_within_raw, na.rm = TRUE)

if (is.na(between_sd) || between_sd == 0) {
  stop(
    "The between-person NLR component has zero or undefined variance.",
    call. = FALSE
  )
}

if (is.na(within_sd) || within_sd == 0) {
  stop(
    "The within-person NLR component has zero or undefined variance.",
    call. = FALSE
  )
}

dat0 <- dat0 %>%
  dplyr::mutate(
    NLR_between_z = as.numeric(
      (NLR_between_raw - between_mean) / between_sd
    ),
    NLR_within_z = as.numeric(
      (NLR_within_raw - within_mean) / within_sd
    )
  )

decomposition_summary <- tibble(
  component = c("Between-person NLR", "Within-person NLR"),
  reference_mean = c(between_mean, within_mean),
  reference_sd = c(between_sd, within_sd),
  n_available = c(
    sum(!is.na(person_nlr$NLR_between_raw)),
    sum(!is.na(dat0$NLR_within_raw))
  )
)

write.csv(
  decomposition_summary,
  file.path(out_dir, "NLR_between_within_standardization_parameters.csv"),
  row.names = FALSE
)

participant_nlr_availability <- person_nlr %>%
  mutate(
    contributes_within_person_information = n_NLR_available >= 2
  )

write.csv(
  participant_nlr_availability,
  file.path(out_dir, "Participant_NLR_availability_for_decomposition.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Check for duplicate baseline records
# ---------------------------------------------------------
duplicate_baseline_records <- dat0 %>%
  filter(EVENT_ID == "BL") %>%
  count(PATNO, name = "n_baseline_records") %>%
  filter(n_baseline_records > 1)

if (nrow(duplicate_baseline_records) > 0) {
  warning(
    paste0(
      nrow(duplicate_baseline_records),
      " participant(s) had duplicate baseline records. ",
      "The first baseline record per participant was retained."
    ),
    call. = FALSE
  )

  write.csv(
    duplicate_baseline_records,
    file.path(out_dir, "Duplicate_baseline_records.csv"),
    row.names = FALSE
  )
}

# ---------------------------------------------------------
# Baseline UPDRS-III extraction for optional baseline-adjusted follow-up model
# ---------------------------------------------------------
baseline_updrs <- dat0 %>%
  dplyr::filter(EVENT_ID == "BL") %>%
  dplyr::select(
    PATNO,
    baseline_updrs3_score_on = updrs3_score_on
  ) %>%
  dplyr::distinct(PATNO, .keep_all = TRUE)

dat0 <- dat0 %>%
  dplyr::left_join(baseline_updrs, by = "PATNO")

# ---------------------------------------------------------
# Missingness summaries before complete-case filtering
# ---------------------------------------------------------
model_vars_main <- c(
  "updrs3_score_on",
  "NLR_between_z",
  "NLR_within_z",
  "year_c",
  "age",
  "sex",
  "bmi",
  "duration_yrs",
  "LEDD",
  "DOMSIDE",
  "SITE",
  "PATNO"
)

missingness_by_variable_main <- dat0 %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(model_vars_main),
      list(
        n_missing = ~sum(is.na(.x)),
        pct_missing = ~100 * mean(is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "variable_stat",
    values_to = "value"
  ) %>%
  tidyr::separate(
    variable_stat,
    into = c("variable", "stat"),
    sep = "_(?=n_missing|pct_missing)",
    remove = TRUE
  ) %>%
  tidyr::pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  dplyr::arrange(variable)

missingness_by_visit_main <- dat0 %>%
  dplyr::group_by(EVENT_ID, year_c) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    n_UPDRS_missing = sum(is.na(updrs3_score_on)),
    pct_UPDRS_missing = 100 * mean(is.na(updrs3_score_on)),
    n_NLR_missing = sum(is.na(NLR)),
    pct_NLR_missing = 100 * mean(is.na(NLR)),
    n_NLR_within_missing = sum(is.na(NLR_within_z)),
    pct_NLR_within_missing = 100 * mean(is.na(NLR_within_z)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year_c)

participant_missingness_main <- dat0 %>%
  dplyr::group_by(PATNO) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    any_model_var_missing = any(
      is.na(updrs3_score_on) |
        is.na(NLR_between_z) |
        is.na(NLR_within_z) |
        is.na(year_c) |
        is.na(age) |
        is.na(sex) |
        is.na(bmi) |
        is.na(duration_yrs) |
        is.na(LEDD) |
        is.na(DOMSIDE) |
        is.na(SITE)
    ),
    .groups = "drop"
  )

participant_missingness_summary_main <- participant_missingness_main %>%
  dplyr::summarise(
    n_participants = dplyr::n(),
    n_with_at_least_one_missing_model_variable = sum(any_model_var_missing),
    pct_with_at_least_one_missing_model_variable = 100 * mean(any_model_var_missing)
  )

# ---------------------------------------------------------
# Complete-case dataset for main model
# ---------------------------------------------------------
dat <- dat0 %>%
  tidyr::drop_na(dplyr::all_of(model_vars_main))

dat <- droplevels(dat)

cat("Main time-varying NLR model complete-case dataset:\n")
cat("Rows:", nrow(dat), "\n")
cat("Unique PD subjects:", dplyr::n_distinct(dat$PATNO), "\n")
cat("Unique sites:", dplyr::n_distinct(dat$SITE), "\n\n")

if (nrow(dat) < 10) {
  stop("Fewer than 10 complete observations were available for the main model.", call. = FALSE)
}

# ---------------------------------------------------------
# Save datasets and missingness outputs
# ---------------------------------------------------------
write.csv(
  dat0,
  file.path(out_dir, "PD_timevarying_NLR_UPDRS_dataset_before_complete_case_filter.csv"),
  row.names = FALSE
)

write.csv(
  dat,
  file.path(out_dir, "PD_timevarying_NLR_UPDRS_complete_case_dataset.csv"),
  row.names = FALSE
)

write.csv(
  person_nlr,
  file.path(out_dir, "Participant_mean_NLR_between_component.csv"),
  row.names = FALSE
)

write.csv(
  missingness_by_variable_main,
  file.path(out_dir, "Missingness_by_model_variable_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  missingness_by_visit_main,
  file.path(out_dir, "Missingness_by_visit_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  participant_missingness_main,
  file.path(out_dir, "Participant_level_missingness_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  participant_missingness_summary_main,
  file.path(out_dir, "Participant_level_missingness_summary_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Descriptives
# ---------------------------------------------------------
desc <- dat %>%
  dplyr::group_by(EVENT_ID, year_c) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    mean_UPDRS = mean(updrs3_score_on, na.rm = TRUE),
    sd_UPDRS = sd(updrs3_score_on, na.rm = TRUE),
    median_UPDRS = median(updrs3_score_on, na.rm = TRUE),
    IQR_UPDRS = IQR(updrs3_score_on, na.rm = TRUE),
    mean_NLR = mean(NLR, na.rm = TRUE),
    sd_NLR = sd(NLR, na.rm = TRUE),
    mean_NLR_between_z = mean(NLR_between_z, na.rm = TRUE),
    sd_NLR_between_z = sd(NLR_between_z, na.rm = TRUE),
    mean_NLR_within_z = mean(NLR_within_z, na.rm = TRUE),
    sd_NLR_within_z = sd(NLR_within_z, na.rm = TRUE),
    mean_LEDD = mean(LEDD, na.rm = TRUE),
    sd_LEDD = sd(LEDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year_c)

write.csv(
  desc,
  file.path(out_dir, "PD_timevarying_NLR_UPDRS_descriptives.csv"),
  row.names = FALSE
)

visit_counts <- dat %>%
  dplyr::group_by(EVENT_ID, year_c) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year_c)

write.csv(
  visit_counts,
  file.path(out_dir, "Visit_counts_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

site_counts <- dat %>%
  dplyr::group_by(SITE) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_rows))

write.csv(
  site_counts,
  file.path(out_dir, "SITE_distribution_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Main model:
# Between-person and within-person time-varying NLR components
# ---------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))

lmm_main <- lmer(
  updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c +
    age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
    (1 | SITE) + (1 | PATNO),
  data = dat,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

main_convergence_messages <- lmm_main@optinfo$conv$lme4$messages

main_convergence_summary <- tibble(
  model = "main_between_within_NLR",
  converged_without_lme4_message = is.null(main_convergence_messages),
  convergence_message = if (is.null(main_convergence_messages)) {
    NA_character_
  } else {
    paste(main_convergence_messages, collapse = " | ")
  },
  singular_fit = lme4::isSingular(lmm_main, tol = 1e-4)
)

write.csv(
  main_convergence_summary,
  file.path(out_dir, "LMM_convergence_main_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

fixef_main <- make_fixef_table(lmm_main)
anova_main <- make_anova_table(lmm_main)
randef_main <- as.data.frame(VarCorr(lmm_main))

write.csv(
  fixef_main,
  file.path(out_dir, "LMM_fixed_effects_timevarying_NLR_UPDRS_main.csv"),
  row.names = FALSE
)

write.csv(
  anova_main,
  file.path(out_dir, "LMM_TypeIII_ANOVA_timevarying_NLR_UPDRS_main.csv"),
  row.names = FALSE
)

write.csv(
  randef_main,
  file.path(out_dir, "LMM_random_effects_variance_timevarying_NLR_UPDRS_main.csv"),
  row.names = FALSE
)

model_fit_main <- tibble(
  model = "main_between_within_NLR",
  AIC = AIC(lmm_main),
  BIC = BIC(lmm_main),
  logLik = as.numeric(logLik(lmm_main)),
  deviance = deviance(lmm_main),
  sigma = sigma(lmm_main),
  n_obs = nobs(lmm_main),
  n_subjects = dplyr::n_distinct(dat$PATNO),
  n_sites = dplyr::n_distinct(dat$SITE)
)

write.csv(
  model_fit_main,
  file.path(out_dir, "LMM_model_fit_indices_timevarying_NLR_UPDRS_main.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Secondary model:
# Does the within-person NLR-UPDRS association vary over time?
# ---------------------------------------------------------
lmm_within_interaction <- lmer(
  updrs3_score_on ~ NLR_between_z + NLR_within_z * year_c +
    age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
    (1 | SITE) + (1 | PATNO),
  data = dat,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

secondary_convergence_messages <-
  lmm_within_interaction@optinfo$conv$lme4$messages

secondary_convergence_summary <- tibble(
  model = "secondary_within_NLR_x_time",
  converged_without_lme4_message = is.null(
    secondary_convergence_messages
  ),
  convergence_message = if (is.null(
    secondary_convergence_messages
  )) {
    NA_character_
  } else {
    paste(secondary_convergence_messages, collapse = " | ")
  },
  singular_fit = lme4::isSingular(
    lmm_within_interaction,
    tol = 1e-4
  )
)

write.csv(
  secondary_convergence_summary,
  file.path(
    out_dir,
    "LMM_convergence_secondary_withinNLRxTime_UPDRS.csv"
  ),
  row.names = FALSE
)

fixef_within_interaction <- make_fixef_table(lmm_within_interaction)
anova_within_interaction <- make_anova_table(lmm_within_interaction)
randef_within_interaction <- as.data.frame(VarCorr(lmm_within_interaction))

write.csv(
  fixef_within_interaction,
  file.path(out_dir, "Secondary_withinNLRxTime_LMM_fixed_effects_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  anova_within_interaction,
  file.path(out_dir, "Secondary_withinNLRxTime_LMM_TypeIII_ANOVA_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  randef_within_interaction,
  file.path(out_dir, "Secondary_withinNLRxTime_random_effects_UPDRS.csv"),
  row.names = FALSE
)

model_fit_within_interaction <- tibble(
  model = "secondary_within_NLR_x_time",
  AIC = AIC(lmm_within_interaction),
  BIC = BIC(lmm_within_interaction),
  logLik = as.numeric(logLik(lmm_within_interaction)),
  deviance = deviance(lmm_within_interaction),
  sigma = sigma(lmm_within_interaction),
  n_obs = nobs(lmm_within_interaction),
  n_subjects = dplyr::n_distinct(dat$PATNO),
  n_sites = dplyr::n_distinct(dat$SITE)
)

write.csv(
  model_fit_within_interaction,
  file.path(out_dir, "Secondary_withinNLRxTime_model_fit_indices_UPDRS.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Simple slopes for within-person NLR effect at each time point
# ---------------------------------------------------------
time_points <- sort(unique(dat$year_c))

within_slopes_by_time <- emtrends(
  lmm_within_interaction,
  ~ year_c,
  var = "NLR_within_z",
  at = list(year_c = time_points)
)

within_slopes_by_time_df <- as.data.frame(summary(within_slopes_by_time, infer = c(TRUE, TRUE))) %>%
  standardize_ci_names() %>%
  dplyr::mutate(
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

write.csv(
  within_slopes_by_time_df,
  file.path(out_dir, "Within_person_NLR_effect_on_UPDRS_by_timepoint.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Sensitivity model:
# Baseline outcome-adjusted follow-up-only model
# ---------------------------------------------------------
model_vars_baseline_adj <- c(
  "updrs3_score_on",
  "baseline_updrs3_score_on",
  "NLR_between_z",
  "NLR_within_z",
  "year_c",
  "age",
  "sex",
  "bmi",
  "duration_yrs",
  "LEDD",
  "DOMSIDE",
  "SITE",
  "PATNO"
)

dat_baseline_adj <- dat0 %>%
  dplyr::filter(year_c > 0) %>%
  tidyr::drop_na(dplyr::all_of(model_vars_baseline_adj))

dat_baseline_adj <- droplevels(dat_baseline_adj)

if (nrow(dat_baseline_adj) < 10) {
  stop(
    "Fewer than 10 complete follow-up observations were available for the baseline-adjusted sensitivity model.",
    call. = FALSE
  )
}

lmm_baseline_adj <- lmer(
  updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c +
    baseline_updrs3_score_on +
    age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
    (1 | SITE) + (1 | PATNO),
  data = dat_baseline_adj,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

fixef_baseline_adj <- make_fixef_table(lmm_baseline_adj)
anova_baseline_adj <- make_anova_table(lmm_baseline_adj)
randef_baseline_adj <- as.data.frame(VarCorr(lmm_baseline_adj))

write.csv(
  dat_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_followup_only_dataset_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  fixef_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_LMM_fixed_effects_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  anova_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_LMM_TypeIII_ANOVA_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  randef_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_random_effects_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

model_fit_baseline_adj <- tibble(
  model = "sensitivity_baseline_adjusted_followup_only",
  AIC = AIC(lmm_baseline_adj),
  BIC = BIC(lmm_baseline_adj),
  logLik = as.numeric(logLik(lmm_baseline_adj)),
  deviance = deviance(lmm_baseline_adj),
  sigma = sigma(lmm_baseline_adj),
  n_obs = nobs(lmm_baseline_adj),
  n_subjects = dplyr::n_distinct(dat_baseline_adj$PATNO),
  n_sites = dplyr::n_distinct(dat_baseline_adj$SITE)
)

write.csv(
  model_fit_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_model_fit_indices_timevarying_NLR_UPDRS.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Sensitivity model:
# Log-NLR between/within decomposition
# ---------------------------------------------------------
person_log_nlr <- dat0 %>%
  dplyr::group_by(PATNO) %>%
  dplyr::summarise(
    log_NLR_between_raw = mean(log_NLR, na.rm = TRUE),
    n_log_NLR_available = sum(!is.na(log_NLR)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    log_NLR_between_raw = ifelse(is.nan(log_NLR_between_raw), NA_real_, log_NLR_between_raw)
  )

dat_log0 <- dat0 %>%
  dplyr::left_join(person_log_nlr, by = "PATNO") %>%
  dplyr::mutate(
    log_NLR_within_raw = log_NLR - log_NLR_between_raw
  )

log_between_mean <- mean(person_log_nlr$log_NLR_between_raw, na.rm = TRUE)
log_between_sd   <- sd(person_log_nlr$log_NLR_between_raw, na.rm = TRUE)

log_within_mean <- mean(dat_log0$log_NLR_within_raw, na.rm = TRUE)
log_within_sd   <- sd(dat_log0$log_NLR_within_raw, na.rm = TRUE)

dat_log0 <- dat_log0 %>%
  dplyr::mutate(
    log_NLR_between_z = as.numeric((log_NLR_between_raw - log_between_mean) / log_between_sd),
    log_NLR_within_z  = as.numeric((log_NLR_within_raw  - log_within_mean) / log_within_sd)
  )

model_vars_log <- c(
  "updrs3_score_on",
  "log_NLR_between_z",
  "log_NLR_within_z",
  "year_c",
  "age",
  "sex",
  "bmi",
  "duration_yrs",
  "LEDD",
  "DOMSIDE",
  "SITE",
  "PATNO"
)

dat_log <- dat_log0 %>%
  tidyr::drop_na(dplyr::all_of(model_vars_log))

dat_log <- droplevels(dat_log)

if (nrow(dat_log) < 10) {
  stop(
    "Fewer than 10 complete observations were available for the log-NLR sensitivity model.",
    call. = FALSE
  )
}

lmm_log <- lmer(
  updrs3_score_on ~ log_NLR_between_z + log_NLR_within_z + year_c +
    age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
    (1 | SITE) + (1 | PATNO),
  data = dat_log,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

fixef_log <- make_fixef_table(lmm_log)
anova_log <- make_anova_table(lmm_log)

write.csv(
  dat_log,
  file.path(out_dir, "Sensitivity_logNLR_between_within_dataset_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  fixef_log,
  file.path(out_dir, "Sensitivity_logNLR_between_within_LMM_fixed_effects_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  anova_log,
  file.path(out_dir, "Sensitivity_logNLR_between_within_LMM_TypeIII_ANOVA_UPDRS.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Model diagnostics for main model
# ---------------------------------------------------------
diagnostic_df <- dat %>%
  dplyr::mutate(
    fitted_value = fitted(lmm_main),
    residual = resid(lmm_main),
    pearson_residual = residual / sigma(lmm_main)
  )

write.csv(
  diagnostic_df,
  file.path(out_dir, "Diagnostic_values_timevarying_NLR_UPDRS_main.csv"),
  row.names = FALSE
)

p_resid_fitted <- ggplot(
  diagnostic_df,
  aes(x = fitted_value, y = residual)
) +
  geom_point(alpha = 0.45, size = 1.6) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  theme_classic(base_size = 13) +
  labs(
    title = "Model diagnostics: residuals vs fitted values",
    x = "Fitted values",
    y = "Residuals"
  )

ggsave(
  filename = file.path(out_dir, "Diagnostic_residuals_vs_fitted_timevarying_NLR_UPDRS.png"),
  plot = p_resid_fitted,
  width = 7,
  height = 5,
  dpi = 600
)

p_qq_resid <- ggplot(
  diagnostic_df,
  aes(sample = residual)
) +
  stat_qq(alpha = 0.55, size = 1.6) +
  stat_qq_line(linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Model diagnostics: Q-Q plot of residuals",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  )

ggsave(
  filename = file.path(out_dir, "Diagnostic_QQ_residuals_timevarying_NLR_UPDRS.png"),
  plot = p_qq_resid,
  width = 7,
  height = 5,
  dpi = 600
)

p_obs_fit <- ggplot(
  diagnostic_df,
  aes(x = fitted_value, y = updrs3_score_on)
) +
  geom_point(alpha = 0.45, size = 1.6) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Model diagnostics: observed vs fitted values",
    x = "Fitted values",
    y = "Observed UPDRS-III"
  )

ggsave(
  filename = file.path(out_dir, "Diagnostic_observed_vs_fitted_timevarying_NLR_UPDRS.png"),
  plot = p_obs_fit,
  width = 7,
  height = 5,
  dpi = 600
)

ranef_patno <- ranef(lmm_main)$PATNO %>%
  tibble::rownames_to_column("PATNO") %>%
  dplyr::rename(random_intercept = `(Intercept)`)

p_qq_ranef_patno <- ggplot(
  ranef_patno,
  aes(sample = random_intercept)
) +
  stat_qq(alpha = 0.55, size = 1.6) +
  stat_qq_line(linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Model diagnostics: Q-Q plot of participant random intercepts",
    x = "Theoretical quantiles",
    y = "Participant random intercepts"
  )

ggsave(
  filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR_UPDRS.png"),
  plot = p_qq_ranef_patno,
  width = 7,
  height = 5,
  dpi = 600
)

ranef_site <- ranef(lmm_main)$SITE %>%
  tibble::rownames_to_column("SITE") %>%
  dplyr::rename(random_intercept = `(Intercept)`)

p_qq_ranef_site <- ggplot(
  ranef_site,
  aes(sample = random_intercept)
) +
  stat_qq(alpha = 0.55, size = 1.6) +
  stat_qq_line(linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Model diagnostics: Q-Q plot of site random intercepts",
    x = "Theoretical quantiles",
    y = "Site random intercepts"
  )

ggsave(
  filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR_UPDRS.png"),
  plot = p_qq_ranef_site,
  width = 7,
  height = 5,
  dpi = 600
)

tryCatch(
  {
    diagnostic_check <- performance::check_model(lmm_main)

    png(
      filename = file.path(
        out_dir,
        "Performance_check_model_timevarying_NLR_UPDRS.png"
      ),
      width = 2400,
      height = 1800,
      res = 300
    )
    print(diagnostic_check)
    dev.off()
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }

    warning(
      paste0(
        "Could not save performance::check_model output: ",
        conditionMessage(e)
      ),
      call. = FALSE
    )
  }
)

# ---------------------------------------------------------
# Visualization 1:
# Relationship between within-person NLR deviation and UPDRS residualized trend
# This is descriptive, not the main inference.
# ---------------------------------------------------------
p_within_scatter <- ggplot(
  dat,
  aes(x = NLR_within_z, y = updrs3_score_on)
) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  theme_classic(base_size = 13) +
  labs(
    title = "Visit-wise within-person NLR deviation and ON-medication UPDRS-III",
    subtitle = "Descriptive plot; inference based on the mixed-effects model",
    x = "Within-person NLR deviation, standardized",
    y = "ON-medication UPDRS-III score"
  )

ggsave(
  filename = file.path(out_dir, "Descriptive_within_person_NLR_vs_UPDRS.png"),
  plot = p_within_scatter,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Descriptive_within_person_NLR_vs_UPDRS.pdf"),
  plot = p_within_scatter,
  width = 7,
  height = 5
)

# ---------------------------------------------------------
# Visualization 2:
# Predicted UPDRS at low / mean / high within-person NLR deviation
# holding between-person NLR at the mean
# ---------------------------------------------------------
within_levels <- c(-1.5, 0, 1.5)
within_labels <- c(
  "Lower-than-usual NLR (-1.5 SD)",
  "Usual NLR",
  "Higher-than-usual NLR (+1.5 SD)"
)

emm_within <- emmeans(
  lmm_main,
  ~ year_c | NLR_within_z,
  at = list(
    year_c = time_points,
    NLR_within_z = within_levels,
    NLR_between_z = 0,
    age = mean(dat$age, na.rm = TRUE),
    bmi = mean(dat$bmi, na.rm = TRUE),
    duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
    LEDD = mean(dat$LEDD, na.rm = TRUE)
  )
)

emm_within_df <- as.data.frame(summary(emm_within, infer = c(TRUE, TRUE))) %>%
  standardize_ci_names() %>%
  dplyr::mutate(
    within_NLR_level = factor(
      NLR_within_z,
      levels = within_levels,
      labels = within_labels
    ),
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

write.csv(
  emm_within_df,
  file.path(out_dir, "Predicted_UPDRS_by_within_person_NLR_levels.csv"),
  row.names = FALSE
)

p_within_pred <- ggplot(
  emm_within_df,
  aes(
    x = year_c,
    y = emmean,
    color = within_NLR_level,
    fill = within_NLR_level,
    group = within_NLR_level
  )
) +
  geom_ribbon(
    aes(ymin = lower.CL, ymax = upper.CL),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  labs(
    title = "Time-varying NLR and ON-medication UPDRS-III",
    subtitle = "Predicted UPDRS-III at lower-than-usual, usual, and higher-than-usual NLR",
    x = "Years from baseline",
    y = "Predicted ON-medication UPDRS-III score",
    color = "Visit-wise NLR",
    fill = "Visit-wise NLR"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "Figure_timevarying_withinNLR_UPDRS_predicted.png"),
  plot = p_within_pred,
  width = 9.5,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Figure_timevarying_withinNLR_UPDRS_predicted.pdf"),
  plot = p_within_pred,
  width = 9.5,
  height = 6.8
)

# ---------------------------------------------------------
# Full text summary
# ---------------------------------------------------------
sink(file.path(out_dir, "LMM_timevarying_NLR_UPDRS_summary.txt"))

cat("============================================================\n")
cat("Aim 2B / Reviewer 1 Comment 2:\n")
cat("Does time-varying NLR track ON-medication UPDRS-III?\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(file_path, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Main model formula:\n")
cat("updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Secondary model formula:\n")
cat("updrs3_score_on ~ NLR_between_z + NLR_within_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Time variable:\n")
cat("year_c is coded in years from baseline; baseline = 0.\n\n")

cat("Outcome:\n")
cat("ON-medication UPDRS-III score.\n\n")

cat("Main predictors:\n")
cat("NLR_between_z = participant mean NLR across available visits, standardized.\n")
cat("NLR_within_z  = visit-wise deviation from participant mean NLR, standardized.\n\n")

cat("Interpretation:\n")
cat("NLR_between_z tests whether individuals with generally higher NLR have higher UPDRS-III across follow-up.\n")
cat("NLR_within_z tests whether visits with higher-than-usual NLR are accompanied by higher UPDRS-III.\n\n")

cat("Detected covariate columns:\n")
cat("age       =", age_col, "\n")
cat("sex       =", sex_col, "\n")
cat("bmi       =", bmi_col, "\n")
cat("site      =", site_col, "\n")
cat("duration  =", duration_col, "\n")
cat("LEDD      =", ledd_col, "\n")
cat("DOMSIDE   =", domside_col, "\n")
cat("UPDRS-III =", updrs_col, "\n\n")

cat("Sample size before complete-case filtering:\n")
cat("Rows =", nrow(dat0), "\n")
cat("Unique PD subjects =", dplyr::n_distinct(dat0$PATNO), "\n")
cat("Unique sites =", dplyr::n_distinct(dat0$SITE), "\n\n")

cat("Sample size in main complete-case model:\n")
cat("Rows =", nrow(dat), "\n")
cat("Unique PD subjects =", dplyr::n_distinct(dat$PATNO), "\n")
cat("Unique sites =", dplyr::n_distinct(dat$SITE), "\n\n")

cat("Visit counts in main model:\n")
print(table(dat$EVENT_ID))
cat("\n\n")

cat("Missingness by model variable, main model:\n")
print(missingness_by_variable_main, row.names = FALSE)
cat("\n\n")

cat("Participant-level missingness summary, main model:\n")
print(participant_missingness_summary_main, row.names = FALSE)
cat("\n\n")

cat("Descriptives by visit:\n")
print(desc, row.names = FALSE)
cat("\n\n")

cat("Random-effects variance, main model:\n")
print(VarCorr(lmm_main), comp = c("Variance", "Std.Dev."))
cat("\n\n")

cat("Model fit indices, main model:\n")
print(model_fit_main, row.names = FALSE)
cat("\n\n")

cat("Fixed effects, main model:\n")
print(fixef_main, row.names = FALSE)
cat("\n\n")

cat("Type III ANOVA, main model:\n")
print(anova_main, row.names = FALSE)
cat("\n\n")

cat("Secondary model: within-person NLR x time\n")
cat("Random-effects variance:\n")
print(VarCorr(lmm_within_interaction), comp = c("Variance", "Std.Dev."))
cat("\n\n")

cat("Model fit indices:\n")
print(model_fit_within_interaction, row.names = FALSE)
cat("\n\n")

cat("Fixed effects:\n")
print(fixef_within_interaction, row.names = FALSE)
cat("\n\n")

cat("Type III ANOVA:\n")
print(anova_within_interaction, row.names = FALSE)
cat("\n\n")

cat("Within-person NLR simple slopes by time point:\n")
print(within_slopes_by_time_df, row.names = FALSE)
cat("\n\n")

cat("Sensitivity model: baseline outcome-adjusted follow-up-only model\n")
cat("Formula:\n")
cat("updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c + baseline_updrs3_score_on + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
cat("Rows =", nrow(dat_baseline_adj), "\n")
cat("Unique PD subjects =", dplyr::n_distinct(dat_baseline_adj$PATNO), "\n")
cat("Unique sites =", dplyr::n_distinct(dat_baseline_adj$SITE), "\n\n")
cat("Fixed effects:\n")
print(fixef_baseline_adj, row.names = FALSE)
cat("\n\n")
cat("Type III ANOVA:\n")
print(anova_baseline_adj, row.names = FALSE)
cat("\n\n")

cat("Sensitivity model: log-NLR between/within decomposition\n")
cat("Formula:\n")
cat("updrs3_score_on ~ log_NLR_between_z + log_NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
cat("Rows =", nrow(dat_log), "\n")
cat("Unique PD subjects =", dplyr::n_distinct(dat_log$PATNO), "\n")
cat("Unique sites =", dplyr::n_distinct(dat_log$SITE), "\n\n")
cat("Fixed effects:\n")
print(fixef_log, row.names = FALSE)
cat("\n\n")
cat("Type III ANOVA:\n")
print(anova_log, row.names = FALSE)
cat("\n\n")

cat("Model assumption checks saved as diagnostic plots:\n")
cat("- Diagnostic_residuals_vs_fitted_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_QQ_residuals_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_observed_vs_fitted_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR_UPDRS.png\n")
cat("- Performance_check_model_timevarying_NLR_UPDRS.png\n\n")

cat("Figures saved:\n")
cat("- Descriptive_within_person_NLR_vs_UPDRS.png/pdf\n")
cat("- Figure_timevarying_withinNLR_UPDRS_predicted.png/pdf\n")

sink()

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("Time-varying NLR -> ON-medication UPDRS-III analysis completed.\n")
cat("Outputs saved to:\n")
cat(
  normalizePath(out_dir, winslash = "/", mustWork = FALSE),
  "\n"
)
cat("============================================================\n")