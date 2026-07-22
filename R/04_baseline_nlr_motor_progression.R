# =========================================================
# 04_baseline_nlr_motor_progression.R
#
# Aim 2A: Does baseline NLR predict motor progression?
# PD-only longitudinal mixed-effects model
#
# Outcome:
#   ON-medication MDS-UPDRS Part III / UPDRS-III
#
# Main predictor:
#   baseline_NLR_z
#   Interpretation: per 1-SD higher baseline NLR
#
# Main model:
#   updrs3_score_on ~ baseline_NLR_z * year_c
#                     + age + sex + bmi + duration_yrs + LEDD + DOMSIDE
#                     + (1 | SITE) + (1 | PATNO)
#
# Sensitivity model for baseline outcome adjustment:
#   Follow-up UPDRS-III only, additionally adjusted for baseline UPDRS-III
#
# Input:
#   Revised Excel file with NLR already calculated
#
# Expected input:
#   data/derived/PPMI_with_NLR_all_visits_updated.xlsx
#
# Default output directory:
#   outputs/04_baseline_nlr_motor_progression
#
# Optional command-line usage:
#   Rscript R/04_baseline_nlr_motor_progression.R \
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
  file.path("outputs", "04_baseline_nlr_motor_progression")
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
  
  # lmerTest output usually contains df and Pr(>|t|)
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
      p_formatted = fmt_p(p_value)
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
      p_formatted = fmt_p(p_value)
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
  "NLR",
  "baseline_NLR",
  "baseline_NLR_z"
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
    
    # Main predictor
    baseline_NLR_z = as.numeric(baseline_NLR_z),
    
    # Optional skewness-sensitive predictor
    log_baseline_NLR = ifelse(
      !is.na(baseline_NLR) & baseline_NLR > 0,
      log(baseline_NLR),
      NA_real_
    )
  ) %>%
  dplyr::select(
    PATNO,
    SITE,
    EVENT_ID,
    year,
    year_c,
    updrs3_score_on,
    baseline_NLR,
    baseline_NLR_z,
    log_baseline_NLR,
    NLR,
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
# Baseline UPDRS-III extraction for baseline-adjusted sensitivity model
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
  "baseline_NLR_z",
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
    n_baseline_NLR_z_missing = sum(is.na(baseline_NLR_z)),
    pct_baseline_NLR_z_missing = 100 * mean(is.na(baseline_NLR_z)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year_c)

participant_missingness_main <- dat0 %>%
  dplyr::group_by(PATNO) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    any_model_var_missing = any(
      is.na(updrs3_score_on) |
        is.na(baseline_NLR_z) |
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

cat("Main model complete-case dataset:\n")
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
  file.path(out_dir, "PD_NLR_motor_progression_dataset_before_complete_case_filter.csv"),
  row.names = FALSE
)

write.csv(
  dat,
  file.path(out_dir, "PD_NLR_motor_progression_complete_case_dataset.csv"),
  row.names = FALSE
)

write.csv(
  missingness_by_variable_main,
  file.path(out_dir, "Missingness_by_model_variable_main_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  missingness_by_visit_main,
  file.path(out_dir, "Missingness_by_visit_main_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  participant_missingness_main,
  file.path(out_dir, "Participant_level_missingness_main_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  participant_missingness_summary_main,
  file.path(out_dir, "Participant_level_missingness_summary_main_UPDRS.csv"),
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
    mean_baseline_NLR = mean(baseline_NLR, na.rm = TRUE),
    sd_baseline_NLR = sd(baseline_NLR, na.rm = TRUE),
    mean_baseline_NLR_z = mean(baseline_NLR_z, na.rm = TRUE),
    sd_baseline_NLR_z = sd(baseline_NLR_z, na.rm = TRUE),
    mean_LEDD = mean(LEDD, na.rm = TRUE),
    sd_LEDD = sd(LEDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year_c)

write.csv(
  desc,
  file.path(out_dir, "PD_NLR_motor_progression_descriptives.csv"),
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
  file.path(out_dir, "Visit_counts_UPDRS_main_model.csv"),
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
  file.path(out_dir, "SITE_distribution_UPDRS_main_model.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Main mixed model
# ---------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))

lmm <- lmer(
  updrs3_score_on ~ baseline_NLR_z * year_c +
    age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
    (1 | SITE) + (1 | PATNO),
  data = dat,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

model_convergence_messages <- lmm@optinfo$conv$lme4$messages

convergence_summary <- tibble(
  converged_without_lme4_message = is.null(model_convergence_messages),
  convergence_message = if (is.null(model_convergence_messages)) {
    NA_character_
  } else {
    paste(model_convergence_messages, collapse = " | ")
  },
  singular_fit = lme4::isSingular(lmm, tol = 1e-4)
)

write.csv(
  convergence_summary,
  file.path(out_dir, "LMM_convergence_and_singularity_UPDRS_baseline_NLR.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Main model outputs
# ---------------------------------------------------------
fixef_tab <- make_fixef_table(lmm)

write.csv(
  fixef_tab,
  file.path(out_dir, "LMM_fixed_effects_UPDRS_baseline_NLR.csv"),
  row.names = FALSE
)

anova_tab <- make_anova_table(lmm)

write.csv(
  anova_tab,
  file.path(out_dir, "LMM_TypeIII_ANOVA_UPDRS_baseline_NLR.csv"),
  row.names = FALSE
)

randef_var <- as.data.frame(VarCorr(lmm))

write.csv(
  randef_var,
  file.path(out_dir, "LMM_random_effects_variance_UPDRS_baseline_NLR.csv"),
  row.names = FALSE
)

model_fit_indices <- tibble(
  AIC = AIC(lmm),
  BIC = BIC(lmm),
  logLik = as.numeric(logLik(lmm)),
  deviance = deviance(lmm),
  sigma = sigma(lmm),
  n_obs = nobs(lmm),
  n_subjects = dplyr::n_distinct(dat$PATNO),
  n_sites = dplyr::n_distinct(dat$SITE)
)

write.csv(
  model_fit_indices,
  file.path(out_dir, "LMM_model_fit_indices_UPDRS_baseline_NLR.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Simple slopes of baseline NLR effect at each time point
# ---------------------------------------------------------
time_points <- sort(unique(dat$year_c))

slopes_by_time <- emtrends(
  lmm,
  ~ year_c,
  var = "baseline_NLR_z",
  at = list(year_c = time_points)
)

slopes_by_time_df <- as.data.frame(summary(slopes_by_time, infer = c(TRUE, TRUE))) %>%
  standardize_ci_names() %>%
  dplyr::mutate(
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

write.csv(
  slopes_by_time_df,
  file.path(out_dir, "Baseline_NLR_effect_on_UPDRS_by_timepoint.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Predicted trajectories for low / mean / high baseline NLR
# Using baseline_NLR_z = -1.5, 0, +1.5
# ---------------------------------------------------------
nlr_z_levels <- c(-1.5, 0, 1.5)
nlr_labels <- c(
  "Low baseline NLR (-1.5 SD)",
  "Mean baseline NLR",
  "High baseline NLR (+1.5 SD)"
)

emm_cont <- emmeans(
  lmm,
  ~ year_c | baseline_NLR_z,
  at = list(
    year_c = time_points,
    baseline_NLR_z = nlr_z_levels,
    age = mean(dat$age, na.rm = TRUE),
    bmi = mean(dat$bmi, na.rm = TRUE),
    duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
    LEDD = mean(dat$LEDD, na.rm = TRUE)
  )
)

emm_cont_df <- as.data.frame(summary(emm_cont, infer = c(TRUE, TRUE))) %>%
  standardize_ci_names() %>%
  dplyr::mutate(
    NLR_group = factor(
      baseline_NLR_z,
      levels = nlr_z_levels,
      labels = nlr_labels
    ),
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

write.csv(
  emm_cont_df,
  file.path(out_dir, "Predicted_UPDRS_trajectories_by_baseline_NLR_continuous.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Tertile-based model-estimated trajectories
# ---------------------------------------------------------
bl_nlr_tertiles <- dat %>%
  dplyr::filter(EVENT_ID == "BL") %>%
  dplyr::select(PATNO, baseline_NLR_z, baseline_NLR) %>%
  dplyr::distinct(PATNO, .keep_all = TRUE) %>%
  dplyr::mutate(
    NLR_tertile_number = ntile(baseline_NLR_z, 3),
    NLR_tertile = factor(
      NLR_tertile_number,
      levels = c(1, 2, 3),
      labels = c("Low NLR tertile", "Middle NLR tertile", "High NLR tertile")
    )
  ) %>%
  dplyr::filter(!is.na(NLR_tertile))

dat_tert <- dat %>%
  dplyr::left_join(
    bl_nlr_tertiles %>% dplyr::select(PATNO, NLR_tertile),
    by = "PATNO"
  ) %>%
  dplyr::filter(!is.na(NLR_tertile)) %>%
  dplyr::mutate(
    NLR_tertile = factor(
      NLR_tertile,
      levels = c("Low NLR tertile", "Middle NLR tertile", "High NLR tertile")
    )
  )

write.csv(
  dat_tert,
  file.path(out_dir, "PD_NLR_motor_progression_dataset_with_tertiles.csv"),
  row.names = FALSE
)

tert_summary <- bl_nlr_tertiles %>%
  dplyr::group_by(NLR_tertile) %>%
  dplyr::summarise(
    n_subjects = dplyr::n(),
    baseline_NLR_z_mean = mean(baseline_NLR_z, na.rm = TRUE),
    baseline_NLR_z_sd = sd(baseline_NLR_z, na.rm = TRUE),
    baseline_NLR_mean = mean(baseline_NLR, na.rm = TRUE),
    baseline_NLR_sd = sd(baseline_NLR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(match(
    NLR_tertile,
    c("Low NLR tertile", "Middle NLR tertile", "High NLR tertile")
  ))

write.csv(
  tert_summary,
  file.path(out_dir, "Baseline_NLR_tertile_summary.csv"),
  row.names = FALSE
)

tert_means <- tert_summary$baseline_NLR_z_mean
tert_labels <- as.character(tert_summary$NLR_tertile)

emm_tert <- emmeans(
  lmm,
  ~ year_c | baseline_NLR_z,
  at = list(
    year_c = time_points,
    baseline_NLR_z = tert_means,
    age = mean(dat$age, na.rm = TRUE),
    bmi = mean(dat$bmi, na.rm = TRUE),
    duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
    LEDD = mean(dat$LEDD, na.rm = TRUE)
  )
)

emm_tert_df <- as.data.frame(summary(emm_tert, infer = c(TRUE, TRUE))) %>%
  standardize_ci_names() %>%
  dplyr::mutate(
    NLR_tertile = factor(
      baseline_NLR_z,
      levels = tert_means,
      labels = tert_labels
    ),
    p_formatted = vapply(p.value, fmt_p, character(1))
  ) %>%
  dplyr::filter(!is.na(NLR_tertile)) %>%
  dplyr::mutate(
    NLR_tertile = droplevels(NLR_tertile)
  )

write.csv(
  emm_tert_df,
  file.path(out_dir, "Predicted_UPDRS_trajectories_by_NLR_tertile.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Raw descriptive plot by baseline NLR tertiles
# ---------------------------------------------------------
raw_tert_df <- dat_tert %>%
  dplyr::group_by(NLR_tertile, year_c, EVENT_ID) %>%
  dplyr::summarise(
    mean_UPDRS = mean(updrs3_score_on, na.rm = TRUE),
    se_UPDRS = sd(updrs3_score_on, na.rm = TRUE) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(NLR_tertile))

write.csv(
  raw_tert_df,
  file.path(out_dir, "Raw_UPDRS_by_baseline_NLR_tertile.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Main model diagnostic outputs
# ---------------------------------------------------------
diagnostic_df <- dat %>%
  dplyr::mutate(
    fitted_value = fitted(lmm),
    residual = resid(lmm),
    pearson_residual = residual / sigma(lmm)
  )

write.csv(
  diagnostic_df,
  file.path(out_dir, "Diagnostic_values_UPDRS_baseline_NLR.csv"),
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
  filename = file.path(out_dir, "Diagnostic_residuals_vs_fitted_UPDRS_baseline_NLR.png"),
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
  filename = file.path(out_dir, "Diagnostic_QQ_residuals_UPDRS_baseline_NLR.png"),
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
  filename = file.path(out_dir, "Diagnostic_observed_vs_fitted_UPDRS_baseline_NLR.png"),
  plot = p_obs_fit,
  width = 7,
  height = 5,
  dpi = 600
)

ranef_patno <- ranef(lmm)$PATNO %>%
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
  filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_PATNO_UPDRS_baseline_NLR.png"),
  plot = p_qq_ranef_patno,
  width = 7,
  height = 5,
  dpi = 600
)

ranef_site <- ranef(lmm)$SITE %>%
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
  filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_SITE_UPDRS_baseline_NLR.png"),
  plot = p_qq_ranef_site,
  width = 7,
  height = 5,
  dpi = 600
)

tryCatch(
  {
    diagnostic_check <- performance::check_model(lmm)

    png(
      filename = file.path(
        out_dir,
        "Performance_check_model_UPDRS_baseline_NLR.png"
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
# Figure 1: continuous baseline NLR interaction plot
# ---------------------------------------------------------
p1 <- ggplot(
  emm_cont_df,
  aes(
    x = year_c,
    y = emmean,
    color = NLR_group,
    fill = NLR_group,
    group = NLR_group
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
    title = "Baseline NLR predicts longitudinal motor severity in Parkinson's disease",
    subtitle = "Model-estimated ON-medication UPDRS-III trajectories at low, mean, and high baseline NLR",
    x = "Years from baseline",
    y = "Predicted ON-medication UPDRS-III score",
    color = "Baseline NLR",
    fill = "Baseline NLR"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "Figure1_NLRxTime_UPDRS_continuous.png"),
  plot = p1,
  width = 9.5,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Figure1_NLRxTime_UPDRS_continuous.pdf"),
  plot = p1,
  width = 9.5,
  height = 6.8
)

# ---------------------------------------------------------
# Figure 2: tertile-based model-estimated trajectories
# ---------------------------------------------------------
p2 <- ggplot(
  emm_tert_df,
  aes(
    x = year_c,
    y = emmean,
    color = NLR_tertile,
    fill = NLR_tertile,
    group = NLR_tertile
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
  scale_color_discrete(drop = TRUE, na.translate = FALSE) +
  scale_fill_discrete(drop = TRUE, na.translate = FALSE) +
  labs(
    title = "Longitudinal motor trajectories stratified by baseline NLR tertiles",
    subtitle = "Model-estimated ON-medication UPDRS-III trajectories in Parkinson's disease",
    x = "Years from baseline",
    y = "Predicted ON-medication UPDRS-III score",
    color = "Baseline NLR tertile",
    fill = "Baseline NLR tertile"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "Figure2_NLR_tertiles_UPDRS_model_estimated.png"),
  plot = p2,
  width = 9.5,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Figure2_NLR_tertiles_UPDRS_model_estimated.pdf"),
  plot = p2,
  width = 9.5,
  height = 6.8
)

# ---------------------------------------------------------
# Figure 3: raw descriptive trajectories by tertile
# ---------------------------------------------------------
p3 <- ggplot(
  raw_tert_df,
  aes(
    x = year_c,
    y = mean_UPDRS,
    color = NLR_tertile,
    group = NLR_tertile
  )
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(
      ymin = mean_UPDRS - se_UPDRS,
      ymax = mean_UPDRS + se_UPDRS
    ),
    width = 0.12
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  scale_color_discrete(drop = TRUE, na.translate = FALSE) +
  labs(
    title = "Observed UPDRS-III means by baseline NLR tertile",
    subtitle = "Raw descriptive trajectories",
    x = "Years from baseline",
    y = "Observed mean ON-medication UPDRS-III",
    color = "Baseline NLR tertile"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "Figure3_raw_UPDRS_by_NLR_tertile.png"),
  plot = p3,
  width = 9.5,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Figure3_raw_UPDRS_by_NLR_tertile.pdf"),
  plot = p3,
  width = 9.5,
  height = 6.8
)

# ---------------------------------------------------------
# Sensitivity model 1:
# Baseline outcome-adjusted follow-up-only model
# ---------------------------------------------------------
model_vars_sens_baseline_adj <- c(
  "updrs3_score_on",
  "baseline_updrs3_score_on",
  "baseline_NLR_z",
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
  tidyr::drop_na(dplyr::all_of(model_vars_sens_baseline_adj))

dat_baseline_adj <- droplevels(dat_baseline_adj)

if (nrow(dat_baseline_adj) < 10) {
  stop(
    "Fewer than 10 complete follow-up observations were available for the baseline-adjusted sensitivity model.",
    call. = FALSE
  )
}

lmm_baseline_adj <- lmer(
  updrs3_score_on ~ baseline_NLR_z * year_c +
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
  file.path(out_dir, "Sensitivity_baseline_adjusted_followup_only_dataset.csv"),
  row.names = FALSE
)

write.csv(
  fixef_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_LMM_fixed_effects_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  anova_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_LMM_TypeIII_ANOVA_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  randef_baseline_adj,
  file.path(out_dir, "Sensitivity_baseline_adjusted_random_effects_UPDRS.csv"),
  row.names = FALSE
)

baseline_adj_fit_indices <- tibble(
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
  baseline_adj_fit_indices,
  file.path(out_dir, "Sensitivity_baseline_adjusted_model_fit_indices_UPDRS.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Sensitivity model 2:
# log baseline NLR predictor
# ---------------------------------------------------------
model_vars_log <- c(
  "updrs3_score_on",
  "log_baseline_NLR",
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

dat_log <- dat0 %>%
  tidyr::drop_na(dplyr::all_of(model_vars_log))

dat_log <- droplevels(dat_log)

if (nrow(dat_log) < 10) {
  stop(
    "Fewer than 10 complete observations were available for the log-NLR sensitivity model.",
    call. = FALSE
  )
}

lmm_log <- lmer(
  updrs3_score_on ~ log_baseline_NLR * year_c +
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
  file.path(out_dir, "Sensitivity_log_baseline_NLR_dataset.csv"),
  row.names = FALSE
)

write.csv(
  fixef_log,
  file.path(out_dir, "Sensitivity_log_baseline_NLR_LMM_fixed_effects_UPDRS.csv"),
  row.names = FALSE
)

write.csv(
  anova_log,
  file.path(out_dir, "Sensitivity_log_baseline_NLR_LMM_TypeIII_ANOVA_UPDRS.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Full text summary
# ---------------------------------------------------------
sink(file.path(out_dir, "LMM_UPDRS_baseline_NLR_summary.txt"))

cat("============================================================\n")
cat("Aim 2A: Does baseline NLR predict motor progression?\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(file_path, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Main model formula:\n")
cat("updrs3_score_on ~ baseline_NLR_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Time variable:\n")
cat("year_c is coded in years from baseline; baseline = 0.\n\n")

cat("Outcome:\n")
cat("ON-medication UPDRS-III score.\n\n")

cat("Main predictor:\n")
cat("baseline_NLR_z; effect estimates reflect 1-SD higher baseline NLR.\n\n")

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
print(VarCorr(lmm), comp = c("Variance", "Std.Dev."))
cat("\n\n")

cat("Model fit indices, main model:\n")
print(model_fit_indices, row.names = FALSE)
cat("\n\n")

cat("Fixed effects, main model:\n")
print(fixef_tab, row.names = FALSE)
cat("\n\n")

cat("Type III ANOVA, main model:\n")
print(anova_tab, row.names = FALSE)
cat("\n\n")

cat("Simple slopes of baseline NLR effect at each time point:\n")
print(slopes_by_time_df, row.names = FALSE)
cat("\n\n")

cat("NLR levels used in continuous interaction figure:\n")
print(data.frame(
  label = nlr_labels,
  baseline_NLR_z = nlr_z_levels
), row.names = FALSE)
cat("\n\n")

cat("Baseline NLR tertile summary:\n")
print(tert_summary, row.names = FALSE)
cat("\n\n")

cat("Sensitivity model: baseline outcome-adjusted follow-up-only model\n")
cat("Formula:\n")
cat("updrs3_score_on ~ baseline_NLR_z * year_c + baseline_updrs3_score_on + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
cat("Rows =", nrow(dat_baseline_adj), "\n")
cat("Unique PD subjects =", dplyr::n_distinct(dat_baseline_adj$PATNO), "\n")
cat("Unique sites =", dplyr::n_distinct(dat_baseline_adj$SITE), "\n\n")
cat("Fixed effects:\n")
print(fixef_baseline_adj, row.names = FALSE)
cat("\n\n")
cat("Type III ANOVA:\n")
print(anova_baseline_adj, row.names = FALSE)
cat("\n\n")

cat("Sensitivity model: log baseline NLR predictor\n")
cat("Formula:\n")
cat("updrs3_score_on ~ log_baseline_NLR * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
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
cat("- Diagnostic_residuals_vs_fitted_UPDRS_baseline_NLR.png\n")
cat("- Diagnostic_QQ_residuals_UPDRS_baseline_NLR.png\n")
cat("- Diagnostic_observed_vs_fitted_UPDRS_baseline_NLR.png\n")
cat("- Diagnostic_QQ_random_intercepts_PATNO_UPDRS_baseline_NLR.png\n")
cat("- Diagnostic_QQ_random_intercepts_SITE_UPDRS_baseline_NLR.png\n")
cat("- Performance_check_model_UPDRS_baseline_NLR.png\n")

sink()

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("Baseline NLR -> UPDRS-III motor progression analysis completed.\n")
cat("Outputs saved to:\n")
cat(
  normalizePath(out_dir, winslash = "/", mustWork = FALSE),
  "\n"
)
cat("============================================================\n")