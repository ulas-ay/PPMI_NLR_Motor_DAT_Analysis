# =========================================================
# 07_timevarying_nlr_dat_binding.R
#
# Aim 3B / Reviewer 1 Comment 2:
# Does time-varying NLR track longitudinal DAT binding?
#
# PD-only longitudinal mixed-effects models
#
# Outcomes:
#   DAT_putamen_R
#   DAT_putamen_L
#   DAT_caudate_R
#   DAT_caudate_L
#
# Visits included:
#   BL, V04, V06, V10
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
#   across follow-up have lower/higher DAT binding.
#
#   NLR_within_z tests whether visits with higher-than-usual NLR
#   are accompanied by lower/higher DAT binding.
#
# Main model:
#   DAT outcome ~ NLR_between_z + NLR_within_z + year_c
#                 + age + sex + bmi + duration_yrs + LEDD + DOMSIDE
#                 + (1 | SITE) + (1 | PATNO)
#
# Secondary model:
#   DAT outcome ~ NLR_between_z + NLR_within_z * year_c
#                 + age + sex + bmi + duration_yrs + LEDD + DOMSIDE
#                 + (1 | SITE) + (1 | PATNO)
#
# Sensitivity model:
#   Follow-up DAT only, additionally adjusted for baseline DAT
#
# Input:
#   Revised Excel file with NLR already calculated
#
# Expected input:
#   data/derived/PPMI_with_NLR_all_visits_updated.xlsx
#
# Default output directory:
#   outputs/07_timevarying_nlr_dat_binding
#
# Optional command-line usage:
#   Rscript R/07_timevarying_nlr_dat_binding.R \
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
  file.path("outputs", "07_timevarying_nlr_dat_binding")
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
# Required variables
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
# Detect columns
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

dat_put_r_col <- find_col(df, c("MIA_PUTAMEN_R", "mia_putamen_r"))
dat_put_l_col <- find_col(df, c("MIA_PUTAMEN_L", "mia_putamen_l"))
dat_cau_r_col <- find_col(df, c("MIA_CAUDATE_R", "mia_caudate_r"))
dat_cau_l_col <- find_col(df, c("MIA_CAUDATE_L", "mia_caudate_l"))

cat("Detected columns:\n")
cat("age_col       =", age_col, "\n")
cat("sex_col       =", sex_col, "\n")
cat("bmi_col       =", bmi_col, "\n")
cat("site_col      =", site_col, "\n")
cat("duration_col  =", duration_col, "\n")
cat("LEDD_col      =", ledd_col, "\n")
cat("DOMSIDE_col   =", domside_col, "\n")
cat("dat_put_r_col =", dat_put_r_col, "\n")
cat("dat_put_l_col =", dat_put_l_col, "\n")
cat("dat_cau_r_col =", dat_cau_r_col, "\n")
cat("dat_cau_l_col =", dat_cau_l_col, "\n\n")

if (any(is.na(c(
  age_col, sex_col, bmi_col, site_col, duration_col,
  ledd_col, domside_col,
  dat_put_r_col, dat_put_l_col, dat_cau_r_col, dat_cau_l_col
)))) {
  stop("Could not detect one or more required columns. Check names(df) and update candidate names.")
}

# ---------------------------------------------------------
# Prepare PD-only DAT dataset
# DAT imaging visits: BL, V04, V06, V10
# ---------------------------------------------------------
dat_all0 <- df %>%
  dplyr::filter(
    PRIMDIAG == 1,
    EVENT_ID %in% c("BL", "V04", "V06", "V10")
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
    DAT_putamen_R = !!dat_put_r_col,
    DAT_putamen_L = !!dat_put_l_col,
    DAT_caudate_R = !!dat_cau_r_col,
    DAT_caudate_L = !!dat_cau_l_col
  ) %>%
  dplyr::mutate(
    PATNO = as.factor(PATNO),
    SITE = as.factor(SITE),
    sex = as.factor(sex),
    DOMSIDE = as.factor(DOMSIDE),
    EVENT_ID = factor(
      EVENT_ID,
      levels = c("BL", "V04", "V06", "V10")
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
    DAT_putamen_R,
    DAT_putamen_L,
    DAT_caudate_R,
    DAT_caudate_L,
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

dat_all0 <- droplevels(dat_all0)

write.csv(
  dat_all0,
  file.path(out_dir, "PD_timevarying_NLR_DAT_dataset_before_outcome_specific_filter.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Decompose visit-wise NLR into between-person and within-person components
# ---------------------------------------------------------
# NLR_between_raw:
#   Participant-specific mean NLR across available DAT imaging visits.
#
# NLR_within_raw:
#   Visit-wise deviation from participant-specific mean NLR.
#
# Standardization:
#   NLR_between_z is standardized across participants.
#   NLR_within_z is standardized across observations.
# ---------------------------------------------------------

person_nlr <- dat_all0 %>%
  dplyr::group_by(PATNO) %>%
  dplyr::summarise(
    NLR_between_raw = mean(NLR, na.rm = TRUE),
    n_NLR_available = sum(!is.na(NLR)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    NLR_between_raw = ifelse(is.nan(NLR_between_raw), NA_real_, NLR_between_raw)
  )

dat_all0 <- dat_all0 %>%
  dplyr::left_join(person_nlr, by = "PATNO") %>%
  dplyr::mutate(
    NLR_within_raw = NLR - NLR_between_raw
  )

between_mean <- mean(person_nlr$NLR_between_raw, na.rm = TRUE)
between_sd <- sd(person_nlr$NLR_between_raw, na.rm = TRUE)

within_mean <- mean(dat_all0$NLR_within_raw, na.rm = TRUE)
within_sd <- sd(dat_all0$NLR_within_raw, na.rm = TRUE)

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

dat_all0 <- dat_all0 %>%
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
    sum(!is.na(dat_all0$NLR_within_raw))
  )
)

write.csv(
  decomposition_summary,
  file.path(out_dir, "NLR_between_within_standardization_parameters_DAT.csv"),
  row.names = FALSE
)

participant_nlr_availability <- person_nlr %>%
  mutate(
    contributes_within_person_information = n_NLR_available >= 2
  )

write.csv(
  participant_nlr_availability,
  file.path(out_dir, "Participant_NLR_availability_DAT_visits.csv"),
  row.names = FALSE
)

write.csv(
  person_nlr,
  file.path(out_dir, "Participant_mean_NLR_between_component_DAT_visits.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Check for duplicate baseline records
# ---------------------------------------------------------
duplicate_baseline_records <- dat_all0 %>%
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
# Baseline DAT extraction for baseline-adjusted sensitivity models
# ---------------------------------------------------------
baseline_dat <- dat_all0 %>%
  dplyr::filter(EVENT_ID == "BL") %>%
  dplyr::select(
    PATNO,
    baseline_DAT_putamen_R = DAT_putamen_R,
    baseline_DAT_putamen_L = DAT_putamen_L,
    baseline_DAT_caudate_R = DAT_caudate_R,
    baseline_DAT_caudate_L = DAT_caudate_L
  ) %>%
  dplyr::distinct(PATNO, .keep_all = TRUE)

dat_all0 <- dat_all0 %>%
  dplyr::left_join(baseline_dat, by = "PATNO")

# ---------------------------------------------------------
# Log-NLR between/within decomposition
# ---------------------------------------------------------
person_log_nlr <- dat_all0 %>%
  dplyr::group_by(PATNO) %>%
  dplyr::summarise(
    log_NLR_between_raw = mean(log_NLR, na.rm = TRUE),
    n_log_NLR_available = sum(!is.na(log_NLR)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    log_NLR_between_raw = ifelse(is.nan(log_NLR_between_raw), NA_real_, log_NLR_between_raw)
  )

dat_all0 <- dat_all0 %>%
  dplyr::left_join(person_log_nlr, by = "PATNO") %>%
  dplyr::mutate(
    log_NLR_within_raw = log_NLR - log_NLR_between_raw
  )

log_between_mean <- mean(
  person_log_nlr$log_NLR_between_raw,
  na.rm = TRUE
)
log_between_sd <- sd(
  person_log_nlr$log_NLR_between_raw,
  na.rm = TRUE
)

log_within_mean <- mean(
  dat_all0$log_NLR_within_raw,
  na.rm = TRUE
)
log_within_sd <- sd(
  dat_all0$log_NLR_within_raw,
  na.rm = TRUE
)

if (is.na(log_between_sd) || log_between_sd == 0) {
  stop(
    "The between-person log-NLR component has zero or undefined variance.",
    call. = FALSE
  )
}

if (is.na(log_within_sd) || log_within_sd == 0) {
  stop(
    "The within-person log-NLR component has zero or undefined variance.",
    call. = FALSE
  )
}

dat_all0 <- dat_all0 %>%
  dplyr::mutate(
    log_NLR_between_z = as.numeric(
      (log_NLR_between_raw - log_between_mean) / log_between_sd
    ),
    log_NLR_within_z = as.numeric(
      (log_NLR_within_raw - log_within_mean) / log_within_sd
    )
  )

log_decomposition_summary <- tibble(
  component = c(
    "Between-person log-NLR",
    "Within-person log-NLR"
  ),
  reference_mean = c(log_between_mean, log_within_mean),
  reference_sd = c(log_between_sd, log_within_sd),
  n_available = c(
    sum(!is.na(person_log_nlr$log_NLR_between_raw)),
    sum(!is.na(dat_all0$log_NLR_within_raw))
  )
)

write.csv(
  log_decomposition_summary,
  file.path(
    out_dir,
    "LogNLR_between_within_standardization_parameters_DAT.csv"
  ),
  row.names = FALSE
)

# ---------------------------------------------------------
# Shared settings
# ---------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))
time_points <- c(0, 1, 2, 4)

# ---------------------------------------------------------
# Main function to run one time-varying NLR -> DAT model
# ---------------------------------------------------------
run_timevarying_dat_model <- function(
    data,
    outcome_var,
    baseline_outcome_var,
    outcome_label,
    file_stub
) {
  
  cat("\n============================================================\n")
  cat("Running time-varying NLR DAT model:", file_stub, "\n")
  cat("Outcome:", outcome_var, "\n")
  cat("============================================================\n")
  
  # -----------------------------------------------------
  # Model variables and missingness before complete-case filtering
  # -----------------------------------------------------
  model_vars_main <- c(
    outcome_var,
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
  
  missingness_by_variable <- data %>%
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
  
  missingness_by_visit <- data %>%
    dplyr::group_by(EVENT_ID, year_c) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_subjects = dplyr::n_distinct(PATNO),
      n_outcome_missing = sum(is.na(.data[[outcome_var]])),
      pct_outcome_missing = 100 * mean(is.na(.data[[outcome_var]])),
      n_NLR_missing = sum(is.na(NLR)),
      pct_NLR_missing = 100 * mean(is.na(NLR)),
      n_NLR_within_missing = sum(is.na(NLR_within_z)),
      pct_NLR_within_missing = 100 * mean(is.na(NLR_within_z)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(year_c)
  
  participant_missingness <- data %>%
    dplyr::group_by(PATNO) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      any_model_var_missing = any(
        is.na(.data[[outcome_var]]) |
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
  
  participant_missingness_summary <- participant_missingness %>%
    dplyr::summarise(
      n_participants = dplyr::n(),
      n_with_at_least_one_missing_model_variable = sum(any_model_var_missing),
      pct_with_at_least_one_missing_model_variable = 100 * mean(any_model_var_missing)
    )
  
  write.csv(
    missingness_by_variable,
    file.path(out_dir, paste0(file_stub, "_Missingness_by_model_variable_timevarying_NLR.csv")),
    row.names = FALSE
  )
  
  write.csv(
    missingness_by_visit,
    file.path(out_dir, paste0(file_stub, "_Missingness_by_visit_timevarying_NLR.csv")),
    row.names = FALSE
  )
  
  write.csv(
    participant_missingness,
    file.path(out_dir, paste0(file_stub, "_Participant_level_missingness_timevarying_NLR.csv")),
    row.names = FALSE
  )
  
  write.csv(
    participant_missingness_summary,
    file.path(out_dir, paste0(file_stub, "_Participant_level_missingness_summary_timevarying_NLR.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Complete-case dataset for main model
  # -----------------------------------------------------
  dat <- data %>%
    tidyr::drop_na(dplyr::all_of(model_vars_main)) %>%
    droplevels()
  
  if (nrow(dat) < 10) {
    stop(paste0("Not enough complete-case observations for ", file_stub))
  }
  
  write.csv(
    dat,
    file.path(out_dir, paste0(file_stub, "_timevarying_NLR_dataset_complete_case.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Descriptives
  # -----------------------------------------------------
  desc <- dat %>%
    dplyr::group_by(EVENT_ID, year_c) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_subjects = dplyr::n_distinct(PATNO),
      mean_DAT = mean(.data[[outcome_var]], na.rm = TRUE),
      sd_DAT = sd(.data[[outcome_var]], na.rm = TRUE),
      median_DAT = median(.data[[outcome_var]], na.rm = TRUE),
      IQR_DAT = IQR(.data[[outcome_var]], na.rm = TRUE),
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
    file.path(out_dir, paste0(file_stub, "_timevarying_NLR_descriptives.csv")),
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
    file.path(out_dir, paste0(file_stub, "_timevarying_NLR_visit_counts.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Main model
  # -----------------------------------------------------
  form_main <- as.formula(
    paste0(
      outcome_var,
      " ~ NLR_between_z + NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )
  
  lmm_main <- lmer(
    form_main,
    data = dat,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )

  main_convergence_messages <- lmm_main@optinfo$conv$lme4$messages

  main_convergence_summary <- tibble(
    outcome = outcome_var,
    model = "main_timevarying_between_within_NLR",
    converged_without_lme4_message = is.null(
      main_convergence_messages
    ),
    convergence_message = if (is.null(
      main_convergence_messages
    )) {
      NA_character_
    } else {
      paste(main_convergence_messages, collapse = " | ")
    },
    singular_fit = lme4::isSingular(lmm_main, tol = 1e-4)
  )

  write.csv(
    main_convergence_summary,
    file.path(
      out_dir,
      paste0(file_stub, "_Main_LMM_convergence_and_singularity.csv")
    ),
    row.names = FALSE
  )
  
  fixef_main <- make_fixef_table(lmm_main)
  anova_main <- make_anova_table(lmm_main)
  randef_main <- as.data.frame(VarCorr(lmm_main))
  
  model_fit_main <- tibble(
    outcome = outcome_var,
    model = "main_timevarying_between_within_NLR",
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
    fixef_main,
    file.path(out_dir, paste0(file_stub, "_LMM_fixed_effects_timevarying_NLR_main.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_main,
    file.path(out_dir, paste0(file_stub, "_LMM_TypeIII_ANOVA_timevarying_NLR_main.csv")),
    row.names = FALSE
  )
  
  write.csv(
    randef_main,
    file.path(out_dir, paste0(file_stub, "_LMM_random_effects_variance_timevarying_NLR_main.csv")),
    row.names = FALSE
  )
  
  write.csv(
    model_fit_main,
    file.path(out_dir, paste0(file_stub, "_LMM_model_fit_indices_timevarying_NLR_main.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Secondary model:
  # Does within-person NLR-DAT association vary over time?
  # -----------------------------------------------------
  form_within_interaction <- as.formula(
    paste0(
      outcome_var,
      " ~ NLR_between_z + NLR_within_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )
  
  lmm_within_interaction <- lmer(
    form_within_interaction,
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
    outcome = outcome_var,
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
      paste0(
        file_stub,
        "_Secondary_withinNLRxTime_convergence_and_singularity.csv"
      )
    ),
    row.names = FALSE
  )
  
  fixef_within_interaction <- make_fixef_table(lmm_within_interaction)
  anova_within_interaction <- make_anova_table(lmm_within_interaction)
  randef_within_interaction <- as.data.frame(VarCorr(lmm_within_interaction))
  
  model_fit_within_interaction <- tibble(
    outcome = outcome_var,
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
    fixef_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_LMM_fixed_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )
  
  write.csv(
    randef_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_random_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    model_fit_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_model_fit_indices.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Simple slopes for within-person NLR effect at each time point
  # -----------------------------------------------------
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
    file.path(out_dir, paste0(file_stub, "_Within_person_NLR_effect_on_DAT_by_timepoint.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Sensitivity model:
  # Baseline DAT-adjusted follow-up-only model
  # -----------------------------------------------------
  model_vars_baseline_adj <- c(
    outcome_var,
    baseline_outcome_var,
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
  
  dat_baseline_adj <- data %>%
    dplyr::filter(year_c > 0) %>%
    tidyr::drop_na(dplyr::all_of(model_vars_baseline_adj)) %>%
    droplevels()
  
  form_baseline_adj <- as.formula(
    paste0(
      outcome_var,
      " ~ NLR_between_z + NLR_within_z + year_c + ",
      baseline_outcome_var,
      " + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )
  
  lmm_baseline_adj <- lmer(
    form_baseline_adj,
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
  
  model_fit_baseline_adj <- tibble(
    outcome = outcome_var,
    model = "sensitivity_baseline_DAT_adjusted_followup_only",
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
    dat_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_dataset.csv")),
    row.names = FALSE
  )
  
  write.csv(
    fixef_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_LMM_fixed_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )
  
  write.csv(
    randef_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_random_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    model_fit_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_model_fit_indices.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Sensitivity model:
  # Log-NLR between/within decomposition
  # -----------------------------------------------------
  model_vars_log <- c(
    outcome_var,
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
  
  dat_log <- data %>%
    tidyr::drop_na(dplyr::all_of(model_vars_log)) %>%
    droplevels()
  
  form_log <- as.formula(
    paste0(
      outcome_var,
      " ~ log_NLR_between_z + log_NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )
  
  lmm_log <- lmer(
    form_log,
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
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_dataset.csv")),
    row.names = FALSE
  )
  
  write.csv(
    fixef_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_LMM_fixed_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Diagnostics for main model
  # -----------------------------------------------------
  diagnostic_df <- dat %>%
    dplyr::mutate(
      fitted_value = fitted(lmm_main),
      residual = resid(lmm_main),
      pearson_residual = residual / sigma(lmm_main)
    )
  
  write.csv(
    diagnostic_df,
    file.path(out_dir, paste0(file_stub, "_Diagnostic_values_timevarying_NLR.csv")),
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
      title = paste0(file_stub, ": residuals vs fitted values"),
      x = "Fitted values",
      y = "Residuals"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_residuals_vs_fitted_timevarying_NLR.png")),
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
      title = paste0(file_stub, ": Q-Q plot of residuals"),
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_residuals_timevarying_NLR.png")),
    plot = p_qq_resid,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  p_obs_fit <- ggplot(
    diagnostic_df,
    aes(x = fitted_value, y = .data[[outcome_var]])
  ) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.7) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(file_stub, ": observed vs fitted values"),
      x = "Fitted values",
      y = paste0("Observed ", outcome_label)
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_observed_vs_fitted_timevarying_NLR.png")),
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
      title = paste0(file_stub, ": Q-Q plot of participant random intercepts"),
      x = "Theoretical quantiles",
      y = "Participant random intercepts"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR.png")),
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
      title = paste0(file_stub, ": Q-Q plot of site random intercepts"),
      x = "Theoretical quantiles",
      y = "Site random intercepts"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR.png")),
    plot = p_qq_ranef_site,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  diagnostic_check <- performance::check_model(lmm_main)
  
  png(
    filename = file.path(out_dir, paste0(file_stub, "_Performance_check_model_timevarying_NLR.png")),
    width = 2400,
    height = 1800,
    res = 300
  )
  print(diagnostic_check)
  dev.off()
  
  # -----------------------------------------------------
  # Visualization:
  # predicted DAT at lower-than-usual, usual, higher-than-usual NLR
  # -----------------------------------------------------
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
    file.path(out_dir, paste0(file_stub, "_Predicted_DAT_by_within_person_NLR_levels.csv")),
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
      breaks = c(0, 1, 2, 4),
      labels = c("BL", "1", "2", "4")
    ) +
    labs(
      title = paste0("Time-varying NLR and ", outcome_label),
      subtitle = "Predicted DAT binding at lower-than-usual, usual, and higher-than-usual NLR",
      x = "Years from baseline",
      y = outcome_label,
      color = "Visit-wise NLR",
      fill = "Visit-wise NLR"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure_timevarying_withinNLR_DAT_predicted.png")),
    plot = p_within_pred,
    width = 8.5,
    height = 5.8,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure_timevarying_withinNLR_DAT_predicted.pdf")),
    plot = p_within_pred,
    width = 8.5,
    height = 5.8
  )
  
  # -----------------------------------------------------
  # Summary text
  # -----------------------------------------------------
  sink(file.path(out_dir, paste0(file_stub, "_LMM_timevarying_NLR_DAT_summary.txt")))
  
  cat("============================================================\n")
  cat("Aim 3B / Reviewer 1 Comment 2:\n")
  cat("Does time-varying NLR track longitudinal DAT binding?\n")
  cat("============================================================\n\n")
  
  cat("Outcome:\n")
  cat(outcome_var, "(", outcome_label, ")\n\n")
  
  cat("Input file:\n")
  cat(file_path, "\n\n")
  
  cat("Output directory:\n")
  cat(out_dir, "\n\n")
  
  cat("Main model formula:\n")
  cat(paste0(
    outcome_var,
    " ~ NLR_between_z + NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))
  
  cat("Secondary model formula:\n")
  cat(paste0(
    outcome_var,
    " ~ NLR_between_z + NLR_within_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))
  
  cat("Time variable:\n")
  cat("year_c is coded in years from baseline; baseline = 0.\n\n")
  
  cat("Main predictors:\n")
  cat("NLR_between_z = participant mean NLR across available DAT imaging visits, standardized.\n")
  cat("NLR_within_z  = visit-wise deviation from participant mean NLR, standardized.\n\n")
  
  cat("Interpretation:\n")
  cat("NLR_between_z tests whether individuals with generally higher NLR have different DAT binding across follow-up.\n")
  cat("NLR_within_z tests whether visits with higher-than-usual NLR are accompanied by different DAT binding.\n\n")
  
  cat("Covariates:\n")
  cat("age, sex, BMI, disease duration, LEDD, DOMSIDE.\n\n")
  
  cat("Sample size before outcome-specific complete-case filtering:\n")
  cat("Rows =", nrow(data), "\n")
  cat("Unique PD subjects =", dplyr::n_distinct(data$PATNO), "\n")
  cat("Unique sites =", dplyr::n_distinct(data$SITE), "\n\n")
  
  cat("Sample size in main complete-case model:\n")
  cat("Rows =", nrow(dat), "\n")
  cat("Unique PD subjects =", dplyr::n_distinct(dat$PATNO), "\n")
  cat("Unique sites =", dplyr::n_distinct(dat$SITE), "\n\n")
  
  cat("Visit counts:\n")
  print(table(dat$EVENT_ID))
  cat("\n\n")
  
  cat("Missingness by model variable:\n")
  print(missingness_by_variable, row.names = FALSE)
  cat("\n\n")
  
  cat("Participant-level missingness summary:\n")
  print(participant_missingness_summary, row.names = FALSE)
  cat("\n\n")
  
  cat("Descriptives by visit:\n")
  print(desc, row.names = FALSE)
  cat("\n\n")
  
  cat("Main model: time-varying NLR between/within components\n")
  cat("Random-effects variance:\n")
  print(VarCorr(lmm_main), comp = c("Variance", "Std.Dev."))
  cat("\n\n")
  cat("Model fit indices:\n")
  print(model_fit_main, row.names = FALSE)
  cat("\n\n")
  cat("Fixed effects:\n")
  print(fixef_main, row.names = FALSE)
  cat("\n\n")
  cat("Type III ANOVA:\n")
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
  
  cat("Sensitivity model: baseline DAT-adjusted follow-up-only model\n")
  cat("Formula:\n")
  cat(paste0(
    outcome_var,
    " ~ NLR_between_z + NLR_within_z + year_c + ",
    baseline_outcome_var,
    " + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))
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
  cat(paste0(
    outcome_var,
    " ~ log_NLR_between_z + log_NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))
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
  cat(paste0("- ", file_stub, "_Diagnostic_residuals_vs_fitted_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_residuals_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_observed_vs_fitted_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Performance_check_model_timevarying_NLR.png\n\n"))
  
  cat("Figure saved:\n")
  cat(paste0("- ", file_stub, "_Figure_timevarying_withinNLR_DAT_predicted.png/pdf\n"))
  
  sink()
  
  return(
    list(
      outcome = outcome_var,
      model_main = lmm_main,
      fixef_main = fixef_main,
      anova_main = anova_main,
      model_within_interaction = lmm_within_interaction,
      fixef_within_interaction = fixef_within_interaction,
      anova_within_interaction = anova_within_interaction,
      within_slopes = within_slopes_by_time_df,
      sensitivity_baseline_adj_fixef = fixef_baseline_adj,
      sensitivity_baseline_adj_anova = anova_baseline_adj,
      sensitivity_log_fixef = fixef_log,
      sensitivity_log_anova = anova_log,
      descriptives = desc,
      missingness = missingness_by_variable,
      participant_missingness_summary = participant_missingness_summary
    )
  )
}

# ---------------------------------------------------------
# Run models
# ---------------------------------------------------------

res_put_r <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_putamen_R",
  baseline_outcome_var = "baseline_DAT_putamen_R",
  outcome_label = "Putaminal DAT binding (right)",
  file_stub = "PUTAMEN_R"
)

res_put_l <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_putamen_L",
  baseline_outcome_var = "baseline_DAT_putamen_L",
  outcome_label = "Putaminal DAT binding (left)",
  file_stub = "PUTAMEN_L"
)

res_cau_r <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_caudate_R",
  baseline_outcome_var = "baseline_DAT_caudate_R",
  outcome_label = "Caudate DAT binding (right)",
  file_stub = "CAUDATE_R"
)

res_cau_l <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_caudate_L",
  baseline_outcome_var = "baseline_DAT_caudate_L",
  outcome_label = "Caudate DAT binding (left)",
  file_stub = "CAUDATE_L"
)

# ---------------------------------------------------------
# Combined tables across outcomes
# ---------------------------------------------------------

combined_fixef_main <- bind_rows(
  res_put_r$fixef_main %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$fixef_main %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$fixef_main %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$fixef_main %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_anova_main <- bind_rows(
  res_put_r$anova_main %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$anova_main %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$anova_main %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$anova_main %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_within_interaction_fixef <- bind_rows(
  res_put_r$fixef_within_interaction %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$fixef_within_interaction %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$fixef_within_interaction %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$fixef_within_interaction %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_within_interaction_anova <- bind_rows(
  res_put_r$anova_within_interaction %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$anova_within_interaction %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$anova_within_interaction %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$anova_within_interaction %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_within_slopes <- bind_rows(
  res_put_r$within_slopes %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$within_slopes %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$within_slopes %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$within_slopes %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_sensitivity_baseline_adj <- bind_rows(
  res_put_r$sensitivity_baseline_adj_fixef %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$sensitivity_baseline_adj_fixef %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$sensitivity_baseline_adj_fixef %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$sensitivity_baseline_adj_fixef %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_sensitivity_log <- bind_rows(
  res_put_r$sensitivity_log_fixef %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$sensitivity_log_fixef %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$sensitivity_log_fixef %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$sensitivity_log_fixef %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

write.csv(
  combined_fixef_main,
  file.path(out_dir, "COMBINED_timevarying_NLR_DAT_main_LMM_fixed_effects_all_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_anova_main,
  file.path(out_dir, "COMBINED_timevarying_NLR_DAT_main_LMM_TypeIII_ANOVA_all_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_within_interaction_fixef,
  file.path(out_dir, "COMBINED_withinNLRxTime_LMM_fixed_effects_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_within_interaction_anova,
  file.path(out_dir, "COMBINED_withinNLRxTime_LMM_TypeIII_ANOVA_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_within_slopes,
  file.path(out_dir, "COMBINED_within_person_NLR_effect_on_DAT_by_timepoint_all_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_baseline_adj,
  file.path(out_dir, "COMBINED_sensitivity_baseline_DAT_adjusted_fixed_effects_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_log,
  file.path(out_dir, "COMBINED_sensitivity_logNLR_between_within_fixed_effects_all_DAT_outcomes.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("Time-varying NLR -> longitudinal DAT binding analyses completed.\n")
cat("All outputs were saved to:\n")
cat(
  normalizePath(out_dir, winslash = "/", mustWork = FALSE),
  "\n"
)
cat("============================================================\n")