# =========================================================
# 06_baseline_nlr_dat_decline.R
#
# Aim 3A: Does baseline NLR predict longitudinal DAT decline?
# PD-only longitudinal mixed-effects models
#
# Outcomes:
#   Primary/Regional:
#     MIA_PUTAMEN_R
#     MIA_PUTAMEN_L
#     MIA_CAUDATE_R
#     MIA_CAUDATE_L
#
# Visits included:
#   BL, V04, V06, V10
#
# Main predictor:
#   baseline_NLR_z
#   Interpretation: per 1-SD higher baseline NLR
#
# Main model:
#   DAT outcome ~ baseline_NLR_z * year_c
#                 + age + sex + bmi + duration_yrs + LEDD + DOMSIDE
#                 + (1 | SITE) + (1 | PATNO)
#
# Sensitivity model:
#   Follow-up DAT only, additionally adjusted for baseline DAT of that region
#
# Input:
#   Revised Excel file with NLR already calculated
#
# Expected input:
#   data/derived/PPMI_with_NLR_all_visits_updated.xlsx
#
# Default output directory:
#   outputs/06_baseline_nlr_dat_decline
#
# Optional command-line usage:
#   Rscript R/06_baseline_nlr_dat_decline.R \
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
  file.path("outputs", "06_baseline_nlr_dat_decline")
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
# Prepare PD-only dataset
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
    
    # Main predictor
    baseline_NLR_z = as.numeric(baseline_NLR_z),
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
    DAT_putamen_R,
    DAT_putamen_L,
    DAT_caudate_R,
    DAT_caudate_L,
    baseline_NLR,
    baseline_NLR_z,
    log_baseline_NLR,
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
  file.path(out_dir, "PD_DAT_analysis_dataset_before_outcome_specific_filter.csv"),
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
# Shared settings
# ---------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))
time_points <- c(0, 1, 2, 4)

# ---------------------------------------------------------
# Main function to run one DAT model
# ---------------------------------------------------------
run_dat_model <- function(
    data,
    outcome_var,
    baseline_outcome_var,
    outcome_label,
    file_stub
) {
  
  cat("\n============================================================\n")
  cat("Running DAT model:", file_stub, "\n")
  cat("Outcome:", outcome_var, "\n")
  cat("============================================================\n")
  
  # -----------------------------------------------------
  # Model variables and missingness before complete-case filtering
  # -----------------------------------------------------
  model_vars_main <- c(
    outcome_var,
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
      n_baseline_NLR_z_missing = sum(is.na(baseline_NLR_z)),
      pct_baseline_NLR_z_missing = 100 * mean(is.na(baseline_NLR_z)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(year_c)
  
  participant_missingness <- data %>%
    dplyr::group_by(PATNO) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      any_model_var_missing = any(
        is.na(.data[[outcome_var]]) |
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
  
  participant_missingness_summary <- participant_missingness %>%
    dplyr::summarise(
      n_participants = dplyr::n(),
      n_with_at_least_one_missing_model_variable = sum(any_model_var_missing),
      pct_with_at_least_one_missing_model_variable = 100 * mean(any_model_var_missing)
    )
  
  write.csv(
    missingness_by_variable,
    file.path(out_dir, paste0(file_stub, "_Missingness_by_model_variable.csv")),
    row.names = FALSE
  )
  
  write.csv(
    missingness_by_visit,
    file.path(out_dir, paste0(file_stub, "_Missingness_by_visit.csv")),
    row.names = FALSE
  )
  
  write.csv(
    participant_missingness,
    file.path(out_dir, paste0(file_stub, "_Participant_level_missingness.csv")),
    row.names = FALSE
  )
  
  write.csv(
    participant_missingness_summary,
    file.path(out_dir, paste0(file_stub, "_Participant_level_missingness_summary.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Complete-case dataset for this outcome
  # -----------------------------------------------------
  dat <- data %>%
    tidyr::drop_na(dplyr::all_of(model_vars_main)) %>%
    droplevels()
  
  if (nrow(dat) < 10) {
    stop(paste0("Not enough complete-case observations for ", file_stub))
  }
  
  write.csv(
    dat,
    file.path(out_dir, paste0(file_stub, "_dataset_complete_case.csv")),
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
    file.path(out_dir, paste0(file_stub, "_descriptives.csv")),
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
    file.path(out_dir, paste0(file_stub, "_visit_counts.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Main model
  # -----------------------------------------------------
  form_main <- as.formula(
    paste0(
      outcome_var,
      " ~ baseline_NLR_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )
  
  lmm <- lmer(
    form_main,
    data = dat,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )

  convergence_messages <- lmm@optinfo$conv$lme4$messages

  convergence_summary <- tibble(
    outcome = outcome_var,
    converged_without_lme4_message = is.null(convergence_messages),
    convergence_message = if (is.null(convergence_messages)) {
      NA_character_
    } else {
      paste(convergence_messages, collapse = " | ")
    },
    singular_fit = lme4::isSingular(lmm, tol = 1e-4)
  )

  write.csv(
    convergence_summary,
    file.path(
      out_dir,
      paste0(file_stub, "_LMM_convergence_and_singularity.csv")
    ),
    row.names = FALSE
  )
  
  fixef_tab <- make_fixef_table(lmm)
  anova_tab <- make_anova_table(lmm)
  randef_var <- as.data.frame(VarCorr(lmm))
  
  model_fit_indices <- tibble(
    outcome = outcome_var,
    model = "main_baseline_NLR_x_time",
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
    fixef_tab,
    file.path(out_dir, paste0(file_stub, "_LMM_fixed_effects_baseline_NLR.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_tab,
    file.path(out_dir, paste0(file_stub, "_LMM_TypeIII_ANOVA_baseline_NLR.csv")),
    row.names = FALSE
  )
  
  write.csv(
    randef_var,
    file.path(out_dir, paste0(file_stub, "_LMM_random_effects_variance_baseline_NLR.csv")),
    row.names = FALSE
  )
  
  write.csv(
    model_fit_indices,
    file.path(out_dir, paste0(file_stub, "_LMM_model_fit_indices_baseline_NLR.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Simple slopes of baseline NLR effect at each time point
  # -----------------------------------------------------
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
    file.path(out_dir, paste0(file_stub, "_Baseline_NLR_effect_on_DAT_by_timepoint.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Continuous plot at baseline_NLR_z = -1.5, 0, +1.5
  # -----------------------------------------------------
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
    file.path(out_dir, paste0(file_stub, "_Predicted_DAT_trajectories_by_baseline_NLR_continuous.csv")),
    row.names = FALSE
  )
  
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
      breaks = c(0, 1, 2, 4),
      labels = c("BL", "1", "2", "4")
    ) +
    labs(
      title = paste0("Baseline NLR and longitudinal decline in ", outcome_label),
      subtitle = "Model-estimated DAT trajectories at low, mean, and high baseline NLR",
      x = "Years from baseline",
      y = outcome_label,
      color = "Baseline NLR",
      fill = "Baseline NLR"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure1_baseline_NLR_continuous.png")),
    plot = p1,
    width = 8.5,
    height = 5.8,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure1_baseline_NLR_continuous.pdf")),
    plot = p1,
    width = 8.5,
    height = 5.8
  )
  
  # -----------------------------------------------------
  # Tertile-based model-estimated trajectories
  # -----------------------------------------------------
  bl_tertiles <- dat %>%
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
      bl_tertiles %>% dplyr::select(PATNO, NLR_tertile),
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
    file.path(out_dir, paste0(file_stub, "_dataset_with_baseline_NLR_tertiles.csv")),
    row.names = FALSE
  )
  
  tert_summary <- bl_tertiles %>%
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
    file.path(out_dir, paste0(file_stub, "_Baseline_NLR_tertile_summary.csv")),
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
    file.path(out_dir, paste0(file_stub, "_Predicted_DAT_trajectories_by_baseline_NLR_tertile.csv")),
    row.names = FALSE
  )
  
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
      breaks = c(0, 1, 2, 4),
      labels = c("BL", "1", "2", "4")
    ) +
    scale_color_discrete(drop = TRUE, na.translate = FALSE) +
    scale_fill_discrete(drop = TRUE, na.translate = FALSE) +
    labs(
      title = paste0("Longitudinal ", outcome_label, " trajectories by baseline NLR tertile"),
      subtitle = "Model-estimated DAT trajectories in Parkinson's disease",
      x = "Years from baseline",
      y = outcome_label,
      color = "Baseline NLR tertile",
      fill = "Baseline NLR tertile"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure2_baseline_NLR_tertiles_model_estimated.png")),
    plot = p2,
    width = 8.5,
    height = 5.8,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure2_baseline_NLR_tertiles_model_estimated.pdf")),
    plot = p2,
    width = 8.5,
    height = 5.8
  )
  
  # -----------------------------------------------------
  # Raw descriptive tertile plot
  # -----------------------------------------------------
  raw_tert_df <- dat_tert %>%
    dplyr::group_by(NLR_tertile, year_c, EVENT_ID) %>%
    dplyr::summarise(
      mean_DAT = mean(.data[[outcome_var]], na.rm = TRUE),
      se_DAT = sd(.data[[outcome_var]], na.rm = TRUE) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(NLR_tertile))
  
  write.csv(
    raw_tert_df,
    file.path(out_dir, paste0(file_stub, "_Raw_DAT_by_baseline_NLR_tertile.csv")),
    row.names = FALSE
  )
  
  p3 <- ggplot(
    raw_tert_df,
    aes(
      x = year_c,
      y = mean_DAT,
      color = NLR_tertile,
      group = NLR_tertile
    )
  ) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(ymin = mean_DAT - se_DAT, ymax = mean_DAT + se_DAT),
      width = 0.12
    ) +
    scale_x_continuous(
      breaks = c(0, 1, 2, 4),
      labels = c("BL", "1", "2", "4")
    ) +
    scale_color_discrete(drop = TRUE, na.translate = FALSE) +
    labs(
      title = paste0("Observed ", outcome_label, " means by baseline NLR tertile"),
      subtitle = "Raw descriptive trajectories",
      x = "Years from baseline",
      y = outcome_label,
      color = "Baseline NLR tertile"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure3_raw_baseline_NLR_tertiles.png")),
    plot = p3,
    width = 8.5,
    height = 5.8,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure3_raw_baseline_NLR_tertiles.pdf")),
    plot = p3,
    width = 8.5,
    height = 5.8
  )
  
  # -----------------------------------------------------
  # Sensitivity model 1:
  # Baseline outcome-adjusted follow-up-only model
  # -----------------------------------------------------
  model_vars_baseline_adj <- c(
    outcome_var,
    baseline_outcome_var,
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
  
  dat_baseline_adj <- data %>%
    dplyr::filter(year_c > 0) %>%
    tidyr::drop_na(dplyr::all_of(model_vars_baseline_adj)) %>%
    droplevels()
  
  form_baseline_adj <- as.formula(
    paste0(
      outcome_var,
      " ~ baseline_NLR_z * year_c + ",
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
  
  baseline_adj_fit_indices <- tibble(
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
    baseline_adj_fit_indices,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_model_fit_indices.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Sensitivity model 2:
  # log baseline NLR predictor
  # -----------------------------------------------------
  model_vars_log <- c(
    outcome_var,
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
  
  dat_log <- data %>%
    tidyr::drop_na(dplyr::all_of(model_vars_log)) %>%
    droplevels()
  
  form_log <- as.formula(
    paste0(
      outcome_var,
      " ~ log_baseline_NLR * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
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
    file.path(out_dir, paste0(file_stub, "_Sensitivity_log_baseline_NLR_dataset.csv")),
    row.names = FALSE
  )
  
  write.csv(
    fixef_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_log_baseline_NLR_LMM_fixed_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_log_baseline_NLR_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Diagnostics for main model
  # -----------------------------------------------------
  diagnostic_df <- dat %>%
    dplyr::mutate(
      fitted_value = fitted(lmm),
      residual = resid(lmm),
      pearson_residual = residual / sigma(lmm)
    )
  
  write.csv(
    diagnostic_df,
    file.path(out_dir, paste0(file_stub, "_Diagnostic_values.csv")),
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
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_residuals_vs_fitted.png")),
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
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_residuals.png")),
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
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_observed_vs_fitted.png")),
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
      title = paste0(file_stub, ": Q-Q plot of participant random intercepts"),
      x = "Theoretical quantiles",
      y = "Participant random intercepts"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_random_intercepts_PATNO.png")),
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
      title = paste0(file_stub, ": Q-Q plot of site random intercepts"),
      x = "Theoretical quantiles",
      y = "Site random intercepts"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_random_intercepts_SITE.png")),
    plot = p_qq_ranef_site,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  diagnostic_check <- performance::check_model(lmm)
  
  png(
    filename = file.path(out_dir, paste0(file_stub, "_Performance_check_model.png")),
    width = 2400,
    height = 1800,
    res = 300
  )
  print(diagnostic_check)
  dev.off()
  
  # -----------------------------------------------------
  # Summary text
  # -----------------------------------------------------
  sink(file.path(out_dir, paste0(file_stub, "_LMM_DAT_baseline_NLR_summary.txt")))
  
  cat("============================================================\n")
  cat("Aim 3A: Does baseline NLR predict longitudinal DAT decline?\n")
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
    " ~ baseline_NLR_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))
  
  cat("Time variable:\n")
  cat("year_c is coded in years from baseline; baseline = 0.\n\n")
  
  cat("Main predictor:\n")
  cat("baseline_NLR_z; effect estimates reflect 1-SD higher baseline NLR.\n\n")
  
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
  
  cat("NLR levels used in Figure 1:\n")
  print(data.frame(
    label = nlr_labels,
    baseline_NLR_z = nlr_z_levels
  ), row.names = FALSE)
  cat("\n\n")
  
  cat("Baseline NLR tertile summary:\n")
  print(tert_summary, row.names = FALSE)
  cat("\n\n")
  
  cat("Sensitivity model: baseline DAT-adjusted follow-up-only model\n")
  cat("Formula:\n")
  cat(paste0(
    outcome_var,
    " ~ baseline_NLR_z * year_c + ",
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
  
  cat("Sensitivity model: log baseline NLR predictor\n")
  cat("Formula:\n")
  cat(paste0(
    outcome_var,
    " ~ log_baseline_NLR * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
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
  cat(paste0("- ", file_stub, "_Diagnostic_residuals_vs_fitted.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_residuals.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_observed_vs_fitted.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_PATNO.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_SITE.png\n"))
  cat(paste0("- ", file_stub, "_Performance_check_model.png\n"))
  
  sink()
  
  return(
    list(
      outcome = outcome_var,
      model = lmm,
      fixef = fixef_tab,
      anova = anova_tab,
      slopes = slopes_by_time_df,
      descriptives = desc,
      missingness = missingness_by_variable,
      participant_missingness_summary = participant_missingness_summary,
      sensitivity_baseline_adj_fixef = fixef_baseline_adj,
      sensitivity_baseline_adj_anova = anova_baseline_adj,
      sensitivity_log_fixef = fixef_log,
      sensitivity_log_anova = anova_log
    )
  )
}

# ---------------------------------------------------------
# Run models
# ---------------------------------------------------------

res_put_r <- run_dat_model(
  data = dat_all0,
  outcome_var = "DAT_putamen_R",
  baseline_outcome_var = "baseline_DAT_putamen_R",
  outcome_label = "Putaminal DAT binding (right)",
  file_stub = "PUTAMEN_R"
)

res_put_l <- run_dat_model(
  data = dat_all0,
  outcome_var = "DAT_putamen_L",
  baseline_outcome_var = "baseline_DAT_putamen_L",
  outcome_label = "Putaminal DAT binding (left)",
  file_stub = "PUTAMEN_L"
)

res_cau_r <- run_dat_model(
  data = dat_all0,
  outcome_var = "DAT_caudate_R",
  baseline_outcome_var = "baseline_DAT_caudate_R",
  outcome_label = "Caudate DAT binding (right)",
  file_stub = "CAUDATE_R"
)

res_cau_l <- run_dat_model(
  data = dat_all0,
  outcome_var = "DAT_caudate_L",
  baseline_outcome_var = "baseline_DAT_caudate_L",
  outcome_label = "Caudate DAT binding (left)",
  file_stub = "CAUDATE_L"
)

# ---------------------------------------------------------
# Combined fixed-effect and ANOVA tables across outcomes
# ---------------------------------------------------------
combined_fixef_main <- bind_rows(
  res_put_r$fixef %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$fixef %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$fixef %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$fixef %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_anova_main <- bind_rows(
  res_put_r$anova %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$anova %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$anova %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$anova %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_slopes <- bind_rows(
  res_put_r$slopes %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$slopes %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$slopes %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$slopes %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

write.csv(
  combined_fixef_main,
  file.path(out_dir, "COMBINED_main_LMM_fixed_effects_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_anova_main,
  file.path(out_dir, "COMBINED_main_LMM_TypeIII_ANOVA_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_slopes,
  file.path(out_dir, "COMBINED_baseline_NLR_effect_on_DAT_by_timepoint_all_outcomes.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("Baseline NLR -> longitudinal DAT decline analyses completed.\n")
cat("All outputs were saved to:\n")
cat(
  normalizePath(out_dir, winslash = "/", mustWork = FALSE),
  "\n"
)
cat("============================================================\n")