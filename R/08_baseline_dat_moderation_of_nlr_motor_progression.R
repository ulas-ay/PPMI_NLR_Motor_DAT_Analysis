# =========================================================
# 08_baseline_dat_moderation_of_nlr_motor_progression.R
#
# Aim 4: Does baseline striatal DAT moderate the association
# between baseline NLR and longitudinal motor progression?
#
# PD-only longitudinal mixed-effects models
#
# Outcome:
#   ON-medication UPDRS-III longitudinal trajectory
#
# Predictor:
#   baseline_NLR_z
#
# Moderators:
#   baseline striatal DAT binding
#   1) MIA_CAUDATE_L
#   2) MIA_CAUDATE_R
#   3) MIA_PUTAMEN_L
#   4) MIA_PUTAMEN_R
#
# Main moderation term:
#   baseline_NLR_z * DAT_BL_z * year_c
#
# Covariates:
#   age, sex, BMI, duration_yrs, LEDD, DOMSIDE
#
# Random intercepts:
#   SITE, PATNO
#
# Sensitivity:
#   Follow-up-only model additionally adjusted for baseline UPDRS-III
#
# Expected input:
#   data/derived/PPMI_with_NLR_all_visits_updated.xlsx
#
# Default output directory:
#   outputs/08_baseline_dat_moderation
#
# Optional command-line usage:
#   Rscript R/08_baseline_dat_moderation_of_nlr_motor_progression.R \
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
  file.path("outputs", "08_baseline_dat_moderation")
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
# Helpers
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
# Detect shared columns
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

caudate_l_col <- find_col(df, c("MIA_CAUDATE_L", "mia_caudate_l"))
caudate_r_col <- find_col(df, c("MIA_CAUDATE_R", "mia_caudate_r"))
putamen_l_col <- find_col(df, c("MIA_PUTAMEN_L", "mia_putamen_l"))
putamen_r_col <- find_col(df, c("MIA_PUTAMEN_R", "mia_putamen_r"))

cat("Detected shared columns:\n")
cat("age_col       =", age_col, "\n")
cat("sex_col       =", sex_col, "\n")
cat("bmi_col       =", bmi_col, "\n")
cat("site_col      =", site_col, "\n")
cat("duration_col  =", duration_col, "\n")
cat("updrs_col     =", updrs_col, "\n")
cat("ledd_col      =", ledd_col, "\n")
cat("domside_col   =", domside_col, "\n")
cat("caudate_l_col =", caudate_l_col, "\n")
cat("caudate_r_col =", caudate_r_col, "\n")
cat("putamen_l_col =", putamen_l_col, "\n")
cat("putamen_r_col =", putamen_r_col, "\n\n")

if (any(is.na(c(
  age_col, sex_col, bmi_col, site_col, duration_col,
  updrs_col, ledd_col, domside_col,
  caudate_l_col, caudate_r_col,
  putamen_l_col, putamen_r_col
)))) {
  stop("Could not detect one or more required columns. Check names(df).")
}

# ---------------------------------------------------------
# Prepare PD-only longitudinal UPDRS dataset
# ---------------------------------------------------------
dat_long0 <- df %>%
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
    UPDRS = !!updrs_col,
    DAT_caudate_L = !!caudate_l_col,
    DAT_caudate_R = !!caudate_r_col,
    DAT_putamen_L = !!putamen_l_col,
    DAT_putamen_R = !!putamen_r_col
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
    
    # Baseline participants are expected to be drug-naive.
    # To keep consistency with previous models, force baseline LEDD to zero.
    LEDD = ifelse(year_c == 0, 0, LEDD_raw),
    
    # Main immune predictor
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
    UPDRS,
    baseline_NLR,
    baseline_NLR_z,
    log_baseline_NLR,
    DAT_caudate_L,
    DAT_caudate_R,
    DAT_putamen_L,
    DAT_putamen_R,
    age,
    sex,
    bmi,
    duration_yrs,
    LEDD,
    LEDD_raw,
    DOMSIDE,
    dplyr::everything()
  )

dat_long0 <- droplevels(dat_long0)

write.csv(
  dat_long0,
  file.path(out_dir, "Shared_longitudinal_UPDRS_dataset_before_moderator_filter.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Check for duplicate baseline records
# ---------------------------------------------------------
duplicate_baseline_records <- dat_long0 %>%
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
# Baseline UPDRS for sensitivity models
# ---------------------------------------------------------
baseline_updrs <- dat_long0 %>%
  dplyr::filter(EVENT_ID == "BL") %>%
  dplyr::select(
    PATNO,
    baseline_UPDRS = UPDRS
  ) %>%
  dplyr::distinct(PATNO, .keep_all = TRUE)

dat_long0 <- dat_long0 %>%
  dplyr::left_join(baseline_updrs, by = "PATNO")

# ---------------------------------------------------------
# Main moderation function
# ---------------------------------------------------------
run_striatal_moderation <- function(dat_var, dat_label, file_stub = dat_var) {
  
  cat("\n============================================================\n")
  cat("Running moderation model:", file_stub, "\n")
  cat("DAT moderator:", dat_var, "\n")
  cat("============================================================\n")
  
  # -----------------------------------------------------
  # Baseline DAT extraction
  # -----------------------------------------------------
  bl_dat <- dat_long0 %>%
    dplyr::filter(EVENT_ID == "BL") %>%
    dplyr::select(PATNO, DAT_BL = all_of(dat_var)) %>%
    dplyr::distinct(PATNO, .keep_all = TRUE)
  
  write.csv(
    bl_dat,
    file.path(out_dir, paste0(file_stub, "_baseline_DAT.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Merge baseline DAT into longitudinal UPDRS data
  # -----------------------------------------------------
  dat_bl_mean <- mean(bl_dat$DAT_BL, na.rm = TRUE)
  dat_bl_sd <- sd(bl_dat$DAT_BL, na.rm = TRUE)

  if (is.na(dat_bl_sd) || dat_bl_sd == 0) {
    stop(
      paste0(
        "Baseline DAT moderator has zero or undefined variance for ",
        file_stub,
        "."
      ),
      call. = FALSE
    )
  }

  dat0 <- dat_long0 %>%
    dplyr::left_join(bl_dat, by = "PATNO") %>%
    dplyr::mutate(
      DAT_BL_z = as.numeric(
        (DAT_BL - dat_bl_mean) / dat_bl_sd
      )
    )

  dat_standardization <- tibble(
    moderator = dat_var,
    baseline_DAT_mean = dat_bl_mean,
    baseline_DAT_sd = dat_bl_sd,
    n_available = sum(!is.na(bl_dat$DAT_BL))
  )

  write.csv(
    dat_standardization,
    file.path(
      out_dir,
      paste0(file_stub, "_baseline_DAT_standardization_parameters.csv")
    ),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Missingness before complete-case filtering
  # -----------------------------------------------------
  model_vars_main <- c(
    "UPDRS",
    "baseline_NLR_z",
    "DAT_BL_z",
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
  
  missingness_by_variable <- dat0 %>%
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
  
  missingness_by_visit <- dat0 %>%
    dplyr::group_by(EVENT_ID, year_c) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_subjects = dplyr::n_distinct(PATNO),
      n_UPDRS_missing = sum(is.na(UPDRS)),
      pct_UPDRS_missing = 100 * mean(is.na(UPDRS)),
      n_baseline_NLR_z_missing = sum(is.na(baseline_NLR_z)),
      pct_baseline_NLR_z_missing = 100 * mean(is.na(baseline_NLR_z)),
      n_DAT_BL_z_missing = sum(is.na(DAT_BL_z)),
      pct_DAT_BL_z_missing = 100 * mean(is.na(DAT_BL_z)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(year_c)
  
  participant_missingness <- dat0 %>%
    dplyr::group_by(PATNO) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      any_model_var_missing = any(
        is.na(UPDRS) |
          is.na(baseline_NLR_z) |
          is.na(DAT_BL_z) |
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
  # Complete-case dataset
  # -----------------------------------------------------
  dat <- dat0 %>%
    tidyr::drop_na(dplyr::all_of(model_vars_main)) %>%
    droplevels()
  
  if (nrow(dat) < 10) {
    stop(paste0("Not enough complete-case observations for ", file_stub))
  }
  
  write.csv(
    dat,
    file.path(out_dir, paste0(file_stub, "_moderation_dataset_complete_case.csv")),
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
      mean_UPDRS = mean(UPDRS, na.rm = TRUE),
      sd_UPDRS = sd(UPDRS, na.rm = TRUE),
      median_UPDRS = median(UPDRS, na.rm = TRUE),
      IQR_UPDRS = IQR(UPDRS, na.rm = TRUE),
      mean_baseline_NLR = mean(baseline_NLR, na.rm = TRUE),
      sd_baseline_NLR = sd(baseline_NLR, na.rm = TRUE),
      mean_baseline_NLR_z = mean(baseline_NLR_z, na.rm = TRUE),
      sd_baseline_NLR_z = sd(baseline_NLR_z, na.rm = TRUE),
      mean_DAT_BL = mean(DAT_BL, na.rm = TRUE),
      sd_DAT_BL = sd(DAT_BL, na.rm = TRUE),
      mean_DAT_BL_z = mean(DAT_BL_z, na.rm = TRUE),
      sd_DAT_BL_z = sd(DAT_BL_z, na.rm = TRUE),
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
  
  # -----------------------------------------------------
  # Main model
  # -----------------------------------------------------
  options(contrasts = c("contr.sum", "contr.poly"))
  
  lmm <- lmer(
    UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c +
      age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
      (1 | SITE) + (1 | PATNO),
    data = dat,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )

  convergence_messages <- lmm@optinfo$conv$lme4$messages

  convergence_summary <- tibble(
    moderator = dat_var,
    model = "main_NLR_DAT_time_moderation",
    converged_without_lme4_message = is.null(
      convergence_messages
    ),
    convergence_message = if (is.null(
      convergence_messages
    )) {
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
  
  model_fit <- tibble(
    moderator = dat_var,
    model = "main_NLR_DAT_time_moderation",
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
    file.path(out_dir, paste0(file_stub, "_LMM_fixed_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_tab,
    file.path(out_dir, paste0(file_stub, "_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )
  
  write.csv(
    randef_var,
    file.path(out_dir, paste0(file_stub, "_LMM_random_effects_variance.csv")),
    row.names = FALSE
  )
  
  write.csv(
    model_fit,
    file.path(out_dir, paste0(file_stub, "_LMM_model_fit_indices.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Simple slopes:
  # effect of baseline NLR at each time point,
  # stratified by baseline DAT level
  # -----------------------------------------------------
  time_points <- sort(unique(dat$year_c))
  dat_levels <- c(-1.5, 0, 1.5)
  
  slopes_nlr_by_dat_time <- emtrends(
    lmm,
    ~ year_c | DAT_BL_z,
    var = "baseline_NLR_z",
    at = list(
      year_c = time_points,
      DAT_BL_z = dat_levels
    )
  )
  
  slopes_nlr_by_dat_time_df <- as.data.frame(summary(slopes_nlr_by_dat_time, infer = c(TRUE, TRUE))) %>%
    standardize_ci_names() %>%
    dplyr::mutate(
      DAT_group = factor(
        DAT_BL_z,
        levels = c(-1.5, 0, 1.5),
        labels = c(
          "Low baseline DAT (-1.5 SD)",
          "Mean baseline DAT",
          "High baseline DAT (+1.5 SD)"
        )
      ),
      p_formatted = vapply(p.value, fmt_p, character(1))
    )
  
  write.csv(
    slopes_nlr_by_dat_time_df,
    file.path(out_dir, paste0(file_stub, "_simple_slopes_NLR_by_DAT_and_time.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Model-derived trajectories:
  # NLR = -1.5, 0, +1.5 SD
  # DAT = -1.5, 0, +1.5 SD
  # -----------------------------------------------------
  nlr_levels <- c(-1.5, 0, 1.5)
  
  emm_grid <- emmeans(
    lmm,
    ~ year_c | baseline_NLR_z * DAT_BL_z,
    at = list(
      year_c = time_points,
      baseline_NLR_z = nlr_levels,
      DAT_BL_z = dat_levels,
      age = mean(dat$age, na.rm = TRUE),
      bmi = mean(dat$bmi, na.rm = TRUE),
      duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
      LEDD = mean(dat$LEDD, na.rm = TRUE)
    )
  )
  
  emm_df <- as.data.frame(summary(emm_grid, infer = c(TRUE, TRUE))) %>%
    standardize_ci_names() %>%
    dplyr::mutate(
      NLR_group = factor(
        baseline_NLR_z,
        levels = nlr_levels,
        labels = c(
          "Low baseline NLR (-1.5 SD)",
          "Mean baseline NLR",
          "High baseline NLR (+1.5 SD)"
        )
      ),
      DAT_group = factor(
        DAT_BL_z,
        levels = dat_levels,
        labels = c(
          "Low baseline DAT (-1.5 SD)",
          "Mean baseline DAT",
          "High baseline DAT (+1.5 SD)"
        )
      ),
      p_formatted = vapply(p.value, fmt_p, character(1))
    )
  
  write.csv(
    emm_df,
    file.path(out_dir, paste0(file_stub, "_predicted_UPDRS_trajectories_1p5SD.csv")),
    row.names = FALSE
  )
  
  p1 <- ggplot(
    emm_df,
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
      alpha = 0.16,
      color = NA
    ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.6) +
    facet_wrap(~ DAT_group) +
    scale_x_continuous(
      breaks = time_points,
      labels = ifelse(time_points == 0, "BL", as.character(time_points))
    ) +
    scale_color_discrete(drop = TRUE, na.translate = FALSE) +
    scale_fill_discrete(drop = TRUE, na.translate = FALSE) +
    labs(
      title = paste0(dat_label, " moderates the association between baseline NLR and motor progression"),
      subtitle = "Model-estimated UPDRS-III trajectories at ±1.5 SD levels of baseline NLR and DAT",
      x = "Years from baseline",
      y = "Predicted ON-medication UPDRS-III score",
      color = "Baseline NLR",
      fill = "Baseline NLR"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure1_moderation_1p5SD.png")),
    plot = p1,
    width = 11,
    height = 5.8,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure1_moderation_1p5SD.pdf")),
    plot = p1,
    width = 11,
    height = 5.8
  )
  
  # -----------------------------------------------------
  # Raw descriptive tertile plot
  # -----------------------------------------------------
  dat_plot <- dat %>%
    dplyr::filter(EVENT_ID == "BL") %>%
    dplyr::select(PATNO, baseline_NLR_z, DAT_BL) %>%
    dplyr::distinct(PATNO, .keep_all = TRUE) %>%
    dplyr::mutate(
      NLR_tertile = ntile(baseline_NLR_z, 3),
      DAT_tertile = ntile(DAT_BL, 3),
      NLR_tertile = factor(
        NLR_tertile,
        levels = c(1, 2, 3),
        labels = c("Low NLR tertile", "Middle NLR tertile", "High NLR tertile")
      ),
      DAT_tertile = factor(
        DAT_tertile,
        levels = c(1, 2, 3),
        labels = c("Low DAT tertile", "Middle DAT tertile", "High DAT tertile")
      )
    )
  
  raw_df <- dat %>%
    dplyr::left_join(
      dat_plot %>% dplyr::select(PATNO, NLR_tertile, DAT_tertile),
      by = "PATNO"
    ) %>%
    dplyr::group_by(DAT_tertile, NLR_tertile, year_c, EVENT_ID) %>%
    dplyr::summarise(
      mean_UPDRS = mean(UPDRS, na.rm = TRUE),
      se_UPDRS = sd(UPDRS, na.rm = TRUE) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  write.csv(
    raw_df,
    file.path(out_dir, paste0(file_stub, "_raw_UPDRS_by_NLR_and_DAT_tertiles.csv")),
    row.names = FALSE
  )
  
  p2 <- ggplot(
    raw_df,
    aes(
      x = year_c,
      y = mean_UPDRS,
      color = NLR_tertile,
      group = NLR_tertile
    )
  ) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.4) +
    geom_errorbar(
      aes(ymin = mean_UPDRS - se_UPDRS, ymax = mean_UPDRS + se_UPDRS),
      width = 0.12
    ) +
    facet_wrap(~ DAT_tertile) +
    scale_x_continuous(
      breaks = time_points,
      labels = ifelse(time_points == 0, "BL", as.character(time_points))
    ) +
    scale_color_discrete(drop = TRUE, na.translate = FALSE) +
    labs(
      title = paste0("Observed motor trajectories by baseline NLR and ", dat_label, " tertiles"),
      subtitle = "Raw descriptive means",
      x = "Years from baseline",
      y = "Observed mean ON-medication UPDRS-III",
      color = "Baseline NLR tertile"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure2_moderation_raw_tertiles.png")),
    plot = p2,
    width = 11,
    height = 5.8,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure2_moderation_raw_tertiles.pdf")),
    plot = p2,
    width = 11,
    height = 5.8
  )
  
  # -----------------------------------------------------
  # Sensitivity model:
  # follow-up-only model adjusted for baseline UPDRS
  # -----------------------------------------------------
  model_vars_sens <- c(
    "UPDRS",
    "baseline_UPDRS",
    "baseline_NLR_z",
    "DAT_BL_z",
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
  
  dat_sens <- dat0 %>%
    dplyr::filter(year_c > 0) %>%
    tidyr::drop_na(dplyr::all_of(model_vars_sens)) %>%
    droplevels()
  
  lmm_sens <- lmer(
    UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c +
      baseline_UPDRS +
      age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
      (1 | SITE) + (1 | PATNO),
    data = dat_sens,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
  
  fixef_sens <- make_fixef_table(lmm_sens)
  anova_sens <- make_anova_table(lmm_sens)
  randef_sens <- as.data.frame(VarCorr(lmm_sens))
  
  model_fit_sens <- tibble(
    moderator = dat_var,
    model = "sensitivity_baseline_UPDRS_adjusted_followup_only",
    AIC = AIC(lmm_sens),
    BIC = BIC(lmm_sens),
    logLik = as.numeric(logLik(lmm_sens)),
    deviance = deviance(lmm_sens),
    sigma = sigma(lmm_sens),
    n_obs = nobs(lmm_sens),
    n_subjects = dplyr::n_distinct(dat_sens$PATNO),
    n_sites = dplyr::n_distinct(dat_sens$SITE)
  )
  
  write.csv(
    dat_sens,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_dataset.csv")),
    row.names = FALSE
  )
  
  write.csv(
    fixef_sens,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_LMM_fixed_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    anova_sens,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )
  
  write.csv(
    randef_sens,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_random_effects.csv")),
    row.names = FALSE
  )
  
  write.csv(
    model_fit_sens,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_model_fit_indices.csv")),
    row.names = FALSE
  )
  
  # -----------------------------------------------------
  # Diagnostic plots for main model
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
    aes(x = fitted_value, y = UPDRS)
  ) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.7) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(file_stub, ": observed vs fitted values"),
      x = "Fitted values",
      y = "Observed ON-medication UPDRS-III"
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
  sink(file.path(out_dir, paste0(file_stub, "_moderation_summary.txt")))
  
  cat("============================================================\n")
  cat("STRIATAL DAT moderation analysis\n")
  cat("============================================================\n\n")
  
  cat("DAT moderator:\n")
  cat(dat_var, "(", dat_label, ")\n\n")
  
  cat("Input file:\n")
  cat(file_path, "\n\n")
  
  cat("Output directory:\n")
  cat(out_dir, "\n\n")
  
  cat("Main model formula:\n")
  cat("UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
  
  cat("Sensitivity model formula:\n")
  cat("UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c + baseline_UPDRS + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
  
  cat("Outcome:\n")
  cat("ON-medication UPDRS-III longitudinal trajectory.\n\n")
  
  cat("Main predictor:\n")
  cat("baseline_NLR_z; effect estimates reflect 1-SD higher baseline NLR.\n\n")
  
  cat("Moderator:\n")
  cat("DAT_BL_z; baseline DAT binding standardized within the analytic sample.\n\n")
  
  cat("Time variable:\n")
  cat("year_c is coded in years from baseline; baseline = 0.\n\n")
  
  cat("Covariates:\n")
  cat("age, sex, BMI, disease duration, LEDD, DOMSIDE.\n\n")
  
  cat("Sample size before complete-case filtering:\n")
  cat("Rows =", nrow(dat0), "\n")
  cat("Unique PD subjects =", dplyr::n_distinct(dat0$PATNO), "\n")
  cat("Unique sites =", dplyr::n_distinct(dat0$SITE), "\n\n")
  
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
  print(model_fit, row.names = FALSE)
  cat("\n\n")
  
  cat("Fixed effects, main model:\n")
  print(fixef_tab, row.names = FALSE)
  cat("\n\n")
  
  cat("Type III ANOVA, main model:\n")
  print(anova_tab, row.names = FALSE)
  cat("\n\n")
  
  cat("Simple slopes of baseline NLR at each time point, stratified by baseline DAT level:\n")
  print(slopes_nlr_by_dat_time_df, row.names = FALSE)
  cat("\n\n")
  
  cat("Sensitivity model: baseline UPDRS-adjusted follow-up-only model\n")
  cat("Rows =", nrow(dat_sens), "\n")
  cat("Unique PD subjects =", dplyr::n_distinct(dat_sens$PATNO), "\n")
  cat("Unique sites =", dplyr::n_distinct(dat_sens$SITE), "\n\n")
  
  cat("Fixed effects, sensitivity model:\n")
  print(fixef_sens, row.names = FALSE)
  cat("\n\n")
  
  cat("Type III ANOVA, sensitivity model:\n")
  print(anova_sens, row.names = FALSE)
  cat("\n\n")
  
  cat("Figures saved:\n")
  cat(paste0("- ", file_stub, "_Figure1_moderation_1p5SD.png/pdf\n"))
  cat(paste0("- ", file_stub, "_Figure2_moderation_raw_tertiles.png/pdf\n\n"))
  
  cat("Model assumption checks saved as diagnostic plots:\n")
  cat(paste0("- ", file_stub, "_Diagnostic_residuals_vs_fitted.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_residuals.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_observed_vs_fitted.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_PATNO.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_SITE.png\n"))
  cat(paste0("- ", file_stub, "_Performance_check_model.png\n"))
  
  sink()
  
  invisible(
    list(
      moderator = dat_var,
      model = lmm,
      fixef = fixef_tab,
      anova = anova_tab,
      slopes = slopes_nlr_by_dat_time_df,
      sensitivity_fixef = fixef_sens,
      sensitivity_anova = anova_sens,
      descriptives = desc,
      missingness = missingness_by_variable,
      participant_missingness_summary = participant_missingness_summary
    )
  )
}

# ---------------------------------------------------------
# Run all striatal moderation models
# ---------------------------------------------------------
res_caudate_l <- run_striatal_moderation(
  dat_var = "DAT_caudate_L",
  dat_label = "Left caudate DAT",
  file_stub = "CAUDATE_L"
)

res_caudate_r <- run_striatal_moderation(
  dat_var = "DAT_caudate_R",
  dat_label = "Right caudate DAT",
  file_stub = "CAUDATE_R"
)

res_putamen_l <- run_striatal_moderation(
  dat_var = "DAT_putamen_L",
  dat_label = "Left putaminal DAT",
  file_stub = "PUTAMEN_L"
)

res_putamen_r <- run_striatal_moderation(
  dat_var = "DAT_putamen_R",
  dat_label = "Right putaminal DAT",
  file_stub = "PUTAMEN_R"
)

# ---------------------------------------------------------
# Combined tables across moderators
# ---------------------------------------------------------
combined_fixef <- bind_rows(
  res_caudate_l$fixef %>% mutate(moderator = "CAUDATE_L"),
  res_caudate_r$fixef %>% mutate(moderator = "CAUDATE_R"),
  res_putamen_l$fixef %>% mutate(moderator = "PUTAMEN_L"),
  res_putamen_r$fixef %>% mutate(moderator = "PUTAMEN_R")
) %>%
  select(moderator, everything())

combined_anova <- bind_rows(
  res_caudate_l$anova %>% mutate(moderator = "CAUDATE_L"),
  res_caudate_r$anova %>% mutate(moderator = "CAUDATE_R"),
  res_putamen_l$anova %>% mutate(moderator = "PUTAMEN_L"),
  res_putamen_r$anova %>% mutate(moderator = "PUTAMEN_R")
) %>%
  select(moderator, everything())

combined_slopes <- bind_rows(
  res_caudate_l$slopes %>% mutate(moderator = "CAUDATE_L"),
  res_caudate_r$slopes %>% mutate(moderator = "CAUDATE_R"),
  res_putamen_l$slopes %>% mutate(moderator = "PUTAMEN_L"),
  res_putamen_r$slopes %>% mutate(moderator = "PUTAMEN_R")
) %>%
  select(moderator, everything())

combined_sensitivity_fixef <- bind_rows(
  res_caudate_l$sensitivity_fixef %>% mutate(moderator = "CAUDATE_L"),
  res_caudate_r$sensitivity_fixef %>% mutate(moderator = "CAUDATE_R"),
  res_putamen_l$sensitivity_fixef %>% mutate(moderator = "PUTAMEN_L"),
  res_putamen_r$sensitivity_fixef %>% mutate(moderator = "PUTAMEN_R")
) %>%
  select(moderator, everything())

combined_sensitivity_anova <- bind_rows(
  res_caudate_l$sensitivity_anova %>% mutate(moderator = "CAUDATE_L"),
  res_caudate_r$sensitivity_anova %>% mutate(moderator = "CAUDATE_R"),
  res_putamen_l$sensitivity_anova %>% mutate(moderator = "PUTAMEN_L"),
  res_putamen_r$sensitivity_anova %>% mutate(moderator = "PUTAMEN_R")
) %>%
  select(moderator, everything())

write.csv(
  combined_fixef,
  file.path(out_dir, "COMBINED_moderation_LMM_fixed_effects_all_DAT_moderators.csv"),
  row.names = FALSE
)

write.csv(
  combined_anova,
  file.path(out_dir, "COMBINED_moderation_LMM_TypeIII_ANOVA_all_DAT_moderators.csv"),
  row.names = FALSE
)

write.csv(
  combined_slopes,
  file.path(out_dir, "COMBINED_simple_slopes_NLR_by_DAT_and_time_all_moderators.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_fixef,
  file.path(out_dir, "COMBINED_sensitivity_baseline_UPDRS_adjusted_fixed_effects_all_DAT_moderators.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_anova,
  file.path(out_dir, "COMBINED_sensitivity_baseline_UPDRS_adjusted_TypeIII_ANOVA_all_DAT_moderators.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("NLR x baseline DAT moderation analyses completed.\n")
cat("All outputs were saved to:\n")
cat(
  normalizePath(out_dir, winslash = "/", mustWork = FALSE),
  "\n"
)
cat("============================================================\n")