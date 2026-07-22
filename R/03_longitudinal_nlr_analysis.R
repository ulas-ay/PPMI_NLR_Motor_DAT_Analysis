# ==============================================================================
# 03_longitudinal_nlr_analysis.R
#
# Purpose
# -------
# Compare longitudinal NLR trajectories between healthy controls (HC) and
# participants with Parkinson's disease (PD) using a linear mixed-effects model.
#
# Primary model
# -------------
# NLR ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO)
#
# Random intercepts
# -----------------
# - Participant (PATNO)
# - Study site (SITE)
#
# Expected input
# --------------
# data/derived/PPMI_with_NLR_all_visits_updated.xlsx
# Sheet: All_data_with_NLR
#
# Default output directory
# ------------------------
# outputs/03_longitudinal_nlr_analysis
#
# Optional command-line usage
# ---------------------------
# Rscript R/03_longitudinal_nlr_analysis.R \
#   path/to/PPMI_with_NLR_all_visits_updated.xlsx \
#   path/to/output_directory
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# Package checks
# ------------------------------------------------------------------------------
required_packages <- c(
  "dplyr",
  "readxl",
  "tidyr",
  "ggplot2",
  "lme4",
  "lmerTest",
  "car",
  "emmeans",
  "tibble",
  "performance"
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
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
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

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(
    "data",
    "derived",
    "PPMI_with_NLR_all_visits_updated.xlsx"
  )
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("outputs", "03_longitudinal_nlr_analysis")
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      normalizePath(input_file, winslash = "/", mustWork = FALSE),
      "\n\nRun 01_calculate_nlr_and_prepare_datasets.R first, ",
      "or provide the derived workbook as the first command-line argument."
    ),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Input file: ", normalizePath(input_file, winslash = "/"))
message(
  "Output directory: ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
find_col <- function(data, candidates) {
  detected <- candidates[candidates %in% names(data)]

  if (length(detected) == 0) {
    return(NA_character_)
  }

  detected[[1]]
}

assert_columns <- function(data, columns, object_name = "data") {
  missing_columns <- setdiff(columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing required column(s) in ", object_name, ": ",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

fmt_p <- function(p) {
  if (length(p) == 0 || is.na(p)) {
    return("\u2014")
  }

  if (p < 0.001) {
    return("<0.001")
  }

  sprintf("%.3f", p)
}

standardize_ci_names <- function(data) {
  rename_map <- c(
    "asymp.LCL" = "lower.CL",
    "asymp.UCL" = "upper.CL",
    "lowerCL" = "lower.CL",
    "upperCL" = "upper.CL"
  )

  for (old_name in names(rename_map)) {
    new_name <- rename_map[[old_name]]

    if (
      !new_name %in% names(data) &&
      old_name %in% names(data)
    ) {
      names(data)[names(data) == old_name] <- new_name
    }
  }

  data
}

create_missingness_table <- function(data, variables) {
  bind_rows(
    lapply(
      variables,
      function(variable) {
        tibble(
          variable = variable,
          n_missing = sum(is.na(data[[variable]])),
          pct_missing = 100 * mean(is.na(data[[variable]]))
        )
      }
    )
  )
}

safe_write_csv <- function(data, filename) {
  full_path <- file.path(output_dir, filename)
  write.csv(data, full_path, row.names = FALSE)
  message("Saved: ", full_path)
}

safe_save_check_model <- function(model, filename) {
  output_path <- file.path(output_dir, filename)

  tryCatch(
    {
      diagnostic_check <- performance::check_model(model)

      png(
        filename = output_path,
        width = 2400,
        height = 1800,
        res = 300
      )
      print(diagnostic_check)
      dev.off()

      message("Saved: ", output_path)
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
}

# ------------------------------------------------------------------------------
# Read derived data
# ------------------------------------------------------------------------------
sheet_name <- "All_data_with_NLR"
available_sheets <- excel_sheets(input_file)

if (!sheet_name %in% available_sheets) {
  stop(
    paste0(
      "Worksheet '", sheet_name, "' was not found. Available sheets: ",
      paste(available_sheets, collapse = ", ")
    ),
    call. = FALSE
  )
}

raw_data <- read_excel(
  input_file,
  sheet = sheet_name
)

message("Rows loaded: ", nrow(raw_data))
message("Columns loaded: ", ncol(raw_data))

required_variables <- c(
  "PATNO",
  "EVENT_ID",
  "PRIMDIAG",
  "GROUP",
  "NLR",
  "year_from_baseline"
)

assert_columns(
  raw_data,
  required_variables,
  object_name = "derived NLR dataset"
)

# ------------------------------------------------------------------------------
# Detect covariates
# ------------------------------------------------------------------------------
age_col <- find_col(
  raw_data,
  c(
    "age", "AGE", "Age", "age_bl",
    "AGE_AT_VISIT", "AGE_AT_BASELINE", "age_at_baseline"
  )
)

sex_col <- find_col(
  raw_data,
  c("sex", "SEX", "Sex", "gender", "GENDER")
)

bmi_col <- find_col(
  raw_data,
  c("BMI", "bmi", "Bmi", "body_mass_index", "BodyMassIndex")
)

site_col <- find_col(
  raw_data,
  c("SITE", "site", "Site", "siteid", "SITEID", "site_id")
)

detected_columns <- tibble(
  variable = c("Age", "Sex", "BMI", "Site"),
  detected_column = c(age_col, sex_col, bmi_col, site_col)
)

safe_write_csv(
  detected_columns,
  "Detected_input_columns.csv"
)

if (any(is.na(c(age_col, sex_col, bmi_col, site_col)))) {
  stop(
    paste0(
      "Could not detect one or more required covariate columns: ",
      "age, sex, BMI, or site. See Detected_input_columns.csv."
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Prepare analysis data
# ------------------------------------------------------------------------------
visit_order <- c(
  "BL", "V04", "V06", "V08",
  "V10", "V12", "V13", "V14"
)

analysis_data_before_complete_case <- raw_data %>%
  filter(
    PRIMDIAG %in% c(1, 17),
    EVENT_ID %in% visit_order
  ) %>%
  transmute(
    PATNO = factor(PATNO),
    SITE = factor(.data[[site_col]]),
    GROUP = factor(GROUP, levels = c("HC", "PD")),
    PRIMDIAG = PRIMDIAG,
    EVENT_ID = factor(EVENT_ID, levels = visit_order),
    year_c = suppressWarnings(as.numeric(year_from_baseline)),
    NLR = suppressWarnings(as.numeric(NLR)),
    log_NLR_model = if_else(
      !is.na(NLR) & NLR > 0,
      log(NLR),
      NA_real_
    ),
    age = suppressWarnings(as.numeric(.data[[age_col]])),
    sex = factor(.data[[sex_col]]),
    bmi = suppressWarnings(as.numeric(.data[[bmi_col]]))
  ) %>%
  droplevels()

# ------------------------------------------------------------------------------
# Missingness before complete-case filtering
# ------------------------------------------------------------------------------
model_variables <- c(
  "NLR",
  "GROUP",
  "year_c",
  "age",
  "sex",
  "bmi",
  "SITE",
  "PATNO"
)

missingness_by_variable <- create_missingness_table(
  analysis_data_before_complete_case,
  model_variables
)

missingness_by_visit <- analysis_data_before_complete_case %>%
  group_by(GROUP, EVENT_ID, year_c) %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(PATNO),
    n_NLR_missing = sum(is.na(NLR)),
    pct_NLR_missing = 100 * mean(is.na(NLR)),
    .groups = "drop"
  ) %>%
  arrange(GROUP, year_c)

participant_missingness <- analysis_data_before_complete_case %>%
  group_by(PATNO, GROUP) %>%
  summarise(
    n_rows = n(),
    any_model_variable_missing = any(
      is.na(NLR) |
        is.na(year_c) |
        is.na(age) |
        is.na(sex) |
        is.na(bmi) |
        is.na(SITE)
    ),
    .groups = "drop"
  )

participant_missingness_summary <- participant_missingness %>%
  group_by(GROUP) %>%
  summarise(
    n_participants = n(),
    n_with_at_least_one_missing_model_variable = sum(
      any_model_variable_missing
    ),
    pct_with_at_least_one_missing_model_variable =
      100 * mean(any_model_variable_missing),
    .groups = "drop"
  ) %>%
  bind_rows(
    participant_missingness %>%
      summarise(
        GROUP = factor(
          "Overall",
          levels = c("HC", "PD", "Overall")
        ),
        n_participants = n(),
        n_with_at_least_one_missing_model_variable = sum(
          any_model_variable_missing
        ),
        pct_with_at_least_one_missing_model_variable =
          100 * mean(any_model_variable_missing)
      )
  )

# ------------------------------------------------------------------------------
# Complete-case model dataset
# ------------------------------------------------------------------------------
analysis_data <- analysis_data_before_complete_case %>%
  drop_na(all_of(model_variables)) %>%
  droplevels()

if (nrow(analysis_data) < 10) {
  stop(
    "Fewer than 10 complete observations were available for the model.",
    call. = FALSE
  )
}

if (n_distinct(analysis_data$GROUP) != 2) {
  stop(
    "Both HC and PD groups must be represented in the complete-case dataset.",
    call. = FALSE
  )
}

if (n_distinct(analysis_data$PATNO) < 2) {
  stop(
    "The complete-case dataset contains fewer than two participants.",
    call. = FALSE
  )
}

if (n_distinct(analysis_data$SITE) < 2) {
  stop(
    "The complete-case dataset contains fewer than two study sites.",
    call. = FALSE
  )
}

message("Complete-case rows: ", nrow(analysis_data))
message(
  "Unique participants: ",
  n_distinct(analysis_data$PATNO)
)
message(
  "Unique sites: ",
  n_distinct(analysis_data$SITE)
)

# ------------------------------------------------------------------------------
# Save analysis and missingness datasets
# ------------------------------------------------------------------------------
safe_write_csv(
  analysis_data_before_complete_case,
  "NLR_longitudinal_dataset_before_complete_case_filter.csv"
)

safe_write_csv(
  analysis_data,
  "NLR_longitudinal_complete_case_dataset.csv"
)

safe_write_csv(
  missingness_by_variable,
  "Missingness_by_model_variable.csv"
)

safe_write_csv(
  missingness_by_visit,
  "Missingness_NLR_by_visit_and_group.csv"
)

safe_write_csv(
  participant_missingness,
  "Participant_level_missingness.csv"
)

safe_write_csv(
  participant_missingness_summary,
  "Participant_level_missingness_summary.csv"
)

# ------------------------------------------------------------------------------
# Descriptive summaries
# ------------------------------------------------------------------------------
nlr_descriptives <- analysis_data %>%
  group_by(GROUP, EVENT_ID, year_c) %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(PATNO),
    n_NLR_available = sum(!is.na(NLR)),
    mean_NLR = mean(NLR),
    sd_NLR = sd(NLR),
    median_NLR = median(NLR),
    IQR_NLR = IQR(NLR),
    min_NLR = min(NLR),
    max_NLR = max(NLR),
    .groups = "drop"
  ) %>%
  arrange(GROUP, year_c)

site_counts <- analysis_data %>%
  group_by(SITE) %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(PATNO),
    .groups = "drop"
  ) %>%
  arrange(desc(n_rows))

visit_counts <- analysis_data %>%
  group_by(EVENT_ID, year_c, GROUP) %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(PATNO),
    .groups = "drop"
  ) %>%
  arrange(year_c, GROUP)

safe_write_csv(
  nlr_descriptives,
  "NLR_longitudinal_descriptives.csv"
)

safe_write_csv(
  site_counts,
  "Site_distribution.csv"
)

safe_write_csv(
  visit_counts,
  "Visit_counts_by_group.csv"
)

# ------------------------------------------------------------------------------
# Primary mixed-effects model
#
# Sum-to-zero contrasts are used so Type III tests from car::Anova() are
# interpretable for the GROUP × year interaction in the presence of covariates.
# ------------------------------------------------------------------------------
old_contrasts <- options("contrasts")
on.exit(options(old_contrasts), add = TRUE)

options(contrasts = c("contr.sum", "contr.poly"))

primary_model <- lmer(
  NLR ~ GROUP * year_c + age + sex + bmi +
    (1 | SITE) + (1 | PATNO),
  data = analysis_data,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

model_convergence_messages <- primary_model@optinfo$conv$lme4$messages

convergence_summary <- tibble(
  converged_without_lme4_message = is.null(
    model_convergence_messages
  ),
  convergence_message = if (is.null(model_convergence_messages)) {
    NA_character_
  } else {
    paste(model_convergence_messages, collapse = " | ")
  },
  singular_fit = isSingular(primary_model, tol = 1e-4)
)

safe_write_csv(
  convergence_summary,
  "LMM_convergence_and_singularity_NLR.csv"
)

# ------------------------------------------------------------------------------
# Fixed effects
# ------------------------------------------------------------------------------
fixed_effects <- as.data.frame(
  coef(summary(primary_model))
) %>%
  rownames_to_column("term") %>%
  rename(
    estimate = Estimate,
    SE = `Std. Error`,
    df = df,
    t_value = `t value`,
    p_value = `Pr(>|t|)`
  ) %>%
  mutate(
    lower_95_CI = estimate - 1.96 * SE,
    upper_95_CI = estimate + 1.96 * SE,
    p_formatted = vapply(p_value, fmt_p, character(1))
  ) %>%
  select(
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

safe_write_csv(
  fixed_effects,
  "LMM_fixed_effects_NLR.csv"
)

# ------------------------------------------------------------------------------
# Type III tests
# ------------------------------------------------------------------------------
type3_anova <- as.data.frame(
  car::Anova(
    primary_model,
    type = 3,
    test.statistic = "Chisq"
  )
) %>%
  rownames_to_column("term") %>%
  rename(
    chi_square = Chisq,
    df = Df,
    p_value = `Pr(>Chisq)`
  ) %>%
  mutate(
    p_formatted = vapply(p_value, fmt_p, character(1))
  )

safe_write_csv(
  type3_anova,
  "LMM_TypeIII_ANOVA_NLR.csv"
)

# ------------------------------------------------------------------------------
# Random effects and model fit
# ------------------------------------------------------------------------------
random_effects_variance <- as.data.frame(
  VarCorr(primary_model)
)

model_fit_indices <- tibble(
  AIC = AIC(primary_model),
  BIC = BIC(primary_model),
  logLik = as.numeric(logLik(primary_model)),
  deviance = deviance(primary_model),
  residual_sigma = sigma(primary_model),
  n_observations = nobs(primary_model),
  n_participants = n_distinct(analysis_data$PATNO),
  n_sites = n_distinct(analysis_data$SITE)
)

safe_write_csv(
  random_effects_variance,
  "LMM_random_effects_variance_NLR.csv"
)

safe_write_csv(
  model_fit_indices,
  "LMM_model_fit_indices_NLR.csv"
)

# ------------------------------------------------------------------------------
# Estimated marginal means at observed study years
# ------------------------------------------------------------------------------
time_points <- sort(
  unique(analysis_data$year_c)
)

estimated_marginal_means <- emmeans(
  primary_model,
  ~ GROUP | year_c,
  at = list(year_c = time_points),
  lmer.df = "satterthwaite"
)

estimated_marginal_means_df <- as.data.frame(
  summary(
    estimated_marginal_means,
    infer = c(TRUE, TRUE)
  )
) %>%
  standardize_ci_names() %>%
  mutate(
    p_formatted = if ("p.value" %in% names(.)) {
      vapply(p.value, fmt_p, character(1))
    } else {
      "\u2014"
    }
  )

safe_write_csv(
  estimated_marginal_means_df,
  "EMMeans_NLR_by_group_and_time.csv"
)

# ------------------------------------------------------------------------------
# HC versus PD contrast at each time point
# ------------------------------------------------------------------------------
group_contrasts_each_time <- as.data.frame(
  pairs(
    estimated_marginal_means,
    adjust = "bonferroni"
  )
) %>%
  mutate(
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

safe_write_csv(
  group_contrasts_each_time,
  "Group_comparisons_NLR_at_each_timepoint.csv"
)

# ------------------------------------------------------------------------------
# Group-specific longitudinal slopes
# ------------------------------------------------------------------------------
group_slopes <- emtrends(
  primary_model,
  ~ GROUP,
  var = "year_c",
  lmer.df = "satterthwaite"
)

group_slopes_df <- as.data.frame(
  summary(
    group_slopes,
    infer = c(TRUE, TRUE)
  )
) %>%
  standardize_ci_names() %>%
  mutate(
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

slope_difference <- as.data.frame(
  pairs(
    group_slopes,
    adjust = "none"
  )
) %>%
  mutate(
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

safe_write_csv(
  group_slopes_df,
  "Group_specific_NLR_slopes.csv"
)

safe_write_csv(
  slope_difference,
  "Slope_difference_NLR_between_groups.csv"
)

# ------------------------------------------------------------------------------
# Baseline adjusted group difference
# ------------------------------------------------------------------------------
baseline_emmeans <- emmeans(
  primary_model,
  ~ GROUP,
  at = list(year_c = 0),
  lmer.df = "satterthwaite"
)

baseline_emmeans_df <- as.data.frame(
  summary(
    baseline_emmeans,
    infer = c(TRUE, TRUE)
  )
) %>%
  standardize_ci_names() %>%
  mutate(
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

baseline_group_contrast <- as.data.frame(
  pairs(
    baseline_emmeans,
    adjust = "none"
  )
) %>%
  mutate(
    p_formatted = vapply(p.value, fmt_p, character(1))
  )

safe_write_csv(
  baseline_emmeans_df,
  "Baseline_EMMeans_NLR.csv"
)

safe_write_csv(
  baseline_group_contrast,
  "Baseline_group_difference_NLR.csv"
)

# ------------------------------------------------------------------------------
# Model diagnostics
# ------------------------------------------------------------------------------
diagnostic_data <- analysis_data %>%
  mutate(
    fitted_value = fitted(primary_model),
    residual = resid(primary_model),
    pearson_residual = residual / sigma(primary_model)
  )

safe_write_csv(
  diagnostic_data,
  "LMM_NLR_diagnostic_values.csv"
)

residuals_fitted_plot <- ggplot(
  diagnostic_data,
  aes(x = fitted_value, y = residual)
) +
  geom_point(alpha = 0.45, size = 1.6) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 0.8
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Residuals versus fitted values",
    x = "Fitted values",
    y = "Residuals"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Diagnostic_residuals_vs_fitted_NLR.png"
  ),
  plot = residuals_fitted_plot,
  width = 7,
  height = 5,
  dpi = 600
)

residual_qq_plot <- ggplot(
  diagnostic_data,
  aes(sample = residual)
) +
  stat_qq(alpha = 0.55, size = 1.6) +
  stat_qq_line(linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Q\u2013Q plot of model residuals",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Diagnostic_QQ_residuals_NLR.png"
  ),
  plot = residual_qq_plot,
  width = 7,
  height = 5,
  dpi = 600
)

observed_fitted_plot <- ggplot(
  diagnostic_data,
  aes(x = fitted_value, y = NLR)
) +
  geom_point(alpha = 0.45, size = 1.6) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.7
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Observed versus fitted NLR",
    x = "Fitted values",
    y = "Observed NLR"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Diagnostic_observed_vs_fitted_NLR.png"
  ),
  plot = observed_fitted_plot,
  width = 7,
  height = 5,
  dpi = 600
)

participant_random_effects <- ranef(primary_model)$PATNO %>%
  rownames_to_column("PATNO") %>%
  rename(random_intercept = `(Intercept)`)

participant_random_effects_qq <- ggplot(
  participant_random_effects,
  aes(sample = random_intercept)
) +
  stat_qq(alpha = 0.55, size = 1.6) +
  stat_qq_line(linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Q\u2013Q plot of participant random intercepts",
    x = "Theoretical quantiles",
    y = "Participant random intercepts"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Diagnostic_QQ_random_intercepts_PATNO_NLR.png"
  ),
  plot = participant_random_effects_qq,
  width = 7,
  height = 5,
  dpi = 600
)

site_random_effects <- ranef(primary_model)$SITE %>%
  rownames_to_column("SITE") %>%
  rename(random_intercept = `(Intercept)`)

site_random_effects_qq <- ggplot(
  site_random_effects,
  aes(sample = random_intercept)
) +
  stat_qq(alpha = 0.55, size = 1.6) +
  stat_qq_line(linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Q\u2013Q plot of site random intercepts",
    x = "Theoretical quantiles",
    y = "Site random intercepts"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Diagnostic_QQ_random_intercepts_SITE_NLR.png"
  ),
  plot = site_random_effects_qq,
  width = 7,
  height = 5,
  dpi = 600
)

safe_save_check_model(
  primary_model,
  "Performance_check_model_NLR.png"
)

# ------------------------------------------------------------------------------
# Longitudinal estimated-marginal-means plot
# ------------------------------------------------------------------------------
plot_data <- estimated_marginal_means_df %>%
  mutate(
    year_label = case_when(
      year_c == 0 ~ "BL",
      TRUE ~ as.character(year_c)
    )
  )

observed_summary <- analysis_data %>%
  group_by(GROUP, year_c) %>%
  summarise(
    observed_mean_NLR = mean(NLR),
    observed_se_NLR = sd(NLR) / sqrt(n()),
    n_observations = n(),
    n_participants = n_distinct(PATNO),
    .groups = "drop"
  )

longitudinal_plot <- ggplot(
  plot_data,
  aes(
    x = year_c,
    y = emmean,
    group = GROUP,
    color = GROUP,
    fill = GROUP
  )
) +
  geom_ribbon(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_point(
    data = observed_summary,
    aes(
      x = year_c,
      y = observed_mean_NLR,
      color = GROUP
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    stroke = 1,
    size = 2.2
  ) +
  scale_x_continuous(
    breaks = 0:7,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  labs(
    title = "Longitudinal trajectory of NLR",
    subtitle = paste(
      "Adjusted estimated marginal means with 95% confidence intervals;",
      "open points indicate observed means"
    ),
    x = "Years from baseline",
    y = "Neutrophil-to-lymphocyte ratio",
    color = "Diagnostic group",
    fill = "Diagnostic group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    output_dir,
    "NLR_longitudinal_plot.png"
  ),
  plot = longitudinal_plot,
  width = 8,
  height = 5.5,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "NLR_longitudinal_plot.pdf"
  ),
  plot = longitudinal_plot,
  width = 8,
  height = 5.5
)

# ------------------------------------------------------------------------------
# Observed individual trajectories
# ------------------------------------------------------------------------------
spaghetti_plot <- ggplot(
  analysis_data,
  aes(
    x = year_c,
    y = NLR,
    group = PATNO,
    color = GROUP
  )
) +
  geom_line(alpha = 0.08) +
  geom_smooth(
    aes(group = GROUP, fill = GROUP),
    method = "loess",
    se = TRUE,
    linewidth = 1.2,
    alpha = 0.18
  ) +
  scale_x_continuous(
    breaks = 0:7,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  labs(
    title = "Observed longitudinal NLR trajectories",
    subtitle = paste(
      "Thin lines indicate individual participants;",
      "thick lines indicate smoothed group trends"
    ),
    x = "Years from baseline",
    y = "Neutrophil-to-lymphocyte ratio",
    color = "Diagnostic group",
    fill = "Diagnostic group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    output_dir,
    "NLR_spaghetti_plot.png"
  ),
  plot = spaghetti_plot,
  width = 8,
  height = 5.5,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "NLR_spaghetti_plot.pdf"
  ),
  plot = spaghetti_plot,
  width = 8,
  height = 5.5
)

# ------------------------------------------------------------------------------
# Log-NLR sensitivity model
# ------------------------------------------------------------------------------
log_analysis_data <- analysis_data %>%
  drop_na(log_NLR_model) %>%
  droplevels()

log_nlr_model <- lmer(
  log_NLR_model ~ GROUP * year_c + age + sex + bmi +
    (1 | SITE) + (1 | PATNO),
  data = log_analysis_data,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

log_nlr_fixed_effects <- as.data.frame(
  coef(summary(log_nlr_model))
) %>%
  rownames_to_column("term") %>%
  rename(
    estimate = Estimate,
    SE = `Std. Error`,
    df = df,
    t_value = `t value`,
    p_value = `Pr(>|t|)`
  ) %>%
  mutate(
    lower_95_CI = estimate - 1.96 * SE,
    upper_95_CI = estimate + 1.96 * SE,
    p_formatted = vapply(p_value, fmt_p, character(1))
  )

log_nlr_type3 <- as.data.frame(
  car::Anova(
    log_nlr_model,
    type = 3,
    test.statistic = "Chisq"
  )
) %>%
  rownames_to_column("term") %>%
  rename(
    chi_square = Chisq,
    df = Df,
    p_value = `Pr(>Chisq)`
  ) %>%
  mutate(
    p_formatted = vapply(p_value, fmt_p, character(1))
  )

safe_write_csv(
  log_nlr_fixed_effects,
  "Sensitivity_logNLR_LMM_fixed_effects.csv"
)

safe_write_csv(
  log_nlr_type3,
  "Sensitivity_logNLR_LMM_TypeIII_ANOVA.csv"
)

# ------------------------------------------------------------------------------
# Full text summary
# ------------------------------------------------------------------------------
summary_file <- file.path(
  output_dir,
  "LMM_NLR_summary.txt"
)

summary_connection <- file(summary_file, open = "wt")
sink(summary_connection)

cat("Longitudinal analysis of NLR\n")
cat("============================\n\n")

cat("Input file:\n")
cat(normalizePath(input_file, winslash = "/"), "\n\n")

cat("Output directory:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)

cat("Primary model formula:\n")
cat(
  paste0(
    "NLR ~ GROUP * year_c + age + sex + bmi + ",
    "(1 | SITE) + (1 | PATNO)\n\n"
  )
)

cat("Time coding:\n")
cat("year_c represents years from baseline; baseline = 0.\n\n")

cat("Detected covariate columns:\n")
print(detected_columns, row.names = FALSE)

cat("\nSample before complete-case filtering:\n")
cat(
  "Rows:",
  nrow(analysis_data_before_complete_case),
  "\n"
)
cat(
  "Unique participants:",
  n_distinct(analysis_data_before_complete_case$PATNO),
  "\n"
)
cat(
  "Unique sites:",
  n_distinct(analysis_data_before_complete_case$SITE),
  "\n\n"
)

cat("Complete-case model sample:\n")
cat("Rows:", nrow(analysis_data), "\n")
cat(
  "Unique participants:",
  n_distinct(analysis_data$PATNO),
  "\n"
)
cat(
  "Unique sites:",
  n_distinct(analysis_data$SITE),
  "\n\n"
)

cat("Convergence and singularity:\n")
print(convergence_summary, row.names = FALSE)

cat("\nMissingness by model variable:\n")
print(missingness_by_variable, row.names = FALSE)

cat("\nParticipant-level missingness summary:\n")
print(participant_missingness_summary, row.names = FALSE)

cat("\nNLR descriptives by group and visit:\n")
print(nlr_descriptives, row.names = FALSE)

cat("\nRandom-effects variance:\n")
print(VarCorr(primary_model), comp = c("Variance", "Std.Dev."))

cat("\nModel fit indices:\n")
print(model_fit_indices, row.names = FALSE)

cat("\nFixed effects:\n")
print(fixed_effects, row.names = FALSE)

cat("\nType III tests:\n")
print(type3_anova, row.names = FALSE)

cat("\nEstimated marginal means by group and time:\n")
print(estimated_marginal_means_df, row.names = FALSE)

cat("\nHC versus PD contrasts at each time point:\n")
print(group_contrasts_each_time, row.names = FALSE)

cat("\nGroup-specific NLR slopes:\n")
print(group_slopes_df, row.names = FALSE)

cat("\nSlope difference between groups:\n")
print(slope_difference, row.names = FALSE)

cat("\nBaseline adjusted group difference:\n")
print(baseline_group_contrast, row.names = FALSE)

cat("\nLog-NLR sensitivity model fixed effects:\n")
print(log_nlr_fixed_effects, row.names = FALSE)

cat("\nLog-NLR sensitivity model Type III tests:\n")
print(log_nlr_type3, row.names = FALSE)

sink()
close(summary_connection)

# ------------------------------------------------------------------------------
# Reproducibility information
# ------------------------------------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

# ------------------------------------------------------------------------------
# Final console message
# ------------------------------------------------------------------------------
message("")
message("Longitudinal NLR analysis completed successfully.")
message(
  "Outputs saved to: ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)
