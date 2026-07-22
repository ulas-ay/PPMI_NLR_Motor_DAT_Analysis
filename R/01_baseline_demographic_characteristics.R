# ==============================================================================
# 02_baseline_demographic_characteristics.R
#
# Purpose
# -------
# Generate Table 1 and supporting baseline/follow-up outputs for the PPMI
# NLR, motor progression, and DAT imaging analyses.
#
# Main outputs
# ------------
# - Baseline demographic, clinical, imaging, and inflammatory characteristics
# - Welch t-test and Wilcoxon sensitivity results
# - Normality and chi-square assumption checks
# - Visit-level follow-up/attrition tables
# - Follow-up flow chart and attrition line plot
# - Manuscript-ready sample-size and baseline-summary text
#
# Expected input
# --------------
# data/PPMI_with_NLR_all_visits_updated.xlsx
# Sheet: All_data_with_NLR
#
# Default output directory
# ------------------------
# outputs/02_baseline_demographic_characteristics
#
# Optional command-line usage
# ---------------------------
# Rscript R/02_baseline_demographic_characteristics.R \
#   path/to/input.xlsx path/to/output_directory
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# Package checks
# ------------------------------------------------------------------------------
required_packages <- c(
  "dplyr",
  "readxl",
  "tidyr",
  "stringr",
  "tibble",
  "broom",
  "ggplot2",
  "gt"
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
library(stringr)
library(tibble)
library(broom)
library(ggplot2)
library(gt)

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
  file.path("outputs", "02_baseline_demographic_characteristics")
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      normalizePath(input_file, winslash = "/", mustWork = FALSE),
      "\n\nProvide the file as the first command-line argument or place it at:\n",
      file.path("data", "PPMI_with_NLR_all_visits_updated.xlsx")
    ),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Input file: ", normalizePath(input_file, winslash = "/"))
message("Output directory: ", normalizePath(output_dir, winslash = "/"))

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
  missing <- setdiff(columns, names(data))

  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required column(s) in ", object_name, ": ",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

fmt_mean_sd <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return("\u2014")
  }

  sprintf(
    paste0("%.", digits, "f \u00b1 %.", digits, "f"),
    mean(x),
    stats::sd(x)
  )
}

fmt_n_pct <- function(n, denominator, digits = 1) {
  if (is.na(n) || is.na(denominator) || denominator == 0) {
    return("\u2014")
  }

  sprintf(
    paste0("%d (%.", digits, "f%%)"),
    n,
    100 * n / denominator
  )
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

safe_welch <- function(data, variable) {
  analysis_data <- data %>%
    filter(!is.na(.data[[variable]]), !is.na(GROUP)) %>%
    transmute(
      GROUP = droplevels(GROUP),
      value = suppressWarnings(as.numeric(.data[[variable]]))
    ) %>%
    filter(!is.na(value))

  if (
    n_distinct(analysis_data$GROUP) != 2 ||
    any(table(analysis_data$GROUP) < 2)
  ) {
    return(
      tibble(
        variable = variable,
        test = "Welch t-test",
        statistic = NA_real_,
        degrees_freedom = NA_real_,
        p_value = NA_real_
      )
    )
  }

  test_result <- tryCatch(
    t.test(value ~ GROUP, data = analysis_data, var.equal = FALSE),
    error = function(e) NULL
  )

  if (is.null(test_result)) {
    return(
      tibble(
        variable = variable,
        test = "Welch t-test",
        statistic = NA_real_,
        degrees_freedom = NA_real_,
        p_value = NA_real_
      )
    )
  }

  tibble(
    variable = variable,
    test = "Welch t-test",
    statistic = unname(test_result$statistic),
    degrees_freedom = unname(test_result$parameter),
    p_value = test_result$p.value
  )
}

safe_wilcox <- function(data, variable) {
  analysis_data <- data %>%
    filter(!is.na(.data[[variable]]), !is.na(GROUP)) %>%
    transmute(
      GROUP = droplevels(GROUP),
      value = suppressWarnings(as.numeric(.data[[variable]]))
    ) %>%
    filter(!is.na(value))

  if (
    n_distinct(analysis_data$GROUP) != 2 ||
    any(table(analysis_data$GROUP) < 1)
  ) {
    return(NA_real_)
  }

  tryCatch(
    wilcox.test(value ~ GROUP, data = analysis_data, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_shapiro <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]

  if (length(x) < 3 || length(x) > 5000 || length(unique(x)) < 3) {
    return(c(W = NA_real_, p = NA_real_))
  }

  result <- tryCatch(shapiro.test(x), error = function(e) NULL)

  if (is.null(result)) {
    return(c(W = NA_real_, p = NA_real_))
  }

  c(
    W = unname(result$statistic),
    p = result$p.value
  )
}

normality_checks <- function(data, variables, variable_labels, output_directory) {
  results <- vector("list", length(variables))
  names(results) <- variables

  for (variable in variables) {
    if (!variable %in% names(data)) {
      next
    }

    summary_table <- data %>%
      filter(!is.na(.data[[variable]]), !is.na(GROUP)) %>%
      group_by(GROUP) %>%
      group_modify(
        ~ {
          values <- suppressWarnings(as.numeric(.x[[variable]]))
          values <- values[!is.na(values)]
          shapiro_result <- safe_shapiro(values)

          tibble(
            n = length(values),
            mean = mean(values),
            sd = sd(values),
            median = median(values),
            IQR = IQR(values),
            shapiro_W = shapiro_result[["W"]],
            shapiro_p = shapiro_result[["p"]]
          )
        }
      ) %>%
      ungroup() %>%
      mutate(
        variable = variable,
        label = variable_labels[[variable]],
        shapiro_p_formatted = vapply(shapiro_p, fmt_p, character(1))
      ) %>%
      select(variable, label, everything())

    results[[variable]] <- summary_table

    qq_plot <- ggplot(
      data %>% filter(!is.na(.data[[variable]]), !is.na(GROUP)),
      aes(sample = .data[[variable]])
    ) +
      stat_qq(alpha = 0.55, size = 1.6) +
      stat_qq_line(linewidth = 0.7) +
      facet_wrap(~GROUP, scales = "free") +
      theme_classic(base_size = 13) +
      labs(
        title = paste0("Q\u2013Q plot: ", variable_labels[[variable]]),
        x = "Theoretical quantiles",
        y = "Sample quantiles"
      )

    ggsave(
      filename = file.path(output_directory, paste0("QQ_", variable, ".png")),
      plot = qq_plot,
      width = 7,
      height = 4.5,
      dpi = 600
    )
  }

  bind_rows(results)
}

chisq_check <- function(data, variable, variable_label) {
  analysis_data <- data %>%
    filter(!is.na(GROUP), !is.na(.data[[variable]]))

  contingency_table <- table(
    analysis_data$GROUP,
    analysis_data[[variable]]
  )

  if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
    return(
      tibble(
        variable = variable,
        label = variable_label,
        test_used = "Not tested",
        min_expected = NA_real_,
        pct_expected_less_5 = NA_real_,
        any_expected_less_1 = NA,
        assumption_satisfied = NA,
        chi_square_p = NA_real_,
        fisher_p = NA_real_,
        selected_p = NA_real_,
        selected_p_formatted = "\u2014"
      )
    )
  }

  chi_result <- suppressWarnings(
    chisq.test(contingency_table, correct = FALSE)
  )

  expected_counts <- chi_result$expected
  min_expected <- min(expected_counts)
  pct_expected_less_5 <- 100 * mean(expected_counts < 5)
  any_expected_less_1 <- any(expected_counts < 1)

  assumption_satisfied <-
    min_expected >= 1 &&
    pct_expected_less_5 <= 20

  fisher_p <- tryCatch(
    fisher.test(contingency_table)$p.value,
    error = function(e) NA_real_
  )

  selected_p <- if (assumption_satisfied) {
    chi_result$p.value
  } else {
    fisher_p
  }

  selected_test <- if (assumption_satisfied) {
    "Chi-square"
  } else {
    "Fisher exact"
  }

  tibble(
    variable = variable,
    label = variable_label,
    test_used = selected_test,
    min_expected = min_expected,
    pct_expected_less_5 = pct_expected_less_5,
    any_expected_less_1 = any_expected_less_1,
    assumption_satisfied = assumption_satisfied,
    chi_square_p = chi_result$p.value,
    fisher_p = fisher_p,
    selected_p = selected_p,
    selected_p_formatted = fmt_p(selected_p)
  )
}

safe_gt_save <- function(table_object, filename) {
  tryCatch(
    {
      gt::gtsave(table_object, filename)
      message("Saved: ", filename)
    },
    error = function(e) {
      warning(
        paste0(
          "Could not save ", filename, ": ", conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
}

# ------------------------------------------------------------------------------
# Read data
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

df <- read_excel(input_file, sheet = sheet_name)

message("Rows loaded: ", nrow(df))
message("Columns loaded: ", ncol(df))

assert_columns(
  df,
  c("PATNO", "PRIMDIAG", "EVENT_ID"),
  object_name = "input dataset"
)

# ------------------------------------------------------------------------------
# Detect variable names
# ------------------------------------------------------------------------------
age_col <- find_col(
  df,
  c(
    "age", "AGE", "Age", "age_bl", "AGE_AT_VISIT",
    "AGE_AT_BASELINE", "age_at_baseline"
  )
)

sex_col <- find_col(
  df,
  c("sex", "SEX", "Sex", "gender", "GENDER")
)

bmi_col <- find_col(
  df,
  c("BMI", "bmi", "Bmi", "body_mass_index", "BodyMassIndex")
)

duration_col <- find_col(
  df,
  c(
    "duration_yrs", "disease_duration",
    "disease_duration_years", "DURATION_YRS"
  )
)

updrs_on_col <- find_col(
  df,
  c(
    "updrs3_score_on", "UPDRS3_score_on",
    "updrs3_on", "UPDRS3_ON"
  )
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

education_col <- find_col(
  df,
  c(
    "EDUCYRS", "education", "Education",
    "education_years", "educyrs", "EDUCATION"
  )
)

moca_col <- find_col(
  df,
  c("MOCA", "MoCA", "moca", "moca_score", "MCA")
)

upsit_col <- find_col(
  df,
  c("UPSIT", "upsit", "UPSIT_score", "upsit_score")
)

caudate_l_col <- find_col(
  df,
  c("MIA_CAUDATE_L", "mia_caudate_l")
)

caudate_r_col <- find_col(
  df,
  c("MIA_CAUDATE_R", "mia_caudate_r")
)

putamen_l_col <- find_col(
  df,
  c("MIA_PUTAMEN_L", "mia_putamen_l")
)

putamen_r_col <- find_col(
  df,
  c("MIA_PUTAMEN_R", "mia_putamen_r")
)

detected_columns <- tibble(
  variable = c(
    "Age", "Sex", "BMI", "Disease duration", "UPDRS-III ON",
    "LEDD", "Dominant side", "Education", "MoCA", "UPSIT",
    "Left caudate DAT", "Right caudate DAT",
    "Left putamen DAT", "Right putamen DAT"
  ),
  detected_column = c(
    age_col, sex_col, bmi_col, duration_col, updrs_on_col,
    ledd_col, domside_col, education_col, moca_col, upsit_col,
    caudate_l_col, caudate_r_col, putamen_l_col, putamen_r_col
  )
)

write.csv(
  detected_columns,
  file.path(output_dir, "Detected_input_columns.csv"),
  row.names = FALSE
)

if (any(is.na(c(age_col, sex_col, bmi_col)))) {
  stop(
    paste0(
      "Could not detect the required age, sex, or BMI column. ",
      "See Detected_input_columns.csv and inspect names(df)."
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Keep PD and healthy controls
# PRIMDIAG: 1 = PD; 17 = healthy control
# ------------------------------------------------------------------------------
df_analysis <- df %>%
  filter(PRIMDIAG %in% c(1, 17)) %>%
  mutate(
    GROUP = factor(
      if_else(PRIMDIAG == 1, "PD", "HC"),
      levels = c("HC", "PD")
    ),
    EVENT_ID = as.character(EVENT_ID)
  )

# ------------------------------------------------------------------------------
# Visit coding
# ------------------------------------------------------------------------------
visit_order <- c("BL", "V04", "V06", "V08", "V10", "V12", "V13", "V14")

visit_dictionary <- tibble(
  EVENT_ID = visit_order,
  visit_label = c(
    "Baseline", "Year 1", "Year 2", "Year 3",
    "Year 4", "Year 5", "Year 6", "Year 7"
  ),
  year_c = 0:7
)

df_analysis <- df_analysis %>%
  left_join(visit_dictionary, by = "EVENT_ID")

# ------------------------------------------------------------------------------
# Baseline dataset
# ------------------------------------------------------------------------------
baseline <- df_analysis %>%
  filter(EVENT_ID == "BL") %>%
  distinct(PATNO, .keep_all = TRUE) %>%
  rename(
    age = all_of(age_col),
    SEX_raw = all_of(sex_col),
    BMI = all_of(bmi_col)
  )

if (nrow(baseline) == 0) {
  stop("No baseline records with EVENT_ID == 'BL' were found.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Recode categorical variables
# ------------------------------------------------------------------------------
baseline <- baseline %>%
  mutate(
    SEX = case_when(
      as.character(SEX_raw) %in% c(
        "0", "Female", "F", "female", "f"
      ) ~ "Female",
      as.character(SEX_raw) %in% c(
        "1", "Male", "M", "male", "m"
      ) ~ "Male",
      TRUE ~ as.character(SEX_raw)
    ),
    SEX = factor(SEX, levels = c("Female", "Male"))
  )

if (!is.na(domside_col)) {
  baseline <- baseline %>%
    rename(DOMSIDE_raw = all_of(domside_col)) %>%
    mutate(
      DOMSIDE = case_when(
        as.character(DOMSIDE_raw) %in% c(
          "1", "Left", "left", "L"
        ) ~ "Left",
        as.character(DOMSIDE_raw) %in% c(
          "2", "Right", "right", "R"
        ) ~ "Right",
        as.character(DOMSIDE_raw) %in% c(
          "3", "Symmetric", "symmetric", "Bilateral"
        ) ~ "Symmetric",
        TRUE ~ as.character(DOMSIDE_raw)
      ),
      DOMSIDE = factor(
        DOMSIDE,
        levels = c("Left", "Right", "Symmetric")
      )
    )
}

# ------------------------------------------------------------------------------
# Baseline NLR
# ------------------------------------------------------------------------------
if ("baseline_NLR" %in% names(baseline)) {
  baseline <- baseline %>%
    mutate(NLR_table = suppressWarnings(as.numeric(baseline_NLR)))
} else if ("NLR" %in% names(baseline)) {
  baseline <- baseline %>%
    mutate(NLR_table = suppressWarnings(as.numeric(NLR)))
} else if (all(c("Neutrophils", "Lymphocytes") %in% names(baseline))) {
  baseline <- baseline %>%
    mutate(
      NLR_table = if_else(
        !is.na(Lymphocytes) & Lymphocytes > 0,
        as.numeric(Neutrophils) / as.numeric(Lymphocytes),
        NA_real_
      )
    )
} else {
  warning(
    "No baseline NLR variable could be identified or calculated.",
    call. = FALSE
  )
  baseline$NLR_table <- NA_real_
}

# ------------------------------------------------------------------------------
# Baseline LEDD
#
# The cohort is de novo and drug-naive at baseline. Therefore, baseline LEDD is
# set to 0 for PD participants and left missing for controls.
# ------------------------------------------------------------------------------
if (!is.na(ledd_col)) {
  baseline <- baseline %>%
    rename(LEDD_raw = all_of(ledd_col)) %>%
    mutate(
      LEDD_table = if_else(GROUP == "PD", 0, NA_real_)
    )
}

# ------------------------------------------------------------------------------
# Rename optional variables
# ------------------------------------------------------------------------------
optional_rename_map <- c(
  duration_yrs = duration_col,
  UPDRSIII_ON = updrs_on_col,
  education_years = education_col,
  MoCA = moca_col,
  UPSIT = upsit_col,
  DAT_caudate_L = caudate_l_col,
  DAT_caudate_R = caudate_r_col,
  DAT_putamen_L = putamen_l_col,
  DAT_putamen_R = putamen_r_col
)

for (new_name in names(optional_rename_map)) {
  old_name <- optional_rename_map[[new_name]]

  if (!is.na(old_name) && old_name %in% names(baseline)) {
    names(baseline)[names(baseline) == old_name] <- new_name
  }
}

# ------------------------------------------------------------------------------
# Save baseline analysis dataset
# ------------------------------------------------------------------------------
write.csv(
  baseline,
  file.path(output_dir, "Baseline_dataset_for_Table1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Variables and labels
# ------------------------------------------------------------------------------
comparison_continuous_variables <- intersect(
  c(
    "age",
    "BMI",
    "education_years",
    "MoCA",
    "UPSIT",
    "Neutrophils",
    "Lymphocytes",
    "Monocytes",
    "NLR_table"
  ),
  names(baseline)
)

comparison_categorical_variables <- intersect(
  "SEX",
  names(baseline)
)

pd_only_continuous_variables <- intersect(
  c(
    "duration_yrs",
    "UPDRSIII_ON",
    "LEDD_table",
    "DAT_caudate_L",
    "DAT_caudate_R",
    "DAT_putamen_L",
    "DAT_putamen_R"
  ),
  names(baseline)
)

pd_only_categorical_variables <- intersect(
  "DOMSIDE",
  names(baseline)
)

variable_labels <- list(
  age = "Age, years",
  BMI = "Body mass index, kg/m\u00b2",
  education_years = "Education, years",
  MoCA = "MoCA score",
  UPSIT = "UPSIT score",
  Neutrophils = "Neutrophils",
  Lymphocytes = "Lymphocytes",
  Monocytes = "Monocytes",
  NLR_table = "Neutrophil-to-lymphocyte ratio",
  SEX = "Male sex, n (%)",
  duration_yrs = "Disease duration, years",
  UPDRSIII_ON = "ON-medication MDS-UPDRS Part III score",
  LEDD_table = "LEDD at baseline",
  DOMSIDE = "Dominant side of motor involvement",
  DAT_caudate_L = "Left caudate DAT binding",
  DAT_caudate_R = "Right caudate DAT binding",
  DAT_putamen_L = "Left putamen DAT binding",
  DAT_putamen_R = "Right putamen DAT binding"
)

# ------------------------------------------------------------------------------
# Sample sizes
# ------------------------------------------------------------------------------
n_hc <- sum(baseline$GROUP == "HC", na.rm = TRUE)
n_pd <- sum(baseline$GROUP == "PD", na.rm = TRUE)
n_overall <- nrow(baseline)

# ------------------------------------------------------------------------------
# Continuous comparison rows
# ------------------------------------------------------------------------------
continuous_rows <- lapply(
  comparison_continuous_variables,
  function(variable) {
    welch_result <- safe_welch(baseline, variable)
    wilcoxon_p <- safe_wilcox(baseline, variable)

    tibble(
      Section = paste(
        "Baseline demographic and peripheral inflammatory",
        "characteristics"
      ),
      Variable = variable_labels[[variable]],
      Overall = fmt_mean_sd(baseline[[variable]]),
      HC = fmt_mean_sd(
        baseline %>%
          filter(GROUP == "HC") %>%
          pull(all_of(variable))
      ),
      PD = fmt_mean_sd(
        baseline %>%
          filter(GROUP == "PD") %>%
          pull(all_of(variable))
      ),
      p = fmt_p(welch_result$p_value),
      Test = "Welch t-test",
      Wilcoxon_sensitivity_p = fmt_p(wilcoxon_p)
    )
  }
) %>%
  bind_rows()

# ------------------------------------------------------------------------------
# Categorical comparison rows
# ------------------------------------------------------------------------------
categorical_rows <- tibble()
chisq_outputs <- tibble()

if ("SEX" %in% comparison_categorical_variables) {
  sex_test <- chisq_check(
    baseline,
    variable = "SEX",
    variable_label = variable_labels[["SEX"]]
  )

  chisq_outputs <- sex_test

  denominator_overall <- sum(!is.na(baseline$SEX))
  denominator_hc <- sum(
    baseline$GROUP == "HC" & !is.na(baseline$SEX)
  )
  denominator_pd <- sum(
    baseline$GROUP == "PD" & !is.na(baseline$SEX)
  )

  n_male_overall <- sum(baseline$SEX == "Male", na.rm = TRUE)
  n_male_hc <- sum(
    baseline$GROUP == "HC" & baseline$SEX == "Male",
    na.rm = TRUE
  )
  n_male_pd <- sum(
    baseline$GROUP == "PD" & baseline$SEX == "Male",
    na.rm = TRUE
  )

  categorical_rows <- tibble(
    Section = paste(
      "Baseline demographic and peripheral inflammatory",
      "characteristics"
    ),
    Variable = variable_labels[["SEX"]],
    Overall = fmt_n_pct(n_male_overall, denominator_overall),
    HC = fmt_n_pct(n_male_hc, denominator_hc),
    PD = fmt_n_pct(n_male_pd, denominator_pd),
    p = fmt_p(sex_test$selected_p),
    Test = sex_test$test_used,
    Wilcoxon_sensitivity_p = "\u2014"
  )
}

# ------------------------------------------------------------------------------
# PD-only continuous rows
# ------------------------------------------------------------------------------
pd_continuous_rows <- lapply(
  pd_only_continuous_variables,
  function(variable) {
    tibble(
      Section = "PD-specific baseline clinical and imaging characteristics",
      Variable = variable_labels[[variable]],
      Overall = "\u2014",
      HC = "\u2014",
      PD = fmt_mean_sd(
        baseline %>%
          filter(GROUP == "PD") %>%
          pull(all_of(variable))
      ),
      p = "\u2014",
      Test = "PD only",
      Wilcoxon_sensitivity_p = "\u2014"
    )
  }
) %>%
  bind_rows()

# ------------------------------------------------------------------------------
# PD-only dominant-side row
# ------------------------------------------------------------------------------
pd_categorical_rows <- tibble()

if ("DOMSIDE" %in% pd_only_categorical_variables) {
  pd_dominant_side <- baseline %>%
    filter(GROUP == "PD", !is.na(DOMSIDE)) %>%
    count(DOMSIDE, name = "n") %>%
    mutate(
      percentage = 100 * n / sum(n),
      text = paste0(
        as.character(DOMSIDE),
        ": ",
        n,
        " (",
        sprintf("%.1f", percentage),
        "%)"
      )
    )

  pd_categorical_rows <- tibble(
    Section = "PD-specific baseline clinical and imaging characteristics",
    Variable = variable_labels[["DOMSIDE"]],
    Overall = "\u2014",
    HC = "\u2014",
    PD = paste(pd_dominant_side$text, collapse = "; "),
    p = "\u2014",
    Test = "PD only",
    Wilcoxon_sensitivity_p = "\u2014"
  )
}

# ------------------------------------------------------------------------------
# Final Table 1
# ------------------------------------------------------------------------------
sample_size_row <- tibble(
  Section = paste(
    "Baseline demographic and peripheral inflammatory",
    "characteristics"
  ),
  Variable = "N",
  Overall = as.character(n_overall),
  HC = as.character(n_hc),
  PD = as.character(n_pd),
  p = "\u2014",
  Test = "\u2014",
  Wilcoxon_sensitivity_p = "\u2014"
)

table1_data <- bind_rows(
  sample_size_row,
  continuous_rows,
  categorical_rows,
  pd_continuous_rows,
  pd_categorical_rows
)

write.csv(
  table1_data,
  file.path(
    output_dir,
    "Table1_Baseline_Demographic_Clinical_Characteristics.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Formatted Table 1
# ------------------------------------------------------------------------------
table1_gt <- table1_data %>%
  select(Section, Variable, Overall, HC, PD, p) %>%
  gt(groupname_col = "Section") %>%
  tab_header(
    title = md(
      paste0(
        "**Table 1. Baseline demographic, clinical, imaging, ",
        "and peripheral inflammatory characteristics**"
      )
    )
  ) %>%
  cols_label(
    Variable = md("**Variable**"),
    Overall = md("**Overall**"),
    HC = md("**HC**"),
    PD = md("**PD**"),
    p = md("**p**")
  ) %>%
  tab_spanner(
    label = md("**Diagnostic group**"),
    columns = c(HC, PD)
  ) %>%
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(4),
    row_group.font.weight = "bold"
  ) %>%
  tab_source_note(
    source_note = md(
      paste0(
        "Continuous variables are presented as mean \u00b1 SD, and categorical ",
        "variables as n (%). Between-group comparisons used Welch's t-tests ",
        "for continuous variables. Categorical variables were compared using ",
        "Pearson's chi-square test when expected-cell assumptions were met and ",
        "Fisher's exact test otherwise. PD-specific clinical and imaging ",
        "variables are summarized for the PD group only."
      )
    )
  )

safe_gt_save(
  table1_gt,
  file.path(
    output_dir,
    "Table1_Baseline_Demographic_Clinical_Characteristics.html"
  )
)

# DOCX export requires a sufficiently recent gt installation and Pandoc.
safe_gt_save(
  table1_gt,
  file.path(
    output_dir,
    "Table1_Baseline_Demographic_Clinical_Characteristics.docx"
  )
)

# ------------------------------------------------------------------------------
# Distributional checks and sensitivity analyses
# ------------------------------------------------------------------------------
normality_results <- normality_checks(
  data = baseline,
  variables = comparison_continuous_variables,
  variable_labels = variable_labels,
  output_directory = output_dir
)

write.csv(
  normality_results,
  file.path(
    output_dir,
    "Welch_ttest_normality_assumption_checks.csv"
  ),
  row.names = FALSE
)

welch_summary <- lapply(
  comparison_continuous_variables,
  function(variable) {
    welch_result <- safe_welch(baseline, variable)
    wilcoxon_p <- safe_wilcox(baseline, variable)

    tibble(
      variable = variable,
      label = variable_labels[[variable]],
      welch_t = welch_result$statistic,
      welch_df = welch_result$degrees_freedom,
      welch_p = welch_result$p_value,
      welch_p_formatted = fmt_p(welch_result$p_value),
      wilcoxon_sensitivity_p = wilcoxon_p,
      wilcoxon_sensitivity_p_formatted = fmt_p(wilcoxon_p)
    )
  }
) %>%
  bind_rows()

write.csv(
  welch_summary,
  file.path(
    output_dir,
    "Welch_ttest_pvalues_with_Wilcoxon_sensitivity.csv"
  ),
  row.names = FALSE
)

write.csv(
  chisq_outputs,
  file.path(
    output_dir,
    "Chi_square_expected_cell_assumption_checks.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Follow-up and attrition counts
# ------------------------------------------------------------------------------
followup_data <- df_analysis %>%
  filter(EVENT_ID %in% visit_order) %>%
  mutate(
    EVENT_ID = factor(EVENT_ID, levels = visit_order),
    GROUP = factor(GROUP, levels = c("HC", "PD"))
  )

if ("NLR" %in% names(followup_data)) {
  followup_data <- followup_data %>%
    mutate(NLR_visit = suppressWarnings(as.numeric(NLR)))
} else if (
  all(c("Neutrophils", "Lymphocytes") %in% names(followup_data))
) {
  followup_data <- followup_data %>%
    mutate(
      NLR_visit = if_else(
        !is.na(Lymphocytes) & Lymphocytes > 0,
        as.numeric(Neutrophils) / as.numeric(Lymphocytes),
        NA_real_
      )
    )
} else {
  followup_data <- followup_data %>%
    mutate(NLR_visit = NA_real_)
}

if (!is.na(updrs_on_col) && updrs_on_col %in% names(followup_data)) {
  followup_data <- followup_data %>%
    mutate(
      UPDRSIII_ON_visit = suppressWarnings(
        as.numeric(.data[[updrs_on_col]])
      )
    )
} else {
  followup_data <- followup_data %>%
    mutate(UPDRSIII_ON_visit = NA_real_)
}

dat_columns_available <- na.omit(
  c(
    caudate_l_col,
    caudate_r_col,
    putamen_l_col,
    putamen_r_col
  )
) %>%
  as.character()

if (length(dat_columns_available) > 0) {
  followup_data <- followup_data %>%
    mutate(
      DAT_any_available = rowSums(
        !is.na(pick(all_of(dat_columns_available)))
      ) > 0
    )
} else {
  followup_data <- followup_data %>%
    mutate(DAT_any_available = FALSE)
}

visit_counts <- followup_data %>%
  group_by(EVENT_ID, visit_label, year_c, GROUP) %>%
  summarise(
    n_subjects_any_record = n_distinct(PATNO),
    n_subjects_with_NLR = n_distinct(PATNO[!is.na(NLR_visit)]),
    n_subjects_with_ON_UPDRSIII = n_distinct(
      PATNO[!is.na(UPDRSIII_ON_visit)]
    ),
    n_subjects_with_any_DAT = n_distinct(
      PATNO[DAT_any_available]
    ),
    .groups = "drop"
  ) %>%
  arrange(GROUP, year_c)

write.csv(
  visit_counts,
  file.path(output_dir, "Followup_attrition_counts_long.csv"),
  row.names = FALSE
)

visit_counts_wide <- visit_counts %>%
  select(
    EVENT_ID,
    visit_label,
    year_c,
    GROUP,
    n_subjects_any_record
  ) %>%
  pivot_wider(
    names_from = GROUP,
    values_from = n_subjects_any_record
  ) %>%
  arrange(year_c)

write.csv(
  visit_counts_wide,
  file.path(
    output_dir,
    "Followup_attrition_counts_wide_any_record.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Follow-up flow chart
# ------------------------------------------------------------------------------
flow_data <- visit_counts %>%
  select(
    EVENT_ID,
    visit_label,
    year_c,
    GROUP,
    n_subjects_any_record
  ) %>%
  mutate(
    x = as.numeric(factor(EVENT_ID, levels = visit_order)),
    y = if_else(GROUP == "HC", 2, 1),
    box_label = paste0(
      visit_label,
      "\n",
      GROUP,
      " n = ",
      n_subjects_any_record
    )
  )

arrow_data <- flow_data %>%
  arrange(GROUP, x) %>%
  group_by(GROUP) %>%
  mutate(
    xend = lead(x),
    yend = lead(y)
  ) %>%
  filter(!is.na(xend)) %>%
  ungroup()

flow_plot <- ggplot() +
  geom_segment(
    data = arrow_data,
    aes(
      x = x + 0.38,
      xend = xend - 0.38,
      y = y,
      yend = yend
    ),
    arrow = grid::arrow(length = grid::unit(0.15, "inches")),
    linewidth = 0.5
  ) +
  geom_label(
    data = flow_data,
    aes(x = x, y = y, label = box_label),
    size = 3.2,
    linewidth = 0.35,
    fill = "white"
  ) +
  scale_x_continuous(
    breaks = seq_along(visit_order),
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7"),
    limits = c(0.5, length(visit_order) + 0.5)
  ) +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("PD", "HC"),
    limits = c(0.4, 2.6)
  ) +
  labs(
    title = "Follow-up and attrition across longitudinal visits",
    subtitle = "Participants with any available record at each visit",
    x = "Years from baseline",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  file.path(output_dir, "Followup_attrition_flowchart.png"),
  plot = flow_plot,
  width = 13,
  height = 4.8,
  dpi = 600
)

ggsave(
  file.path(output_dir, "Followup_attrition_flowchart.pdf"),
  plot = flow_plot,
  width = 13,
  height = 4.8
)

# ------------------------------------------------------------------------------
# Attrition line plot
# ------------------------------------------------------------------------------
attrition_plot <- ggplot(
  visit_counts,
  aes(
    x = year_c,
    y = n_subjects_any_record,
    group = GROUP,
    color = GROUP
  )
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = n_subjects_any_record),
    vjust = -0.8,
    size = 3.3,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = 0:7,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  labs(
    title = "Longitudinal follow-up availability",
    subtitle = "Participants with any available record at each visit",
    x = "Years from baseline",
    y = "Number of participants",
    color = "Group"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "Followup_attrition_lineplot.png"),
  plot = attrition_plot,
  width = 8.5,
  height = 5.5,
  dpi = 600
)

ggsave(
  file.path(output_dir, "Followup_attrition_lineplot.pdf"),
  plot = attrition_plot,
  width = 8.5,
  height = 5.5
)

# ------------------------------------------------------------------------------
# Manuscript-ready text
# ------------------------------------------------------------------------------
sample_size_text <- c(
  "Sample size justification",
  "=========================",
  "",
  paste0(
    "Baseline analyses included all eligible PPMI participants with a ",
    "diagnostic classification of PD or healthy control and an available ",
    "baseline visit record. This yielded ",
    n_pd,
    " participants with PD and ",
    n_hc,
    " healthy controls."
  ),
  "",
  paste0(
    "No formal a priori power calculation was performed because this was a ",
    "secondary analysis of an existing longitudinal observational cohort. ",
    "The analytic sample size was therefore determined by the number of ",
    "eligible participants with available blood-count data, clinical ",
    "assessments, imaging measures, and covariates. Effect estimates should ",
    "be reported with 95% confidence intervals and exact p values, while ",
    "missingness and follow-up availability are summarized at both the ",
    "observation and participant levels."
  )
)

writeLines(
  sample_size_text,
  file.path(output_dir, "Sample_size_justification_text.txt")
)

baseline_summary_text <- c(
  "Baseline demographic and clinical characteristics",
  "=================================================",
  "",
  paste0(
    "Baseline analyses included ",
    n_pd,
    " participants with PD and ",
    n_hc,
    " healthy controls."
  ),
  "",
  paste0(
    "Baseline demographic, clinical, imaging, and peripheral inflammatory ",
    "characteristics are summarized in Table 1. Continuous variables were ",
    "compared using Welch's t-tests. Distributional characteristics were ",
    "examined using Q\u2013Q plots and Shapiro\u2013Wilk tests, and Wilcoxon rank-sum ",
    "tests were conducted as sensitivity analyses. Categorical variables ",
    "were compared using Pearson's chi-square test when expected-cell ",
    "assumptions were met and Fisher's exact test otherwise. Follow-up ",
    "availability and attrition are summarized in the accompanying tables ",
    "and figures."
  )
)

writeLines(
  baseline_summary_text,
  file.path(output_dir, "Manuscript_baseline_summary_starter.txt")
)

# ------------------------------------------------------------------------------
# Reproducibility information
# ------------------------------------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

# ------------------------------------------------------------------------------
# Console summary
# ------------------------------------------------------------------------------
message("")
message("Baseline sample sizes")
message("---------------------")
message("HC: ", n_hc)
message("PD: ", n_pd)
message("Overall: ", n_overall)

message("")
message(
  "Continuous comparison variables: ",
  paste(comparison_continuous_variables, collapse = ", ")
)

message(
  "PD-only continuous variables: ",
  paste(pd_only_continuous_variables, collapse = ", ")
)

message("")
message("Analysis completed successfully.")
message(
  "Outputs saved to: ",
  normalizePath(output_dir, winslash = "/")
)
