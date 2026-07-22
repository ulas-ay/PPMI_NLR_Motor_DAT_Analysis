# ==============================================================================
# 01_calculate_nlr_and_prepare_datasets.R
#
# Purpose
# -------
# Calculate visit-wise and baseline inflammatory ratios from merged PPMI blood
# data and create analysis-ready datasets for subsequent scripts.
#
# Derived variables
# -----------------
# - NLR: neutrophil-to-lymphocyte ratio
# - MLR: monocyte-to-lymphocyte ratio
# - PLR: platelet-to-lymphocyte ratio
# - Baseline versions of each ratio and blood-cell measure
# - Visit-wise and baseline-standardized NLR variables
# - Natural-log-transformed NLR variables
#
# Expected input
# --------------
# data/raw/PPMI_with_serum_all_visits_all_blood_merged.xlsx
#
# Main derived dataset
# --------------------
# data/derived/PPMI_with_NLR_all_visits_updated.xlsx
#
# Supporting outputs
# ------------------
# outputs/01_nlr_calculation
#
# Optional command-line usage
# ---------------------------
# Rscript R/01_calculate_nlr_and_prepare_datasets.R \
#   path/to/raw_input.xlsx \
#   path/to/derived_output.xlsx \
#   path/to/supporting_output_directory
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# Package checks
# ------------------------------------------------------------------------------
required_packages <- c(
  "readxl",
  "openxlsx",
  "dplyr",
  "tidyr",
  "ggplot2",
  "tibble"
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

library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(
    "data",
    "raw",
    "PPMI_with_serum_all_visits_all_blood_merged.xlsx"
  )
}

derived_file <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(
    "data",
    "derived",
    "PPMI_with_NLR_all_visits_updated.xlsx"
  )
}

output_dir <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path("outputs", "01_nlr_calculation")
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      normalizePath(input_file, winslash = "/", mustWork = FALSE),
      "\n\nProvide the raw Excel file as the first command-line argument or ",
      "place it at the default location shown above."
    ),
    call. = FALSE
  )
}

dir.create(dirname(derived_file), recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Input file: ", normalizePath(input_file, winslash = "/"))
message(
  "Derived workbook: ",
  normalizePath(derived_file, winslash = "/", mustWork = FALSE)
)
message(
  "Supporting outputs: ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
assert_columns <- function(data, required_columns, object_name = "data") {
  missing_columns <- setdiff(required_columns, names(data))

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

safe_ratio <- function(numerator, denominator) {
  numerator <- suppressWarnings(as.numeric(numerator))
  denominator <- suppressWarnings(as.numeric(denominator))

  ratio <- ifelse(
    !is.na(numerator) &
      !is.na(denominator) &
      denominator > 0,
    numerator / denominator,
    NA_real_
  )

  ratio[!is.finite(ratio)] <- NA_real_
  ratio
}

safe_standardize <- function(x, reference_mean = NULL, reference_sd = NULL) {
  x <- suppressWarnings(as.numeric(x))

  if (is.null(reference_mean)) {
    reference_mean <- mean(x, na.rm = TRUE)
  }

  if (is.null(reference_sd)) {
    reference_sd <- sd(x, na.rm = TRUE)
  }

  if (
    is.na(reference_mean) ||
    is.na(reference_sd) ||
    reference_sd == 0
  ) {
    return(rep(NA_real_, length(x)))
  }

  as.numeric((x - reference_mean) / reference_sd)
}

safe_log <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(x) & x > 0, log(x), NA_real_)
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

safe_welch <- function(data, value_variable, group_variable = "GROUP") {
  analysis_data <- data %>%
    select(all_of(c(value_variable, group_variable))) %>%
    drop_na() %>%
    mutate(
      "{value_variable}" := suppressWarnings(
        as.numeric(.data[[value_variable]])
      ),
      "{group_variable}" := droplevels(
        factor(.data[[group_variable]])
      )
    ) %>%
    filter(!is.na(.data[[value_variable]]))

  group_counts <- table(analysis_data[[group_variable]])

  if (
    length(group_counts) != 2 ||
    any(group_counts < 2)
  ) {
    return(
      tibble(
        variable = value_variable,
        n = nrow(analysis_data),
        statistic = NA_real_,
        df = NA_real_,
        p = NA_real_
      )
    )
  }

  formula_object <- reformulate(
    group_variable,
    response = value_variable
  )

  test_result <- tryCatch(
    t.test(
      formula_object,
      data = analysis_data,
      var.equal = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(test_result)) {
    return(
      tibble(
        variable = value_variable,
        n = nrow(analysis_data),
        statistic = NA_real_,
        df = NA_real_,
        p = NA_real_
      )
    )
  }

  tibble(
    variable = value_variable,
    n = nrow(analysis_data),
    statistic = unname(test_result$statistic),
    df = unname(test_result$parameter),
    p = test_result$p.value
  )
}

safe_descriptive_summary <- function(data, variable) {
  values <- suppressWarnings(as.numeric(data[[variable]]))
  available_values <- values[!is.na(values)]

  if (length(available_values) == 0) {
    return(
      tibble(
        n_rows = nrow(data),
        n_available = 0L,
        mean = NA_real_,
        sd = NA_real_,
        median = NA_real_,
        IQR = NA_real_,
        min = NA_real_,
        max = NA_real_
      )
    )
  }

  tibble(
    n_rows = nrow(data),
    n_available = length(available_values),
    mean = mean(available_values),
    sd = sd(available_values),
    median = median(available_values),
    IQR = IQR(available_values),
    min = min(available_values),
    max = max(available_values)
  )
}

create_missingness_table <- function(data, variables) {
  variables <- intersect(variables, names(data))

  if (length(variables) == 0) {
    return(
      tibble(
        variable = character(),
        n_missing = integer(),
        pct_missing = numeric()
      )
    )
  }

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

write_csv_output <- function(data, filename) {
  full_path <- file.path(output_dir, filename)
  write.csv(data, full_path, row.names = FALSE)
  message("Saved: ", full_path)
}

# ------------------------------------------------------------------------------
# Read raw data
# ------------------------------------------------------------------------------
raw_data <- read_excel(input_file)

message("Rows loaded: ", nrow(raw_data))
message("Columns loaded: ", ncol(raw_data))

required_variables <- c(
  "PATNO",
  "EVENT_ID",
  "PRIMDIAG",
  "Neutrophils",
  "Lymphocytes",
  "Monocytes",
  "Platelets"
)

assert_columns(
  raw_data,
  required_variables,
  object_name = "raw input dataset"
)

# ------------------------------------------------------------------------------
# Visit dictionary
# ------------------------------------------------------------------------------
visit_dictionary <- tribble(
  ~EVENT_ID, ~year_from_baseline,
  "BL", 0,
  "V04", 1,
  "V06", 2,
  "V08", 3,
  "V10", 4,
  "V12", 5,
  "V13", 6,
  "V14", 7
)

# ------------------------------------------------------------------------------
# Calculate visit-wise inflammatory indices
# ------------------------------------------------------------------------------
data_with_ratios <- raw_data %>%
  mutate(
    PATNO = as.character(PATNO),
    EVENT_ID = as.character(EVENT_ID),
    GROUP = case_when(
      PRIMDIAG == 1 ~ "PD",
      PRIMDIAG == 17 ~ "HC",
      TRUE ~ NA_character_
    ),
    GROUP = factor(GROUP, levels = c("HC", "PD")),
    NLR = safe_ratio(Neutrophils, Lymphocytes),
    MLR = safe_ratio(Monocytes, Lymphocytes),
    PLR = safe_ratio(Platelets, Lymphocytes)
  ) %>%
  left_join(visit_dictionary, by = "EVENT_ID")

message("Visit-wise NLR, MLR, and PLR calculated.")

# ------------------------------------------------------------------------------
# Create one baseline record per participant
# ------------------------------------------------------------------------------
duplicate_baseline_records <- data_with_ratios %>%
  filter(EVENT_ID == "BL") %>%
  count(PATNO, name = "n_baseline_records") %>%
  filter(n_baseline_records > 1)

if (nrow(duplicate_baseline_records) > 0) {
  warning(
    paste0(
      nrow(duplicate_baseline_records),
      " participant(s) had duplicate baseline records. ",
      "The first baseline record per participant was retained. ",
      "See Duplicate_baseline_records.csv."
    ),
    call. = FALSE
  )

  write_csv_output(
    duplicate_baseline_records,
    "Duplicate_baseline_records.csv"
  )
}

baseline_ratios <- data_with_ratios %>%
  filter(EVENT_ID == "BL") %>%
  arrange(PATNO) %>%
  distinct(PATNO, .keep_all = TRUE) %>%
  select(
    PATNO,
    baseline_NLR = NLR,
    baseline_MLR = MLR,
    baseline_PLR = PLR,
    baseline_Neutrophils = Neutrophils,
    baseline_Lymphocytes = Lymphocytes,
    baseline_Monocytes = Monocytes,
    baseline_Platelets = Platelets
  )

data_with_ratios <- data_with_ratios %>%
  left_join(
    baseline_ratios,
    by = "PATNO",
    relationship = "many-to-one"
  )

message("Baseline inflammatory measures merged into all visit rows.")

# ------------------------------------------------------------------------------
# Standardized and log-transformed NLR variables
#
# NLR_z:
#   Standardized across all available visit-level observations.
#
# baseline_NLR_z:
#   Baseline NLR standardized across participants with an available baseline
#   NLR value, then propagated across all visits for each participant.
# ------------------------------------------------------------------------------
visit_nlr_mean <- mean(data_with_ratios$NLR, na.rm = TRUE)
visit_nlr_sd <- sd(data_with_ratios$NLR, na.rm = TRUE)

baseline_nlr_mean <- mean(
  baseline_ratios$baseline_NLR,
  na.rm = TRUE
)
baseline_nlr_sd <- sd(
  baseline_ratios$baseline_NLR,
  na.rm = TRUE
)

data_with_ratios <- data_with_ratios %>%
  mutate(
    NLR_z = safe_standardize(
      NLR,
      reference_mean = visit_nlr_mean,
      reference_sd = visit_nlr_sd
    ),
    baseline_NLR_z = safe_standardize(
      baseline_NLR,
      reference_mean = baseline_nlr_mean,
      reference_sd = baseline_nlr_sd
    ),
    log_NLR = safe_log(NLR),
    log_baseline_NLR = safe_log(baseline_NLR)
  )

message("Standardized and log-transformed NLR variables created.")

# ------------------------------------------------------------------------------
# Reorder key variables
# ------------------------------------------------------------------------------
key_variables <- c(
  "PATNO",
  "EVENT_ID",
  "year_from_baseline",
  "PRIMDIAG",
  "GROUP",
  "NLR",
  "NLR_z",
  "log_NLR",
  "baseline_NLR",
  "baseline_NLR_z",
  "log_baseline_NLR",
  "MLR",
  "PLR",
  "baseline_MLR",
  "baseline_PLR",
  "Neutrophils",
  "Lymphocytes",
  "Monocytes",
  "Platelets",
  "baseline_Neutrophils",
  "baseline_Lymphocytes",
  "baseline_Monocytes",
  "baseline_Platelets"
)

data_with_ratios <- data_with_ratios %>%
  select(
    all_of(intersect(key_variables, names(data_with_ratios))),
    everything()
  )

# ------------------------------------------------------------------------------
# Analysis datasets
# ------------------------------------------------------------------------------
analysis_data <- data_with_ratios %>%
  filter(
    PRIMDIAG %in% c(1, 17),
    EVENT_ID %in% visit_dictionary$EVENT_ID
  )

pd_data <- analysis_data %>%
  filter(GROUP == "PD")

hc_data <- analysis_data %>%
  filter(GROUP == "HC")

baseline_data <- analysis_data %>%
  filter(EVENT_ID == "BL") %>%
  distinct(PATNO, .keep_all = TRUE)

pd_baseline_data <- baseline_data %>%
  filter(GROUP == "PD")

hc_baseline_data <- baseline_data %>%
  filter(GROUP == "HC")

# ------------------------------------------------------------------------------
# Descriptive summaries
# ------------------------------------------------------------------------------
overall_nlr_summary <- safe_descriptive_summary(
  analysis_data,
  "NLR"
) %>%
  mutate(
    n_participants = n_distinct(analysis_data$PATNO),
    .before = 1
  ) %>%
  rename(
    n_NLR_available = n_available,
    NLR_mean = mean,
    NLR_sd = sd,
    NLR_median = median,
    NLR_IQR = IQR,
    NLR_min = min,
    NLR_max = max
  )

baseline_nlr_by_group <- baseline_data %>%
  group_by(GROUP) %>%
  group_modify(
    ~ {
      summary_row <- safe_descriptive_summary(.x, "NLR")

      summary_row %>%
        mutate(
          n_participants = n_distinct(.x$PATNO),
          .before = 1
        )
    }
  ) %>%
  ungroup() %>%
  rename(
    n_NLR_available = n_available,
    NLR_mean = mean,
    NLR_sd = sd,
    NLR_median = median,
    NLR_IQR = IQR,
    NLR_min = min,
    NLR_max = max
  )

visitwise_nlr_by_group <- analysis_data %>%
  group_by(EVENT_ID, year_from_baseline, GROUP) %>%
  group_modify(
    ~ {
      summary_row <- safe_descriptive_summary(.x, "NLR")

      summary_row %>%
        mutate(
          n_participants = n_distinct(.x$PATNO),
          .before = 1
        )
    }
  ) %>%
  ungroup() %>%
  arrange(year_from_baseline, GROUP) %>%
  rename(
    n_NLR_available = n_available,
    NLR_mean = mean,
    NLR_sd = sd,
    NLR_median = median,
    NLR_IQR = IQR,
    NLR_min = min,
    NLR_max = max
  )

baseline_welch_nlr <- safe_welch(
  baseline_data,
  value_variable = "NLR",
  group_variable = "GROUP"
) %>%
  mutate(
    p_formatted = vapply(p, fmt_p, character(1))
  )

# ------------------------------------------------------------------------------
# Missingness
# ------------------------------------------------------------------------------
nlr_related_variables <- c(
  "Neutrophils",
  "Lymphocytes",
  "Monocytes",
  "Platelets",
  "NLR",
  "MLR",
  "PLR",
  "baseline_NLR",
  "baseline_NLR_z"
)

missingness_all_visits <- create_missingness_table(
  analysis_data,
  nlr_related_variables
)

missingness_baseline <- create_missingness_table(
  baseline_data,
  nlr_related_variables
)

# ------------------------------------------------------------------------------
# Save CSV outputs
# ------------------------------------------------------------------------------
write_csv_output(
  data_with_ratios,
  "PPMI_all_visits_with_NLR.csv"
)

write_csv_output(
  analysis_data,
  "PPMI_PD_HC_analysis_dataset_with_NLR.csv"
)

write_csv_output(
  pd_data,
  "PPMI_PD_only_all_visits_with_NLR.csv"
)

write_csv_output(
  hc_data,
  "PPMI_HC_only_all_visits_with_NLR.csv"
)

write_csv_output(
  baseline_data,
  "PPMI_baseline_PD_HC_with_NLR.csv"
)

write_csv_output(
  pd_baseline_data,
  "PPMI_baseline_PD_only_with_NLR.csv"
)

write_csv_output(
  hc_baseline_data,
  "PPMI_baseline_HC_only_with_NLR.csv"
)

write_csv_output(
  overall_nlr_summary,
  "NLR_overall_summary.csv"
)

write_csv_output(
  baseline_nlr_by_group,
  "Baseline_NLR_by_group_summary.csv"
)

write_csv_output(
  visitwise_nlr_by_group,
  "Visitwise_NLR_by_group_summary.csv"
)

write_csv_output(
  baseline_welch_nlr,
  "Baseline_NLR_Welch_test_PD_vs_HC.csv"
)

write_csv_output(
  missingness_all_visits,
  "NLR_missingness_all_visits.csv"
)

write_csv_output(
  missingness_baseline,
  "NLR_missingness_baseline.csv"
)

# ------------------------------------------------------------------------------
# Save derived Excel workbook
# ------------------------------------------------------------------------------
workbook <- createWorkbook()

workbook_sheets <- list(
  All_data_with_NLR = data_with_ratios,
  PD_HC_analysis_dataset = analysis_data,
  PD_only_all_visits = pd_data,
  HC_only_all_visits = hc_data,
  Baseline_PD_HC = baseline_data,
  Baseline_NLR_by_group = baseline_nlr_by_group,
  Visitwise_NLR_by_group = visitwise_nlr_by_group,
  Baseline_Welch_NLR = baseline_welch_nlr,
  Missingness_all_visits = missingness_all_visits,
  Missingness_baseline = missingness_baseline
)

for (sheet_name in names(workbook_sheets)) {
  addWorksheet(workbook, sheet_name)
  writeData(
    workbook,
    sheet = sheet_name,
    x = workbook_sheets[[sheet_name]],
    keepNA = TRUE
  )
  freezePane(
    workbook,
    sheet = sheet_name,
    firstRow = TRUE
  )
  setColWidths(
    workbook,
    sheet = sheet_name,
    cols = 1:ncol(workbook_sheets[[sheet_name]]),
    widths = "auto"
  )
}

saveWorkbook(
  workbook,
  derived_file,
  overwrite = TRUE
)

message("Derived Excel workbook saved: ", derived_file)

# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------
baseline_nlr_plot <- ggplot(
  baseline_data,
  aes(x = GROUP, y = NLR, color = GROUP)
) +
  geom_jitter(
    width = 0.15,
    alpha = 0.65,
    size = 2,
    na.rm = TRUE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3,
    color = "black",
    na.rm = TRUE
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    color = "black",
    na.rm = TRUE
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Baseline NLR by diagnostic group",
    x = "Diagnostic group",
    y = "Neutrophil-to-lymphocyte ratio"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(
    output_dir,
    "Baseline_NLR_by_group.png"
  ),
  plot = baseline_nlr_plot,
  width = 6,
  height = 5,
  dpi = 600
)

visit_summary_for_plot <- visitwise_nlr_by_group %>%
  filter(
    !is.na(year_from_baseline),
    !is.na(NLR_mean)
  )

visitwise_nlr_plot <- ggplot(
  visit_summary_for_plot,
  aes(
    x = year_from_baseline,
    y = NLR_mean,
    group = GROUP,
    color = GROUP
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(
      ymin = NLR_mean - NLR_sd / sqrt(n_NLR_available),
      ymax = NLR_mean + NLR_sd / sqrt(n_NLR_available)
    ),
    width = 0.12,
    na.rm = TRUE
  ) +
  scale_x_continuous(breaks = 0:7) +
  theme_classic(base_size = 13) +
  labs(
    title = "Visit-wise NLR by diagnostic group",
    x = "Years from baseline",
    y = "Mean neutrophil-to-lymphocyte ratio",
    color = "Diagnostic group"
  ) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(
    output_dir,
    "Visitwise_NLR_by_group.png"
  ),
  plot = visitwise_nlr_plot,
  width = 7,
  height = 5,
  dpi = 600
)

baseline_nlr_distribution_plot <- ggplot(
  baseline_data,
  aes(x = NLR, fill = GROUP)
) +
  geom_histogram(
    alpha = 0.55,
    bins = 30,
    position = "identity",
    na.rm = TRUE
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Distribution of baseline NLR",
    x = "Baseline neutrophil-to-lymphocyte ratio",
    y = "Count",
    fill = "Diagnostic group"
  ) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(
    output_dir,
    "Baseline_NLR_distribution.png"
  ),
  plot = baseline_nlr_distribution_plot,
  width = 7,
  height = 5,
  dpi = 600
)

# ------------------------------------------------------------------------------
# Summary text
# ------------------------------------------------------------------------------
summary_file <- file.path(output_dir, "NLR_calculation_summary.txt")

summary_connection <- file(summary_file, open = "wt")
sink(summary_connection)

cat("NLR calculation and data-preparation summary\n")
cat("============================================\n\n")

cat("Input file:\n")
cat(normalizePath(input_file, winslash = "/"), "\n\n")

cat("Derived workbook:\n")
cat(
  normalizePath(
    derived_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)

cat("Raw dataset dimensions:\n")
cat("Rows:", nrow(raw_data), "\n")
cat("Columns:", ncol(raw_data), "\n\n")

cat("Dataset dimensions after ratio calculation:\n")
cat("Rows:", nrow(data_with_ratios), "\n")
cat("Columns:", ncol(data_with_ratios), "\n\n")

cat("PD/HC analysis dataset:\n")
cat("Rows:", nrow(analysis_data), "\n")
cat(
  "Unique participants:",
  n_distinct(analysis_data$PATNO),
  "\n\n"
)

cat("Baseline PD/HC dataset:\n")
cat("Rows:", nrow(baseline_data), "\n")
cat(
  "Unique participants:",
  n_distinct(baseline_data$PATNO),
  "\n\n"
)

cat("NLR standardization parameters:\n")
cat("Visit-wise NLR mean:", visit_nlr_mean, "\n")
cat("Visit-wise NLR SD:", visit_nlr_sd, "\n")
cat("Baseline NLR mean:", baseline_nlr_mean, "\n")
cat("Baseline NLR SD:", baseline_nlr_sd, "\n\n")

cat("Visit distribution:\n")
print(table(analysis_data$EVENT_ID, useNA = "ifany"))

cat("\nGroup distribution across all analysis rows:\n")
print(table(analysis_data$GROUP, useNA = "ifany"))

cat("\nBaseline group distribution:\n")
print(table(baseline_data$GROUP, useNA = "ifany"))

cat("\nOverall NLR summary:\n")
print(overall_nlr_summary)

cat("\nBaseline NLR by group:\n")
print(baseline_nlr_by_group)

cat("\nBaseline Welch test, PD versus HC:\n")
print(baseline_welch_nlr)

cat("\nVisit-wise NLR by group:\n")
print(visitwise_nlr_by_group)

cat("\nNLR-related missingness across all visits:\n")
print(missingness_all_visits)

cat("\nNLR-related missingness at baseline:\n")
print(missingness_baseline)

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
message("NLR calculation and dataset preparation completed successfully.")
message(
  "Main derived workbook: ",
  normalizePath(derived_file, winslash = "/", mustWork = FALSE)
)
message(
  "Supporting outputs: ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)
