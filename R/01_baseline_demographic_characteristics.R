# =========================================================
# Table 1. Baseline Demographic and Clinical Characteristics
# PD vs HC
#
# This script:
#   1. Restricts the dataset to participants with Parkinson's disease
#      and healthy controls based on PRIMDIAG coding
#   2. Creates a baseline-only dataset with one row per participant
#   3. Generates Table 1 for demographic, clinical, and laboratory variables
#   4. Exports longitudinal visit counts by diagnostic group
#
# Required input:
#   data/PPMI_analysis_dataset.xlsx
#
# Outputs:
#   outputs/tables/Table1_Demographic_Clinical_Characteristics.html
#   outputs/tables/Table1_Demographic_Clinical_Characteristics.rtf
#   outputs/tables/Table1_Demographic_Clinical_Characteristics.csv
#   outputs/tables/Visit_counts_by_group.csv
#   outputs/tables/Visit_counts_by_group_wide.csv
#   outputs/tables/Demographic_summary_text.txt
# =========================================================


# ---------------------------------------------------------
# Load required packages
# ---------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readxl",
  "tidyr",
  "stringr",
  "gtsummary",
  "gt",
  "broom",
  "tibble"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following required packages are not installed: ",
      paste(missing_packages, collapse = ", "),
      ". Please install them before running this script."
    )
  )
}

invisible(
  lapply(required_packages, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  })
)


# ---------------------------------------------------------
# Define paths
#
# Note:
# The individual-level PPMI dataset is not included in this
# repository. After obtaining access to the PPMI data and
# preparing the analysis dataset, place the file in the data/
# folder or update the path below accordingly.
# ---------------------------------------------------------

data_path <- file.path("data", "PPMI_analysis_dataset.xlsx")
out_dir   <- file.path("outputs", "tables")

if (!file.exists(data_path)) {
  stop(
    paste0(
      "Input file not found: ", data_path, "\n",
      "Please place the analysis dataset in the data/ folder ",
      "or update data_path in this script."
    )
  )
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------

df <- readxl::read_excel(data_path)


# ---------------------------------------------------------
# Define diagnostic groups
#
# The analysis is restricted to:
#   PRIMDIAG = 1   Parkinson's disease
#   PRIMDIAG = 17  Healthy control
#
# Other diagnostic categories, if present, are not included
# in this Table 1 analysis.
# ---------------------------------------------------------

df <- df %>%
  dplyr::filter(PRIMDIAG %in% c(1, 17)) %>%
  dplyr::mutate(
    GROUP = factor(
      dplyr::case_when(
        PRIMDIAG == 1  ~ "PD",
        PRIMDIAG == 17 ~ "HC",
        TRUE ~ NA_character_
      ),
      levels = c("HC", "PD")
    )
  )


# ---------------------------------------------------------
# Create baseline-only dataset
#
# Baseline analyses use one row per participant to avoid
# inflating the sample size due to repeated longitudinal records.
# ---------------------------------------------------------

bl <- df %>%
  dplyr::filter(EVENT_ID == "BL") %>%
  dplyr::distinct(PATNO, .keep_all = TRUE)


# ---------------------------------------------------------
# Recode categorical variables
#
# SEX:
#   0 = Female
#   1 = Male
#
# DOMSIDE:
#   1 = Left
#   2 = Right
#   3 = Symmetric
# ---------------------------------------------------------

bl <- bl %>%
  dplyr::mutate(
    SEX = factor(
      SEX,
      levels = c(0, 1),
      labels = c("Female", "Male")
    ),
    DOMSIDE = factor(
      DOMSIDE,
      levels = c(1, 2, 3),
      labels = c("Left", "Right", "Symmetric")
    )
  )


# ---------------------------------------------------------
# Define variables for Table 1
# ---------------------------------------------------------

vars_cont <- c(
  "age",
  "BMI",
  "Monocytes",
  "Lymphocytes",
  "Neutrophils"
)

vars_cat <- c(
  "SEX",
  "DOMSIDE"
)

table1_vars <- c(vars_cont, vars_cat)


# ---------------------------------------------------------
# Check whether all required variables exist
# ---------------------------------------------------------

required_columns <- c(
  "PATNO",
  "EVENT_ID",
  "PRIMDIAG",
  "GROUP",
  table1_vars
)

missing_columns <- setdiff(required_columns, names(bl))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the dataset: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# ---------------------------------------------------------
# Generate Table 1
#
# Continuous variables:
#   Mean ± SD, Welch two-sample t-test
#
# Categorical variables:
#   n (%), chi-square test
# ---------------------------------------------------------

tbl1 <- bl %>%
  dplyr::select(
    GROUP,
    dplyr::all_of(vars_cont),
    dplyr::all_of(vars_cat)
  ) %>%
  gtsummary::tbl_summary(
    by = GROUP,
    statistic = list(
      gtsummary::all_continuous() ~ "{mean} ± {sd}",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = gtsummary::all_continuous() ~ 2,
    missing = "ifany",
    label = list(
      age ~ "Age, years",
      BMI ~ "Body mass index",
      Monocytes ~ "Monocytes",
      Lymphocytes ~ "Lymphocytes",
      Neutrophils ~ "Neutrophils",
      SEX ~ "Sex",
      DOMSIDE ~ "Dominant side"
    )
  ) %>%
  gtsummary::add_p(
    test = list(
      gtsummary::all_continuous() ~ "t.test",
      gtsummary::all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = ~ gtsummary::style_pvalue(.x, digits = 3)
  ) %>%
  gtsummary::add_overall(last = TRUE) %>%
  gtsummary::modify_header(
    label ~ "**Variable**",
    stat_1 ~ "**HC**",
    stat_2 ~ "**PD**",
    stat_0 ~ "**Overall**",
    p.value ~ "**p**"
  ) %>%
  gtsummary::bold_labels()


# ---------------------------------------------------------
# Export formatted Table 1
# ---------------------------------------------------------

tbl_gt <- gtsummary::as_gt(tbl1)

gt::gtsave(
  tbl_gt,
  file.path(out_dir, "Table1_Demographic_Clinical_Characteristics.html")
)

gt::gtsave(
  tbl_gt,
  file.path(out_dir, "Table1_Demographic_Clinical_Characteristics.rtf")
)


# ---------------------------------------------------------
# Export Table 1 as raw dataframe
# ---------------------------------------------------------

tbl_df <- tibble::as_tibble(tbl1, col_labels = FALSE)

write.csv(
  tbl_df,
  file.path(out_dir, "Table1_Demographic_Clinical_Characteristics.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Visit counts across longitudinal follow-up
#
# n_rows:
#   Number of rows available at each visit
#
# n_subjects:
#   Number of unique participants available at each visit
# ---------------------------------------------------------

visit_counts <- df %>%
  dplyr::group_by(EVENT_ID, GROUP) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    .groups = "drop"
  ) %>%
  dplyr::arrange(EVENT_ID, GROUP)

write.csv(
  visit_counts,
  file.path(out_dir, "Visit_counts_by_group.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Wide-format visit-count table
# ---------------------------------------------------------

visit_counts_wide <- visit_counts %>%
  dplyr::select(EVENT_ID, GROUP, n_subjects) %>%
  tidyr::pivot_wider(
    names_from = GROUP,
    values_from = n_subjects
  ) %>%
  dplyr::arrange(EVENT_ID)

write.csv(
  visit_counts_wide,
  file.path(out_dir, "Visit_counts_by_group_wide.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Text summary for manuscript
# ---------------------------------------------------------

n_hc <- bl %>%
  dplyr::filter(GROUP == "HC") %>%
  nrow()

n_pd <- bl %>%
  dplyr::filter(GROUP == "PD") %>%
  nrow()

txt <- c(
  paste0(
    "Baseline analyses included ",
    n_pd,
    " participants with Parkinson's disease and ",
    n_hc,
    " healthy controls."
  ),
  "",
  paste0(
    "Demographic and laboratory characteristics were summarized ",
    "using baseline values only to avoid inflation of sample size ",
    "due to repeated observations in the longitudinal dataset."
  )
)

writeLines(
  txt,
  file.path(out_dir, "Demographic_summary_text.txt")
)


# ---------------------------------------------------------
# Console output
# ---------------------------------------------------------

cat("\nBaseline sample sizes:\n")
cat("HC =", n_hc, "\n")
cat("PD =", n_pd, "\n")

cat("\nVisit counts by group, wide format:\n")
print(visit_counts_wide)

cat("\nAll outputs saved to:\n", out_dir, "\n")