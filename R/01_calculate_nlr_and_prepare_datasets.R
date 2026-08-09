# =========================================================
# 01_calculate_nlr_and_prepare_datasets.R
#
# Calculate neutrophil-to-lymphocyte ratio (NLR) and prepare
# participant-visit-level datasets for downstream analyses.
#
# This script:
#   - calculates NLR, MLR, and PLR from absolute blood counts;
#   - audits duplicate PATNO-EVENT_ID records;
#   - collapses the source data to one row per participant-visit;
#   - derives participant-level baseline NLR variables;
#   - prevents row multiplication during baseline joins;
#   - creates analysis-ready PD/HC and PD-only datasets;
#   - exports descriptive, missingness, and audit summaries.
#
# Input/output paths are intentionally not hard-coded.
# Provide them either as command-line arguments:
#   Rscript R/01_calculate_nlr_and_prepare_datasets.R <input.xlsx> <output_root>
# or through environment variables:
#   PPMI_INPUT_FILE=/path/to/input.xlsx
#   PPMI_OUTPUT_ROOT=/path/to/output
#
# The individual-level PPMI data used by this script are not
# distributed with the repository.
# =========================================================

# ---------------------------------------------------------
# Packages
# ---------------------------------------------------------
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
            "Missing required R package(s): ",
            paste(missing_packages, collapse = ", "),
            ".\nInstall them before running this script, for example with:\n",
            "install.packages(c(",
            paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
            "))"
        ),
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(readxl)
    library(openxlsx)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tibble)
})

# ---------------------------------------------------------
# Input file and output directory
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

file_path <- if (length(args) >= 1 && nzchar(args[1])) {
    args[1]
} else {
    Sys.getenv("PPMI_INPUT_FILE", unset = "")
}

output_root <- if (length(args) >= 2 && nzchar(args[2])) {
    args[2]
} else {
    Sys.getenv("PPMI_OUTPUT_ROOT", unset = "")
}

if (!nzchar(file_path)) {
    stop(
        paste0(
            "No input file was specified.\n",
            "Pass the source Excel file as the first command-line argument ",
            "or set the PPMI_INPUT_FILE environment variable."
        ),
        call. = FALSE
    )
}

if (!nzchar(output_root)) {
    stop(
        paste0(
            "No output directory was specified.\n",
            "Pass the output root as the second command-line argument ",
            "or set the PPMI_OUTPUT_ROOT environment variable."
        ),
        call. = FALSE
    )
}

file_path <- normalizePath(file_path, mustWork = FALSE)
output_root <- normalizePath(output_root, mustWork = FALSE)
out_dir <- file.path(output_root, "01_NLR")

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

if (!file.exists(file_path)) {
    stop(
        paste0(
            "Input file not found:\n",
            file_path,
            "\n\nCheck the supplied input path before running this script."
        ),
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------
safe_ratio <- function(num, den) {
    out <- num / den
    out[is.infinite(out)] <- NA_real_
    out[den == 0] <- NA_real_
    out
}

fmt_p <- function(p) {
    ifelse(
        is.na(p),
        NA_character_,
        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    )
}

safe_first_nonmissing <- function(x) {
    nonmiss <- x[!is.na(x)]
    if (length(nonmiss) == 0) return(x[1])
    nonmiss[1]
}

collapse_vector <- function(x) {
    nonmiss <- x[!is.na(x)]
    
    if (length(nonmiss) == 0) {
        return(x[1])
    }
    
    if (inherits(x, "Date")) {
        return(safe_first_nonmissing(x))
    }
    
    if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) {
        return(safe_first_nonmissing(x))
    }
    
    if (is.numeric(x)) {
        if (dplyr::n_distinct(nonmiss) <= 1) {
            return(nonmiss[1])
        } else {
            return(mean(nonmiss, na.rm = TRUE))
        }
    }
    
    safe_first_nonmissing(x)
}

collapse_to_participant_visit <- function(data) {
    
    if (!all(c("PATNO", "EVENT_ID") %in% names(data))) {
        stop(
            "collapse_to_participant_visit() requires PATNO and EVENT_ID columns.",
            call. = FALSE
        )
    }
    
    data %>%
        arrange(PATNO, EVENT_ID) %>%
        group_by(PATNO, EVENT_ID) %>%
        summarise(
            across(
                .cols = everything(),
                .fns = collapse_vector
            ),
            .groups = "drop"
        )
}

safe_welch <- function(data, value_var, group_var = "GROUP") {
    
    dd <- data %>%
        dplyr::select(dplyr::all_of(c(value_var, group_var))) %>%
        tidyr::drop_na() %>%
        dplyr::mutate(
            GROUP_test = factor(.data[[group_var]], levels = c("PD", "HC"))
        )
    
    if (length(unique(dd$GROUP_test)) != 2 || nrow(dd) < 3) {
        return(tibble(
            variable = value_var,
            n = nrow(dd),
            estimate_PD_minus_HC = NA_real_,
            statistic = NA_real_,
            df = NA_real_,
            p = NA_real_
        ))
    }
    
    tt <- tryCatch(
        t.test(dd[[value_var]] ~ dd$GROUP_test, var.equal = FALSE),
        error = function(e) NULL
    )
    
    if (is.null(tt)) {
        return(tibble(
            variable = value_var,
            n = nrow(dd),
            estimate_PD_minus_HC = NA_real_,
            statistic = NA_real_,
            df = NA_real_,
            p = NA_real_
        ))
    }
    
    tibble(
        variable = value_var,
        n = nrow(dd),
        estimate_PD_minus_HC = unname(diff(rev(tt$estimate))),
        statistic = unname(tt$statistic),
        df = unname(tt$parameter),
        p = tt$p.value
    )
}

make_missingness_table <- function(data, vars) {
    
    vars <- vars[vars %in% names(data)]
    
    data %>%
        summarise(
            across(
                all_of(vars),
                list(
                    n_missing = ~sum(is.na(.x)),
                    pct_missing = ~100 * mean(is.na(.x))
                ),
                .names = "{.col}_{.fn}"
            )
        ) %>%
        pivot_longer(
            cols = everything(),
            names_to = "variable_stat",
            values_to = "value"
        ) %>%
        separate(
            variable_stat,
            into = c("variable", "stat"),
            sep = "_(?=n_missing|pct_missing)",
            remove = TRUE
        ) %>%
        pivot_wider(
            names_from = stat,
            values_from = value
        ) %>%
        arrange(variable)
}

# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------
df_raw <- readxl::read_excel(file_path)

cat("\nData loaded successfully.\n")
cat("Raw rows:", nrow(df_raw), "\n")
cat("Raw columns:", ncol(df_raw), "\n\n")

# ---------------------------------------------------------
# Required column check
# ---------------------------------------------------------
required_vars <- c(
    "PATNO",
    "EVENT_ID",
    "PRIMDIAG",
    "Neutrophils",
    "Lymphocytes",
    "Monocytes",
    "Platelets"
)

missing_vars <- setdiff(required_vars, names(df_raw))

if (length(missing_vars) > 0) {
    stop(
        paste0(
            "The following required columns were not found in the input dataset:\n",
            paste(missing_vars, collapse = ", "),
            "\n\nRun names(df_raw) to inspect the available column names."
        ),
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Event-year mapping
# ---------------------------------------------------------
event_map <- tibble::tribble(
    ~EVENT_ID, ~year_from_baseline,
    "BL",  0,
    "V04", 1,
    "V06", 2,
    "V08", 3,
    "V10", 4,
    "V12", 5,
    "V13", 6,
    "V14", 7
)

# ---------------------------------------------------------
# Basic formatting and initial NLR calculation
# ---------------------------------------------------------
df_nlr_raw <- df_raw %>%
    mutate(
        EVENT_ID = as.character(EVENT_ID),
        PATNO = as.character(PATNO),
        
        GROUP = case_when(
            PRIMDIAG == 1  ~ "PD",
            PRIMDIAG == 17 ~ "HC",
            TRUE ~ NA_character_
        ),
        GROUP = factor(GROUP, levels = c("HC", "PD")),
        
        NLR = safe_ratio(Neutrophils, Lymphocytes),
        MLR = safe_ratio(Monocytes, Lymphocytes),
        PLR = safe_ratio(Platelets, Lymphocytes)
    ) %>%
    left_join(event_map, by = "EVENT_ID")

cat("NLR, MLR, and PLR were calculated for the raw rows.\n")

# ---------------------------------------------------------
# Duplicate audit before collapsing
# ---------------------------------------------------------
df_pd_hc_raw <- df_nlr_raw %>%
    filter(
        PRIMDIAG %in% c(1, 17),
        EVENT_ID %in% event_map$EVENT_ID
    )

duplicate_keys <- df_pd_hc_raw %>%
    count(PATNO, EVENT_ID, name = "n_rows") %>%
    filter(n_rows > 1) %>%
    arrange(desc(n_rows), PATNO, EVENT_ID)

duplicate_rows <- df_pd_hc_raw %>%
    semi_join(duplicate_keys, by = c("PATNO", "EVENT_ID")) %>%
    arrange(PATNO, EVENT_ID)

if (nrow(duplicate_rows) > 0) {
    
    duplicate_summary_by_group_visit <- duplicate_rows %>%
        count(GROUP, EVENT_ID, year_from_baseline, PATNO, name = "n_rows") %>%
        group_by(GROUP, EVENT_ID, year_from_baseline) %>%
        summarise(
            n_duplicate_participant_visits = n(),
            n_rows_in_duplicate_participant_visits = sum(n_rows),
            n_extra_rows_due_to_duplicates = sum(n_rows - 1),
            max_rows_per_participant_visit = max(n_rows),
            .groups = "drop"
        ) %>%
        arrange(GROUP, year_from_baseline)
    
    variables_to_check <- setdiff(names(df_pd_hc_raw), c("PATNO", "EVENT_ID"))
    
    duplicate_discrepancy_by_variable <- bind_rows(
        lapply(
            variables_to_check,
            function(v) {
                duplicate_rows %>%
                    group_by(PATNO, EVENT_ID) %>%
                    summarise(
                        n_distinct_nonmissing = n_distinct(.data[[v]], na.rm = TRUE),
                        .groups = "drop"
                    ) %>%
                    summarise(
                        variable = v,
                        n_duplicate_participant_visits_checked = n(),
                        n_with_more_than_one_distinct_nonmissing_value =
                            sum(n_distinct_nonmissing > 1, na.rm = TRUE),
                        pct_with_more_than_one_distinct_nonmissing_value =
                            100 * mean(n_distinct_nonmissing > 1, na.rm = TRUE),
                        .groups = "drop"
                    )
            }
        )
    ) %>%
        arrange(desc(n_with_more_than_one_distinct_nonmissing_value), variable)
    
} else {
    
    duplicate_summary_by_group_visit <- tibble(
        GROUP = character(),
        EVENT_ID = character(),
        year_from_baseline = numeric(),
        n_duplicate_participant_visits = integer(),
        n_rows_in_duplicate_participant_visits = integer(),
        n_extra_rows_due_to_duplicates = integer(),
        max_rows_per_participant_visit = integer()
    )
    
    duplicate_discrepancy_by_variable <- tibble(
        variable = character(),
        n_duplicate_participant_visits_checked = integer(),
        n_with_more_than_one_distinct_nonmissing_value = integer(),
        pct_with_more_than_one_distinct_nonmissing_value = numeric()
    )
}

duplicate_audit_summary <- tibble(
    raw_rows_PD_HC = nrow(df_pd_hc_raw),
    raw_unique_participants = n_distinct(df_pd_hc_raw$PATNO),
    raw_unique_participant_visits = n_distinct(paste(df_pd_hc_raw$PATNO, df_pd_hc_raw$EVENT_ID)),
    n_duplicate_participant_visit_keys = nrow(duplicate_keys),
    n_rows_in_duplicate_participant_visits = nrow(duplicate_rows),
    n_extra_rows_due_to_duplicates = nrow(df_pd_hc_raw) -
        n_distinct(paste(df_pd_hc_raw$PATNO, df_pd_hc_raw$EVENT_ID))
)

write.csv(
    duplicate_keys,
    file.path(out_dir, "Duplicate_PATNO_EVENT_keys_before_collapsing.csv"),
    row.names = FALSE
)

write.csv(
    duplicate_rows,
    file.path(out_dir, "Duplicate_PATNO_EVENT_rows_before_collapsing.csv"),
    row.names = FALSE
)

write.csv(
    duplicate_summary_by_group_visit,
    file.path(out_dir, "Duplicate_PATNO_EVENT_summary_by_group_visit.csv"),
    row.names = FALSE
)

write.csv(
    duplicate_discrepancy_by_variable,
    file.path(out_dir, "Duplicate_PATNO_EVENT_discrepancy_by_variable.csv"),
    row.names = FALSE
)

write.csv(
    duplicate_audit_summary,
    file.path(out_dir, "Duplicate_PATNO_EVENT_audit_summary.csv"),
    row.names = FALSE
)

cat("\nDuplicate audit summary:\n")
print(duplicate_audit_summary)

# ---------------------------------------------------------
# Collapse full dataset to participant-visit level
# ---------------------------------------------------------
df_nlr_visit <- collapse_to_participant_visit(df_nlr_raw)

df_nlr_visit <- df_nlr_visit %>%
    mutate(
        PATNO = as.character(PATNO),
        EVENT_ID = as.character(EVENT_ID),
        GROUP = factor(as.character(GROUP), levels = c("HC", "PD"))
    )

cleaning_summary <- tibble(
    rows_before_collapsing_all = nrow(df_nlr_raw),
    rows_after_collapsing_all = nrow(df_nlr_visit),
    rows_removed_by_collapsing_all = nrow(df_nlr_raw) - nrow(df_nlr_visit),
    unique_participants_before_all = n_distinct(df_nlr_raw$PATNO),
    unique_participants_after_all = n_distinct(df_nlr_visit$PATNO),
    unique_participant_visits_before_all = n_distinct(paste(df_nlr_raw$PATNO, df_nlr_raw$EVENT_ID)),
    unique_participant_visits_after_all = n_distinct(paste(df_nlr_visit$PATNO, df_nlr_visit$EVENT_ID))
)

write.csv(
    cleaning_summary,
    file.path(out_dir, "Participant_visit_cleaning_summary_all_rows.csv"),
    row.names = FALSE
)

cat("\nParticipant-visit cleaning summary:\n")
print(cleaning_summary)

# ---------------------------------------------------------
# Extract baseline NLR from participant-visit-level data
# IMPORTANT: one row per PATNO
# ---------------------------------------------------------
baseline_nlr <- df_nlr_visit %>%
    filter(EVENT_ID == "BL") %>%
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

baseline_nlr_audit <- baseline_nlr %>%
    summarise(
        n_baseline_rows = n(),
        n_unique_PATNO = n_distinct(PATNO),
        n_missing_baseline_NLR = sum(is.na(baseline_NLR)),
        pct_missing_baseline_NLR = 100 * mean(is.na(baseline_NLR))
    )

write.csv(
    baseline_nlr_audit,
    file.path(out_dir, "Baseline_NLR_participant_level_audit.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Merge baseline NLR into all participant-visit rows
# This join should NOT multiply rows
# ---------------------------------------------------------
rows_before_join <- nrow(df_nlr_visit)

df_nlr <- df_nlr_visit %>%
    left_join(baseline_nlr, by = "PATNO")

rows_after_join <- nrow(df_nlr)

baseline_join_audit <- tibble(
    rows_before_baseline_join = rows_before_join,
    rows_after_baseline_join = rows_after_join,
    rows_added_by_join = rows_after_join - rows_before_join,
    baseline_nlr_rows = nrow(baseline_nlr),
    baseline_nlr_unique_PATNO = n_distinct(baseline_nlr$PATNO)
)

write.csv(
    baseline_join_audit,
    file.path(out_dir, "Baseline_NLR_join_audit.csv"),
    row.names = FALSE
)

if (rows_after_join != rows_before_join) {
    stop(
        paste0(
            "Baseline NLR join multiplied rows. Rows before join: ",
            rows_before_join,
            "; rows after join: ",
            rows_after_join,
            ". Check baseline_nlr uniqueness."
        ),
        call. = FALSE
    )
}

cat("\nBaseline NLR was merged into all participant-visit rows.\n")
print(baseline_join_audit)

# ---------------------------------------------------------
# Z-standardized NLR variables
# ---------------------------------------------------------
nlr_mean_all <- mean(df_nlr$NLR, na.rm = TRUE)
nlr_sd_all   <- sd(df_nlr$NLR, na.rm = TRUE)

baseline_nlr_mean <- mean(baseline_nlr$baseline_NLR, na.rm = TRUE)
baseline_nlr_sd   <- sd(baseline_nlr$baseline_NLR, na.rm = TRUE)

df_nlr <- df_nlr %>%
    mutate(
        NLR_z = as.numeric((NLR - nlr_mean_all) / nlr_sd_all),
        baseline_NLR_z = as.numeric((baseline_NLR - baseline_nlr_mean) / baseline_nlr_sd),
        
        log_NLR = ifelse(!is.na(NLR) & NLR > 0, log(NLR), NA_real_),
        log_baseline_NLR = ifelse(!is.na(baseline_NLR) & baseline_NLR > 0, log(baseline_NLR), NA_real_)
    )

cat("Z-standardized and log-transformed NLR variables were added.\n")

standardization_audit <- tibble(
    nlr_mean_all = nlr_mean_all,
    nlr_sd_all = nlr_sd_all,
    baseline_nlr_mean = baseline_nlr_mean,
    baseline_nlr_sd = baseline_nlr_sd
)

write.csv(
    standardization_audit,
    file.path(out_dir, "NLR_standardization_audit.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Reorder columns: key variables first, then everything else
# ---------------------------------------------------------
key_cols <- c(
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

key_cols <- key_cols[key_cols %in% names(df_nlr)]

df_nlr <- df_nlr %>%
    select(
        all_of(key_cols),
        everything()
    )

# ---------------------------------------------------------
# Analysis subset: PD and HC only, mapped visits only
# ---------------------------------------------------------
analysis_df <- df_nlr %>%
    filter(
        PRIMDIAG %in% c(1, 17),
        EVENT_ID %in% event_map$EVENT_ID
    )

pd_df <- analysis_df %>%
    filter(GROUP == "PD")

hc_df <- analysis_df %>%
    filter(GROUP == "HC")

bl_df <- analysis_df %>%
    filter(EVENT_ID == "BL") %>%
    distinct(PATNO, .keep_all = TRUE)

pd_bl_df <- bl_df %>%
    filter(GROUP == "PD")

hc_bl_df <- bl_df %>%
    filter(GROUP == "HC")

analysis_dataset_audit <- tibble(
    all_rows_after_collapsing = nrow(df_nlr),
    analysis_rows_PD_HC_mapped_visits = nrow(analysis_df),
    analysis_unique_participants = n_distinct(analysis_df$PATNO),
    analysis_unique_participant_visits = n_distinct(paste(analysis_df$PATNO, analysis_df$EVENT_ID)),
    analysis_duplicate_PATNO_EVENT_after_cleaning =
        nrow(analysis_df) - n_distinct(paste(analysis_df$PATNO, analysis_df$EVENT_ID)),
    baseline_rows_PD_HC = nrow(bl_df),
    baseline_unique_participants_PD_HC = n_distinct(bl_df$PATNO),
    pd_rows_all_visits = nrow(pd_df),
    pd_unique_participants = n_distinct(pd_df$PATNO),
    hc_rows_all_visits = nrow(hc_df),
    hc_unique_participants = n_distinct(hc_df$PATNO),
    pd_baseline_rows = nrow(pd_bl_df),
    pd_baseline_unique_participants = n_distinct(pd_bl_df$PATNO),
    hc_baseline_rows = nrow(hc_bl_df),
    hc_baseline_unique_participants = n_distinct(hc_bl_df$PATNO)
)

write.csv(
    analysis_dataset_audit,
    file.path(out_dir, "Analysis_dataset_participant_visit_audit.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Descriptive summaries
# ---------------------------------------------------------
overall_nlr_summary <- analysis_df %>%
    summarise(
        n_rows = n(),
        n_participants = n_distinct(PATNO),
        n_participant_visits = n_distinct(paste(PATNO, EVENT_ID)),
        n_NLR_available = sum(!is.na(NLR)),
        NLR_mean = mean(NLR, na.rm = TRUE),
        NLR_sd = sd(NLR, na.rm = TRUE),
        NLR_median = median(NLR, na.rm = TRUE),
        NLR_IQR = IQR(NLR, na.rm = TRUE),
        NLR_min = min(NLR, na.rm = TRUE),
        NLR_max = max(NLR, na.rm = TRUE)
    )

baseline_nlr_by_group <- bl_df %>%
    group_by(GROUP) %>%
    summarise(
        n = n(),
        n_participants = n_distinct(PATNO),
        n_NLR_available = sum(!is.na(NLR)),
        NLR_mean = mean(NLR, na.rm = TRUE),
        NLR_sd = sd(NLR, na.rm = TRUE),
        NLR_median = median(NLR, na.rm = TRUE),
        NLR_IQR = IQR(NLR, na.rm = TRUE),
        NLR_min = min(NLR, na.rm = TRUE),
        NLR_max = max(NLR, na.rm = TRUE),
        .groups = "drop"
    )

visitwise_nlr_by_group <- analysis_df %>%
    group_by(EVENT_ID, year_from_baseline, GROUP) %>%
    summarise(
        n = n(),
        n_participants = n_distinct(PATNO),
        n_NLR_available = sum(!is.na(NLR)),
        NLR_mean = mean(NLR, na.rm = TRUE),
        NLR_sd = sd(NLR, na.rm = TRUE),
        NLR_median = median(NLR, na.rm = TRUE),
        NLR_IQR = IQR(NLR, na.rm = TRUE),
        NLR_min = min(NLR, na.rm = TRUE),
        NLR_max = max(NLR, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(year_from_baseline, GROUP)

baseline_welch_nlr <- safe_welch(bl_df, "NLR", "GROUP") %>%
    mutate(
        p_formatted = fmt_p(p)
    )

# ---------------------------------------------------------
# Missingness summaries for NLR-related variables
# ---------------------------------------------------------
nlr_vars <- c(
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

missingness_all <- make_missingness_table(analysis_df, nlr_vars)

missingness_baseline <- make_missingness_table(bl_df, nlr_vars)

# ---------------------------------------------------------
# Save CSV outputs
# ---------------------------------------------------------
write.csv(
    df_nlr,
    file.path(out_dir, "PPMI_all_visits_with_NLR_participant_visit_level.csv"),
    row.names = FALSE
)

write.csv(
    analysis_df,
    file.path(out_dir, "PPMI_PD_HC_analysis_dataset_with_NLR.csv"),
    row.names = FALSE
)

write.csv(
    pd_df,
    file.path(out_dir, "PPMI_PD_only_all_visits_with_NLR.csv"),
    row.names = FALSE
)

write.csv(
    bl_df,
    file.path(out_dir, "PPMI_baseline_PD_HC_with_NLR.csv"),
    row.names = FALSE
)

write.csv(
    overall_nlr_summary,
    file.path(out_dir, "NLR_overall_summary.csv"),
    row.names = FALSE
)

write.csv(
    baseline_nlr_by_group,
    file.path(out_dir, "Baseline_NLR_by_group_summary.csv"),
    row.names = FALSE
)

write.csv(
    visitwise_nlr_by_group,
    file.path(out_dir, "Visitwise_NLR_by_group_summary.csv"),
    row.names = FALSE
)

write.csv(
    baseline_welch_nlr,
    file.path(out_dir, "Baseline_NLR_Welch_test_PD_vs_HC.csv"),
    row.names = FALSE
)

write.csv(
    missingness_all,
    file.path(out_dir, "NLR_missingness_all_visits.csv"),
    row.names = FALSE
)

write.csv(
    missingness_baseline,
    file.path(out_dir, "NLR_missingness_baseline.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Save updated Excel workbook
# ---------------------------------------------------------
xlsx_out <- file.path(out_dir, "PPMI_with_NLR_all_visits_updated.xlsx")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "All_data_with_NLR")
openxlsx::writeData(wb, "All_data_with_NLR", df_nlr)

openxlsx::addWorksheet(wb, "PD_HC_analysis_dataset")
openxlsx::writeData(wb, "PD_HC_analysis_dataset", analysis_df)

openxlsx::addWorksheet(wb, "PD_only_all_visits")
openxlsx::writeData(wb, "PD_only_all_visits", pd_df)

openxlsx::addWorksheet(wb, "Baseline_PD_HC")
openxlsx::writeData(wb, "Baseline_PD_HC", bl_df)

openxlsx::addWorksheet(wb, "Baseline_NLR_by_group")
openxlsx::writeData(wb, "Baseline_NLR_by_group", baseline_nlr_by_group)

openxlsx::addWorksheet(wb, "Visitwise_NLR_by_group")
openxlsx::writeData(wb, "Visitwise_NLR_by_group", visitwise_nlr_by_group)

openxlsx::addWorksheet(wb, "Baseline_Welch_NLR")
openxlsx::writeData(wb, "Baseline_Welch_NLR", baseline_welch_nlr)

openxlsx::addWorksheet(wb, "Missingness_all_visits")
openxlsx::writeData(wb, "Missingness_all_visits", missingness_all)

openxlsx::addWorksheet(wb, "Missingness_baseline")
openxlsx::writeData(wb, "Missingness_baseline", missingness_baseline)

openxlsx::addWorksheet(wb, "Duplicate_audit_summary")
openxlsx::writeData(wb, "Duplicate_audit_summary", duplicate_audit_summary)

openxlsx::addWorksheet(wb, "Duplicate_by_visit")
openxlsx::writeData(wb, "Duplicate_by_visit", duplicate_summary_by_group_visit)

openxlsx::addWorksheet(wb, "Duplicate_discrepancy")
openxlsx::writeData(wb, "Duplicate_discrepancy", duplicate_discrepancy_by_variable)

openxlsx::addWorksheet(wb, "Cleaning_summary")
openxlsx::writeData(wb, "Cleaning_summary", cleaning_summary)

openxlsx::addWorksheet(wb, "Baseline_join_audit")
openxlsx::writeData(wb, "Baseline_join_audit", baseline_join_audit)

openxlsx::addWorksheet(wb, "Analysis_dataset_audit")
openxlsx::writeData(wb, "Analysis_dataset_audit", analysis_dataset_audit)

openxlsx::addWorksheet(wb, "Standardization_audit")
openxlsx::writeData(wb, "Standardization_audit", standardization_audit)

openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)

cat("\nUpdated Excel workbook saved to:\n")
cat(xlsx_out, "\n")

# ---------------------------------------------------------
# Plots
# ---------------------------------------------------------

# Baseline NLR by group
p_bl_nlr <- ggplot(
    bl_df,
    aes(x = GROUP, y = NLR, color = GROUP)
) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 2) +
    stat_summary(
        fun = mean,
        geom = "point",
        size = 3,
        color = "black"
    ) +
    stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.15,
        color = "black"
    ) +
    theme_classic(base_size = 13) +
    labs(
        title = "Baseline NLR by group",
        x = "Group",
        y = "Neutrophil-to-lymphocyte ratio, NLR"
    ) +
    theme(
        legend.position = "none"
    )

ggsave(
    filename = file.path(out_dir, "Baseline_NLR_by_group.png"),
    plot = p_bl_nlr,
    width = 6,
    height = 5,
    dpi = 600
)

# Visit-wise NLR trajectory by group
p_visit_nlr <- ggplot(
    analysis_df,
    aes(x = year_from_baseline, y = NLR, color = GROUP)
) +
    geom_point(alpha = 0.25, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    theme_classic(base_size = 13) +
    scale_x_continuous(
        breaks = sort(unique(analysis_df$year_from_baseline)),
        labels = sort(unique(analysis_df$year_from_baseline))
    ) +
    labs(
        title = "NLR across visits",
        x = "Years from baseline",
        y = "Neutrophil-to-lymphocyte ratio, NLR",
        color = "Group"
    )

ggsave(
    filename = file.path(out_dir, "Visitwise_NLR_by_group.png"),
    plot = p_visit_nlr,
    width = 7,
    height = 5,
    dpi = 600
)

# Baseline NLR histogram
p_hist_nlr <- ggplot(
    bl_df,
    aes(x = NLR, fill = GROUP)
) +
    geom_histogram(alpha = 0.55, bins = 30, position = "identity") +
    theme_classic(base_size = 13) +
    labs(
        title = "Distribution of baseline NLR",
        x = "Baseline NLR",
        y = "Count",
        fill = "Group"
    )

ggsave(
    filename = file.path(out_dir, "Baseline_NLR_distribution.png"),
    plot = p_hist_nlr,
    width = 7,
    height = 5,
    dpi = 600
)

# ---------------------------------------------------------
# Summary text file
# ---------------------------------------------------------
sink(file.path(out_dir, "NLR_summary.txt"))

cat("NLR calculation summary\n")
cat("=======================\n\n")

cat("Input file:\n")
cat(file_path, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Total rows in original dataset:", nrow(df_raw), "\n")
cat("Total columns in original dataset:", ncol(df_raw), "\n\n")

cat("Rows after raw NLR calculation:", nrow(df_nlr_raw), "\n")
cat("Rows after participant-visit collapsing:", nrow(df_nlr), "\n")
cat("Rows removed by collapsing:", nrow(df_nlr_raw) - nrow(df_nlr), "\n\n")

cat("Duplicate audit summary among PD/HC mapped visits:\n")
print(duplicate_audit_summary)

cat("\n\nParticipant-visit cleaning summary, all rows:\n")
print(cleaning_summary)

cat("\n\nBaseline NLR join audit:\n")
print(baseline_join_audit)

cat("\n\nAnalysis dataset audit:\n")
print(analysis_dataset_audit)

cat("\n\nStandardization audit:\n")
print(standardization_audit)

cat("\n\nPD/HC analysis rows:", nrow(analysis_df), "\n")
cat("PD/HC unique participants:", n_distinct(analysis_df$PATNO), "\n")
cat("PD/HC unique participant-visits:", n_distinct(paste(analysis_df$PATNO, analysis_df$EVENT_ID)), "\n\n")

cat("Baseline PD/HC rows:", nrow(bl_df), "\n")
cat("Baseline PD/HC unique participants:", n_distinct(bl_df$PATNO), "\n\n")

cat("Visit distribution in PD/HC analysis dataset:\n")
print(table(analysis_df$EVENT_ID, useNA = "ifany"))

cat("\n\nGroup distribution in PD/HC analysis dataset:\n")
print(table(analysis_df$GROUP, useNA = "ifany"))

cat("\n\nBaseline group distribution:\n")
print(table(bl_df$GROUP, useNA = "ifany"))

cat("\n\nOverall NLR summary:\n")
print(overall_nlr_summary)

cat("\n\nBaseline NLR by group:\n")
print(baseline_nlr_by_group)

cat("\n\nBaseline Welch test for NLR, PD vs HC:\n")
print(baseline_welch_nlr)

cat("\n\nVisitwise NLR by group:\n")
print(visitwise_nlr_by_group)

cat("\n\nNLR missingness, all visits:\n")
print(missingness_all)

cat("\n\nNLR missingness, baseline:\n")
print(missingness_baseline)

cat("\n\nDuplicate discrepancy by variable:\n")
print(duplicate_discrepancy_by_variable)

sink()

# ---------------------------------------------------------
# Final message
# ---------------------------------------------------------
cat("\n============================================================\n")
cat("NLR datasets and outputs were saved successfully.\n")
cat("Output directory:\n")
cat(out_dir, "\n\n")
cat("Main updated Excel workbook:\n")
cat(xlsx_out, "\n\n")
cat("Use the following workbook for downstream analyses:\n")
cat(xlsx_out, "\n\n")
cat("KEY AUDIT FILES:\n")
cat(file.path(out_dir, "Duplicate_PATNO_EVENT_audit_summary.csv"), "\n")
cat(file.path(out_dir, "Duplicate_PATNO_EVENT_summary_by_group_visit.csv"), "\n")
cat(file.path(out_dir, "Duplicate_PATNO_EVENT_discrepancy_by_variable.csv"), "\n")
cat(file.path(out_dir, "Baseline_NLR_join_audit.csv"), "\n")
cat("============================================================\n")