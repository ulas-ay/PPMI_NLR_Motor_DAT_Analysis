# =========================================================
# 03_longitudinal_nlr_analysis.R
#
# Longitudinal analysis of neutrophil-to-lymphocyte ratio
# in Parkinson's disease (PD) and healthy controls (HC).
#
# Primary model:
#   NLR ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO)
#
# This script:
#   - uses the participant-visit-level dataset created by Script 01;
#   - audits PATNO-EVENT_ID duplicates before modelling;
#   - collapses residual duplicates to one row per participant-visit;
#   - reports outcome-inclusive and predictor/covariate-only missingness;
#   - reports participant-level missingness;
#   - fits the primary longitudinal mixed-effects model;
#   - estimates marginal means, group contrasts, and longitudinal slopes;
#   - performs model diagnostics; and
#   - fits a log-transformed NLR sensitivity model.
#
# Default input:
#   <PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
#
# Output:
#   <PPMI_OUTPUT_ROOT>/03_NLR_LONGITUDINAL
#
# Paths can be supplied using command-line arguments or environment
# variables. See the "Paths" section below.
# =========================================================

# ---------------------------------------------------------
# Packages
# ---------------------------------------------------------
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
    "performance",
    "openxlsx"
)

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
    stop(
        paste0(
            "Missing required R package(s): ",
            paste(missing_packages, collapse = ", "),
            ".\nInstall them before running this script."
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
library(openxlsx)

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

output_root <- if (length(args) >= 2 && nzchar(args[2])) {
    args[2]
} else {
    Sys.getenv("PPMI_OUTPUT_ROOT", unset = file.path(getwd(), "outputs"))
}

file_path <- if (length(args) >= 1 && nzchar(args[1])) {
    args[1]
} else {
    env_input <- Sys.getenv("PPMI_NLR_DATA_FILE", unset = "")
    if (nzchar(env_input)) {
        env_input
    } else {
        file.path(
            output_root,
            "01_NLR",
            "PPMI_with_NLR_all_visits_updated.xlsx"
        )
    }
}

out_dir <- file.path(output_root, "03_NLR_LONGITUDINAL")

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

if (!file.exists(file_path)) {
    stop(
        paste0(
            "Input file not found:\n",
            file_path,
            "\n\nRun 01_calculate_nlr_and_prepare_datasets.R first, ",
            "or provide the input file explicitly."
        ),
        call. = FALSE
    )
}

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
        tibble::rownames_to_column("term") %>%
        dplyr::rename(
            estimate = Estimate,
            SE = `Std. Error`,
            t_value = `t value`
        )
    
    if (!"df" %in% names(out)) {
        out$df <- NA_real_
    }
    
    if ("Pr(>|t|)" %in% names(out)) {
        out <- out %>% dplyr::rename(p_value = `Pr(>|t|)`)
    } else if ("Pr(>|z|)" %in% names(out)) {
        out <- out %>% dplyr::rename(p_value = `Pr(>|z|)`)
    } else {
        out$p_value <- NA_real_
    }
    
    out %>%
        mutate(
            lower_95_CI = estimate - 1.96 * SE,
            upper_95_CI = estimate + 1.96 * SE,
            p_formatted = fmt_p(p_value)
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
}

make_anova_table <- function(model) {
    out <- as.data.frame(car::Anova(model, type = 3, test.statistic = "Chisq")) %>%
        tibble::rownames_to_column("term")
    
    if ("Chisq" %in% names(out)) {
        out <- out %>% rename(test_statistic = Chisq)
    } else if ("F" %in% names(out)) {
        out <- out %>% rename(test_statistic = F)
    }
    
    if ("Df" %in% names(out)) {
        out <- out %>% rename(df = Df)
    }
    
    if ("Pr(>Chisq)" %in% names(out)) {
        out <- out %>% rename(p_value = `Pr(>Chisq)`)
    } else if ("Pr(>F)" %in% names(out)) {
        out <- out %>% rename(p_value = `Pr(>F)`)
    } else {
        out$p_value <- NA_real_
    }
    
    out %>%
        mutate(
            p_formatted = fmt_p(p_value)
        )
}

make_missingness_table <- function(data, vars) {
    
    vars <- vars[vars %in% names(data)]
    
    if (length(vars) == 0) {
        return(tibble(
            variable = character(),
            n_missing = numeric(),
            pct_missing = numeric()
        ))
    }
    
    data %>%
        ungroup() %>%
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

make_row_missing_flag <- function(data, vars, flag_name) {
    
    vars <- vars[vars %in% names(data)]
    
    if (length(vars) == 0) {
        data[[flag_name]] <- FALSE
        return(data)
    }
    
    data %>%
        ungroup() %>%
        mutate(
            "{flag_name}" := if_any(all_of(vars), is.na)
        )
}

make_participant_missingness <- function(data, vars, label) {
    
    flag_data <- make_row_missing_flag(
        data = data,
        vars = vars,
        flag_name = "row_has_missing_for_this_definition"
    )
    
    participant_level <- flag_data %>%
        group_by(PATNO, GROUP) %>%
        summarise(
            n_rows = n(),
            any_missing = any(row_has_missing_for_this_definition),
            .groups = "drop"
        ) %>%
        mutate(
            missingness_definition = label
        )
    
    summary <- participant_level %>%
        summarise(
            missingness_definition = label,
            n_participants = n(),
            n_with_at_least_one_missing = sum(any_missing),
            pct_with_at_least_one_missing = 100 * mean(any_missing),
            .groups = "drop"
        )
    
    list(
        participant_level = participant_level,
        summary = summary
    )
}

safe_check_model_png <- function(model, filename) {
    out_file <- file.path(out_dir, filename)
    
    tryCatch(
        {
            diagnostic_check <- performance::check_model(model)
            
            png(
                filename = out_file,
                width = 2400,
                height = 1800,
                res = 300
            )
            print(diagnostic_check)
            dev.off()
        },
        error = function(e) {
            if (dev.cur() > 1) dev.off()
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

# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------
df_raw <- readxl::read_excel(file_path, sheet = "All_data_with_NLR")

cat("\nData loaded successfully.\n")
cat("Rows:", nrow(df_raw), "\n")
cat("Columns:", ncol(df_raw), "\n\n")

# ---------------------------------------------------------
# Required variables
# ---------------------------------------------------------
required_vars <- c(
    "PATNO",
    "EVENT_ID",
    "PRIMDIAG",
    "GROUP",
    "NLR",
    "year_from_baseline"
)

missing_vars <- setdiff(required_vars, names(df_raw))

if (length(missing_vars) > 0) {
    stop(
        paste0(
            "The following required column(s) were not found in the dataset:\n",
            paste(missing_vars, collapse = ", "),
            "\n\nInspect names(df_raw) and update the input data or script if needed."
        ),
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Detect covariate columns
# ---------------------------------------------------------
age_col <- find_col(
    df_raw,
    c("age", "AGE", "Age", "age_bl", "AGE_AT_VISIT", "AGE_AT_BASELINE", "age_at_baseline")
)

sex_col <- find_col(
    df_raw,
    c("sex", "SEX", "Sex", "gender", "GENDER")
)

bmi_col <- find_col(
    df_raw,
    c("BMI", "bmi", "Bmi", "body_mass_index", "BodyMassIndex")
)

site_col <- find_col(
    df_raw,
    c("SITE", "site", "Site", "siteid", "SITEID", "site_id")
)

detected_columns <- tibble(
    variable = c("age", "sex", "bmi", "SITE"),
    detected_column = c(age_col, sex_col, bmi_col, site_col)
)

write.csv(
    detected_columns,
    file.path(out_dir, "Detected_input_columns.csv"),
    row.names = FALSE
)

cat("Detected columns:\n")
print(detected_columns)
cat("\n")

if (any(is.na(c(age_col, sex_col, bmi_col, site_col)))) {
    stop(
        "Could not detect one or more required columns: age, sex, BMI, SITE. Check names(df_raw) and update candidate names.",
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Prepare raw analysis dataset
# ---------------------------------------------------------
visit_order <- c("BL", "V04", "V06", "V08", "V10", "V12", "V13", "V14")

dat_raw <- df_raw %>%
    filter(
        PRIMDIAG %in% c(1, 17),
        EVENT_ID %in% visit_order
    ) %>%
    mutate(
        PATNO = as.character(PATNO),
        EVENT_ID = as.character(EVENT_ID),
        GROUP = factor(as.character(GROUP), levels = c("HC", "PD"))
    )

# ---------------------------------------------------------
# Duplicate audit before modelling
# ---------------------------------------------------------
duplicate_keys <- dat_raw %>%
    count(PATNO, EVENT_ID, name = "n_rows") %>%
    filter(n_rows > 1) %>%
    arrange(desc(n_rows), PATNO, EVENT_ID)

duplicate_rows <- dat_raw %>%
    semi_join(duplicate_keys, by = c("PATNO", "EVENT_ID")) %>%
    arrange(PATNO, EVENT_ID)

if (nrow(duplicate_rows) > 0) {
    
    duplicate_summary_by_group_visit <- duplicate_rows %>%
        count(GROUP, EVENT_ID, PATNO, name = "n_rows") %>%
        group_by(GROUP, EVENT_ID) %>%
        summarise(
            n_duplicate_participant_visits = n(),
            n_rows_in_duplicate_participant_visits = sum(n_rows),
            n_extra_rows_due_to_duplicates = sum(n_rows - 1),
            max_rows_per_participant_visit = max(n_rows),
            .groups = "drop"
        ) %>%
        arrange(GROUP, EVENT_ID)
    
    variables_to_check <- setdiff(names(dat_raw), c("PATNO", "EVENT_ID"))
    
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
    rows_before_cleaning = nrow(dat_raw),
    unique_participants_before_cleaning = n_distinct(dat_raw$PATNO),
    unique_participant_visits_before_cleaning =
        n_distinct(paste(dat_raw$PATNO, dat_raw$EVENT_ID)),
    n_duplicate_participant_visit_keys = nrow(duplicate_keys),
    n_rows_in_duplicate_participant_visits = nrow(duplicate_rows),
    n_extra_rows_due_to_duplicates =
        nrow(dat_raw) - n_distinct(paste(dat_raw$PATNO, dat_raw$EVENT_ID))
)

write.csv(
    duplicate_keys,
    file.path(out_dir, "Duplicate_PATNO_EVENT_keys_before_model_cleaning.csv"),
    row.names = FALSE
)

write.csv(
    duplicate_rows,
    file.path(out_dir, "Duplicate_PATNO_EVENT_rows_before_model_cleaning.csv"),
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

# ---------------------------------------------------------
# Collapse to participant-visit level
# ---------------------------------------------------------
dat_visit <- collapse_to_participant_visit(dat_raw)

dat_visit <- dat_visit %>%
    mutate(
        PATNO = as.character(PATNO),
        EVENT_ID = as.character(EVENT_ID),
        GROUP = factor(as.character(GROUP), levels = c("HC", "PD"))
    )

cleaning_summary <- tibble(
    rows_before_cleaning = nrow(dat_raw),
    rows_after_participant_visit_cleaning = nrow(dat_visit),
    rows_removed_by_cleaning = nrow(dat_raw) - nrow(dat_visit),
    unique_participants_before_cleaning = n_distinct(dat_raw$PATNO),
    unique_participants_after_cleaning = n_distinct(dat_visit$PATNO),
    unique_participant_visits_before_cleaning =
        n_distinct(paste(dat_raw$PATNO, dat_raw$EVENT_ID)),
    unique_participant_visits_after_cleaning =
        n_distinct(paste(dat_visit$PATNO, dat_visit$EVENT_ID)),
    duplicate_PATNO_EVENT_after_cleaning =
        nrow(dat_visit) - n_distinct(paste(dat_visit$PATNO, dat_visit$EVENT_ID))
)

write.csv(
    cleaning_summary,
    file.path(out_dir, "Participant_visit_cleaning_summary.csv"),
    row.names = FALSE
)

if (cleaning_summary$duplicate_PATNO_EVENT_after_cleaning != 0) {
    stop(
        "Duplicate PATNO-EVENT_ID rows remain after participant-visit cleaning.",
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Prepare analysis dataset
# IMPORTANT:
# We explicitly keep only modelling columns here.
# This prevents duplicate column-name errors.
# ---------------------------------------------------------
dat0 <- dat_visit %>%
    mutate(
        PATNO_model = factor(PATNO),
        SITE_model = factor(.data[[site_col]]),
        GROUP_model = factor(as.character(GROUP), levels = c("HC", "PD")),
        EVENT_ID_model = factor(EVENT_ID, levels = visit_order),
        
        year_model = year_from_baseline,
        year_c_model = as.numeric(year_from_baseline),
        
        NLR_model = as.numeric(NLR),
        log_NLR_model = ifelse(!is.na(NLR_model) & NLR_model > 0, log(NLR_model), NA_real_),
        
        age_model = as.numeric(.data[[age_col]]),
        sex_model = factor(.data[[sex_col]]),
        bmi_model = as.numeric(.data[[bmi_col]])
    ) %>%
    select(
        PATNO = PATNO_model,
        SITE = SITE_model,
        GROUP = GROUP_model,
        PRIMDIAG,
        EVENT_ID = EVENT_ID_model,
        year = year_model,
        year_c = year_c_model,
        NLR = NLR_model,
        log_NLR_model,
        age = age_model,
        sex = sex_model,
        bmi = bmi_model
    )

dat0 <- droplevels(dat0)

# ---------------------------------------------------------
# Missingness before complete-case filtering
# ---------------------------------------------------------
model_vars <- c(
    "NLR",
    "GROUP",
    "year_c",
    "age",
    "sex",
    "bmi",
    "SITE",
    "PATNO"
)

predictor_covariate_vars <- c(
    "GROUP",
    "year_c",
    "age",
    "sex",
    "bmi",
    "SITE",
    "PATNO"
)

missingness_by_model_variable <- make_missingness_table(dat0, model_vars)

missingness_by_predictor_covariate <- make_missingness_table(
    dat0,
    predictor_covariate_vars
)

# Create row-level missingness flags BEFORE grouping.
dat0_missing_flags <- dat0 %>%
    ungroup() %>%
    mutate(
        row_has_missing_model_variable = if_any(all_of(model_vars), is.na),
        row_has_missing_predictor_covariate = if_any(all_of(predictor_covariate_vars), is.na)
    )

missingness_by_visit <- dat0_missing_flags %>%
    group_by(GROUP, EVENT_ID, year_c) %>%
    summarise(
        n_rows = n(),
        n_subjects = n_distinct(PATNO),
        
        n_NLR_missing = sum(is.na(NLR)),
        pct_NLR_missing = 100 * mean(is.na(NLR)),
        
        n_model_variable_missing_rows = sum(row_has_missing_model_variable),
        pct_model_variable_missing_rows = 100 * mean(row_has_missing_model_variable),
        
        n_predictor_or_covariate_missing_rows = sum(row_has_missing_predictor_covariate),
        pct_predictor_or_covariate_missing_rows = 100 * mean(row_has_missing_predictor_covariate),
        
        .groups = "drop"
    ) %>%
    arrange(GROUP, year_c)

participant_missing_model_vars <- make_participant_missingness(
    data = dat0,
    vars = model_vars,
    label = "Model variables including outcome NLR"
)

participant_missing_predictors <- make_participant_missingness(
    data = dat0,
    vars = predictor_covariate_vars,
    label = "Predictors/covariates only; outcome NLR excluded"
)

participant_level_missingness <- bind_rows(
    participant_missing_model_vars$participant_level,
    participant_missing_predictors$participant_level
)

participant_level_missingness_summary <- bind_rows(
    participant_missing_model_vars$summary,
    participant_missing_predictors$summary
)

# ---------------------------------------------------------
# Complete-case dataset for the main model
# ---------------------------------------------------------
dat <- dat0 %>%
    drop_na(all_of(model_vars))

dat <- droplevels(dat)

cat("Complete-case model dataset:\n")
cat("Rows:", nrow(dat), "\n")
cat("Unique participants:", n_distinct(dat$PATNO), "\n")
cat("Unique sites:", n_distinct(dat$SITE), "\n\n")

if (nrow(dat) < 10) {
    stop("Insufficient complete observations for model fitting.", call. = FALSE)
}

# ---------------------------------------------------------
# Save analysis datasets and missingness outputs
# ---------------------------------------------------------
write.csv(
    dat0,
    file.path(out_dir, "NLR_longitudinal_dataset_before_complete_case_filter.csv"),
    row.names = FALSE
)

write.csv(
    dat,
    file.path(out_dir, "NLR_longitudinal_complete_case_dataset.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_model_variable,
    file.path(out_dir, "Missingness_by_model_variable_including_outcome.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_predictor_covariate,
    file.path(out_dir, "Missingness_by_predictor_covariate_only.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_visit,
    file.path(out_dir, "Missingness_NLR_by_visit_and_group.csv"),
    row.names = FALSE
)

write.csv(
    participant_level_missingness,
    file.path(out_dir, "Participant_level_missingness_long.csv"),
    row.names = FALSE
)

write.csv(
    participant_level_missingness_summary,
    file.path(out_dir, "Participant_level_missingness_summary.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Descriptive summaries
# ---------------------------------------------------------
desc_nlr <- dat %>%
    group_by(GROUP, EVENT_ID, year_c) %>%
    summarise(
        n_rows = n(),
        n_subjects = n_distinct(PATNO),
        mean_NLR = mean(NLR, na.rm = TRUE),
        sd_NLR = sd(NLR, na.rm = TRUE),
        median_NLR = median(NLR, na.rm = TRUE),
        IQR_NLR = IQR(NLR, na.rm = TRUE),
        min_NLR = min(NLR, na.rm = TRUE),
        max_NLR = max(NLR, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(GROUP, year_c)

write.csv(
    desc_nlr,
    file.path(out_dir, "NLR_longitudinal_descriptives.csv"),
    row.names = FALSE
)

site_counts <- dat %>%
    group_by(SITE) %>%
    summarise(
        n_rows = n(),
        n_subjects = n_distinct(PATNO),
        .groups = "drop"
    ) %>%
    arrange(desc(n_rows))

write.csv(
    site_counts,
    file.path(out_dir, "SITE_distribution.csv"),
    row.names = FALSE
)

visit_counts <- dat %>%
    group_by(EVENT_ID, year_c, GROUP) %>%
    summarise(
        n_rows = n(),
        n_subjects = n_distinct(PATNO),
        .groups = "drop"
    ) %>%
    arrange(year_c, GROUP)

write.csv(
    visit_counts,
    file.path(out_dir, "Visit_counts_by_group.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Mixed model
# ---------------------------------------------------------
# Sum-to-zero contrasts are used for Type III tests.
options(contrasts = c("contr.sum", "contr.poly"))

lmm <- lmer(
    NLR ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO),
    data = dat,
    REML = FALSE
)

# ---------------------------------------------------------
# Model outputs
# ---------------------------------------------------------
anova_type3 <- car::Anova(lmm, type = 3, test.statistic = "Chisq")

fixef_tab <- make_fixef_table(lmm)

write.csv(
    fixef_tab,
    file.path(out_dir, "LMM_fixed_effects_NLR.csv"),
    row.names = FALSE
)

anova_tab <- make_anova_table(lmm)

write.csv(
    anova_tab,
    file.path(out_dir, "LMM_TypeIII_ANOVA_NLR.csv"),
    row.names = FALSE
)

randef_var <- as.data.frame(VarCorr(lmm))

write.csv(
    randef_var,
    file.path(out_dir, "LMM_random_effects_variance_NLR.csv"),
    row.names = FALSE
)

model_fit_indices <- tibble(
    AIC = AIC(lmm),
    BIC = BIC(lmm),
    logLik = as.numeric(logLik(lmm)),
    deviance = deviance(lmm),
    sigma = sigma(lmm),
    n_obs = nobs(lmm),
    n_subjects = n_distinct(dat$PATNO),
    n_sites = n_distinct(dat$SITE)
)

write.csv(
    model_fit_indices,
    file.path(out_dir, "LMM_model_fit_indices_NLR.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Estimated marginal means
# ---------------------------------------------------------
time_points <- sort(unique(dat$year_c))

emm <- emmeans(
    lmm,
    ~ GROUP | year_c,
    at = list(year_c = time_points)
)

emm_df <- as.data.frame(summary(emm, infer = c(TRUE, TRUE))) %>%
    standardize_ci_names() %>%
    mutate(
        p_formatted = if ("p.value" %in% names(.)) fmt_p(p.value) else NA_character_
    )

write.csv(
    emm_df,
    file.path(out_dir, "EMMeans_NLR_by_group_and_time.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Group contrasts at each time point
# ---------------------------------------------------------
group_contrasts_each_time <- as.data.frame(
    contrast(emm, method = "pairwise", adjust = "bonferroni")
) %>%
    mutate(
        p_formatted = fmt_p(p.value)
    )

write.csv(
    group_contrasts_each_time,
    file.path(out_dir, "Group_comparisons_NLR_at_each_timepoint.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Slopes by group
# ---------------------------------------------------------
emtr <- emtrends(lmm, ~ GROUP, var = "year_c")

emtr_df <- as.data.frame(summary(emtr, infer = c(TRUE, TRUE))) %>%
    standardize_ci_names() %>%
    mutate(
        p_formatted = if ("p.value" %in% names(.)) fmt_p(p.value) else NA_character_
    )

write.csv(
    emtr_df,
    file.path(out_dir, "Group_specific_NLR_slopes.csv"),
    row.names = FALSE
)

slope_contrast <- as.data.frame(
    pairs(emtr, adjust = "bonferroni")
) %>%
    mutate(
        p_formatted = fmt_p(p.value)
    )

write.csv(
    slope_contrast,
    file.path(out_dir, "Slope_difference_NLR_between_groups.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Baseline group difference
# ---------------------------------------------------------
emm_baseline <- emmeans(lmm, ~ GROUP, at = list(year_c = 0))

baseline_emm_df <- as.data.frame(summary(emm_baseline, infer = c(TRUE, TRUE))) %>%
    standardize_ci_names() %>%
    mutate(
        p_formatted = if ("p.value" %in% names(.)) fmt_p(p.value) else NA_character_
    )

baseline_contrast <- as.data.frame(
    pairs(emm_baseline, adjust = "bonferroni")
) %>%
    mutate(
        p_formatted = fmt_p(p.value)
    )

write.csv(
    baseline_emm_df,
    file.path(out_dir, "Baseline_EMMeans_NLR.csv"),
    row.names = FALSE
)

write.csv(
    baseline_contrast,
    file.path(out_dir, "Baseline_group_difference_NLR.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Model diagnostics
# ---------------------------------------------------------
diagnostic_df <- dat %>%
    mutate(
        fitted_value = fitted(lmm),
        residual = resid(lmm),
        pearson_residual = residual / sigma(lmm)
    )

write.csv(
    diagnostic_df,
    file.path(out_dir, "LMM_NLR_diagnostic_values.csv"),
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
    filename = file.path(out_dir, "Diagnostic_residuals_vs_fitted_NLR.png"),
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
    filename = file.path(out_dir, "Diagnostic_QQ_residuals_NLR.png"),
    plot = p_qq_resid,
    width = 7,
    height = 5,
    dpi = 600
)

p_obs_fit <- ggplot(
    diagnostic_df,
    aes(x = fitted_value, y = NLR)
) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.7) +
    theme_classic(base_size = 13) +
    labs(
        title = "Model diagnostics: observed vs fitted values",
        x = "Fitted values",
        y = "Observed NLR"
    )

ggsave(
    filename = file.path(out_dir, "Diagnostic_observed_vs_fitted_NLR.png"),
    plot = p_obs_fit,
    width = 7,
    height = 5,
    dpi = 600
)

ranef_patno <- ranef(lmm)$PATNO %>%
    tibble::rownames_to_column("PATNO") %>%
    rename(random_intercept = `(Intercept)`)

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
    filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_PATNO_NLR.png"),
    plot = p_qq_ranef_patno,
    width = 7,
    height = 5,
    dpi = 600
)

ranef_site <- ranef(lmm)$SITE %>%
    tibble::rownames_to_column("SITE") %>%
    rename(random_intercept = `(Intercept)`)

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
    filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_SITE_NLR.png"),
    plot = p_qq_ranef_site,
    width = 7,
    height = 5,
    dpi = 600
)

safe_check_model_png(lmm, "Performance_check_model_NLR.png")

# ---------------------------------------------------------
# Plotting data
# ---------------------------------------------------------
plot_df <- emm_df %>%
    mutate(
        EVENT_ID_label = case_when(
            year_c == 0 ~ "BL",
            year_c == 1 ~ "1",
            year_c == 2 ~ "2",
            year_c == 3 ~ "3",
            year_c == 4 ~ "4",
            year_c == 5 ~ "5",
            year_c == 6 ~ "6",
            year_c == 7 ~ "7",
            TRUE ~ as.character(year_c)
        )
    )

raw_df <- dat %>%
    group_by(GROUP, year_c, EVENT_ID) %>%
    summarise(
        mean_NLR = mean(NLR, na.rm = TRUE),
        se_NLR = sd(NLR, na.rm = TRUE) / sqrt(n()),
        n = n(),
        n_subjects = n_distinct(PATNO),
        .groups = "drop"
    )

p <- ggplot(
    plot_df,
    aes(x = year_c, y = emmean, group = GROUP, color = GROUP, fill = GROUP)
) +
    geom_ribbon(
        aes(ymin = lower.CL, ymax = upper.CL),
        alpha = 0.18,
        color = NA
    ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_point(
        data = raw_df,
        aes(x = year_c, y = mean_NLR, group = GROUP, color = GROUP),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        stroke = 1,
        size = 2.2
    ) +
    scale_x_continuous(
        breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
        labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
    ) +
    labs(
        title = "Longitudinal trajectory of NLR",
        subtitle = "Estimated marginal means from the mixed-effects model with 95% CI",
        x = "Years from baseline",
        y = "Neutrophil-to-lymphocyte ratio, NLR",
        color = "Group",
        fill = "Group"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

ggsave(
    filename = file.path(out_dir, "NLR_longitudinal_plot.png"),
    plot = p,
    width = 8,
    height = 5.5,
    dpi = 600
)

ggsave(
    filename = file.path(out_dir, "NLR_longitudinal_plot.pdf"),
    plot = p,
    width = 8,
    height = 5.5
)

p_spaghetti <- ggplot(
    dat,
    aes(x = year_c, y = NLR, group = PATNO, color = GROUP)
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
        breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
        labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
    ) +
    labs(
        title = "Observed longitudinal trajectories of NLR",
        subtitle = "Thin lines indicate individual participants; thick lines indicate smoothed group trends",
        x = "Years from baseline",
        y = "Neutrophil-to-lymphocyte ratio, NLR",
        color = "Group",
        fill = "Group"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold")
    )

ggsave(
    filename = file.path(out_dir, "NLR_spaghetti_plot.png"),
    plot = p_spaghetti,
    width = 8,
    height = 5.5,
    dpi = 600
)

ggsave(
    filename = file.path(out_dir, "NLR_spaghetti_plot.pdf"),
    plot = p_spaghetti,
    width = 8,
    height = 5.5
)

# ---------------------------------------------------------
# Log-NLR sensitivity model
# ---------------------------------------------------------
log_model_vars <- c(
    "log_NLR_model",
    "GROUP",
    "year_c",
    "age",
    "sex",
    "bmi",
    "SITE",
    "PATNO"
)

dat_log <- dat0 %>%
    drop_na(all_of(log_model_vars))

dat_log <- droplevels(dat_log)

lmm_log <- lmer(
    log_NLR_model ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO),
    data = dat_log,
    REML = FALSE
)

fixef_log_tab <- make_fixef_table(lmm_log)

write.csv(
    fixef_log_tab,
    file.path(out_dir, "LogNLR_LMM_fixed_effects.csv"),
    row.names = FALSE
)

anova_log_tab <- make_anova_table(lmm_log)

write.csv(
    anova_log_tab,
    file.path(out_dir, "LogNLR_LMM_TypeIII_ANOVA.csv"),
    row.names = FALSE
)

log_model_fit_indices <- tibble(
    AIC = AIC(lmm_log),
    BIC = BIC(lmm_log),
    logLik = as.numeric(logLik(lmm_log)),
    deviance = deviance(lmm_log),
    sigma = sigma(lmm_log),
    n_obs = nobs(lmm_log),
    n_subjects = n_distinct(dat_log$PATNO),
    n_sites = n_distinct(dat_log$SITE)
)

write.csv(
    log_model_fit_indices,
    file.path(out_dir, "LogNLR_LMM_model_fit_indices.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Workbook output
# ---------------------------------------------------------
wb <- openxlsx::createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
    sheet_name <- substr(sheet_name, 1, 31)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = data)
}

add_sheet(wb, "Detected_columns", detected_columns)
add_sheet(wb, "Duplicate_summary", duplicate_audit_summary)
add_sheet(wb, "Duplicate_by_visit", duplicate_summary_by_group_visit)
add_sheet(wb, "Duplicate_discrepancy", duplicate_discrepancy_by_variable)
add_sheet(wb, "Cleaning_summary", cleaning_summary)
add_sheet(wb, "Missing_model_vars", missingness_by_model_variable)
add_sheet(wb, "Missing_predictors", missingness_by_predictor_covariate)
add_sheet(wb, "Participant_missing", participant_level_missingness_summary)
add_sheet(wb, "Missing_by_visit", missingness_by_visit)
add_sheet(wb, "Complete_case_counts", model_fit_indices)
add_sheet(wb, "Descriptives", desc_nlr)
add_sheet(wb, "Visit_counts", visit_counts)
add_sheet(wb, "Site_counts", site_counts)
add_sheet(wb, "Fixed_effects", fixef_tab)
add_sheet(wb, "TypeIII_ANOVA", anova_tab)
add_sheet(wb, "Random_effects", randef_var)
add_sheet(wb, "Model_fit", model_fit_indices)
add_sheet(wb, "EMMeans", emm_df)
add_sheet(wb, "Group_contrasts", group_contrasts_each_time)
add_sheet(wb, "Slopes", emtr_df)
add_sheet(wb, "Slope_contrast", slope_contrast)
add_sheet(wb, "Baseline_EMM", baseline_emm_df)
add_sheet(wb, "Baseline_contrast", baseline_contrast)
add_sheet(wb, "Log_fixed_effects", fixef_log_tab)
add_sheet(wb, "Log_TypeIII_ANOVA", anova_log_tab)
add_sheet(wb, "Log_model_fit", log_model_fit_indices)

openxlsx::saveWorkbook(
    wb,
    file = file.path(out_dir, "NLR_longitudinal_outputs.xlsx"),
    overwrite = TRUE
)

# ---------------------------------------------------------
# Full text summary
# ---------------------------------------------------------
sink(file.path(out_dir, "LMM_NLR_summary.txt"))

cat("============================================================\n")
cat("Longitudinal analysis of NLR\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(file_path, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Model formula:\n")
cat("NLR ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO)\n\n")

cat("Time variable:\n")
cat("year_c is coded in years from baseline; baseline = 0.\n\n")

cat("Detected covariate columns:\n")
print(detected_columns, row.names = FALSE)
cat("\n\n")

cat("Duplicate audit before model cleaning:\n")
print(duplicate_audit_summary, row.names = FALSE)
cat("\n\n")

cat("Participant-visit cleaning summary:\n")
print(cleaning_summary, row.names = FALSE)
cat("\n\n")

cat("Sample size before complete-case filtering:\n")
cat("Rows =", nrow(dat0), "\n")
cat("Unique subjects =", n_distinct(dat0$PATNO), "\n")
cat("Unique sites =", n_distinct(dat0$SITE), "\n")
cat("Unique participant-visits =", n_distinct(paste(dat0$PATNO, dat0$EVENT_ID)), "\n\n")

cat("Sample size in complete-case model:\n")
cat("Rows =", nrow(dat), "\n")
cat("Unique subjects =", n_distinct(dat$PATNO), "\n")
cat("Unique sites =", n_distinct(dat$SITE), "\n")
cat("Unique participant-visits =", n_distinct(paste(dat$PATNO, dat$EVENT_ID)), "\n\n")

cat("Group counts, complete-case rows:\n")
print(table(dat$GROUP))
cat("\n")

cat("Visit counts, complete-case rows:\n")
print(table(dat$EVENT_ID))
cat("\n\n")

cat("Missingness by model variable, including outcome NLR:\n")
print(missingness_by_model_variable, row.names = FALSE)
cat("\n\n")

cat("Missingness by predictors/covariates only, excluding outcome NLR:\n")
print(missingness_by_predictor_covariate, row.names = FALSE)
cat("\n\n")

cat("Missingness by visit and group:\n")
print(missingness_by_visit, row.names = FALSE)
cat("\n\n")

cat("Participant-level missingness summary:\n")
print(participant_level_missingness_summary, row.names = FALSE)
cat("\n\n")

cat("NLR descriptives by group and visit, complete-case data:\n")
print(desc_nlr, row.names = FALSE)
cat("\n\n")

cat("Random-effects variance:\n")
print(VarCorr(lmm), comp = c("Variance", "Std.Dev."))
cat("\n\n")

cat("Model fit indices:\n")
print(model_fit_indices, row.names = FALSE)
cat("\n\n")

cat("Fixed effects summary:\n")
print(fixef_tab, row.names = FALSE)
cat("\n\n")

cat("Type III ANOVA:\n")
print(anova_tab, row.names = FALSE)
cat("\n\n")

cat("Estimated marginal means by group and time:\n")
print(emm_df, row.names = FALSE)
cat("\n\n")

cat("Group comparisons at each timepoint:\n")
print(group_contrasts_each_time, row.names = FALSE)
cat("\n\n")

cat("Group-specific NLR slopes:\n")
print(emtr_df, row.names = FALSE)
cat("\n\n")

cat("Slope difference between groups:\n")
print(slope_contrast, row.names = FALSE)
cat("\n\n")

cat("Baseline group difference:\n")
print(baseline_contrast, row.names = FALSE)
cat("\n\n")

cat("Log-NLR model fixed effects:\n")
print(fixef_log_tab, row.names = FALSE)
cat("\n\n")

cat("Log-NLR Type III ANOVA:\n")
print(anova_log_tab, row.names = FALSE)
cat("\n\n")

cat("Log-NLR model fit indices:\n")
print(log_model_fit_indices, row.names = FALSE)
cat("\n\n")

cat("Model assumption checks saved as diagnostic plots:\n")
cat("- Diagnostic_residuals_vs_fitted_NLR.png\n")
cat("- Diagnostic_QQ_residuals_NLR.png\n")
cat("- Diagnostic_observed_vs_fitted_NLR.png\n")
cat("- Diagnostic_QQ_random_intercepts_PATNO_NLR.png\n")
cat("- Diagnostic_QQ_random_intercepts_SITE_NLR.png\n")
cat("- Performance_check_model_NLR.png\n")

sink()

# ---------------------------------------------------------
# Final message
# ---------------------------------------------------------
cat("\n============================================================\n")
cat("Longitudinal NLR analysis completed.\n")
cat("All outputs were saved to:\n")
cat(out_dir, "\n\n")
cat("Main summary file:\n")
cat(file.path(out_dir, "LMM_NLR_summary.txt"), "\n\n")
cat("Excel workbook:\n")
cat(file.path(out_dir, "NLR_longitudinal_outputs.xlsx"), "\n")
cat("============================================================\n")