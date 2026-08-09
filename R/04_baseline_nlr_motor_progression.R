# =========================================================
# 04_baseline_nlr_motor_progression.R
#
# Baseline neutrophil-to-lymphocyte ratio (NLR) and
# longitudinal motor progression in Parkinson's disease.
#
# Primary outcome:
#   ON-medication MDS-UPDRS Part III / UPDRS-III
#
# Primary predictor:
#   baseline_NLR_z
#   Effect estimates correspond to a 1-SD higher baseline NLR.
#
# Primary model:
#   updrs3_score_on ~ baseline_NLR_z * year_c +
#       age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
#       (1 | SITE) + (1 | PATNO)
#
# Sensitivity analyses:
#   1) Follow-up-only model additionally adjusted for baseline UPDRS-III.
#   2) Model using log-transformed baseline NLR.
#
# Predicted trajectories and simple slopes are calculated manually
# from the fixed-effect estimates and model covariance matrix; this
# script intentionally does not rely on emmeans() or emtrends().
#
# Default input:
#   <PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
#
# Output:
#   <PPMI_OUTPUT_ROOT>/04_NLR_MOTOR_PROGRESSION
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

out_dir <- file.path(output_root, "04_NLR_MOTOR_PROGRESSION")

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
    
    if (length(nonmiss) == 0) return(x[1])
    
    if (inherits(x, "Date")) return(safe_first_nonmissing(x))
    if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(safe_first_nonmissing(x))
    
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
        group_by(PATNO) %>%
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

safe_lmer <- function(formula, data, model_name) {
    if (nrow(data) < 10) {
        stop(
            paste0(model_name, " has fewer than 10 complete-case observations."),
            call. = FALSE
        )
    }
    
    lmer(
        formula = formula,
        data = data,
        REML = FALSE
    )
}

get_interaction_term <- function(term_names, term_a, term_b) {
    candidates <- c(
        paste0(term_a, ":", term_b),
        paste0(term_b, ":", term_a)
    )
    
    hit <- candidates[candidates %in% term_names]
    
    if (length(hit) == 0) {
        stop(
            paste0(
                "Could not find interaction term for ",
                term_a,
                " and ",
                term_b,
                ". Available fixed-effect terms:\n",
                paste(term_names, collapse = "\n")
            ),
            call. = FALSE
        )
    }
    
    hit[1]
}

manual_simple_slopes_baseline_nlr <- function(model, time_points) {
    
    beta <- lme4::fixef(model)
    vc <- as.matrix(vcov(model))
    
    term_main <- "baseline_NLR_z"
    term_int <- get_interaction_term(
        term_names = names(beta),
        term_a = "baseline_NLR_z",
        term_b = "year_c"
    )
    
    if (!term_main %in% names(beta)) {
        stop("baseline_NLR_z fixed-effect term not found.", call. = FALSE)
    }
    
    out <- lapply(
        time_points,
        function(t) {
            
            L <- rep(0, length(beta))
            names(L) <- names(beta)
            
            L[term_main] <- 1
            L[term_int] <- t
            
            estimate <- sum(L * beta)
            SE <- sqrt(as.numeric(t(L) %*% vc %*% L))
            
            z_value <- estimate / SE
            p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
            
            tibble(
                year_c = t,
                baseline_NLR_z_effect = estimate,
                SE = SE,
                lower_95_CI = estimate - 1.96 * SE,
                upper_95_CI = estimate + 1.96 * SE,
                z_value = z_value,
                p_value = p_value,
                p_formatted = fmt_p(p_value)
            )
        }
    )
    
    bind_rows(out)
}

manual_predictions_from_lmer <- function(model, newdata) {
    
    beta <- lme4::fixef(model)
    vc <- as.matrix(vcov(model))
    
    fixed_formula <- lme4::nobars(formula(model))
    fixed_terms <- delete.response(terms(fixed_formula))
    
    contrasts_arg <- attr(model.matrix(model), "contrasts")
    
    X <- model.matrix(
        fixed_terms,
        data = newdata,
        contrasts.arg = contrasts_arg
    )
    
    missing_cols <- setdiff(names(beta), colnames(X))
    extra_cols <- setdiff(colnames(X), names(beta))
    
    if (length(missing_cols) > 0) {
        stop(
            paste0(
                "Prediction model matrix is missing fixed-effect columns:\n",
                paste(missing_cols, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    
    if (length(extra_cols) > 0) {
        X <- X[, names(beta), drop = FALSE]
    } else {
        X <- X[, names(beta), drop = FALSE]
    }
    
    estimate <- as.numeric(X %*% beta)
    SE <- sqrt(diag(X %*% vc %*% t(X)))
    
    newdata %>%
        mutate(
            emmean = estimate,
            SE = SE,
            lower.CL = estimate - 1.96 * SE,
            upper.CL = estimate + 1.96 * SE
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

missing_vars <- setdiff(required_vars, names(df_raw))

if (length(missing_vars) > 0) {
    stop(
        paste0(
            "The following required column(s) were not found in the dataset:\n",
            paste(missing_vars, collapse = ", "),
            "\n\nRun Script 01 first or verify the supplied NLR dataset."
        ),
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Detect key columns
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

duration_col <- find_col(
    df_raw,
    c("duration_yrs", "disease_duration", "disease_duration_years", "DURATION_YRS")
)

updrs_col <- find_col(
    df_raw,
    c("updrs3_score_on", "UPDRS3_score_on", "updrs3_on", "UPDRS3_ON")
)

ledd_col <- find_col(
    df_raw,
    c(
        "LEDD", "ledd", "LEDDTOT", "LEDD_total", "ledd_total",
        "total_ledd", "LEDD_Total", "LEDD_TOT"
    )
)

domside_col <- find_col(
    df_raw,
    c("DOMSIDE", "domside", "DomSide", "dominant_side")
)

detected_columns <- tibble(
    variable = c(
        "age",
        "sex",
        "bmi",
        "SITE",
        "duration_yrs",
        "updrs3_score_on",
        "LEDD",
        "DOMSIDE"
    ),
    detected_column = c(
        age_col,
        sex_col,
        bmi_col,
        site_col,
        duration_col,
        updrs_col,
        ledd_col,
        domside_col
    )
)

write.csv(
    detected_columns,
    file.path(out_dir, "Detected_input_columns.csv"),
    row.names = FALSE
)

cat("Detected columns:\n")
print(detected_columns)
cat("\n")

if (any(is.na(c(age_col, sex_col, bmi_col, site_col, duration_col, updrs_col, ledd_col, domside_col)))) {
    stop(
        "Could not detect one or more required columns. Check names(df_raw) and update candidate names.",
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Prepare raw PD-only longitudinal dataset
# ---------------------------------------------------------
visit_order <- c("BL", "V04", "V06", "V08", "V10", "V12", "V13", "V14")

dat_raw <- df_raw %>%
    filter(
        PRIMDIAG == 1,
        EVENT_ID %in% visit_order
    ) %>%
    mutate(
        PATNO = as.character(PATNO),
        EVENT_ID = as.character(EVENT_ID)
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
    
    duplicate_summary_by_visit <- duplicate_rows %>%
        count(EVENT_ID, PATNO, name = "n_rows") %>%
        group_by(EVENT_ID) %>%
        summarise(
            n_duplicate_participant_visits = n(),
            n_rows_in_duplicate_participant_visits = sum(n_rows),
            n_extra_rows_due_to_duplicates = sum(n_rows - 1),
            max_rows_per_participant_visit = max(n_rows),
            .groups = "drop"
        ) %>%
        arrange(EVENT_ID)
    
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
    
    duplicate_summary_by_visit <- tibble(
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
    duplicate_summary_by_visit,
    file.path(out_dir, "Duplicate_PATNO_EVENT_summary_by_visit.csv"),
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
        EVENT_ID = as.character(EVENT_ID)
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
# Prepare explicit modelling dataset
# ---------------------------------------------------------
dat0 <- dat_visit %>%
    mutate(
        PATNO_model = factor(PATNO),
        SITE_model = factor(.data[[site_col]]),
        EVENT_ID_model = factor(EVENT_ID, levels = visit_order),
        
        year_model = year_from_baseline,
        year_c_model = as.numeric(year_from_baseline),
        
        updrs3_score_on_model = as.numeric(.data[[updrs_col]]),
        
        baseline_NLR_model = as.numeric(baseline_NLR),
        baseline_NLR_z_model = as.numeric(baseline_NLR_z),
        log_baseline_NLR_model = ifelse(
            !is.na(baseline_NLR_model) & baseline_NLR_model > 0,
            log(baseline_NLR_model),
            NA_real_
        ),
        
        NLR_model = as.numeric(NLR),
        
        age_model = as.numeric(.data[[age_col]]),
        sex_model = factor(.data[[sex_col]]),
        bmi_model = as.numeric(.data[[bmi_col]]),
        duration_yrs_model = as.numeric(.data[[duration_col]]),
        
        LEDD_raw_model = as.numeric(.data[[ledd_col]]),
        LEDD_model = ifelse(year_c_model == 0, 0, LEDD_raw_model),
        
        DOMSIDE_model = factor(.data[[domside_col]])
    ) %>%
    select(
        PATNO = PATNO_model,
        SITE = SITE_model,
        EVENT_ID = EVENT_ID_model,
        PRIMDIAG,
        year = year_model,
        year_c = year_c_model,
        updrs3_score_on = updrs3_score_on_model,
        baseline_NLR = baseline_NLR_model,
        baseline_NLR_z = baseline_NLR_z_model,
        log_baseline_NLR = log_baseline_NLR_model,
        NLR = NLR_model,
        age = age_model,
        sex = sex_model,
        bmi = bmi_model,
        duration_yrs = duration_yrs_model,
        LEDD = LEDD_model,
        LEDD_raw = LEDD_raw_model,
        DOMSIDE = DOMSIDE_model
    )

dat0 <- droplevels(dat0)

# ---------------------------------------------------------
# Baseline UPDRS-III extraction and safe merge
# ---------------------------------------------------------
baseline_updrs <- dat0 %>%
    filter(EVENT_ID == "BL") %>%
    select(
        PATNO,
        baseline_updrs3_score_on = updrs3_score_on
    ) %>%
    distinct(PATNO, .keep_all = TRUE)

baseline_updrs_audit <- tibble(
    baseline_updrs_rows = nrow(baseline_updrs),
    baseline_updrs_unique_PATNO = n_distinct(baseline_updrs$PATNO),
    n_missing_baseline_updrs = sum(is.na(baseline_updrs$baseline_updrs3_score_on)),
    pct_missing_baseline_updrs =
        100 * mean(is.na(baseline_updrs$baseline_updrs3_score_on))
)

rows_before_baseline_updrs_join <- nrow(dat0)

dat0 <- dat0 %>%
    left_join(baseline_updrs, by = "PATNO")

rows_after_baseline_updrs_join <- nrow(dat0)

baseline_updrs_join_audit <- tibble(
    rows_before_baseline_updrs_join = rows_before_baseline_updrs_join,
    rows_after_baseline_updrs_join = rows_after_baseline_updrs_join,
    rows_added_by_join = rows_after_baseline_updrs_join - rows_before_baseline_updrs_join,
    baseline_updrs_rows = nrow(baseline_updrs),
    baseline_updrs_unique_PATNO = n_distinct(baseline_updrs$PATNO)
)

write.csv(
    baseline_updrs_audit,
    file.path(out_dir, "Baseline_UPDRS_audit.csv"),
    row.names = FALSE
)

write.csv(
    baseline_updrs_join_audit,
    file.path(out_dir, "Baseline_UPDRS_join_audit.csv"),
    row.names = FALSE
)

if (rows_after_baseline_updrs_join != rows_before_baseline_updrs_join) {
    stop(
        "Baseline UPDRS join multiplied rows. Check baseline_updrs uniqueness.",
        call. = FALSE
    )
}

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

predictor_covariate_vars_main <- c(
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

missingness_by_model_variable_main <- make_missingness_table(
    dat0,
    model_vars_main
)

missingness_by_predictor_covariate_main <- make_missingness_table(
    dat0,
    predictor_covariate_vars_main
)

dat0_missing_flags <- dat0 %>%
    ungroup() %>%
    mutate(
        row_has_missing_model_variable =
            if_any(all_of(model_vars_main), is.na),
        row_has_missing_predictor_covariate =
            if_any(all_of(predictor_covariate_vars_main), is.na)
    )

missingness_by_visit_main <- dat0_missing_flags %>%
    group_by(EVENT_ID, year_c) %>%
    summarise(
        n_rows = n(),
        n_subjects = n_distinct(PATNO),
        
        n_UPDRS_missing = sum(is.na(updrs3_score_on)),
        pct_UPDRS_missing = 100 * mean(is.na(updrs3_score_on)),
        
        n_baseline_NLR_z_missing = sum(is.na(baseline_NLR_z)),
        pct_baseline_NLR_z_missing = 100 * mean(is.na(baseline_NLR_z)),
        
        n_model_variable_missing_rows = sum(row_has_missing_model_variable),
        pct_model_variable_missing_rows = 100 * mean(row_has_missing_model_variable),
        
        n_predictor_or_covariate_missing_rows = sum(row_has_missing_predictor_covariate),
        pct_predictor_or_covariate_missing_rows = 100 * mean(row_has_missing_predictor_covariate),
        
        .groups = "drop"
    ) %>%
    arrange(year_c)

participant_missing_model_vars_main <- make_participant_missingness(
    data = dat0,
    vars = model_vars_main,
    label = "Model variables including outcome UPDRS-III"
)

participant_missing_predictors_main <- make_participant_missingness(
    data = dat0,
    vars = predictor_covariate_vars_main,
    label = "Predictors/covariates only; outcome UPDRS-III excluded"
)

participant_level_missingness_main <- bind_rows(
    participant_missing_model_vars_main$participant_level,
    participant_missing_predictors_main$participant_level
)

participant_level_missingness_summary_main <- bind_rows(
    participant_missing_model_vars_main$summary,
    participant_missing_predictors_main$summary
)

# ---------------------------------------------------------
# Complete-case dataset for main model
# ---------------------------------------------------------
dat <- dat0 %>%
    drop_na(all_of(model_vars_main))

dat <- droplevels(dat)

cat("Main model complete-case dataset:\n")
cat("Rows:", nrow(dat), "\n")
cat("Unique PD subjects:", n_distinct(dat$PATNO), "\n")
cat("Unique sites:", n_distinct(dat$SITE), "\n\n")

if (nrow(dat) < 10) {
    stop("Insufficient complete observations for the main model.", call. = FALSE)
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
    missingness_by_model_variable_main,
    file.path(out_dir, "Missingness_by_model_variable_main_UPDRS_including_outcome.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_predictor_covariate_main,
    file.path(out_dir, "Missingness_by_predictor_covariate_main_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_visit_main,
    file.path(out_dir, "Missingness_by_visit_main_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    participant_level_missingness_main,
    file.path(out_dir, "Participant_level_missingness_main_UPDRS_long.csv"),
    row.names = FALSE
)

write.csv(
    participant_level_missingness_summary_main,
    file.path(out_dir, "Participant_level_missingness_summary_main_UPDRS.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Descriptives
# ---------------------------------------------------------
desc <- dat %>%
    group_by(EVENT_ID, year_c) %>%
    summarise(
        n_rows = n(),
        n_subjects = n_distinct(PATNO),
        mean_UPDRS = mean(updrs3_score_on, na.rm = TRUE),
        sd_UPDRS = sd(updrs3_score_on, na.rm = TRUE),
        median_UPDRS = median(updrs3_score_on, na.rm = TRUE),
        IQR_UPDRS = IQR(updrs3_score_on, na.rm = TRUE),
        min_UPDRS = min(updrs3_score_on, na.rm = TRUE),
        max_UPDRS = max(updrs3_score_on, na.rm = TRUE),
        mean_baseline_NLR = mean(baseline_NLR, na.rm = TRUE),
        sd_baseline_NLR = sd(baseline_NLR, na.rm = TRUE),
        mean_baseline_NLR_z = mean(baseline_NLR_z, na.rm = TRUE),
        sd_baseline_NLR_z = sd(baseline_NLR_z, na.rm = TRUE),
        mean_LEDD = mean(LEDD, na.rm = TRUE),
        sd_LEDD = sd(LEDD, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(year_c)

write.csv(
    desc,
    file.path(out_dir, "PD_NLR_motor_progression_descriptives.csv"),
    row.names = FALSE
)

visit_counts <- dat %>%
    group_by(EVENT_ID, year_c) %>%
    summarise(
        n_rows = n(),
        n_subjects = n_distinct(PATNO),
        .groups = "drop"
    ) %>%
    arrange(year_c)

write.csv(
    visit_counts,
    file.path(out_dir, "Visit_counts_UPDRS_main_model.csv"),
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
    file.path(out_dir, "SITE_distribution_UPDRS_main_model.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Main mixed model
# ---------------------------------------------------------
# Sum-to-zero contrasts are used for Type III tests.
options(contrasts = c("contr.sum", "contr.poly"))

lmm <- safe_lmer(
    updrs3_score_on ~ baseline_NLR_z * year_c +
        age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
        (1 | SITE) + (1 | PATNO),
    data = dat,
    model_name = "Main UPDRS baseline NLR model"
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
    n_subjects = n_distinct(dat$PATNO),
    n_sites = n_distinct(dat$SITE)
)

write.csv(
    model_fit_indices,
    file.path(out_dir, "LMM_model_fit_indices_UPDRS_baseline_NLR.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Manual simple slopes of baseline NLR effect at each time point
# ---------------------------------------------------------
time_points <- sort(unique(dat$year_c))

slopes_by_time_df <- manual_simple_slopes_baseline_nlr(
    model = lmm,
    time_points = time_points
)

write.csv(
    slopes_by_time_df,
    file.path(out_dir, "Baseline_NLR_effect_on_UPDRS_by_timepoint.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Manual predicted trajectories for low / mean / high baseline NLR
# ---------------------------------------------------------
nlr_z_levels <- c(-1.5, 0, 1.5)

nlr_labels <- c(
    "Low baseline NLR (-1.5 SD)",
    "Mean baseline NLR",
    "High baseline NLR (+1.5 SD)"
)

prediction_reference <- tibble(
    age = mean(dat$age, na.rm = TRUE),
    sex = factor(levels(dat$sex)[1], levels = levels(dat$sex)),
    bmi = mean(dat$bmi, na.rm = TRUE),
    duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
    LEDD = mean(dat$LEDD, na.rm = TRUE),
    DOMSIDE = factor(levels(dat$DOMSIDE)[1], levels = levels(dat$DOMSIDE))
)

pred_grid_cont <- expand.grid(
    year_c = time_points,
    baseline_NLR_z = nlr_z_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
) %>%
    as_tibble() %>%
    mutate(
        age = prediction_reference$age[1],
        sex = prediction_reference$sex[1],
        bmi = prediction_reference$bmi[1],
        duration_yrs = prediction_reference$duration_yrs[1],
        LEDD = prediction_reference$LEDD[1],
        DOMSIDE = prediction_reference$DOMSIDE[1]
    ) %>%
    mutate(
        sex = factor(sex, levels = levels(dat$sex)),
        DOMSIDE = factor(DOMSIDE, levels = levels(dat$DOMSIDE))
    )

emm_cont_df <- manual_predictions_from_lmer(
    model = lmm,
    newdata = pred_grid_cont
) %>%
    mutate(
        NLR_group = factor(
            baseline_NLR_z,
            levels = nlr_z_levels,
            labels = nlr_labels
        )
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
    filter(EVENT_ID == "BL") %>%
    select(PATNO, baseline_NLR_z, baseline_NLR) %>%
    distinct(PATNO, .keep_all = TRUE) %>%
    filter(!is.na(baseline_NLR_z)) %>%
    mutate(
        NLR_tertile_number = ntile(baseline_NLR_z, 3),
        NLR_tertile = factor(
            NLR_tertile_number,
            levels = c(1, 2, 3),
            labels = c("Low NLR tertile", "Middle NLR tertile", "High NLR tertile")
        )
    )

dat_tert <- dat %>%
    left_join(
        bl_nlr_tertiles %>% select(PATNO, NLR_tertile),
        by = "PATNO"
    ) %>%
    filter(!is.na(NLR_tertile)) %>%
    mutate(
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
    group_by(NLR_tertile) %>%
    summarise(
        n_subjects = n(),
        baseline_NLR_z_mean = mean(baseline_NLR_z, na.rm = TRUE),
        baseline_NLR_z_sd = sd(baseline_NLR_z, na.rm = TRUE),
        baseline_NLR_mean = mean(baseline_NLR, na.rm = TRUE),
        baseline_NLR_sd = sd(baseline_NLR, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(match(
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

if (length(tert_means) > 0) {
    
    pred_grid_tert <- expand.grid(
        year_c = time_points,
        baseline_NLR_z = tert_means,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    ) %>%
        as_tibble() %>%
        mutate(
            age = prediction_reference$age[1],
            sex = prediction_reference$sex[1],
            bmi = prediction_reference$bmi[1],
            duration_yrs = prediction_reference$duration_yrs[1],
            LEDD = prediction_reference$LEDD[1],
            DOMSIDE = prediction_reference$DOMSIDE[1]
        ) %>%
        mutate(
            sex = factor(sex, levels = levels(dat$sex)),
            DOMSIDE = factor(DOMSIDE, levels = levels(dat$DOMSIDE))
        )
    
    emm_tert_df <- manual_predictions_from_lmer(
        model = lmm,
        newdata = pred_grid_tert
    ) %>%
        mutate(
            NLR_tertile = factor(
                baseline_NLR_z,
                levels = tert_means,
                labels = tert_labels
            )
        ) %>%
        filter(!is.na(NLR_tertile)) %>%
        mutate(
            NLR_tertile = droplevels(NLR_tertile)
        )
    
} else {
    emm_tert_df <- tibble()
}

write.csv(
    emm_tert_df,
    file.path(out_dir, "Predicted_UPDRS_trajectories_by_NLR_tertile.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Raw descriptive plot by baseline NLR tertiles
# ---------------------------------------------------------
raw_tert_df <- dat_tert %>%
    group_by(NLR_tertile, year_c, EVENT_ID) %>%
    summarise(
        mean_UPDRS = mean(updrs3_score_on, na.rm = TRUE),
        se_UPDRS = sd(updrs3_score_on, na.rm = TRUE) / sqrt(n()),
        n = n(),
        n_subjects = n_distinct(PATNO),
        .groups = "drop"
    ) %>%
    filter(!is.na(NLR_tertile))

write.csv(
    raw_tert_df,
    file.path(out_dir, "Raw_UPDRS_by_baseline_NLR_tertile.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Main model diagnostic outputs
# ---------------------------------------------------------
diagnostic_df <- dat %>%
    mutate(
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
    filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_PATNO_UPDRS_baseline_NLR.png"),
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
    filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_SITE_UPDRS_baseline_NLR.png"),
    plot = p_qq_ranef_site,
    width = 7,
    height = 5,
    dpi = 600
)

safe_check_model_png(lmm, "Performance_check_model_UPDRS_baseline_NLR.png")

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
if (nrow(emm_tert_df) > 0) {
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
}

# ---------------------------------------------------------
# Figure 3: raw descriptive trajectories by tertile
# ---------------------------------------------------------
if (nrow(raw_tert_df) > 0) {
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
}

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

predictor_covariate_vars_sens_baseline_adj <- c(
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

dat_baseline_adj0 <- dat0 %>%
    filter(year_c > 0)

missingness_by_model_variable_baseline_adj <- make_missingness_table(
    dat_baseline_adj0,
    model_vars_sens_baseline_adj
)

missingness_by_predictor_covariate_baseline_adj <- make_missingness_table(
    dat_baseline_adj0,
    predictor_covariate_vars_sens_baseline_adj
)

dat_baseline_adj <- dat_baseline_adj0 %>%
    drop_na(all_of(model_vars_sens_baseline_adj))

dat_baseline_adj <- droplevels(dat_baseline_adj)

lmm_baseline_adj <- safe_lmer(
    updrs3_score_on ~ baseline_NLR_z * year_c +
        baseline_updrs3_score_on +
        age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
        (1 | SITE) + (1 | PATNO),
    data = dat_baseline_adj,
    model_name = "Baseline UPDRS-adjusted follow-up-only model"
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
    missingness_by_model_variable_baseline_adj,
    file.path(out_dir, "Sensitivity_baseline_adjusted_missingness_by_model_variable_including_outcome.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_predictor_covariate_baseline_adj,
    file.path(out_dir, "Sensitivity_baseline_adjusted_missingness_by_predictor_covariate.csv"),
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
    n_subjects = n_distinct(dat_baseline_adj$PATNO),
    n_sites = n_distinct(dat_baseline_adj$SITE)
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

predictor_covariate_vars_log <- c(
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

missingness_by_model_variable_log <- make_missingness_table(
    dat0,
    model_vars_log
)

missingness_by_predictor_covariate_log <- make_missingness_table(
    dat0,
    predictor_covariate_vars_log
)

dat_log <- dat0 %>%
    drop_na(all_of(model_vars_log))

dat_log <- droplevels(dat_log)

lmm_log <- safe_lmer(
    updrs3_score_on ~ log_baseline_NLR * year_c +
        age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
        (1 | SITE) + (1 | PATNO),
    data = dat_log,
    model_name = "Log baseline NLR UPDRS model"
)

fixef_log <- make_fixef_table(lmm_log)
anova_log <- make_anova_table(lmm_log)

log_fit_indices <- tibble(
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
    dat_log,
    file.path(out_dir, "Sensitivity_log_baseline_NLR_dataset.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_model_variable_log,
    file.path(out_dir, "Sensitivity_log_baseline_NLR_missingness_by_model_variable_including_outcome.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_predictor_covariate_log,
    file.path(out_dir, "Sensitivity_log_baseline_NLR_missingness_by_predictor_covariate.csv"),
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

write.csv(
    log_fit_indices,
    file.path(out_dir, "Sensitivity_log_baseline_NLR_model_fit_indices_UPDRS.csv"),
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
add_sheet(wb, "Duplicate_by_visit", duplicate_summary_by_visit)
add_sheet(wb, "Duplicate_discrepancy", duplicate_discrepancy_by_variable)
add_sheet(wb, "Cleaning_summary", cleaning_summary)
add_sheet(wb, "Baseline_UPDRS_audit", baseline_updrs_audit)
add_sheet(wb, "Baseline_UPDRS_join", baseline_updrs_join_audit)

add_sheet(wb, "Missing_model_main", missingness_by_model_variable_main)
add_sheet(wb, "Missing_predictors_main", missingness_by_predictor_covariate_main)
add_sheet(wb, "Missing_by_visit_main", missingness_by_visit_main)
add_sheet(wb, "Participant_missing", participant_level_missingness_summary_main)

add_sheet(wb, "Main_fit", model_fit_indices)
add_sheet(wb, "Main_fixed_effects", fixef_tab)
add_sheet(wb, "Main_TypeIII_ANOVA", anova_tab)
add_sheet(wb, "Main_random_effects", randef_var)
add_sheet(wb, "Main_descriptives", desc)
add_sheet(wb, "Main_visit_counts", visit_counts)
add_sheet(wb, "Main_site_counts", site_counts)
add_sheet(wb, "Slopes_by_time", slopes_by_time_df)
add_sheet(wb, "Pred_continuous", emm_cont_df)
add_sheet(wb, "Tertile_summary", tert_summary)
add_sheet(wb, "Pred_tertiles", emm_tert_df)
add_sheet(wb, "Raw_tertiles", raw_tert_df)

add_sheet(wb, "Sens_BLadj_fit", baseline_adj_fit_indices)
add_sheet(wb, "Sens_BLadj_fixed", fixef_baseline_adj)
add_sheet(wb, "Sens_BLadj_ANOVA", anova_baseline_adj)
add_sheet(wb, "Sens_BLadj_miss_model", missingness_by_model_variable_baseline_adj)
add_sheet(wb, "Sens_BLadj_miss_pred", missingness_by_predictor_covariate_baseline_adj)

add_sheet(wb, "Sens_log_fit", log_fit_indices)
add_sheet(wb, "Sens_log_fixed", fixef_log)
add_sheet(wb, "Sens_log_ANOVA", anova_log)
add_sheet(wb, "Sens_log_miss_model", missingness_by_model_variable_log)
add_sheet(wb, "Sens_log_miss_pred", missingness_by_predictor_covariate_log)

openxlsx::saveWorkbook(
    wb,
    file = file.path(out_dir, "UPDRS_baseline_NLR_outputs.xlsx"),
    overwrite = TRUE
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
print(detected_columns, row.names = FALSE)
cat("\n\n")

cat("Duplicate audit before model cleaning:\n")
print(duplicate_audit_summary, row.names = FALSE)
cat("\n\n")

cat("Participant-visit cleaning summary:\n")
print(cleaning_summary, row.names = FALSE)
cat("\n\n")

cat("Baseline UPDRS audit:\n")
print(baseline_updrs_audit, row.names = FALSE)
cat("\n\n")

cat("Baseline UPDRS join audit:\n")
print(baseline_updrs_join_audit, row.names = FALSE)
cat("\n\n")

cat("Sample size before complete-case filtering:\n")
cat("Rows =", nrow(dat0), "\n")
cat("Unique PD subjects =", n_distinct(dat0$PATNO), "\n")
cat("Unique sites =", n_distinct(dat0$SITE), "\n")
cat("Unique participant-visits =", n_distinct(paste(dat0$PATNO, dat0$EVENT_ID)), "\n\n")

cat("Sample size in main complete-case model:\n")
cat("Rows =", nrow(dat), "\n")
cat("Unique PD subjects =", n_distinct(dat$PATNO), "\n")
cat("Unique sites =", n_distinct(dat$SITE), "\n")
cat("Unique participant-visits =", n_distinct(paste(dat$PATNO, dat$EVENT_ID)), "\n\n")

cat("Visit counts in main model:\n")
print(table(dat$EVENT_ID))
cat("\n\n")

cat("Missingness by model variable, main model, including outcome UPDRS-III:\n")
print(missingness_by_model_variable_main, row.names = FALSE)
cat("\n\n")

cat("Missingness by predictors/covariates only, main model, excluding outcome UPDRS-III:\n")
print(missingness_by_predictor_covariate_main, row.names = FALSE)
cat("\n\n")

cat("Missingness by visit, main model:\n")
print(missingness_by_visit_main, row.names = FALSE)
cat("\n\n")

cat("Participant-level missingness summary, main model:\n")
print(participant_level_missingness_summary_main, row.names = FALSE)
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
cat("Manual linear combinations of fixed effects: baseline_NLR_z + baseline_NLR_z:year_c * year_c.\n")
print(slopes_by_time_df, row.names = FALSE)
cat("\n\n")

cat("NLR levels used in continuous interaction figure:\n")
print(data.frame(
    label = nlr_labels,
    baseline_NLR_z = nlr_z_levels
), row.names = FALSE)
cat("\n\n")

cat("Prediction reference values for manual predicted trajectories:\n")
print(prediction_reference, row.names = FALSE)
cat("\n\n")

cat("Baseline NLR tertile summary:\n")
print(tert_summary, row.names = FALSE)
cat("\n\n")

cat("Sensitivity model: baseline outcome-adjusted follow-up-only model\n")
cat("Formula:\n")
cat("updrs3_score_on ~ baseline_NLR_z * year_c + baseline_updrs3_score_on + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
cat("Rows =", nrow(dat_baseline_adj), "\n")
cat("Unique PD subjects =", n_distinct(dat_baseline_adj$PATNO), "\n")
cat("Unique sites =", n_distinct(dat_baseline_adj$SITE), "\n\n")
cat("Missingness by model variable, including outcome:\n")
print(missingness_by_model_variable_baseline_adj, row.names = FALSE)
cat("\n\n")
cat("Missingness by predictors/covariates only:\n")
print(missingness_by_predictor_covariate_baseline_adj, row.names = FALSE)
cat("\n\n")
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
cat("Unique PD subjects =", n_distinct(dat_log$PATNO), "\n")
cat("Unique sites =", n_distinct(dat_log$SITE), "\n\n")
cat("Missingness by model variable, including outcome:\n")
print(missingness_by_model_variable_log, row.names = FALSE)
cat("\n\n")
cat("Missingness by predictors/covariates only:\n")
print(missingness_by_predictor_covariate_log, row.names = FALSE)
cat("\n\n")
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

# ---------------------------------------------------------
# Final message
# ---------------------------------------------------------
cat("\n============================================================\n")
cat("Baseline NLR -> UPDRS-III motor progression analysis completed.\n")
cat("All outputs were saved to:\n")
cat(out_dir, "\n\n")
cat("Main summary file:\n")
cat(file.path(out_dir, "LMM_UPDRS_baseline_NLR_summary.txt"), "\n\n")
cat("Excel workbook:\n")
cat(file.path(out_dir, "UPDRS_baseline_NLR_outputs.xlsx"), "\n")
cat("============================================================\n")