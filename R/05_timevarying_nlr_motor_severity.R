# =========================================================
# 05_timevarying_nlr_motor_severity.R
#
# Time-varying neutrophil-to-lymphocyte ratio (NLR) and
# longitudinal motor severity in Parkinson's disease.
#
# Outcome:
#   ON-medication MDS-UPDRS Part III / UPDRS-III
#
# Main NLR components:
#   NLR_between_z
#       Between-person component: participant-specific mean NLR
#       across available visits, standardized.
#
#   NLR_within_z
#       Within-person component: visit-wise deviation from the
#       participant-specific mean NLR, standardized.
#
# Primary model:
#   updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c +
#       age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
#       (1 | SITE) + (1 | PATNO)
#
# Secondary model:
#   updrs3_score_on ~ NLR_between_z + NLR_within_z * year_c +
#       age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
#       (1 | SITE) + (1 | PATNO)
#
# Sensitivity analyses:
#   1) Follow-up-only model adjusted for baseline UPDRS-III.
#   2) Log-transformed NLR between-/within-person decomposition.
#
# Predicted trajectories and within-person simple slopes are
# calculated manually from the fixed-effect estimates and model
# covariance matrix; this script intentionally does not rely on
# emmeans() or emtrends().
#
# Default input:
#   <PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
#
# Output:
#   <PPMI_OUTPUT_ROOT>/05_TIMEVARYING_NLR_UPDRS
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

out_dir <- file.path(output_root, "05_TIMEVARYING_NLR_UPDRS")

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

manual_simple_slopes <- function(model, focal_term, moderator_term, moderator_values) {
    
    beta <- lme4::fixef(model)
    vc <- as.matrix(vcov(model))
    
    term_int <- get_interaction_term(
        term_names = names(beta),
        term_a = focal_term,
        term_b = moderator_term
    )
    
    if (!focal_term %in% names(beta)) {
        stop(
            paste0(focal_term, " fixed-effect term not found."),
            call. = FALSE
        )
    }
    
    out <- lapply(
        moderator_values,
        function(t) {
            
            L <- rep(0, length(beta))
            names(L) <- names(beta)
            
            L[focal_term] <- 1
            L[term_int] <- t
            
            estimate <- sum(L * beta)
            SE <- sqrt(as.numeric(t(L) %*% vc %*% L))
            z_value <- estimate / SE
            p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
            
            tibble(
                year_c = t,
                effect = estimate,
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
    
    if (length(missing_cols) > 0) {
        stop(
            paste0(
                "Prediction model matrix is missing fixed-effect columns:\n",
                paste(missing_cols, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    
    X <- X[, names(beta), drop = FALSE]
    
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
# Required variables
# ---------------------------------------------------------
required_vars <- c(
    "PATNO",
    "EVENT_ID",
    "PRIMDIAG",
    "year_from_baseline",
    "NLR"
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
# Prepare raw PD-only dataset
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
        
        NLR_model = as.numeric(NLR),
        log_NLR_model = ifelse(!is.na(NLR_model) & NLR_model > 0, log(NLR_model), NA_real_),
        
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
        NLR = NLR_model,
        log_NLR = log_NLR_model,
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
# Decompose visit-wise NLR into between-person and within-person components
# ---------------------------------------------------------
person_nlr <- dat0 %>%
    group_by(PATNO) %>%
    summarise(
        NLR_between_raw = mean(NLR, na.rm = TRUE),
        n_NLR_available = sum(!is.na(NLR)),
        .groups = "drop"
    ) %>%
    mutate(
        NLR_between_raw = ifelse(is.nan(NLR_between_raw), NA_real_, NLR_between_raw)
    )

between_mean <- mean(person_nlr$NLR_between_raw, na.rm = TRUE)
between_sd <- sd(person_nlr$NLR_between_raw, na.rm = TRUE)

dat0 <- dat0 %>%
    left_join(person_nlr, by = "PATNO") %>%
    mutate(
        NLR_within_raw = NLR - NLR_between_raw
    )

within_mean <- mean(dat0$NLR_within_raw, na.rm = TRUE)
within_sd <- sd(dat0$NLR_within_raw, na.rm = TRUE)

dat0 <- dat0 %>%
    mutate(
        NLR_between_z = as.numeric((NLR_between_raw - between_mean) / between_sd),
        NLR_within_z = as.numeric((NLR_within_raw - within_mean) / within_sd)
    )

nlr_decomposition_audit <- tibble(
    n_participants_in_person_nlr = nrow(person_nlr),
    n_participants_with_at_least_one_NLR = sum(person_nlr$n_NLR_available > 0),
    between_mean = between_mean,
    between_sd = between_sd,
    within_mean = within_mean,
    within_sd = within_sd,
    n_missing_NLR_between_z = sum(is.na(dat0$NLR_between_z)),
    pct_missing_NLR_between_z = 100 * mean(is.na(dat0$NLR_between_z)),
    n_missing_NLR_within_z = sum(is.na(dat0$NLR_within_z)),
    pct_missing_NLR_within_z = 100 * mean(is.na(dat0$NLR_within_z))
)

write.csv(
    person_nlr,
    file.path(out_dir, "Participant_mean_NLR_between_component.csv"),
    row.names = FALSE
)

write.csv(
    nlr_decomposition_audit,
    file.path(out_dir, "NLR_between_within_decomposition_audit.csv"),
    row.names = FALSE
)

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

predictor_covariate_vars_main <- c(
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
        
        n_NLR_missing = sum(is.na(NLR)),
        pct_NLR_missing = 100 * mean(is.na(NLR)),
        
        n_NLR_between_z_missing = sum(is.na(NLR_between_z)),
        pct_NLR_between_z_missing = 100 * mean(is.na(NLR_between_z)),
        
        n_NLR_within_z_missing = sum(is.na(NLR_within_z)),
        pct_NLR_within_z_missing = 100 * mean(is.na(NLR_within_z)),
        
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

cat("Main time-varying NLR model complete-case dataset:\n")
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
    file.path(out_dir, "PD_timevarying_NLR_UPDRS_dataset_before_complete_case_filter.csv"),
    row.names = FALSE
)

write.csv(
    dat,
    file.path(out_dir, "PD_timevarying_NLR_UPDRS_complete_case_dataset.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_model_variable_main,
    file.path(out_dir, "Missingness_by_model_variable_timevarying_NLR_UPDRS_including_outcome.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_predictor_covariate_main,
    file.path(out_dir, "Missingness_by_predictor_covariate_timevarying_NLR_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_visit_main,
    file.path(out_dir, "Missingness_by_visit_timevarying_NLR_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    participant_level_missingness_main,
    file.path(out_dir, "Participant_level_missingness_timevarying_NLR_UPDRS_long.csv"),
    row.names = FALSE
)

write.csv(
    participant_level_missingness_summary_main,
    file.path(out_dir, "Participant_level_missingness_summary_timevarying_NLR_UPDRS.csv"),
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
    arrange(year_c)

write.csv(
    desc,
    file.path(out_dir, "PD_timevarying_NLR_UPDRS_descriptives.csv"),
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
    file.path(out_dir, "Visit_counts_timevarying_NLR_UPDRS.csv"),
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
    file.path(out_dir, "SITE_distribution_timevarying_NLR_UPDRS.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Main model:
# Between-person and within-person time-varying NLR components
# ---------------------------------------------------------
# Sum-to-zero contrasts are used for Type III tests.
options(contrasts = c("contr.sum", "contr.poly"))

lmm_main <- safe_lmer(
    updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c +
        age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
        (1 | SITE) + (1 | PATNO),
    data = dat,
    model_name = "Main time-varying NLR UPDRS model"
)

fixef_main <- make_fixef_table(lmm_main)
anova_main <- make_anova_table(lmm_main)
randef_main <- as.data.frame(VarCorr(lmm_main))

write.csv(
    fixef_main,
    file.path(out_dir, "LMM_fixed_effects_timevarying_NLR_UPDRS_main.csv"),
    row.names = FALSE
)

write.csv(
    anova_main,
    file.path(out_dir, "LMM_TypeIII_ANOVA_timevarying_NLR_UPDRS_main.csv"),
    row.names = FALSE
)

write.csv(
    randef_main,
    file.path(out_dir, "LMM_random_effects_variance_timevarying_NLR_UPDRS_main.csv"),
    row.names = FALSE
)

model_fit_main <- tibble(
    model = "main_between_within_NLR",
    AIC = AIC(lmm_main),
    BIC = BIC(lmm_main),
    logLik = as.numeric(logLik(lmm_main)),
    deviance = deviance(lmm_main),
    sigma = sigma(lmm_main),
    n_obs = nobs(lmm_main),
    n_subjects = n_distinct(dat$PATNO),
    n_sites = n_distinct(dat$SITE)
)

write.csv(
    model_fit_main,
    file.path(out_dir, "LMM_model_fit_indices_timevarying_NLR_UPDRS_main.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Secondary model:
# Does within-person NLR-UPDRS association vary over time?
# ---------------------------------------------------------
lmm_within_interaction <- safe_lmer(
    updrs3_score_on ~ NLR_between_z + NLR_within_z * year_c +
        age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
        (1 | SITE) + (1 | PATNO),
    data = dat,
    model_name = "Secondary within-person NLR by time UPDRS model"
)

fixef_within_interaction <- make_fixef_table(lmm_within_interaction)
anova_within_interaction <- make_anova_table(lmm_within_interaction)
randef_within_interaction <- as.data.frame(VarCorr(lmm_within_interaction))

write.csv(
    fixef_within_interaction,
    file.path(out_dir, "Secondary_withinNLRxTime_LMM_fixed_effects_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    anova_within_interaction,
    file.path(out_dir, "Secondary_withinNLRxTime_LMM_TypeIII_ANOVA_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    randef_within_interaction,
    file.path(out_dir, "Secondary_withinNLRxTime_random_effects_UPDRS.csv"),
    row.names = FALSE
)

model_fit_within_interaction <- tibble(
    model = "secondary_within_NLR_x_time",
    AIC = AIC(lmm_within_interaction),
    BIC = BIC(lmm_within_interaction),
    logLik = as.numeric(logLik(lmm_within_interaction)),
    deviance = deviance(lmm_within_interaction),
    sigma = sigma(lmm_within_interaction),
    n_obs = nobs(lmm_within_interaction),
    n_subjects = n_distinct(dat$PATNO),
    n_sites = n_distinct(dat$SITE)
)

write.csv(
    model_fit_within_interaction,
    file.path(out_dir, "Secondary_withinNLRxTime_model_fit_indices_UPDRS.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Manual simple slopes for within-person NLR effect at each time point
# ---------------------------------------------------------
time_points <- sort(unique(dat$year_c))

within_slopes_by_time_df <- manual_simple_slopes(
    model = lmm_within_interaction,
    focal_term = "NLR_within_z",
    moderator_term = "year_c",
    moderator_values = time_points
) %>%
    rename(
        NLR_within_z_effect = effect
    )

write.csv(
    within_slopes_by_time_df,
    file.path(out_dir, "Within_person_NLR_effect_on_UPDRS_by_timepoint.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Sensitivity model:
# Baseline outcome-adjusted follow-up-only model
# ---------------------------------------------------------
model_vars_baseline_adj <- c(
    "updrs3_score_on",
    "baseline_updrs3_score_on",
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

predictor_covariate_vars_baseline_adj <- c(
    "baseline_updrs3_score_on",
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

dat_baseline_adj0 <- dat0 %>%
    filter(year_c > 0)

missingness_by_model_variable_baseline_adj <- make_missingness_table(
    dat_baseline_adj0,
    model_vars_baseline_adj
)

missingness_by_predictor_covariate_baseline_adj <- make_missingness_table(
    dat_baseline_adj0,
    predictor_covariate_vars_baseline_adj
)

dat_baseline_adj <- dat_baseline_adj0 %>%
    drop_na(all_of(model_vars_baseline_adj))

dat_baseline_adj <- droplevels(dat_baseline_adj)

lmm_baseline_adj <- safe_lmer(
    updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c +
        baseline_updrs3_score_on +
        age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
        (1 | SITE) + (1 | PATNO),
    data = dat_baseline_adj,
    model_name = "Baseline UPDRS-adjusted follow-up-only time-varying NLR model"
)

fixef_baseline_adj <- make_fixef_table(lmm_baseline_adj)
anova_baseline_adj <- make_anova_table(lmm_baseline_adj)
randef_baseline_adj <- as.data.frame(VarCorr(lmm_baseline_adj))

write.csv(
    dat_baseline_adj,
    file.path(out_dir, "Sensitivity_baseline_adjusted_followup_only_dataset_timevarying_NLR_UPDRS.csv"),
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
    file.path(out_dir, "Sensitivity_baseline_adjusted_LMM_fixed_effects_timevarying_NLR_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    anova_baseline_adj,
    file.path(out_dir, "Sensitivity_baseline_adjusted_LMM_TypeIII_ANOVA_timevarying_NLR_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    randef_baseline_adj,
    file.path(out_dir, "Sensitivity_baseline_adjusted_random_effects_timevarying_NLR_UPDRS.csv"),
    row.names = FALSE
)

model_fit_baseline_adj <- tibble(
    model = "sensitivity_baseline_adjusted_followup_only",
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
    model_fit_baseline_adj,
    file.path(out_dir, "Sensitivity_baseline_adjusted_model_fit_indices_timevarying_NLR_UPDRS.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Sensitivity model:
# Log-NLR between/within decomposition
# ---------------------------------------------------------
person_log_nlr <- dat0 %>%
    group_by(PATNO) %>%
    summarise(
        log_NLR_between_raw = mean(log_NLR, na.rm = TRUE),
        n_log_NLR_available = sum(!is.na(log_NLR)),
        .groups = "drop"
    ) %>%
    mutate(
        log_NLR_between_raw = ifelse(is.nan(log_NLR_between_raw), NA_real_, log_NLR_between_raw)
    )

dat_log0 <- dat0 %>%
    left_join(person_log_nlr, by = "PATNO") %>%
    mutate(
        log_NLR_within_raw = log_NLR - log_NLR_between_raw
    )

log_between_mean <- mean(person_log_nlr$log_NLR_between_raw, na.rm = TRUE)
log_between_sd <- sd(person_log_nlr$log_NLR_between_raw, na.rm = TRUE)

log_within_mean <- mean(dat_log0$log_NLR_within_raw, na.rm = TRUE)
log_within_sd <- sd(dat_log0$log_NLR_within_raw, na.rm = TRUE)

dat_log0 <- dat_log0 %>%
    mutate(
        log_NLR_between_z = as.numeric((log_NLR_between_raw - log_between_mean) / log_between_sd),
        log_NLR_within_z = as.numeric((log_NLR_within_raw - log_within_mean) / log_within_sd)
    )

log_nlr_decomposition_audit <- tibble(
    n_participants_in_person_log_nlr = nrow(person_log_nlr),
    n_participants_with_at_least_one_log_NLR = sum(person_log_nlr$n_log_NLR_available > 0),
    log_between_mean = log_between_mean,
    log_between_sd = log_between_sd,
    log_within_mean = log_within_mean,
    log_within_sd = log_within_sd,
    n_missing_log_NLR_between_z = sum(is.na(dat_log0$log_NLR_between_z)),
    pct_missing_log_NLR_between_z = 100 * mean(is.na(dat_log0$log_NLR_between_z)),
    n_missing_log_NLR_within_z = sum(is.na(dat_log0$log_NLR_within_z)),
    pct_missing_log_NLR_within_z = 100 * mean(is.na(dat_log0$log_NLR_within_z))
)

write.csv(
    person_log_nlr,
    file.path(out_dir, "Participant_mean_logNLR_between_component.csv"),
    row.names = FALSE
)

write.csv(
    log_nlr_decomposition_audit,
    file.path(out_dir, "LogNLR_between_within_decomposition_audit.csv"),
    row.names = FALSE
)

model_vars_log <- c(
    "updrs3_score_on",
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

predictor_covariate_vars_log <- c(
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

missingness_by_model_variable_log <- make_missingness_table(
    dat_log0,
    model_vars_log
)

missingness_by_predictor_covariate_log <- make_missingness_table(
    dat_log0,
    predictor_covariate_vars_log
)

dat_log <- dat_log0 %>%
    drop_na(all_of(model_vars_log))

dat_log <- droplevels(dat_log)

lmm_log <- safe_lmer(
    updrs3_score_on ~ log_NLR_between_z + log_NLR_within_z + year_c +
        age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
        (1 | SITE) + (1 | PATNO),
    data = dat_log,
    model_name = "Log-NLR between-within UPDRS model"
)

fixef_log <- make_fixef_table(lmm_log)
anova_log <- make_anova_table(lmm_log)

model_fit_log <- tibble(
    model = "sensitivity_logNLR_between_within",
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
    file.path(out_dir, "Sensitivity_logNLR_between_within_dataset_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_model_variable_log,
    file.path(out_dir, "Sensitivity_logNLR_missingness_by_model_variable_including_outcome.csv"),
    row.names = FALSE
)

write.csv(
    missingness_by_predictor_covariate_log,
    file.path(out_dir, "Sensitivity_logNLR_missingness_by_predictor_covariate.csv"),
    row.names = FALSE
)

write.csv(
    fixef_log,
    file.path(out_dir, "Sensitivity_logNLR_between_within_LMM_fixed_effects_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    anova_log,
    file.path(out_dir, "Sensitivity_logNLR_between_within_LMM_TypeIII_ANOVA_UPDRS.csv"),
    row.names = FALSE
)

write.csv(
    model_fit_log,
    file.path(out_dir, "Sensitivity_logNLR_between_within_model_fit_indices_UPDRS.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Model diagnostics for main model
# ---------------------------------------------------------
diagnostic_df <- dat %>%
    mutate(
        fitted_value = fitted(lmm_main),
        residual = resid(lmm_main),
        pearson_residual = residual / sigma(lmm_main)
    )

write.csv(
    diagnostic_df,
    file.path(out_dir, "Diagnostic_values_timevarying_NLR_UPDRS_main.csv"),
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
    filename = file.path(out_dir, "Diagnostic_residuals_vs_fitted_timevarying_NLR_UPDRS.png"),
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
    filename = file.path(out_dir, "Diagnostic_QQ_residuals_timevarying_NLR_UPDRS.png"),
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
    filename = file.path(out_dir, "Diagnostic_observed_vs_fitted_timevarying_NLR_UPDRS.png"),
    plot = p_obs_fit,
    width = 7,
    height = 5,
    dpi = 600
)

ranef_patno <- ranef(lmm_main)$PATNO %>%
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
    filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR_UPDRS.png"),
    plot = p_qq_ranef_patno,
    width = 7,
    height = 5,
    dpi = 600
)

ranef_site <- ranef(lmm_main)$SITE %>%
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
    filename = file.path(out_dir, "Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR_UPDRS.png"),
    plot = p_qq_ranef_site,
    width = 7,
    height = 5,
    dpi = 600
)

safe_check_model_png(lmm_main, "Performance_check_model_timevarying_NLR_UPDRS.png")

# ---------------------------------------------------------
# Visualization 1:
# Descriptive within-person NLR deviation and UPDRS
# ---------------------------------------------------------
p_within_scatter <- ggplot(
    dat,
    aes(x = NLR_within_z, y = updrs3_score_on)
) +
    geom_point(alpha = 0.35, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
    theme_classic(base_size = 13) +
    labs(
        title = "Visit-wise within-person NLR deviation and ON-medication UPDRS-III",
        subtitle = "Descriptive plot; inference based on the mixed-effects model",
        x = "Within-person NLR deviation, standardized",
        y = "ON-medication UPDRS-III score"
    )

ggsave(
    filename = file.path(out_dir, "Descriptive_within_person_NLR_vs_UPDRS.png"),
    plot = p_within_scatter,
    width = 7,
    height = 5,
    dpi = 600
)

ggsave(
    filename = file.path(out_dir, "Descriptive_within_person_NLR_vs_UPDRS.pdf"),
    plot = p_within_scatter,
    width = 7,
    height = 5
)

# ---------------------------------------------------------
# Visualization 2:
# Manual predicted UPDRS at low / mean / high within-person NLR
# from the main model, holding between-person NLR at mean
# ---------------------------------------------------------
within_levels <- c(-1.5, 0, 1.5)

within_labels <- c(
    "Lower-than-usual NLR (-1.5 SD)",
    "Usual NLR",
    "Higher-than-usual NLR (+1.5 SD)"
)

prediction_reference <- tibble(
    age = mean(dat$age, na.rm = TRUE),
    sex = factor(levels(dat$sex)[1], levels = levels(dat$sex)),
    bmi = mean(dat$bmi, na.rm = TRUE),
    duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
    LEDD = mean(dat$LEDD, na.rm = TRUE),
    DOMSIDE = factor(levels(dat$DOMSIDE)[1], levels = levels(dat$DOMSIDE))
)

pred_grid_within <- expand.grid(
    year_c = time_points,
    NLR_within_z = within_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
) %>%
    as_tibble() %>%
    mutate(
        NLR_between_z = 0,
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

emm_within_df <- manual_predictions_from_lmer(
    model = lmm_main,
    newdata = pred_grid_within
) %>%
    mutate(
        within_NLR_level = factor(
            NLR_within_z,
            levels = within_levels,
            labels = within_labels
        )
    )

write.csv(
    emm_within_df,
    file.path(out_dir, "Predicted_UPDRS_by_within_person_NLR_levels.csv"),
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
        breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
        labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
    ) +
    labs(
        title = "Time-varying NLR and ON-medication UPDRS-III",
        subtitle = "Predicted UPDRS-III at lower-than-usual, usual, and higher-than-usual NLR",
        x = "Years from baseline",
        y = "Predicted ON-medication UPDRS-III score",
        color = "Visit-wise NLR",
        fill = "Visit-wise NLR"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

ggsave(
    filename = file.path(out_dir, "Figure_timevarying_withinNLR_UPDRS_predicted.png"),
    plot = p_within_pred,
    width = 9.5,
    height = 6.8,
    dpi = 600
)

ggsave(
    filename = file.path(out_dir, "Figure_timevarying_withinNLR_UPDRS_predicted.pdf"),
    plot = p_within_pred,
    width = 9.5,
    height = 6.8
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
add_sheet(wb, "NLR_decomp_audit", nlr_decomposition_audit)
add_sheet(wb, "Person_NLR", person_nlr)
add_sheet(wb, "Baseline_UPDRS_audit", baseline_updrs_audit)
add_sheet(wb, "Baseline_UPDRS_join", baseline_updrs_join_audit)

add_sheet(wb, "Missing_model_main", missingness_by_model_variable_main)
add_sheet(wb, "Missing_predictors_main", missingness_by_predictor_covariate_main)
add_sheet(wb, "Missing_by_visit_main", missingness_by_visit_main)
add_sheet(wb, "Participant_missing", participant_level_missingness_summary_main)

add_sheet(wb, "Main_fit", model_fit_main)
add_sheet(wb, "Main_fixed_effects", fixef_main)
add_sheet(wb, "Main_TypeIII_ANOVA", anova_main)
add_sheet(wb, "Main_random_effects", randef_main)
add_sheet(wb, "Main_descriptives", desc)
add_sheet(wb, "Main_visit_counts", visit_counts)
add_sheet(wb, "Main_site_counts", site_counts)

add_sheet(wb, "Secondary_fit", model_fit_within_interaction)
add_sheet(wb, "Secondary_fixed", fixef_within_interaction)
add_sheet(wb, "Secondary_ANOVA", anova_within_interaction)
add_sheet(wb, "Within_slopes_time", within_slopes_by_time_df)

add_sheet(wb, "Sens_BLadj_fit", model_fit_baseline_adj)
add_sheet(wb, "Sens_BLadj_fixed", fixef_baseline_adj)
add_sheet(wb, "Sens_BLadj_ANOVA", anova_baseline_adj)
add_sheet(wb, "Sens_BLadj_miss_model", missingness_by_model_variable_baseline_adj)
add_sheet(wb, "Sens_BLadj_miss_pred", missingness_by_predictor_covariate_baseline_adj)

add_sheet(wb, "LogNLR_decomp_audit", log_nlr_decomposition_audit)
add_sheet(wb, "Sens_log_fit", model_fit_log)
add_sheet(wb, "Sens_log_fixed", fixef_log)
add_sheet(wb, "Sens_log_ANOVA", anova_log)
add_sheet(wb, "Sens_log_miss_model", missingness_by_model_variable_log)
add_sheet(wb, "Sens_log_miss_pred", missingness_by_predictor_covariate_log)

add_sheet(wb, "Pred_within_NLR", emm_within_df)

openxlsx::saveWorkbook(
    wb,
    file = file.path(out_dir, "Timevarying_NLR_UPDRS_outputs.xlsx"),
    overwrite = TRUE
)

# ---------------------------------------------------------
# Full text summary
# ---------------------------------------------------------
sink(file.path(out_dir, "LMM_timevarying_NLR_UPDRS_summary.txt"))

cat("============================================================\n")
cat("Aim 2B / Reviewer 1 Comment 2:\n")
cat("Does time-varying NLR track ON-medication UPDRS-III?\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(file_path, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Main model formula:\n")
cat("updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Secondary model formula:\n")
cat("updrs3_score_on ~ NLR_between_z + NLR_within_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Time variable:\n")
cat("year_c is coded in years from baseline; baseline = 0.\n\n")

cat("Outcome:\n")
cat("ON-medication UPDRS-III score.\n\n")

cat("Main predictors:\n")
cat("NLR_between_z = participant mean NLR across available visits, standardized.\n")
cat("NLR_within_z  = visit-wise deviation from participant mean NLR, standardized.\n\n")

cat("Interpretation:\n")
cat("NLR_between_z tests whether individuals with generally higher NLR have higher UPDRS-III across follow-up.\n")
cat("NLR_within_z tests whether visits with higher-than-usual NLR are accompanied by higher UPDRS-III.\n\n")

cat("Detected covariate columns:\n")
print(detected_columns, row.names = FALSE)
cat("\n\n")

cat("Duplicate audit before model cleaning:\n")
print(duplicate_audit_summary, row.names = FALSE)
cat("\n\n")

cat("Participant-visit cleaning summary:\n")
print(cleaning_summary, row.names = FALSE)
cat("\n\n")

cat("NLR between/within decomposition audit:\n")
print(nlr_decomposition_audit, row.names = FALSE)
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
print(VarCorr(lmm_main), comp = c("Variance", "Std.Dev."))
cat("\n\n")

cat("Model fit indices, main model:\n")
print(model_fit_main, row.names = FALSE)
cat("\n\n")

cat("Fixed effects, main model:\n")
print(fixef_main, row.names = FALSE)
cat("\n\n")

cat("Type III ANOVA, main model:\n")
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
cat("Manual linear combinations of fixed effects: NLR_within_z + NLR_within_z:year_c * year_c.\n")
print(within_slopes_by_time_df, row.names = FALSE)
cat("\n\n")

cat("Sensitivity model: baseline outcome-adjusted follow-up-only model\n")
cat("Formula:\n")
cat("updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c + baseline_updrs3_score_on + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
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

cat("Sensitivity model: log-NLR between/within decomposition\n")
cat("Formula:\n")
cat("updrs3_score_on ~ log_NLR_between_z + log_NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")
cat("Rows =", nrow(dat_log), "\n")
cat("Unique PD subjects =", n_distinct(dat_log$PATNO), "\n")
cat("Unique sites =", n_distinct(dat_log$SITE), "\n\n")
cat("Log-NLR decomposition audit:\n")
print(log_nlr_decomposition_audit, row.names = FALSE)
cat("\n\n")
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

cat("Prediction reference values for manual predicted trajectories:\n")
print(prediction_reference, row.names = FALSE)
cat("\n\n")

cat("Model assumption checks saved as diagnostic plots:\n")
cat("- Diagnostic_residuals_vs_fitted_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_QQ_residuals_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_observed_vs_fitted_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR_UPDRS.png\n")
cat("- Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR_UPDRS.png\n")
cat("- Performance_check_model_timevarying_NLR_UPDRS.png\n\n")

cat("Figures saved:\n")
cat("- Descriptive_within_person_NLR_vs_UPDRS.png/pdf\n")
cat("- Figure_timevarying_withinNLR_UPDRS_predicted.png/pdf\n\n")

cat("Excel workbook:\n")
cat(file.path(out_dir, "Timevarying_NLR_UPDRS_outputs.xlsx"), "\n")

sink()

# ---------------------------------------------------------
# Final message
# ---------------------------------------------------------
cat("\n============================================================\n")
cat("Time-varying NLR -> ON-medication UPDRS-III analysis tamamlandı.\n")
cat("All outputs were saved to:\n")
cat(out_dir, "\n\n")
cat("Main summary file:\n")
cat(file.path(out_dir, "LMM_timevarying_NLR_UPDRS_summary.txt"), "\n\n")
cat("Excel workbook:\n")
cat(file.path(out_dir, "Timevarying_NLR_UPDRS_outputs.xlsx"), "\n")
cat("============================================================\n")