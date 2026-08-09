# =========================================================
# 08_baseline_dat_moderation_of_nlr_motor_progression.R
#
# Moderation of the association between baseline
# neutrophil-to-lymphocyte ratio (NLR) and longitudinal motor
# progression by baseline striatal dopamine transporter (DAT)
# binding in Parkinson's disease.
#
# Outcome:
#   ON-medication MDS-UPDRS Part III / UPDRS-III
#
# Predictor:
#   baseline_NLR_z
#
# Baseline DAT moderators:
#   - Left caudate
#   - Right caudate
#   - Left putamen
#   - Right putamen
#
# Primary moderation model:
#   UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c +
#       age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
#       (1 | SITE) + (1 | PATNO)
#
# Sensitivity analysis:
#   Follow-up-only model additionally adjusted for baseline UPDRS-III.
#
# Simple slopes and predicted trajectories are calculated manually
# from fixed-effect estimates and the model covariance matrix; this
# script intentionally does not rely on emmeans() or emtrends().
#
# Default input:
#   <PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
#
# Output:
#   <PPMI_OUTPUT_ROOT>/08_NLR_MODERATION_STRIATAL_DAT
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

out_dir <- file.path(output_root, "08_NLR_MODERATION_STRIATAL_DAT")

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

get_interaction_term_by_components <- function(term_names, components) {
    
    split_terms <- strsplit(term_names, ":", fixed = TRUE)
    
    hit <- term_names[
        vapply(
            split_terms,
            function(x) {
                identical(sort(x), sort(components))
            },
            logical(1)
        )
    ]
    
    if (length(hit) == 0) {
        stop(
            paste0(
                "Could not find interaction term for components:\n",
                paste(components, collapse = " : "),
                "\n\nAvailable fixed-effect terms:\n",
                paste(term_names, collapse = "\n")
            ),
            call. = FALSE
        )
    }
    
    hit[1]
}

manual_nlr_slopes_by_dat_and_time <- function(
        model,
        dat_values,
        time_values
) {
    
    beta <- lme4::fixef(model)
    vc <- as.matrix(vcov(model))
    term_names <- names(beta)
    
    term_nlr <- "baseline_NLR_z"
    
    term_nlr_dat <- get_interaction_term_by_components(
        term_names = term_names,
        components = c("baseline_NLR_z", "DAT_BL_z")
    )
    
    term_nlr_time <- get_interaction_term_by_components(
        term_names = term_names,
        components = c("baseline_NLR_z", "year_c")
    )
    
    term_nlr_dat_time <- get_interaction_term_by_components(
        term_names = term_names,
        components = c("baseline_NLR_z", "DAT_BL_z", "year_c")
    )
    
    if (!term_nlr %in% term_names) {
        stop("baseline_NLR_z fixed-effect term not found.", call. = FALSE)
    }
    
    out <- expand.grid(
        DAT_BL_z = dat_values,
        year_c = time_values,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    ) %>%
        as_tibble()
    
    bind_rows(
        lapply(
            seq_len(nrow(out)),
            function(i) {
                
                d <- out$DAT_BL_z[i]
                t <- out$year_c[i]
                
                L <- rep(0, length(beta))
                names(L) <- names(beta)
                
                L[term_nlr] <- 1
                L[term_nlr_dat] <- d
                L[term_nlr_time] <- t
                L[term_nlr_dat_time] <- d * t
                
                estimate <- sum(L * beta)
                SE <- sqrt(as.numeric(t(L) %*% vc %*% L))
                z_value <- estimate / SE
                p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
                
                tibble(
                    DAT_BL_z = d,
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
    )
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

add_sheet <- function(wb, sheet_name, data) {
    sheet_name <- substr(sheet_name, 1, 31)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = data)
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
# Detect columns
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

dat_caudate_l_col <- find_col(df_raw, c("MIA_CAUDATE_L", "mia_caudate_l"))
dat_caudate_r_col <- find_col(df_raw, c("MIA_CAUDATE_R", "mia_caudate_r"))
dat_putamen_l_col <- find_col(df_raw, c("MIA_PUTAMEN_L", "mia_putamen_l"))
dat_putamen_r_col <- find_col(df_raw, c("MIA_PUTAMEN_R", "mia_putamen_r"))

detected_columns <- tibble(
    variable = c(
        "age",
        "sex",
        "bmi",
        "SITE",
        "duration_yrs",
        "UPDRS",
        "LEDD",
        "DOMSIDE",
        "DAT_caudate_L",
        "DAT_caudate_R",
        "DAT_putamen_L",
        "DAT_putamen_R"
    ),
    detected_column = c(
        age_col,
        sex_col,
        bmi_col,
        site_col,
        duration_col,
        updrs_col,
        ledd_col,
        domside_col,
        dat_caudate_l_col,
        dat_caudate_r_col,
        dat_putamen_l_col,
        dat_putamen_r_col
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

if (any(is.na(c(
    age_col,
    sex_col,
    bmi_col,
    site_col,
    duration_col,
    updrs_col,
    ledd_col,
    domside_col,
    dat_caudate_l_col,
    dat_caudate_r_col,
    dat_putamen_l_col,
    dat_putamen_r_col
)))) {
    stop(
        "Could not detect one or more required columns. Check names(df_raw) and update candidate names.",
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Prepare raw PD-only longitudinal UPDRS dataset
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
# Prepare explicit longitudinal modelling dataset
# ---------------------------------------------------------
dat_long0 <- dat_visit %>%
    mutate(
        PATNO_model = factor(PATNO),
        SITE_model = factor(.data[[site_col]]),
        EVENT_ID_model = factor(EVENT_ID, levels = visit_order),
        
        year_model = year_from_baseline,
        year_c_model = as.numeric(year_from_baseline),
        
        UPDRS_model = as.numeric(.data[[updrs_col]]),
        
        baseline_NLR_model = as.numeric(baseline_NLR),
        baseline_NLR_z_model = as.numeric(baseline_NLR_z),
        log_baseline_NLR_model = ifelse(
            !is.na(baseline_NLR_model) & baseline_NLR_model > 0,
            log(baseline_NLR_model),
            NA_real_
        ),
        
        DAT_caudate_L_model = as.numeric(.data[[dat_caudate_l_col]]),
        DAT_caudate_R_model = as.numeric(.data[[dat_caudate_r_col]]),
        DAT_putamen_L_model = as.numeric(.data[[dat_putamen_l_col]]),
        DAT_putamen_R_model = as.numeric(.data[[dat_putamen_r_col]]),
        
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
        UPDRS = UPDRS_model,
        baseline_NLR = baseline_NLR_model,
        baseline_NLR_z = baseline_NLR_z_model,
        log_baseline_NLR = log_baseline_NLR_model,
        DAT_caudate_L = DAT_caudate_L_model,
        DAT_caudate_R = DAT_caudate_R_model,
        DAT_putamen_L = DAT_putamen_L_model,
        DAT_putamen_R = DAT_putamen_R_model,
        age = age_model,
        sex = sex_model,
        bmi = bmi_model,
        duration_yrs = duration_yrs_model,
        LEDD = LEDD_model,
        LEDD_raw = LEDD_raw_model,
        DOMSIDE = DOMSIDE_model
    )

dat_long0 <- droplevels(dat_long0)

write.csv(
    dat_long0,
    file.path(out_dir, "Shared_longitudinal_UPDRS_dataset_before_moderator_filter.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Baseline UPDRS extraction and safe merge
# ---------------------------------------------------------
baseline_updrs <- dat_long0 %>%
    filter(EVENT_ID == "BL") %>%
    select(
        PATNO,
        baseline_UPDRS = UPDRS
    ) %>%
    distinct(PATNO, .keep_all = TRUE)

baseline_updrs_audit <- tibble(
    baseline_updrs_rows = nrow(baseline_updrs),
    baseline_updrs_unique_PATNO = n_distinct(baseline_updrs$PATNO),
    n_missing_baseline_UPDRS = sum(is.na(baseline_updrs$baseline_UPDRS)),
    pct_missing_baseline_UPDRS =
        100 * mean(is.na(baseline_updrs$baseline_UPDRS))
)

rows_before_baseline_updrs_join <- nrow(dat_long0)

dat_long0 <- dat_long0 %>%
    left_join(baseline_updrs, by = "PATNO")

rows_after_baseline_updrs_join <- nrow(dat_long0)

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
# Shared settings
# ---------------------------------------------------------
# Sum-to-zero contrasts are used for Type III tests.
options(contrasts = c("contr.sum", "contr.poly"))

# ---------------------------------------------------------
# Main moderation function
# ---------------------------------------------------------
run_striatal_moderation <- function(
        data,
        dat_var,
        dat_label,
        file_stub
) {
    
    cat("\n============================================================\n")
    cat("Running moderation model:", file_stub, "\n")
    cat("DAT moderator:", dat_var, "\n")
    cat("============================================================\n")
    
    # -----------------------------------------------------
    # Baseline DAT extraction
    # -----------------------------------------------------
    bl_dat_raw <- data %>%
        filter(EVENT_ID == "BL") %>%
        select(
            PATNO,
            DAT_BL = all_of(dat_var)
        ) %>%
        distinct(PATNO, .keep_all = TRUE)
    
    dat_bl_mean <- mean(bl_dat_raw$DAT_BL, na.rm = TRUE)
    dat_bl_sd <- sd(bl_dat_raw$DAT_BL, na.rm = TRUE)
    
    bl_dat <- bl_dat_raw %>%
        mutate(
            DAT_BL_z = as.numeric((DAT_BL - dat_bl_mean) / dat_bl_sd)
        )
    
    baseline_dat_audit <- tibble(
        moderator = dat_var,
        dat_label = dat_label,
        baseline_DAT_rows = nrow(bl_dat),
        baseline_DAT_unique_PATNO = n_distinct(bl_dat$PATNO),
        n_missing_DAT_BL = sum(is.na(bl_dat$DAT_BL)),
        pct_missing_DAT_BL = 100 * mean(is.na(bl_dat$DAT_BL)),
        DAT_BL_mean_for_standardization = dat_bl_mean,
        DAT_BL_sd_for_standardization = dat_bl_sd,
        n_missing_DAT_BL_z = sum(is.na(bl_dat$DAT_BL_z)),
        pct_missing_DAT_BL_z = 100 * mean(is.na(bl_dat$DAT_BL_z))
    )
    
    write.csv(
        bl_dat,
        file.path(out_dir, paste0(file_stub, "_baseline_DAT.csv")),
        row.names = FALSE
    )
    
    write.csv(
        baseline_dat_audit,
        file.path(out_dir, paste0(file_stub, "_baseline_DAT_audit.csv")),
        row.names = FALSE
    )
    
    # -----------------------------------------------------
    # Merge baseline DAT into longitudinal UPDRS data
    # -----------------------------------------------------
    rows_before_dat_join <- nrow(data)
    
    dat0 <- data %>%
        left_join(bl_dat, by = "PATNO")
    
    rows_after_dat_join <- nrow(dat0)
    
    baseline_dat_join_audit <- tibble(
        moderator = dat_var,
        rows_before_DAT_join = rows_before_dat_join,
        rows_after_DAT_join = rows_after_dat_join,
        rows_added_by_join = rows_after_dat_join - rows_before_dat_join,
        baseline_DAT_rows = nrow(bl_dat),
        baseline_DAT_unique_PATNO = n_distinct(bl_dat$PATNO)
    )
    
    write.csv(
        baseline_dat_join_audit,
        file.path(out_dir, paste0(file_stub, "_baseline_DAT_join_audit.csv")),
        row.names = FALSE
    )
    
    if (rows_after_dat_join != rows_before_dat_join) {
        stop(
            paste0(
                file_stub,
                ": Baseline DAT join multiplied rows. Check baseline DAT uniqueness."
            ),
            call. = FALSE
        )
    }
    
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
    
    predictor_covariate_vars_main <- c(
        "baseline_NLR_z",
        "DAT_BL_z",
        "year_c",
        "age",
        "sex",
        "bmi",
        "duration_yrs",
        "LEDD",
        "DOMSIDE",
        "SITE"
    )
    
    missingness_by_model_variable <- make_missingness_table(
        dat0,
        model_vars_main
    )
    
    missingness_by_predictor_covariate <- make_missingness_table(
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
    
    missingness_by_visit <- dat0_missing_flags %>%
        group_by(EVENT_ID, year_c) %>%
        summarise(
            n_rows = n(),
            n_subjects = n_distinct(PATNO),
            
            n_UPDRS_missing = sum(is.na(UPDRS)),
            pct_UPDRS_missing = 100 * mean(is.na(UPDRS)),
            
            n_baseline_NLR_z_missing = sum(is.na(baseline_NLR_z)),
            pct_baseline_NLR_z_missing = 100 * mean(is.na(baseline_NLR_z)),
            
            n_DAT_BL_z_missing = sum(is.na(DAT_BL_z)),
            pct_DAT_BL_z_missing = 100 * mean(is.na(DAT_BL_z)),
            
            n_model_variable_missing_rows = sum(row_has_missing_model_variable),
            pct_model_variable_missing_rows = 100 * mean(row_has_missing_model_variable),
            
            n_predictor_or_covariate_missing_rows = sum(row_has_missing_predictor_covariate),
            pct_predictor_or_covariate_missing_rows = 100 * mean(row_has_missing_predictor_covariate),
            
            .groups = "drop"
        ) %>%
        arrange(year_c)
    
    participant_missing_model_vars <- make_participant_missingness(
        data = dat0,
        vars = model_vars_main,
        label = paste0(file_stub, ": model variables including outcome UPDRS-III")
    )
    
    participant_missing_predictors <- make_participant_missingness(
        data = dat0,
        vars = predictor_covariate_vars_main,
        label = paste0(file_stub, ": predictors/covariates only; outcome UPDRS-III excluded")
    )
    
    participant_level_missingness <- bind_rows(
        participant_missing_model_vars$participant_level,
        participant_missing_predictors$participant_level
    )
    
    participant_missingness_summary <- bind_rows(
        participant_missing_model_vars$summary,
        participant_missing_predictors$summary
    )
    
    write.csv(
        missingness_by_model_variable,
        file.path(out_dir, paste0(file_stub, "_Missingness_by_model_variable_including_outcome.csv")),
        row.names = FALSE
    )
    
    write.csv(
        missingness_by_predictor_covariate,
        file.path(out_dir, paste0(file_stub, "_Missingness_by_predictor_covariate.csv")),
        row.names = FALSE
    )
    
    write.csv(
        missingness_by_visit,
        file.path(out_dir, paste0(file_stub, "_Missingness_by_visit.csv")),
        row.names = FALSE
    )
    
    write.csv(
        participant_level_missingness,
        file.path(out_dir, paste0(file_stub, "_Participant_level_missingness_long.csv")),
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
        drop_na(all_of(model_vars_main)) %>%
        droplevels()
    
    if (nrow(dat) < 10) {
        stop(
            paste0("Not enough complete-case observations for ", file_stub),
            call. = FALSE
        )
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
        group_by(EVENT_ID, year_c) %>%
        summarise(
            n_rows = n(),
            n_subjects = n_distinct(PATNO),
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
        arrange(year_c)
    
    visit_counts <- dat %>%
        group_by(EVENT_ID, year_c) %>%
        summarise(
            n_rows = n(),
            n_subjects = n_distinct(PATNO),
            .groups = "drop"
        ) %>%
        arrange(year_c)
    
    site_counts <- dat %>%
        group_by(SITE) %>%
        summarise(
            n_rows = n(),
            n_subjects = n_distinct(PATNO),
            .groups = "drop"
        ) %>%
        arrange(desc(n_rows))
    
    write.csv(
        desc,
        file.path(out_dir, paste0(file_stub, "_descriptives.csv")),
        row.names = FALSE
    )
    
    write.csv(
        visit_counts,
        file.path(out_dir, paste0(file_stub, "_visit_counts.csv")),
        row.names = FALSE
    )
    
    write.csv(
        site_counts,
        file.path(out_dir, paste0(file_stub, "_site_counts.csv")),
        row.names = FALSE
    )
    
    # -----------------------------------------------------
    # Main model
    # -----------------------------------------------------
    lmm <- safe_lmer(
        UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c +
            age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
            (1 | SITE) + (1 | PATNO),
        data = dat,
        model_name = paste0(file_stub, " main NLR x DAT x time moderation model")
    )
    
    fixef_tab <- make_fixef_table(lmm)
    anova_tab <- make_anova_table(lmm)
    randef_var <- as.data.frame(VarCorr(lmm))
    
    model_fit <- tibble(
        moderator = dat_var,
        dat_label = dat_label,
        model = "main_NLR_DAT_time_moderation",
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
    # Manual simple slopes:
    # effect of baseline NLR at each time point,
    # stratified by baseline DAT level
    # -----------------------------------------------------
    time_points <- sort(unique(dat$year_c))
    dat_levels <- c(-1.5, 0, 1.5)
    
    slopes_nlr_by_dat_time_df <- manual_nlr_slopes_by_dat_and_time(
        model = lmm,
        dat_values = dat_levels,
        time_values = time_points
    ) %>%
        mutate(
            DAT_group = factor(
                DAT_BL_z,
                levels = dat_levels,
                labels = c(
                    "Low baseline DAT (-1.5 SD)",
                    "Mean baseline DAT",
                    "High baseline DAT (+1.5 SD)"
                )
            )
        )
    
    write.csv(
        slopes_nlr_by_dat_time_df,
        file.path(out_dir, paste0(file_stub, "_simple_slopes_NLR_by_DAT_and_time.csv")),
        row.names = FALSE
    )
    
    # -----------------------------------------------------
    # Manual predicted trajectories:
    # NLR = -1.5, 0, +1.5 SD
    # DAT = -1.5, 0, +1.5 SD
    # -----------------------------------------------------
    nlr_levels <- c(-1.5, 0, 1.5)
    
    prediction_reference <- tibble(
        age = mean(dat$age, na.rm = TRUE),
        sex = factor(levels(dat$sex)[1], levels = levels(dat$sex)),
        bmi = mean(dat$bmi, na.rm = TRUE),
        duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
        LEDD = mean(dat$LEDD, na.rm = TRUE),
        DOMSIDE = factor(levels(dat$DOMSIDE)[1], levels = levels(dat$DOMSIDE))
    )
    
    pred_grid <- expand.grid(
        year_c = time_points,
        baseline_NLR_z = nlr_levels,
        DAT_BL_z = dat_levels,
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
    
    emm_df <- manual_predictions_from_lmer(
        model = lmm,
        newdata = pred_grid
    ) %>%
        mutate(
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
            )
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
        filter(EVENT_ID == "BL") %>%
        select(PATNO, baseline_NLR_z, DAT_BL) %>%
        distinct(PATNO, .keep_all = TRUE) %>%
        filter(!is.na(baseline_NLR_z), !is.na(DAT_BL)) %>%
        mutate(
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
        left_join(
            dat_plot %>% select(PATNO, NLR_tertile, DAT_tertile),
            by = "PATNO"
        ) %>%
        filter(!is.na(NLR_tertile), !is.na(DAT_tertile)) %>%
        group_by(DAT_tertile, NLR_tertile, year_c, EVENT_ID) %>%
        summarise(
            mean_UPDRS = mean(UPDRS, na.rm = TRUE),
            se_UPDRS = sd(UPDRS, na.rm = TRUE) / sqrt(n()),
            n = n(),
            n_subjects = n_distinct(PATNO),
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
    
    predictor_covariate_vars_sens <- c(
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
        "SITE"
    )
    
    dat_sens0 <- dat0 %>%
        filter(year_c > 0)
    
    missingness_by_model_variable_sens <- make_missingness_table(
        dat_sens0,
        model_vars_sens
    )
    
    missingness_by_predictor_covariate_sens <- make_missingness_table(
        dat_sens0,
        predictor_covariate_vars_sens
    )
    
    dat_sens <- dat_sens0 %>%
        drop_na(all_of(model_vars_sens)) %>%
        droplevels()
    
    lmm_sens <- safe_lmer(
        UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c +
            baseline_UPDRS +
            age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
            (1 | SITE) + (1 | PATNO),
        data = dat_sens,
        model_name = paste0(file_stub, " sensitivity baseline UPDRS-adjusted moderation model")
    )
    
    fixef_sens <- make_fixef_table(lmm_sens)
    anova_sens <- make_anova_table(lmm_sens)
    randef_sens <- as.data.frame(VarCorr(lmm_sens))
    
    model_fit_sens <- tibble(
        moderator = dat_var,
        dat_label = dat_label,
        model = "sensitivity_baseline_UPDRS_adjusted_followup_only",
        AIC = AIC(lmm_sens),
        BIC = BIC(lmm_sens),
        logLik = as.numeric(logLik(lmm_sens)),
        deviance = deviance(lmm_sens),
        sigma = sigma(lmm_sens),
        n_obs = nobs(lmm_sens),
        n_subjects = n_distinct(dat_sens$PATNO),
        n_sites = n_distinct(dat_sens$SITE)
    )
    
    write.csv(
        dat_sens,
        file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_dataset.csv")),
        row.names = FALSE
    )
    
    write.csv(
        missingness_by_model_variable_sens,
        file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_missingness_by_model_variable_including_outcome.csv")),
        row.names = FALSE
    )
    
    write.csv(
        missingness_by_predictor_covariate_sens,
        file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_UPDRS_adjusted_missingness_by_predictor_covariate.csv")),
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
    # Diagnostics for main model
    # -----------------------------------------------------
    diagnostic_df <- dat %>%
        mutate(
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
        rename(random_intercept = `(Intercept)`)
    
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
        rename(random_intercept = `(Intercept)`)
    
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
    
    safe_check_model_png(
        model = lmm,
        filename = paste0(file_stub, "_Performance_check_model.png")
    )
    
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
    cat("DAT_BL_z; baseline DAT binding standardized across baseline participants with available DAT for this region.\n\n")
    
    cat("Time variable:\n")
    cat("year_c is coded in years from baseline; baseline = 0.\n\n")
    
    cat("Covariates:\n")
    cat("age, sex, BMI, disease duration, LEDD, DOMSIDE.\n\n")
    
    cat("Baseline DAT audit:\n")
    print(baseline_dat_audit, row.names = FALSE)
    cat("\n\n")
    
    cat("Baseline DAT join audit:\n")
    print(baseline_dat_join_audit, row.names = FALSE)
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
    
    cat("Visit counts:\n")
    print(table(dat$EVENT_ID))
    cat("\n\n")
    
    cat("Missingness by model variable, including outcome UPDRS-III:\n")
    print(missingness_by_model_variable, row.names = FALSE)
    cat("\n\n")
    
    cat("Missingness by predictors/covariates only, excluding outcome UPDRS-III:\n")
    print(missingness_by_predictor_covariate, row.names = FALSE)
    cat("\n\n")
    
    cat("Missingness by visit:\n")
    print(missingness_by_visit, row.names = FALSE)
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
    cat("Manual linear combinations of fixed effects:\n")
    cat("NLR effect = b_NLR + b_NLR:DAT * DAT + b_NLR:time * time + b_NLR:DAT:time * DAT * time\n\n")
    print(slopes_nlr_by_dat_time_df, row.names = FALSE)
    cat("\n\n")
    
    cat("Prediction reference values:\n")
    print(prediction_reference, row.names = FALSE)
    cat("\n\n")
    
    cat("Sensitivity model: baseline UPDRS-adjusted follow-up-only model\n")
    cat("Rows =", nrow(dat_sens), "\n")
    cat("Unique PD subjects =", n_distinct(dat_sens$PATNO), "\n")
    cat("Unique sites =", n_distinct(dat_sens$SITE), "\n\n")
    
    cat("Sensitivity missingness by model variable, including outcome:\n")
    print(missingness_by_model_variable_sens, row.names = FALSE)
    cat("\n\n")
    
    cat("Sensitivity missingness by predictors/covariates only:\n")
    print(missingness_by_predictor_covariate_sens, row.names = FALSE)
    cat("\n\n")
    
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
    
    return(
        list(
            moderator = dat_var,
            dat_label = dat_label,
            
            model = lmm,
            fixef = fixef_tab,
            anova = anova_tab,
            randef = randef_var,
            model_fit = model_fit,
            
            slopes = slopes_nlr_by_dat_time_df,
            predictions = emm_df,
            
            sensitivity_model = lmm_sens,
            sensitivity_fixef = fixef_sens,
            sensitivity_anova = anova_sens,
            sensitivity_randef = randef_sens,
            sensitivity_model_fit = model_fit_sens,
            
            descriptives = desc,
            visit_counts = visit_counts,
            site_counts = site_counts,
            
            missingness_model = missingness_by_model_variable,
            missingness_predictor = missingness_by_predictor_covariate,
            missingness_visit = missingness_by_visit,
            participant_missingness_summary = participant_missingness_summary,
            
            sensitivity_missingness_model = missingness_by_model_variable_sens,
            sensitivity_missingness_predictor = missingness_by_predictor_covariate_sens,
            
            baseline_dat_audit = baseline_dat_audit,
            baseline_dat_join_audit = baseline_dat_join_audit
        )
    )
}

# ---------------------------------------------------------
# Run all striatal moderation models
# ---------------------------------------------------------
res_caudate_l <- run_striatal_moderation(
    data = dat_long0,
    dat_var = "DAT_caudate_L",
    dat_label = "Left caudate DAT",
    file_stub = "CAUDATE_L"
)

res_caudate_r <- run_striatal_moderation(
    data = dat_long0,
    dat_var = "DAT_caudate_R",
    dat_label = "Right caudate DAT",
    file_stub = "CAUDATE_R"
)

res_putamen_l <- run_striatal_moderation(
    data = dat_long0,
    dat_var = "DAT_putamen_L",
    dat_label = "Left putaminal DAT",
    file_stub = "PUTAMEN_L"
)

res_putamen_r <- run_striatal_moderation(
    data = dat_long0,
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

combined_model_fit <- bind_rows(
    res_caudate_l$model_fit,
    res_caudate_r$model_fit,
    res_putamen_l$model_fit,
    res_putamen_r$model_fit
)

combined_slopes <- bind_rows(
    res_caudate_l$slopes %>% mutate(moderator = "CAUDATE_L"),
    res_caudate_r$slopes %>% mutate(moderator = "CAUDATE_R"),
    res_putamen_l$slopes %>% mutate(moderator = "PUTAMEN_L"),
    res_putamen_r$slopes %>% mutate(moderator = "PUTAMEN_R")
) %>%
    select(moderator, everything())

combined_predictions <- bind_rows(
    res_caudate_l$predictions %>% mutate(moderator = "CAUDATE_L"),
    res_caudate_r$predictions %>% mutate(moderator = "CAUDATE_R"),
    res_putamen_l$predictions %>% mutate(moderator = "PUTAMEN_L"),
    res_putamen_r$predictions %>% mutate(moderator = "PUTAMEN_R")
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

combined_sensitivity_model_fit <- bind_rows(
    res_caudate_l$sensitivity_model_fit,
    res_caudate_r$sensitivity_model_fit,
    res_putamen_l$sensitivity_model_fit,
    res_putamen_r$sensitivity_model_fit
)

combined_missingness_predictor <- bind_rows(
    res_caudate_l$missingness_predictor %>% mutate(moderator = "CAUDATE_L"),
    res_caudate_r$missingness_predictor %>% mutate(moderator = "CAUDATE_R"),
    res_putamen_l$missingness_predictor %>% mutate(moderator = "PUTAMEN_L"),
    res_putamen_r$missingness_predictor %>% mutate(moderator = "PUTAMEN_R")
) %>%
    select(moderator, everything())

combined_missingness_model <- bind_rows(
    res_caudate_l$missingness_model %>% mutate(moderator = "CAUDATE_L"),
    res_caudate_r$missingness_model %>% mutate(moderator = "CAUDATE_R"),
    res_putamen_l$missingness_model %>% mutate(moderator = "PUTAMEN_L"),
    res_putamen_r$missingness_model %>% mutate(moderator = "PUTAMEN_R")
) %>%
    select(moderator, everything())

combined_participant_missingness <- bind_rows(
    res_caudate_l$participant_missingness_summary %>% mutate(moderator = "CAUDATE_L"),
    res_caudate_r$participant_missingness_summary %>% mutate(moderator = "CAUDATE_R"),
    res_putamen_l$participant_missingness_summary %>% mutate(moderator = "PUTAMEN_L"),
    res_putamen_r$participant_missingness_summary %>% mutate(moderator = "PUTAMEN_R")
) %>%
    select(moderator, everything())

combined_baseline_dat_audit <- bind_rows(
    res_caudate_l$baseline_dat_audit,
    res_caudate_r$baseline_dat_audit,
    res_putamen_l$baseline_dat_audit,
    res_putamen_r$baseline_dat_audit
)

combined_baseline_dat_join_audit <- bind_rows(
    res_caudate_l$baseline_dat_join_audit,
    res_caudate_r$baseline_dat_join_audit,
    res_putamen_l$baseline_dat_join_audit,
    res_putamen_r$baseline_dat_join_audit
)

# ---------------------------------------------------------
# Save combined CSV files
# ---------------------------------------------------------
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
    combined_model_fit,
    file.path(out_dir, "COMBINED_moderation_LMM_model_fit_indices_all_DAT_moderators.csv"),
    row.names = FALSE
)

write.csv(
    combined_slopes,
    file.path(out_dir, "COMBINED_simple_slopes_NLR_by_DAT_and_time_all_moderators.csv"),
    row.names = FALSE
)

write.csv(
    combined_predictions,
    file.path(out_dir, "COMBINED_predicted_UPDRS_trajectories_1p5SD_all_DAT_moderators.csv"),
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

write.csv(
    combined_sensitivity_model_fit,
    file.path(out_dir, "COMBINED_sensitivity_baseline_UPDRS_adjusted_model_fit_indices_all_DAT_moderators.csv"),
    row.names = FALSE
)

write.csv(
    combined_missingness_predictor,
    file.path(out_dir, "COMBINED_missingness_by_predictor_covariate_all_DAT_moderators.csv"),
    row.names = FALSE
)

write.csv(
    combined_missingness_model,
    file.path(out_dir, "COMBINED_missingness_by_model_variable_including_outcome_all_DAT_moderators.csv"),
    row.names = FALSE
)

write.csv(
    combined_participant_missingness,
    file.path(out_dir, "COMBINED_participant_level_missingness_summary_all_DAT_moderators.csv"),
    row.names = FALSE
)

write.csv(
    combined_baseline_dat_audit,
    file.path(out_dir, "COMBINED_baseline_DAT_audit_all_DAT_moderators.csv"),
    row.names = FALSE
)

write.csv(
    combined_baseline_dat_join_audit,
    file.path(out_dir, "COMBINED_baseline_DAT_join_audit_all_DAT_moderators.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Workbook output
# ---------------------------------------------------------
wb <- openxlsx::createWorkbook()

add_sheet(wb, "Detected_columns", detected_columns)
add_sheet(wb, "Duplicate_summary", duplicate_audit_summary)
add_sheet(wb, "Duplicate_by_visit", duplicate_summary_by_visit)
add_sheet(wb, "Duplicate_discrepancy", duplicate_discrepancy_by_variable)
add_sheet(wb, "Cleaning_summary", cleaning_summary)
add_sheet(wb, "Baseline_UPDRS_audit", baseline_updrs_audit)
add_sheet(wb, "Baseline_UPDRS_join", baseline_updrs_join_audit)
add_sheet(wb, "Baseline_DAT_audit", combined_baseline_dat_audit)
add_sheet(wb, "Baseline_DAT_join", combined_baseline_dat_join_audit)

add_sheet(wb, "Combined_fixef", combined_fixef)
add_sheet(wb, "Combined_ANOVA", combined_anova)
add_sheet(wb, "Combined_fit", combined_model_fit)
add_sheet(wb, "Combined_slopes", combined_slopes)
add_sheet(wb, "Combined_predictions", combined_predictions)

add_sheet(wb, "Sens_fixef", combined_sensitivity_fixef)
add_sheet(wb, "Sens_ANOVA", combined_sensitivity_anova)
add_sheet(wb, "Sens_fit", combined_sensitivity_model_fit)

add_sheet(wb, "Missing_predictors", combined_missingness_predictor)
add_sheet(wb, "Missing_model", combined_missingness_model)
add_sheet(wb, "Participant_missing", combined_participant_missingness)

add_sheet(wb, "CAU_L_desc", res_caudate_l$descriptives)
add_sheet(wb, "CAU_R_desc", res_caudate_r$descriptives)
add_sheet(wb, "PUT_L_desc", res_putamen_l$descriptives)
add_sheet(wb, "PUT_R_desc", res_putamen_r$descriptives)

openxlsx::saveWorkbook(
    wb,
    file = file.path(out_dir, "NLR_striatal_DAT_moderation_outputs.xlsx"),
    overwrite = TRUE
)

# ---------------------------------------------------------
# Combined summary text
# ---------------------------------------------------------
sink(file.path(out_dir, "COMBINED_NLR_striatal_DAT_moderation_summary.txt"))

cat("============================================================\n")
cat("Aim 4: Does baseline striatal DAT moderate the association\n")
cat("between baseline NLR and longitudinal motor progression?\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(file_path, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Main model formula for each moderator:\n")
cat("UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Sensitivity model formula for each moderator:\n")
cat("UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c + baseline_UPDRS + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("DAT moderators:\n")
cat("- DAT_caudate_L\n")
cat("- DAT_caudate_R\n")
cat("- DAT_putamen_L\n")
cat("- DAT_putamen_R\n\n")

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

cat("Baseline DAT audit across moderators:\n")
print(combined_baseline_dat_audit, row.names = FALSE)
cat("\n\n")

cat("Baseline DAT join audit across moderators:\n")
print(combined_baseline_dat_join_audit, row.names = FALSE)
cat("\n\n")

cat("Detected columns:\n")
print(detected_columns, row.names = FALSE)
cat("\n\n")

cat("Combined fixed effects, main models:\n")
print(combined_fixef, row.names = FALSE)
cat("\n\n")

cat("Combined Type III ANOVA, main models:\n")
print(combined_anova, row.names = FALSE)
cat("\n\n")

cat("Combined model fit indices, main models:\n")
print(combined_model_fit, row.names = FALSE)
cat("\n\n")

cat("Combined simple slopes of baseline NLR by DAT level and time:\n")
print(combined_slopes, row.names = FALSE)
cat("\n\n")

cat("Combined predictor/covariate missingness:\n")
print(combined_missingness_predictor, row.names = FALSE)
cat("\n\n")

cat("Combined model-variable missingness including UPDRS outcome:\n")
print(combined_missingness_model, row.names = FALSE)
cat("\n\n")

cat("Combined participant-level missingness summary:\n")
print(combined_participant_missingness, row.names = FALSE)
cat("\n\n")

cat("Combined sensitivity fixed effects:\n")
print(combined_sensitivity_fixef, row.names = FALSE)
cat("\n\n")

cat("Combined sensitivity Type III ANOVA:\n")
print(combined_sensitivity_anova, row.names = FALSE)
cat("\n\n")

cat("Excel workbook:\n")
cat(file.path(out_dir, "NLR_striatal_DAT_moderation_outputs.xlsx"), "\n\n")

sink()

# ---------------------------------------------------------
# Final message
# ---------------------------------------------------------
cat("\n============================================================\n")
cat("NLR x baseline striatal DAT moderation analyses completed.\n")
cat("All outputs were saved to:\n")
cat(out_dir, "\n\n")
cat("Main combined summary:\n")
cat(file.path(out_dir, "COMBINED_NLR_striatal_DAT_moderation_summary.txt"), "\n\n")
cat("Excel workbook:\n")
cat(file.path(out_dir, "NLR_striatal_DAT_moderation_outputs.xlsx"), "\n")
cat("============================================================\n")