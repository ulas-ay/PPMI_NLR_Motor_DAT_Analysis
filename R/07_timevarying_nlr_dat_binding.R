# =========================================================
# 07_timevarying_nlr_dat_binding.R
#
# Time-varying neutrophil-to-lymphocyte ratio (NLR) and
# longitudinal striatal dopamine transporter (DAT) binding
# in Parkinson's disease.
#
# Regional DAT outcomes:
#   - Right putamen
#   - Left putamen
#   - Right caudate
#   - Left caudate
#
# DAT imaging visits:
#   Baseline, year 1, year 2, and year 4
#   (EVENT_ID: BL, V04, V06, V10)
#
# NLR components:
#   NLR_between_z
#       Between-person component: participant-specific mean NLR
#       across available DAT imaging visits, standardized.
#
#   NLR_within_z
#       Within-person component: visit-wise deviation from the
#       participant-specific mean NLR, standardized.
#
# Primary regional model:
#   DAT outcome ~ NLR_between_z + NLR_within_z + year_c +
#       age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
#       (1 | SITE) + (1 | PATNO)
#
# Secondary regional model:
#   DAT outcome ~ NLR_between_z + NLR_within_z * year_c +
#       age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
#       (1 | SITE) + (1 | PATNO)
#
# Sensitivity analyses:
#   1) Follow-up-only model adjusted for baseline DAT in the
#      corresponding striatal region.
#   2) Log-transformed NLR between-/within-person decomposition.
#
# Predicted trajectories and simple slopes are calculated manually
# from fixed-effect estimates and the model covariance matrix; this
# script intentionally does not rely on emmeans() or emtrends().
#
# Default input:
#   <PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
#
# Output:
#   <PPMI_OUTPUT_ROOT>/07_TIMEVARYING_NLR_DAT
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

out_dir <- file.path(output_root, "07_TIMEVARYING_NLR_DAT")

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

dat_put_r_col <- find_col(df_raw, c("MIA_PUTAMEN_R", "mia_putamen_r"))
dat_put_l_col <- find_col(df_raw, c("MIA_PUTAMEN_L", "mia_putamen_l"))
dat_cau_r_col <- find_col(df_raw, c("MIA_CAUDATE_R", "mia_caudate_r"))
dat_cau_l_col <- find_col(df_raw, c("MIA_CAUDATE_L", "mia_caudate_l"))

detected_columns <- tibble(
  variable = c(
    "age",
    "sex",
    "bmi",
    "SITE",
    "duration_yrs",
    "LEDD",
    "DOMSIDE",
    "DAT_putamen_R",
    "DAT_putamen_L",
    "DAT_caudate_R",
    "DAT_caudate_L"
  ),
  detected_column = c(
    age_col,
    sex_col,
    bmi_col,
    site_col,
    duration_col,
    ledd_col,
    domside_col,
    dat_put_r_col,
    dat_put_l_col,
    dat_cau_r_col,
    dat_cau_l_col
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
  ledd_col,
  domside_col,
  dat_put_r_col,
  dat_put_l_col,
  dat_cau_r_col,
  dat_cau_l_col
)))) {
  stop(
    "Could not detect one or more required columns. Check names(df_raw) and update candidate names.",
    call. = FALSE
  )
}

# ---------------------------------------------------------
# Prepare raw PD-only DAT dataset
# DAT imaging visits: BL, V04, V06, V10
# ---------------------------------------------------------
visit_order <- c("BL", "V04", "V06", "V10")

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
dat_all0 <- dat_visit %>%
  mutate(
    PATNO_model = factor(PATNO),
    SITE_model = factor(.data[[site_col]]),
    EVENT_ID_model = factor(EVENT_ID, levels = visit_order),

    year_model = year_from_baseline,
    year_c_model = as.numeric(year_from_baseline),

    DAT_putamen_R_model = as.numeric(.data[[dat_put_r_col]]),
    DAT_putamen_L_model = as.numeric(.data[[dat_put_l_col]]),
    DAT_caudate_R_model = as.numeric(.data[[dat_cau_r_col]]),
    DAT_caudate_L_model = as.numeric(.data[[dat_cau_l_col]]),

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
    DAT_putamen_R = DAT_putamen_R_model,
    DAT_putamen_L = DAT_putamen_L_model,
    DAT_caudate_R = DAT_caudate_R_model,
    DAT_caudate_L = DAT_caudate_L_model,
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

dat_all0 <- droplevels(dat_all0)

write.csv(
  dat_all0,
  file.path(out_dir, "PD_timevarying_NLR_DAT_dataset_before_outcome_specific_filter.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# NLR between/within decomposition
# ---------------------------------------------------------
person_nlr <- dat_all0 %>%
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

dat_all0 <- dat_all0 %>%
  left_join(person_nlr, by = "PATNO") %>%
  mutate(
    NLR_within_raw = NLR - NLR_between_raw
  )

within_mean <- mean(dat_all0$NLR_within_raw, na.rm = TRUE)
within_sd <- sd(dat_all0$NLR_within_raw, na.rm = TRUE)

dat_all0 <- dat_all0 %>%
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
  n_missing_NLR_between_z = sum(is.na(dat_all0$NLR_between_z)),
  pct_missing_NLR_between_z = 100 * mean(is.na(dat_all0$NLR_between_z)),
  n_missing_NLR_within_z = sum(is.na(dat_all0$NLR_within_z)),
  pct_missing_NLR_within_z = 100 * mean(is.na(dat_all0$NLR_within_z))
)

write.csv(
  person_nlr,
  file.path(out_dir, "Participant_mean_NLR_between_component_DAT_visits.csv"),
  row.names = FALSE
)

write.csv(
  nlr_decomposition_audit,
  file.path(out_dir, "NLR_between_within_decomposition_audit_DAT.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Baseline DAT extraction and safe merge
# ---------------------------------------------------------
baseline_dat <- dat_all0 %>%
  filter(EVENT_ID == "BL") %>%
  select(
    PATNO,
    baseline_DAT_putamen_R = DAT_putamen_R,
    baseline_DAT_putamen_L = DAT_putamen_L,
    baseline_DAT_caudate_R = DAT_caudate_R,
    baseline_DAT_caudate_L = DAT_caudate_L
  ) %>%
  distinct(PATNO, .keep_all = TRUE)

baseline_dat_audit <- tibble(
  baseline_dat_rows = nrow(baseline_dat),
  baseline_dat_unique_PATNO = n_distinct(baseline_dat$PATNO),
  n_missing_baseline_DAT_putamen_R = sum(is.na(baseline_dat$baseline_DAT_putamen_R)),
  pct_missing_baseline_DAT_putamen_R = 100 * mean(is.na(baseline_dat$baseline_DAT_putamen_R)),
  n_missing_baseline_DAT_putamen_L = sum(is.na(baseline_dat$baseline_DAT_putamen_L)),
  pct_missing_baseline_DAT_putamen_L = 100 * mean(is.na(baseline_dat$baseline_DAT_putamen_L)),
  n_missing_baseline_DAT_caudate_R = sum(is.na(baseline_dat$baseline_DAT_caudate_R)),
  pct_missing_baseline_DAT_caudate_R = 100 * mean(is.na(baseline_dat$baseline_DAT_caudate_R)),
  n_missing_baseline_DAT_caudate_L = sum(is.na(baseline_dat$baseline_DAT_caudate_L)),
  pct_missing_baseline_DAT_caudate_L = 100 * mean(is.na(baseline_dat$baseline_DAT_caudate_L))
)

rows_before_baseline_dat_join <- nrow(dat_all0)

dat_all0 <- dat_all0 %>%
  left_join(baseline_dat, by = "PATNO")

rows_after_baseline_dat_join <- nrow(dat_all0)

baseline_dat_join_audit <- tibble(
  rows_before_baseline_dat_join = rows_before_baseline_dat_join,
  rows_after_baseline_dat_join = rows_after_baseline_dat_join,
  rows_added_by_join = rows_after_baseline_dat_join - rows_before_baseline_dat_join,
  baseline_dat_rows = nrow(baseline_dat),
  baseline_dat_unique_PATNO = n_distinct(baseline_dat$PATNO)
)

write.csv(
  baseline_dat_audit,
  file.path(out_dir, "Baseline_DAT_audit.csv"),
  row.names = FALSE
)

write.csv(
  baseline_dat_join_audit,
  file.path(out_dir, "Baseline_DAT_join_audit.csv"),
  row.names = FALSE
)

if (rows_after_baseline_dat_join != rows_before_baseline_dat_join) {
  stop(
    "Baseline DAT join multiplied rows. Check baseline_dat uniqueness.",
    call. = FALSE
  )
}

# ---------------------------------------------------------
# Log-NLR between/within decomposition
# ---------------------------------------------------------
person_log_nlr <- dat_all0 %>%
  group_by(PATNO) %>%
  summarise(
    log_NLR_between_raw = mean(log_NLR, na.rm = TRUE),
    n_log_NLR_available = sum(!is.na(log_NLR)),
    .groups = "drop"
  ) %>%
  mutate(
    log_NLR_between_raw = ifelse(is.nan(log_NLR_between_raw), NA_real_, log_NLR_between_raw)
  )

dat_all0 <- dat_all0 %>%
  left_join(person_log_nlr, by = "PATNO") %>%
  mutate(
    log_NLR_within_raw = log_NLR - log_NLR_between_raw
  )

log_between_mean <- mean(person_log_nlr$log_NLR_between_raw, na.rm = TRUE)
log_between_sd <- sd(person_log_nlr$log_NLR_between_raw, na.rm = TRUE)

log_within_mean <- mean(dat_all0$log_NLR_within_raw, na.rm = TRUE)
log_within_sd <- sd(dat_all0$log_NLR_within_raw, na.rm = TRUE)

dat_all0 <- dat_all0 %>%
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
  n_missing_log_NLR_between_z = sum(is.na(dat_all0$log_NLR_between_z)),
  pct_missing_log_NLR_between_z = 100 * mean(is.na(dat_all0$log_NLR_between_z)),
  n_missing_log_NLR_within_z = sum(is.na(dat_all0$log_NLR_within_z)),
  pct_missing_log_NLR_within_z = 100 * mean(is.na(dat_all0$log_NLR_within_z))
)

write.csv(
  person_log_nlr,
  file.path(out_dir, "Participant_mean_logNLR_between_component_DAT_visits.csv"),
  row.names = FALSE
)

write.csv(
  log_nlr_decomposition_audit,
  file.path(out_dir, "LogNLR_between_within_decomposition_audit_DAT.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# Shared settings
# ---------------------------------------------------------
# Sum-to-zero contrasts are used for Type III tests.
options(contrasts = c("contr.sum", "contr.poly"))
time_points <- c(0, 1, 2, 4)

# ---------------------------------------------------------
# Main function to run one time-varying NLR -> DAT model
# ---------------------------------------------------------
run_timevarying_dat_model <- function(
    data,
    outcome_var,
    baseline_outcome_var,
    outcome_label,
    file_stub
) {

  cat("\n============================================================\n")
  cat("Running time-varying NLR DAT model:", file_stub, "\n")
  cat("Outcome:", outcome_var, "\n")
  cat("============================================================\n")

  # -----------------------------------------------------
  # Missingness before complete-case filtering
  # -----------------------------------------------------
  model_vars_main <- c(
    outcome_var,
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

  missingness_by_model_variable <- make_missingness_table(
    data,
    model_vars_main
  )

  missingness_by_predictor_covariate <- make_missingness_table(
    data,
    predictor_covariate_vars_main
  )

  data_missing_flags <- data %>%
    ungroup() %>%
    mutate(
      row_has_missing_model_variable =
        if_any(all_of(model_vars_main), is.na),
      row_has_missing_predictor_covariate =
        if_any(all_of(predictor_covariate_vars_main), is.na)
    )

  missingness_by_visit <- data_missing_flags %>%
    group_by(EVENT_ID, year_c) %>%
    summarise(
      n_rows = n(),
      n_subjects = n_distinct(PATNO),

      n_outcome_missing = sum(is.na(.data[[outcome_var]])),
      pct_outcome_missing = 100 * mean(is.na(.data[[outcome_var]])),

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

  participant_missing_model_vars <- make_participant_missingness(
    data = data,
    vars = model_vars_main,
    label = paste0(file_stub, ": model variables including DAT outcome")
  )

  participant_missing_predictors <- make_participant_missingness(
    data = data,
    vars = predictor_covariate_vars_main,
    label = paste0(file_stub, ": predictors/covariates only; DAT outcome excluded")
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
    file.path(out_dir, paste0(file_stub, "_Missingness_by_model_variable_timevarying_NLR_including_outcome.csv")),
    row.names = FALSE
  )

  write.csv(
    missingness_by_predictor_covariate,
    file.path(out_dir, paste0(file_stub, "_Missingness_by_predictor_covariate_timevarying_NLR.csv")),
    row.names = FALSE
  )

  write.csv(
    missingness_by_visit,
    file.path(out_dir, paste0(file_stub, "_Missingness_by_visit_timevarying_NLR.csv")),
    row.names = FALSE
  )

  write.csv(
    participant_level_missingness,
    file.path(out_dir, paste0(file_stub, "_Participant_level_missingness_timevarying_NLR_long.csv")),
    row.names = FALSE
  )

  write.csv(
    participant_missingness_summary,
    file.path(out_dir, paste0(file_stub, "_Participant_level_missingness_summary_timevarying_NLR.csv")),
    row.names = FALSE
  )

  # -----------------------------------------------------
  # Complete-case dataset
  # -----------------------------------------------------
  dat <- data %>%
    drop_na(all_of(model_vars_main)) %>%
    droplevels()

  dat <- droplevels(dat)

  if (nrow(dat) < 10) {
    stop(paste0("Not enough complete-case observations for ", file_stub), call. = FALSE)
  }

  write.csv(
    dat,
    file.path(out_dir, paste0(file_stub, "_timevarying_NLR_dataset_complete_case.csv")),
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
      mean_DAT = mean(.data[[outcome_var]], na.rm = TRUE),
      sd_DAT = sd(.data[[outcome_var]], na.rm = TRUE),
      median_DAT = median(.data[[outcome_var]], na.rm = TRUE),
      IQR_DAT = IQR(.data[[outcome_var]], na.rm = TRUE),
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
    file.path(out_dir, paste0(file_stub, "_timevarying_NLR_descriptives.csv")),
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
    file.path(out_dir, paste0(file_stub, "_timevarying_NLR_visit_counts.csv")),
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
    file.path(out_dir, paste0(file_stub, "_timevarying_NLR_site_counts.csv")),
    row.names = FALSE
  )

  # -----------------------------------------------------
  # Main model
  # -----------------------------------------------------
  form_main <- as.formula(
    paste0(
      outcome_var,
      " ~ NLR_between_z + NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )

  lmm_main <- safe_lmer(
    formula = form_main,
    data = dat,
    model_name = paste0(file_stub, " main time-varying NLR DAT model")
  )

  fixef_main <- make_fixef_table(lmm_main)
  anova_main <- make_anova_table(lmm_main)
  randef_main <- as.data.frame(VarCorr(lmm_main))

  model_fit_main <- tibble(
    outcome = outcome_var,
    model = "main_timevarying_between_within_NLR",
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
    fixef_main,
    file.path(out_dir, paste0(file_stub, "_LMM_fixed_effects_timevarying_NLR_main.csv")),
    row.names = FALSE
  )

  write.csv(
    anova_main,
    file.path(out_dir, paste0(file_stub, "_LMM_TypeIII_ANOVA_timevarying_NLR_main.csv")),
    row.names = FALSE
  )

  write.csv(
    randef_main,
    file.path(out_dir, paste0(file_stub, "_LMM_random_effects_variance_timevarying_NLR_main.csv")),
    row.names = FALSE
  )

  write.csv(
    model_fit_main,
    file.path(out_dir, paste0(file_stub, "_LMM_model_fit_indices_timevarying_NLR_main.csv")),
    row.names = FALSE
  )

  # -----------------------------------------------------
  # Secondary model: within-person NLR x time
  # -----------------------------------------------------
  form_within_interaction <- as.formula(
    paste0(
      outcome_var,
      " ~ NLR_between_z + NLR_within_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )

  lmm_within_interaction <- safe_lmer(
    formula = form_within_interaction,
    data = dat,
    model_name = paste0(file_stub, " secondary within-person NLR by time DAT model")
  )

  fixef_within_interaction <- make_fixef_table(lmm_within_interaction)
  anova_within_interaction <- make_anova_table(lmm_within_interaction)
  randef_within_interaction <- as.data.frame(VarCorr(lmm_within_interaction))

  model_fit_within_interaction <- tibble(
    outcome = outcome_var,
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
    fixef_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_LMM_fixed_effects.csv")),
    row.names = FALSE
  )

  write.csv(
    anova_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )

  write.csv(
    randef_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_random_effects.csv")),
    row.names = FALSE
  )

  write.csv(
    model_fit_within_interaction,
    file.path(out_dir, paste0(file_stub, "_Secondary_withinNLRxTime_model_fit_indices.csv")),
    row.names = FALSE
  )

  # -----------------------------------------------------
  # Manual simple slopes for within-person NLR effect by time
  # -----------------------------------------------------
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
    file.path(out_dir, paste0(file_stub, "_Within_person_NLR_effect_on_DAT_by_timepoint.csv")),
    row.names = FALSE
  )

  # -----------------------------------------------------
  # Sensitivity model 1:
  # Baseline DAT-adjusted follow-up-only model
  # -----------------------------------------------------
  model_vars_baseline_adj <- c(
    outcome_var,
    baseline_outcome_var,
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
    baseline_outcome_var,
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

  dat_baseline_adj0 <- data %>%
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
    drop_na(all_of(model_vars_baseline_adj)) %>%
    droplevels()

  form_baseline_adj <- as.formula(
    paste0(
      outcome_var,
      " ~ NLR_between_z + NLR_within_z + year_c + ",
      baseline_outcome_var,
      " + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )

  lmm_baseline_adj <- safe_lmer(
    formula = form_baseline_adj,
    data = dat_baseline_adj,
    model_name = paste0(file_stub, " baseline DAT-adjusted follow-up model")
  )

  fixef_baseline_adj <- make_fixef_table(lmm_baseline_adj)
  anova_baseline_adj <- make_anova_table(lmm_baseline_adj)
  randef_baseline_adj <- as.data.frame(VarCorr(lmm_baseline_adj))

  model_fit_baseline_adj <- tibble(
    outcome = outcome_var,
    model = "sensitivity_baseline_DAT_adjusted_followup_only",
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
    dat_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_dataset.csv")),
    row.names = FALSE
  )

  write.csv(
    missingness_by_model_variable_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_missingness_by_model_variable_including_outcome.csv")),
    row.names = FALSE
  )

  write.csv(
    missingness_by_predictor_covariate_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_missingness_by_predictor_covariate.csv")),
    row.names = FALSE
  )

  write.csv(
    fixef_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_LMM_fixed_effects.csv")),
    row.names = FALSE
  )

  write.csv(
    anova_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )

  write.csv(
    randef_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_random_effects.csv")),
    row.names = FALSE
  )

  write.csv(
    model_fit_baseline_adj,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_baseline_DAT_adjusted_model_fit_indices.csv")),
    row.names = FALSE
  )

  # -----------------------------------------------------
  # Sensitivity model 2:
  # Log-NLR between/within decomposition
  # -----------------------------------------------------
  model_vars_log <- c(
    outcome_var,
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
    data,
    model_vars_log
  )

  missingness_by_predictor_covariate_log <- make_missingness_table(
    data,
    predictor_covariate_vars_log
  )

  dat_log <- data %>%
    drop_na(all_of(model_vars_log)) %>%
    droplevels()

  form_log <- as.formula(
    paste0(
      outcome_var,
      " ~ log_NLR_between_z + log_NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + ",
      "(1 | SITE) + (1 | PATNO)"
    )
  )

  lmm_log <- safe_lmer(
    formula = form_log,
    data = dat_log,
    model_name = paste0(file_stub, " log-NLR between-within DAT model")
  )

  fixef_log <- make_fixef_table(lmm_log)
  anova_log <- make_anova_table(lmm_log)

  model_fit_log <- tibble(
    outcome = outcome_var,
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
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_dataset.csv")),
    row.names = FALSE
  )

  write.csv(
    missingness_by_model_variable_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_missingness_by_model_variable_including_outcome.csv")),
    row.names = FALSE
  )

  write.csv(
    missingness_by_predictor_covariate_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_missingness_by_predictor_covariate.csv")),
    row.names = FALSE
  )

  write.csv(
    fixef_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_LMM_fixed_effects.csv")),
    row.names = FALSE
  )

  write.csv(
    anova_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )

  write.csv(
    model_fit_log,
    file.path(out_dir, paste0(file_stub, "_Sensitivity_logNLR_between_within_model_fit_indices.csv")),
    row.names = FALSE
  )

  # -----------------------------------------------------
  # Diagnostics for main model
  # -----------------------------------------------------
  diagnostic_df <- dat %>%
    mutate(
      fitted_value = fitted(lmm_main),
      residual = resid(lmm_main),
      pearson_residual = residual / sigma(lmm_main)
    )

  write.csv(
    diagnostic_df,
    file.path(out_dir, paste0(file_stub, "_Diagnostic_values_timevarying_NLR.csv")),
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
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_residuals_vs_fitted_timevarying_NLR.png")),
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
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_residuals_timevarying_NLR.png")),
    plot = p_qq_resid,
    width = 7,
    height = 5,
    dpi = 600
  )

  p_obs_fit <- ggplot(
    diagnostic_df,
    aes(x = fitted_value, y = .data[[outcome_var]])
  ) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.7) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(file_stub, ": observed vs fitted values"),
      x = "Fitted values",
      y = paste0("Observed ", outcome_label)
    )

  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_observed_vs_fitted_timevarying_NLR.png")),
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
      title = paste0(file_stub, ": Q-Q plot of participant random intercepts"),
      x = "Theoretical quantiles",
      y = "Participant random intercepts"
    )

  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR.png")),
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
      title = paste0(file_stub, ": Q-Q plot of site random intercepts"),
      x = "Theoretical quantiles",
      y = "Site random intercepts"
    )

  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR.png")),
    plot = p_qq_ranef_site,
    width = 7,
    height = 5,
    dpi = 600
  )

  safe_check_model_png(
    model = lmm_main,
    filename = paste0(file_stub, "_Performance_check_model_timevarying_NLR.png")
  )

  # -----------------------------------------------------
  # Manual prediction figure:
  # Predicted DAT at lower-than-usual, usual, higher-than-usual NLR
  # -----------------------------------------------------
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
    file.path(out_dir, paste0(file_stub, "_Predicted_DAT_by_within_person_NLR_levels.csv")),
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
      breaks = c(0, 1, 2, 4),
      labels = c("BL", "1", "2", "4")
    ) +
    labs(
      title = paste0("Time-varying NLR and ", outcome_label),
      subtitle = "Predicted DAT binding at lower-than-usual, usual, and higher-than-usual NLR",
      x = "Years from baseline",
      y = outcome_label,
      color = "Visit-wise NLR",
      fill = "Visit-wise NLR"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure_timevarying_withinNLR_DAT_predicted.png")),
    plot = p_within_pred,
    width = 8.5,
    height = 5.8,
    dpi = 600
  )

  ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure_timevarying_withinNLR_DAT_predicted.pdf")),
    plot = p_within_pred,
    width = 8.5,
    height = 5.8
  )

  # -----------------------------------------------------
  # Summary text
  # -----------------------------------------------------
  sink(file.path(out_dir, paste0(file_stub, "_LMM_timevarying_NLR_DAT_summary.txt")))

  cat("============================================================\n")
  cat("Aim 3B / Reviewer 1 Comment 2:\n")
  cat("Does time-varying NLR track longitudinal DAT binding?\n")
  cat("============================================================\n\n")

  cat("Outcome:\n")
  cat(outcome_var, "(", outcome_label, ")\n\n")

  cat("Input file:\n")
  cat(file_path, "\n\n")

  cat("Output directory:\n")
  cat(out_dir, "\n\n")

  cat("Main model formula:\n")
  cat(paste0(
    outcome_var,
    " ~ NLR_between_z + NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))

  cat("Secondary model formula:\n")
  cat(paste0(
    outcome_var,
    " ~ NLR_between_z + NLR_within_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))

  cat("Time variable:\n")
  cat("year_c is coded in years from baseline; baseline = 0.\n\n")

  cat("Main predictors:\n")
  cat("NLR_between_z = participant mean NLR across available DAT imaging visits, standardized.\n")
  cat("NLR_within_z  = visit-wise deviation from participant mean NLR, standardized.\n\n")

  cat("Interpretation:\n")
  cat("NLR_between_z tests whether individuals with generally higher NLR have different DAT binding across follow-up.\n")
  cat("NLR_within_z tests whether visits with higher-than-usual NLR are accompanied by different DAT binding.\n\n")

  cat("Covariates:\n")
  cat("age, sex, BMI, disease duration, LEDD, DOMSIDE.\n\n")

  cat("Duplicate audit before model cleaning:\n")
  print(duplicate_audit_summary, row.names = FALSE)
  cat("\n\n")

  cat("Participant-visit cleaning summary:\n")
  print(cleaning_summary, row.names = FALSE)
  cat("\n\n")

  cat("Sample size before outcome-specific complete-case filtering:\n")
  cat("Rows =", nrow(data), "\n")
  cat("Unique PD subjects =", n_distinct(data$PATNO), "\n")
  cat("Unique sites =", n_distinct(data$SITE), "\n")
  cat("Unique participant-visits =", n_distinct(paste(data$PATNO, data$EVENT_ID)), "\n\n")

  cat("Sample size in main complete-case model:\n")
  cat("Rows =", nrow(dat), "\n")
  cat("Unique PD subjects =", n_distinct(dat$PATNO), "\n")
  cat("Unique sites =", n_distinct(dat$SITE), "\n")
  cat("Unique participant-visits =", n_distinct(paste(dat$PATNO, dat$EVENT_ID)), "\n\n")

  cat("Visit counts:\n")
  print(table(dat$EVENT_ID))
  cat("\n\n")

  cat("Missingness by model variable, including DAT outcome:\n")
  print(missingness_by_model_variable, row.names = FALSE)
  cat("\n\n")

  cat("Missingness by predictors/covariates only, DAT outcome excluded:\n")
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

  cat("Main model: time-varying NLR between/within components\n")
  cat("Random-effects variance:\n")
  print(VarCorr(lmm_main), comp = c("Variance", "Std.Dev."))
  cat("\n\n")

  cat("Model fit indices:\n")
  print(model_fit_main, row.names = FALSE)
  cat("\n\n")

  cat("Fixed effects:\n")
  print(fixef_main, row.names = FALSE)
  cat("\n\n")

  cat("Type III ANOVA:\n")
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

  cat("Sensitivity model: baseline DAT-adjusted follow-up-only model\n")
  cat("Formula:\n")
  cat(paste0(
    outcome_var,
    " ~ NLR_between_z + NLR_within_z + year_c + ",
    baseline_outcome_var,
    " + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))
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
  cat(paste0(
    outcome_var,
    " ~ log_NLR_between_z + log_NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
  ))
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

  cat("Prediction reference values:\n")
  print(prediction_reference, row.names = FALSE)
  cat("\n\n")

  cat("Model assumption checks saved as diagnostic plots:\n")
  cat(paste0("- ", file_stub, "_Diagnostic_residuals_vs_fitted_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_residuals_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_observed_vs_fitted_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_PATNO_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Diagnostic_QQ_random_intercepts_SITE_timevarying_NLR.png\n"))
  cat(paste0("- ", file_stub, "_Performance_check_model_timevarying_NLR.png\n\n"))

  cat("Figure saved:\n")
  cat(paste0("- ", file_stub, "_Figure_timevarying_withinNLR_DAT_predicted.png/pdf\n"))

  sink()

  return(
    list(
      outcome = outcome_var,
      file_stub = file_stub,

      model_main = lmm_main,
      fixef_main = fixef_main,
      anova_main = anova_main,
      randef_main = randef_main,
      model_fit_main = model_fit_main,

      model_within_interaction = lmm_within_interaction,
      fixef_within_interaction = fixef_within_interaction,
      anova_within_interaction = anova_within_interaction,
      randef_within_interaction = randef_within_interaction,
      model_fit_within_interaction = model_fit_within_interaction,
      within_slopes = within_slopes_by_time_df,

      sensitivity_baseline_adj_fixef = fixef_baseline_adj,
      sensitivity_baseline_adj_anova = anova_baseline_adj,
      sensitivity_baseline_adj_fit = model_fit_baseline_adj,

      sensitivity_log_fixef = fixef_log,
      sensitivity_log_anova = anova_log,
      sensitivity_log_fit = model_fit_log,

      descriptives = desc,
      visit_counts = visit_counts,
      site_counts = site_counts,

      missingness_model = missingness_by_model_variable,
      missingness_predictor = missingness_by_predictor_covariate,
      missingness_visit = missingness_by_visit,
      participant_missingness_summary = participant_missingness_summary,

      prediction_within = emm_within_df
    )
  )
}

# ---------------------------------------------------------
# Run models
# ---------------------------------------------------------
res_put_r <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_putamen_R",
  baseline_outcome_var = "baseline_DAT_putamen_R",
  outcome_label = "Putaminal DAT binding (right)",
  file_stub = "PUTAMEN_R"
)

res_put_l <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_putamen_L",
  baseline_outcome_var = "baseline_DAT_putamen_L",
  outcome_label = "Putaminal DAT binding (left)",
  file_stub = "PUTAMEN_L"
)

res_cau_r <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_caudate_R",
  baseline_outcome_var = "baseline_DAT_caudate_R",
  outcome_label = "Caudate DAT binding (right)",
  file_stub = "CAUDATE_R"
)

res_cau_l <- run_timevarying_dat_model(
  data = dat_all0,
  outcome_var = "DAT_caudate_L",
  baseline_outcome_var = "baseline_DAT_caudate_L",
  outcome_label = "Caudate DAT binding (left)",
  file_stub = "CAUDATE_L"
)

# ---------------------------------------------------------
# Combined tables across outcomes
# ---------------------------------------------------------
combined_fixef_main <- bind_rows(
  res_put_r$fixef_main %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$fixef_main %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$fixef_main %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$fixef_main %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_anova_main <- bind_rows(
  res_put_r$anova_main %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$anova_main %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$anova_main %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$anova_main %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_model_fit_main <- bind_rows(
  res_put_r$model_fit_main,
  res_put_l$model_fit_main,
  res_cau_r$model_fit_main,
  res_cau_l$model_fit_main
)

combined_within_interaction_fixef <- bind_rows(
  res_put_r$fixef_within_interaction %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$fixef_within_interaction %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$fixef_within_interaction %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$fixef_within_interaction %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_within_interaction_anova <- bind_rows(
  res_put_r$anova_within_interaction %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$anova_within_interaction %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$anova_within_interaction %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$anova_within_interaction %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_model_fit_within_interaction <- bind_rows(
  res_put_r$model_fit_within_interaction,
  res_put_l$model_fit_within_interaction,
  res_cau_r$model_fit_within_interaction,
  res_cau_l$model_fit_within_interaction
)

combined_within_slopes <- bind_rows(
  res_put_r$within_slopes %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$within_slopes %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$within_slopes %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$within_slopes %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_sensitivity_baseline_adj_fixef <- bind_rows(
  res_put_r$sensitivity_baseline_adj_fixef %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$sensitivity_baseline_adj_fixef %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$sensitivity_baseline_adj_fixef %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$sensitivity_baseline_adj_fixef %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_sensitivity_baseline_adj_anova <- bind_rows(
  res_put_r$sensitivity_baseline_adj_anova %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$sensitivity_baseline_adj_anova %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$sensitivity_baseline_adj_anova %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$sensitivity_baseline_adj_anova %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_sensitivity_baseline_adj_fit <- bind_rows(
  res_put_r$sensitivity_baseline_adj_fit,
  res_put_l$sensitivity_baseline_adj_fit,
  res_cau_r$sensitivity_baseline_adj_fit,
  res_cau_l$sensitivity_baseline_adj_fit
)

combined_sensitivity_log_fixef <- bind_rows(
  res_put_r$sensitivity_log_fixef %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$sensitivity_log_fixef %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$sensitivity_log_fixef %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$sensitivity_log_fixef %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_sensitivity_log_anova <- bind_rows(
  res_put_r$sensitivity_log_anova %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$sensitivity_log_anova %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$sensitivity_log_anova %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$sensitivity_log_anova %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_sensitivity_log_fit <- bind_rows(
  res_put_r$sensitivity_log_fit,
  res_put_l$sensitivity_log_fit,
  res_cau_r$sensitivity_log_fit,
  res_cau_l$sensitivity_log_fit
)

combined_missingness_predictor <- bind_rows(
  res_put_r$missingness_predictor %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$missingness_predictor %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$missingness_predictor %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$missingness_predictor %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_missingness_model <- bind_rows(
  res_put_r$missingness_model %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$missingness_model %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$missingness_model %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$missingness_model %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_participant_missingness <- bind_rows(
  res_put_r$participant_missingness_summary %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$participant_missingness_summary %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$participant_missingness_summary %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$participant_missingness_summary %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

combined_prediction_within <- bind_rows(
  res_put_r$prediction_within %>% mutate(outcome = "PUTAMEN_R"),
  res_put_l$prediction_within %>% mutate(outcome = "PUTAMEN_L"),
  res_cau_r$prediction_within %>% mutate(outcome = "CAUDATE_R"),
  res_cau_l$prediction_within %>% mutate(outcome = "CAUDATE_L")
) %>%
  select(outcome, everything())

# ---------------------------------------------------------
# Save combined CSVs
# ---------------------------------------------------------
write.csv(
  combined_fixef_main,
  file.path(out_dir, "COMBINED_timevarying_NLR_DAT_main_LMM_fixed_effects_all_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_anova_main,
  file.path(out_dir, "COMBINED_timevarying_NLR_DAT_main_LMM_TypeIII_ANOVA_all_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_model_fit_main,
  file.path(out_dir, "COMBINED_timevarying_NLR_DAT_main_model_fit_indices_all_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_within_interaction_fixef,
  file.path(out_dir, "COMBINED_withinNLRxTime_LMM_fixed_effects_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_within_interaction_anova,
  file.path(out_dir, "COMBINED_withinNLRxTime_LMM_TypeIII_ANOVA_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_model_fit_within_interaction,
  file.path(out_dir, "COMBINED_withinNLRxTime_model_fit_indices_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_within_slopes,
  file.path(out_dir, "COMBINED_within_person_NLR_effect_on_DAT_by_timepoint_all_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_baseline_adj_fixef,
  file.path(out_dir, "COMBINED_sensitivity_baseline_DAT_adjusted_fixed_effects_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_baseline_adj_anova,
  file.path(out_dir, "COMBINED_sensitivity_baseline_DAT_adjusted_ANOVA_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_baseline_adj_fit,
  file.path(out_dir, "COMBINED_sensitivity_baseline_DAT_adjusted_model_fit_indices_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_log_fixef,
  file.path(out_dir, "COMBINED_sensitivity_logNLR_between_within_fixed_effects_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_log_anova,
  file.path(out_dir, "COMBINED_sensitivity_logNLR_between_within_ANOVA_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_sensitivity_log_fit,
  file.path(out_dir, "COMBINED_sensitivity_logNLR_between_within_model_fit_indices_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_missingness_predictor,
  file.path(out_dir, "COMBINED_missingness_by_predictor_covariate_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_missingness_model,
  file.path(out_dir, "COMBINED_missingness_by_model_variable_including_outcome_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_participant_missingness,
  file.path(out_dir, "COMBINED_participant_level_missingness_summary_all_DAT_outcomes.csv"),
  row.names = FALSE
)

write.csv(
  combined_prediction_within,
  file.path(out_dir, "COMBINED_predicted_DAT_by_within_person_NLR_levels_all_outcomes.csv"),
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
add_sheet(wb, "NLR_decomp_audit", nlr_decomposition_audit)
add_sheet(wb, "LogNLR_decomp_audit", log_nlr_decomposition_audit)
add_sheet(wb, "Baseline_DAT_audit", baseline_dat_audit)
add_sheet(wb, "Baseline_DAT_join", baseline_dat_join_audit)

add_sheet(wb, "Combined_main_fixef", combined_fixef_main)
add_sheet(wb, "Combined_main_ANOVA", combined_anova_main)
add_sheet(wb, "Combined_main_fit", combined_model_fit_main)

add_sheet(wb, "Combined_within_fixef", combined_within_interaction_fixef)
add_sheet(wb, "Combined_within_ANOVA", combined_within_interaction_anova)
add_sheet(wb, "Combined_within_fit", combined_model_fit_within_interaction)
add_sheet(wb, "Combined_within_slopes", combined_within_slopes)

add_sheet(wb, "Combined_BLadj_fixef", combined_sensitivity_baseline_adj_fixef)
add_sheet(wb, "Combined_BLadj_ANOVA", combined_sensitivity_baseline_adj_anova)
add_sheet(wb, "Combined_BLadj_fit", combined_sensitivity_baseline_adj_fit)

add_sheet(wb, "Combined_log_fixef", combined_sensitivity_log_fixef)
add_sheet(wb, "Combined_log_ANOVA", combined_sensitivity_log_anova)
add_sheet(wb, "Combined_log_fit", combined_sensitivity_log_fit)

add_sheet(wb, "Combined_miss_pred", combined_missingness_predictor)
add_sheet(wb, "Combined_miss_model", combined_missingness_model)
add_sheet(wb, "Combined_part_miss", combined_participant_missingness)
add_sheet(wb, "Combined_prediction", combined_prediction_within)

add_sheet(wb, "PUT_R_desc", res_put_r$descriptives)
add_sheet(wb, "PUT_L_desc", res_put_l$descriptives)
add_sheet(wb, "CAU_R_desc", res_cau_r$descriptives)
add_sheet(wb, "CAU_L_desc", res_cau_l$descriptives)

openxlsx::saveWorkbook(
  wb,
  file = file.path(out_dir, "Timevarying_NLR_DAT_outputs.xlsx"),
  overwrite = TRUE
)

# ---------------------------------------------------------
# Overall summary text
# ---------------------------------------------------------
sink(file.path(out_dir, "COMBINED_timevarying_NLR_DAT_summary.txt"))

cat("============================================================\n")
cat("Aim 3B / Reviewer 1 Comment 2:\n")
cat("Does time-varying NLR track longitudinal DAT binding?\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(file_path, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("DAT outcomes:\n")
cat("- DAT_putamen_R\n")
cat("- DAT_putamen_L\n")
cat("- DAT_caudate_R\n")
cat("- DAT_caudate_L\n\n")

cat("DAT visits included:\n")
cat("BL, V04, V06, V10\n\n")

cat("Main model formula for each outcome:\n")
cat("DAT outcome ~ NLR_between_z + NLR_within_z + year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Secondary model formula for each outcome:\n")
cat("DAT outcome ~ NLR_between_z + NLR_within_z * year_c + age + sex + bmi + duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n")

cat("Sensitivity model 1:\n")
cat("Follow-up DAT only, additionally adjusted for baseline DAT of the same region.\n\n")

cat("Sensitivity model 2:\n")
cat("Log-NLR between/within decomposition model.\n\n")

cat("Duplicate audit before model cleaning:\n")
print(duplicate_audit_summary, row.names = FALSE)
cat("\n\n")

cat("Participant-visit cleaning summary:\n")
print(cleaning_summary, row.names = FALSE)
cat("\n\n")

cat("NLR between/within decomposition audit:\n")
print(nlr_decomposition_audit, row.names = FALSE)
cat("\n\n")

cat("Log-NLR between/within decomposition audit:\n")
print(log_nlr_decomposition_audit, row.names = FALSE)
cat("\n\n")

cat("Baseline DAT audit:\n")
print(baseline_dat_audit, row.names = FALSE)
cat("\n\n")

cat("Baseline DAT join audit:\n")
print(baseline_dat_join_audit, row.names = FALSE)
cat("\n\n")

cat("Detected columns:\n")
print(detected_columns, row.names = FALSE)
cat("\n\n")

cat("Combined fixed effects, main models:\n")
print(combined_fixef_main, row.names = FALSE)
cat("\n\n")

cat("Combined Type III ANOVA, main models:\n")
print(combined_anova_main, row.names = FALSE)
cat("\n\n")

cat("Combined model fit indices, main models:\n")
print(combined_model_fit_main, row.names = FALSE)
cat("\n\n")

cat("Combined within-person NLR x time fixed effects:\n")
print(combined_within_interaction_fixef, row.names = FALSE)
cat("\n\n")

cat("Combined within-person NLR x time ANOVA:\n")
print(combined_within_interaction_anova, row.names = FALSE)
cat("\n\n")

cat("Combined within-person NLR simple slopes by time point:\n")
print(combined_within_slopes, row.names = FALSE)
cat("\n\n")

cat("Combined predictor/covariate missingness:\n")
print(combined_missingness_predictor, row.names = FALSE)
cat("\n\n")

cat("Combined model-variable missingness including DAT outcome:\n")
print(combined_missingness_model, row.names = FALSE)
cat("\n\n")

cat("Combined participant-level missingness summary:\n")
print(combined_participant_missingness, row.names = FALSE)
cat("\n\n")

cat("Combined sensitivity: baseline DAT-adjusted fixed effects:\n")
print(combined_sensitivity_baseline_adj_fixef, row.names = FALSE)
cat("\n\n")

cat("Combined sensitivity: baseline DAT-adjusted ANOVA:\n")
print(combined_sensitivity_baseline_adj_anova, row.names = FALSE)
cat("\n\n")

cat("Combined sensitivity: log-NLR between/within fixed effects:\n")
print(combined_sensitivity_log_fixef, row.names = FALSE)
cat("\n\n")

cat("Combined sensitivity: log-NLR between/within ANOVA:\n")
print(combined_sensitivity_log_anova, row.names = FALSE)
cat("\n\n")

cat("Excel workbook:\n")
cat(file.path(out_dir, "Timevarying_NLR_DAT_outputs.xlsx"), "\n\n")

sink()

# ---------------------------------------------------------
# Final message
# ---------------------------------------------------------
cat("\n============================================================\n")
cat("Time-varying NLR -> longitudinal DAT binding analyses completed.\n")
cat("All outputs were saved to:\n")
cat(out_dir, "\n\n")
cat("Main combined summary:\n")
cat(file.path(out_dir, "COMBINED_timevarying_NLR_DAT_summary.txt"), "\n\n")
cat("Excel workbook:\n")
cat(file.path(out_dir, "Timevarying_NLR_DAT_outputs.xlsx"), "\n")
cat("============================================================\n")