# =========================================================
# 02_baseline_demographic_characteristics.R
#
# Generate baseline demographic, clinical, imaging, and
# peripheral inflammatory characteristics (Table 1) and
# longitudinal follow-up availability summaries.
#
# This script:
#   - reads the participant-visit-level dataset created by Script 01;
#   - audits PATNO-EVENT_ID duplicates before analysis;
#   - collapses data to one row per participant-visit if needed;
#   - generates baseline PD/HC descriptive and group-comparison outputs;
#   - performs Welch t-tests with Wilcoxon sensitivity analyses;
#   - checks continuous-variable normality and chi-square assumptions;
#   - summarizes PD-specific clinical and regional DAT characteristics;
#   - exports follow-up/attrition tables and figures.
#
# Input/output paths are intentionally not hard-coded.
# Provide them either as command-line arguments:
#   Rscript R/02_baseline_demographic_characteristics.R <nlr_dataset.xlsx> <output_root>
# or through environment variables:
#   PPMI_NLR_DATA_FILE=/path/to/PPMI_with_NLR_all_visits_updated.xlsx
#   PPMI_OUTPUT_ROOT=/path/to/output
#
# If PPMI_NLR_DATA_FILE is not supplied, the script looks for:
#   <PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
#
# The individual-level PPMI data used by this script are not
# distributed with the repository.
# =========================================================

# ---------------------------------------------------------
# Packages
# ---------------------------------------------------------
required_packages <- c(
    "dplyr",
    "readxl",
    "tidyr",
    "tibble",
    "broom",
    "ggplot2",
    "gt",
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
            ".\nInstall them before running this script, for example with:\n",
            "install.packages(c(",
            paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
            "))"
        ),
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(dplyr)
    library(readxl)
    library(tidyr)
    library(tibble)
    library(broom)
    library(ggplot2)
    library(gt)
    library(openxlsx)
})

# ---------------------------------------------------------
# Input file and output directory
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

output_root <- if (length(args) >= 2 && nzchar(args[2])) {
    args[2]
} else {
    Sys.getenv("PPMI_OUTPUT_ROOT", unset = "")
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

output_root <- normalizePath(output_root, mustWork = FALSE)

input_file <- if (length(args) >= 1 && nzchar(args[1])) {
    args[1]
} else {
    Sys.getenv("PPMI_NLR_DATA_FILE", unset = "")
}

if (!nzchar(input_file)) {
    input_file <- file.path(
        output_root,
        "01_NLR",
        "PPMI_with_NLR_all_visits_updated.xlsx"
    )
}

input_file <- normalizePath(input_file, mustWork = FALSE)
out_dir <- file.path(output_root, "02_BASELINE_CHARACTERISTICS")

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

if (!file.exists(input_file)) {
    stop(
        paste0(
            "Input file not found:\n",
            input_file,
            "\n\nRun Script 01 first, provide the NLR dataset as the first ",
            "command-line argument, or set PPMI_NLR_DATA_FILE."
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

fmt_mean_sd <- function(x, digits = 2) {
    x <- as.numeric(x)
    if (all(is.na(x))) return("—")
    sprintf(
        paste0("%.", digits, "f ± %.", digits, "f"),
        mean(x, na.rm = TRUE),
        sd(x, na.rm = TRUE)
    )
}

fmt_median_iqr <- function(x, digits = 2) {
    x <- as.numeric(x)
    if (all(is.na(x))) return("—")
    q <- quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
    sprintf(
        paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
        q[2], q[1], q[3]
    )
}

fmt_n_pct <- function(n, denom, digits = 1) {
    if (is.na(n) || is.na(denom) || denom == 0) return("—")
    sprintf(
        paste0("%d (%.", digits, "f%%)"),
        n,
        100 * n / denom
    )
}

fmt_p <- function(p) {
    ifelse(
        is.na(p),
        "—",
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

safe_welch <- function(data, var) {
    d <- data %>%
        filter(!is.na(.data[[var]]), !is.na(GROUP)) %>%
        mutate(GROUP = factor(GROUP, levels = c("PD", "HC"))) %>%
        select(GROUP, value = all_of(var))
    
    if (n_distinct(d$GROUP) < 2) {
        return(tibble(
            variable = var,
            test = "Welch t-test",
            estimate_PD_minus_HC = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_
        ))
    }
    
    out <- tryCatch(
        broom::tidy(t.test(value ~ GROUP, data = d, var.equal = FALSE)),
        error = function(e) NULL
    )
    
    if (is.null(out)) {
        return(tibble(
            variable = var,
            test = "Welch t-test",
            estimate_PD_minus_HC = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_
        ))
    }
    
    tibble(
        variable = var,
        test = "Welch t-test",
        estimate_PD_minus_HC = out$estimate,
        statistic = out$statistic,
        p_value = out$p.value
    )
}

safe_wilcox <- function(data, var) {
    d <- data %>%
        filter(!is.na(.data[[var]]), !is.na(GROUP)) %>%
        mutate(GROUP = factor(GROUP, levels = c("PD", "HC"))) %>%
        select(GROUP, value = all_of(var))
    
    if (n_distinct(d$GROUP) < 2) {
        return(NA_real_)
    }
    
    tryCatch(
        wilcox.test(value ~ GROUP, data = d, exact = FALSE)$p.value,
        error = function(e) NA_real_
    )
}

normality_checks <- function(data, vars, labels, out_dir) {
    out <- list()
    
    for (v in vars) {
        if (!v %in% names(data)) next
        
        tmp <- data %>%
            filter(!is.na(.data[[v]]), !is.na(GROUP)) %>%
            group_by(GROUP) %>%
            summarise(
                n = n(),
                mean = mean(.data[[v]], na.rm = TRUE),
                sd = sd(.data[[v]], na.rm = TRUE),
                median = median(.data[[v]], na.rm = TRUE),
                IQR = IQR(.data[[v]], na.rm = TRUE),
                shapiro_W = ifelse(
                    n() >= 3 & n() <= 5000,
                    tryCatch(shapiro.test(.data[[v]])$statistic, error = function(e) NA_real_),
                    NA_real_
                ),
                shapiro_p = ifelse(
                    n() >= 3 & n() <= 5000,
                    tryCatch(shapiro.test(.data[[v]])$p.value, error = function(e) NA_real_),
                    NA_real_
                ),
                .groups = "drop"
            ) %>%
            mutate(
                variable = v,
                label = labels[[v]],
                shapiro_p_formatted = fmt_p(shapiro_p)
            ) %>%
            select(variable, label, everything())
        
        out[[v]] <- tmp
        
        p <- ggplot(
            data %>% filter(!is.na(.data[[v]]), !is.na(GROUP)),
            aes(sample = .data[[v]])
        ) +
            stat_qq(alpha = 0.55, size = 1.6) +
            stat_qq_line(linewidth = 0.7) +
            facet_wrap(~ GROUP, scales = "free") +
            theme_classic(base_size = 13) +
            labs(
                title = paste0("Q-Q plot: ", labels[[v]]),
                x = "Theoretical quantiles",
                y = "Sample quantiles"
            )
        
        ggsave(
            filename = file.path(out_dir, paste0("QQ_", v, ".png")),
            plot = p,
            width = 7,
            height = 4.5,
            dpi = 600
        )
    }
    
    bind_rows(out)
}

chisq_check <- function(data, var, label) {
    d <- data %>%
        filter(!is.na(GROUP), !is.na(.data[[var]]))
    
    tab <- table(d$GROUP, d[[var]])
    
    if (nrow(tab) < 2 || ncol(tab) < 2) {
        return(tibble(
            variable = var,
            label = label,
            test_used = "Not tested",
            min_expected = NA_real_,
            pct_expected_less_5 = NA_real_,
            any_expected_less_1 = NA,
            assumption_satisfied = NA,
            chi_square_p = NA_real_,
            fisher_p = NA_real_,
            selected_p = NA_real_,
            selected_p_formatted = "—"
        ))
    }
    
    chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
    expected <- chi$expected
    
    min_expected <- min(expected)
    pct_expected_less_5 <- 100 * mean(expected < 5)
    any_expected_less_1 <- any(expected < 1)
    
    assumption_satisfied <- min_expected >= 1 && pct_expected_less_5 <= 20
    
    fisher_p <- tryCatch(
        fisher.test(tab)$p.value,
        error = function(e) NA_real_
    )
    
    selected_p <- if (assumption_satisfied) chi$p.value else fisher_p
    test_used <- if (assumption_satisfied) "Chi-square" else "Fisher exact"
    
    tibble(
        variable = var,
        label = label,
        test_used = test_used,
        min_expected = min_expected,
        pct_expected_less_5 = pct_expected_less_5,
        any_expected_less_1 = any_expected_less_1,
        assumption_satisfied = assumption_satisfied,
        chi_square_p = chi$p.value,
        fisher_p = fisher_p,
        selected_p = selected_p,
        selected_p_formatted = fmt_p(selected_p)
    )
}

# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------
df_raw <- readxl::read_excel(input_file, sheet = "All_data_with_NLR")

cat("\nData loaded:\n")
cat("Rows:", nrow(df_raw), "\n")
cat("Columns:", ncol(df_raw), "\n\n")

required_core <- c("PATNO", "EVENT_ID", "PRIMDIAG")
missing_core <- setdiff(required_core, names(df_raw))

if (length(missing_core) > 0) {
    stop(
        paste0(
            "Missing required core column(s): ",
            paste(missing_core, collapse = ", ")
        ),
        call. = FALSE
    )
}

# ---------------------------------------------------------
# Detect variables
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

duration_col <- find_col(
    df_raw,
    c("duration_yrs", "disease_duration", "disease_duration_years", "DURATION_YRS")
)

updrs_on_col <- find_col(
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

education_col <- find_col(
    df_raw,
    c("EDUCYRS", "education", "Education", "education_years", "educyrs", "EDUCATION")
)

moca_col <- find_col(
    df_raw,
    c("MOCA", "MoCA", "moca", "moca_score", "MCA")
)

upsit_col <- find_col(
    df_raw,
    c("UPSIT", "upsit", "UPSIT_score", "upsit_score")
)

caudate_l_col <- find_col(df_raw, c("MIA_CAUDATE_L", "mia_caudate_l"))
caudate_r_col <- find_col(df_raw, c("MIA_CAUDATE_R", "mia_caudate_r"))
putamen_l_col <- find_col(df_raw, c("MIA_PUTAMEN_L", "mia_putamen_l"))
putamen_r_col <- find_col(df_raw, c("MIA_PUTAMEN_R", "mia_putamen_r"))

site_col <- find_col(
    df_raw,
    c("SITE", "site", "Site", "siteid", "SITEID", "site_id")
)

detected_columns <- tibble(
    variable = c(
        "age_col",
        "sex_col",
        "bmi_col",
        "duration_col",
        "updrs_on_col",
        "ledd_col",
        "domside_col",
        "education_col",
        "moca_col",
        "upsit_col",
        "site_col",
        "caudate_l_col",
        "caudate_r_col",
        "putamen_l_col",
        "putamen_r_col"
    ),
    detected_column = c(
        age_col,
        sex_col,
        bmi_col,
        duration_col,
        updrs_on_col,
        ledd_col,
        domside_col,
        education_col,
        moca_col,
        upsit_col,
        site_col,
        caudate_l_col,
        caudate_r_col,
        putamen_l_col,
        putamen_r_col
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

required_detected <- c(age_col, sex_col, bmi_col)

if (any(is.na(required_detected))) {
    stop("Could not detect age, sex, or BMI columns. Check names(df_raw).")
}

# ---------------------------------------------------------
# Keep PD and HC
# PRIMDIAG: 1 = PD, 17 = HC
# ---------------------------------------------------------
df2_raw <- df_raw %>%
    filter(PRIMDIAG %in% c(1, 17)) %>%
    mutate(
        GROUP = factor(
            ifelse(PRIMDIAG == 1, "PD", "HC"),
            levels = c("HC", "PD")
        ),
        EVENT_ID = as.character(EVENT_ID)
    )

# ---------------------------------------------------------
# Visit coding
# ---------------------------------------------------------
visit_order <- c("BL", "V04", "V06", "V08", "V10", "V12", "V13", "V14")

visit_labels <- tibble(
    EVENT_ID = visit_order,
    visit_label = c(
        "Baseline",
        "Year 1",
        "Year 2",
        "Year 3",
        "Year 4",
        "Year 5",
        "Year 6",
        "Year 7"
    ),
    year_c = c(0, 1, 2, 3, 4, 5, 6, 7)
)

df2_raw <- df2_raw %>%
    left_join(visit_labels, by = "EVENT_ID")

# ---------------------------------------------------------
# PATNO-EVENT_ID duplicate audit
# ---------------------------------------------------------
duplicate_keys <- df2_raw %>%
    count(PATNO, EVENT_ID, name = "n_rows") %>%
    filter(n_rows > 1) %>%
    arrange(desc(n_rows), PATNO, EVENT_ID)

write.csv(
    duplicate_keys,
    file.path(out_dir, "Duplicate_PATNO_EVENT_keys_before_collapsing.csv"),
    row.names = FALSE
)

duplicate_rows <- df2_raw %>%
    semi_join(duplicate_keys, by = c("PATNO", "EVENT_ID")) %>%
    arrange(PATNO, EVENT_ID)

write.csv(
    duplicate_rows,
    file.path(out_dir, "Duplicate_PATNO_EVENT_rows_before_collapsing.csv"),
    row.names = FALSE
)

if (nrow(duplicate_rows) > 0) {
    
    duplicate_summary_by_group_visit <- duplicate_rows %>%
        count(GROUP, EVENT_ID, visit_label, year_c, PATNO, name = "n_rows") %>%
        group_by(GROUP, EVENT_ID, visit_label, year_c) %>%
        summarise(
            n_duplicate_participant_visits = n(),
            n_rows_in_duplicate_participant_visits = sum(n_rows),
            n_extra_rows_due_to_duplicates = sum(n_rows - 1),
            max_rows_per_participant_visit = max(n_rows),
            .groups = "drop"
        ) %>%
        arrange(GROUP, year_c)
    
} else {
    
    duplicate_summary_by_group_visit <- tibble(
        GROUP = character(),
        EVENT_ID = character(),
        visit_label = character(),
        year_c = numeric(),
        n_duplicate_participant_visits = integer(),
        n_rows_in_duplicate_participant_visits = integer(),
        n_extra_rows_due_to_duplicates = integer(),
        max_rows_per_participant_visit = integer()
    )
}

write.csv(
    duplicate_summary_by_group_visit,
    file.path(out_dir, "Duplicate_PATNO_EVENT_summary_by_group_visit.csv"),
    row.names = FALSE
)

if (nrow(duplicate_rows) > 0) {
    
    variables_to_check <- setdiff(names(df2_raw), c("PATNO", "EVENT_ID"))
    
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
    
    duplicate_discrepancy_by_variable <- tibble(
        variable = character(),
        n_duplicate_participant_visits_checked = integer(),
        n_with_more_than_one_distinct_nonmissing_value = integer(),
        pct_with_more_than_one_distinct_nonmissing_value = numeric()
    )
}

write.csv(
    duplicate_discrepancy_by_variable,
    file.path(out_dir, "Duplicate_PATNO_EVENT_discrepancy_by_variable.csv"),
    row.names = FALSE
)

duplicate_audit_summary <- tibble(
    raw_rows_PD_HC = nrow(df2_raw),
    raw_unique_participants = n_distinct(df2_raw$PATNO),
    raw_unique_participant_visits = n_distinct(paste(df2_raw$PATNO, df2_raw$EVENT_ID)),
    n_duplicate_participant_visit_keys = nrow(duplicate_keys),
    n_rows_in_duplicate_participant_visits = nrow(duplicate_rows),
    n_extra_rows_due_to_duplicates = nrow(df2_raw) -
        n_distinct(paste(df2_raw$PATNO, df2_raw$EVENT_ID))
)

write.csv(
    duplicate_audit_summary,
    file.path(out_dir, "Duplicate_PATNO_EVENT_audit_summary.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Collapse to participant-visit level
# ---------------------------------------------------------
df2 <- collapse_to_participant_visit(df2_raw)

df2 <- df2 %>%
    mutate(
        GROUP = factor(as.character(GROUP), levels = c("HC", "PD")),
        EVENT_ID = as.character(EVENT_ID)
    )

participant_visit_cleaning_summary <- tibble(
    rows_before_collapsing = nrow(df2_raw),
    rows_after_collapsing = nrow(df2),
    rows_removed_by_collapsing = nrow(df2_raw) - nrow(df2),
    unique_participants_before = n_distinct(df2_raw$PATNO),
    unique_participants_after = n_distinct(df2$PATNO),
    unique_participant_visits_before = n_distinct(paste(df2_raw$PATNO, df2_raw$EVENT_ID)),
    unique_participant_visits_after = n_distinct(paste(df2$PATNO, df2$EVENT_ID))
)

write.csv(
    participant_visit_cleaning_summary,
    file.path(out_dir, "Participant_visit_cleaning_summary.csv"),
    row.names = FALSE
)

write.csv(
    df2,
    file.path(out_dir, "Participant_visit_level_dataset_used_for_table_and_attrition.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Baseline-only dataset
# ---------------------------------------------------------
bl <- df2 %>%
    filter(EVENT_ID == "BL") %>%
    distinct(PATNO, .keep_all = TRUE) %>%
    rename(
        age = !!age_col,
        SEX_raw = !!sex_col,
        BMI = !!bmi_col
    )

# ---------------------------------------------------------
# Recode categorical variables
# ---------------------------------------------------------
bl <- bl %>%
    mutate(
        SEX = case_when(
            as.character(SEX_raw) %in% c("0", "Female", "F", "female", "f") ~ "Female",
            as.character(SEX_raw) %in% c("1", "Male", "M", "male", "m") ~ "Male",
            TRUE ~ as.character(SEX_raw)
        ),
        SEX = factor(SEX, levels = c("Female", "Male"))
    )

if (!is.na(domside_col) && domside_col %in% names(bl)) {
    bl <- bl %>%
        rename(DOMSIDE_raw = !!domside_col) %>%
        mutate(
            DOMSIDE = case_when(
                as.character(DOMSIDE_raw) %in% c("1", "Left", "left", "L") ~ "Left",
                as.character(DOMSIDE_raw) %in% c("2", "Right", "right", "R") ~ "Right",
                as.character(DOMSIDE_raw) %in% c("3", "Symmetric", "symmetric", "Bilateral") ~ "Symmetric",
                TRUE ~ as.character(DOMSIDE_raw)
            ),
            DOMSIDE = factor(DOMSIDE, levels = c("Left", "Right", "Symmetric"))
        )
}

# ---------------------------------------------------------
# Create baseline NLR variable for Table 1
# ---------------------------------------------------------
if ("baseline_NLR" %in% names(bl)) {
    bl <- bl %>%
        mutate(NLR_table = baseline_NLR)
} else if ("NLR" %in% names(bl)) {
    bl <- bl %>%
        mutate(NLR_table = NLR)
} else if (all(c("Neutrophils", "Lymphocytes") %in% names(bl))) {
    bl <- bl %>%
        mutate(NLR_table = Neutrophils / Lymphocytes)
} else {
    warning("No NLR variable could be generated.")
    bl$NLR_table <- NA_real_
}

# ---------------------------------------------------------
# Force baseline LEDD to 0 if present
# ---------------------------------------------------------
if (!is.na(ledd_col) && ledd_col %in% names(bl)) {
    bl <- bl %>%
        rename(LEDD_raw = !!ledd_col) %>%
        mutate(LEDD_table = ifelse(GROUP == "PD", 0, NA_real_))
}

# ---------------------------------------------------------
# Rename optional clinical/imaging variables
# ---------------------------------------------------------
if (!is.na(duration_col) && duration_col %in% names(bl)) {
    bl <- bl %>% rename(duration_yrs = !!duration_col)
}

if (!is.na(updrs_on_col) && updrs_on_col %in% names(bl)) {
    bl <- bl %>% rename(UPDRSIII_ON = !!updrs_on_col)
}

if (!is.na(education_col) && education_col %in% names(bl)) {
    bl <- bl %>% rename(education_years = !!education_col)
}

if (!is.na(moca_col) && moca_col %in% names(bl)) {
    bl <- bl %>% rename(MoCA = !!moca_col)
}

if (!is.na(upsit_col) && upsit_col %in% names(bl)) {
    bl <- bl %>% rename(UPSIT = !!upsit_col)
}

if (!is.na(caudate_l_col) && caudate_l_col %in% names(bl)) {
    bl <- bl %>% rename(DAT_caudate_L = !!caudate_l_col)
}

if (!is.na(caudate_r_col) && caudate_r_col %in% names(bl)) {
    bl <- bl %>% rename(DAT_caudate_R = !!caudate_r_col)
}

if (!is.na(putamen_l_col) && putamen_l_col %in% names(bl)) {
    bl <- bl %>% rename(DAT_putamen_L = !!putamen_l_col)
}

if (!is.na(putamen_r_col) && putamen_r_col %in% names(bl)) {
    bl <- bl %>% rename(DAT_putamen_R = !!putamen_r_col)
}

# ---------------------------------------------------------
# Save baseline analysis dataset
# ---------------------------------------------------------
write.csv(
    bl,
    file.path(out_dir, "Baseline_dataset_for_Table1.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Table variables
# ---------------------------------------------------------
comparison_cont_vars <- c(
    "age",
    "BMI",
    "education_years",
    "MoCA",
    "UPSIT",
    "Neutrophils",
    "Lymphocytes",
    "Monocytes",
    "NLR_table"
)

comparison_cont_vars <- comparison_cont_vars[comparison_cont_vars %in% names(bl)]

comparison_cat_vars <- c("SEX")
comparison_cat_vars <- comparison_cat_vars[comparison_cat_vars %in% names(bl)]

pd_only_cont_vars <- c(
    "duration_yrs",
    "UPDRSIII_ON",
    "LEDD_table",
    "DAT_caudate_L",
    "DAT_caudate_R",
    "DAT_putamen_L",
    "DAT_putamen_R"
)

pd_only_cont_vars <- pd_only_cont_vars[pd_only_cont_vars %in% names(bl)]

pd_only_cat_vars <- c("DOMSIDE")
pd_only_cat_vars <- pd_only_cat_vars[pd_only_cat_vars %in% names(bl)]

labels <- list(
    age = "Age, years",
    BMI = "Body mass index, kg/m²",
    education_years = "Education, years",
    MoCA = "MoCA score",
    UPSIT = "UPSIT score",
    Neutrophils = "Absolute neutrophil count, ×10³/µL",
    Lymphocytes = "Absolute lymphocyte count, ×10³/µL",
    Monocytes = "Absolute monocyte count, ×10³/µL",
    NLR_table = "Neutrophil-to-lymphocyte ratio",
    SEX = "Male sex, n (%)",
    duration_yrs = "Disease duration, years",
    UPDRSIII_ON = "ON-medication UPDRS-III score",
    LEDD_table = "LEDD at baseline",
    DOMSIDE = "Dominant side of motor involvement",
    DAT_caudate_L = "Left caudate DAT binding",
    DAT_caudate_R = "Right caudate DAT binding",
    DAT_putamen_L = "Left putamen DAT binding",
    DAT_putamen_R = "Right putamen DAT binding"
)

# ---------------------------------------------------------
# Sample sizes
# ---------------------------------------------------------
n_hc <- bl %>% filter(GROUP == "HC") %>% nrow()
n_pd <- bl %>% filter(GROUP == "PD") %>% nrow()
n_overall <- nrow(bl)

# ---------------------------------------------------------
# Build comparison rows: continuous variables
# ---------------------------------------------------------
continuous_rows <- list()
welch_outputs <- list()

for (v in comparison_cont_vars) {
    
    test <- safe_welch(bl, v)
    wilcox_p <- safe_wilcox(bl, v)
    
    welch_outputs[[v]] <- test %>%
        mutate(
            label = labels[[v]],
            wilcoxon_sensitivity_p = wilcox_p,
            welch_p_formatted = fmt_p(p_value),
            wilcoxon_sensitivity_p_formatted = fmt_p(wilcox_p)
        )
    
    continuous_rows[[v]] <- tibble(
        Section = "Baseline demographic and peripheral inflammatory characteristics",
        Variable = labels[[v]],
        Overall = fmt_mean_sd(bl[[v]]),
        HC = fmt_mean_sd(bl %>% filter(GROUP == "HC") %>% pull(all_of(v))),
        PD = fmt_mean_sd(bl %>% filter(GROUP == "PD") %>% pull(all_of(v))),
        p = fmt_p(test$p_value),
        Test = "Welch t-test",
        Wilcoxon_sensitivity_p = fmt_p(wilcox_p)
    )
}

continuous_rows <- bind_rows(continuous_rows)
welch_outputs <- bind_rows(welch_outputs)

# ---------------------------------------------------------
# Build comparison rows: categorical variables
# Currently: sex summarized as male n (%)
# ---------------------------------------------------------
categorical_rows <- list()
chisq_outputs <- list()

if ("SEX" %in% comparison_cat_vars) {
    
    chi <- chisq_check(bl, "SEX", labels[["SEX"]])
    chisq_outputs[["SEX"]] <- chi
    
    denom_overall <- bl %>% filter(!is.na(SEX)) %>% nrow()
    denom_hc <- bl %>% filter(GROUP == "HC", !is.na(SEX)) %>% nrow()
    denom_pd <- bl %>% filter(GROUP == "PD", !is.na(SEX)) %>% nrow()
    
    n_male_overall <- bl %>% filter(SEX == "Male") %>% nrow()
    n_male_hc <- bl %>% filter(GROUP == "HC", SEX == "Male") %>% nrow()
    n_male_pd <- bl %>% filter(GROUP == "PD", SEX == "Male") %>% nrow()
    
    categorical_rows[["SEX"]] <- tibble(
        Section = "Baseline demographic and peripheral inflammatory characteristics",
        Variable = labels[["SEX"]],
        Overall = fmt_n_pct(n_male_overall, denom_overall),
        HC = fmt_n_pct(n_male_hc, denom_hc),
        PD = fmt_n_pct(n_male_pd, denom_pd),
        p = fmt_p(chi$selected_p),
        Test = chi$test_used,
        Wilcoxon_sensitivity_p = "—"
    )
}

categorical_rows <- bind_rows(categorical_rows)
chisq_outputs <- bind_rows(chisq_outputs)

# ---------------------------------------------------------
# Build PD-only rows: continuous clinical and imaging features
# ---------------------------------------------------------
pd_rows_cont <- list()

for (v in pd_only_cont_vars) {
    
    pd_rows_cont[[v]] <- tibble(
        Section = "PD-specific baseline clinical and imaging characteristics",
        Variable = labels[[v]],
        Overall = "—",
        HC = "—",
        PD = fmt_mean_sd(bl %>% filter(GROUP == "PD") %>% pull(all_of(v))),
        p = "—",
        Test = "PD only",
        Wilcoxon_sensitivity_p = "—"
    )
}

pd_rows_cont <- bind_rows(pd_rows_cont)

# ---------------------------------------------------------
# Build PD-only rows: dominant side
# ---------------------------------------------------------
pd_rows_cat <- list()

if ("DOMSIDE" %in% pd_only_cat_vars) {
    
    pd_dom <- bl %>%
        filter(GROUP == "PD", !is.na(DOMSIDE)) %>%
        count(DOMSIDE) %>%
        mutate(
            pct = 100 * n / sum(n),
            text = paste0(as.character(DOMSIDE), ": ", n, " (", sprintf("%.1f", pct), "%)")
        )
    
    pd_rows_cat[["DOMSIDE"]] <- tibble(
        Section = "PD-specific baseline clinical and imaging characteristics",
        Variable = labels[["DOMSIDE"]],
        Overall = "—",
        HC = "—",
        PD = paste(pd_dom$text, collapse = "; "),
        p = "—",
        Test = "PD only",
        Wilcoxon_sensitivity_p = "—"
    )
}

pd_rows_cat <- bind_rows(pd_rows_cat)

# ---------------------------------------------------------
# Final Table 1
# ---------------------------------------------------------
table1_df <- bind_rows(
    continuous_rows,
    categorical_rows,
    pd_rows_cont,
    pd_rows_cat
)

n_row <- tibble(
    Section = "Baseline demographic and peripheral inflammatory characteristics",
    Variable = "N",
    Overall = as.character(n_overall),
    HC = as.character(n_hc),
    PD = as.character(n_pd),
    p = "—",
    Test = "—",
    Wilcoxon_sensitivity_p = "—"
)

table1_df <- bind_rows(n_row, table1_df)

write.csv(
    table1_df,
    file.path(out_dir, "Table1_Baseline_Demographic_Clinical_Characteristics.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Formatted GT table
# ---------------------------------------------------------
gt_tbl <- table1_df %>%
    select(Section, Variable, Overall, HC, PD, p) %>%
    gt(groupname_col = "Section") %>%
    tab_header(
        title = md("**Table 1. Baseline demographic, clinical, imaging, and peripheral inflammatory characteristics**")
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
            "Continuous variables are shown as mean ± SD. Categorical variables are shown as n (%). Group comparisons were performed using Welch's t-tests for continuous variables and chi-square tests for categorical variables when expected-cell assumptions were satisfied; Fisher's exact test was used otherwise. PD-specific clinical and imaging variables are summarized for the PD group only."
        )
    )

gtsave(
    gt_tbl,
    file.path(out_dir, "Table1_Baseline_Demographic_Clinical_Characteristics.html")
)

gtsave(
    gt_tbl,
    file.path(out_dir, "Table1_Baseline_Demographic_Clinical_Characteristics.rtf")
)

# ---------------------------------------------------------
# Assumption checks for Welch t-tests
# ---------------------------------------------------------
norm_checks <- normality_checks(
    data = bl,
    vars = comparison_cont_vars,
    labels = labels,
    out_dir = out_dir
)

write.csv(
    norm_checks,
    file.path(out_dir, "Welch_ttest_normality_assumption_checks.csv"),
    row.names = FALSE
)

write.csv(
    welch_outputs,
    file.path(out_dir, "Welch_ttest_pvalues_with_Wilcoxon_sensitivity.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Assumption checks for chi-square tests
# ---------------------------------------------------------
write.csv(
    chisq_outputs,
    file.path(out_dir, "Chi_square_expected_cell_assumption_checks.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Follow-up / attrition counts
# ---------------------------------------------------------
df_follow <- df2 %>%
    filter(EVENT_ID %in% visit_order) %>%
    mutate(
        EVENT_ID = factor(EVENT_ID, levels = visit_order),
        GROUP = factor(GROUP, levels = c("HC", "PD"))
    )

# Create visit-wise NLR variable
if ("NLR" %in% names(df_follow)) {
    df_follow <- df_follow %>% mutate(NLR_visit = NLR)
} else if (all(c("Neutrophils", "Lymphocytes") %in% names(df_follow))) {
    df_follow <- df_follow %>% mutate(NLR_visit = Neutrophils / Lymphocytes)
} else {
    df_follow <- df_follow %>% mutate(NLR_visit = NA_real_)
}

# Create UPDRS variable in long data
if (!is.na(updrs_on_col) && updrs_on_col %in% names(df_follow)) {
    df_follow <- df_follow %>%
        mutate(UPDRSIII_ON_visit = .data[[updrs_on_col]])
} else {
    df_follow <- df_follow %>%
        mutate(UPDRSIII_ON_visit = NA_real_)
}

# Create DAT availability indicator
dat_cols_available <- c(
    caudate_l_col,
    caudate_r_col,
    putamen_l_col,
    putamen_r_col
)

dat_cols_available <- dat_cols_available[!is.na(dat_cols_available)]
dat_cols_available <- dat_cols_available[dat_cols_available %in% names(df_follow)]

if (length(dat_cols_available) > 0) {
    df_follow <- df_follow %>%
        mutate(
            DAT_any_available = ifelse(
                rowSums(!is.na(across(all_of(dat_cols_available)))) > 0,
                TRUE,
                FALSE
            )
        )
} else {
    df_follow <- df_follow %>%
        mutate(DAT_any_available = FALSE)
}

visit_counts <- df_follow %>%
    group_by(EVENT_ID, visit_label, year_c, GROUP) %>%
    summarise(
        n_subjects_any_record = n_distinct(PATNO),
        n_subjects_with_NLR = n_distinct(PATNO[!is.na(NLR_visit)]),
        n_subjects_with_ON_UPDRSIII = n_distinct(PATNO[!is.na(UPDRSIII_ON_visit)]),
        n_subjects_with_any_DAT = n_distinct(PATNO[DAT_any_available]),
        .groups = "drop"
    ) %>%
    arrange(GROUP, year_c)

write.csv(
    visit_counts,
    file.path(out_dir, "Followup_attrition_counts_long.csv"),
    row.names = FALSE
)

visit_counts_wide <- visit_counts %>%
    select(EVENT_ID, visit_label, year_c, GROUP, n_subjects_any_record) %>%
    pivot_wider(
        names_from = GROUP,
        values_from = n_subjects_any_record
    ) %>%
    arrange(year_c)

write.csv(
    visit_counts_wide,
    file.path(out_dir, "Followup_attrition_counts_wide_any_record.csv"),
    row.names = FALSE
)

visit_counts_wide_all_availability <- visit_counts %>%
    select(
        EVENT_ID,
        visit_label,
        year_c,
        GROUP,
        n_subjects_any_record,
        n_subjects_with_NLR,
        n_subjects_with_ON_UPDRSIII,
        n_subjects_with_any_DAT
    ) %>%
    pivot_wider(
        names_from = GROUP,
        values_from = c(
            n_subjects_any_record,
            n_subjects_with_NLR,
            n_subjects_with_ON_UPDRSIII,
            n_subjects_with_any_DAT
        )
    ) %>%
    arrange(year_c)

write.csv(
    visit_counts_wide_all_availability,
    file.path(out_dir, "Followup_attrition_counts_wide_all_availability.csv"),
    row.names = FALSE
)

# ---------------------------------------------------------
# Attrition flow chart
# ---------------------------------------------------------
flow_df <- visit_counts %>%
    select(EVENT_ID, visit_label, year_c, GROUP, n_subjects_any_record) %>%
    mutate(
        x = as.numeric(factor(EVENT_ID, levels = visit_order)),
        y = ifelse(GROUP == "HC", 2, 1),
        box_label = paste0(
            visit_label,
            "\n",
            GROUP,
            " n = ",
            n_subjects_any_record
        )
    )

arrow_df <- flow_df %>%
    arrange(GROUP, x) %>%
    group_by(GROUP) %>%
    mutate(
        xend = lead(x),
        yend = lead(y)
    ) %>%
    filter(!is.na(xend)) %>%
    ungroup()

p_flow <- ggplot() +
    geom_segment(
        data = arrow_df,
        aes(x = x + 0.38, xend = xend - 0.38, y = y, yend = yend),
        arrow = arrow(length = grid::unit(0.15, "inches")),
        linewidth = 0.5
    ) +
    geom_label(
        data = flow_df,
        aes(x = x, y = y, label = box_label),
        size = 3.2,
        label.size = 0.35,
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
        subtitle = "Number of participants with any available record at each visit",
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
    filename = file.path(out_dir, "Followup_attrition_flowchart.png"),
    plot = p_flow,
    width = 13,
    height = 4.8,
    dpi = 600
)

ggsave(
    filename = file.path(out_dir, "Followup_attrition_flowchart.pdf"),
    plot = p_flow,
    width = 13,
    height = 4.8
)

# ---------------------------------------------------------
# Attrition line plot
# ---------------------------------------------------------
p_attrition <- ggplot(
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
    filename = file.path(out_dir, "Followup_attrition_lineplot.png"),
    plot = p_attrition,
    width = 8.5,
    height = 5.5,
    dpi = 600
)

ggsave(
    filename = file.path(out_dir, "Followup_attrition_lineplot.pdf"),
    plot = p_attrition,
    width = 8.5,
    height = 5.5
)

# ---------------------------------------------------------
# Sample size justification text
# ---------------------------------------------------------
sample_text <- c(
    "Sample size justification",
    "=========================",
    "",
    paste0(
        "Baseline analyses included all eligible PPMI participants with diagnostic classification as PD or HC and an available baseline visit record. This yielded ",
        n_pd,
        " participants with PD and ",
        n_hc,
        " healthy controls."
    ),
    "",
    "No formal a priori power calculation was performed because this was a secondary analysis of an existing longitudinal observational cohort. The analytic sample size was therefore determined by the number of eligible participants with available baseline blood count data, clinical assessments, imaging measures, and covariates. To support interpretation, effect estimates are reported with 95% confidence intervals and precise p values, and missingness and follow-up availability are summarized at both the observation and participant levels."
)

writeLines(
    sample_text,
    file.path(out_dir, "Sample_size_justification_text.txt")
)

# ---------------------------------------------------------
# Manuscript-ready baseline summary starter
# ---------------------------------------------------------
summary_text <- c(
    "Baseline demographic and clinical characteristics",
    "===============================================",
    "",
    paste0(
        "Baseline analyses included ",
        n_pd,
        " participants with PD and ",
        n_hc,
        " healthy controls."
    ),
    "",
    "Baseline demographic, clinical, imaging, and peripheral inflammatory characteristics are summarized in Table 1. Continuous variables were compared using Welch's t-tests after assessment of distributional assumptions using Q-Q plots and Shapiro-Wilk tests. Categorical variables were compared using chi-square tests when expected-cell assumptions were satisfied; otherwise, Fisher's exact test was used. Follow-up availability and attrition are summarized in the accompanying attrition table and flow chart."
)

writeLines(
    summary_text,
    file.path(out_dir, "Manuscript_baseline_summary_starter.txt")
)

# ---------------------------------------------------------
# Excel workbook with key outputs
# ---------------------------------------------------------
wb <- openxlsx::createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
    sheet_name <- substr(sheet_name, 1, 31)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = data)
}

add_sheet(wb, "Table1", table1_df)
add_sheet(wb, "Baseline_dataset", bl)
add_sheet(wb, "Detected_columns", detected_columns)
add_sheet(wb, "Duplicate_summary", duplicate_audit_summary)
add_sheet(wb, "Duplicate_by_visit", duplicate_summary_by_group_visit)
add_sheet(wb, "Duplicate_discrepancy", duplicate_discrepancy_by_variable)
add_sheet(wb, "Cleaning_summary", participant_visit_cleaning_summary)
add_sheet(wb, "Followup_long", visit_counts)
add_sheet(wb, "Followup_wide", visit_counts_wide_all_availability)
add_sheet(wb, "Welch_Wilcoxon", welch_outputs)
add_sheet(wb, "Normality_checks", norm_checks)
add_sheet(wb, "Chi_square_checks", chisq_outputs)

openxlsx::saveWorkbook(
    wb,
    file = file.path(out_dir, "Table1_and_attrition_outputs.xlsx"),
    overwrite = TRUE
)

# ---------------------------------------------------------
# Console output
# ---------------------------------------------------------
cat("\nBaseline sample sizes:\n")
cat("HC =", n_hc, "\n")
cat("PD =", n_pd, "\n")
cat("Overall =", n_overall, "\n")

cat("\nPATNO-EVENT_ID duplicate audit:\n")
print(duplicate_audit_summary)

cat("\nParticipant-visit cleaning summary:\n")
print(participant_visit_cleaning_summary)

cat("\nComparison continuous variables included:\n")
print(comparison_cont_vars)

cat("\nPD-only continuous variables included:\n")
print(pd_only_cont_vars)

cat("\nVisit counts by group:\n")
print(visit_counts_wide)

cat("\nAll outputs saved to:\n", out_dir, "\n")