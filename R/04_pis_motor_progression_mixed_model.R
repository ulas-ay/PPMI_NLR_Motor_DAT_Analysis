# =========================================================
# Baseline PIS and Motor Progression in Parkinson's Disease
#
# This script:
#   1. Restricts the dataset to Parkinson's disease participants
#   2. Reconstructs baseline PCA-derived Peripheral Immune Score (PIS)
#      using neutrophil, lymphocyte, and monocyte counts
#   3. Carries baseline PIS forward across all longitudinal visits
#   4. Tests whether baseline PIS predicts longitudinal motor severity
#      using a linear mixed-effects model
#   5. Exports model summaries, simple slopes, predicted trajectories,
#      tertile-based summaries, and figures
#
# Required input:
#   data/PPMI_analysis_dataset.xlsx
#
# Main model:
#   updrs3_score_on ~ PIS_BL * year_c + age + sex + bmi +
#                     duration_yrs + LEDD + DOMSIDE +
#                     (1 | SITE) + (1 | PATNO)
#
# Outputs:
#   outputs/pis_motor_progression/PD_PIS_motor_progression_dataset.csv
#   outputs/pis_motor_progression/PD_PIS_motor_progression_descriptives.csv
#   outputs/pis_motor_progression/LMM_fixed_effects.csv
#   outputs/pis_motor_progression/LMM_TypeIII_ANOVA.csv
#   outputs/pis_motor_progression/LMM_random_effects_variance.csv
#   outputs/pis_motor_progression/PIS_effect_on_UPDRS_by_timepoint.csv
#   outputs/pis_motor_progression/Predicted_UPDRS_trajectories_by_baseline_PIS_continuous.csv
#   outputs/pis_motor_progression/PD_PIS_motor_progression_dataset_with_tertiles.csv
#   outputs/pis_motor_progression/Baseline_PIS_tertile_means.csv
#   outputs/pis_motor_progression/Predicted_UPDRS_trajectories_by_PIS_tertile.csv
#   outputs/pis_motor_progression/Figure1_PISxTime_UPDRS_continuous.png
#   outputs/pis_motor_progression/Figure1_PISxTime_UPDRS_continuous.pdf
#   outputs/pis_motor_progression/Figure2_PIS_tertiles_UPDRS.png
#   outputs/pis_motor_progression/Figure2_PIS_tertiles_UPDRS.pdf
#   outputs/pis_motor_progression/Figure3_raw_UPDRS_by_PIS_tertile.png
#   outputs/pis_motor_progression/LMM_summary.txt
# =========================================================


# ---------------------------------------------------------
# Load required packages
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
out_dir   <- file.path("outputs", "pis_motor_progression")

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
# Helper functions
# ---------------------------------------------------------

safe_ratio <- function(num, den) {
  out <- num / den
  out[is.infinite(out)] <- NA_real_
  out
}


find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}


standardize_ci_names <- function(x) {
  if (!"lower.CL" %in% names(x) && "asymp.LCL" %in% names(x)) {
    x <- dplyr::rename(x, lower.CL = asymp.LCL)
  }
  
  if (!"upper.CL" %in% names(x) && "asymp.UCL" %in% names(x)) {
    x <- dplyr::rename(x, upper.CL = asymp.UCL)
  }
  
  if (!"lower.CL" %in% names(x) && "lowerCL" %in% names(x)) {
    x <- dplyr::rename(x, lower.CL = lowerCL)
  }
  
  if (!"upper.CL" %in% names(x) && "upperCL" %in% names(x)) {
    x <- dplyr::rename(x, upper.CL = upperCL)
  }
  
  if (!all(c("lower.CL", "upper.CL") %in% names(x))) {
    print(names(x))
    stop("Confidence interval columns were not found in the emmeans output.")
  }
  
  x
}


# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------

df <- readxl::read_excel(data_path)

cat("\nData successfully loaded.\n")
cat("Number of rows:", nrow(df), "\n")
cat("Number of columns:", ncol(df), "\n\n")


# ---------------------------------------------------------
# Define visit-year mapping
# ---------------------------------------------------------

event_map <- tibble::tribble(
  ~EVENT_ID, ~year,
  "BL",   0,
  "V04",  1,
  "V06",  2,
  "V08",  3,
  "V10",  4,
  "V12",  5,
  "V13",  6,
  "V14",  7
)

time_points <- event_map$year


# ---------------------------------------------------------
# Detect key columns robustly
# ---------------------------------------------------------

age_col <- find_col(
  df,
  c("age")
)

sex_col <- find_col(
  df,
  c("SEX")
)

bmi_col <- find_col(
  df,
  c("BMI")
)

site_col <- find_col(
  df,
  c("SITE")
)

duration_col <- find_col(
  df,
  c("duration_yrs")
)

updrs_col <- find_col(
  df,
  c("updrs3_score_on")
)

ledd_col <- find_col(
  df,
  c("LEDD")
)

domside_col <- find_col(
  df,
  c("DOMSIDE")
)

cat("Detected columns:\n")
cat("age_col      =", age_col, "\n")
cat("sex_col      =", sex_col, "\n")
cat("bmi_col      =", bmi_col, "\n")
cat("site_col     =", site_col, "\n")
cat("duration_col =", duration_col, "\n")
cat("updrs_col    =", updrs_col, "\n")
cat("ledd_col     =", ledd_col, "\n")
cat("domside_col  =", domside_col, "\n\n")

detected_cols <- c(
  age_col,
  sex_col,
  bmi_col,
  site_col,
  duration_col,
  updrs_col,
  ledd_col,
  domside_col
)

if (any(is.na(detected_cols))) {
  stop(
    paste0(
      "Could not detect one or more required columns. ",
      "Please check names(df) and update the candidate names in the script."
    )
  )
}


# ---------------------------------------------------------
# Check required columns
# ---------------------------------------------------------

immune_vars <- c(
  "Neutrophils",
  "Lymphocytes",
  "Monocytes"
)

required_vars <- c(
  "PATNO",
  "EVENT_ID",
  "PRIMDIAG",
  immune_vars,
  detected_cols
)

missing_vars <- setdiff(required_vars, names(df))

if (length(missing_vars) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the dataset:\n",
      paste(missing_vars, collapse = ", "),
      "\n\nPlease check the column names in the input dataset."
    )
  )
}


# ---------------------------------------------------------
# Prepare PD-only longitudinal dataset
#
# PRIMDIAG:
#   1 = Parkinson's disease
# ---------------------------------------------------------

df2 <- df %>%
  dplyr::mutate(EVENT_ID = as.character(EVENT_ID)) %>%
  dplyr::inner_join(event_map, by = "EVENT_ID") %>%
  dplyr::filter(PRIMDIAG == 1) %>%
  dplyr::mutate(GROUP = "PD")

cat("Number of PD rows after visit restriction:", nrow(df2), "\n")

cat("\nRemaining EVENT_ID distribution:\n")
print(table(df2$EVENT_ID, useNA = "ifany"))


# ---------------------------------------------------------
# Compute baseline PCA-derived PIS
#
# PIS is derived from baseline neutrophil, lymphocyte, and
# monocyte counts. The PC1 direction is aligned so that higher
# PIS corresponds to higher NLR-related immune activation.
# ---------------------------------------------------------

bl <- df2 %>%
  dplyr::filter(EVENT_ID == "BL") %>%
  dplyr::mutate(
    NLR = safe_ratio(Neutrophils, Lymphocytes),
    MLR = safe_ratio(Monocytes, Lymphocytes)
  ) %>%
  dplyr::select(
    PATNO,
    EVENT_ID,
    year,
    dplyr::all_of(immune_vars),
    NLR,
    MLR
  )

pca_input_bl <- bl %>%
  dplyr::select(PATNO, dplyr::all_of(immune_vars)) %>%
  tidyr::drop_na()

cat("\nNumber of PD participants with complete baseline immune data:", nrow(pca_input_bl), "\n")

if (nrow(pca_input_bl) < 5) {
  stop(
    paste0(
      "Insufficient complete baseline observations for PCA. ",
      "At least 5 complete rows are recommended."
    )
  )
}

pca_fit <- stats::prcomp(
  pca_input_bl[, immune_vars],
  scale. = TRUE
)

bl_scores_raw <- tibble::tibble(
  PATNO = pca_input_bl$PATNO,
  PIS_raw = as.numeric(pca_fit$x[, 1])
)

bl_scores <- bl %>%
  dplyr::left_join(bl_scores_raw, by = "PATNO")


# ---------------------------------------------------------
# Orient PC1 to correlate positively with NLR
# ---------------------------------------------------------

tmp_cor_nlr <- suppressWarnings(
  stats::cor(
    bl_scores$PIS_raw,
    bl_scores$NLR,
    use = "pairwise.complete.obs",
    method = "spearman"
  )
)

if (!is.na(tmp_cor_nlr) && tmp_cor_nlr < 0) {
  rotation_pc1 <- -pca_fit$rotation[, 1]
} else {
  rotation_pc1 <- pca_fit$rotation[, 1]
}

cat("\nPC1 orientation check: Spearman correlation between PIS_raw and NLR =\n")
print(tmp_cor_nlr)


# ---------------------------------------------------------
# Final baseline PIS
# ---------------------------------------------------------

scaled_bl <- scale(
  bl_scores[, immune_vars],
  center = pca_fit$center,
  scale = pca_fit$scale
)

bl_pis <- tibble::tibble(
  PATNO = bl_scores$PATNO,
  PIS_BL = as.numeric(scaled_bl %*% rotation_pc1)
) %>%
  dplyr::distinct(PATNO, .keep_all = TRUE)


# ---------------------------------------------------------
# Build PD longitudinal dataset
#
# Baseline PIS is carried forward to all visits.
# Baseline LEDD is set to 0 to represent the drug-naive
# baseline state when applicable.
# ---------------------------------------------------------

dat <- df2 %>%
  dplyr::left_join(bl_pis, by = "PATNO") %>%
  dplyr::select(
    PATNO,
    EVENT_ID,
    year,
    dplyr::all_of(site_col),
    dplyr::all_of(age_col),
    dplyr::all_of(sex_col),
    dplyr::all_of(bmi_col),
    dplyr::all_of(duration_col),
    dplyr::all_of(ledd_col),
    dplyr::all_of(domside_col),
    dplyr::all_of(updrs_col),
    PIS_BL
  ) %>%
  dplyr::rename(
    SITE = dplyr::all_of(site_col),
    age = dplyr::all_of(age_col),
    sex = dplyr::all_of(sex_col),
    bmi = dplyr::all_of(bmi_col),
    duration_yrs = dplyr::all_of(duration_col),
    LEDD_raw = dplyr::all_of(ledd_col),
    DOMSIDE = dplyr::all_of(domside_col),
    updrs3_score_on = dplyr::all_of(updrs_col)
  ) %>%
  dplyr::mutate(
    SITE = as.factor(SITE),
    sex = as.factor(sex),
    DOMSIDE = as.factor(DOMSIDE),
    year_c = year,
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    duration_yrs = as.numeric(duration_yrs),
    LEDD_raw = as.numeric(LEDD_raw),
    updrs3_score_on = as.numeric(updrs3_score_on),
    LEDD = ifelse(year_c == 0, 0, LEDD_raw)
  ) %>%
  tidyr::drop_na(
    PIS_BL,
    age,
    sex,
    bmi,
    duration_yrs,
    LEDD,
    DOMSIDE,
    updrs3_score_on,
    SITE
  ) %>%
  droplevels()

if (dplyr::n_distinct(dat$PATNO) < 3) {
  stop("Too few unique PD participants are available for mixed-effects modeling.")
}


# ---------------------------------------------------------
# Save analysis dataset
# ---------------------------------------------------------

write.csv(
  dat,
  file.path(out_dir, "PD_PIS_motor_progression_dataset.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Descriptive summaries
# ---------------------------------------------------------

desc <- dat %>%
  dplyr::group_by(EVENT_ID, year) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    mean_UPDRS = mean(updrs3_score_on, na.rm = TRUE),
    sd_UPDRS = stats::sd(updrs3_score_on, na.rm = TRUE),
    mean_PIS_BL = mean(PIS_BL, na.rm = TRUE),
    sd_PIS_BL = stats::sd(PIS_BL, na.rm = TRUE),
    mean_LEDD = mean(LEDD, na.rm = TRUE),
    sd_LEDD = stats::sd(LEDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year)

write.csv(
  desc,
  file.path(out_dir, "PD_PIS_motor_progression_descriptives.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Linear mixed-effects model
#
# Sum-to-zero contrasts are used for Type III tests.
# ---------------------------------------------------------

old_contrasts <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))

lmm <- lmerTest::lmer(
  updrs3_score_on ~ PIS_BL * year_c +
    age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
    (1 | SITE) + (1 | PATNO),
  data = dat,
  REML = FALSE,
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 100000)
  )
)

model_summary <- summary(lmm)
anova_type3 <- car::Anova(lmm, type = 3)


# ---------------------------------------------------------
# Export model outputs
# ---------------------------------------------------------

fixef_tab <- as.data.frame(coef(summary(lmm)))
fixef_tab$term <- rownames(fixef_tab)

fixef_tab <- fixef_tab %>%
  dplyr::relocate(term)

write.csv(
  fixef_tab,
  file.path(out_dir, "LMM_fixed_effects.csv"),
  row.names = FALSE
)

anova_tab <- as.data.frame(anova_type3)
anova_tab$term <- rownames(anova_tab)

anova_tab <- anova_tab %>%
  dplyr::relocate(term)

write.csv(
  anova_tab,
  file.path(out_dir, "LMM_TypeIII_ANOVA.csv"),
  row.names = FALSE
)

randef_var <- as.data.frame(lme4::VarCorr(lmm))

write.csv(
  randef_var,
  file.path(out_dir, "LMM_random_effects_variance.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Simple slopes of baseline PIS effect at each time point
# ---------------------------------------------------------

slopes_by_time <- emmeans::emtrends(
  lmm,
  ~ year_c,
  var = "PIS_BL",
  at = list(year_c = time_points)
)

slopes_by_time_df <- as.data.frame(
  summary(slopes_by_time, infer = c(TRUE, TRUE))
)

slopes_by_time_df <- standardize_ci_names(slopes_by_time_df)

write.csv(
  slopes_by_time_df,
  file.path(out_dir, "PIS_effect_on_UPDRS_by_timepoint.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Predicted trajectories for low, mean, and high baseline PIS
#
# Low and high PIS are defined as mean ± 1.5 SD.
# ---------------------------------------------------------

pis_mean <- mean(dat$PIS_BL, na.rm = TRUE)
pis_sd <- stats::sd(dat$PIS_BL, na.rm = TRUE)

pis_levels <- c(
  pis_mean - 1.5 * pis_sd,
  pis_mean,
  pis_mean + 1.5 * pis_sd
)

pis_labels <- c(
  "Low baseline PIS (-1.5 SD)",
  "Mean baseline PIS",
  "High baseline PIS (+1.5 SD)"
)

emm_cont <- emmeans::emmeans(
  lmm,
  ~ year_c | PIS_BL,
  at = list(
    year_c = time_points,
    PIS_BL = pis_levels,
    age = mean(dat$age, na.rm = TRUE),
    bmi = mean(dat$bmi, na.rm = TRUE),
    duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
    LEDD = mean(dat$LEDD, na.rm = TRUE)
  )
)

emm_cont_df <- as.data.frame(
  summary(emm_cont, infer = c(TRUE, TRUE))
)

emm_cont_df <- standardize_ci_names(emm_cont_df)

emm_cont_df <- emm_cont_df %>%
  dplyr::mutate(
    PIS_group = factor(
      PIS_BL,
      levels = pis_levels,
      labels = pis_labels
    )
  )

write.csv(
  emm_cont_df,
  file.path(out_dir, "Predicted_UPDRS_trajectories_by_baseline_PIS_continuous.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Tertile-based model-estimated trajectories
#
# Tertiles are used for visualization only. The primary model
# treats PIS_BL as a continuous predictor.
# ---------------------------------------------------------

bl_tertiles <- bl_pis %>%
  dplyr::mutate(
    PIS_tertile = dplyr::ntile(PIS_BL, 3),
    PIS_tertile = factor(
      PIS_tertile,
      levels = c(1, 2, 3),
      labels = c(
        "Low PIS tertile",
        "Middle PIS tertile",
        "High PIS tertile"
      )
    )
  ) %>%
  dplyr::filter(!is.na(PIS_tertile))

dat_tert <- dat %>%
  dplyr::left_join(bl_tertiles, by = "PATNO") %>%
  dplyr::filter(!is.na(PIS_tertile)) %>%
  dplyr::mutate(
    PIS_tertile = factor(
      PIS_tertile,
      levels = c(
        "Low PIS tertile",
        "Middle PIS tertile",
        "High PIS tertile"
      )
    )
  )

write.csv(
  dat_tert,
  file.path(out_dir, "PD_PIS_motor_progression_dataset_with_tertiles.csv"),
  row.names = FALSE
)

tert_summary <- bl_tertiles %>%
  dplyr::group_by(PIS_tertile) %>%
  dplyr::summarise(
    PIS_BL_mean = mean(PIS_BL, na.rm = TRUE),
    n_subjects = dplyr::n_distinct(PATNO),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    match(
      PIS_tertile,
      c(
        "Low PIS tertile",
        "Middle PIS tertile",
        "High PIS tertile"
      )
    )
  )

write.csv(
  tert_summary,
  file.path(out_dir, "Baseline_PIS_tertile_means.csv"),
  row.names = FALSE
)

tert_means <- tert_summary$PIS_BL_mean
tert_labels <- as.character(tert_summary$PIS_tertile)

emm_tert <- emmeans::emmeans(
  lmm,
  ~ year_c | PIS_BL,
  at = list(
    year_c = time_points,
    PIS_BL = tert_means,
    age = mean(dat$age, na.rm = TRUE),
    bmi = mean(dat$bmi, na.rm = TRUE),
    duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
    LEDD = mean(dat$LEDD, na.rm = TRUE)
  )
)

emm_tert_df <- as.data.frame(
  summary(emm_tert, infer = c(TRUE, TRUE))
)

emm_tert_df <- standardize_ci_names(emm_tert_df)

emm_tert_df <- emm_tert_df %>%
  dplyr::mutate(
    PIS_tertile = factor(
      PIS_BL,
      levels = tert_means,
      labels = tert_labels
    )
  ) %>%
  dplyr::filter(!is.na(PIS_tertile)) %>%
  dplyr::mutate(
    PIS_tertile = droplevels(PIS_tertile)
  )

write.csv(
  emm_tert_df,
  file.path(out_dir, "Predicted_UPDRS_trajectories_by_PIS_tertile.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Figure 1: continuous baseline PIS interaction plot
# ---------------------------------------------------------

p1 <- ggplot2::ggplot(
  emm_cont_df,
  ggplot2::aes(
    x = year_c,
    y = emmean,
    color = PIS_group,
    fill = PIS_group,
    group = PIS_group
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
    alpha = 0.18,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::scale_x_continuous(
    breaks = event_map$year,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  ggplot2::labs(
    title = "Baseline PIS predicts longitudinal motor severity in Parkinson's disease",
    subtitle = "Model-estimated UPDRS-III trajectories at low, mean, and high baseline PIS",
    x = "Years",
    y = "Predicted UPDRS-III score",
    color = "Baseline PIS",
    fill = "Baseline PIS"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "Figure1_PISxTime_UPDRS_continuous.png"),
  plot = p1,
  width = 9.5,
  height = 6.8,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Figure1_PISxTime_UPDRS_continuous.pdf"),
  plot = p1,
  width = 9.5,
  height = 6.8
)


# ---------------------------------------------------------
# Figure 2: model-estimated trajectories by PIS tertile
# ---------------------------------------------------------

p2 <- ggplot2::ggplot(
  emm_tert_df,
  ggplot2::aes(
    x = year_c,
    y = emmean,
    color = PIS_tertile,
    fill = PIS_tertile,
    group = PIS_tertile
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
    alpha = 0.18,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::scale_x_continuous(
    breaks = event_map$year,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  ggplot2::scale_color_discrete(drop = TRUE, na.translate = FALSE) +
  ggplot2::scale_fill_discrete(drop = TRUE, na.translate = FALSE) +
  ggplot2::labs(
    title = "Longitudinal motor trajectories stratified by baseline PIS tertiles",
    subtitle = "Model-estimated UPDRS-III trajectories in Parkinson's disease",
    x = "Years",
    y = "Predicted UPDRS-III score",
    color = "Baseline PIS tertile",
    fill = "Baseline PIS tertile"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "Figure2_PIS_tertiles_UPDRS.png"),
  plot = p2,
  width = 9.5,
  height = 6.8,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Figure2_PIS_tertiles_UPDRS.pdf"),
  plot = p2,
  width = 9.5,
  height = 6.8
)


# ---------------------------------------------------------
# Figure 3: raw descriptive trajectories by tertile
# ---------------------------------------------------------

raw_tert_df <- dat_tert %>%
  dplyr::group_by(PIS_tertile, year_c, EVENT_ID) %>%
  dplyr::summarise(
    mean_UPDRS = mean(updrs3_score_on, na.rm = TRUE),
    se_UPDRS = stats::sd(updrs3_score_on, na.rm = TRUE) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(PIS_tertile))

p3 <- ggplot2::ggplot(
  raw_tert_df,
  ggplot2::aes(
    x = year_c,
    y = mean_UPDRS,
    color = PIS_tertile,
    group = PIS_tertile
  )
) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_UPDRS - se_UPDRS,
      ymax = mean_UPDRS + se_UPDRS
    ),
    width = 0.12
  ) +
  ggplot2::scale_x_continuous(
    breaks = event_map$year,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  ggplot2::scale_color_discrete(drop = TRUE, na.translate = FALSE) +
  ggplot2::labs(
    title = "Observed UPDRS-III means by baseline PIS tertile",
    subtitle = "Raw descriptive trajectories",
    x = "Years",
    y = "Observed mean UPDRS-III",
    color = "Baseline PIS tertile"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "Figure3_raw_UPDRS_by_PIS_tertile.png"),
  plot = p3,
  width = 9.5,
  height = 6.8,
  dpi = 600
)


# ---------------------------------------------------------
# Save full text summary
# ---------------------------------------------------------

sink(file.path(out_dir, "LMM_summary.txt"))

cat("============================================================\n")
cat("Baseline PIS and motor progression in Parkinson's disease\n")
cat("============================================================\n\n")

cat("Model formula:\n")
cat(
  "updrs3_score_on ~ PIS_BL * year_c + age + sex + bmi + ",
  "duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n",
  sep = ""
)

cat("Detected covariate columns:\n")
cat("age       =", age_col, "\n")
cat("sex       =", sex_col, "\n")
cat("bmi       =", bmi_col, "\n")
cat("site      =", site_col, "\n")
cat("duration  =", duration_col, "\n")
cat("LEDD      =", ledd_col, "\n")
cat("DOMSIDE   =", domside_col, "\n")
cat("UPDRS-III =", updrs_col, "\n\n")

cat("Sample size:\n")
cat("Rows =", nrow(dat), "\n")
cat("Unique PD subjects =", dplyr::n_distinct(dat$PATNO), "\n")
cat("Unique sites =", dplyr::n_distinct(dat$SITE), "\n\n")

cat("Visit counts:\n")
print(table(dat$EVENT_ID))
cat("\n")

cat("Random-effects variance:\n")
print(lme4::VarCorr(lmm), comp = c("Variance", "Std.Dev."))
cat("\n\n")

cat("Fixed effects summary:\n")
print(model_summary)
cat("\n\n")

cat("Type III ANOVA:\n")
print(anova_type3)
cat("\n\n")

cat("Simple slopes of baseline PIS effect at each time point:\n")
print(slopes_by_time_df, row.names = FALSE)
cat("\n\n")

cat("PIS levels used in Figure 1, mean ± 1.5 SD:\n")
print(
  data.frame(
    label = pis_labels,
    PIS_value = pis_levels
  ),
  row.names = FALSE
)
cat("\n\n")

cat("Baseline PIS tertile means used in Figure 2:\n")
print(tert_summary, row.names = FALSE)
cat("\n\n")

sink()


# ---------------------------------------------------------
# Restore contrast options
# ---------------------------------------------------------

options(old_contrasts)


# ---------------------------------------------------------
# Console output
# ---------------------------------------------------------

cat("\n============================================================\n")
cat("PIS-motor progression mixed-effects analysis completed.\n")
cat("Outputs saved to:\n")
cat(out_dir, "\n")
cat("============================================================\n")