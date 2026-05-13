# =========================================================
# Longitudinal Analysis of Peripheral Immune Score (PIS)
# HC vs PD
#
# This script:
#   1. Restricts the dataset to Parkinson's disease participants
#      and healthy controls based on PRIMDIAG coding
#   2. Reconstructs the baseline PCA-derived Peripheral Immune Score
#      using neutrophil, lymphocyte, and monocyte counts
#   3. Projects all available visits onto the baseline PCA axis
#   4. Tests longitudinal group differences in PIS using a linear
#      mixed-effects model
#   5. Estimates group-specific trajectories, group contrasts at each
#      time point, and group-specific slopes
#   6. Exports model tables, descriptive summaries, and figures
#
# Required input:
#   data/PPMI_analysis_dataset.xlsx
#
# Main model:
#   PIS ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO)
#
# Outputs:
#   outputs/pis_longitudinal/PIS_longitudinal_dataset.csv
#   outputs/pis_longitudinal/PIS_longitudinal_descriptives.csv
#   outputs/pis_longitudinal/SITE_distribution.csv
#   outputs/pis_longitudinal/LMM_fixed_effects.csv
#   outputs/pis_longitudinal/LMM_TypeIII_ANOVA.csv
#   outputs/pis_longitudinal/LMM_random_effects_variance.csv
#   outputs/pis_longitudinal/EMMeans_PIS_by_group_and_time.csv
#   outputs/pis_longitudinal/Group_comparisons_at_each_timepoint.csv
#   outputs/pis_longitudinal/Group_specific_slopes.csv
#   outputs/pis_longitudinal/Slope_difference_between_groups.csv
#   outputs/pis_longitudinal/Baseline_EMMeans.csv
#   outputs/pis_longitudinal/Baseline_group_difference.csv
#   outputs/pis_longitudinal/PIS_longitudinal_plot.png
#   outputs/pis_longitudinal/PIS_longitudinal_plot.pdf
#   outputs/pis_longitudinal/PIS_spaghetti_plot.png
#   outputs/pis_longitudinal/LMM_summary.txt
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
out_dir   <- file.path("outputs", "pis_longitudinal")

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


standardize_emmeans_ci <- function(x) {
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
# Detect covariate columns robustly
#
# This allows the script to work if the same variables have
# slightly different names across curated versions of the dataset.
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

cat("Detected covariate columns:\n")
cat("age_col  =", age_col, "\n")
cat("sex_col  =", sex_col, "\n")
cat("bmi_col  =", bmi_col, "\n")
cat("site_col =", site_col, "\n\n")

if (is.na(age_col) || is.na(sex_col) || is.na(bmi_col) || is.na(site_col)) {
  stop(
    paste0(
      "Could not detect one or more required covariate columns: ",
      "age, sex, BMI, or SITE. Please check names(df) and update ",
      "the candidate names in the script."
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
  age_col,
  sex_col,
  bmi_col,
  site_col
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
# Define diagnostic groups and retain selected visits
#
# PRIMDIAG:
#   1  = Parkinson's disease
#   17 = Healthy control
# ---------------------------------------------------------

df2 <- df %>%
  dplyr::mutate(EVENT_ID = as.character(EVENT_ID)) %>%
  dplyr::inner_join(event_map, by = "EVENT_ID") %>%
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

cat("Number of rows after visit and diagnostic-group restriction:", nrow(df2), "\n")

cat("\nRemaining EVENT_ID distribution:\n")
print(table(df2$EVENT_ID, useNA = "ifany"))

cat("\nRemaining GROUP distribution:\n")
print(table(df2$GROUP, useNA = "ifany"))


# ---------------------------------------------------------
# Baseline PCA to define PIS
#
# PIS is derived from baseline neutrophil, lymphocyte, and
# monocyte counts. The baseline PCA solution is then used as
# a common axis to project immune profiles at all visits.
# ---------------------------------------------------------

bl <- df2 %>%
  dplyr::filter(EVENT_ID == "BL") %>%
  dplyr::mutate(
    NLR = safe_ratio(Neutrophils, Lymphocytes),
    MLR = safe_ratio(Monocytes, Lymphocytes)
  ) %>%
  dplyr::select(
    PATNO,
    GROUP,
    EVENT_ID,
    year,
    dplyr::all_of(immune_vars),
    NLR,
    MLR
  )

pca_input_bl <- bl %>%
  dplyr::select(PATNO, dplyr::all_of(immune_vars)) %>%
  tidyr::drop_na()

cat("\nNumber of participants with complete baseline immune data:", nrow(pca_input_bl), "\n")

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
# Orient PC1 for biological interpretability
#
# Goal:
#   Higher PIS should correspond to a more NLR-related,
#   neutrophil-dominant immune profile.
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
# Project all visits into the baseline PCA space
# ---------------------------------------------------------

all_proj <- df2 %>%
  dplyr::mutate(
    NLR = safe_ratio(Neutrophils, Lymphocytes),
    MLR = safe_ratio(Monocytes, Lymphocytes)
  ) %>%
  dplyr::select(
    PATNO,
    GROUP,
    EVENT_ID,
    year,
    dplyr::all_of(site_col),
    dplyr::all_of(immune_vars),
    dplyr::all_of(age_col),
    dplyr::all_of(sex_col),
    dplyr::all_of(bmi_col),
    NLR,
    MLR
  )

proj_complete <- all_proj %>%
  tidyr::drop_na(dplyr::all_of(immune_vars)) %>%
  tidyr::drop_na(dplyr::all_of(c(age_col, sex_col, bmi_col, site_col)))

if (nrow(proj_complete) == 0) {
  stop("No complete observations are available for longitudinal PIS modeling.")
}

scaled_mat <- scale(
  proj_complete[, immune_vars],
  center = pca_fit$center,
  scale = pca_fit$scale
)

proj_complete$PIS <- as.numeric(
  scaled_mat %*% rotation_pc1
)


# ---------------------------------------------------------
# Rename covariates to standard names
# ---------------------------------------------------------

dat <- proj_complete %>%
  dplyr::rename(
    age = dplyr::all_of(age_col),
    sex = dplyr::all_of(sex_col),
    bmi = dplyr::all_of(bmi_col),
    SITE = dplyr::all_of(site_col)
  ) %>%
  dplyr::mutate(
    year_c = year,
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    sex = as.factor(sex),
    SITE = as.factor(SITE),
    GROUP = factor(GROUP, levels = c("HC", "PD"))
  ) %>%
  droplevels()

if (nlevels(dat$GROUP) < 2) {
  stop("The analysis dataset does not contain both HC and PD groups.")
}

if (dplyr::n_distinct(dat$PATNO) < 3) {
  stop("Too few unique participants are available for mixed-effects modeling.")
}


# ---------------------------------------------------------
# Save analysis dataset
# ---------------------------------------------------------

write.csv(
  dat,
  file.path(out_dir, "PIS_longitudinal_dataset.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Descriptive summaries
# ---------------------------------------------------------

desc_n <- dat %>%
  dplyr::group_by(GROUP, EVENT_ID, year) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    mean_PIS = mean(PIS, na.rm = TRUE),
    sd_PIS = stats::sd(PIS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(GROUP, year)

write.csv(
  desc_n,
  file.path(out_dir, "PIS_longitudinal_descriptives.csv"),
  row.names = FALSE
)

site_counts <- dat %>%
  dplyr::group_by(SITE) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(PATNO),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_rows))

write.csv(
  site_counts,
  file.path(out_dir, "SITE_distribution.csv"),
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
  PIS ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO),
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
# Export fixed effects and Type III ANOVA
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
# Estimated marginal means by group and time
# ---------------------------------------------------------

emm <- emmeans::emmeans(
  lmm,
  ~ GROUP | year_c,
  at = list(year_c = time_points)
)

emm_df <- as.data.frame(
  summary(emm, infer = c(TRUE, TRUE))
)

emm_df <- standardize_emmeans_ci(emm_df)

write.csv(
  emm_df,
  file.path(out_dir, "EMMeans_PIS_by_group_and_time.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Group contrasts at each time point
# ---------------------------------------------------------

group_contrasts_each_time <- as.data.frame(
  emmeans::contrast(
    emm,
    method = "pairwise",
    adjust = "bonferroni"
  )
)

write.csv(
  group_contrasts_each_time,
  file.path(out_dir, "Group_comparisons_at_each_timepoint.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Group-specific slopes and slope difference
# ---------------------------------------------------------

emtr <- emmeans::emtrends(
  lmm,
  ~ GROUP,
  var = "year_c"
)

emtr_df <- as.data.frame(emtr)

write.csv(
  emtr_df,
  file.path(out_dir, "Group_specific_slopes.csv"),
  row.names = FALSE
)

slope_contrast <- as.data.frame(
  pairs(
    emtr,
    adjust = "bonferroni"
  )
)

write.csv(
  slope_contrast,
  file.path(out_dir, "Slope_difference_between_groups.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Baseline estimated marginal means and baseline contrast
# ---------------------------------------------------------

emm_baseline <- emmeans::emmeans(
  lmm,
  ~ GROUP,
  at = list(year_c = 0)
)

baseline_emm_df <- as.data.frame(
  summary(emm_baseline, infer = c(TRUE, TRUE))
)

baseline_emm_df <- standardize_emmeans_ci(baseline_emm_df)

baseline_contrast <- as.data.frame(
  pairs(
    emm_baseline,
    adjust = "bonferroni"
  )
)

write.csv(
  baseline_emm_df,
  file.path(out_dir, "Baseline_EMMeans.csv"),
  row.names = FALSE
)

write.csv(
  baseline_contrast,
  file.path(out_dir, "Baseline_group_difference.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Prepare plotting data
# ---------------------------------------------------------

plot_df <- emm_df %>%
  dplyr::mutate(
    EVENT_ID = dplyr::case_when(
      year_c == 0 ~ "BL",
      year_c == 1 ~ "V04",
      year_c == 2 ~ "V06",
      year_c == 3 ~ "V08",
      year_c == 4 ~ "V10",
      year_c == 5 ~ "V12",
      year_c == 6 ~ "V13",
      year_c == 7 ~ "V14",
      TRUE ~ as.character(year_c)
    )
  )

raw_df <- dat %>%
  dplyr::group_by(GROUP, year_c, EVENT_ID) %>%
  dplyr::summarise(
    mean_PIS = mean(PIS, na.rm = TRUE),
    se_PIS = stats::sd(PIS, na.rm = TRUE) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )


# ---------------------------------------------------------
# Longitudinal estimated marginal means plot
# ---------------------------------------------------------

p <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = year_c,
    y = emmean,
    group = GROUP,
    color = GROUP,
    fill = GROUP
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
    alpha = 0.18,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_point(
    data = raw_df,
    ggplot2::aes(
      x = year_c,
      y = mean_PIS,
      group = GROUP,
      color = GROUP
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    stroke = 1,
    size = 2.2
  ) +
  ggplot2::scale_x_continuous(
    breaks = event_map$year,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  ggplot2::labs(
    title = "Longitudinal trajectory of Peripheral Immune Score (PIS)",
    subtitle = "Estimated marginal means from the mixed-effects model with 95% CI",
    x = "Years",
    y = "Peripheral Immune Score (PIS)",
    color = "Group",
    fill = "Group"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "PIS_longitudinal_plot.png"),
  plot = p,
  width = 8,
  height = 5.5,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_dir, "PIS_longitudinal_plot.pdf"),
  plot = p,
  width = 8,
  height = 5.5
)


# ---------------------------------------------------------
# Optional spaghetti plot of observed values
# ---------------------------------------------------------

p_spaghetti <- ggplot2::ggplot(
  dat,
  ggplot2::aes(
    x = year_c,
    y = PIS,
    group = PATNO,
    color = GROUP
  )
) +
  ggplot2::geom_line(alpha = 0.08) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = GROUP,
      fill = GROUP
    ),
    method = "loess",
    se = TRUE,
    linewidth = 1.2,
    alpha = 0.18
  ) +
  ggplot2::scale_x_continuous(
    breaks = event_map$year,
    labels = c("BL", "1", "2", "3", "4", "5", "6", "7")
  ) +
  ggplot2::labs(
    title = "Observed longitudinal trajectories of PIS",
    subtitle = "Thin lines indicate individual participants; thick lines indicate smoothed group trends",
    x = "Years",
    y = "Peripheral Immune Score (PIS)",
    color = "Group",
    fill = "Group"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "PIS_spaghetti_plot.png"),
  plot = p_spaghetti,
  width = 8,
  height = 5.5,
  dpi = 600
)


# ---------------------------------------------------------
# Save full text summary
# ---------------------------------------------------------

sink(file.path(out_dir, "LMM_summary.txt"))

cat("============================================================\n")
cat("Longitudinal analysis of Peripheral Immune Score (PIS)\n")
cat("============================================================\n\n")

cat("Model formula:\n")
cat("PIS ~ GROUP * year_c + age + sex + bmi + (1 | SITE) + (1 | PATNO)\n\n")

cat("Detected covariate columns:\n")
cat("age  =", age_col, "\n")
cat("sex  =", sex_col, "\n")
cat("bmi  =", bmi_col, "\n")
cat("site =", site_col, "\n\n")

cat("Sample size:\n")
cat("Rows =", nrow(dat), "\n")
cat("Unique subjects =", dplyr::n_distinct(dat$PATNO), "\n")
cat("Unique sites =", dplyr::n_distinct(dat$SITE), "\n\n")

cat("Group counts:\n")
print(table(dat$GROUP))
cat("\n")

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

cat("Estimated marginal means by group and time:\n")
print(emm_df, row.names = FALSE)
cat("\n\n")

cat("Group comparisons at each time point:\n")
print(group_contrasts_each_time, row.names = FALSE)
cat("\n\n")

cat("Group-specific slopes:\n")
print(emtr_df, row.names = FALSE)
cat("\n\n")

cat("Slope difference between groups:\n")
print(slope_contrast, row.names = FALSE)
cat("\n\n")

cat("Baseline group difference:\n")
print(baseline_contrast, row.names = FALSE)
cat("\n")

sink()


# ---------------------------------------------------------
# Restore contrast options
# ---------------------------------------------------------

options(old_contrasts)


# ---------------------------------------------------------
# Console output
# ---------------------------------------------------------

cat("\n============================================================\n")
cat("Longitudinal PIS mixed-effects analysis completed.\n")
cat("Outputs saved to:\n")
cat(out_dir, "\n")
cat("============================================================\n")