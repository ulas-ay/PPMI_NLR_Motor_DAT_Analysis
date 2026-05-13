# =========================================================
# Striatal DAT Moderation Models
#
# This script:
#   1. Restricts the dataset to Parkinson's disease participants
#   2. Reconstructs baseline PCA-derived Peripheral Immune Score (PIS)
#      using neutrophil, lymphocyte, and monocyte counts
#   3. Extracts baseline striatal DAT binding measures
#   4. Tests whether baseline striatal DAT binding moderates the
#      association between baseline PIS and longitudinal motor progression
#   5. Runs separate moderation models for:
#        - Left caudate DAT binding
#        - Right caudate DAT binding
#        - Left putaminal DAT binding
#        - Right putaminal DAT binding
#   6. Exports model summaries, simple slopes, predicted trajectories,
#      raw descriptive plots, and combined moderation summary tables
#
# Required input:
#   data/PPMI_analysis_dataset.xlsx
#
# Main model:
#   UPDRS ~ PIS_BL_z * DAT_BL_z * year_c +
#           age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
#           (1 | SITE) + (1 | PATNO)
#
# Outputs:
#   outputs/striatal_dat_moderation/
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
out_dir   <- file.path("outputs", "striatal_dat_moderation")

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
#
# Longitudinal motor trajectories include visits from baseline
# through V19 when available.
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
  "V14",  7,
)


# ---------------------------------------------------------
# Detect shared columns robustly
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

caudate_l_col <- find_col(
  df,
  c("MIA_CAUDATE_L")
)

caudate_r_col <- find_col(
  df,
  c("MIA_CAUDATE_R")
)

putamen_l_col <- find_col(
  df,
  c("MIA_PUTAMEN_L")
)

putamen_r_col <- find_col(
  df,
  c("MIA_PUTAMEN_R")
)

cat("Detected shared columns:\n")
cat("age_col       =", age_col, "\n")
cat("sex_col       =", sex_col, "\n")
cat("bmi_col       =", bmi_col, "\n")
cat("site_col      =", site_col, "\n")
cat("duration_col  =", duration_col, "\n")
cat("updrs_col     =", updrs_col, "\n")
cat("ledd_col      =", ledd_col, "\n")
cat("domside_col   =", domside_col, "\n")
cat("caudate_l_col =", caudate_l_col, "\n")
cat("caudate_r_col =", caudate_r_col, "\n")
cat("putamen_l_col =", putamen_l_col, "\n")
cat("putamen_r_col =", putamen_r_col, "\n\n")

detected_cols <- c(
  age_col,
  sex_col,
  bmi_col,
  site_col,
  duration_col,
  updrs_col,
  ledd_col,
  domside_col,
  caudate_l_col,
  caudate_r_col,
  putamen_l_col,
  putamen_r_col
)

if (any(is.na(detected_cols))) {
  stop(
    paste0(
      "Could not detect one or more required columns. ",
      "Please check names(df) and update the candidate names in this script."
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
  dplyr::filter(PRIMDIAG == 1)

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

rotation_pc1 <- if (!is.na(tmp_cor_nlr) && tmp_cor_nlr < 0) {
  -pca_fit$rotation[, 1]
} else {
  pca_fit$rotation[, 1]
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

write.csv(
  bl_pis,
  file.path(out_dir, "Baseline_PIS_PD.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Shared UPDRS longitudinal dataset
#
# Baseline PIS is carried forward to all visits.
# Baseline LEDD is set to 0 to represent the drug-naive
# baseline state when applicable.
# ---------------------------------------------------------

dat_long <- df2 %>%
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
    dplyr::everything()
  ) %>%
  dplyr::rename(
    SITE = dplyr::all_of(site_col),
    age = dplyr::all_of(age_col),
    sex = dplyr::all_of(sex_col),
    bmi = dplyr::all_of(bmi_col),
    duration_yrs = dplyr::all_of(duration_col),
    LEDD_raw = dplyr::all_of(ledd_col),
    DOMSIDE = dplyr::all_of(domside_col),
    UPDRS = dplyr::all_of(updrs_col)
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
    UPDRS = as.numeric(UPDRS),
    LEDD = ifelse(year_c == 0, 0, LEDD_raw)
  ) %>%
  tidyr::drop_na(
    PIS_BL,
    age,
    sex,
    bmi,
    duration_yrs,
    DOMSIDE,
    UPDRS,
    SITE
  ) %>%
  droplevels()

if (dplyr::n_distinct(dat_long$PATNO) < 3) {
  stop("Too few unique PD participants are available for moderation modeling.")
}

write.csv(
  dat_long,
  file.path(out_dir, "Shared_longitudinal_UPDRS_dataset.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Set contrasts for Type III tests
# ---------------------------------------------------------

old_contrasts <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))


# ---------------------------------------------------------
# Main moderation function
# ---------------------------------------------------------

run_striatal_moderation <- function(dat_var, dat_label, file_stub = dat_var) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Running striatal DAT moderation model for:", dat_label, "\n")
  cat("DAT moderator:", dat_var, "\n")
  cat("------------------------------------------------------------\n")
  
  dat_col <- find_col(df2, c(dat_var))
  
  if (is.na(dat_col)) {
    stop(paste("Could not find DAT column:", dat_var))
  }
  
  
  # -------------------------------------------------------
  # Extract baseline DAT
  # -------------------------------------------------------
  
  bl_dat <- df2 %>%
    dplyr::filter(EVENT_ID == "BL") %>%
    dplyr::select(PATNO, dplyr::all_of(dat_col)) %>%
    dplyr::rename(DAT_BL = dplyr::all_of(dat_col)) %>%
    dplyr::mutate(DAT_BL = as.numeric(DAT_BL)) %>%
    dplyr::distinct(PATNO, .keep_all = TRUE)
  
  write.csv(
    bl_dat,
    file.path(out_dir, paste0(file_stub, "_baseline_DAT.csv")),
    row.names = FALSE
  )
  
  
  # -------------------------------------------------------
  # Merge baseline DAT into longitudinal UPDRS dataset
  # -------------------------------------------------------
  
  dat <- dat_long %>%
    dplyr::left_join(bl_dat, by = "PATNO") %>%
    tidyr::drop_na(DAT_BL) %>%
    dplyr::mutate(
      PIS_BL_z = as.numeric(scale(PIS_BL)),
      DAT_BL_z = as.numeric(scale(DAT_BL))
    ) %>%
    droplevels()
  
  if (dplyr::n_distinct(dat$PATNO) < 3) {
    stop(
      paste0(
        "Too few unique participants for moderation model: ",
        dat_label
      )
    )
  }
  
  write.csv(
    dat,
    file.path(out_dir, paste0(file_stub, "_moderation_dataset.csv")),
    row.names = FALSE
  )
  
  
  # -------------------------------------------------------
  # Descriptive statistics
  # -------------------------------------------------------
  
  desc <- dat %>%
    dplyr::group_by(EVENT_ID, year) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_subjects = dplyr::n_distinct(PATNO),
      mean_UPDRS = mean(UPDRS, na.rm = TRUE),
      sd_UPDRS = stats::sd(UPDRS, na.rm = TRUE),
      mean_PIS_BL = mean(PIS_BL, na.rm = TRUE),
      sd_PIS_BL = stats::sd(PIS_BL, na.rm = TRUE),
      mean_DAT_BL = mean(DAT_BL, na.rm = TRUE),
      sd_DAT_BL = stats::sd(DAT_BL, na.rm = TRUE),
      mean_LEDD = mean(LEDD, na.rm = TRUE),
      sd_LEDD = stats::sd(LEDD, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(year)
  
  write.csv(
    desc,
    file.path(out_dir, paste0(file_stub, "_descriptives.csv")),
    row.names = FALSE
  )
  
  
  # -------------------------------------------------------
  # Linear mixed-effects moderation model
  # -------------------------------------------------------
  
  lmm <- lmerTest::lmer(
    UPDRS ~ PIS_BL_z * DAT_BL_z * year_c +
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
  
  
  # -------------------------------------------------------
  # Export model outputs
  # -------------------------------------------------------
  
  fixef_tab <- as.data.frame(coef(summary(lmm)))
  fixef_tab$term <- rownames(fixef_tab)
  
  fixef_tab <- fixef_tab %>%
    dplyr::relocate(term)
  
  write.csv(
    fixef_tab,
    file.path(out_dir, paste0(file_stub, "_LMM_fixed_effects.csv")),
    row.names = FALSE
  )
  
  anova_tab <- as.data.frame(anova_type3)
  anova_tab$term <- rownames(anova_tab)
  
  anova_tab <- anova_tab %>%
    dplyr::relocate(term)
  
  write.csv(
    anova_tab,
    file.path(out_dir, paste0(file_stub, "_LMM_TypeIII_ANOVA.csv")),
    row.names = FALSE
  )
  
  randef_var <- as.data.frame(lme4::VarCorr(lmm))
  
  write.csv(
    randef_var,
    file.path(out_dir, paste0(file_stub, "_LMM_random_effects_variance.csv")),
    row.names = FALSE
  )
  
  
  # -------------------------------------------------------
  # Simple slopes of PIS at each time point by DAT level
  #
  # DAT levels use standardized values:
  #   -1.5 SD, mean, +1.5 SD
  # -------------------------------------------------------
  
  time_points <- sort(unique(dat$year_c))
  dat_levels <- c(-1.5, 0, 1.5)
  
  slopes_pis_by_dat_time <- emmeans::emtrends(
    lmm,
    ~ year_c | DAT_BL_z,
    var = "PIS_BL_z",
    at = list(
      year_c = time_points,
      DAT_BL_z = dat_levels
    )
  )
  
  slopes_pis_by_dat_time_df <- as.data.frame(
    summary(slopes_pis_by_dat_time, infer = c(TRUE, TRUE))
  )
  
  slopes_pis_by_dat_time_df <- standardize_ci_names(slopes_pis_by_dat_time_df)
  
  slopes_pis_by_dat_time_df <- slopes_pis_by_dat_time_df %>%
    dplyr::mutate(
      DAT_group = factor(
        DAT_BL_z,
        levels = c(-1.5, 0, 1.5),
        labels = c(
          "Low baseline DAT (-1.5 SD)",
          "Mean baseline DAT",
          "High baseline DAT (+1.5 SD)"
        )
      )
    )
  
  write.csv(
    slopes_pis_by_dat_time_df,
    file.path(out_dir, paste0(file_stub, "_simple_slopes_PIS_by_DAT_and_time.csv")),
    row.names = FALSE
  )
  
  
  # -------------------------------------------------------
  # Predicted trajectories at PIS and DAT values:
  #   PIS = -1.5 SD, mean, +1.5 SD
  #   DAT = -1.5 SD, mean, +1.5 SD
  # -------------------------------------------------------
  
  emm_grid <- emmeans::emmeans(
    lmm,
    ~ year_c | PIS_BL_z * DAT_BL_z,
    at = list(
      year_c = time_points,
      PIS_BL_z = c(-1.5, 0, 1.5),
      DAT_BL_z = c(-1.5, 0, 1.5),
      age = mean(dat$age, na.rm = TRUE),
      bmi = mean(dat$bmi, na.rm = TRUE),
      duration_yrs = mean(dat$duration_yrs, na.rm = TRUE),
      LEDD = mean(dat$LEDD, na.rm = TRUE)
    )
  )
  
  emm_df <- as.data.frame(
    summary(emm_grid, infer = c(TRUE, TRUE))
  )
  
  emm_df <- standardize_ci_names(emm_df)
  
  emm_df <- emm_df %>%
    dplyr::mutate(
      PIS_group = factor(
        PIS_BL_z,
        levels = c(-1.5, 0, 1.5),
        labels = c(
          "Low PIS (-1.5 SD)",
          "Mean PIS",
          "High PIS (+1.5 SD)"
        )
      ),
      DAT_group = factor(
        DAT_BL_z,
        levels = c(-1.5, 0, 1.5),
        labels = c(
          "Low DAT (-1.5 SD)",
          "Mean DAT",
          "High DAT (+1.5 SD)"
        )
      )
    )
  
  write.csv(
    emm_df,
    file.path(out_dir, paste0(file_stub, "_predicted_UPDRS_trajectories_1p5SD.csv")),
    row.names = FALSE
  )
  
  
  # -------------------------------------------------------
  # Figure 1: model-estimated moderation trajectories
  # -------------------------------------------------------
  
  p1 <- ggplot2::ggplot(
    emm_df,
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
      alpha = 0.16,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::facet_wrap(~ DAT_group) +
    ggplot2::scale_x_continuous(
      breaks = time_points,
      labels = ifelse(time_points == 0, "BL", as.character(time_points))
    ) +
    ggplot2::scale_color_discrete(drop = TRUE, na.translate = FALSE) +
    ggplot2::scale_fill_discrete(drop = TRUE, na.translate = FALSE) +
    ggplot2::labs(
      title = paste0(
        dat_label,
        " moderates the association between baseline PIS and motor progression"
      ),
      subtitle = paste0(
        "Model-estimated UPDRS-III trajectories at ±1.5 SD levels ",
        "of baseline PIS and DAT"
      ),
      x = "Years",
      y = "Predicted UPDRS-III score",
      color = "Baseline PIS",
      fill = "Baseline PIS"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold")
    )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure1_moderation_1p5SD.png")),
    plot = p1,
    width = 11,
    height = 5.8,
    dpi = 600
  )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure1_moderation_1p5SD.pdf")),
    plot = p1,
    width = 11,
    height = 5.8
  )
  
  
  # -------------------------------------------------------
  # Figure 2: tertile-based raw descriptive plot
  #
  # Tertiles are used for descriptive visualization only.
  # The primary moderation model uses standardized continuous
  # PIS and DAT values.
  # -------------------------------------------------------
  
  dat_plot <- dat %>%
    dplyr::distinct(PATNO, .keep_all = TRUE) %>%
    dplyr::mutate(
      PIS_tertile = dplyr::ntile(PIS_BL, 3),
      DAT_tertile = dplyr::ntile(DAT_BL, 3)
    ) %>%
    dplyr::select(PATNO, PIS_tertile, DAT_tertile) %>%
    dplyr::mutate(
      PIS_tertile = factor(
        PIS_tertile,
        levels = c(1, 2, 3),
        labels = c(
          "Low PIS tertile",
          "Middle PIS tertile",
          "High PIS tertile"
        )
      ),
      DAT_tertile = factor(
        DAT_tertile,
        levels = c(1, 2, 3),
        labels = c(
          "Low DAT tertile",
          "Middle DAT tertile",
          "High DAT tertile"
        )
      )
    )
  
  raw_df <- dat %>%
    dplyr::left_join(dat_plot, by = "PATNO") %>%
    dplyr::group_by(DAT_tertile, PIS_tertile, year_c, EVENT_ID) %>%
    dplyr::summarise(
      mean_UPDRS = mean(UPDRS, na.rm = TRUE),
      se_UPDRS = stats::sd(UPDRS, na.rm = TRUE) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  write.csv(
    raw_df,
    file.path(out_dir, paste0(file_stub, "_raw_descriptive_UPDRS_by_PIS_DAT_tertiles.csv")),
    row.names = FALSE
  )
  
  p2 <- ggplot2::ggplot(
    raw_df,
    ggplot2::aes(
      x = year_c,
      y = mean_UPDRS,
      color = PIS_tertile,
      group = PIS_tertile
    )
  ) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_UPDRS - se_UPDRS,
        ymax = mean_UPDRS + se_UPDRS
      ),
      width = 0.12
    ) +
    ggplot2::facet_wrap(~ DAT_tertile) +
    ggplot2::scale_x_continuous(
      breaks = time_points,
      labels = ifelse(time_points == 0, "BL", as.character(time_points))
    ) +
    ggplot2::scale_color_discrete(drop = TRUE, na.translate = FALSE) +
    ggplot2::labs(
      title = paste0(
        "Observed motor trajectories by baseline PIS and ",
        dat_label,
        " tertiles"
      ),
      subtitle = "Raw descriptive means",
      x = "Years",
      y = "Observed mean UPDRS-III",
      color = "Baseline PIS tertile"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure2_moderation_raw.png")),
    plot = p2,
    width = 11,
    height = 5.8,
    dpi = 600
  )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(file_stub, "_Figure2_moderation_raw.pdf")),
    plot = p2,
    width = 11,
    height = 5.8
  )
  
  
  # -------------------------------------------------------
  # Save text summary
  # -------------------------------------------------------
  
  sink(file.path(out_dir, paste0(file_stub, "_moderation_summary.txt")))
  
  cat("============================================================\n")
  cat("Striatal DAT moderation analysis\n")
  cat("============================================================\n\n")
  
  cat("DAT moderator:\n")
  cat(dat_var, "(", dat_label, ")\n\n")
  
  cat("Model formula:\n")
  cat(
    paste0(
      "UPDRS ~ PIS_BL_z * DAT_BL_z * year_c + age + sex + bmi + ",
      "duration_yrs + LEDD + DOMSIDE + (1 | SITE) + (1 | PATNO)\n\n"
    )
  )
  
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
  
  cat("Simple slopes of PIS at each time point, stratified by baseline DAT ±1.5 SD:\n")
  print(slopes_pis_by_dat_time_df, row.names = FALSE)
  cat("\n\n")
  
  sink()
  
  
  # -------------------------------------------------------
  # Return model objects
  # -------------------------------------------------------
  
  invisible(
    list(
      model = lmm,
      summary = model_summary,
      anova = anova_type3,
      fixef = fixef_tab,
      slopes = slopes_pis_by_dat_time_df,
      descriptives = desc
    )
  )
}


# ---------------------------------------------------------
# Run all striatal DAT moderation models
# ---------------------------------------------------------

res_caudate_l <- run_striatal_moderation(
  dat_var = "MIA_CAUDATE_L",
  dat_label = "Left caudate DAT",
  file_stub = "CAUDATE_L"
)

res_caudate_r <- run_striatal_moderation(
  dat_var = "MIA_CAUDATE_R",
  dat_label = "Right caudate DAT",
  file_stub = "CAUDATE_R"
)

res_putamen_l <- run_striatal_moderation(
  dat_var = "MIA_PUTAMEN_L",
  dat_label = "Left putaminal DAT",
  file_stub = "PUTAMEN_L"
)

res_putamen_r <- run_striatal_moderation(
  dat_var = "MIA_PUTAMEN_R",
  dat_label = "Right putaminal DAT",
  file_stub = "PUTAMEN_R"
)


# ---------------------------------------------------------
# Combined overview tables across moderation models
# ---------------------------------------------------------

combined_type3 <- dplyr::bind_rows(
  as.data.frame(res_caudate_l$anova) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(moderator = "CAUDATE_L"),
  
  as.data.frame(res_caudate_r$anova) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(moderator = "CAUDATE_R"),
  
  as.data.frame(res_putamen_l$anova) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(moderator = "PUTAMEN_L"),
  
  as.data.frame(res_putamen_r$anova) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(moderator = "PUTAMEN_R")
) %>%
  dplyr::select(moderator, term, dplyr::everything())

write.csv(
  combined_type3,
  file.path(out_dir, "Combined_moderation_TypeIII_ANOVA.csv"),
  row.names = FALSE
)

combined_slopes <- dplyr::bind_rows(
  res_caudate_l$slopes %>% dplyr::mutate(moderator = "CAUDATE_L"),
  res_caudate_r$slopes %>% dplyr::mutate(moderator = "CAUDATE_R"),
  res_putamen_l$slopes %>% dplyr::mutate(moderator = "PUTAMEN_L"),
  res_putamen_r$slopes %>% dplyr::mutate(moderator = "PUTAMEN_R")
) %>%
  dplyr::select(moderator, dplyr::everything())

write.csv(
  combined_slopes,
  file.path(out_dir, "Combined_simple_slopes_PIS_by_DAT_and_time.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Restore contrast options
# ---------------------------------------------------------

options(old_contrasts)


# ---------------------------------------------------------
# Console output
# ---------------------------------------------------------

cat("\n============================================================\n")
cat("Striatal DAT moderation analyses completed.\n")
cat("Outputs saved to:\n")
cat(out_dir, "\n")
cat("============================================================\n")