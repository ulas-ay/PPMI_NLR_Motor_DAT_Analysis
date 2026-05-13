# =========================================================
# Baseline Immune PCA and Longitudinal Projection
#
# This script:
#   1. Restricts the dataset to Parkinson's disease participants
#      and healthy controls based on PRIMDIAG coding
#   2. Creates a baseline immune score using PCA on peripheral
#      immune cell counts
#   3. Aligns the first principal component so that higher values
#      indicate higher NLR-related immune activation
#   4. Evaluates correlations between the immune score and NLR/MLR
#   5. Performs visitwise PCA as a sensitivity/stability check
#   6. Projects all available visits onto the baseline PCA axis
#      to derive longitudinal immune scores
#
# Required input:
#   data/PPMI_analysis_dataset.xlsx
#
# Outputs:
#   outputs/pca/BL_PCA_variance_explained.csv
#   outputs/pca/BL_PCA_loadings_PC1.csv
#   outputs/pca/BL_subject_immune_scores.csv
#   outputs/pca/BL_PCA_scree_plot.png
#   outputs/pca/BL_immuneScore_correlations_with_NLR_MLR.csv
#   outputs/pca/BL_immuneScore_vs_NLR.png
#   outputs/pca/BL_immuneScore_vs_MLR.png
#   outputs/pca/Visitwise_PCA_variance_explained.csv
#   outputs/pca/Visitwise_PCA_loadings_PC1.csv
#   outputs/pca/AllVisits_projected_immune_scores_from_BL_PCA.csv
#   outputs/pca/Visitwise_correlations_projectedImmuneScore_with_NLR_MLR.csv
#   outputs/pca/AllVisits_projected_immune_score_by_group.png
#   outputs/pca/PCA_summary.txt
# =========================================================


# ---------------------------------------------------------
# Load required packages
# ---------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readxl",
  "tidyr",
  "ggplot2",
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
out_dir   <- file.path("outputs", "pca")

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


safe_cor_test <- function(data, x, y, method = "spearman") {
  dd <- data %>%
    dplyr::select(dplyr::all_of(c(x, y))) %>%
    tidyr::drop_na()
  
  if (nrow(dd) < 3) {
    return(
      tibble::tibble(
        n = nrow(dd),
        rho = NA_real_,
        p = NA_real_
      )
    )
  }
  
  ct <- suppressWarnings(
    cor.test(
      dd[[x]],
      dd[[y]],
      method = method,
      exact = FALSE
    )
  )
  
  tibble::tibble(
    n = nrow(dd),
    rho = unname(ct$estimate),
    p = ct$p.value
  )
}


make_scatter <- function(data, yvar, ylab, out_name) {
  ann <- safe_cor_test(data, "immune_score", yvar)
  
  ann_lab <- paste0(
    "Spearman rho = ",
    ifelse(is.na(ann$rho), "NA", sprintf("%.3f", ann$rho)),
    "\np = ",
    ifelse(is.na(ann$p), "NA", formatC(ann$p, format = "e", digits = 2)),
    "\nn = ",
    ann$n
  )
  
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = immune_score,
      y = .data[[yvar]],
      color = GROUP
    )
  ) +
    ggplot2::geom_point(size = 2.3, alpha = 0.8) +
    ggplot2::geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title = paste0("Baseline immune score vs ", yvar),
      x = "Immune score, PC1",
      y = ylab,
      color = "Group"
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      label = ann_lab,
      hjust = 1.05,
      vjust = 1.1,
      size = 4
    )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, out_name),
    plot = p,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  return(p)
}


# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------

df <- readxl::read_excel(data_path)

cat("\nData successfully loaded.\n")
cat("Number of rows:", nrow(df), "\n")
cat("Number of columns:", ncol(df), "\n\n")


# ---------------------------------------------------------
# Check required columns
# ---------------------------------------------------------

required_vars <- c(
  "PATNO",
  "EVENT_ID",
  "PRIMDIAG",
  "Neutrophils",
  "Lymphocytes",
  "Monocytes"
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
# Define visit-year mapping
#
# The selected visits correspond to the longitudinal follow-up
# used for the immune PCA projection analysis.
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
# Create baseline dataset
# ---------------------------------------------------------

immune_vars <- c(
  "Neutrophils",
  "Lymphocytes",
  "Monocytes"
)

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
    MLR,
    dplyr::everything()
  )

cat("\nNumber of baseline rows:", nrow(bl), "\n")


# ---------------------------------------------------------
# Prepare baseline PCA input
# ---------------------------------------------------------

pca_input_bl <- bl %>%
  dplyr::select(PATNO, dplyr::all_of(immune_vars)) %>%
  tidyr::drop_na()

cat("Number of participants with complete baseline immune data:", nrow(pca_input_bl), "\n")

if (nrow(pca_input_bl) < 5) {
  stop(
    paste0(
      "Insufficient complete baseline observations for PCA. ",
      "At least 5 complete rows are recommended."
    )
  )
}


# ---------------------------------------------------------
# Run baseline PCA
# ---------------------------------------------------------

pca_fit <- stats::prcomp(
  pca_input_bl[, immune_vars],
  scale. = TRUE
)


# ---------------------------------------------------------
# Extract subject-level raw PC1 scores
# ---------------------------------------------------------

bl_scores_raw <- tibble::tibble(
  PATNO = pca_input_bl$PATNO,
  immune_score_raw = as.numeric(pca_fit$x[, 1])
)

bl_scores <- bl %>%
  dplyr::left_join(bl_scores_raw, by = "PATNO")


# ---------------------------------------------------------
# Orient PC1 for biological interpretability
#
# Goal:
#   Higher immune_score should correspond to higher NLR.
#   If the raw PC1 is negatively correlated with NLR, the sign
#   of PC1 and its loadings is reversed.
# ---------------------------------------------------------

tmp_cor_nlr <- suppressWarnings(
  stats::cor(
    bl_scores$immune_score_raw,
    bl_scores$NLR,
    use = "pairwise.complete.obs",
    method = "spearman"
  )
)

if (!is.na(tmp_cor_nlr) && tmp_cor_nlr < 0) {
  bl_scores <- bl_scores %>%
    dplyr::mutate(immune_score = -immune_score_raw)
  
  rotation_pc1 <- -pca_fit$rotation[, 1]
  
} else {
  bl_scores <- bl_scores %>%
    dplyr::mutate(immune_score = immune_score_raw)
  
  rotation_pc1 <- pca_fit$rotation[, 1]
}

cat("\nSpearman correlation between oriented immune_score and NLR:\n")
print(
  suppressWarnings(
    stats::cor(
      bl_scores$immune_score,
      bl_scores$NLR,
      use = "pairwise.complete.obs",
      method = "spearman"
    )
  )
)


# ---------------------------------------------------------
# PCA summary tables
# ---------------------------------------------------------

eig <- pca_fit$sdev^2
var_explained <- eig / sum(eig)
cum_var <- cumsum(var_explained)

pca_variance_table <- tibble::tibble(
  PC = paste0("PC", seq_along(eig)),
  eigenvalue = eig,
  proportion_variance = var_explained,
  cumulative_variance = cum_var
)

pca_loadings_table <- tibble::tibble(
  variable = rownames(pca_fit$rotation),
  PC1_loading = as.numeric(rotation_pc1),
  abs_PC1_loading = abs(as.numeric(rotation_pc1))
) %>%
  dplyr::arrange(dplyr::desc(abs_PC1_loading))

subject_scores_table <- bl_scores %>%
  dplyr::select(
    PATNO,
    GROUP,
    immune_score,
    immune_score_raw,
    NLR,
    MLR,
    Neutrophils,
    Lymphocytes,
    Monocytes
  )


# ---------------------------------------------------------
# Export baseline PCA tables
# ---------------------------------------------------------

write.csv(
  pca_variance_table,
  file.path(out_dir, "BL_PCA_variance_explained.csv"),
  row.names = FALSE
)

write.csv(
  pca_loadings_table,
  file.path(out_dir, "BL_PCA_loadings_PC1.csv"),
  row.names = FALSE
)

write.csv(
  subject_scores_table,
  file.path(out_dir, "BL_subject_immune_scores.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Scree plot
# ---------------------------------------------------------

scree_df <- pca_variance_table

p_scree <- ggplot2::ggplot(
  scree_df,
  ggplot2::aes(x = PC, y = proportion_variance)
) +
  ggplot2::geom_col() +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", 100 * proportion_variance)),
    vjust = -0.3,
    size = 4
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::labs(
    title = "Baseline PCA: explained variance",
    x = "Principal component",
    y = "Proportion of variance explained"
  ) +
  ggplot2::ylim(
    0,
    max(scree_df$proportion_variance, na.rm = TRUE) * 1.15
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "BL_PCA_scree_plot.png"),
  plot = p_scree,
  width = 7,
  height = 5,
  dpi = 600
)


# ---------------------------------------------------------
# Correlations: immune_score vs NLR / MLR
# ---------------------------------------------------------

cor_overall_nlr <- safe_cor_test(
  bl_scores,
  "immune_score",
  "NLR"
) %>%
  dplyr::mutate(group = "Overall", marker = "NLR")

cor_overall_mlr <- safe_cor_test(
  bl_scores,
  "immune_score",
  "MLR"
) %>%
  dplyr::mutate(group = "Overall", marker = "MLR")

cor_by_group <- dplyr::bind_rows(
  bl_scores %>%
    dplyr::group_by(GROUP) %>%
    dplyr::group_modify(~ safe_cor_test(.x, "immune_score", "NLR")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(marker = "NLR"),
  
  bl_scores %>%
    dplyr::group_by(GROUP) %>%
    dplyr::group_modify(~ safe_cor_test(.x, "immune_score", "MLR")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(marker = "MLR")
) %>%
  dplyr::rename(group = GROUP)

cor_table <- dplyr::bind_rows(
  cor_overall_nlr,
  cor_overall_mlr,
  cor_by_group
) %>%
  dplyr::select(group, marker, n, rho, p)

write.csv(
  cor_table,
  file.path(out_dir, "BL_immuneScore_correlations_with_NLR_MLR.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Scatter plots
# ---------------------------------------------------------

p_nlr <- make_scatter(
  data = bl_scores,
  yvar = "NLR",
  ylab = "Neutrophil-to-lymphocyte ratio, NLR",
  out_name = "BL_immuneScore_vs_NLR.png"
)

p_mlr <- make_scatter(
  data = bl_scores,
  yvar = "MLR",
  ylab = "Monocyte-to-lymphocyte ratio, MLR",
  out_name = "BL_immuneScore_vs_MLR.png"
)


# ---------------------------------------------------------
# Visitwise PCA stability check
#
# Purpose:
#   Run PCA separately at each visit to evaluate whether the
#   structure of PC1 is broadly stable across follow-up visits.
# ---------------------------------------------------------

visit_list <- event_map$EVENT_ID

visit_pca_results <- list()
visit_loading_results <- list()

for (ev in visit_list) {
  
  d_ev <- df2 %>%
    dplyr::filter(EVENT_ID == ev) %>%
    dplyr::mutate(
      NLR = safe_ratio(Neutrophils, Lymphocytes),
      MLR = safe_ratio(Monocytes, Lymphocytes)
    ) %>%
    dplyr::select(
      PATNO,
      GROUP,
      EVENT_ID,
      dplyr::all_of(immune_vars),
      NLR,
      MLR
    ) %>%
    tidyr::drop_na(dplyr::all_of(immune_vars))
  
  if (nrow(d_ev) < 5) {
    cat(
      "\n",
      ev,
      ": insufficient complete observations for PCA; skipped. n =",
      nrow(d_ev),
      "\n"
    )
    next
  }
  
  fit_ev <- stats::prcomp(
    d_ev[, immune_vars],
    scale. = TRUE
  )
  
  eig_ev <- fit_ev$sdev^2
  var_ev <- eig_ev / sum(eig_ev)
  
  load_ev <- fit_ev$rotation[, 1]
  
  # Align visitwise PC1 direction with baseline PC1
  align_sign <- sign(sum(load_ev * rotation_pc1))
  
  if (is.na(align_sign) || align_sign == 0) {
    align_sign <- 1
  }
  
  load_ev <- load_ev * align_sign
  
  visit_pca_results[[ev]] <- tibble::tibble(
    EVENT_ID = ev,
    n_complete = nrow(d_ev),
    PC = paste0("PC", seq_along(eig_ev)),
    eigenvalue = eig_ev,
    proportion_variance = var_ev,
    cumulative_variance = cumsum(var_ev)
  )
  
  visit_loading_results[[ev]] <- tibble::tibble(
    EVENT_ID = ev,
    variable = names(load_ev),
    PC1_loading = as.numeric(load_ev)
  )
}

visit_variance_table <- dplyr::bind_rows(visit_pca_results)
visit_loadings_table <- dplyr::bind_rows(visit_loading_results)

write.csv(
  visit_variance_table,
  file.path(out_dir, "Visitwise_PCA_variance_explained.csv"),
  row.names = FALSE
)

write.csv(
  visit_loadings_table,
  file.path(out_dir, "Visitwise_PCA_loadings_PC1.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Project all visits onto the baseline PCA axis
#
# This creates a longitudinal immune score on a common PCA
# scale derived from the baseline PCA solution.
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
    dplyr::all_of(immune_vars),
    NLR,
    MLR
  )

proj_complete <- all_proj %>%
  tidyr::drop_na(dplyr::all_of(immune_vars))

if (nrow(proj_complete) == 0) {
  stop("No complete observations are available for all-visit PCA projection.")
}


# ---------------------------------------------------------
# Standardize using baseline PCA center and scale values
# ---------------------------------------------------------

scaled_mat <- scale(
  proj_complete[, immune_vars],
  center = pca_fit$center,
  scale = pca_fit$scale
)

proj_complete$immune_score_projected <- as.numeric(
  scaled_mat %*% rotation_pc1
)

write.csv(
  proj_complete,
  file.path(out_dir, "AllVisits_projected_immune_scores_from_BL_PCA.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Visitwise correlations of projected immune score with NLR/MLR
# ---------------------------------------------------------

visit_corrs <- dplyr::bind_rows(
  proj_complete %>%
    dplyr::group_by(EVENT_ID) %>%
    dplyr::group_modify(
      ~ safe_cor_test(.x, "immune_score_projected", "NLR")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(marker = "NLR"),
  
  proj_complete %>%
    dplyr::group_by(EVENT_ID) %>%
    dplyr::group_modify(
      ~ safe_cor_test(.x, "immune_score_projected", "MLR")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(marker = "MLR")
)

write.csv(
  visit_corrs,
  file.path(out_dir, "Visitwise_correlations_projectedImmuneScore_with_NLR_MLR.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# Longitudinal plot of projected immune score
# ---------------------------------------------------------

p_proj <- ggplot2::ggplot(
  proj_complete,
  ggplot2::aes(
    x = year,
    y = immune_score_projected,
    color = GROUP
  )
) +
  ggplot2::geom_point(alpha = 0.35, size = 1.8) +
  ggplot2::geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(proj_complete$year)),
    labels = sort(unique(proj_complete$year))
  ) +
  ggplot2::labs(
    title = "Projected immune score across visits",
    x = "Year",
    y = "Projected immune score from baseline PCA",
    color = "Group"
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "AllVisits_projected_immune_score_by_group.png"),
  plot = p_proj,
  width = 7,
  height = 5,
  dpi = 600
)


# ---------------------------------------------------------
# Save PCA summary
# ---------------------------------------------------------

sink(file.path(out_dir, "PCA_summary.txt"))

cat("Baseline PCA summary\n")
cat("====================\n\n")
print(summary(pca_fit))

cat("\n\nBaseline PCA variance table\n")
cat("===========================\n")
print(pca_variance_table)

cat("\n\nOriented PC1 loadings\n")
cat("====================\n")
print(pca_loadings_table)

cat("\n\nBaseline correlations with NLR/MLR\n")
cat("==================================\n")
print(cor_table)

cat("\n\nVisitwise variance explained\n")
cat("===========================\n")
print(visit_variance_table)

cat("\n\nVisitwise PC1 loadings\n")
cat("=====================\n")
print(visit_loadings_table)

cat("\n\nVisitwise correlations of projected immune score with NLR/MLR\n")
cat("=============================================================\n")
print(visit_corrs)

sink()


# ---------------------------------------------------------
# Console output
# ---------------------------------------------------------

cat("\n============================================================\n")
cat("Baseline immune PCA and longitudinal projection completed.\n")
cat("Outputs saved to:\n")
cat(out_dir, "\n")
cat("============================================================\n")