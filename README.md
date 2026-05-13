# PPMI Peripheral Immune Score, DAT Binding, and Motor Progression Analysis

This repository contains the R code used for the manuscript examining the association between peripheral immune activation, striatal DAT binding, and longitudinal motor progression in Parkinson's disease using data from the Parkinson's Progression Markers Initiative (PPMI).

## Data availability

The individual-level data used in this study are not included in this repository.

Data used in the preparation of this article were obtained on 2025-12-11 from the Parkinson's Progression Markers Initiative (PPMI) database, RRID:SCR_006431.

Researchers can request access to the PPMI dataset through:

www.ppmi-info.org/access-data-specimens/download-data

For up-to-date information on the study, visit:

www.ppmi-info.org

## Analysis scripts

### `R/01_baseline_demographic_characteristics.R`

Generates baseline demographic, clinical, and laboratory comparison tables for Parkinson's disease participants and healthy controls.

### `R/02_baseline_immune_pca_projection.R`

Performs baseline PCA using peripheral immune markers and derives the Peripheral Immune Score (PIS). The baseline PCA solution is also projected onto longitudinal visits.

### `R/03_longitudinal_pis_mixed_model.R`

Tests longitudinal group differences in PIS between Parkinson's disease participants and healthy controls.

### `R/04_pis_motor_progression_mixed_model.R`

Tests whether baseline PIS predicts longitudinal motor progression in Parkinson's disease using UPDRS-III as the outcome.

### `R/05_pis_dat_decline_mixed_models.R`

Tests whether baseline PIS predicts longitudinal DAT-SPECT decline in putamen and caudate regions.

### `R/06_striatal_dat_moderation_models.R`

Tests whether baseline striatal DAT binding moderates the association between baseline PIS and longitudinal motor progression.

### `R/07_striatal_dat_moderation_sensitivity_4y.R`

Repeats the striatal DAT moderation analyses as a sensitivity analysis restricted to the first 4 years of follow-up.

## How to run

Place the curated local PPMI analysis dataset at the file path specified in each R script, or update the `data_path` variable in the scripts.

Then run:

```r
source("run_all_analyses.R")