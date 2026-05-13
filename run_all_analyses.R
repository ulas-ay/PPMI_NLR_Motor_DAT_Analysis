# =========================================================
# Run all manuscript analyses
# =========================================================

source("R/01_baseline_demographic_characteristics.R")
source("R/02_baseline_immune_pca_projection.R")
source("R/03_longitudinal_pis_mixed_model.R")
source("R/04_pis_motor_progression_mixed_model.R")
source("R/05_pis_dat_decline_mixed_models.R")
source("R/06_striatal_dat_moderation_models.R")
source("R/07_striatal_dat_moderation_sensitivity_4y.R")

cat("\nAll analyses completed.\n")