# PPMI NLR, DAT Binding, and Motor Progression Analyses

This repository contains the R code used for analyses examining neutrophil-to-lymphocyte ratio (NLR), striatal dopamine transporter (DAT) binding, and longitudinal motor progression in Parkinson's disease using data from the Parkinson's Progression Markers Initiative (PPMI).

The workflow includes:

- calculation and preparation of visit-wise and baseline NLR variables;
- baseline comparisons between Parkinson's disease participants and healthy controls;
- longitudinal NLR analyses;
- baseline NLR models of motor progression and DAT decline;
- between-person and within-person time-varying NLR analyses;
- moderation analyses testing whether baseline striatal DAT binding modifies the association between baseline NLR and motor progression;
- a sensitivity analysis restricted to the DAT imaging interval.

## Data availability

Individual-level PPMI data are not included in this repository.

Data used in the preparation of this study were obtained on **2025-12-11** from the Parkinson's Progression Markers Initiative database (**PPMI; RRID:SCR_006431**).

Researchers may request access to PPMI data through:

<https://www.ppmi-info.org/access-data-specimens/download-data>

For current information about the PPMI study, visit:

<https://www.ppmi-info.org>

Users are responsible for complying with the PPMI data-use agreement. Restricted participant-level data and generated participant-level outputs must not be committed to a public repository.

## Repository structure

```text
.
├── R/
│   ├── 01_calculate_nlr_and_prepare_datasets.R
│   ├── 02_baseline_demographic_characteristics.R
│   ├── 03_longitudinal_nlr_analysis.R
│   ├── 04_baseline_nlr_motor_progression.R
│   ├── 05_timevarying_nlr_motor_severity.R
│   ├── 06_baseline_nlr_dat_decline.R
│   ├── 07_timevarying_nlr_dat_binding.R
│   ├── 08_baseline_dat_moderation_of_nlr_motor_progression.R
│   └── 09_baseline_dat_moderation_4year_sensitivity.R
├── run_all_analyses.R
├── README.md
└── LICENSE
```

No data or output directories are included in the repository.

All input datasets and generated outputs are stored outside the repository. Paths are supplied through command-line arguments, environment variables, or `run_all_analyses.R`.

## Analysis workflow

The scripts are intended to be run in numerical order.

### `R/01_calculate_nlr_and_prepare_datasets.R`

Calculates visit-wise and baseline inflammatory ratios from the merged PPMI blood dataset.

Derived variables include:

- neutrophil-to-lymphocyte ratio (NLR);
- monocyte-to-lymphocyte ratio (MLR);
- platelet-to-lymphocyte ratio (PLR);
- baseline blood-cell and inflammatory-ratio variables;
- standardized visit-wise and baseline NLR variables;
- log-transformed NLR variables.

The script requires:

1. the raw merged PPMI Excel workbook;
2. a common external output root.

Script 01 automatically creates its output directory and writes the main derived workbook to:

```text
<PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
```

The derived workbook is then used by Scripts 02-09.

### `R/02_baseline_demographic_characteristics.R`

Generates baseline demographic, clinical, imaging, and inflammatory comparisons between Parkinson's disease participants and healthy controls.

The script also produces:

- Welch t-tests and Wilcoxon sensitivity analyses;
- normality checks;
- chi-square assumption checks;
- visit-level follow-up and attrition tables;
- follow-up flow and attrition plots;
- manuscript-ready sample-size and baseline-summary text.

### `R/03_longitudinal_nlr_analysis.R`

Compares longitudinal NLR trajectories between Parkinson's disease participants and healthy controls using a linear mixed-effects model.

Primary model:

```r
NLR ~ GROUP * year_c + age + sex + bmi +
  (1 | SITE) + (1 | PATNO)
```

### `R/04_baseline_nlr_motor_progression.R`

Tests whether baseline NLR predicts longitudinal motor progression in participants with Parkinson's disease.

Outcome:

```text
ON-medication MDS-UPDRS Part III / UPDRS-III
```

Primary model:

```r
updrs3_score_on ~ baseline_NLR_z * year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

The script also includes:

- simple-slope analyses;
- model-estimated trajectories at low, mean, and high baseline NLR;
- tertile-based descriptive and model-estimated trajectories;
- a follow-up-only model adjusted for baseline UPDRS-III;
- a log-baseline-NLR sensitivity model.

### `R/05_timevarying_nlr_motor_severity.R`

Tests whether time-varying NLR tracks longitudinal motor severity in Parkinson's disease.

Visit-wise NLR is decomposed into:

- `NLR_between_z`: between-person differences in mean NLR;
- `NLR_within_z`: within-person deviations from each participant's own mean NLR.

Primary model:

```r
updrs3_score_on ~ NLR_between_z + NLR_within_z + year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

Secondary interaction model:

```r
updrs3_score_on ~ NLR_between_z + NLR_within_z * year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

### `R/06_baseline_nlr_dat_decline.R`

Tests whether baseline NLR predicts longitudinal DAT decline in the bilateral putamen and caudate.

DAT imaging visits:

```text
BL, V04, V06, and V10
```

Regional outcomes:

- right putamen;
- left putamen;
- right caudate;
- left caudate.

Primary model for each regional outcome:

```r
DAT_outcome ~ baseline_NLR_z * year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

The script also includes baseline-DAT-adjusted follow-up models and log-baseline-NLR sensitivity models.

### `R/07_timevarying_nlr_dat_binding.R`

Tests whether between-person and within-person variation in NLR tracks longitudinal DAT binding in the bilateral putamen and caudate.

Primary model for each regional outcome:

```r
DAT_outcome ~ NLR_between_z + NLR_within_z + year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

Secondary interaction model:

```r
DAT_outcome ~ NLR_between_z + NLR_within_z * year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

### `R/08_baseline_dat_moderation_of_nlr_motor_progression.R`

Tests whether baseline striatal DAT binding moderates the association between baseline NLR and longitudinal motor progression.

Moderators:

- left caudate DAT;
- right caudate DAT;
- left putamen DAT;
- right putamen DAT.

Primary moderation model:

```r
UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

The script includes simple slopes, model-estimated trajectories, raw tertile-based descriptive plots, and a follow-up-only sensitivity model adjusted for baseline UPDRS-III.

### `R/09_baseline_dat_moderation_4year_sensitivity.R`

Repeats the baseline DAT moderation analyses as a sensitivity analysis restricted to the DAT imaging interval.

Included visits:

```text
BL, V04, V06, and V10
```

Corresponding follow-up years:

```text
0, 1, 2, and 4
```

Primary model:

```r
UPDRS ~ baseline_NLR_z * DAT_BL_z * year_c +
  age + sex + bmi + duration_yrs + LEDD + DOMSIDE +
  (1 | SITE) + (1 | PATNO)
```

## Software requirements

The analyses were written in R and use packages available from CRAN.

Across the full workflow, the required packages include:

```r
c(
  "readxl",
  "openxlsx",
  "dplyr",
  "tidyr",
  "tibble",
  "broom",
  "ggplot2",
  "gt",
  "lme4",
  "lmerTest",
  "car",
  "emmeans",
  "performance"
)
```

Install the required packages before running the workflow:

```r
install.packages(
  c(
    "readxl",
    "openxlsx",
    "dplyr",
    "tidyr",
      "tibble",
    "broom",
    "ggplot2",
    "gt",
    "lme4",
    "lmerTest",
    "car",
    "emmeans",
    "performance"
  )
)
```

Each script checks for required packages and stops with an informative message when a package is missing. The scripts do not install packages automatically.

## Input and output locations

The repository does not contain data or generated outputs.

The raw PPMI workbook and all generated analysis outputs must be stored outside the repository.

All scripts use a common external output root (`PPMI_OUTPUT_ROOT`). Script 01 creates the derived NLR workbook automatically at:

```text
<PPMI_OUTPUT_ROOT>/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
```

Script 01 therefore requires two paths:

```text
1. Raw merged PPMI Excel workbook
2. Common output root
```

Scripts 02-09 also require two paths:

```text
1. Derived NLR workbook created by Script 01
2. The same common output root
```

Each script creates its own subdirectory under the common output root. The current workflow uses:

```text
01_NLR/
02_BASELINE_CHARACTERISTICS/
03_NLR_LONGITUDINAL/
04_NLR_MOTOR_PROGRESSION/
05_TIMEVARYING_NLR_UPDRS/
06_NLR_DAT_DECLINE/
07_TIMEVARYING_NLR_DAT/
08_NLR_MODERATION_STRIATAL_DAT/
09_NLR_MODERATION_STRIATAL_DAT_4Y/
workflow_logs/
```

The derived workbook created by Script 01 must contain the worksheet:

```text
All_data_with_NLR
```

Scripts perform column checks and attempt to identify common alternative names for variables such as age, sex, BMI, study site, disease duration, LEDD, dominant side, UPDRS-III, and regional DAT measures.

Paths may be supplied through command-line arguments. The analysis scripts also support the `PPMI_OUTPUT_ROOT` environment variable; Scripts 02-09 may additionally use `PPMI_NLR_DATA_FILE` to specify the derived NLR workbook.

## How to run

Run commands from the repository root.

Because data and output directories are external to the repository, explicit file paths should be supplied or configured through the workflow wrapper.

### Run Script 01

```bash
Rscript R/01_calculate_nlr_and_prepare_datasets.R \
  /path/to/PPMI_with_serum_all_visits_all_blood_merged.xlsx \
  /path/to/analysis_outputs
```

Script 01 will create:

```text
/path/to/analysis_outputs/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx
```

### Run Scripts 02-09 individually

Example:

```bash
Rscript R/04_baseline_nlr_motor_progression.R \
  /path/to/analysis_outputs/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx \
  /path/to/analysis_outputs
```

The same structure applies to Scripts 02-09:

```bash
Rscript R/<script_name>.R \
  /path/to/analysis_outputs/01_NLR/PPMI_with_NLR_all_visits_updated.xlsx \
  /path/to/analysis_outputs
```

Each script creates its own output subdirectory under the common output root.

### Run the complete workflow

The complete workflow can be launched through:

```bash
Rscript run_all_analyses.R
```

Before running it, define the external raw-data path and the main output directory in `run_all_analyses.R`.

Alternatively, provide both paths on the command line:

```bash
Rscript run_all_analyses.R \
  /path/to/PPMI_with_serum_all_visits_all_blood_merged.xlsx \
  /path/to/analysis_outputs
```

The wrapper derives the Script 01 output workbook path automatically and executes the scripts in this order:

```text
01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09
```

If any script fails, the workflow stops and the corresponding log file is retained under:

```text
<PPMI_OUTPUT_ROOT>/workflow_logs/
```

## Outputs

Each analysis writes its results to a script-specific subdirectory created automatically under the common external output root.

Outputs include, where applicable:

- complete-case and pre-filtering datasets;
- missingness summaries;
- descriptive statistics;
- fixed-effect tables;
- Type III model tests;
- random-effect variance estimates;
- model-fit indices;
- estimated marginal means and/or model-derived predictions;
- simple-slope analyses;
- sensitivity-model results;
- model-diagnostic plots;
- publication-quality PNG and PDF figures;
- text summaries;
- `sessionInfo.txt` files documenting the R environment.

Because some outputs may contain participant-level information, they should remain outside the public repository and should be handled according to the PPMI data-use agreement.

## Statistical framework

Longitudinal analyses use linear mixed-effects models with:

- participant-level random intercepts;
- study-site random intercepts;
- maximum-likelihood estimation;
- Type III tests;
- covariate adjustment specified within each script;
- model diagnostics and convergence checks;
- model-specific exclusion of observations missing variables required for that analysis, while retaining available longitudinal observations from participants with partial follow-up.

Time is coded as years from baseline, with baseline equal to zero.

Motor outcomes refer to ON-medication UPDRS-III values as represented in the curated analysis dataset.

Time-varying NLR analyses separate:

- between-person differences in mean NLR;
- within-person visit-level deviations from each participant's mean NLR.

## Reproducibility

Each script supports external input and output paths through command-line arguments; the workflow also supports a common `PPMI_OUTPUT_ROOT`, and Scripts 02-09 may use `PPMI_NLR_DATA_FILE` for the derived NLR workbook.

Each output directory contains a `sessionInfo.txt` file recording the R version, platform, and loaded package versions used during execution.

For stronger package-version reproducibility, the repository may additionally be managed with `renv`.

## Citation

When using this code, cite the associated manuscript and acknowledge the Parkinson's Progression Markers Initiative according to the current PPMI data-use and publication requirements.

## License

See the `LICENSE` file for the terms governing use and redistribution of the code in this repository.

## Disclaimer

This repository provides analysis code only. It does not redistribute PPMI data and does not replace the official PPMI documentation, data dictionary, or data-use agreement.