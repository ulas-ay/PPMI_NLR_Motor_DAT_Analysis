# ==============================================================================
# run_all_analyses.R
#
# Purpose
# -------
# Run the complete PPMI NLR, DAT binding, and motor progression workflow in
# numerical order.
#
# Repository contents expected
# ----------------------------
# R/
#   01_calculate_nlr_and_prepare_datasets.R
#   02_baseline_demographic_characteristics.R
#   03_longitudinal_nlr_analysis.R
#   04_baseline_nlr_motor_progression.R
#   05_timevarying_nlr_motor_severity.R
#   06_baseline_nlr_dat_decline.R
#   07_timevarying_nlr_dat_binding.R
#   08_baseline_dat_moderation_of_nlr_motor_progression.R
#   09_baseline_dat_moderation_4year_sensitivity.R
#
# Data and outputs are intentionally stored outside the repository.
#
# Usage
# -----
# 1. Edit the three paths in the USER CONFIGURATION section below.
# 2. From the repository root, run:
#
#      Rscript run_all_analyses.R
#
# Optional command-line usage
# ---------------------------
# Paths can also be supplied without editing this file:
#
# Rscript run_all_analyses.R \
#   /path/to/raw_input.xlsx \
#   /path/to/derived_output.xlsx \
#   /path/to/main_output_directory
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# USER CONFIGURATION
# ------------------------------------------------------------------------------

# Raw merged PPMI workbook used by Script 01.
raw_input_file <- "/path/to/PPMI_with_serum_all_visits_all_blood_merged.xlsx"

# Derived workbook created by Script 01 and used by Scripts 02-09.
derived_output_file <- "/path/to/PPMI_with_NLR_all_visits_updated.xlsx"

# Parent directory in which script-specific output folders will be created.
main_output_directory <- "/path/to/analysis_outputs"

# ------------------------------------------------------------------------------
# Optional command-line path overrides
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  raw_input_file <- args[[1]]
}

if (length(args) >= 2) {
  derived_output_file <- args[[2]]
}

if (length(args) >= 3) {
  main_output_directory <- args[[3]]
}

if (length(args) > 3) {
  warning(
    "Only the first three command-line arguments are used.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

normalize_existing_file <- function(path, label) {
  if (
    is.na(path) ||
    !nzchar(path) ||
    grepl("^/path/to/", path)
  ) {
    stop(
      paste0(
        label,
        " has not been configured.\n",
        "Edit run_all_analyses.R or provide the path as a command-line argument."
      ),
      call. = FALSE
    )
  }

  if (!file.exists(path)) {
    stop(
      paste0(
        label,
        " was not found:\n",
        normalizePath(path, winslash = "/", mustWork = FALSE)
      ),
      call. = FALSE
    )
  }

  normalizePath(path, winslash = "/", mustWork = TRUE)
}

normalize_output_path <- function(path, label) {
  if (
    is.na(path) ||
    !nzchar(path) ||
    grepl("^/path/to/", path)
  ) {
    stop(
      paste0(
        label,
        " has not been configured.\n",
        "Edit run_all_analyses.R or provide the path as a command-line argument."
      ),
      call. = FALSE
    )
  }

  normalizePath(path, winslash = "/", mustWork = FALSE)
}

get_script_directory <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_args, value = TRUE)

  if (length(file_argument) > 0) {
    script_path <- sub("^--file=", "", file_argument[[1]])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

run_r_script <- function(script_path, script_args, log_file) {
  command <- file.path(R.home("bin"), "Rscript")

  message("\n", strrep("=", 78))
  message("Running: ", basename(script_path))
  message(strrep("=", 78))

  start_time <- Sys.time()

  output <- system2(
    command = command,
    args = c(
      "--vanilla",
      shQuote(script_path),
      vapply(script_args, shQuote, character(1))
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  exit_status <- attr(output, "status")

  if (is.null(exit_status)) {
    exit_status <- 0L
  }

  writeLines(
    c(
      paste0("Script: ", script_path),
      paste0("Started: ", format(start_time, "%Y-%m-%d %H:%M:%S")),
      paste0("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("Exit status: ", exit_status),
      "",
      output
    ),
    con = log_file
  )

  if (length(output) > 0) {
    cat(paste(output, collapse = "\n"), "\n")
  }

  if (exit_status != 0L) {
    stop(
      paste0(
        "Analysis failed in ",
        basename(script_path),
        ".\nSee the log file:\n",
        normalizePath(log_file, winslash = "/", mustWork = FALSE)
      ),
      call. = FALSE
    )
  }

  message(
    "Completed: ",
    basename(script_path),
    " (",
    round(
      as.numeric(
        difftime(Sys.time(), start_time, units = "mins")
      ),
      2
    ),
    " minutes)"
  )

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Resolve repository paths
# ------------------------------------------------------------------------------

repository_root <- get_script_directory()
r_directory <- file.path(repository_root, "R")

if (!dir.exists(r_directory)) {
  stop(
    paste0(
      "The R directory was not found:\n",
      normalizePath(r_directory, winslash = "/", mustWork = FALSE)
    ),
    call. = FALSE
  )
}

analysis_scripts <- c(
  "01_calculate_nlr_and_prepare_datasets.R",
  "02_baseline_demographic_characteristics.R",
  "03_longitudinal_nlr_analysis.R",
  "04_baseline_nlr_motor_progression.R",
  "05_timevarying_nlr_motor_severity.R",
  "06_baseline_nlr_dat_decline.R",
  "07_timevarying_nlr_dat_binding.R",
  "08_baseline_dat_moderation_of_nlr_motor_progression.R",
  "09_baseline_dat_moderation_4year_sensitivity.R"
)

script_paths <- file.path(r_directory, analysis_scripts)
missing_scripts <- analysis_scripts[!file.exists(script_paths)]

if (length(missing_scripts) > 0) {
  stop(
    paste0(
      "The following analysis scripts were not found in the R directory:\n",
      paste(missing_scripts, collapse = "\n")
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Validate external paths
# ------------------------------------------------------------------------------

raw_input_file <- normalize_existing_file(
  raw_input_file,
  "Raw input workbook"
)

derived_output_file <- normalize_output_path(
  derived_output_file,
  "Derived output workbook"
)

main_output_directory <- normalize_output_path(
  main_output_directory,
  "Main output directory"
)

dir.create(
  dirname(derived_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  main_output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(dirname(derived_output_file))) {
  stop(
    paste0(
      "Could not create the derived-data directory:\n",
      dirname(derived_output_file)
    ),
    call. = FALSE
  )
}

if (!dir.exists(main_output_directory)) {
  stop(
    paste0(
      "Could not create the main output directory:\n",
      main_output_directory
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Script-specific output directories
# ------------------------------------------------------------------------------

output_directories <- c(
  file.path(main_output_directory, "01_nlr_calculation"),
  file.path(
    main_output_directory,
    "02_baseline_demographic_characteristics"
  ),
  file.path(main_output_directory, "03_longitudinal_nlr_analysis"),
  file.path(main_output_directory, "04_baseline_nlr_motor_progression"),
  file.path(main_output_directory, "05_timevarying_nlr_motor_severity"),
  file.path(main_output_directory, "06_baseline_nlr_dat_decline"),
  file.path(main_output_directory, "07_timevarying_nlr_dat_binding"),
  file.path(main_output_directory, "08_baseline_dat_moderation"),
  file.path(
    main_output_directory,
    "09_baseline_dat_moderation_4year_sensitivity"
  )
)

invisible(
  lapply(
    output_directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

# ------------------------------------------------------------------------------
# Analysis plan
# ------------------------------------------------------------------------------

analysis_plan <- list(
  list(
    script = script_paths[[1]],
    args = c(
      raw_input_file,
      derived_output_file,
      output_directories[[1]]
    )
  ),
  list(
    script = script_paths[[2]],
    args = c(
      derived_output_file,
      output_directories[[2]]
    )
  ),
  list(
    script = script_paths[[3]],
    args = c(
      derived_output_file,
      output_directories[[3]]
    )
  ),
  list(
    script = script_paths[[4]],
    args = c(
      derived_output_file,
      output_directories[[4]]
    )
  ),
  list(
    script = script_paths[[5]],
    args = c(
      derived_output_file,
      output_directories[[5]]
    )
  ),
  list(
    script = script_paths[[6]],
    args = c(
      derived_output_file,
      output_directories[[6]]
    )
  ),
  list(
    script = script_paths[[7]],
    args = c(
      derived_output_file,
      output_directories[[7]]
    )
  ),
  list(
    script = script_paths[[8]],
    args = c(
      derived_output_file,
      output_directories[[8]]
    )
  ),
  list(
    script = script_paths[[9]],
    args = c(
      derived_output_file,
      output_directories[[9]]
    )
  )
)

# ------------------------------------------------------------------------------
# Save workflow configuration
# ------------------------------------------------------------------------------

workflow_log_directory <- file.path(
  main_output_directory,
  "workflow_logs"
)

dir.create(
  workflow_log_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

configuration <- data.frame(
  item = c(
    "repository_root",
    "raw_input_file",
    "derived_output_file",
    "main_output_directory",
    "R_version",
    "platform",
    "workflow_started"
  ),
  value = c(
    repository_root,
    raw_input_file,
    derived_output_file,
    main_output_directory,
    R.version.string,
    R.version$platform,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  configuration,
  file.path(workflow_log_directory, "workflow_configuration.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Run analyses
# ------------------------------------------------------------------------------

workflow_start <- Sys.time()

message("\nPPMI NLR analysis workflow")
message("Repository root: ", repository_root)
message("Raw input: ", raw_input_file)
message("Derived workbook: ", derived_output_file)
message("Output root: ", main_output_directory)
message("Scripts to run: ", length(analysis_plan))

completed_scripts <- character(0)

for (index in seq_along(analysis_plan)) {
  plan_item <- analysis_plan[[index]]
  script_name <- basename(plan_item$script)

  log_file <- file.path(
    workflow_log_directory,
    paste0(
      sprintf("%02d", index),
      "_",
      tools::file_path_sans_ext(script_name),
      ".log"
    )
  )

  run_r_script(
    script_path = plan_item$script,
    script_args = plan_item$args,
    log_file = log_file
  )

  completed_scripts <- c(completed_scripts, script_name)

  progress_table <- data.frame(
    order = seq_along(completed_scripts),
    script = completed_scripts,
    status = "completed",
    stringsAsFactors = FALSE
  )

  write.csv(
    progress_table,
    file.path(workflow_log_directory, "workflow_progress.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# Final workflow summary
# ------------------------------------------------------------------------------

workflow_end <- Sys.time()
elapsed_minutes <- as.numeric(
  difftime(workflow_end, workflow_start, units = "mins")
)

workflow_summary <- data.frame(
  workflow_status = "completed",
  scripts_completed = length(completed_scripts),
  workflow_started = format(
    workflow_start,
    "%Y-%m-%d %H:%M:%S"
  ),
  workflow_finished = format(
    workflow_end,
    "%Y-%m-%d %H:%M:%S"
  ),
  elapsed_minutes = elapsed_minutes,
  derived_workbook = derived_output_file,
  output_root = main_output_directory,
  stringsAsFactors = FALSE
)

write.csv(
  workflow_summary,
  file.path(workflow_log_directory, "workflow_summary.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(workflow_log_directory, "workflow_sessionInfo.txt")
)

message("\n", strrep("=", 78))
message("All nine analyses completed successfully.")
message(
  "Total elapsed time: ",
  round(elapsed_minutes, 2),
  " minutes"
)
message("Derived workbook: ", derived_output_file)
message("Output root: ", main_output_directory)
message("Workflow logs: ", workflow_log_directory)
message(strrep("=", 78))
