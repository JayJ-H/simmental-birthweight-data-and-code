# Run from the extracted release root.
if (.Platform$OS.type == "windows") {
  Sys.setlocale("LC_ALL", "English_United States.utf8")
}
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
stopifnot(dir.exists("original_analysis"), dir.exists("revision_analysis"))
if ("--verify-only" %in% args) {
  d <- read.csv("original_analysis/data/processed/singleton_analysis_with_tail_phenotypes.csv",
                colClasses = c(dam_id = "character", calf_id = "character"))
  stopifnot(nrow(d) == 9384, length(unique(d$dam_id)) == 4624,
            sum(d$lower_tail) == 948, sum(d$upper_tail) == 955)
  cat("Included analytical cohort verified. No models were refitted.\n")
  quit(status = 0)
}
required <- c("readxl", "dplyr", "tidyr", "stringr", "lubridate", "forcats",
              "purrr", "tibble", "openxlsx", "ggplot2", "patchwork", "cowplot",
              "scales", "emmeans", "car", "lme4", "lmerTest", "broom",
              "broom.mixed", "quantreg", "sandwich", "lmtest", "mgcv",
              "scatterplot3d", "svglite", "ragg", "ggrepel", "multcompView")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing R packages: ", paste(missing, collapse = ", "))
run_stage <- function(directory, scripts) {
  setwd(file.path(root, directory))
  on.exit(setwd(root))
  for (script in scripts) {
    message("Running ", directory, "/", script)
    source(script, local = .GlobalEnv, encoding = "UTF-8")
  }
}
input <- "original_analysis/data/processed/singleton_analysis_with_tail_phenotypes.csv"
before <- read.csv(input, colClasses = c(dam_id = "character", calf_id = "character",
                    service_sire = "character", maternal_grandsire = "character"),
                   check.names = FALSE)
run_stage("original_analysis", c("scripts/01_data_cleaning.R", "scripts/06_tail_phenotypes.R"))
after <- read.csv(input, colClasses = c(dam_id = "character", calf_id = "character",
                   service_sire = "character", maternal_grandsire = "character"),
                  check.names = FALSE)
stopifnot(nrow(after) == 9384, length(unique(after$dam_id)) == 4624,
          sum(after$lower_tail) == 948, sum(after$upper_tail) == 955)
stopifnot(isTRUE(all.equal(before, after, check.attributes = FALSE, tolerance = 1e-10)))
cat("Raw-to-analysis reconstruction matched all included columns and rows.\n")
if ("--clean-only" %in% args) quit(status = 0)

run_stage("original_analysis", c(
  "scripts/02_qc_summary.R", "scripts/03_descriptive_stats.R",
  "scripts/05_lmm_birth_weight.R", "scripts/07_glmm_tail_models.R",
  "scripts/08_quantile_regression.R", "scripts/09_service_sire_mgs_exploratory.R",
  "scripts/10_sensitivity_analysis.R", "scripts/13_make_enhanced_figures.R",
  "scripts/14_make_3d_figures.R", "scripts/15_make_advanced_candidate_figures.R",
  "scripts/16_make_reference_style_triptych.R", "scripts/17_make_final_grouped_figures.R"))
run_stage("revision_analysis", c("scripts/01_revision_analysis.R",
                                "scripts/03_check_sparse_fit.R",
                                "scripts/01_revision_analysis.R"))
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("revision_analysis/scripts/14_restore_main_figures.R", shQuote(root)))
stopifnot(status == 0)
key_files <- c("annual_variability.csv", "thresholds.csv",
               "table3_lmm_means.csv", "table4_primary_OR.csv",
               "tail_model_effects.csv", "effect_comparisons.csv",
               "quantile_regression_cluster_bootstrap.csv")
for (name in key_files) {
  expected <- read.csv(file.path(root, "expected_results/revision", name), check.names = FALSE)
  observed <- read.csv(file.path(root, "revision_analysis/results", name), check.names = FALSE)
  stopifnot(identical(names(expected), names(observed)), nrow(expected) == nrow(observed))
  for (column in names(expected)) {
    if (is.numeric(expected[[column]])) {
      stopifnot(isTRUE(all.equal(expected[[column]], observed[[column]],
                                check.attributes = FALSE, tolerance = 0.005)))
    } else {
      stopifnot(identical(as.character(expected[[column]]), as.character(observed[[column]])))
    }
  }
}
cat("Final numerical estimates matched reference outputs within tolerance 0.005.\n")
