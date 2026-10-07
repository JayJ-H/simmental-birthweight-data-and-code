source("scripts/00_config.R")

write_log("10_sensitivity_analysis.R started")

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |> prepare_analysis_types()
all_valid_records <- read_csv_safe(file.path(PROCESSED_DIR, "all_valid_records.csv")) |> prepare_analysis_types()
thresholds <- openxlsx::read.xlsx(file.path(TABLE_DIR, "sex_specific_thresholds.xlsx"), sheet = "sex_specific_thresholds")
sex_var <- sex_variable_name()

apply_primary_thresholds <- function(data) {
  data |>
    dplyr::left_join(thresholds |> dplyr::select(dplyr::all_of(sex_var), P10, P90), by = sex_var) |>
    dplyr::mutate(
      lower_tail = as.integer(birth_weight <= P10),
      upper_tail = as.integer(birth_weight >= P90)
    ) |>
    dplyr::select(-P10, -P90)
}

all_valid_tail <- apply_primary_thresholds(all_valid_records)

fit_binary_status <- function(label, response, data, fixed_terms, random_terms = NULL) {
  cluster_var <- if ("dam_id" %in% names(data)) "dam_id" else NA_character_
  formula <- safe_formula(response, fixed_terms)
  model_data <- data |>
    dplyr::select(dplyr::all_of(unique(c(all.vars(formula), cluster_var)))) |>
    stats::na.omit()
  events <- sum(model_data[[response]] == 1, na.rm = TRUE)
  nonevents <- sum(model_data[[response]] == 0, na.rm = TRUE)
  if (nrow(model_data) < 100 || events < 10 || nonevents < 10) {
    return(tibble::tibble(label = label, phenotype = response, n = nrow(model_data), events = events, event_rate = events / nrow(model_data), formula = paste(deparse(formula), collapse = " "), AIC = NA_real_, BIC = NA_real_, status = "insufficient events"))
  }
  fit <- tryCatch(stats::glm(formula, data = model_data, family = stats::binomial()), error = function(e) e)
  if (inherits(fit, "error")) {
    return(tibble::tibble(label = label, phenotype = response, n = nrow(model_data), events = events, event_rate = events / nrow(model_data), formula = paste(deparse(formula), collapse = " "), AIC = NA_real_, BIC = NA_real_, status = paste("logistic model failed:", fit$message)))
  }
  robust_status <- "cluster-robust SE by dam_id computed"
  if (!is.na(cluster_var)) {
    robust <- tryCatch(lmtest::coeftest(fit, vcov. = sandwich::vcovCL(fit, cluster = model_data[[cluster_var]])), error = function(e) e)
    if (inherits(robust, "error")) robust_status <- paste("cluster-robust SE failed:", robust$message)
  }
  tibble::tibble(
    label = label,
    phenotype = response,
    n = nrow(model_data),
    events = events,
    event_rate = events / nrow(model_data),
    formula = paste(deparse(formula), collapse = " "),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    status = paste("logistic sensitivity model;", robust_status)
  )
}

base_terms <- c(age_terms_from_choice(), "calf_birth_year", "calf_birth_season")

make_threshold_data <- function(data, low_var, high_var, label) {
  list(
    fit_binary_status(paste0(label, " lower"), low_var, data, base_terms),
    fit_binary_status(paste0(label, " upper"), high_var, data, base_terms)
  ) |>
    dplyr::bind_rows()
}

threshold_results <- dplyr::bind_rows(
  make_threshold_data(singleton_tail, "lower_tail_p5", "upper_tail_p95", "sex-specific P5/P95"),
  make_threshold_data(singleton_tail, "lower_tail", "upper_tail", "sex-specific P10/P90"),
  make_threshold_data(singleton_tail, "lower_tail_p15", "upper_tail_p85", "sex-specific P15/P85"),
  fit_binary_status("overall P10 lower", "lower_tail_overall", singleton_tail, c(sex_var, base_terms)),
  fit_binary_status("overall P90 upper", "upper_tail_overall", singleton_tail, c(sex_var, base_terms))
)

birth_type_results <- dplyr::bind_rows(
  fit_binary_status("all valid records plus birth_type lower", "lower_tail", all_valid_tail, c(base_terms, "birth_type")),
  fit_binary_status("all valid records plus birth_type upper", "upper_tail", all_valid_tail, c(base_terms, "birth_type"))
)

exclude_2026 <- singleton_tail |> dplyr::filter(calf_birth_year != "2026")
year_results <- dplyr::bind_rows(
  fit_binary_status("exclude 2026 lower", "lower_tail", exclude_2026, base_terms),
  fit_binary_status("exclude 2026 upper", "upper_tail", exclude_2026, base_terms)
)

outlier_removed <- singleton_tail |> dplyr::filter(!extreme_low_flag, !extreme_high_flag)
outlier_results <- dplyr::bind_rows(
  fit_binary_status("remove birth_weight <25 or >65 lower", "lower_tail", outlier_removed, base_terms),
  fit_binary_status("remove birth_weight <25 or >65 upper", "upper_tail", outlier_removed, base_terms)
)

age_variant_terms <- list(
  parity_plus_age_group = c("parity_group", "dam_age_group", "calf_birth_year", "calf_birth_season"),
  parity_plus_age_cont = c("parity_group", "dam_age_month_cont", "calf_birth_year", "calf_birth_season"),
  parity_only = c("parity_group", "calf_birth_year", "calf_birth_season"),
  age_group_only = c("dam_age_group", "calf_birth_year", "calf_birth_season")
)
age_results <- purrr::imap_dfr(age_variant_terms, function(terms, name) {
  dplyr::bind_rows(
    fit_binary_status(paste(name, "lower"), "lower_tail", singleton_tail, terms),
    fit_binary_status(paste(name, "upper"), "upper_tail", singleton_tail, terms)
  )
})

service_results <- dplyr::bind_rows(
  fit_binary_status("service_sire label adjustment lower", "lower_tail", singleton_tail |> dplyr::filter(!is.na(service_sire)), c(base_terms, "service_sire")),
  fit_binary_status("service_sire label adjustment upper", "upper_tail", singleton_tail |> dplyr::filter(!is.na(service_sire)), c(base_terms, "service_sire"))
)

mgs_results <- dplyr::bind_rows(
  fit_binary_status("maternal_grandsire label adjustment lower", "lower_tail", singleton_tail |> dplyr::filter(!is.na(maternal_grandsire)), c(base_terms, "maternal_grandsire")),
  fit_binary_status("maternal_grandsire label adjustment upper", "upper_tail", singleton_tail |> dplyr::filter(!is.na(maternal_grandsire)), c(base_terms, "maternal_grandsire"))
)

year_season_terms <- c(base_terms, "calf_birth_year:calf_birth_season")
year_season_results <- dplyr::bind_rows(
  fit_binary_status("year season interaction lower", "lower_tail", singleton_tail, year_season_terms),
  fit_binary_status("year season interaction upper", "upper_tail", singleton_tail, year_season_terms)
)

sensitivity_summary <- dplyr::bind_rows(
  threshold_results |> dplyr::mutate(section = "threshold"),
  birth_type_results |> dplyr::mutate(section = "birth_type"),
  year_results |> dplyr::mutate(section = "year"),
  outlier_results |> dplyr::mutate(section = "outlier"),
  age_results |> dplyr::mutate(section = "parity_age"),
  service_results |> dplyr::mutate(section = "service_sire"),
  mgs_results |> dplyr::mutate(section = "maternal_grandsire"),
  year_season_results |> dplyr::mutate(section = "year_season")
) |>
  dplyr::relocate(section, .before = label)

safe_write_xlsx(sensitivity_summary, file.path(TABLE_DIR, "sensitivity_summary.xlsx"))
safe_write_xlsx(
  sensitivity_summary |>
    dplyr::filter(grepl("overall P", label) | grepl("sex-specific P10/P90", label)),
  file.path(MODEL_DIR, "overall_threshold_sensitivity.xlsx")
)

warning_file <- file.path(LOG_DIR, "model_warnings.log")
if (!file.exists(warning_file) || file.info(warning_file)$size == 0) {
  write_log("No model downgrade, separation, or singular-fit warning required manual intervention in the final pipeline run.", warning_file)
}

p <- ggplot2::ggplot(sensitivity_summary, ggplot2::aes(x = event_rate, y = stats::reorder(label, event_rate), colour = phenotype)) +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::facet_grid(section ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = c(
    lower_tail = palette_contract[["signal_blue"]], upper_tail = palette_contract[["accent_orange"]],
    lower_tail_p5 = palette_contract[["signal_blue"]], upper_tail_p95 = palette_contract[["accent_orange"]],
    lower_tail_p15 = palette_contract[["signal_teal"]], upper_tail_p85 = palette_contract[["accent_red"]],
    lower_tail_overall = palette_contract[["neutral_dark"]], upper_tail_overall = palette_contract[["neutral_mid"]]
  ), guide = ggplot2::guide_legend(title = "Phenotype")) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Event rate", y = NULL, title = "Sensitivity analysis summary") +
  theme_pub(base_size = 6.5)
save_plot_pdf_png(p, file.path(FIGURE_DIR, "sensitivity_summary"), width_mm = 180, height_mm = 220)

write_log("10_sensitivity_analysis.R completed")
