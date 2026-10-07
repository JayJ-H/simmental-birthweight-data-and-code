source("scripts/00_config.R")

write_log("09_service_sire_mgs_exploratory.R started")

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |> prepare_analysis_types()
sex_var <- sex_variable_name()

explore_label <- function(data, label_var, min_n, display_name) {
  data_sub <- data |>
    dplyr::filter(!is.na(.data[[label_var]])) |>
    dplyr::add_count(.data[[label_var]], name = "label_n") |>
    dplyr::filter(label_n >= min_n) |>
    droplevels()

  raw_summary <- data_sub |>
    dplyr::group_by(.data[[label_var]]) |>
    dplyr::summarise(
      n = dplyr::n(),
      raw_mean_birth_weight = mean(birth_weight, na.rm = TRUE),
      SD = stats::sd(birth_weight, na.rm = TRUE),
      SE = SD / sqrt(n),
      lower_tail_rate = mean(lower_tail, na.rm = TRUE),
      upper_tail_rate = mean(upper_tail, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n), .data[[label_var]])
  names(raw_summary)[1] <- label_var

  if (nrow(data_sub) == 0 || dplyr::n_distinct(data_sub[[label_var]]) < 2) {
    return(list(raw_summary = raw_summary, adjusted = tibble::tibble(), conditional = tibble::tibble(), data = data_sub))
  }

  fixed_terms <- c(sex_var, age_terms_from_choice(), "calf_birth_year", "calf_birth_season", label_var)
  formula <- safe_formula("birth_weight", fixed_terms)
  model_data <- data_sub |>
    dplyr::select(dplyr::all_of(all.vars(formula))) |>
    stats::na.omit()
  lm_fit <- stats::lm(formula, data = model_data)
  emm <- tryCatch(emmeans::emmeans(lm_fit, stats::as.formula(paste("~", label_var))), error = function(e) e)
  adjusted <- if (inherits(emm, "error")) {
    tibble::tibble(status = emm$message)
  } else {
    as.data.frame(emm) |>
      dplyr::arrange(emmean)
  }

  random_formula <- safe_formula("birth_weight", c(sex_var, age_terms_from_choice(), "calf_birth_year", "calf_birth_season"), random_terms = paste0("(1 | ", label_var, ")"))
  random_fit <- tryCatch(lme4::lmer(random_formula, data = model_data, REML = TRUE), error = function(e) e)
  conditional <- if (inherits(random_fit, "error")) {
    tibble::tibble(status = random_fit$message)
  } else {
    ran <- as.data.frame(lme4::ranef(random_fit)[[label_var]]) |>
      tibble::rownames_to_column(label_var)
    names(ran)[2] <- "conditional_random_intercept"
    ran |>
      dplyr::left_join(raw_summary |> dplyr::select(dplyr::all_of(label_var), n), by = label_var) |>
      dplyr::arrange(conditional_random_intercept)
  }

  list(raw_summary = raw_summary, adjusted = adjusted, conditional = conditional, data = data_sub)
}

service_30 <- explore_label(singleton_tail, "service_sire", 30, "service sire")
service_50 <- explore_label(singleton_tail, "service_sire", 50, "service sire")
mgs_30 <- explore_label(singleton_tail, "maternal_grandsire", MIN_MGS_N, "maternal grandsire")

safe_write_xlsx(
  list(
    n_ge_30_raw_summary = service_30$raw_summary,
    n_ge_30_adjusted_estimates = service_30$adjusted,
    n_ge_30_conditional_predictions = service_30$conditional,
    n_ge_50_raw_summary = service_50$raw_summary,
    n_ge_50_adjusted_estimates = service_50$adjusted,
    n_ge_50_conditional_predictions = service_50$conditional,
    notes = tibble::tibble(note = "Service-sire results are exploratory phenotypic label-associated differences only.")
  ),
  file.path(TABLE_DIR, "service_sire_exploratory.xlsx")
)

safe_write_xlsx(
  list(
    n_ge_30_raw_summary = mgs_30$raw_summary,
    n_ge_30_adjusted_estimates = mgs_30$adjusted,
    n_ge_30_conditional_predictions = mgs_30$conditional,
    notes = tibble::tibble(note = "Maternal-grandsire-label results are exploratory phenotypic label-associated differences only.")
  ),
  file.path(TABLE_DIR, "maternal_grandsire_exploratory.xlsx")
)

plot_exploratory <- function(summary_table, label_var, title, filename) {
  if (nrow(summary_table) == 0) return(invisible(NULL))
  plot_data <- summary_table |>
    dplyr::arrange(raw_mean_birth_weight) |>
    dplyr::mutate(label_plot = factor(.data[[label_var]], levels = .data[[label_var]]))
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = raw_mean_birth_weight, y = label_plot)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = raw_mean_birth_weight - 1.96 * SE, xmax = raw_mean_birth_weight + 1.96 * SE), height = 0.18, linewidth = 0.3, colour = palette_contract[["neutral_mid"]]) +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = palette_contract[["signal_blue"]], alpha = 0.85) +
    ggplot2::scale_size_continuous(name = "Records", range = c(1, 3.5)) +
    ggplot2::labs(x = "Raw mean birth weight (kg)", y = NULL, title = title) +
    theme_pub(base_size = 6.5)
  save_plot_pdf_png(p, file.path(FIGURE_DIR, filename), width_mm = 160, height_mm = 180)
}

plot_exploratory(service_30$raw_summary, "service_sire", "Service-sire-associated phenotypic differences", "service_sire_exploratory")
plot_exploratory(mgs_30$raw_summary, "maternal_grandsire", "Maternal-grandsire-label-associated phenotypic differences", "maternal_grandsire_exploratory")

write_log("09_service_sire_mgs_exploratory.R completed")
