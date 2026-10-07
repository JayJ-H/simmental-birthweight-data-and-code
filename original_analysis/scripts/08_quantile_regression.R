source("scripts/00_config.R")

write_log("08_quantile_regression.R started")

singleton_analysis <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis.csv")) |> prepare_analysis_types()
sex_var <- sex_variable_name()

fixed_terms <- c(sex_var, age_terms_from_choice(), "calf_birth_year", "calf_birth_season")
formula <- safe_formula("birth_weight", fixed_terms)
model_data <- singleton_analysis |>
  dplyr::select(dplyr::all_of(all.vars(formula))) |>
  stats::na.omit()

taus <- c(0.10, 0.50, 0.90)

fit_one_tau <- function(tau) {
  fit <- quantreg::rq(formula, tau = tau, data = model_data)
  s <- summary(fit, se = "nid")
  coefs <- as.data.frame(s$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  names(coefs) <- c("estimate", "std_error", "t_value", "p_value", "term")
  coefs |>
    dplyr::mutate(
      tau = tau,
      CI_low = estimate - 1.96 * std_error,
      CI_high = estimate + 1.96 * std_error,
      .before = 1
    )
}

qr_results <- purrr::map_dfr(taus, fit_one_tau)

safe_write_xlsx(
  list(
    coefficients = qr_results,
    model_info = tibble::tibble(
      formula = paste(deparse(formula), collapse = " "),
      n = nrow(model_data),
      note = "Quantile regression is a supplementary distributional analysis."
    )
  ),
  file.path(MODEL_DIR, "quantile_regression_results.xlsx")
)

plot_data <- qr_results |>
  dplyr::filter(term != "(Intercept)") |>
  dplyr::mutate(term = stringr::str_replace_all(term, "`", ""))

if (nrow(plot_data) > 0) {
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate, y = stats::reorder(term, estimate), colour = factor(tau))) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = palette_contract[["neutral_mid"]]) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = 0.15, linewidth = 0.35, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::geom_point(size = 1.2, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::scale_colour_manual(values = c("0.1" = palette_contract[["signal_blue"]], "0.5" = palette_contract[["neutral_dark"]], "0.9" = palette_contract[["accent_orange"]]), name = "Quantile") +
    ggplot2::labs(x = "Estimated difference in birth weight (kg)", y = NULL, title = "Quantile regression effects") +
    theme_pub()
  save_plot_pdf_png(p, file.path(FIGURE_DIR, "quantile_regression_effects"), width_mm = 160, height_mm = 130)
}

write_log("08_quantile_regression.R completed")
