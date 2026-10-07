source("scripts/00_config.R")

write_log("11_make_figures.R started")

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |> prepare_analysis_types()
sex_var <- sex_variable_name()
sex_title <- if (SEX_RULE_CONFIRMED) "male and female" else "G-coded and non-G-coded"
sex_values <- if (SEX_RULE_CONFIRMED) c(male = palette_contract[["signal_blue"]], female = palette_contract[["accent_orange"]]) else c(`G-coded` = palette_contract[["signal_blue"]], `non-G-coded` = palette_contract[["accent_orange"]])

flow <- openxlsx::read.xlsx(file.path(TABLE_DIR, "analysis_dataset_flow.xlsx"))
flow$step <- factor(flow$step, levels = rev(flow$step))
p1 <- ggplot2::ggplot(flow, ggplot2::aes(x = n, y = step)) +
  ggplot2::geom_col(fill = palette_contract[["signal_blue"]], width = 0.65) +
  ggplot2::geom_text(ggplot2::aes(label = scales::comma(n)), hjust = -0.1, size = 2.2) +
  ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.15))) +
  ggplot2::labs(x = "Records", y = NULL, title = "Data editing and analysis population") +
  theme_pub()
save_plot_pdf_png(p1, file.path(FIGURE_DIR, "Figure1_data_flow"), width_mm = 150, height_mm = 85)

thresholds <- openxlsx::read.xlsx(file.path(TABLE_DIR, "sex_specific_thresholds.xlsx"), sheet = "sex_specific_thresholds")
names(thresholds)[1] <- sex_var
p2 <- ggplot2::ggplot(singleton_tail, ggplot2::aes(x = birth_weight, fill = .data[[sex_var]], colour = .data[[sex_var]])) +
  ggplot2::geom_density(alpha = 0.18, linewidth = 0.45) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = P10, colour = .data[[sex_var]]), linetype = "dashed", linewidth = 0.35, show.legend = FALSE) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = P90, colour = .data[[sex_var]]), linetype = "dashed", linewidth = 0.35, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = sex_values, name = NULL) +
  ggplot2::scale_colour_manual(values = sex_values, name = NULL) +
  ggplot2::labs(x = "Birth weight (kg)", y = "Density", title = paste("Birth-weight distributions of", sex_title, "Simmental calves")) +
  theme_pub()
save_plot_pdf_png(p2, file.path(FIGURE_DIR, "Figure2_birth_weight_distribution_by_sex_or_code"), width_mm = 170, height_mm = 110)

emm_path <- file.path(PROCESSED_DIR, "lmm_emmeans_for_figures.csv")
if (file.exists(emm_path)) {
  emm <- read_csv_safe(emm_path)
  level_cols <- intersect(names(emm), c(sex_var, "parity_group", "calf_birth_year", "calf_birth_season", "dam_age_group"))
  if (length(level_cols) > 0 && nrow(emm) > 0) {
    emm <- emm |>
      dplyr::mutate(level = purrr::map2_chr(term, dplyr::row_number(), function(term_value, row_index) {
        if (term_value %in% names(emm)) as.character(emm[[term_value]][row_index]) else NA_character_
      })) |>
      dplyr::mutate(
        ci_low = dplyr::coalesce(.data[[intersect(c("lower.CL", "asymp.LCL"), names(emm))[1]]], emmean - 1.96 * SE),
        ci_high = dplyr::coalesce(.data[[intersect(c("upper.CL", "asymp.UCL"), names(emm))[1]]], emmean + 1.96 * SE)
      ) |>
      dplyr::filter(!is.na(level), term %in% c(sex_var, "parity_group", "calf_birth_year", "calf_birth_season"))
    p3 <- ggplot2::ggplot(emm, ggplot2::aes(x = level, y = emmean)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.18, linewidth = 0.3, colour = palette_contract[["neutral_mid"]]) +
      ggplot2::geom_point(size = 1.4, colour = palette_contract[["signal_blue"]]) +
      ggplot2::facet_wrap(~ term, scales = "free_x") +
      ggplot2::labs(x = NULL, y = "Estimated marginal mean birth weight (kg)", title = "Core fixed-effect estimates from the LMM") +
      theme_pub(base_size = 6.5) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_plot_pdf_png(p3, file.path(FIGURE_DIR, "Figure3_lmm_emmeans"), width_mm = 183, height_mm = 120)
  }
}

forest_path <- file.path(PROCESSED_DIR, "glmm_tail_forest_data.csv")
if (file.exists(forest_path)) {
  forest <- read_csv_safe(forest_path) |>
    dplyr::filter(is.finite(OR), is.finite(CI_low), is.finite(CI_high)) |>
    dplyr::mutate(term_clean = stringr::str_replace_all(term, c("parity_group" = "Parity ", "calf_birth_year" = "Year ", "calf_birth_season" = "Season ", "dam_age_group" = "Dam age ")))
  if (nrow(forest) > 0) {
    p4 <- ggplot2::ggplot(forest, ggplot2::aes(x = OR, y = stats::reorder(term_clean, OR), colour = phenotype)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = palette_contract[["neutral_mid"]]) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = 0.16, linewidth = 0.35, position = ggplot2::position_dodge(width = 0.55)) +
      ggplot2::geom_point(size = 1.2, position = ggplot2::position_dodge(width = 0.55)) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_colour_manual(values = c(lower_tail = palette_contract[["signal_blue"]], upper_tail = palette_contract[["accent_orange"]]), name = NULL) +
      ggplot2::labs(x = "Odds ratio (log scale)", y = NULL, title = "Lower- and upper-tail birth-weight phenotypes") +
      theme_pub(base_size = 6.5)
    save_plot_pdf_png(p4, file.path(FIGURE_DIR, "Figure4_tail_or_forest"), width_mm = 170, height_mm = 130)
    save_plot_pdf_png(p4, file.path(FIGURE_DIR, "glmm_tail_or_forest"), width_mm = 170, height_mm = 130)
  }
}

hist_plot <- ggplot2::ggplot(singleton_tail, ggplot2::aes(x = birth_weight)) +
  ggplot2::geom_histogram(binwidth = 1, fill = palette_contract[["neutral_light"]], colour = "white", linewidth = 0.1) +
  ggplot2::labs(x = "Birth weight (kg)", y = "Records", title = "Overall singleton birth-weight distribution") +
  theme_pub()
save_plot_pdf_png(hist_plot, file.path(FIGURE_DIR, "Supplementary_overall_birth_weight_histogram"), width_mm = 150, height_mm = 90)

parity_box <- ggplot2::ggplot(singleton_tail, ggplot2::aes(x = parity_group, y = birth_weight)) +
  ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.25, linewidth = 0.3, fill = palette_contract[["neutral_light"]]) +
  ggplot2::labs(x = "Parity group", y = "Birth weight (kg)", title = "Birth weight by parity group") +
  theme_pub()
save_plot_pdf_png(parity_box, file.path(FIGURE_DIR, "Supplementary_parity_boxplot"), width_mm = 135, height_mm = 90)

heat <- singleton_tail |>
  dplyr::group_by(calf_birth_year, calf_birth_season) |>
  dplyr::summarise(mean_birth_weight = mean(birth_weight, na.rm = TRUE), lower_rate = mean(lower_tail, na.rm = TRUE), upper_rate = mean(upper_tail, na.rm = TRUE), .groups = "drop")
heat_plot <- ggplot2::ggplot(heat, ggplot2::aes(x = calf_birth_season, y = calf_birth_year, fill = mean_birth_weight)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::scale_fill_gradient(low = "#E8F1F4", high = palette_contract[["accent_orange"]], name = "Mean kg") +
  ggplot2::labs(x = "Calving season", y = "Calving year", title = "Year-season mean birth weight") +
  theme_pub()
save_plot_pdf_png(heat_plot, file.path(FIGURE_DIR, "Supplementary_year_season_mean_heatmap"), width_mm = 120, height_mm = 90)

if (file.exists(file.path(MODEL_DIR, "lmm_birth_weight_final.rds"))) {
  lmm <- readRDS(file.path(MODEL_DIR, "lmm_birth_weight_final.rds"))
  diag_df <- tibble::tibble(fitted = stats::fitted(lmm), residual = stats::residuals(lmm))
  resid_plot <- ggplot2::ggplot(diag_df, ggplot2::aes(x = fitted, y = residual)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.35, colour = palette_contract[["neutral_mid"]]) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    ggplot2::labs(x = "Fitted birth weight (kg)", y = "Residual", title = "LMM residual diagnostics") +
    theme_pub()
  save_plot_pdf_png(resid_plot, file.path(FIGURE_DIR, "Supplementary_lmm_residual_diagnostics"), width_mm = 135, height_mm = 90)
}

write_log("11_make_figures.R completed")
