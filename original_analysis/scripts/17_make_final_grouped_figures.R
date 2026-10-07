source("scripts/00_config.R")

write_log("17_make_final_grouped_figures.R started")

FINAL_FIGURE_DIR <- file.path(FIGURE_DIR, "final_grouped")
FINAL_SOURCE_DIR <- file.path(FINAL_FIGURE_DIR, "source_data")
create_dir_if_missing(FINAL_FIGURE_DIR)
create_dir_if_missing(FINAL_SOURCE_DIR)

pal <- c(
  dark = "#252A32",
  mid = "#6F7782",
  grid = "#ECEFF2",
  female = "#D48642",
  male = "#2F6DAE",
  lower = "#3B6EA5",
  normal = "#E8E1D3",
  upper = "#C56A3A",
  green = "#6DA58D",
  purple = "#8A75A8",
  grey = "#AEB5BE"
)

theme_final <- function(base_size = 6.4, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.28, colour = pal[["dark"]]),
      axis.ticks = ggplot2::element_line(linewidth = 0.22, colour = pal[["dark"]]),
      axis.title = ggplot2::element_text(size = base_size, colour = pal[["dark"]]),
      axis.text = ggplot2::element_text(size = base_size - 0.8, colour = "#33373D"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 1.0),
      legend.key.size = grid::unit(3.1, "mm"),
      legend.spacing.x = grid::unit(1.4, "mm"),
      strip.background = ggplot2::element_rect(fill = "#F1F3F5", colour = NA),
      strip.text = ggplot2::element_text(size = base_size - 0.4, face = "bold", colour = pal[["dark"]]),
      plot.title = ggplot2::element_text(size = base_size + 0.3, face = "bold", colour = pal[["dark"]]),
      plot.subtitle = ggplot2::element_text(size = base_size - 1.0, colour = pal[["mid"]]),
      panel.grid.major = ggplot2::element_line(linewidth = 0.14, colour = pal[["grid"]]),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 4, 3, 3)
    )
}

save_final <- function(plot, filename, width_mm = 183, height_mm = 170, dpi = 450) {
  path <- file.path(FINAL_FIGURE_DIR, filename)
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggplot2::ggsave(paste0(path, ".pdf"), plot, width = w, height = h, device = grDevices::cairo_pdf)
  ragg::agg_png(paste0(path, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  grDevices::dev.off()
  svglite::svglite(paste0(path, ".svg"), width = w, height = h)
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(path, ".tiff"), width = w, height = h, units = "in", res = 600)
  print(plot)
  grDevices::dev.off()
  invisible(path)
}

sex_cols <- c(female = pal[["female"]], male = pal[["male"]])
tail_cols <- c("Lower tail" = pal[["lower"]], "Reference range" = pal[["normal"]], "Upper tail" = pal[["upper"]])

pct <- function(x, accuracy = 1) scales::percent(x, accuracy = accuracy)
safe_odds <- function(events, total) (events + 0.5) / (total - events + 0.5)
short_label <- function(x, width = 14) {
  x <- as.character(x)
  dplyr::if_else(nchar(x) > width, paste0(substr(x, 1, width - 1), "..."), x)
}

singleton <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |>
  prepare_analysis_types() |>
  dplyr::mutate(
    calf_sex = factor(as.character(calf_sex), levels = c("female", "male")),
    parity_group = factor(as.character(parity_group), levels = c("1", "2", "3", "4", "5", "6", ">=7")),
    calf_birth_year = factor(as.character(calf_birth_year), levels = sort(unique(as.character(calf_birth_year)))),
    calf_birth_season = factor(as.character(calf_birth_season), levels = c("Spring", "Summer", "Autumn", "Winter")),
    tail_class = dplyr::case_when(
      lower_tail == 1 ~ "Lower tail",
      upper_tail == 1 ~ "Upper tail",
      TRUE ~ "Reference range"
    ),
    tail_class = factor(tail_class, levels = names(tail_cols))
  )

all_valid <- read_csv_safe(file.path(PROCESSED_DIR, "all_valid_records.csv")) |>
  prepare_analysis_types()

# Shared summaries ------------------------------------------------------------

monthly <- read_csv_safe(file.path(FIGURE_DIR, "enhanced", "source_data", "figure2_monthly_summary.csv")) |>
  dplyr::mutate(
    calving_month = as.Date(calving_month),
    calf_sex = factor(calf_sex, levels = c("female", "male"))
  )

year_season <- singleton |>
  dplyr::group_by(calf_birth_year, calf_birth_season) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    lower_tail_rate = mean(lower_tail, na.rm = TRUE),
    upper_tail_rate = mean(upper_tail, na.rm = TRUE),
    tail_burden = lower_tail_rate + upper_tail_rate,
    .groups = "drop"
  )

parity_age <- singleton |>
  dplyr::group_by(parity_group, dam_age_group) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    lower_tail_rate = mean(lower_tail, na.rm = TRUE),
    upper_tail_rate = mean(upper_tail, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(monthly, file.path(FINAL_SOURCE_DIR, "shared_monthly_summary.csv"))
write_csv_safe(year_season, file.path(FINAL_SOURCE_DIR, "shared_year_season_summary.csv"))
write_csv_safe(parity_age, file.path(FINAL_SOURCE_DIR, "shared_parity_age_summary.csv"))

# Figure 1: cohort and distribution atlas ------------------------------------

flow <- data.frame(
  step = factor(
    c("Raw records", "Valid records", "Singleton main analysis", "Lower-tail events", "Upper-tail events"),
    levels = c("Raw records", "Valid records", "Singleton main analysis", "Lower-tail events", "Upper-tail events")
  ),
  n = c(nrow(readxl::read_excel(RAW_FILE, guess_max = 10000)), nrow(all_valid), nrow(singleton), sum(singleton$lower_tail == 1), sum(singleton$upper_tail == 1))
)

missing_key <- data.frame(
  field = c("birth_weight", "calf_sex", "dam_age_month", "parity", "service_sire", "maternal_grandsire"),
  missing_rate = c(
    mean(is.na(singleton$birth_weight)),
    mean(is.na(singleton$calf_sex)),
    mean(is.na(singleton$dam_age_month)),
    mean(is.na(singleton$parity)),
    mean(is.na(singleton$service_sire) | singleton$service_sire == ""),
    mean(is.na(singleton$maternal_grandsire) | singleton$maternal_grandsire == "")
  )
)

sex_count <- singleton |>
  dplyr::count(calf_sex, name = "n") |>
  dplyr::mutate(prop = n / sum(n))

age_response_data <- singleton |>
  dplyr::filter(!is.na(dam_age_month), !is.na(birth_weight), !is.na(calf_sex))
age_limits <- as.numeric(stats::quantile(age_response_data$dam_age_month, c(0.03, 0.97), na.rm = TRUE))
age_model_data <- age_response_data |>
  dplyr::filter(dam_age_month >= age_limits[[1]], dam_age_month <= age_limits[[2]])
age_fit <- mgcv::gam(
  birth_weight ~ calf_sex + s(dam_age_month, by = calf_sex, k = 8),
  data = age_model_data,
  method = "REML"
)
age_grid <- expand.grid(
  dam_age_month = seq(age_limits[[1]], age_limits[[2]], length.out = 150),
  calf_sex = levels(age_model_data$calf_sex)
)
age_pred <- stats::predict(age_fit, newdata = age_grid, se.fit = TRUE)
age_grid$pred <- as.numeric(age_pred$fit)
age_grid$lower <- age_grid$pred - 1.96 * as.numeric(age_pred$se.fit)
age_grid$upper <- age_grid$pred + 1.96 * as.numeric(age_pred$se.fit)

p1A <- ggplot2::ggplot(flow, ggplot2::aes(x = n, y = step)) +
  ggplot2::geom_col(fill = "#4F86B8", width = 0.68) +
  ggplot2::geom_text(ggplot2::aes(label = scales::comma(n)), hjust = -0.08, size = 2.0) +
  ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.15))) +
  ggplot2::labs(x = "Records", y = NULL, title = "A  Analysis set") +
  theme_final()

p1B <- ggplot2::ggplot(missing_key, ggplot2::aes(x = stats::reorder(field, missing_rate), y = missing_rate)) +
  ggplot2::geom_col(fill = "#AEB5BE", width = 0.65) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = NULL, y = "Missing", title = "B  Key-field missingness") +
  theme_final()

p1C <- ggplot2::ggplot(sex_count, ggplot2::aes(x = calf_sex, y = prop, fill = calf_sex)) +
  ggplot2::geom_col(width = 0.58, colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(scales::comma(n), "\n", pct(prop, 0.1))), vjust = -0.25, size = 1.9, lineheight = 0.85) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.58), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Proportion", title = "C  Sex composition") +
  theme_final() +
  ggplot2::theme(legend.position = "none")

p1D <- ggplot2::ggplot(singleton, ggplot2::aes(x = birth_weight, colour = calf_sex, fill = calf_sex)) +
  ggplot2::geom_density(alpha = 0.14, linewidth = 0.62) +
  ggplot2::scale_colour_manual(values = sex_cols) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::coord_cartesian(xlim = c(20, 70)) +
  ggplot2::labs(x = "Birth weight (kg)", y = "Density", title = "D  Sex-specific density") +
  theme_final()

ridge_parity <- read_csv_safe(file.path(FIGURE_DIR, "advanced_candidates", "source_data", "advanced_D_ridge_parity_source.csv"))
ridge_season <- read_csv_safe(file.path(FIGURE_DIR, "advanced_candidates", "source_data", "advanced_D_ridge_season_source.csv"))

p1E <- ggplot2::ggplot(ridge_parity, ggplot2::aes(x = x, group = interaction(calf_sex, group_label))) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = group_index, ymax = group_index + density_scaled, fill = calf_sex), alpha = 0.62, colour = "white", linewidth = 0.15) +
  ggplot2::facet_grid(. ~ calf_sex) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::scale_y_continuous(breaks = 1:7, labels = c("1", "2", "3", "4", "5", "6", ">=7")) +
  ggplot2::coord_cartesian(xlim = c(22, 66)) +
  ggplot2::labs(x = "Birth weight (kg)", y = "Parity", title = "E  Parity ridges") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(legend.position = "none")

p1F <- ggplot2::ggplot(ridge_season, ggplot2::aes(x = x, group = interaction(calf_sex, group_label))) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = group_index, ymax = group_index + density_scaled, fill = calf_sex), alpha = 0.62, colour = "white", linewidth = 0.15) +
  ggplot2::facet_grid(. ~ calf_sex) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::scale_y_continuous(breaks = 1:4, labels = c("Spring", "Summer", "Autumn", "Winter")) +
  ggplot2::coord_cartesian(xlim = c(22, 66)) +
  ggplot2::labs(x = "Birth weight (kg)", y = "Season", title = "F  Seasonal ridges") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(legend.position = "none")

p1G <- ggplot2::ggplot(age_response_data, ggplot2::aes(x = dam_age_month, fill = calf_sex)) +
  ggplot2::geom_histogram(bins = 36, alpha = 0.72, position = "identity", colour = "white", linewidth = 0.12) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::coord_cartesian(xlim = c(20, 130)) +
  ggplot2::labs(x = "Dam age (months)", y = "Records", title = "G  Dam-age distribution") +
  theme_final()

p1H <- ggplot2::ggplot(age_grid, ggplot2::aes(x = dam_age_month, y = pred, colour = calf_sex, fill = calf_sex)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.13, linewidth = 0, colour = NA, show.legend = FALSE) +
  ggplot2::geom_line(linewidth = 0.62) +
  ggplot2::scale_colour_manual(values = sex_cols) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::labs(x = "Dam age (months)", y = "Predicted BW (kg)", title = "H  Age response") +
  theme_final()

p1I <- ggplot2::ggplot(parity_age, ggplot2::aes(x = dam_age_group, y = parity_group, fill = mean_birth_weight)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(n >= 20, sprintf("%.1f", mean_birth_weight), "")), size = 1.8) +
  ggplot2::scale_fill_gradientn(colours = c("#EEF3F6", "#B8D7E8", "#5E9BC6", "#244F79")) +
  ggplot2::labs(x = "Dam-age group", y = "Parity", fill = "Mean BW", title = "I  Parity-age mean BW") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "right")

fig1 <- (p1A | p1B | p1C) / (p1D | p1E | p1F) / (p1G | p1H | p1I) +
  patchwork::plot_layout(heights = c(0.9, 1.2, 1.1))
save_final(fig1, "Final_Figure1_cohort_distribution_atlas", height_mm = 190)

# Figure 2: temporal and seasonal tail dynamics --------------------------------

p2A <- ggplot2::ggplot(monthly, ggplot2::aes(x = calving_month, y = mean_birth_weight, colour = calf_sex, fill = calf_sex)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high), alpha = 0.14, colour = NA) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 0.75) +
  ggplot2::scale_colour_manual(values = sex_cols) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ggplot2::labs(x = NULL, y = "Mean BW (kg)", title = "A  Monthly birth-weight trajectory") +
  theme_final()

p2B <- ggplot2::ggplot(monthly, ggplot2::aes(x = calving_month, y = n, fill = calf_sex)) +
  ggplot2::geom_col(width = 23, alpha = 0.85) +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ggplot2::labs(x = NULL, y = "Records", title = "B  Monthly record volume") +
  theme_final()

monthly_long <- monthly |>
  tidyr::pivot_longer(c(lower_tail_rate, upper_tail_rate), names_to = "tail", values_to = "rate") |>
  dplyr::mutate(tail = dplyr::recode(tail, lower_tail_rate = "Lower-tail", upper_tail_rate = "Upper-tail"))

p2C <- ggplot2::ggplot(monthly_long, ggplot2::aes(x = calving_month, y = rate, colour = tail)) +
  ggplot2::geom_line(linewidth = 0.45) +
  ggplot2::facet_grid(. ~ calf_sex) +
  ggplot2::scale_colour_manual(values = c("Lower-tail" = pal[["lower"]], "Upper-tail" = pal[["upper"]])) +
  ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = NULL, y = "Monthly rate", title = "C  Monthly tail frequency") +
  theme_final(base_size = 6.0)

p2D <- ggplot2::ggplot(year_season, ggplot2::aes(x = calf_birth_season, y = calf_birth_year, fill = mean_birth_weight)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", mean_birth_weight)), size = 1.8) +
  ggplot2::scale_fill_gradientn(colours = c("#EEF3F6", "#B8D7E8", "#5E9BC6", "#244F79")) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Mean BW", title = "D  Year-season mean BW") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")

heat_rate <- function(var, title, fill_col) {
  ggplot2::ggplot(year_season, ggplot2::aes(x = calf_birth_season, y = calf_birth_year, fill = .data[[var]])) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = pct(.data[[var]], 0.1)), size = 1.65) +
    ggplot2::scale_fill_gradientn(colours = c("#F7F3EE", "#E8C5A2", fill_col, "#74312B"), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Rate", title = title) +
    theme_final(base_size = 5.9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
}
p2E <- heat_rate("lower_tail_rate", "E  Lower-tail rate", pal[["lower"]])
p2F <- heat_rate("upper_tail_rate", "F  Upper-tail rate", pal[["upper"]])
p2G <- heat_rate("tail_burden", "G  Combined tail burden", pal[["upper"]])

season_year_n <- year_season |>
  dplyr::mutate(label = paste0(scales::comma(n), "\n", pct(tail_burden, 0.1)))
p2H <- ggplot2::ggplot(season_year_n, ggplot2::aes(x = calf_birth_season, y = calf_birth_year, fill = n)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = scales::comma(n)), size = 1.7) +
  ggplot2::scale_fill_gradientn(colours = c("#F4F6F7", "#B9D2C5", "#6DA58D", "#2C6550")) +
  ggplot2::labs(x = NULL, y = NULL, title = "H  Year-season volume") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")

p2I <- ggplot2::ggplot(year_season, ggplot2::aes(x = lower_tail_rate, y = upper_tail_rate, size = n, colour = calf_birth_year)) +
  ggplot2::geom_point(alpha = 0.86) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_size_continuous(range = c(1.2, 4.0)) +
  ggplot2::labs(x = "Lower-tail rate", y = "Upper-tail rate", title = "I  Tail-balance scatter") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(legend.position = "none")

fig2 <- (p2A | p2B | p2C) / (p2D | p2E | p2F) / (p2G | p2H | p2I) +
  patchwork::plot_layout(heights = c(1.05, 1, 1))
save_final(fig2, "Final_Figure2_temporal_tail_dynamics", height_mm = 195)

# Figure 3: adjusted model evidence ------------------------------------------

emm <- read_csv_safe(file.path(PROCESSED_DIR, "lmm_emmeans_for_figures.csv"))
forest <- read_csv_safe(file.path(FIGURE_DIR, "enhanced", "source_data", "figure4_glmm_emmeans_or.csv")) |>
  dplyr::mutate(
    phenotype = factor(phenotype, levels = c("Lower-tail phenotype", "Upper-tail phenotype")),
    term_group = factor(term_group, levels = c("Parity", "Year", "Season")),
    direction = dplyr::case_when(CI_low > 1 ~ "Higher odds", CI_high < 1 ~ "Lower odds", TRUE ~ "Uncertain")
  )

plot_emm <- function(data, xvar, title) {
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[xvar]], y = emmean)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.12, linewidth = 0.35, colour = pal[["mid"]]) +
    ggplot2::geom_point(size = 1.7, colour = pal[["dark"]]) +
    ggplot2::labs(x = NULL, y = "Adjusted mean BW (kg)", title = title) +
    theme_final(base_size = 5.9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}
p3A <- plot_emm(emm |> dplyr::filter(term == "calf_sex"), "calf_sex", "A  Sex")
p3B <- plot_emm(emm |> dplyr::filter(term == "parity_group"), "parity_group", "B  Parity")
p3C <- plot_emm(emm |> dplyr::filter(term == "calf_birth_year"), "calf_birth_year", "C  Year")
p3D <- plot_emm(emm |> dplyr::filter(term == "calf_birth_season"), "calf_birth_season", "D  Season")

forest_panel <- function(phen, title) {
  dat <- forest |>
    dplyr::filter(phenotype == phen) |>
    dplyr::mutate(term_clean = factor(term_clean, levels = rev(unique(term_clean))))
  ggplot2::ggplot(dat, ggplot2::aes(x = OR, y = term_clean)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "22", linewidth = 0.3, colour = pal[["mid"]]) +
    ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = term_clean, colour = direction), linewidth = 0.48) +
    ggplot2::geom_point(ggplot2::aes(colour = direction), size = 1.35) +
    ggplot2::facet_grid(term_group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
    ggplot2::scale_colour_manual(values = c("Lower odds" = pal[["lower"]], "Uncertain" = pal[["grey"]], "Higher odds" = pal[["upper"]])) +
    ggplot2::labs(x = "Odds ratio", y = NULL, title = title) +
    theme_final(base_size = 5.7) +
    ggplot2::theme(legend.position = "none")
}
p3E <- forest_panel("Lower-tail phenotype", "E  Lower-tail GLMM")
p3F <- forest_panel("Upper-tail phenotype", "F  Upper-tail GLMM")

fig3 <- (p3A | p3B | p3C | p3D) / (p3E | p3F) +
  patchwork::plot_layout(heights = c(0.78, 1.8))
save_final(fig3, "Final_Figure3_adjusted_model_evidence", height_mm = 160)

# Figure 4: tail phenotype and robustness ------------------------------------

tail_stack_year <- read_csv_safe(file.path(FIGURE_DIR, "advanced_candidates", "source_data", "advanced_A_year_season_tail_stack.csv")) |>
  dplyr::mutate(tail_class = factor(tail_class, levels = names(tail_cols)))
tail_stack_parity <- read_csv_safe(file.path(FIGURE_DIR, "advanced_candidates", "source_data", "advanced_A_parity_sex_tail_stack.csv")) |>
  dplyr::mutate(tail_class = factor(tail_class, levels = names(tail_cols)))
effect_risk <- read_csv_safe(file.path(FIGURE_DIR, "advanced_candidates", "source_data", "advanced_B_descriptive_effect_risk_map.csv"))

p4A <- ggplot2::ggplot(tail_stack_year, ggplot2::aes(x = calf_birth_season, y = n, fill = tail_class)) +
  ggplot2::geom_col(position = "fill", width = 0.7, colour = "white", linewidth = 0.2) +
  ggplot2::facet_grid(. ~ calf_birth_year) +
  ggplot2::scale_fill_manual(values = tail_cols) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = NULL, y = "Proportion", title = "A  Year-season tail composition") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

p4B <- ggplot2::ggplot(tail_stack_parity, ggplot2::aes(x = parity_group, y = n, fill = tail_class)) +
  ggplot2::geom_col(position = "fill", width = 0.7, colour = "white", linewidth = 0.2) +
  ggplot2::facet_grid(. ~ calf_sex) +
  ggplot2::scale_fill_manual(values = tail_cols) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Parity", y = "Proportion", title = "B  Parity-sex tail composition") +
  theme_final(base_size = 5.9) +
  ggplot2::theme(legend.position = "none")

effect_panel <- function(phen, title) {
  lab <- effect_risk |>
    dplyr::filter(phenotype == phen) |>
    dplyr::slice_max(score, n = 6, with_ties = FALSE)
  ggplot2::ggplot(effect_risk |> dplyr::filter(phenotype == phen), ggplot2::aes(x = mean_delta, y = log2_or)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
    ggplot2::geom_vline(xintercept = 0, linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
    ggplot2::geom_point(ggplot2::aes(size = n, colour = domain), alpha = 0.82) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = label, colour = domain), size = 1.65, show.legend = FALSE, max.overlaps = Inf, segment.size = 0.15) +
    ggplot2::scale_size_continuous(range = c(1.1, 4.2)) +
    ggplot2::labs(x = "Mean BW difference (kg)", y = "log2 odds ratio", title = title) +
    theme_final(base_size = 5.8) +
    ggplot2::theme(legend.position = "none")
}
p4C <- effect_panel("Lower-tail phenotype", "C  Lower-tail effect-risk")
p4D <- effect_panel("Upper-tail phenotype", "D  Upper-tail effect-risk")

qcoef <- openxlsx::read.xlsx(file.path(MODEL_DIR, "quantile_regression_results.xlsx"), sheet = "coefficients") |>
  dplyr::filter(term != "(Intercept)", grepl("calf_sex|parity_group|calf_birth_season", term)) |>
  dplyr::mutate(term = dplyr::case_when(
    term == "calf_sex1" ~ "Sex",
    grepl("parity_group", term) ~ gsub("parity_group", "Parity ", term),
    grepl("calf_birth_season", term) ~ gsub("calf_birth_season", "Season ", term),
    TRUE ~ term
  ))
p4E <- ggplot2::ggplot(qcoef, ggplot2::aes(x = estimate, y = stats::reorder(term, estimate), colour = factor(tau))) +
  ggplot2::geom_vline(xintercept = 0, linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
  ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = term), linewidth = 0.36, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::geom_point(size = 1.1, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::labs(x = "Quantile coefficient (kg)", y = NULL, colour = "Tau", title = "E  Quantile-regression effects") +
  theme_final(base_size = 5.7)

sens <- openxlsx::read.xlsx(file.path(TABLE_DIR, "sensitivity_summary.xlsx"), sheet = 1) |>
  dplyr::mutate(
    tail_direction = dplyr::if_else(grepl("upper", phenotype), "Upper tail", "Lower tail"),
    label_compact = gsub(" lower$| upper$", "", label),
    label_compact = dplyr::recode(
      label_compact,
      "all valid records plus birth_type" = "birth type included",
      "remove birth_weight <25 or >65" = "25-65 kg only",
      "parity_plus_age_group" = "parity + age group",
      "parity_plus_age_cont" = "parity + age spline",
      "parity_only" = "parity only",
      "age_group_only" = "age group only",
      "service_sire label adjustment" = "service sire adjusted",
      "maternal_grandsire label adjustment" = "MGS adjusted",
      "year season interaction" = "year x season"
    ),
    label_compact = stats::reorder(label_compact, event_rate)
  )
p4F <- ggplot2::ggplot(sens, ggplot2::aes(x = event_rate, y = label_compact, colour = tail_direction, shape = tail_direction)) +
  ggplot2::geom_point(size = 1.15, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_colour_manual(values = c("Lower tail" = pal[["lower"]], "Upper tail" = pal[["upper"]])) +
  ggplot2::scale_shape_manual(values = c("Lower tail" = 16, "Upper tail" = 17)) +
  ggplot2::labs(x = "Event rate", y = NULL, title = "F  Sensitivity event rates") +
  theme_final(base_size = 5.4) +
  ggplot2::theme(legend.position = "none")

fig4 <- (p4A | p4B) / (p4C | p4D) / (p4E | p4F) +
  patchwork::plot_layout(heights = c(1.0, 1.0, 1.15))
save_final(fig4, "Final_Figure4_tail_robustness_atlas", height_mm = 185)

# Figure 5: exploratory labels and response surfaces --------------------------

label_volcano <- read_csv_safe(file.path(FIGURE_DIR, "advanced_candidates", "reference_style", "source_data", "panel_B_exploratory_label_volcano.csv")) |>
  dplyr::mutate(
    direction = dplyr::case_when(q_value < 0.05 & delta_kg >= 1 ~ "Higher BW", q_value < 0.05 & delta_kg <= -1 ~ "Lower BW", TRUE ~ "Not highlighted")
  )
volcano_labels <- label_volcano |>
  dplyr::arrange(q_value, dplyr::desc(abs(delta_kg))) |>
  dplyr::slice_head(n = 8)

p5A <- p1H + ggplot2::labs(title = "A  Dam-age response")
p5B <- ggplot2::ggplot(label_volcano, ggplot2::aes(x = delta_kg, y = neg_log10_q)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
  ggplot2::geom_point(ggplot2::aes(fill = direction, shape = label_type, size = n), colour = "white", stroke = 0.25, alpha = 0.82) +
  ggrepel::geom_text_repel(data = volcano_labels, ggplot2::aes(label = label), size = 1.65, max.overlaps = Inf, segment.size = 0.15, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = c("Lower BW" = pal[["lower"]], "Not highlighted" = pal[["grey"]], "Higher BW" = pal[["upper"]])) +
  ggplot2::scale_shape_manual(values = c("Service sire" = 21, "Maternal grandsire" = 24)) +
  ggplot2::scale_size_continuous(range = c(1.1, 3.8)) +
  ggplot2::labs(x = "Mean difference (kg)", y = expression(-log[10]("BH-adjusted P")), title = "B  Exploratory label screen") +
  theme_final(base_size = 5.8) +
  ggplot2::theme(legend.position = "none")

rankings <- read_csv_safe(file.path(FIGURE_DIR, "enhanced", "source_data", "figure5_label_rankings.csv"))
rank_panel <- function(type_pattern, title, point_col) {
  dat <- rankings |>
    dplyr::filter(grepl(type_pattern, label_type, ignore.case = TRUE)) |>
    dplyr::slice_max(raw_mean_birth_weight, n = 12, with_ties = FALSE) |>
    dplyr::mutate(label_short = short_label(label_plot, 12))
  global_mean <- mean(singleton$birth_weight, na.rm = TRUE)
  ggplot2::ggplot(dat, ggplot2::aes(x = raw_mean_birth_weight, y = stats::reorder(label_short, raw_mean_birth_weight), size = n)) +
    ggplot2::geom_vline(xintercept = global_mean, linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
    ggplot2::geom_segment(ggplot2::aes(x = low, xend = high, yend = label_short), linewidth = 0.35, colour = pal[["grey"]]) +
    ggplot2::geom_point(colour = point_col, alpha = 0.88) +
    ggplot2::scale_size_continuous(range = c(1.1, 3.2)) +
    ggplot2::labs(x = "Mean BW (kg)", y = NULL, title = title) +
    theme_final(base_size = 5.4) +
    ggplot2::theme(legend.position = "none")
}
p5C <- rank_panel("service", "C  Service-sire labels", pal[["male"]])
p5D <- rank_panel("maternal", "D  Maternal-grandsire labels", pal[["purple"]])

surface <- read_csv_safe(file.path(FIGURE_DIR, "enhanced", "3d", "source_data", "3d_surface_source.csv")) |>
  dplyr::mutate(calf_sex = factor(calf_sex, levels = c("female", "male")))
surface_panel <- function(sex_value, title) {
  ggplot2::ggplot(surface |> dplyr::filter(calf_sex == sex_value), ggplot2::aes(x = dam_age_month, y = calf_birth_month_num, fill = pred)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = pred), colour = "white", linewidth = 0.18, alpha = 0.65) +
    ggplot2::scale_fill_gradientn(colours = if (sex_value == "female") c("#F6E9DB", "#D48642", "#7E3F24") else c("#E6EEF7", "#2F6DAE", "#17375E")) +
    ggplot2::labs(x = "Dam age (months)", y = "Calving month", fill = "Pred BW", title = title) +
    theme_final(base_size = 5.8) +
    ggplot2::theme(legend.position = "right")
}
p5E <- surface_panel("female", "E  Female response surface")
p5F <- surface_panel("male", "F  Male response surface")

fig5 <- (p5A | p5B) / (p5C | p5D) / (p5E | p5F) +
  patchwork::plot_layout(heights = c(1.05, 1.0, 1.05))
save_final(fig5, "Final_Figure5_exploratory_response_surfaces", height_mm = 180)

write_log("17_make_final_grouped_figures.R completed")
