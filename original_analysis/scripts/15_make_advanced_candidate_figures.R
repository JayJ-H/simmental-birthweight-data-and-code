source("scripts/00_config.R")

write_log("15_make_advanced_candidate_figures.R started")

CANDIDATE_DIR <- file.path(FIGURE_DIR, "advanced_candidates")
CANDIDATE_SOURCE_DIR <- file.path(CANDIDATE_DIR, "source_data")
create_dir_if_missing(CANDIDATE_DIR)
create_dir_if_missing(CANDIDATE_SOURCE_DIR)

candidate_palette <- c(
  neutral_dark = "#252A32",
  neutral_mid = "#6F7782",
  neutral_light = "#E7E9ED",
  lower = "#3B6EA5",
  normal = "#E8E1D3",
  upper = "#C56A3A",
  female = "#D48642",
  male = "#2F6DAE",
  parity = "#6379A8",
  season = "#70A58D",
  year = "#B96E55",
  sex = "#8E78AA"
)

theme_candidate <- function(base_size = 7.2, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.32, colour = candidate_palette[["neutral_dark"]]),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = candidate_palette[["neutral_dark"]]),
      axis.title = ggplot2::element_text(size = base_size, colour = candidate_palette[["neutral_dark"]]),
      axis.text = ggplot2::element_text(size = base_size - 0.7, colour = "#2F3338"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 0.8),
      strip.background = ggplot2::element_rect(fill = "#F1F3F5", colour = NA),
      strip.text = ggplot2::element_text(size = base_size - 0.2, face = "bold", colour = candidate_palette[["neutral_dark"]]),
      plot.title = ggplot2::element_text(size = base_size + 1.3, face = "bold", colour = candidate_palette[["neutral_dark"]]),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.3, colour = candidate_palette[["neutral_mid"]]),
      plot.caption = ggplot2::element_text(size = base_size - 1.2, colour = candidate_palette[["neutral_mid"]], hjust = 0),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.18, colour = "#ECEFF2"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 7, 5, 5)
    )
}

save_candidate <- function(plot, filename, width_mm = 183, height_mm = 130, dpi = 450) {
  path <- file.path(CANDIDATE_DIR, filename)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  ggplot2::ggsave(paste0(path, ".pdf"), plot, width = width_in, height = height_in, device = grDevices::cairo_pdf)
  ragg::agg_png(paste0(path, ".png"), width = width_in, height = height_in, units = "in", res = dpi)
  print(plot)
  grDevices::dev.off()
  svglite::svglite(paste0(path, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(path, ".tiff"), width = width_in, height = height_in, units = "in", res = 600)
  print(plot)
  grDevices::dev.off()

  invisible(path)
}

clean_percent <- function(x, accuracy = 0.1) {
  scales::percent(x, accuracy = accuracy)
}

safe_odds <- function(events, total) {
  (events + 0.5) / (total - events + 0.5)
}

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |>
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
    tail_class = factor(tail_class, levels = c("Lower tail", "Reference range", "Upper tail"))
  )

tail_fill <- c(
  "Lower tail" = candidate_palette[["lower"]],
  "Reference range" = candidate_palette[["normal"]],
  "Upper tail" = candidate_palette[["upper"]]
)

domain_fill <- c(
  Sex = candidate_palette[["sex"]],
  Parity = candidate_palette[["parity"]],
  Season = candidate_palette[["season"]],
  Year = candidate_palette[["year"]]
)

sex_fill <- c(female = candidate_palette[["female"]], male = candidate_palette[["male"]])

# Figure A: tail phenotype landscape -----------------------------------------

year_season_tail <- singleton_tail |>
  dplyr::count(calf_birth_year, calf_birth_season, tail_class, name = "n") |>
  tidyr::complete(
    calf_birth_year,
    calf_birth_season,
    tail_class,
    fill = list(n = 0)
  ) |>
  dplyr::group_by(calf_birth_year, calf_birth_season) |>
  dplyr::mutate(
    total = sum(n),
    proportion = dplyr::if_else(total > 0, n / total, NA_real_),
    label = dplyr::if_else(total > 0 & tail_class != "Reference range" & proportion >= 0.055, clean_percent(proportion), "")
  ) |>
  dplyr::ungroup()

parity_sex_tail <- singleton_tail |>
  dplyr::count(calf_sex, parity_group, tail_class, name = "n") |>
  tidyr::complete(calf_sex, parity_group, tail_class, fill = list(n = 0)) |>
  dplyr::group_by(calf_sex, parity_group) |>
  dplyr::mutate(
    total = sum(n),
    proportion = dplyr::if_else(total > 0, n / total, NA_real_),
    label = dplyr::if_else(tail_class != "Reference range" & proportion >= 0.07, clean_percent(proportion), "")
  ) |>
  dplyr::ungroup()

tail_burden <- singleton_tail |>
  dplyr::group_by(calf_birth_year, calf_birth_season) |>
  dplyr::summarise(
    n = dplyr::n(),
    lower_tail_rate = mean(lower_tail, na.rm = TRUE),
    upper_tail_rate = mean(upper_tail, na.rm = TRUE),
    tail_burden = lower_tail_rate + upper_tail_rate,
    .groups = "drop"
  ) |>
  dplyr::mutate(label = paste0(clean_percent(tail_burden), "\n", "n=", n))

write_csv_safe(year_season_tail, file.path(CANDIDATE_SOURCE_DIR, "advanced_A_year_season_tail_stack.csv"))
write_csv_safe(parity_sex_tail, file.path(CANDIDATE_SOURCE_DIR, "advanced_A_parity_sex_tail_stack.csv"))
write_csv_safe(tail_burden, file.path(CANDIDATE_SOURCE_DIR, "advanced_A_tail_burden_heatmap.csv"))

pA_year <- ggplot2::ggplot(year_season_tail, ggplot2::aes(x = calf_birth_season, y = n, fill = tail_class)) +
  ggplot2::geom_col(position = "fill", width = 0.74, colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    position = ggplot2::position_fill(vjust = 0.5),
    size = 2.0,
    colour = "#1F2328",
    lineheight = 0.85
  ) +
  ggplot2::facet_grid(. ~ calf_birth_year) +
  ggplot2::scale_fill_manual(values = tail_fill) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  ggplot2::labs(
    x = NULL,
    y = "Record proportion",
    title = "Calendar landscape of birth-weight tail phenotypes",
    subtitle = "Bars show lower-tail, reference-range, and upper-tail proportions within each year-season stratum"
  ) +
  theme_candidate(base_size = 7.0) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

pA_parity <- ggplot2::ggplot(parity_sex_tail, ggplot2::aes(x = parity_group, y = n, fill = tail_class)) +
  ggplot2::geom_col(position = "fill", width = 0.74, colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    position = ggplot2::position_fill(vjust = 0.5),
    size = 1.85,
    colour = "#1F2328"
  ) +
  ggplot2::facet_grid(. ~ calf_sex) +
  ggplot2::scale_fill_manual(values = tail_fill) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  ggplot2::labs(x = "Parity group", y = "Record proportion", title = "Tail composition by parity and sex") +
  theme_candidate(base_size = 6.7) +
  ggplot2::theme(legend.position = "none")

pA_heat <- ggplot2::ggplot(tail_burden, ggplot2::aes(x = calf_birth_season, y = calf_birth_year, fill = tail_burden)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 2.0, lineheight = 0.82, colour = "#1F2328") +
  ggplot2::scale_fill_gradientn(colours = c("#F3F0EA", "#DDAE7B", "#B85B35", "#6D2F2A"), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Tail burden", title = "Combined tail burden") +
  theme_candidate(base_size = 6.7) +
  ggplot2::theme(
    legend.position = "right",
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  )

advanced_A <- pA_year / (pA_parity | pA_heat) +
  patchwork::plot_layout(heights = c(1.2, 1), widths = c(1.55, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 8, face = "bold"))

save_candidate(advanced_A, "Advanced_FigureA_tail_phenotype_landscape", width_mm = 183, height_mm = 142)

# Figure B: descriptive effect-risk map ---------------------------------------

overall_n <- nrow(singleton_tail)
overall_mean <- mean(singleton_tail$birth_weight, na.rm = TRUE)
overall_lower_odds <- safe_odds(sum(singleton_tail$lower_tail, na.rm = TRUE), overall_n)
overall_upper_odds <- safe_odds(sum(singleton_tail$upper_tail, na.rm = TRUE), overall_n)

factor_summary <- function(data, var, domain, label_fun = identity) {
  data |>
    dplyr::filter(!is.na(.data[[var]])) |>
    dplyr::group_by(level = .data[[var]]) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_birth_weight = mean(birth_weight, na.rm = TRUE),
      lower_events = sum(lower_tail, na.rm = TRUE),
      upper_events = sum(upper_tail, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      domain = domain,
      label = label_fun(as.character(level)),
      mean_delta = mean_birth_weight - overall_mean,
      lower_log2_or = log2(safe_odds(lower_events, n) / overall_lower_odds),
      upper_log2_or = log2(safe_odds(upper_events, n) / overall_upper_odds)
    )
}

effect_risk <- dplyr::bind_rows(
  factor_summary(singleton_tail, "calf_sex", "Sex", function(x) ifelse(x == "female", "Female", "Male")),
  factor_summary(singleton_tail, "parity_group", "Parity", function(x) paste("Parity", x)),
  factor_summary(singleton_tail, "calf_birth_season", "Season"),
  factor_summary(singleton_tail, "calf_birth_year", "Year")
) |>
  dplyr::mutate(
    label = gsub("Parity >=7", "Parity >=7", label, fixed = TRUE),
    domain = factor(domain, levels = c("Sex", "Parity", "Season", "Year"))
  )

effect_risk_long <- effect_risk |>
  tidyr::pivot_longer(
    c(lower_log2_or, upper_log2_or),
    names_to = "phenotype",
    values_to = "log2_or"
  ) |>
  dplyr::mutate(
    phenotype = dplyr::recode(
      phenotype,
      lower_log2_or = "Lower-tail phenotype",
      upper_log2_or = "Upper-tail phenotype"
    ),
    score = abs(mean_delta) / stats::sd(mean_delta, na.rm = TRUE) + abs(log2_or)
  )

effect_labels <- effect_risk_long |>
  dplyr::group_by(phenotype) |>
  dplyr::slice_max(score, n = 8, with_ties = FALSE) |>
  dplyr::ungroup()

write_csv_safe(effect_risk_long, file.path(CANDIDATE_SOURCE_DIR, "advanced_B_descriptive_effect_risk_map.csv"))

pB <- ggplot2::ggplot(effect_risk_long, ggplot2::aes(x = mean_delta, y = log2_or)) +
  ggplot2::annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#F3F7FB", alpha = 0.75) +
  ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#FBF3EE", alpha = 0.75) +
  ggplot2::geom_hline(yintercept = 0, linetype = "22", linewidth = 0.28, colour = "#7F8790") +
  ggplot2::geom_vline(xintercept = 0, linetype = "22", linewidth = 0.28, colour = "#7F8790") +
  ggplot2::geom_point(
    ggplot2::aes(colour = domain, size = n),
    alpha = 0.86,
    stroke = 0.25
  ) +
  ggplot2::geom_text(
    data = effect_labels,
    ggplot2::aes(label = label, colour = domain),
    size = 2.15,
    hjust = -0.08,
    vjust = -0.25,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  ggplot2::facet_grid(. ~ phenotype) +
  ggplot2::scale_colour_manual(values = domain_fill) +
  ggplot2::scale_size_continuous(range = c(1.7, 6.2), breaks = c(500, 1500, 3000, 5000), labels = scales::comma) +
  ggplot2::coord_cartesian(xlim = c(-2.7, 2.7), ylim = c(-2.2, 2.2), clip = "off") +
  ggplot2::labs(
    x = "Mean birth-weight difference from population mean (kg)",
    y = "Observed log2 odds ratio vs population",
    size = "Records",
    title = "Phenotype effect-risk map",
    subtitle = "Each point is one sex, parity, season, or year stratum; position combines mean-weight shift and tail enrichment"
  ) +
  theme_candidate(base_size = 7.1) +
  ggplot2::theme(
    legend.box = "horizontal",
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.18, colour = "#ECEFF2")
  )

save_candidate(pB, "Advanced_FigureB_effect_risk_map", width_mm = 183, height_mm = 112)

# Figure C: adjusted model evidence forest ------------------------------------

forest_path <- file.path(ENHANCED_SOURCE_DIR <- file.path(FIGURE_DIR, "enhanced", "source_data"), "figure4_glmm_emmeans_or.csv")
if (!file.exists(forest_path)) {
  stop("Missing enhanced GLMM source data: ", forest_path, call. = FALSE)
}

forest <- read_csv_safe(forest_path) |>
  dplyr::mutate(
    phenotype = factor(phenotype, levels = c("Lower-tail phenotype", "Upper-tail phenotype")),
    term_group = factor(term_group, levels = c("Parity", "Year", "Season")),
    direction = dplyr::case_when(
      CI_low > 1 ~ "Higher odds",
      CI_high < 1 ~ "Lower odds",
      TRUE ~ "Uncertain"
    ),
    direction = factor(direction, levels = c("Lower odds", "Uncertain", "Higher odds")),
    evidence = pmin(-log10(p_value), 20)
  )

forest_order <- forest |>
  dplyr::group_by(term_clean) |>
  dplyr::summarise(order_value = max(abs(log(OR)), na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(order_value)

forest <- forest |>
  dplyr::mutate(term_clean = factor(term_clean, levels = forest_order$term_clean))

write_csv_safe(forest, file.path(CANDIDATE_SOURCE_DIR, "advanced_C_adjusted_model_forest.csv"))

pC <- ggplot2::ggplot(forest, ggplot2::aes(x = OR, y = term_clean)) +
  ggplot2::annotate("rect", xmin = 0.2, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#F3F7FB", alpha = 0.75) +
  ggplot2::annotate("rect", xmin = 1, xmax = 8, ymin = -Inf, ymax = Inf, fill = "#FBF3EE", alpha = 0.75) +
  ggplot2::geom_vline(xintercept = 1, linetype = "22", linewidth = 0.32, colour = "#5E6670") +
  ggplot2::geom_segment(
    ggplot2::aes(x = CI_low, xend = CI_high, yend = term_clean, colour = direction),
    linewidth = 0.55,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = direction, size = evidence),
    shape = 21,
    colour = "white",
    stroke = 0.25,
    alpha = 0.94
  ) +
  ggplot2::facet_grid(term_group ~ phenotype, scales = "free_y", space = "free_y") +
  ggplot2::scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  ggplot2::scale_colour_manual(values = c("Lower odds" = candidate_palette[["lower"]], "Uncertain" = "#8C939D", "Higher odds" = candidate_palette[["upper"]])) +
  ggplot2::scale_fill_manual(values = c("Lower odds" = candidate_palette[["lower"]], "Uncertain" = "#8C939D", "Higher odds" = candidate_palette[["upper"]])) +
  ggplot2::scale_size_continuous(range = c(1.4, 4.2), breaks = c(2, 5, 10, 20)) +
  ggplot2::labs(
    x = "Adjusted odds ratio relative to reference level (log scale)",
    y = NULL,
    size = expression(-log[10](P)),
    title = "Adjusted evidence map for birth-weight tail phenotypes",
    subtitle = "GLMM-derived pairwise contrasts; color encodes direction and point size encodes evidence strength"
  ) +
  theme_candidate(base_size = 6.8) +
  ggplot2::theme(
    legend.box = "horizontal",
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.18, colour = "#ECEFF2")
  )

save_candidate(pC, "Advanced_FigureC_model_evidence_forest", width_mm = 183, height_mm = 145)

# Figure D: distribution atlas -------------------------------------------------

make_density_ridges <- function(data, group_var, group_levels, group_title) {
  pieces <- list()
  idx <- 1
  for (sex_value in levels(data$calf_sex)) {
    for (group_value in group_levels) {
      dat <- data |>
        dplyr::filter(calf_sex == sex_value, as.character(.data[[group_var]]) == group_value, !is.na(birth_weight))
      if (nrow(dat) < 20) next
      dens <- stats::density(dat$birth_weight, from = 20, to = 70, n = 260, adjust = 1.05, na.rm = TRUE)
      group_index <- match(group_value, group_levels)
      pieces[[idx]] <- data.frame(
        x = dens$x,
        density = dens$y,
        density_scaled = dens$y / max(dens$y, na.rm = TRUE) * 0.72,
        group_index = group_index,
        group_label = group_value,
        calf_sex = sex_value,
        n = nrow(dat),
        median_birth_weight = stats::median(dat$birth_weight, na.rm = TRUE),
        group_title = group_title,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  dplyr::bind_rows(pieces)
}

ridge_parity <- make_density_ridges(singleton_tail, "parity_group", c("1", "2", "3", "4", "5", "6", ">=7"), "Parity group")
ridge_season <- make_density_ridges(singleton_tail, "calf_birth_season", c("Spring", "Summer", "Autumn", "Winter"), "Calving season")

ridge_medians <- dplyr::bind_rows(ridge_parity, ridge_season) |>
  dplyr::distinct(group_title, calf_sex, group_label, group_index, n, median_birth_weight)

write_csv_safe(ridge_parity, file.path(CANDIDATE_SOURCE_DIR, "advanced_D_ridge_parity_source.csv"))
write_csv_safe(ridge_season, file.path(CANDIDATE_SOURCE_DIR, "advanced_D_ridge_season_source.csv"))
write_csv_safe(ridge_medians, file.path(CANDIDATE_SOURCE_DIR, "advanced_D_ridge_medians.csv"))

plot_ridges <- function(ridge_data, median_data, y_labels, title_text) {
  ggplot2::ggplot(ridge_data, ggplot2::aes(x = x, group = interaction(calf_sex, group_label))) +
    ggplot2::geom_segment(
      data = ridge_data |>
        dplyr::distinct(calf_sex, group_label, group_index),
      ggplot2::aes(x = 20, xend = 70, y = group_index, yend = group_index),
      inherit.aes = FALSE,
      linewidth = 0.18,
      colour = "#D9DEE5"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = group_index, ymax = group_index + density_scaled, fill = calf_sex),
      alpha = 0.62,
      colour = "white",
      linewidth = 0.20
    ) +
    ggplot2::geom_point(
      data = median_data,
      ggplot2::aes(x = median_birth_weight, y = group_index + 0.08, fill = calf_sex),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.7,
      colour = "white",
      stroke = 0.28
    ) +
    ggplot2::facet_grid(. ~ calf_sex) +
    ggplot2::scale_fill_manual(values = sex_fill) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(y_labels),
      labels = y_labels,
      expand = ggplot2::expansion(mult = c(0.02, 0.15))
    ) +
    ggplot2::coord_cartesian(xlim = c(22, 66), clip = "off") +
    ggplot2::labs(x = "Birth weight (kg)", y = NULL, title = title_text) +
    theme_candidate(base_size = 6.8) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_line(linewidth = 0.18, colour = "#ECEFF2")
    )
}

pD_parity <- plot_ridges(
  ridge_parity,
  ridge_medians |> dplyr::filter(group_title == "Parity group"),
  c("1", "2", "3", "4", "5", "6", ">=7"),
  "Parity distribution atlas"
)

pD_season <- plot_ridges(
  ridge_season,
  ridge_medians |> dplyr::filter(group_title == "Calving season"),
  c("Spring", "Summer", "Autumn", "Winter"),
  "Seasonal distribution atlas"
)

advanced_D <- pD_parity / pD_season +
  patchwork::plot_layout(heights = c(1.2, 0.9)) +
  patchwork::plot_annotation(
    title = "Birth-weight distribution atlas",
    subtitle = "Density ridges show the full phenotype distribution; white-centered points mark group medians",
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 8.8, face = "bold", colour = candidate_palette[["neutral_dark"]]),
    plot.subtitle = ggplot2::element_text(size = 6.8, colour = candidate_palette[["neutral_mid"]]),
    plot.tag = ggplot2::element_text(size = 8, face = "bold")
  )

save_candidate(advanced_D, "Advanced_FigureD_distribution_atlas", width_mm = 183, height_mm = 142)

write_log("15_make_advanced_candidate_figures.R completed")
