source("scripts/00_config.R")

write_log("16_make_reference_style_triptych.R started")

TRIPTYCH_DIR <- file.path(FIGURE_DIR, "advanced_candidates", "reference_style")
TRIPTYCH_SOURCE_DIR <- file.path(TRIPTYCH_DIR, "source_data")
create_dir_if_missing(TRIPTYCH_DIR)
create_dir_if_missing(TRIPTYCH_SOURCE_DIR)

pal <- c(
  dark = "#252A32",
  mid = "#6F7782",
  grid = "#ECEFF2",
  female = "#D48642",
  male = "#2F6DAE",
  lower = "#3B6EA5",
  upper = "#C56A3A",
  service = "#8B6FA9",
  mgs = "#5F9D8B"
)

theme_triptych <- function(base_size = 6.7, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.30, colour = pal[["dark"]]),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = pal[["dark"]]),
      axis.title = ggplot2::element_text(size = base_size, colour = pal[["dark"]]),
      axis.text = ggplot2::element_text(size = base_size - 0.7, colour = "#33373D"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 0.8),
      strip.background = ggplot2::element_rect(fill = "#F1F3F5", colour = NA),
      strip.text = ggplot2::element_text(size = base_size - 0.2, face = "bold", colour = pal[["dark"]]),
      plot.title = ggplot2::element_text(size = base_size + 0.4, face = "bold", colour = pal[["dark"]]),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.9, colour = pal[["mid"]]),
      plot.caption = ggplot2::element_text(size = base_size - 1.2, colour = pal[["mid"]], hjust = 0),
      panel.grid.major = ggplot2::element_line(linewidth = 0.16, colour = pal[["grid"]]),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}

save_triptych <- function(plot, filename, width_mm = 183, height_mm = 72, dpi = 500) {
  path <- file.path(TRIPTYCH_DIR, filename)
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

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |>
  prepare_analysis_types() |>
  dplyr::mutate(
    calf_sex = factor(as.character(calf_sex), levels = c("female", "male")),
    calf_birth_year = factor(as.character(calf_birth_year), levels = sort(unique(as.character(calf_birth_year)))),
    calf_birth_season = factor(as.character(calf_birth_season), levels = c("Spring", "Summer", "Autumn", "Winter"))
  )

# Panel A: dam-age response curve --------------------------------------------

age_data <- singleton_tail |>
  dplyr::filter(!is.na(dam_age_month), !is.na(birth_weight), !is.na(calf_sex))

age_limits <- as.numeric(stats::quantile(age_data$dam_age_month, c(0.03, 0.97), na.rm = TRUE))
age_model_data <- age_data |>
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
age_grid$se <- as.numeric(age_pred$se.fit)
age_grid$lower <- age_grid$pred - 1.96 * age_grid$se
age_grid$upper <- age_grid$pred + 1.96 * age_grid$se

age_bins <- age_model_data |>
  dplyr::mutate(age_bin = cut(dam_age_month, breaks = seq(floor(age_limits[[1]]), ceiling(age_limits[[2]]), length.out = 12), include.lowest = TRUE)) |>
  dplyr::group_by(calf_sex, age_bin) |>
  dplyr::summarise(
    n = dplyr::n(),
    age_mid = mean(dam_age_month, na.rm = TRUE),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    se = stats::sd(birth_weight, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) |>
  dplyr::filter(n >= 25)

write_csv_safe(age_grid, file.path(TRIPTYCH_SOURCE_DIR, "panel_A_dam_age_response_predictions.csv"))
write_csv_safe(age_bins, file.path(TRIPTYCH_SOURCE_DIR, "panel_A_dam_age_response_bins.csv"))

pA <- ggplot2::ggplot(age_grid, ggplot2::aes(x = dam_age_month, y = pred, colour = calf_sex, fill = calf_sex)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.13, linewidth = 0, colour = NA, show.legend = FALSE) +
  ggplot2::geom_line(linewidth = 0.72, show.legend = FALSE) +
  ggplot2::geom_point(
    data = age_bins,
    ggplot2::aes(x = age_mid, y = mean_birth_weight, colour = calf_sex),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    stroke = 0.35,
    size = 1.15,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = age_grid |>
      dplyr::group_by(calf_sex) |>
      dplyr::slice_max(dam_age_month, n = 1) |>
      dplyr::ungroup(),
    ggplot2::aes(label = dplyr::recode(as.character(calf_sex), female = "Female", male = "Male")),
    hjust = -0.05,
    size = 2.05,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = c(female = pal[["female"]], male = pal[["male"]]), labels = c(female = "Female", male = "Male")) +
  ggplot2::scale_fill_manual(values = c(female = pal[["female"]], male = pal[["male"]]), labels = c(female = "Female", male = "Male")) +
  ggplot2::coord_cartesian(xlim = age_limits + c(0, 8), ylim = c(40, 50), clip = "off") +
  ggplot2::labs(
    x = "Dam age at calving (months)",
    y = "Predicted birth weight (kg)",
    title = "A  Dam-age response"
  ) +
  theme_triptych(base_size = 6.2) +
  ggplot2::theme(legend.position = "none")

# Panel B: exploratory label volcano -----------------------------------------

make_label_volcano <- function(data, var, label_type, min_n = 30) {
  labels <- data |>
    dplyr::filter(!is.na(.data[[var]]), .data[[var]] != "", !is.na(birth_weight)) |>
    dplyr::count(label = .data[[var]], name = "n") |>
    dplyr::filter(n >= min_n) |>
    dplyr::pull(label)

  rows <- lapply(labels, function(label_value) {
    in_group <- data |>
      dplyr::filter(.data[[var]] == label_value, !is.na(birth_weight)) |>
      dplyr::pull(birth_weight)
    out_group <- data |>
      dplyr::filter(.data[[var]] != label_value | is.na(.data[[var]]), !is.na(birth_weight)) |>
      dplyr::pull(birth_weight)
    test <- stats::t.test(in_group, out_group)
    data.frame(
      label = as.character(label_value),
      label_type = label_type,
      n = length(in_group),
      mean_birth_weight = mean(in_group, na.rm = TRUE),
      delta_kg = mean(in_group, na.rm = TRUE) - mean(out_group, na.rm = TRUE),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

label_volcano <- dplyr::bind_rows(
  make_label_volcano(singleton_tail, "service_sire", "Service sire", min_n = 35),
  make_label_volcano(singleton_tail, "maternal_grandsire", "Maternal grandsire", min_n = 30)
) |>
  dplyr::mutate(
    q_value = stats::p.adjust(p_value, method = "BH"),
    neg_log10_q = pmin(-log10(q_value), 8),
    direction = dplyr::case_when(
      q_value < 0.05 & delta_kg >= 1 ~ "Higher BW",
      q_value < 0.05 & delta_kg <= -1 ~ "Lower BW",
      TRUE ~ "Not highlighted"
    ),
    direction = factor(direction, levels = c("Lower BW", "Not highlighted", "Higher BW")),
    label_type = factor(label_type, levels = c("Service sire", "Maternal grandsire"))
  )

volcano_labels <- label_volcano |>
  dplyr::arrange(q_value, dplyr::desc(abs(delta_kg))) |>
  dplyr::slice_head(n = 4)

write_csv_safe(label_volcano, file.path(TRIPTYCH_SOURCE_DIR, "panel_B_exploratory_label_volcano.csv"))

pB <- ggplot2::ggplot(label_volcano, ggplot2::aes(x = delta_kg, y = neg_log10_q)) +
  ggplot2::annotate("rect", xmin = -Inf, xmax = -1, ymin = -log10(0.05), ymax = Inf, fill = "#EEF5FA", alpha = 0.85) +
  ggplot2::annotate("rect", xmin = 1, xmax = Inf, ymin = -log10(0.05), ymax = Inf, fill = "#FBF0EA", alpha = 0.85) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.25, colour = "#7F8790") +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "22", linewidth = 0.25, colour = "#7F8790") +
  ggplot2::geom_point(
    ggplot2::aes(fill = direction, shape = label_type, size = n),
    alpha = 0.78,
    colour = "white",
    stroke = 0.25,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = volcano_labels,
    ggplot2::aes(label = label),
    size = 1.65,
    box.padding = 0.18,
    point.padding = 0.08,
    min.segment.length = 0,
    segment.size = 0.16,
    segment.colour = "#9AA2AC",
    force = 2.5,
    force_pull = 0.3,
    max.time = 2,
    max.overlaps = Inf,
    colour = pal[["dark"]],
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = c("Lower BW" = pal[["lower"]], "Not highlighted" = "#B7BDC6", "Higher BW" = pal[["upper"]])) +
  ggplot2::scale_shape_manual(values = c("Service sire" = 21, "Maternal grandsire" = 24)) +
  ggplot2::scale_size_continuous(range = c(1.0, 3.6), breaks = c(50, 100, 200, 400)) +
  ggplot2::coord_cartesian(xlim = c(-4.5, 4.5), ylim = c(0, 8.2), clip = "on") +
  ggplot2::labs(
    x = "Mean difference (kg)",
    y = expression(-log[10]("BH-adjusted P")),
    title = "B  Label screen"
  ) +
  theme_triptych(base_size = 6.2) +
  ggplot2::theme(legend.position = "none")

# Panel C: age distribution plus year-season scatter --------------------------

stratum <- singleton_tail |>
  dplyr::filter(!is.na(dam_age_month), !is.na(birth_weight), !is.na(calf_birth_year), !is.na(calf_birth_season)) |>
  dplyr::group_by(calf_birth_year, calf_birth_season) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_age = mean(dam_age_month, na.rm = TRUE),
    se_age = stats::sd(dam_age_month, na.rm = TRUE) / sqrt(n),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    se_birth_weight = stats::sd(birth_weight, na.rm = TRUE) / sqrt(n),
    lower_tail_rate = mean(lower_tail, na.rm = TRUE),
    upper_tail_rate = mean(upper_tail, na.rm = TRUE),
    tail_burden = lower_tail_rate + upper_tail_rate,
    .groups = "drop"
  ) |>
  dplyr::filter(n >= 30) |>
  dplyr::mutate(
    label = paste(calf_birth_year, substr(as.character(calf_birth_season), 1, 3), sep = "-"),
    ci_age = 1.96 * se_age,
    ci_birth_weight = 1.96 * se_birth_weight
  )

stratum_labels <- stratum |>
  dplyr::slice_max(tail_burden, n = 2, with_ties = FALSE)

write_csv_safe(stratum, file.path(TRIPTYCH_SOURCE_DIR, "panel_C_year_season_age_weight_scatter.csv"))

pC_top <- ggplot2::ggplot(age_data, ggplot2::aes(x = dam_age_month)) +
  ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 34, fill = "#E88478", colour = "white", linewidth = 0.12, alpha = 0.72) +
  ggplot2::geom_density(linewidth = 0.55, colour = pal[["dark"]], alpha = 0.85) +
  ggplot2::coord_cartesian(xlim = c(20, 125), clip = "off") +
  ggplot2::labs(x = NULL, y = "Density", title = "C  Age-structured strata") +
  theme_triptych(base_size = 6.3) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(2, 5, 0, 4)
  )

pC_bottom <- ggplot2::ggplot(stratum, ggplot2::aes(x = mean_age, y = mean_birth_weight)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = mean_age - ci_age, xend = mean_age + ci_age, yend = mean_birth_weight),
    linewidth = 0.32,
    colour = "#A9B1BA",
    alpha = 0.90
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_birth_weight - ci_birth_weight, ymax = mean_birth_weight + ci_birth_weight),
    width = 0,
    linewidth = 0.32,
    colour = "#A9B1BA",
    alpha = 0.90
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = tail_burden, size = n),
    shape = 21,
    colour = "white",
    stroke = 0.32,
    alpha = 0.96
  ) +
  ggrepel::geom_text_repel(
    data = stratum_labels,
    ggplot2::aes(label = label),
    size = 1.65,
    box.padding = 0.15,
    point.padding = 0.08,
    min.segment.length = 0,
    segment.size = 0.15,
    segment.colour = "#9AA2AC",
    force = 2,
    force_pull = 0.3,
    max.time = 2,
    colour = pal[["dark"]],
    max.overlaps = Inf
  ) +
  ggplot2::scale_fill_gradientn(colours = c("#EAF1F4", "#E0B37B", "#C56A3A", "#74312B"), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_size_continuous(range = c(1.4, 4.2), breaks = c(100, 500, 900)) +
  ggplot2::guides(
    size = "none",
    fill = ggplot2::guide_colorbar(
      title = "Tail burden",
      title.position = "top",
      barwidth = grid::unit(2.4, "cm"),
      barheight = grid::unit(0.20, "cm")
    )
  ) +
  ggplot2::coord_cartesian(xlim = c(20, 125), ylim = c(40, 48.5), clip = "off") +
  ggplot2::labs(
    x = "Mean dam age at calving (months)",
    y = "Mean birth weight (kg)",
    fill = "Tail burden",
    size = "Records"
  ) +
  theme_triptych(base_size = 6.3) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(0.45, "cm"),
    legend.title = ggplot2::element_text(size = 5.7, colour = pal[["dark"]]),
    plot.margin = ggplot2::margin(0, 5, 4, 4)
  )

pC <- pC_top / pC_bottom +
  patchwork::plot_layout(heights = c(0.42, 1.0))

triptych <- (pA | pB | pC) +
  patchwork::plot_layout(widths = c(1.0, 1.08, 1.02))

save_triptych(triptych, "Advanced_FigureE_reference_style_triptych", width_mm = 183, height_mm = 104)
save_triptych(triptych, "Advanced_FigureE_reference_style_triptych_v2", width_mm = 183, height_mm = 86)
save_triptych(triptych, "Advanced_FigureE_reference_style_triptych_v3", width_mm = 183, height_mm = 78)

write_log("16_make_reference_style_triptych.R completed")
