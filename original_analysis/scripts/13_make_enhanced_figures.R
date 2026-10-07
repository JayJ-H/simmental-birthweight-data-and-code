source("scripts/00_config.R")

write_log("13_make_enhanced_figures.R started")

ENHANCED_FIGURE_DIR <- file.path(FIGURE_DIR, "enhanced")
ENHANCED_SOURCE_DIR <- file.path(ENHANCED_FIGURE_DIR, "source_data")
create_dir_if_missing(ENHANCED_FIGURE_DIR)
create_dir_if_missing(ENHANCED_SOURCE_DIR)

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |>
  prepare_analysis_types()
sex_var <- sex_variable_name()

singleton_tail <- singleton_tail |>
  dplyr::mutate(
    parity_group = factor(as.character(parity_group), levels = c("1", "2", "3", "4", "5", "6", ">=7"), ordered = FALSE),
    dam_age_group = factor(as.character(dam_age_group), levels = c("<=30", "31-42", "43-60", ">60"), ordered = FALSE),
    calf_birth_season = factor(as.character(calf_birth_season), levels = c("Spring", "Summer", "Autumn", "Winter"), ordered = FALSE),
    calf_sex = factor(as.character(calf_sex), levels = c("female", "male"), ordered = FALSE)
  )

theme_enhanced <- function(base_size = 7.4, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#242424"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = "#242424"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.7, colour = "#222222"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 0.8),
      strip.background = ggplot2::element_rect(fill = "#F3F5F7", colour = NA),
      strip.text = ggplot2::element_text(size = base_size - 0.2, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 1.2, face = "bold", colour = "#1F2430"),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.2, colour = "#555B66"),
      plot.caption = ggplot2::element_text(size = base_size - 1.2, colour = "#696F7A", hjust = 0),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.18, colour = "#E9EDF2"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 7, 5, 5)
    )
}

sex_palette <- c(male = "#2F6DAE", female = "#D78642")
tail_palette <- c(lower_tail = "#3C77AF", upper_tail = "#C96D3A")
heat_palette <- c("#F2F5F7", "#B8D7E8", "#5E9BC6", "#244F79")
rate_palette <- c("#F7F3EE", "#E8C5A2", "#C9783E", "#7E3F24")

save_enhanced <- function(plot, filename, width_mm = 183, height_mm = 125, dpi = 600) {
  path <- file.path(ENHANCED_FIGURE_DIR, filename)
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggplot2::ggsave(paste0(path, ".pdf"), plot, width = w, height = h, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(path, ".png"), plot, width = w, height = h, dpi = dpi)
  svglite::svglite(paste0(path, ".svg"), width = w, height = h)
  print(plot)
  grDevices::dev.off()
  invisible(path)
}

clean_phenotype <- function(x) {
  dplyr::recode(x, lower_tail = "Lower-tail phenotype", upper_tail = "Upper-tail phenotype", .default = x)
}

thresholds <- openxlsx::read.xlsx(file.path(TABLE_DIR, "sex_specific_thresholds.xlsx"), sheet = "sex_specific_thresholds")
names(thresholds)[1] <- sex_var
threshold_long <- thresholds |>
  tidyr::pivot_longer(c(P10, P90), names_to = "threshold", values_to = "birth_weight") |>
  dplyr::mutate(label = paste0(threshold, "=", birth_weight, " kg"))

tail_rate <- singleton_tail |>
  dplyr::group_by(.data[[sex_var]]) |>
  dplyr::summarise(
    lower_tail = mean(lower_tail, na.rm = TRUE),
    upper_tail = mean(upper_tail, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(c(lower_tail, upper_tail), names_to = "phenotype", values_to = "rate") |>
  dplyr::mutate(phenotype_label = clean_phenotype(phenotype))
names(tail_rate)[1] <- sex_var

write_csv_safe(threshold_long, file.path(ENHANCED_SOURCE_DIR, "figure1_thresholds.csv"))
write_csv_safe(tail_rate, file.path(ENHANCED_SOURCE_DIR, "figure1_tail_rate.csv"))

p1_density <- ggplot2::ggplot(singleton_tail, ggplot2::aes(x = birth_weight, fill = .data[[sex_var]], colour = .data[[sex_var]])) +
  ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 44, position = "identity", alpha = 0.14, linewidth = 0) +
  ggplot2::geom_density(linewidth = 0.62, alpha = 0.04) +
  ggplot2::geom_vline(data = threshold_long, ggplot2::aes(xintercept = birth_weight, colour = .data[[sex_var]], linetype = threshold), linewidth = 0.38, show.legend = FALSE) +
  ggplot2::geom_text(
    data = threshold_long,
    ggplot2::aes(x = birth_weight, y = 0.061, label = label, colour = .data[[sex_var]]),
    inherit.aes = FALSE,
    angle = 90,
    vjust = -0.35,
    hjust = 1,
    size = 2.1,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = sex_palette) +
  ggplot2::scale_colour_manual(values = sex_palette) +
  ggplot2::scale_linetype_manual(values = c(P10 = "22", P90 = "22")) +
  ggplot2::coord_cartesian(xlim = c(20, 70), ylim = c(0, 0.064), clip = "off") +
  ggplot2::labs(
    x = "Birth weight (kg)",
    y = "Density",
    title = "Sex-specific birth-weight distributions",
    subtitle = "Dashed lines mark within-sex P10 and P90 thresholds"
  ) +
  theme_enhanced()

p1_violin <- ggplot2::ggplot(singleton_tail, ggplot2::aes(x = parity_group, y = birth_weight, fill = .data[[sex_var]])) +
  ggplot2::geom_violin(position = ggplot2::position_dodge(width = 0.76), width = 0.72, alpha = 0.50, trim = TRUE, colour = NA) +
  ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.76), width = 0.13, outlier.shape = NA, alpha = 0.78, linewidth = 0.25) +
  ggplot2::stat_summary(
    ggplot2::aes(group = .data[[sex_var]]),
    fun = median,
    geom = "point",
    position = ggplot2::position_dodge(width = 0.76),
    size = 0.85,
    colour = "#151515",
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = sex_palette) +
  ggplot2::coord_cartesian(ylim = c(25, 65)) +
  ggplot2::labs(x = "Parity group", y = "Birth weight (kg)", title = "Distribution by parity and sex") +
  theme_enhanced() +
  ggplot2::theme(legend.position = "none")

p1_tail <- ggplot2::ggplot(tail_rate, ggplot2::aes(x = .data[[sex_var]], y = rate, fill = phenotype)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.62), width = 0.50, alpha = 0.95) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(rate, accuracy = 0.1)),
    position = ggplot2::position_dodge(width = 0.62),
    vjust = -0.45,
    size = 2.25
  ) +
  ggplot2::scale_fill_manual(values = tail_palette, labels = clean_phenotype) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.13), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Event rate", title = "Tail phenotype frequency") +
  theme_enhanced()

figure1 <- p1_density / (p1_violin | p1_tail) +
  patchwork::plot_layout(heights = c(1.15, 1), widths = c(1.65, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 9, face = "bold"))
save_enhanced(figure1, "Enhanced_Figure1_distribution_tail", width_mm = 183, height_mm = 155)

monthly_summary <- singleton_tail |>
  dplyr::mutate(calving_month = lubridate::floor_date(calving_date, "month")) |>
  dplyr::group_by(calving_month, .data[[sex_var]]) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    SD = stats::sd(birth_weight, na.rm = TRUE),
    SE = SD / sqrt(n),
    lower_tail_rate = mean(lower_tail, na.rm = TRUE),
    upper_tail_rate = mean(upper_tail, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    low = mean_birth_weight - 1.96 * SE,
    high = mean_birth_weight + 1.96 * SE
  )
names(monthly_summary)[2] <- sex_var
write_csv_safe(monthly_summary, file.path(ENHANCED_SOURCE_DIR, "figure2_monthly_summary.csv"))

p2_line <- ggplot2::ggplot(monthly_summary, ggplot2::aes(x = calving_month, y = mean_birth_weight, colour = .data[[sex_var]], fill = .data[[sex_var]])) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high), alpha = 0.14, colour = NA) +
  ggplot2::geom_line(linewidth = 0.62) +
  ggplot2::geom_point(size = 1.25, stroke = 0.2) +
  ggplot2::scale_colour_manual(values = sex_palette) +
  ggplot2::scale_fill_manual(values = sex_palette) +
  ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  ggplot2::labs(x = NULL, y = "Mean birth weight (kg)", title = "Monthly birth-weight trajectory", subtitle = "Ribbons show approximate 95% CI for the monthly mean") +
  theme_enhanced() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

p2_count <- ggplot2::ggplot(monthly_summary, ggplot2::aes(x = calving_month, y = n, fill = .data[[sex_var]])) +
  ggplot2::geom_col(width = 24, alpha = 0.82, position = "stack") +
  ggplot2::scale_fill_manual(values = sex_palette) +
  ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  ggplot2::labs(x = NULL, y = "Singleton records", title = "Monthly record volume") +
  theme_enhanced() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")

tail_month <- monthly_summary |>
  dplyr::select(calving_month, dplyr::all_of(sex_var), lower_tail_rate, upper_tail_rate) |>
  tidyr::pivot_longer(c(lower_tail_rate, upper_tail_rate), names_to = "phenotype", values_to = "rate") |>
  dplyr::mutate(phenotype = dplyr::recode(phenotype, lower_tail_rate = "Lower-tail", upper_tail_rate = "Upper-tail"))

p2_tail <- ggplot2::ggplot(tail_month, ggplot2::aes(x = calving_month, y = rate, colour = phenotype)) +
  ggplot2::geom_line(linewidth = 0.52) +
  ggplot2::facet_wrap(stats::as.formula(paste("~", sex_var)), nrow = 1) +
  ggplot2::scale_colour_manual(values = c(`Lower-tail` = tail_palette[[1]], `Upper-tail` = tail_palette[[2]])) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ggplot2::labs(x = NULL, y = "Monthly event rate", title = "Monthly tail phenotype frequency") +
  theme_enhanced()

figure2 <- p2_line / (p2_count | p2_tail) +
  patchwork::plot_layout(heights = c(1.05, 1), widths = c(1, 1.15)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 9, face = "bold"))
save_enhanced(figure2, "Enhanced_Figure2_temporal_patterns", width_mm = 183, height_mm = 145)

year_season <- singleton_tail |>
  dplyr::group_by(calf_birth_year, calf_birth_season) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    lower_tail_rate = mean(lower_tail, na.rm = TRUE),
    upper_tail_rate = mean(upper_tail, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    year = factor(calf_birth_year),
    season = factor(calf_birth_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
  )
write_csv_safe(year_season, file.path(ENHANCED_SOURCE_DIR, "figure3_year_season.csv"))

heatmap_panel <- function(data, value, label, fill_label, palette, percent = FALSE) {
  text_label <- if (percent) scales::percent(data[[value]], accuracy = 0.1) else sprintf("%.1f", data[[value]])
  ggplot2::ggplot(data, ggplot2::aes(x = season, y = year, fill = .data[[value]])) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
    ggplot2::geom_text(ggplot2::aes(label = text_label), size = 2.1, colour = "#202020") +
    ggplot2::scale_fill_gradientn(colours = palette, name = fill_label) +
    ggplot2::labs(x = NULL, y = NULL, title = label) +
    theme_enhanced(base_size = 6.8) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

p3_mean <- heatmap_panel(year_season, "mean_birth_weight", "Mean birth weight", "kg", heat_palette, FALSE)
p3_lower <- heatmap_panel(year_season, "lower_tail_rate", "Lower-tail frequency", "rate", rate_palette, TRUE)
p3_upper <- heatmap_panel(year_season, "upper_tail_rate", "Upper-tail frequency", "rate", rate_palette, TRUE)
figure3 <- p3_mean | p3_lower | p3_upper
figure3 <- figure3 +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Year-season structure of singleton birth weight",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold", colour = "#1F2430"))
  ) &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 9, face = "bold"))
save_enhanced(figure3, "Enhanced_Figure3_year_season_heatmaps", width_mm = 183, height_mm = 90)

build_emmeans_or <- function(model_path, phenotype_label) {
  if (!file.exists(model_path)) return(tibble::tibble())
  model <- readRDS(model_path)
  terms_to_plot <- c("parity_group", "calf_birth_year", "calf_birth_season")
  purrr::map_dfr(terms_to_plot, function(term_name) {
    emm <- tryCatch(emmeans::emmeans(model, stats::as.formula(paste("~", term_name))), error = function(e) e)
    if (inherits(emm, "error")) return(tibble::tibble())
    emm_grid <- as.data.frame(emm)
    ref_value <- switch(term_name,
                        parity_group = "1",
                        calf_birth_year = "calf_birth_year2023",
                        calf_birth_season = "Spring",
                        NA_character_)
    ref_index <- match(ref_value, as.character(emm_grid[[term_name]]))
    if (is.na(ref_index) && term_name == "calf_birth_year") {
      ref_index <- match("2023", as.character(emm_grid[[term_name]]))
    }
    if (is.na(ref_index)) ref_index <- 1
    ctr <- tryCatch(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_index, adjust = "none"), error = function(e) e)
    if (inherits(ctr, "error")) return(tibble::tibble())
    out <- as.data.frame(summary(ctr, infer = c(TRUE, TRUE), type = "link"))
    low_col <- intersect(c("asymp.LCL", "lower.CL"), names(out))[1]
    high_col <- intersect(c("asymp.UCL", "upper.CL"), names(out))[1]
    contrast_clean <- gsub("calf_birth_year", "", out$contrast)
    contrast_clean <- gsub(" - ", " vs ", contrast_clean, fixed = TRUE)
    term_clean <- if (term_name == "parity_group") {
      paste("Parity", contrast_clean)
    } else if (term_name == "calf_birth_year") {
      paste("Year", contrast_clean)
    } else if (term_name == "calf_birth_season") {
      contrast_clean
    } else {
      contrast_clean
    }
    out |>
      dplyr::transmute(
        phenotype = phenotype_label,
        term_group = dplyr::recode(term_name, parity_group = "Parity", calf_birth_year = "Year", calf_birth_season = "Season"),
        term_clean = term_clean,
        OR = exp(estimate),
        CI_low = exp(.data[[low_col]]),
        CI_high = exp(.data[[high_col]]),
        p_value = p.value
      )
  })
}

forest <- dplyr::bind_rows(
  build_emmeans_or(file.path(MODEL_DIR, "glmm_lower_tail_main.rds"), "Lower-tail phenotype"),
  build_emmeans_or(file.path(MODEL_DIR, "glmm_upper_tail_main.rds"), "Upper-tail phenotype")
) |>
  dplyr::filter(is.finite(OR), is.finite(CI_low), is.finite(CI_high))

if (nrow(forest) > 0) {
  write_csv_safe(forest, file.path(ENHANCED_SOURCE_DIR, "figure4_glmm_emmeans_or.csv"))

  p4 <- ggplot2::ggplot(forest, ggplot2::aes(y = stats::reorder(term_clean, OR), x = OR, colour = phenotype)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.28, colour = "#7A8088") +
    ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = term_clean), linewidth = 0.48, position = ggplot2::position_dodge(width = 0.58)) +
    ggplot2::geom_point(size = 1.45, position = ggplot2::position_dodge(width = 0.58)) +
    ggplot2::facet_grid(term_group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8")) +
    ggplot2::scale_colour_manual(values = c(`Lower-tail phenotype` = tail_palette[[1]], `Upper-tail phenotype` = tail_palette[[2]])) +
    ggplot2::labs(
      x = "Odds ratio relative to reference level (log scale)",
      y = NULL,
      title = "Adjusted contrasts for distribution-defined tail phenotypes",
      subtitle = "Contrasts are derived from GLMM estimated marginal means; reference levels are parity 1, 2023, and Spring"
    ) +
    theme_enhanced(base_size = 7.0)
  save_enhanced(p4, "Enhanced_Figure4_glmm_forest", width_mm = 170, height_mm = 125)
}

read_optional_sheet <- function(path, sheet) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(openxlsx::read.xlsx(path, sheet = sheet), error = function(e) data.frame())
}

service_summary <- read_optional_sheet(file.path(TABLE_DIR, "service_sire_exploratory.xlsx"), "n_ge_50_raw_summary")
mgs_summary <- read_optional_sheet(file.path(TABLE_DIR, "maternal_grandsire_exploratory.xlsx"), "n_ge_30_raw_summary")

rank_plot_data <- dplyr::bind_rows(
  service_summary |>
    dplyr::filter(!is.na(service_sire), n >= 50) |>
    dplyr::mutate(label_type = "Service sire", label = as.character(service_sire)),
  mgs_summary |>
    dplyr::filter(!is.na(maternal_grandsire), n >= 50) |>
    dplyr::mutate(label_type = "Maternal-grandsire label", label = as.character(maternal_grandsire))
) |>
  dplyr::group_by(label_type) |>
  dplyr::arrange(raw_mean_birth_weight, .by_group = TRUE) |>
  dplyr::mutate(rank = dplyr::row_number(), keep = rank <= 8 | rank > dplyr::n() - 8) |>
  dplyr::filter(keep) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    label_plot = factor(label, levels = unique(label[order(label_type, raw_mean_birth_weight)])),
    low = raw_mean_birth_weight - 1.96 * SE,
    high = raw_mean_birth_weight + 1.96 * SE
  )

if (nrow(rank_plot_data) > 0) {
  write_csv_safe(rank_plot_data, file.path(ENHANCED_SOURCE_DIR, "figure5_label_rankings.csv"))
  p5 <- ggplot2::ggplot(rank_plot_data, ggplot2::aes(x = raw_mean_birth_weight, y = label_plot)) +
    ggplot2::geom_segment(ggplot2::aes(x = low, xend = high, yend = label_plot), colour = "#8C929A", linewidth = 0.38) +
    ggplot2::geom_point(ggplot2::aes(size = n, colour = label_type), alpha = 0.90) +
    ggplot2::facet_grid(label_type ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_manual(values = c(`Service sire` = "#2F6DAE", `Maternal-grandsire label` = "#6B8F3A")) +
    ggplot2::scale_size_continuous(range = c(1.2, 3.6), name = "Records") +
    ggplot2::labs(
      x = "Raw mean birth weight (kg)",
      y = NULL,
      title = "Exploratory label-associated phenotypic differences",
      subtitle = "Lowest and highest label-level raw means among labels with at least 50 records"
    ) +
    theme_enhanced(base_size = 6.8)
  save_enhanced(p5, "Enhanced_Figure5_label_rankings", width_mm = 170, height_mm = 150)
}

age_parity <- singleton_tail |>
  dplyr::filter(!is.na(dam_age_group), !is.na(parity_group)) |>
  dplyr::group_by(parity_group, dam_age_group) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    lower_tail_rate = mean(lower_tail, na.rm = TRUE),
    upper_tail_rate = mean(upper_tail, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_safe(age_parity, file.path(ENHANCED_SOURCE_DIR, "figure6_age_parity_summary.csv"))

p6_heat <- ggplot2::ggplot(age_parity, ggplot2::aes(x = dam_age_group, y = parity_group, fill = mean_birth_weight)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(sprintf("%.1f", mean_birth_weight), "\n", "n=", n)), size = 2.0, lineheight = 0.92) +
  ggplot2::scale_fill_gradientn(colours = heat_palette, name = "Mean kg") +
  ggplot2::labs(x = "Dam age group (months)", y = "Parity group", title = "Parity-age structure of birth weight") +
  theme_enhanced(base_size = 6.9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")

scatter_data <- singleton_tail |>
  dplyr::filter(!is.na(dam_age_month), !is.na(birth_weight)) |>
  dplyr::arrange(row_id) |>
  dplyr::slice(unique(round(seq(1, dplyr::n(), length.out = min(1800, dplyr::n())))))

p6_scatter <- ggplot2::ggplot(scatter_data, ggplot2::aes(x = dam_age_month, y = birth_weight, colour = .data[[sex_var]])) +
  ggplot2::geom_point(alpha = 0.18, size = 0.55) +
  ggplot2::geom_smooth(method = "loess", se = TRUE, linewidth = 0.55, span = 0.45) +
  ggplot2::scale_colour_manual(values = sex_palette) +
  ggplot2::coord_cartesian(xlim = c(18, 125), ylim = c(25, 70)) +
  ggplot2::labs(x = "Dam age at calving (months)", y = "Birth weight (kg)", title = "Continuous dam-age pattern") +
  theme_enhanced(base_size = 6.9)

figure6 <- p6_heat | p6_scatter
figure6 <- figure6 + patchwork::plot_layout(widths = c(1, 1.15)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Dam age and parity are strongly coupled in the field records",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold", colour = "#1F2430"))
  ) &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 9, face = "bold"))
save_enhanced(figure6, "Enhanced_Figure6_age_parity_structure", width_mm = 183, height_mm = 95)

write_log("13_make_enhanced_figures.R completed")
