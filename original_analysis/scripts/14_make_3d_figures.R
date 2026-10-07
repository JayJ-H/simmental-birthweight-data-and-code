source("scripts/00_config.R")

write_log("14_make_3d_figures.R started")

THREED_DIR <- file.path(FIGURE_DIR, "enhanced", "3d")
THREED_SOURCE_DIR <- file.path(THREED_DIR, "source_data")
create_dir_if_missing(THREED_DIR)
create_dir_if_missing(THREED_SOURCE_DIR)

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |>
  prepare_analysis_types() |>
  dplyr::mutate(
    calf_sex = factor(as.character(calf_sex), levels = c("female", "male")),
    calf_birth_month_num = as.integer(calf_birth_month),
    parity_group = factor(as.character(parity_group), levels = c("1", "2", "3", "4", "5", "6", ">=7"))
  )

sex_palette <- c(female = "#D78642", male = "#2F6DAE")

save_base_plot <- function(filename, plot_fun, width = 8.2, height = 6.2, pointsize = 10) {
  pdf_path <- file.path(THREED_DIR, paste0(filename, ".pdf"))
  png_path <- file.path(THREED_DIR, paste0(filename, ".png"))
  svg_path <- file.path(THREED_DIR, paste0(filename, ".svg"))
  tiff_path <- file.path(THREED_DIR, paste0(filename, ".tiff"))

  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "sans")
  plot_fun()
  grDevices::dev.off()

  ragg::agg_png(png_path, width = width, height = height, units = "in", res = 400, pointsize = pointsize)
  plot_fun()
  grDevices::dev.off()

  svglite::svglite(svg_path, width = width, height = height, pointsize = pointsize)
  plot_fun()
  grDevices::dev.off()

  ragg::agg_tiff(tiff_path, width = width, height = height, units = "in", res = 600, pointsize = pointsize)
  plot_fun()
  grDevices::dev.off()
}

scatter_data <- singleton_tail |>
  dplyr::filter(!is.na(dam_age_month), !is.na(parity), !is.na(birth_weight), !is.na(calf_sex)) |>
  dplyr::arrange(calf_sex, dam_age_month, parity) |>
  dplyr::slice(unique(round(seq(1, dplyr::n(), length.out = min(2200, dplyr::n())))))

write_csv_safe(scatter_data, file.path(THREED_SOURCE_DIR, "3d_scatter_source.csv"))

scatter_expr <- function() {
  par(mar = c(3.4, 3.6, 2.8, 1.0), family = "sans")
  s3d <- scatterplot3d::scatterplot3d(
    x = scatter_data$dam_age_month,
    y = scatter_data$parity,
    z = scatter_data$birth_weight,
    color = sex_palette[as.character(scatter_data$calf_sex)],
    pch = 16,
    cex.symbols = 0.55,
    angle = 46,
    scale.y = 0.95,
    grid = TRUE,
    box = FALSE,
    xlab = "Dam age at calving (months)",
    ylab = "Parity",
    zlab = "Birth weight (kg)",
    main = "3D scatter of birth weight, parity, and dam age"
  )
  legend("topleft", legend = c("female", "male"), col = sex_palette, pch = 16, bty = "n", cex = 0.9)
  mtext("Supplementary-style visualization; downsampled to improve readability", side = 3, line = 0.2, cex = 0.8, adj = 0)
}
save_base_plot("3D_Scatter_birth_weight_parity_age", scatter_expr, width = 8.6, height = 6.4)

surface_data <- singleton_tail |>
  dplyr::filter(!is.na(dam_age_month), !is.na(calf_birth_month_num), !is.na(birth_weight), !is.na(calf_sex))

surface_age_limits <- as.numeric(stats::quantile(surface_data$dam_age_month, probs = c(0.05, 0.95), na.rm = TRUE))
surface_age_limits <- round(surface_age_limits, 1)

fit_surface <- function(sex_value) {
  dat <- surface_data |>
    dplyr::filter(calf_sex == sex_value)
  gam_fit <- mgcv::bam(
    birth_weight ~
      s(dam_age_month, k = 6) +
      s(calf_birth_month_num, bs = "cc", k = 6) +
      ti(dam_age_month, calf_birth_month_num, bs = c("tp", "cc"), k = c(5, 5)),
    data = dat,
    method = "fREML",
    discrete = TRUE,
    knots = list(calf_birth_month_num = c(0.5, 12.5))
  )
  age_seq <- seq(surface_age_limits[[1]], surface_age_limits[[2]], length.out = 55)
  month_seq <- seq(1, 12, length.out = 34)
  grid <- expand.grid(
    dam_age_month = age_seq,
    calf_birth_month_num = month_seq
  )
  pred <- stats::predict(gam_fit, newdata = grid, se.fit = TRUE)
  grid$pred <- as.numeric(pred$fit)
  grid$pred_se <- as.numeric(pred$se.fit)
  grid$age_support <- sprintf("%.1f-%.1f months (5th-95th percentile)", surface_age_limits[[1]], surface_age_limits[[2]])
  list(
    fit = gam_fit,
    z = matrix(grid$pred, nrow = length(age_seq), ncol = length(month_seq)),
    age_seq = age_seq,
    month_seq = month_seq,
    grid = grid
  )
}

female_surface <- fit_surface("female")
male_surface <- fit_surface("male")
surface_zlim <- range(c(female_surface$z, male_surface$z), na.rm = TRUE)
surface_zlim <- c(floor(surface_zlim[[1]]), ceiling(surface_zlim[[2]]))
write_csv_safe(
  dplyr::bind_rows(
    female_surface$grid |>
      dplyr::mutate(calf_sex = "female"),
    male_surface$grid |>
      dplyr::mutate(calf_sex = "male")
  ),
  file.path(THREED_SOURCE_DIR, "3d_surface_source.csv")
)
write_csv_safe(
  data.frame(
    item = c("dam_age_lower", "dam_age_upper", "z_lower", "z_upper", "model_note"),
    value = c(
      surface_age_limits[[1]],
      surface_age_limits[[2]],
      surface_zlim[[1]],
      surface_zlim[[2]],
      "Sex-stratified GAM with cyclic month smooth; predictions restricted to 5th-95th percentile dam-age support."
    )
  ),
  file.path(THREED_SOURCE_DIR, "3d_surface_metadata.csv")
)

surface_colors <- function(z, cols, zlim) {
  breaks <- seq(zlim[[1]], zlim[[2]], length.out = length(cols) + 1)
  cols[findInterval(z[-1, -1], breaks, all.inside = TRUE)]
}

surface_expr <- function() {
  par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2.1, 0.8), oma = c(0, 0, 1.9, 0), family = "sans")
  cols_f <- grDevices::colorRampPalette(c("#F6E9DB", "#E3B27F", "#B96A2B"))(80)
  cols_m <- grDevices::colorRampPalette(c("#E6EEF7", "#8FB6DC", "#2F6DAE"))(80)
  mesh_col <- grDevices::adjustcolor("white", alpha.f = 0.35)

  persp(
    x = female_surface$age_seq,
    y = female_surface$month_seq,
    z = female_surface$z,
    theta = 36,
    phi = 28,
    expand = 0.62,
    zlim = surface_zlim,
    col = surface_colors(female_surface$z, cols_f, surface_zlim),
    border = mesh_col,
    shade = 0.25,
    ticktype = "detailed",
    xlab = "Dam age (months)",
    ylab = "Calving month",
    zlab = "Predicted BW (kg)",
    main = "Female"
  )
  persp(
    x = male_surface$age_seq,
    y = male_surface$month_seq,
    z = male_surface$z,
    theta = 36,
    phi = 28,
    expand = 0.62,
    zlim = surface_zlim,
    col = surface_colors(male_surface$z, cols_m, surface_zlim),
    border = mesh_col,
    shade = 0.25,
    ticktype = "detailed",
    xlab = "Dam age (months)",
    ylab = "Calving month",
    zlab = "Predicted BW (kg)",
    main = "Male"
  )
  mtext(
    sprintf("Support-limited 3D response surface from cyclic GAM predictions; dam age %.1f-%.1f months",
            surface_age_limits[[1]], surface_age_limits[[2]]),
    side = 3,
    line = 0.3,
    outer = TRUE,
    cex = 0.82
  )
}
save_base_plot("3D_Surface_predicted_birth_weight", surface_expr, width = 10.5, height = 5.6)

write_log("14_make_3d_figures.R completed")
