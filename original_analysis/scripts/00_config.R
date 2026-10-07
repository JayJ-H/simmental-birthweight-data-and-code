options(stringsAsFactors = FALSE)
options(contrasts = c("contr.sum", "contr.poly"))

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

DATA_DIR <- file.path(PROJECT_ROOT, "data")
RAW_DIR <- file.path(DATA_DIR, "raw")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "outputs")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
MODEL_DIR <- file.path(OUTPUT_DIR, "models")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
MANUSCRIPT_DIR <- file.path(PROJECT_ROOT, "manuscript_notes")

SEX_RULE_CONFIRMED <- TRUE
SURVIVAL_STATUS_CONFIRMED <- "survived_day_of_birth"
DATA_EXPORT_DATE <- as.Date("2026-07-02")
DAM_AGE_INTERPRETATION <- "calving_date_based"

PRIMARY_DAM_AGE_TERM <- "parity_only"
MIN_SERVICE_SIRE_N <- c(30, 50)
MIN_MGS_N <- 30

create_dir_if_missing <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

for (path in c(DATA_DIR, RAW_DIR, PROCESSED_DIR, OUTPUT_DIR, TABLE_DIR,
               FIGURE_DIR, MODEL_DIR, LOG_DIR, MANUSCRIPT_DIR)) {
  create_dir_if_missing(path)
}

write_log <- function(message, file = file.path(LOG_DIR, "pipeline.log")) {
  create_dir_if_missing(dirname(file))
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", message)
  cat(line, file = file, append = TRUE, sep = "\n")
  message(message)
}

check_package <- function(pkg, required = TRUE) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok && required) {
    stop("Required R package not installed: ", pkg,
         ". Install it with install.packages('", pkg, "').", call. = FALSE)
  }
  ok
}

required_packages <- c(
  "readxl", "dplyr", "tidyr", "stringr", "lubridate", "forcats", "purrr",
  "tibble", "openxlsx", "ggplot2", "patchwork", "cowplot", "scales",
  "emmeans", "car", "lme4", "lmerTest", "broom", "broom.mixed",
  "quantreg", "sandwich", "lmtest", "mgcv", "scatterplot3d", "svglite", "ragg"
)
optional_packages <- c("glmmTMB", "MuMIn", "performance", "logistf", "brglm2")

invisible(lapply(required_packages, check_package, required = TRUE))
for (pkg in optional_packages) {
  if (!check_package(pkg, required = FALSE)) {
    write_log(paste("Optional package unavailable; related robustness output may be skipped:", pkg))
  }
}

detect_raw_file <- function() {
  env_file <- Sys.getenv("RAW_FILE", unset = NA_character_)
  if (!is.na(env_file) && nzchar(env_file)) {
    if (!file.exists(env_file)) stop("RAW_FILE does not exist: ", env_file, call. = FALSE)
    return(normalizePath(env_file, winslash = "/"))
  }
  files <- list.files(RAW_DIR, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No Excel file found in data/raw/.", call. = FALSE)
  if (length(files) > 1) {
    stop("Multiple Excel files found in data/raw/. Set RAW_FILE explicitly.", call. = FALSE)
  }
  normalizePath(files[[1]], winslash = "/")
}

RAW_FILE <- detect_raw_file()

safe_read_excel <- function(path = RAW_FILE, sheet = 1) {
  readxl::read_excel(path, sheet = sheet, guess_max = 10000) |>
    as.data.frame(check.names = FALSE)
}

sanitize_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  substr(x, 1, 31)
}

safe_write_xlsx <- function(x, path) {
  create_dir_if_missing(dirname(path))
  if (is.list(x) && !is.data.frame(x)) {
    names(x) <- make.unique(sanitize_sheet_name(names(x)))
  }
  openxlsx::write.xlsx(x, file = path, overwrite = TRUE, asTable = TRUE)
  invisible(path)
}

write_csv_safe <- function(x, path) {
  create_dir_if_missing(dirname(path))
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

read_csv_safe <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)
}

clean_id <- function(x) {
  if (is.numeric(x)) {
    y <- format(x, scientific = FALSE, trim = TRUE)
  } else {
    y <- as.character(x)
  }
  y <- stringr::str_trim(y)
  y[y %in% c("", "NA", "NaN", "NULL", "<NA>")] <- NA_character_
  y
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

parse_excel_date <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x_chr <- stringr::str_trim(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL", "<NA>")] <- NA_character_
  parsed <- suppressWarnings(lubridate::ymd(x_chr))
  miss <- is.na(parsed) & !is.na(x_chr)
  if (any(miss)) parsed[miss] <- suppressWarnings(lubridate::ymd_hms(x_chr[miss]))
  miss <- is.na(parsed) & !is.na(x_chr)
  if (any(miss)) parsed[miss] <- suppressWarnings(lubridate::mdy(x_chr[miss]))
  as.Date(parsed)
}

subtract_months <- function(date_value, month_value) {
  date_value <- as.Date(date_value)
  n <- max(length(date_value), length(month_value))
  date_value <- rep(date_value, length.out = n)
  month_value <- rep(month_value, length.out = n)
  out <- rep(as.Date(NA), length(date_value))
  ok <- !is.na(date_value) & !is.na(month_value)
  if (any(ok)) {
    dates <- as.POSIXlt(date_value[ok])
    months_back <- round(month_value[ok])
    for (i in seq_along(months_back)) {
      total_month <- (dates$year[i] + 1900L) * 12L + (dates$mon[i] + 1L) - months_back[i]
      target_year <- (total_month - 1L) %/% 12L
      target_month <- (total_month - 1L) %% 12L + 1L
      next_month <- target_month + 1L
      next_year <- target_year
      if (next_month == 13L) {
        next_month <- 1L
        next_year <- target_year + 1L
      }
      first_next <- as.Date(sprintf("%04d-%02d-01", next_year, next_month))
      last_day <- as.integer(format(first_next - 1L, "%d"))
      target_day <- min(dates$mday[i], last_day)
      out[which(ok)[i]] <- as.Date(sprintf("%04d-%02d-%02d", target_year, target_month, target_day))
    }
  }
  as.Date(out)
}

season_from_month <- function(month_value) {
  out <- dplyr::case_when(
    month_value %in% c(3, 4, 5) ~ "Spring",
    month_value %in% c(6, 7, 8) ~ "Summer",
    month_value %in% c(9, 10, 11) ~ "Autumn",
    month_value %in% c(12, 1, 2) ~ "Winter",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("Spring", "Summer", "Autumn", "Winter"), ordered = FALSE)
}

prepare_analysis_types <- function(df) {
  date_cols <- intersect(names(df), c(
    "calving_date", "dam_birth_date_estimated", "dam_birth_date_export_based"
  ))
  for (col in date_cols) df[[col]] <- as.Date(df[[col]])

  numeric_cols <- intersect(names(df), c(
    "dam_age_month", "parity", "birth_weight", "calf_birth_year",
    "calf_birth_month", "dam_birth_year", "dam_birth_month",
    "dam_age_month_export_based_at_calving", "dam_age_month_cont"
  ))
  for (col in numeric_cols) df[[col]] <- to_numeric_safe(df[[col]])

  factor_cols <- intersect(names(df), c(
    "calf_sex", "calf_sex_code", "birth_type", "parity_group", "dam_age_group",
    "calf_birth_year", "calf_birth_season", "dam_birth_year", "dam_birth_season",
    "service_sire", "maternal_grandsire", "dam_id"
  ))
  for (col in factor_cols) df[[col]] <- as.factor(df[[col]])
  df
}

sex_variable_name <- function() {
  if (isTRUE(SEX_RULE_CONFIRMED)) "calf_sex" else "calf_sex_code"
}

sex_label_title <- function() {
  if (isTRUE(SEX_RULE_CONFIRMED)) "male and female" else "G-coded and non-G-coded"
}

summary_stats <- function(data, value = "birth_weight") {
  x <- data[[value]]
  tibble::tibble(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    SD = stats::sd(x, na.rm = TRUE),
    SE = SD / sqrt(n),
    median = stats::median(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    P5 = stats::quantile(x, 0.05, na.rm = TRUE, names = FALSE),
    P10 = stats::quantile(x, 0.10, na.rm = TRUE, names = FALSE),
    P15 = stats::quantile(x, 0.15, na.rm = TRUE, names = FALSE),
    P25 = stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    P75 = stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    P85 = stats::quantile(x, 0.85, na.rm = TRUE, names = FALSE),
    P90 = stats::quantile(x, 0.90, na.rm = TRUE, names = FALSE),
    P95 = stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE)
  )
}

safe_formula <- function(response, fixed_terms, random_terms = NULL) {
  rhs <- paste(fixed_terms[!is.na(fixed_terms) & nzchar(fixed_terms)], collapse = " + ")
  if (!is.null(random_terms) && length(random_terms) > 0) {
    rhs <- paste(c(rhs, random_terms), collapse = " + ")
  }
  stats::as.formula(paste(response, "~", rhs))
}

age_terms_from_choice <- function(choice = PRIMARY_DAM_AGE_TERM) {
  switch(choice,
         parity_age_group = c("parity_group", "dam_age_group"),
         parity_age_cont = c("parity_group", "dam_age_month_cont"),
         parity_only = c("parity_group"),
         age_group_only = c("dam_age_group"),
         c("parity_group"))
}

tidy_lm_anova <- function(model) {
  out <- as.data.frame(car::Anova(model, type = 3))
  out$term <- rownames(out)
  rownames(out) <- NULL
  dplyr::relocate(out, term)
}

or_table_from_model <- function(model) {
  coefs <- as.data.frame(summary(model)$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  names(coefs) <- gsub("Pr\\(>\\|z\\|\\)", "p_value", names(coefs))
  names(coefs) <- gsub("Pr\\(>z\\)", "p_value", names(coefs))
  names(coefs) <- gsub("Std. Error", "std_error", names(coefs), fixed = TRUE)
  names(coefs) <- gsub("Estimate", "estimate", names(coefs), fixed = TRUE)
  if (!"p_value" %in% names(coefs)) {
    p_col <- grep("Pr", names(coefs), value = TRUE)
    if (length(p_col)) names(coefs)[match(p_col[[1]], names(coefs))] <- "p_value"
  }
  coefs |>
    dplyr::mutate(
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std_error),
      CI_high = exp(estimate + 1.96 * std_error)
    ) |>
    dplyr::select(term, estimate, std_error, dplyr::everything())
}

theme_pub <- function(base_size = 7, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.2),
      legend.text = ggplot2::element_text(size = base_size - 0.7),
      strip.text = ggplot2::element_text(size = base_size - 0.2, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}

palette_contract <- c(
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  signal_blue = "#3182BD",
  signal_teal = "#33B5A5",
  accent_red = "#D24B40",
  accent_orange = "#E28E2C"
)

save_plot_pdf_png <- function(plot, filename, width_mm = 183, height_mm = 120, dpi = 600) {
  create_dir_if_missing(dirname(filename))
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  ggplot2::ggsave(paste0(filename, ".pdf"), plot = plot, width = width, height = height,
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(filename, ".png"), plot = plot, width = width, height = height,
                  dpi = dpi)
  invisible(filename)
}

write_log("00_config.R loaded")
