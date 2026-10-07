source("scripts/00_config.R")

write_log("04_least_squares_fixed_effects.R started")

all_valid_records <- read_csv_safe(file.path(PROCESSED_DIR, "all_valid_records.csv")) |> prepare_analysis_types()
singleton_analysis <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis.csv")) |> prepare_analysis_types()

sex_var <- sex_variable_name()
age_terms <- age_terms_from_choice()

fit_ols <- function(data, model_name, fixed_terms) {
  formula <- safe_formula("birth_weight", fixed_terms)
  model_data <- data |>
    dplyr::select(dplyr::all_of(all.vars(formula))) |>
    stats::na.omit()
  model <- stats::lm(formula, data = model_data)
  list(
    model_name = model_name,
    formula = formula,
    data = model_data,
    model = model,
    tidy = broom::tidy(model) |> dplyr::mutate(model = model_name, .before = 1),
    glance = broom::glance(model) |> dplyr::mutate(model = model_name, formula = paste(deparse(formula), collapse = " "), n = nrow(model_data), .before = 1),
    anova = tidy_lm_anova(model) |> dplyr::mutate(model = model_name, .before = 1)
  )
}

model_1_terms <- c(sex_var, "birth_type", age_terms, "calf_birth_year", "calf_birth_season")
model_2_terms <- c(sex_var, age_terms, "calf_birth_year", "calf_birth_season")
model_3_terms <- c(sex_var, age_terms, "calf_birth_year", "calf_birth_season", "dam_birth_year", "dam_birth_season")

fits <- list(
  all_valid_with_birth_type = fit_ols(all_valid_records, "all_valid_with_birth_type", model_1_terms),
  singleton_main = fit_ols(singleton_analysis, "singleton_main", model_2_terms),
  singleton_supplementary_dam_birth = fit_ols(singleton_analysis, "singleton_supplementary_dam_birth", model_3_terms)
)

emm_for_model <- function(fit) {
  terms <- attr(stats::terms(fit$model), "term.labels")
  terms <- terms[!grepl(":", terms)]
  purrr::map_dfr(terms, function(term) {
    if (is.numeric(fit$data[[term]])) return(tibble::tibble())
    emm <- tryCatch(emmeans::emmeans(fit$model, stats::as.formula(paste("~", term))), error = function(e) e)
    if (inherits(emm, "error")) {
      return(tibble::tibble(model = fit$model_name, term = term, status = emm$message))
    }
    as.data.frame(emm) |>
      dplyr::mutate(model = fit$model_name, term = term, status = "ok", .before = 1)
  })
}

pairs_for_model <- function(fit) {
  terms <- attr(stats::terms(fit$model), "term.labels")
  terms <- terms[!grepl(":", terms)]
  purrr::map_dfr(terms, function(term) {
    if (is.numeric(fit$data[[term]])) return(tibble::tibble())
    if (dplyr::n_distinct(fit$data[[term]]) < 2) return(tibble::tibble())
    emm <- tryCatch(emmeans::emmeans(fit$model, stats::as.formula(paste("~", term))), error = function(e) e)
    if (inherits(emm, "error")) {
      return(tibble::tibble(model = fit$model_name, term = term, status = emm$message))
    }
    pairs <- tryCatch(as.data.frame(emmeans::pairs(emm, adjust = "tukey")), error = function(e) e)
    if (inherits(pairs, "error")) {
      return(tibble::tibble(model = fit$model_name, term = term, status = pairs$message))
    }
    pairs |>
      dplyr::mutate(model = fit$model_name, term = term, status = "ok", .before = 1)
  })
}

anova_all <- purrr::map_dfr(fits, "anova") |>
  dplyr::mutate(
    significance = dplyr::case_when(
      `Pr(>F)` < 0.05 ~ "significant",
      `Pr(>F)` < 0.10 ~ "tendency",
      TRUE ~ "not significant"
    )
  )

safe_write_xlsx(
  list(
    model_summary = purrr::map_dfr(fits, "glance"),
    coefficients = purrr::map_dfr(fits, "tidy"),
    type_III_ANOVA = anova_all,
    estimated_marginal_means = purrr::map_dfr(fits, emm_for_model),
    tukey_comparisons = purrr::map_dfr(fits, pairs_for_model)
  ),
  file.path(TABLE_DIR, "least_squares_fixed_effects.xlsx")
)

saveRDS(fits$singleton_main$model, file.path(MODEL_DIR, "ols_singleton_main.rds"))

write_log("04_least_squares_fixed_effects.R completed")
