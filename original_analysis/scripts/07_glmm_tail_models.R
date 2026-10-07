source("scripts/00_config.R")

log_file <- file.path(LOG_DIR, "glmm_model_log.txt")
warning_file <- file.path(LOG_DIR, "model_warnings.log")
write_log("07_glmm_tail_models.R started", log_file)

singleton_tail <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv")) |> prepare_analysis_types()

main_terms <- c(age_terms_from_choice(), "calf_birth_year", "calf_birth_season")

capture_glmer <- function(formula, data) {
  warnings <- character()
  model <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        formula,
        data = data,
        family = stats::binomial(link = "logit"),
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(model = model, warnings = unique(warnings))
}

tidy_glm_like <- function(model, model_name) {
  coefs <- as.data.frame(summary(model)$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  names(coefs) <- gsub("Estimate", "estimate", names(coefs), fixed = TRUE)
  names(coefs) <- gsub("Std. Error", "std_error", names(coefs), fixed = TRUE)
  z_col <- grep("z value|t value", names(coefs), value = TRUE)
  p_col <- grep("Pr", names(coefs), value = TRUE)
  if (length(z_col)) names(coefs)[match(z_col[[1]], names(coefs))] <- "statistic"
  if (length(p_col)) names(coefs)[match(p_col[[1]], names(coefs))] <- "p_value"
  coefs |>
    dplyr::mutate(
      model = model_name,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std_error),
      CI_high = exp(estimate + 1.96 * std_error),
      .before = 1
    )
}

fit_cluster_glm <- function(response, fixed_terms, data, model_name) {
  formula <- safe_formula(response, fixed_terms)
  model_data <- data |>
    dplyr::select(dplyr::all_of(c(all.vars(formula), "dam_id"))) |>
    stats::na.omit()
  glm_fit <- stats::glm(formula, data = model_data, family = stats::binomial())
  robust <- tryCatch(lmtest::coeftest(glm_fit, vcov. = sandwich::vcovCL(glm_fit, cluster = model_data$dam_id)), error = function(e) e)
  if (inherits(robust, "error")) {
    table <- tidy_glm_like(glm_fit, paste0(model_name, "_glm_fallback")) |>
      dplyr::mutate(note = paste("cluster robust SE failed:", robust$message))
  } else {
    table <- as.data.frame(robust) |>
      tibble::rownames_to_column("term")
    names(table) <- c("term", "estimate", "std_error", "statistic", "p_value")
    table <- table |>
      dplyr::mutate(
        model = paste0(model_name, "_glm_cluster_robust"),
        OR = exp(estimate),
        CI_low = exp(estimate - 1.96 * std_error),
        CI_high = exp(estimate + 1.96 * std_error),
        note = "ordinary logistic regression with cluster-robust SE by dam_id",
        .before = 1
      )
  }
  list(model = glm_fit, data = model_data, table = table)
}

fit_tail_model <- function(response, fixed_terms, random_terms, data, model_name) {
  formula <- safe_formula(response, fixed_terms, random_terms)
  model_data <- data |>
    dplyr::select(dplyr::all_of(all.vars(formula))) |>
    stats::na.omit()
  if (sum(model_data[[response]] == 1) < 10 || sum(model_data[[response]] == 0) < 10) {
    return(list(
      status = "insufficient events",
      model_name = model_name,
      formula = formula,
      data = model_data,
      table = tibble::tibble(model = model_name, status = "insufficient events")
    ))
  }
  fit <- capture_glmer(formula, model_data)
  if (inherits(fit$model, "error")) {
    write_log(paste(model_name, "GLMM failed:", fit$model$message), warning_file)
    fallback <- fit_cluster_glm(response, fixed_terms, model_data, model_name)
    return(list(status = "glm_cluster_fallback", model_name = model_name, formula = formula, data = model_data, model = fallback$model, table = fallback$table, warnings = fit$model$message))
  }
  singular <- lme4::isSingular(fit$model, tol = 1e-4)
  if (length(fit$warnings) > 0) write_log(paste(model_name, paste(fit$warnings, collapse = " | ")), warning_file)
  if (singular) write_log(paste(model_name, "singular fit; ordinary logistic cluster-robust fallback also reported"), warning_file)
  table <- tidy_glm_like(fit$model, model_name) |>
    dplyr::mutate(note = ifelse(singular, "GLMM singular fit", "GLMM"))
  fallback_table <- if (singular) fit_cluster_glm(response, fixed_terms, model_data, model_name)$table else tibble::tibble()
  list(
    status = ifelse(singular, "fit_singular", "fit"),
    model_name = model_name,
    formula = formula,
    data = model_data,
    model = fit$model,
    table = dplyr::bind_rows(table, fallback_table),
    warnings = fit$warnings
  )
}

lower_main <- fit_tail_model("lower_tail", main_terms, "(1 | dam_id)", singleton_tail, "lower_tail_main")
upper_main <- fit_tail_model("upper_tail", main_terms, "(1 | dam_id)", singleton_tail, "upper_tail_main")

service_data <- singleton_tail |> dplyr::filter(!is.na(service_sire))
lower_service <- fit_tail_model("lower_tail", main_terms, c("(1 | dam_id)", "(1 | service_sire)"), service_data, "lower_tail_plus_service_sire")
upper_service <- fit_tail_model("upper_tail", main_terms, c("(1 | dam_id)", "(1 | service_sire)"), service_data, "upper_tail_plus_service_sire")

model_diag <- function(fit) {
  if (is.null(fit$model) || inherits(fit$model, "glm")) {
    return(tibble::tibble(model = fit$model_name, status = fit$status, formula = paste(deparse(fit$formula), collapse = " "), n = nrow(fit$data), AIC = if (!is.null(fit$model)) stats::AIC(fit$model) else NA_real_, BIC = if (!is.null(fit$model)) stats::BIC(fit$model) else NA_real_, random_effect_variance = NA_real_, warnings = paste(fit$warnings, collapse = " | ")))
  }
  var_df <- as.data.frame(lme4::VarCorr(fit$model))
  tibble::tibble(
    model = fit$model_name,
    status = fit$status,
    formula = paste(deparse(fit$formula), collapse = " "),
    n = nrow(fit$data),
    AIC = stats::AIC(fit$model),
    BIC = stats::BIC(fit$model),
    random_effect_variance = paste(var_df$grp, round(var_df$vcov, 6), sep = "=", collapse = "; "),
    warnings = paste(fit$warnings, collapse = " | ")
  )
}

lower_results <- list(lower_main, lower_service)
upper_results <- list(upper_main, upper_service)

safe_write_xlsx(
  list(
    diagnostics = purrr::map_dfr(lower_results, model_diag),
    OR_table = purrr::map_dfr(lower_results, "table")
  ),
  file.path(MODEL_DIR, "glmm_lower_tail_results.xlsx")
)
safe_write_xlsx(
  list(
    diagnostics = purrr::map_dfr(upper_results, model_diag),
    OR_table = purrr::map_dfr(upper_results, "table")
  ),
  file.path(MODEL_DIR, "glmm_upper_tail_results.xlsx")
)

forest_data <- dplyr::bind_rows(
  purrr::map_dfr(lower_results, "table") |> dplyr::mutate(phenotype = "lower_tail"),
  purrr::map_dfr(upper_results, "table") |> dplyr::mutate(phenotype = "upper_tail")
) |>
  dplyr::filter(!term %in% "(Intercept)", grepl("_main$", model))
write_csv_safe(forest_data, file.path(PROCESSED_DIR, "glmm_tail_forest_data.csv"))

saveRDS(lower_main$model, file.path(MODEL_DIR, "glmm_lower_tail_main.rds"))
saveRDS(upper_main$model, file.path(MODEL_DIR, "glmm_upper_tail_main.rds"))

write_log("07_glmm_tail_models.R completed", log_file)
