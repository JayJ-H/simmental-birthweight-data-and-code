try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})
root <- normalizePath("..", winslash = "/", mustWork = TRUE)
old <- file.path(root, "original_analysis")
out <- normalizePath(".", winslash = "/")
for (p in c("results", "results/models")) {
  dir.create(file.path(out, p), recursive = TRUE, showWarnings = FALSE)
}
csv <- function(x, name) {
  write.csv(x, file.path(out, "results", paste0(name, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
d <- read.csv(file.path(old, "data/processed/singleton_analysis_with_tail_phenotypes.csv"),
              fileEncoding = "UTF-8", colClasses = c(dam_id = "character",
              calf_id = "character", service_sire = "character",
              maternal_grandsire = "character"))
stopifnot(nrow(d) == 9384, sum(d$lower_tail) == 948, sum(d$upper_tail) == 955)
d$parity_group <- factor(d$parity_group, levels = c("1", "2", "3", "4", "5", "6", ">=7"))
d$calf_birth_year <- factor(d$calf_birth_year, levels = as.character(2023:2026))
d$calf_birth_season <- factor(d$calf_birth_season,
                             levels = c("Spring", "Summer", "Autumn", "Winter"))
d$calf_sex <- factor(d$calf_sex, levels = c("female", "male"))
for (v in c("dam_id", "service_sire", "maternal_grandsire")) d[[v]] <- factor(d[[v]])
fixed <- c("parity_group", "calf_birth_year", "calf_birth_season")
refs <- c(parity_group = "1", calf_birth_year = "2023",
          calf_birth_season = "Spring", calf_sex = "female")
effects <- function(m, id, outcome) {
  do.call(rbind, lapply(fixed, function(v) {
    e <- emmeans(m, reformulate(v))
    ref <- match(refs[[v]], as.character(as.data.frame(e)[[v]]))
    stopifnot(!is.na(ref))
    z <- as.data.frame(summary(contrast(e, "trt.vs.ctrl", ref = ref, adjust = "none"),
                               infer = c(TRUE, TRUE), type = "link"))
    z$contrast <- gsub(v, "", as.character(z$contrast), fixed = TRUE)
    z$contrast <- gsub("[()]", "", z$contrast)
    lo <- intersect(c("asymp.LCL", "lower.CL"), names(z))[1]
    hi <- intersect(c("asymp.UCL", "upper.CL"), names(z))[1]
    data.frame(model = id, outcome = outcome, factor = v, comparison = z$contrast,
               log_OR = z$estimate, SE = z$SE, OR = exp(z$estimate),
               CI_low = exp(z[[lo]]), CI_high = exp(z[[hi]]), p = z$p.value)
  }))
}
jobs <- list()
add_job <- function(id, outcome, response, subset, extra = NULL, cache = NULL) {
  jobs[[length(jobs) + 1L]] <<- list(id = id, outcome = outcome, response = response,
                                   rows = which(subset), extra = extra, cache = cache)
}
for (o in c("lower", "upper")) {
  y <- paste0(o, "_tail")
  add_job(paste0("primary_", o), o, y, rep(TRUE, nrow(d)),
          cache = file.path(old, "outputs/models", paste0("glmm_", o, "_tail_main.rds")))
  add_job(paste0("p5p95_", o), o, if(o == "lower") "lower_tail_p5" else "upper_tail_p95", rep(TRUE, nrow(d)))
  add_job(paste0("p15p85_", o), o, if(o == "lower") "lower_tail_p15" else "upper_tail_p85", rep(TRUE, nrow(d)))
  add_job(paste0("complete_years_", o), o, y, d$calf_birth_year != "2026")
  add_job(paste0("weight_25_65_", o), o, y, d$birth_weight >= 25 & d$birth_weight <= 65)
  add_job(paste0("sire_subset_base_", o), o, y, !is.na(d$service_sire))
  add_job(paste0("sire_subset_adjusted_", o), o, y, !is.na(d$service_sire), "service_sire")
  add_job(paste0("mgs_subset_base_", o), o, y, !is.na(d$maternal_grandsire))
  add_job(paste0("mgs_subset_adjusted_", o), o, y, !is.na(d$maternal_grandsire), "maternal_grandsire")
}
fit_job <- function(job) {
  try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE)
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(emmeans))
  dd <- droplevels(d[job$rows, ])
  saved <- file.path(out, "results/models", paste0(job$id, ".rds"))
  warnings <- character()
  if (!is.null(job$cache)) {
    m <- readRDS(job$cache)
    origin <- "original saved main model; reference contrasts recomputed"
  } else if (file.exists(saved)) {
    m <- readRDS(saved)
    origin <- "revision cached fit"
  } else {
    f <- reformulate(c(fixed, "(1 | dam_id)",
                       if(!is.null(job$extra)) paste0("(1 | ", job$extra, ")")),
                     response = job$response)
    m <- withCallingHandlers(
      glmer(f, data = dd, family = binomial(),
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))),
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}
    )
    saveRDS(m, saved)
    origin <- "new revision fit"
  }
  stopifnot(nobs(m) == nrow(dd))
  conv <- m@optinfo$conv$lme4$messages
  var <- as.data.frame(VarCorr(m))
  diag <- data.frame(model = job$id, outcome = job$outcome, n = nobs(m),
                     events = sum(dd[[job$response]]),
                     dams = length(unique(dd$dam_id)),
                     extra_groups = if(is.null(job$extra)) NA_integer_ else nlevels(dd[[job$extra]]),
                     nAGQ = unname(getME(m,"devcomp")$dims[["nAGQ"]]),
                     singular = isSingular(m, tol = 1e-4),
                     optimizer_code = paste(m@optinfo$conv$opt, collapse = ";"),
                     messages = paste(unique(c(conv, warnings)), collapse = "; "),
                     variance = paste(var$grp, signif(var$vcov, 6), sep = "=", collapse = "; "),
                     formula = paste(deparse(formula(m)), collapse = " "),
                     origin = origin)
  write.csv(diag, file.path(out, "results/models", paste0(job$id, "_diagnostics.csv")), row.names = FALSE)
  list(diagnostics = diag, coefficients = effects(m, job$id, job$outcome))
}
cl <- parallel::makePSOCKcluster(4)
parallel::clusterExport(cl, c("d", "out", "fixed", "refs", "effects", "fit_job"))
cat("Fitting", length(jobs), "GLMM scenarios on four workers\n")
fits <- tryCatch(parallel::parLapplyLB(cl, jobs, fit_job),
                 finally = parallel::stopCluster(cl))
diag <- do.call(rbind, lapply(fits, `[[`, "diagnostics"))
eff <- do.call(rbind, lapply(fits, `[[`, "coefficients"))
csv(diag, "model_diagnostics")
csv(eff, "tail_model_effects")
stopifnot(all(diag$optimizer_code == "0"), all(diag$n > 0))
cat("GLMM scenarios complete\n")

# Primary LMM table comes entirely from the archived mixed model, not OLS.
lmm <- readRDS(file.path(old, "outputs/models/lmm_birth_weight_final.rds"))
anova <- as.data.frame(anova(lmm, type = 3))
anova$factor <- rownames(anova)
csv(anova, "lmm_overall_tests")
means <- list()
pairs <- list()
for (v in c("calf_sex", fixed)) {
  e <- emmeans(lmm, reformulate(v), lmer.df = "asymptotic")
  z <- as.data.frame(e)
  cc <- as.data.frame(summary(pairs(e, adjust = "tukey"), infer = c(TRUE, TRUE)))
  lv <- as.character(z[[v]])
  mat <- matrix(1, length(lv), length(lv), dimnames = list(lv, lv))
  for(i in seq_len(nrow(cc))) {
    nm <- strsplit(as.character(cc$contrast[i]), " - ", fixed = TRUE)[[1]]
    nm <- gsub("^\\(|\\)$", "", nm)
    nm <- gsub(v, "", nm, fixed = TRUE)
    if(!all(nm %in% lv)) stop("Cannot resolve LMM contrast: ", cc$contrast[i])
    mat[nm[1], nm[2]] <- mat[nm[2], nm[1]] <- cc$p.value[i]
  }
  ord <- order(z$emmean, decreasing = TRUE)
  letters <- multcompView::multcompLetters(mat[ord, ord, drop = FALSE], threshold = .05)$Letters
  pcol <- grep("Pr", names(anova), value = TRUE)[1]
  means[[v]] <- data.frame(factor = v, level = lv, mean = z$emmean,
                          SE = z$SE, CI_low = z$asymp.LCL, CI_high = z$asymp.UCL,
                          letters = unname(letters[lv]),
                          overall_p = anova[v, pcol])
  cc$factor <- v
  pairs[[v]] <- cc
}
csv(do.call(rbind, means), "table3_lmm_means")
csv(do.call(rbind, pairs), "table3_tukey_comparisons")
csv(as.data.frame(VarCorr(lmm)), "lmm_variance")
csv(eff[grepl("^primary_", eff$model), ], "table4_primary_OR")

year_stats <- do.call(rbind, lapply(split(d, d$calf_birth_year), function(x) {
  b <- x$birth_weight
  data.frame(year = as.character(x$calf_birth_year[1]), n = nrow(x), mean = mean(b),
             SD = sd(b), median = median(b), Q1 = unname(quantile(b, .25)),
             Q3 = unname(quantile(b, .75)), IQR = IQR(b),
             min = min(b), max = max(b), CV_percent = 100 * sd(b)/mean(b),
             lower_n = sum(x$lower_tail), upper_n = sum(x$upper_tail),
             lower_percent = 100 * mean(x$lower_tail), upper_percent = 100 * mean(x$upper_tail))
}))
csv(year_stats, "annual_variability")
csv(as.data.frame(table(year = d$calf_birth_year, season = d$calf_birth_season)),
    "year_season_counts")
csv(d[d$birth_weight < 25 | d$birth_weight > 65,
      c("row_id", "dam_id", "calf_id", "birth_weight", "calving_date", "calf_sex",
        "service_sire", "maternal_grandsire")], "extreme_records_83")
thresholds <- do.call(rbind, lapply(split(d, d$calf_sex), function(x) {
  q <- quantile(x$birth_weight, c(.05,.1,.15,.85,.9,.95))
  data.frame(sex = as.character(x$calf_sex[1]), n = nrow(x),
             P5=q[1], P10=q[2], P15=q[3], P85=q[4], P90=q[5], P95=q[6])
}))
csv(thresholds, "thresholds")
missing_summary <- list()
missing_categories <- list()
label_counts <- list()
for (v in c("service_sire", "maternal_grandsire")) {
  for (status in c("available", "missing")) {
    dd <- d[if(status == "available") !is.na(d[[v]]) else is.na(d[[v]]), ]
    missing_summary[[paste(v,status)]] <- data.frame(field = v, status = status,
      n = nrow(dd), mean = mean(dd$birth_weight), SD = sd(dd$birth_weight),
      lower_percent = 100*mean(dd$lower_tail), upper_percent = 100*mean(dd$upper_tail))
    for (f in c("calf_sex", fixed)) {
      tt <- as.data.frame(table(dd[[f]]))
      names(tt) <- c("level", "n")
      missing_categories[[paste(v,status,f)]] <- data.frame(field=v,status=status,
        factor=f, tt, percent=100*tt$n/nrow(dd))
    }
  }
  tab <- as.data.frame(table(label = d[[v]], year = d$calf_birth_year))
  tab$field <- v
  label_counts[[v]] <- tab[tab$Freq > 0, ]
}
csv(do.call(rbind, missing_summary), "ancestry_missingness_summary")
csv(do.call(rbind, missing_categories), "ancestry_missingness_categories")
csv(do.call(rbind, label_counts), "ancestry_labels_by_year")

# Compare effects on the identical sample when assessing ancestry-label adjustment.
comparisons <- list()
for (model in unique(eff$model)) {
  if(grepl("^primary_", model)) next
  sub <- eff[eff$model == model, ]
  target <- if(grepl("_adjusted_", model)) sub("_adjusted_", "_base_", model) else
    paste0("primary_", sub$outcome[1])
  base <- eff[eff$model == target, c("factor", "comparison", "log_OR", "OR")]
  names(base)[3:4] <- c("reference_log_OR","reference_OR")
  z <- merge(sub, base, by=c("factor","comparison"))
  z$reference_model <- target
  z$same_direction <- sign(z$log_OR) == sign(z$reference_log_OR)
  z$OR_ratio <- z$OR/z$reference_OR
  comparisons[[model]] <- z
}
csv(do.call(rbind, comparisons), "effect_comparisons")

qr_file <- file.path(out, "results/qr_bootstrap.rds")
X <- model.matrix(~ calf_sex + parity_group + calf_birth_year + calf_birth_season, d,
                  contrasts.arg = list(calf_sex="contr.treatment", parity_group="contr.treatment",
                    calf_birth_year="contr.treatment", calf_birth_season="contr.treatment"))
Y <- d$birth_weight
taus <- c(.1, .5, .9)
groups <- split(seq_len(nrow(d)), d$dam_id, drop = TRUE)
B <- 999L
point <- sapply(taus, function(tau) quantreg::rq.fit.fnb(X, Y, tau = tau)$coefficients)
boot_one <- function(i) {
  set.seed(20261005L + i)
  ix <- unlist(groups[sample.int(length(groups), length(groups), replace=TRUE)],
               use.names=FALSE)
  warn <- character()
  ans <- tryCatch(withCallingHandlers(
    sapply(taus, function(tau) quantreg::rq.fit.fnb(X[ix,,drop=FALSE], Y[ix], tau=tau)$coefficients),
    warning=function(w){warn <<- c(warn,conditionMessage(w)); invokeRestart("muffleWarning")}),
    error=function(e)e)
  if(inherits(ans,"error")) return(list(i=i, values=NULL, message=ans$message))
  list(i=i, values=ans, message=paste(unique(warn),collapse="; "))
}
if(file.exists(qr_file)) {
  boot <- readRDS(qr_file)
} else {
  cat("Starting", B, "dam-level pairs-bootstrap replicates at three quantiles\n")
  cl <- parallel::makePSOCKcluster(4)
  parallel::clusterExport(cl, c("X","Y","taus","groups","boot_one"))
  boot <- tryCatch(parallel::parLapplyLB(cl, seq_len(B), boot_one),
                   finally=parallel::stopCluster(cl))
  saveRDS(boot, qr_file)
}
needs_refit <- which(vapply(boot, function(z) nzchar(z$message) &&
  is.null(z$refit_method), logical(1)))
if(length(needs_refit)) {
  backup <- file.path(out, "results/qr_bootstrap_initial.rds")
  if(!file.exists(backup)) file.copy(qr_file, backup)
  for(i in needs_refit) {
    set.seed(20261005L + i)
    ix <- unlist(groups[sample.int(length(groups), length(groups), replace=TRUE)],
                 use.names=FALSE)
    stopifnot(qr(X[ix,,drop=FALSE])$rank == ncol(X))
    warning_br <- character()
    refit <- withCallingHandlers(
      sapply(taus, function(tau)
        quantreg::rq.fit.br(X[ix,,drop=FALSE], Y[ix], tau=tau)$coefficients),
      warning=function(w) {
        warning_br <<- c(warning_br, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    stopifnot(all(is.finite(refit)),
              !any(grepl("singular|converg|error", warning_br, ignore.case=TRUE)))
    loss <- function(b, j) {
      r <- Y[ix] - as.vector(X[ix,,drop=FALSE] %*% b)
      sum(r * (taus[j] - (r < 0)))
    }
    old_loss <- vapply(seq_along(taus), function(j) loss(boot[[i]]$values[,j],j), numeric(1))
    new_loss <- vapply(seq_along(taus), function(j) loss(refit[,j],j), numeric(1))
    stopifnot(all(new_loss <= old_loss + 1e-5 * (1 + abs(old_loss))))
    boot[[i]]$initial_message <- boot[[i]]$message
    boot[[i]]$message <- paste(unique(warning_br), collapse="; ")
    boot[[i]]$max_coefficient_change <- max(abs(refit-boot[[i]]$values))
    boot[[i]]$max_loss_change <- max(abs(new_loss-old_loss))
    boot[[i]]$refit_method <- "Barrodale-Roberts; full-rank design; objective checked"
    boot[[i]]$values <- refit
    cat("Verified warned QR bootstrap replicate", i, "\n")
  }
  saveRDS(boot, qr_file)
}
valid <- vapply(boot, function(z) !is.null(z$values) && all(is.finite(z$values)), logical(1))
csv(data.frame(replicate=seq_len(B),success=valid,
               message=vapply(boot,`[[`,character(1),"message"),
               initial_message=vapply(boot,function(z) if(is.null(z$initial_message)) "" else z$initial_message,character(1)),
               refit_method=vapply(boot,function(z) if(is.null(z$refit_method)) "" else z$refit_method,character(1)),
               max_coefficient_change=vapply(boot,function(z) if(is.null(z$max_coefficient_change)) NA_real_ else z$max_coefficient_change,numeric(1)),
               max_loss_change=vapply(boot,function(z) if(is.null(z$max_loss_change)) NA_real_ else z$max_loss_change,numeric(1))),
    "qr_bootstrap_diagnostics")
stopifnot(sum(valid) >= .99*B)
qr <- do.call(rbind, lapply(seq_along(taus), function(j) {
  mat <- do.call(rbind, lapply(boot[valid], function(z) z$values[,j]))
  data.frame(tau=taus[j],term=colnames(X),estimate=point[,j],
             bootstrap_SE=apply(mat,2,sd),
             CI_low=apply(mat,2,quantile,.025),CI_high=apply(mat,2,quantile,.975),
             successful_replicates=sum(valid), n=nrow(d), dams=length(groups),
             method="dam-level pairs bootstrap; empirical percentile 95% CI")
}))
csv(qr, "quantile_regression_cluster_bootstrap")
writeLines(capture.output(sessionInfo()), file.path(out,"results/session_info.txt"))
cat("All analyses complete; QR successes:",sum(valid),"/",B,"\n")
