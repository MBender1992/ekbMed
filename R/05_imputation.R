# Multiple imputation helpers

#' Create reusable clinical imputations
#'
#' Explicitly selects targets and predictors for mice chained-equation imputation.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param impute_vars Names of variables explicitly selected for imputation. Treatment, event and follow-up cannot be imputation targets.
#' @param predictor_vars Additional predictor names for the imputation models.
#' @param auxiliary_vars Names of auxiliary imputation predictors.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param survival_imputation Whether to include the event indicator and Nelson-Aalen cumulative hazard as imputation predictors.
#' @param m Number of completed imputations.
#' @param maxit Number of chained-equation iterations.
#' @param seed Random seed supplied to mice for reproducibility.
#' @param methods Optional named character vector overriding mice methods for variables in impute_vars.
#' @param printFlag Whether mice prints iteration progress.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return An `ekb_imputation` list containing mids, original_data, imputation_data, methods, predictorMatrix and reproducibility settings.
#' @details Imputation is performed once upstream. Targets predict one another; additional explicit predictors are included. With survival_imputation enabled, the event indicator and Nelson-Aalen hazard are predictors; raw survival time is not automatically used. Treatment/time/event must be observed when supplied. Missing predictors must themselves be listed as imputation targets. Additional arguments go to mice::mice. Assess convergence and plausibility of imputations; missing-at-random assumptions are not verified automatically.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("mice", quietly = TRUE)) {
#' d$age[c(3, 12, 25, 40)] <- NA
#' mi <- impute_clinical(d, "age", predictor_vars = "sex",
#'   treatment = "treatment", time = "time", event = "event",
#'   m = 2, maxit = 2, seed = 21)
#'   mi$m
#' }
impute_clinical <- function(data,
                            impute_vars,
                            predictor_vars = NULL,
                            auxiliary_vars = NULL,
                            treatment = NULL,
                            time = NULL,
                            event = NULL,
                            survival_imputation = !is.null(time) && !is.null(event),
                            m = 20,
                            maxit = 10,
                            seed = 1234,
                            methods = NULL,
                            printFlag = FALSE,
                            ...) {
  .ekb_require("mice")
  dat <- as.data.frame(data)
  .ekb_assert_columns(dat, c(impute_vars, predictor_vars, auxiliary_vars, treatment, time, event))
  if (!length(impute_vars)) stop("impute_vars must be specified explicitly.", call. = FALSE)
  if (any(c(time, event, treatment) %in% impute_vars, na.rm = TRUE)) {
    stop("time, event and treatment should not be included in impute_vars in this workflow.", call. = FALSE)
  }
  if (!is.null(event)) dat[[event]] <- .ekb_binary_event(dat[[event]], event)

  checked_predictors <- setdiff(unique(c(auxiliary_vars, predictor_vars, treatment, time, event)), impute_vars)
  na_aux <- checked_predictors[vapply(checked_predictors, function(v) anyNA(dat[[v]]), logical(1))]
  if (length(na_aux)) {
    stop(
      sprintf("Missing predictors must be explicitly imputed; treatment/time/event must be observed: %s", paste(na_aux, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.null(time)) .ekb_validate_time(dat, time)
  hazard_name <- NULL
  if (isTRUE(survival_imputation)) {
    if (is.null(time) || is.null(event)) stop("time and event are required when survival_imputation = TRUE.", call. = FALSE)
    hazard_name <- ".ekb_nelson_aalen"
    if (hazard_name %in% names(dat)) stop("Reserved column name '.ekb_nelson_aalen' already exists.", call. = FALSE)
    dat[[hazard_name]] <- do.call(mice::nelsonaalen, list(data = dat, timevar = time, statusvar = event))
  }

  pred_vars <- unique(c(predictor_vars, auxiliary_vars, treatment))
  if (isTRUE(survival_imputation)) pred_vars <- unique(c(pred_vars, event, hazard_name))
  pred_vars <- setdiff(pred_vars, impute_vars)
  .ekb_assert_columns(dat, pred_vars)

  method <- mice::make.method(dat)
  method[] <- ""
  auto <- mice::make.method(dat)
  method[impute_vars] <- auto[impute_vars]
  if (!is.null(methods)) {
    if (is.null(names(methods))) stop("methods must be a named character vector.", call. = FALSE)
    bad <- setdiff(names(methods), impute_vars)
    if (length(bad)) stop("methods may only override variables listed in impute_vars.", call. = FALSE)
    method[names(methods)] <- methods
  }

  pm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat), dimnames = list(names(dat), names(dat)))
  for (v in impute_vars) {
    preds <- setdiff(unique(c(pred_vars, setdiff(impute_vars, v))), v)
    preds <- preds[preds %in% names(dat)]
    pm[v, preds] <- 1L
  }

  mids <- mice::mice(
    data = dat,
    m = m,
    maxit = maxit,
    method = method,
    predictorMatrix = pm,
    seed = seed,
    printFlag = printFlag,
    ...
  )

  structure(
    list(
      mids = mids,
      original_data = as.data.frame(data),
      imputation_data = dat,
      impute_vars = impute_vars,
      predictor_vars = pred_vars,
      auxiliary_vars = auxiliary_vars,
      treatment = treatment,
      time = time,
      event = event,
      survival_imputation = survival_imputation,
      nelson_aalen_variable = hazard_name,
      method = method,
      predictorMatrix = pm,
      m = m,
      maxit = maxit,
      seed = seed,
      call = match.call()
    ),
    class = "ekb_imputation"
  )
}

.ekb_get_mids <- function(imputation) {
  mi <- if (inherits(imputation, "ekb_imputation")) imputation$mids else imputation
  if (!inherits(mi, "mids")) {
    stop("imputation must be an ekb_imputation or mice::mids object.", call. = FALSE)
  }
  mi
}

.ekb_mi_diagnostic_summary <- function(diagnostics) {
  if (is.null(diagnostics) || !length(diagnostics)) return(NULL)
  out <- do.call(rbind, lapply(seq_along(diagnostics), function(i) {
    d <- diagnostics[[i]]
    data.frame(
      imputation = i,
      max_abs_smd = d$max_abs_smd,
      min_ess = min(d$ess$ess, na.rm = TRUE),
      positivity_flag = d$positivity_flag,
      stringsAsFactors = FALSE
    )
  }))
  .ekb_as_tibble(out)
}

#' Pool Cox models across existing imputations
#'
#' Reuses an upstream imputation, optionally estimating propensity scores and weights separately within each completed dataset.
#'
#' @param imputation An existing `ekb_imputation` or `mice::mids` object. It is reused and never re-imputed downstream.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param covariates Character vector of prespecified outcome-model covariates. Include treatment explicitly for treatment-effect models.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param ps_covariates Character vector of baseline propensity-score predictors. In MI functions, NULL disables weighting; supplying predictors enables within-imputation weighting.
#' @param estimand Target estimand passed to [WeightIt::weightit()], default `"ATE"`.
#' @param weight_method Propensity-score estimation method passed to WeightIt.
#' @param stabilize Whether to request stabilized weights from WeightIt.
#' @param treatment_reference Optional reference level. Otherwise the first existing factor level (or default factor ordering) is used.
#' @param robust NULL delegates variance selection to [survival::coxph()]; non-integer weights normally trigger robust SE. TRUE explicitly requests sandwich variance, including with integer weights; FALSE requests model-based variance.
#' @param conf.level Confidence level between zero and one.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return An `ekb_mi_cox` list with imputation, fits, pooled, pooled_summary, treatment_effect and weighting diagnostics.
#' @details Each completed imputation is analysed separately. Supplying ps_covariates estimates a fresh propensity model and weights in every imputation; otherwise models are unweighted. Additional arguments are passed to WeightIt::weightit. Rubin pooling uses mice::pool, which uses robust.se from broom when present (mice >= 3.13.2). Set robust = TRUE to guarantee robust variance even for integer weights. A completed imputation is not a complete-case analysis. The returned imputation is the original upstream object.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("mice", quietly = TRUE)) {
#' d$age[c(3, 12, 25, 40)] <- NA
#' mi <- impute_clinical(d, "age", predictor_vars = "sex",
#'   treatment = "treatment", time = "time", event = "event",
#'   m = 2, maxit = 2, seed = 21)
#'   fit_mi_cox(mi, "time", "event", c("treatment", "age"))
#' }
fit_mi_cox <- function(imputation,
                       time,
                       event,
                       covariates,
                       treatment = NULL,
                       ps_covariates = NULL,
                       estimand = "ATE",
                       weight_method = "glm",
                       stabilize = FALSE,
                       treatment_reference = NULL,
                       robust = NULL,
                       conf.level = 0.95,
                       ...) {
  .ekb_require("mice")
  .ekb_require("broom")

  mi <- .ekb_get_mids(imputation)
  d1 <- mice::complete(mi, 1)
  .ekb_assert_columns(d1, c(time, event, covariates, treatment, ps_covariates))
  if (!length(covariates)) stop("covariates must contain at least one outcome-model covariate.", call. = FALSE)

  weighted <- length(ps_covariates %||% character()) > 0L

  if (weighted) {
    message(
      "IPTW enabled: propensity scores and weights will be estimated ",
      "separately within each imputed dataset (estimand = ", estimand, ")."
    )
  }

  treatment_reference_final <- NULL
  treatment_comparator_final <- NULL
  if (!is.null(treatment)) {
    trt1 <- .ekb_binary_factor(d1[[treatment]], reference = treatment_reference, name = treatment)
    treatment_reference_final <- levels(trt1)[1]
    treatment_comparator_final <- levels(trt1)[2]
  }
  if (weighted && is.null(treatment)) {
    stop("treatment must be specified when ps_covariates are used.", call. = FALSE)
  }
  if (weighted && !treatment %in% covariates) {
    stop("treatment must also be included in covariates for an IPTW treatment-effect model.", call. = FALSE)
  }

  fits <- vector("list", mi$m)
  weighting <- if (weighted) vector("list", mi$m) else NULL
  diagnostics <- if (weighted) vector("list", mi$m) else NULL

  for (i in seq_len(mi$m)) {
    d <- mice::complete(mi, i)
    if (!is.null(treatment)) {
      d[[treatment]] <- .ekb_binary_factor(d[[treatment]], reference = treatment_reference, name = treatment)
    }

    w <- NULL
    if (weighted) {
      w <- suppressMessages(
        fit_weighting(
          d,
          treatment = treatment,
          ps_covariates = ps_covariates,
          method = weight_method,
          estimand = estimand,
          stabilize = stabilize,
          treatment_reference = treatment_reference,
          ...
        )
      )
      weighting[[i]] <- w
      diagnostics[[i]] <- diagnose_weighting(w)
    }

    fits[[i]] <- fit_cox(
      d,
      time = time,
      event = event,
      covariates = covariates,
      weights = if (weighted) w$weights else NULL,
      robust = robust
    )$fit
  }

  pooled <- mice::pool(mice::as.mira(fits))
  pooled_summary <- .ekb_pool_summary(pooled, conf.level)

  treatment_effect <- NULL
  if (!is.null(treatment) && treatment %in% covariates) {
    treatment_term <- .ekb_treatment_term(fits[[1]], treatment)
    treatment_effect <- pooled_summary[
      pooled_summary$term == treatment_term,
      ,
      drop = FALSE
    ]
  }

  structure(
    list(
      imputation = imputation,
      fits = fits,
      pooled = pooled,
      pooled_summary = .ekb_as_tibble(pooled_summary),
      treatment_effect = if (is.null(treatment_effect)) NULL else .ekb_as_tibble(treatment_effect),
      weighted = weighted,
      weighting = weighting,
      diagnostics = diagnostics,
      diagnostic_summary = .ekb_mi_diagnostic_summary(diagnostics),
      weighting_approach = if (weighted) "within-imputation" else "none",
      time = time,
      event = event,
      covariates = covariates,
      treatment = treatment,
      ps_covariates = ps_covariates,
      estimand = if (weighted) estimand else NA_character_,
      weight_method = if (weighted) weight_method else NA_character_,
      stabilize = if (weighted) stabilize else NA,
      treatment_reference = treatment_reference_final,
      treatment_comparator = treatment_comparator_final,
      call = match.call()
    ),
    class = "ekb_mi_cox"
  )
}
