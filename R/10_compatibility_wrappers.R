# Compatibility wrappers for selected historical emR names.
# These are intentionally thin and should not contain independent statistical implementations.

#' Legacy ATE weighting wrapper
#'
#' Compatibility wrapper delegating to fit_weighting.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param vars Legacy character vector of covariate names.
#' @param prop.var Legacy treatment-column name.
#' @return A numeric vector of unstabilized ATE weights.
#' @family compatibility wrappers
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("WeightIt", quietly = TRUE)) ate_weights(d, "age", "treatment")
ate_weights <- function(data, vars, prop.var) {
  fit_weighting(data, treatment = prop.var, ps_covariates = setdiff(vars, prop.var), method = "glm", estimand = "ATE")$weights
}

#' Legacy survival landmark wrapper
#'
#' Compatibility wrapper delegating to fit_survival and summarize_survival.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param status Legacy event-column name.
#' @param var Legacy grouping or subgroup-column name.
#' @param times Numeric landmark times.
#' @return A landmark summary table with a formatted column.
#' @family compatibility wrappers
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' survival_time(d, "time", "event", "treatment", times = c(12, 24))
survival_time <- function(data, time, status, var = NULL, times) {
  fit <- fit_survival(data, time = time, event = status, group = var)
  sm <- summarize_survival(fit, landmarks = times, extend = FALSE)$landmarks
  if (!nrow(sm)) return(data.frame())
  sm$formatted <- sprintf("%.1f%% (%.1f%%\u2013%.1f%%)", 100 * sm$estimate, 100 * sm$conf.low, 100 * sm$conf.high)
  sm
}

#' Legacy median survival wrapper
#'
#' Compatibility wrapper for median survival and optional log-rank testing.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param status Legacy event-column name.
#' @param var Legacy grouping or subgroup-column name.
#' @param round Decimal places in formatted medians.
#' @param statistics Whether to attach the log-rank p-value.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @param conf.type Confidence-interval transformation passed to [survival::survfit()].
#' @return A median table with formatted estimates and optional p.value attribute.
#' @family compatibility wrappers
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' add_median_survival(d, "time", "event", "treatment")
add_median_survival <- function(data, time, status, var, round = 1, statistics = TRUE, weights = NULL, conf.type = "log-log") {
  fit <- fit_survival(data, time = time, event = status, group = var, weights = weights, conf.type = conf.type)
  sm <- summarize_survival(fit, landmarks = numeric())$median
  sm$formatted <- ifelse(
    is.na(sm$median),
    "Not reached",
    sprintf(paste0("%.", round, "f (%.", round, "f\u2013%.", round, "f)"), sm$median, sm$conf.low, sm$conf.high)
  )
  if (isTRUE(statistics)) {
    lr <- survival_logrank(data, time = time, event = status, group = var, weights = weights)
    attr(sm, "p.value") <- lr$p.value
  }
  sm
}

#' Legacy Cox output wrapper
#'
#' Compatibility wrapper fitting only prespecified full models.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param status Legacy event-column name.
#' @param vars Legacy character vector of covariate names.
#' @param fixed.var Additional outcome covariates always included.
#' @param output Return a formatted table or the underlying Cox fit.
#' @param modeltype Legacy model type. Only full models are implemented; other values warn and fit the full model.
#' @param p.thres Unused compatibility argument; no automatic model selection is performed.
#' @param niter Unused compatibility argument.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @return A formatted coefficient table or a survival::coxph fit.
#' @family compatibility wrappers
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' cox_output(d, "time", "event", c("treatment", "age"))
cox_output <- function(data, time, status, vars, fixed.var = NULL, output = c("table", "fit"), modeltype = "full", p.thres = 0.1, niter = 10, weights = NULL) {
  output <- match.arg(output)
  if (!identical(modeltype, "full")) {
    warning("Backward selection is retained only as legacy behavior and is not implemented in this refactor. Fitting the full model.", call. = FALSE)
  }
  fit <- fit_cox(data, time = time, event = status, covariates = unique(c(vars, fixed.var)), weights = weights)
  if (output == "fit") return(fit$fit)
  format_cox(fit)
}

#' Legacy subgroup Cox wrapper
#'
#' Historical name retained for reproducibility; this is a subgroup analysis, not a meta-analysis.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param status Legacy event-column name.
#' @param vars Legacy character vector of covariate names.
#' @param var Legacy grouping or subgroup-column name.
#' @param meta.group Legacy treatment-column name (despite its historical spelling).
#' @param univariate Whether to omit additional adjustment covariates.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @return The subgroup_cox results tibble.
#' @family compatibility wrappers
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' coxph_meta_analysis(d, "time", "event", "age", "sex", "treatment")
coxph_meta_analysis <- function(data, time, status, vars, var, meta.group, univariate = FALSE, weights = NULL) {
  covs <- if (isTRUE(univariate)) NULL else setdiff(vars, c(var, meta.group))
  subgroup_cox(
    data = data,
    time = time,
    event = status,
    treatment = meta.group,
    subgroups = var,
    covariates = covs,
    weights = weights,
    include_overall = TRUE,
    interaction = FALSE
  )$results
}

# Historical MI+IPTW entry point. New code should use:
# impute_clinical() -> fit_mi_cox(..., ps_covariates = ...)
#' Legacy combined imputation and IPTW wrapper
#'
#' Compatibility wrapper calling impute_clinical and fit_mi_cox; new analyses should use these steps explicitly.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param impute_vars Names of variables explicitly selected for imputation. Treatment, event and follow-up cannot be imputation targets.
#' @param ps_covariates Character vector of baseline propensity-score predictors. In MI functions, NULL disables weighting; supplying predictors enables within-imputation weighting.
#' @param auxiliary_vars Names of auxiliary imputation predictors.
#' @param outcome_covariates Additional covariates for the outcome model in the compatibility wrapper.
#' @param m Number of completed imputations.
#' @param maxit Number of chained-equation iterations.
#' @param seed Random seed supplied to mice for reproducibility.
#' @param estimand Target estimand passed to [WeightIt::weightit()], default `"ATE"`.
#' @param weight_method Propensity-score estimation method passed to WeightIt.
#' @param stabilize Whether to request stabilized weights from WeightIt.
#' @param survival_imputation Whether to include the event indicator and Nelson-Aalen cumulative hazard as imputation predictors.
#' @param robust NULL delegates variance selection to [survival::coxph()]; non-integer weights normally trigger robust SE. TRUE explicitly requests sandwich variance, including with integer weights; FALSE requests model-based variance.
#' @param conf.level Confidence level between zero and one.
#' @param treatment_reference Optional reference level. Otherwise the first existing factor level (or default factor ordering) is used.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return An ekb_mi_cox object.
#' @details Legacy only: warns and delegates to impute_clinical followed by fit_mi_cox. Additional arguments go to impute_clinical. The new API avoids repeating imputation across analyses.
#' @family compatibility wrappers
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' \dontrun{
#' # Legacy only; prefer impute_clinical() followed by fit_mi_cox().
#' fit_mi_iptw_cox(d, "time", "event", "treatment", "age", "age", m = 2)
#' }
fit_mi_iptw_cox <- function(data,
                            time,
                            event,
                            treatment,
                            impute_vars,
                            ps_covariates,
                            auxiliary_vars = NULL,
                            outcome_covariates = NULL,
                            m = 20,
                            maxit = 10,
                            seed = 1234,
                            estimand = "ATE",
                            weight_method = "glm",
                            stabilize = FALSE,
                            survival_imputation = TRUE,
                            robust = NULL,
                            conf.level = 0.95,
                            treatment_reference = NULL,
                            ...) {
  warning(
    "fit_mi_iptw_cox() is retained for compatibility. Prefer impute_clinical() followed by fit_mi_cox(..., ps_covariates = ...).",
    call. = FALSE
  )

  mi <- impute_clinical(
    data = data,
    impute_vars = impute_vars,
    predictor_vars = unique(c(ps_covariates, outcome_covariates)),
    auxiliary_vars = auxiliary_vars,
    treatment = treatment,
    time = time,
    event = event,
    survival_imputation = survival_imputation,
    m = m,
    maxit = maxit,
    seed = seed,
    ...
  )

  fit_mi_cox(
    imputation = mi,
    time = time,
    event = event,
    covariates = unique(c(treatment, setdiff(outcome_covariates %||% character(), treatment))),
    treatment = treatment,
    ps_covariates = ps_covariates,
    estimand = estimand,
    weight_method = weight_method,
    stabilize = stabilize,
    treatment_reference = treatment_reference,
    robust = robust,
    conf.level = conf.level
  )
}

