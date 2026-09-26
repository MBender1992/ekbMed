# Propensity-score weighting / IPTW

#' Estimate reusable propensity-score weights
#'
#' Fits a binary-treatment propensity model using WeightIt.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param ps_covariates Character vector of baseline propensity-score predictors. In MI functions, NULL disables weighting; supplying predictors enables within-imputation weighting.
#' @param method Method passed to the weighting or adjusted-survival backend; see Details.
#' @param estimand Target estimand passed to [WeightIt::weightit()], default `"ATE"`.
#' @param stabilize Whether to request stabilized weights from WeightIt.
#' @param treatment_reference Optional reference level. Otherwise the first existing factor level (or default factor ordering) is used.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return An `ekb_weighting` list containing the WeightIt object, aligned weights, propensity_score and model metadata.
#' @details Requires complete baseline PS predictors and exactly two observed treatment levels. With logistic ATE weighting, the score is the probability of the second treatment level. Default method is glm, default estimand is ATE, and no stabilization or truncation is automatic. Additional arguments are passed to WeightIt::weightit. Inspect diagnose_weighting before interpreting weighted outcomes.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("WeightIt", quietly = TRUE)) {
#'   w <- fit_weighting(d, "treatment", c("age", "sex"))
#'   fit_cox(d, "time", "event", "treatment", weights = w)
#' }
fit_weighting <- function(data,
                          treatment,
                          ps_covariates,
                          method = "glm",
                          estimand = "ATE",
                          stabilize = FALSE,
                          treatment_reference = NULL,
                          ...) {
  .ekb_require("WeightIt")
  .ekb_assert_columns(data, c(treatment, ps_covariates))
  if (!length(ps_covariates)) stop("ps_covariates must contain at least one baseline covariate.", call. = FALSE)

  dat <- as.data.frame(data)
  .ekb_complete_ps(dat, treatment, ps_covariates)
  dat[[treatment]] <- .ekb_binary_factor(dat[[treatment]], reference = treatment_reference, name = treatment)
  formula <- .ekb_model_formula(treatment, ps_covariates)
  wobj <- WeightIt::weightit(
    formula,
    data = dat,
    method = method,
    estimand = estimand,
    stabilize = stabilize,
    ...
  )

  ps <- NULL
  if (!is.null(wobj$ps)) ps <- as.numeric(wobj$ps)
  if (is.null(ps) && !is.null(wobj$distance)) ps <- as.numeric(wobj$distance)

  structure(
    list(
      object = wobj,
      data = dat,
      weights = as.numeric(wobj$weights),
      propensity_score = ps,
      treatment = treatment,
      ps_covariates = ps_covariates,
      method = method,
      estimand = estimand,
      stabilize = stabilize,
      call = match.call()
    ),
    class = "ekb_weighting"
  )
}

#' Calculate logistic ATE weights directly
#'
#' Provides an explicit logistic-regression reference implementation for unstabilized ATE weights.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param ps_covariates Character vector of baseline propensity-score predictors. In MI functions, NULL disables weighting; supplying predictors enables within-imputation weighting.
#' @param treatment_reference Optional reference level. Otherwise the first existing factor level (or default factor ordering) is used.
#' @return A list containing model, propensity_score and weights.
#' @details Uses logistic regression and the second treatment level as the positive outcome. Treated weights are 1/PS and reference weights are 1/(1-PS). Requires complete predictors and does not stabilize or truncate.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' manual_ate_weights(d, "treatment", c("age", "sex"))
manual_ate_weights <- function(data, treatment, ps_covariates, treatment_reference = NULL) {
  .ekb_assert_columns(data, c(treatment, ps_covariates))
  dat <- as.data.frame(data)
  .ekb_complete_ps(dat, treatment, ps_covariates)
  dat[[treatment]] <- .ekb_binary_factor(dat[[treatment]], reference = treatment_reference, name = treatment)
  y <- as.integer(dat[[treatment]] == levels(dat[[treatment]])[2])
  form <- .ekb_model_formula(treatment, ps_covariates)
  mod <- stats::glm(form, data = dat, family = stats::binomial())
  ps <- stats::predict(mod, newdata = dat, type = "response")
  w <- ifelse(y == 1L, 1 / ps, 1 / (1 - ps))
  list(model = mod, propensity_score = ps, weights = w)
}

#' Inspect balance overlap and effective sample size
#'
#' Diagnoses an ekb_weighting object using cobalt.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param smd_threshold Absolute standardized mean difference threshold for diagnostics.
#' @param ps_limits Lower and upper propensity-score thresholds used to flag limited overlap.
#' @param weight_probs Probabilities at which to report weight quantiles.
#' @return A list with balance, balance_object, ess, weight_quantiles, ps_summary, max_abs_smd and positivity_flag.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("WeightIt", quietly = TRUE) &&
#'     requireNamespace("cobalt", quietly = TRUE)) {
#'   diagnose_weighting(fit_weighting(d, "treatment", "age"))
#' }
diagnose_weighting <- function(object,
                               smd_threshold = 0.10,
                               ps_limits = c(0.01, 0.99),
                               weight_probs = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)) {
  .ekb_require("cobalt")
  if (!inherits(object, "ekb_weighting")) stop("object must be created by fit_weighting().", call. = FALSE)

  bt <- cobalt::bal.tab(
    object$object,
    un = TRUE,
    binary = "std",
    thresholds = c(m = smd_threshold),
    quick = FALSE
  )
  balance <- as.data.frame(bt$Balance)
  balance$variable <- rownames(balance)
  rownames(balance) <- NULL
  balance <- balance[, c("variable", setdiff(names(balance), "variable")), drop = FALSE]

  trt <- object$data[[object$treatment]]
  lev <- levels(trt)
  ess <- do.call(rbind, lapply(lev, function(g) {
    idx <- trt == g
    data.frame(group = g, n = sum(idx), ess = .ekb_weighted_ess(object$weights[idx]), stringsAsFactors = FALSE)
  }))

  q <- do.call(rbind, lapply(lev, function(g) {
    idx <- trt == g
    qq <- stats::quantile(object$weights[idx], probs = weight_probs, na.rm = TRUE, names = FALSE)
    data.frame(group = g, probability = weight_probs, weight = qq, stringsAsFactors = FALSE)
  }))

  ps_summary <- NULL
  positivity_flag <- NA
  if (!is.null(object$propensity_score)) {
    ps <- object$propensity_score
    ps_summary <- do.call(rbind, lapply(lev, function(g) {
      idx <- trt == g
      data.frame(
        group = g,
        min = min(ps[idx], na.rm = TRUE),
        q01 = stats::quantile(ps[idx], 0.01, na.rm = TRUE, names = FALSE),
        median = stats::median(ps[idx], na.rm = TRUE),
        q99 = stats::quantile(ps[idx], 0.99, na.rm = TRUE, names = FALSE),
        max = max(ps[idx], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    positivity_flag <- any(ps < ps_limits[1] | ps > ps_limits[2], na.rm = TRUE)
  }

  smd_cols <- grep("Diff.Adj", names(balance), value = TRUE)
  max_smd <- if (length(smd_cols)) max(abs(balance[[smd_cols[1]]]), na.rm = TRUE) else NA_real_

  list(
    balance = .ekb_as_tibble(balance),
    balance_object = bt,
    ess = .ekb_as_tibble(ess),
    weight_quantiles = .ekb_as_tibble(q),
    ps_summary = if (is.null(ps_summary)) NULL else .ekb_as_tibble(ps_summary),
    max_abs_smd = max_smd,
    smd_threshold = smd_threshold,
    positivity_flag = positivity_flag,
    ps_limits = ps_limits
  )
}

#' Plot covariate balance
#'
#' Produces a cobalt love plot for an ekb_weighting object.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param smd_threshold Absolute standardized mean difference threshold for diagnostics.
#' @param var.order Variable ordering passed to [cobalt::love.plot()].
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return A ggplot object produced by cobalt.
#' @details Additional arguments are passed to cobalt::love.plot.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("WeightIt", quietly = TRUE) &&
#'     requireNamespace("cobalt", quietly = TRUE)) {
#'   plot_balance(fit_weighting(d, "treatment", "age"))
#' }
plot_balance <- function(object, smd_threshold = 0.10, var.order = "unadjusted", ...) {
  .ekb_require("cobalt")
  if (!inherits(object, "ekb_weighting")) stop("object must be created by fit_weighting().", call. = FALSE)
  cobalt::love.plot(
    object$object,
    stats = "mean.diffs",
    abs = TRUE,
    binary = "std",
    thresholds = c(m = smd_threshold),
    var.order = var.order,
    line = TRUE,
    stars = "none",
    ...
  )
}

#' Cap weights at prespecified quantiles
#'
#' Explicitly caps numeric weights; fitting functions never apply truncation automatically.
#'
#' @param weights A positive numeric weight vector.
#' @param probs Two increasing quantile probabilities defining weight caps.
#' @return A list with original and capped weights, caps, probabilities and n_changed.
#' @export
#' @examples
#' truncate_weights(c(1, 1.5, 2, 4, 12), probs = c(0.05, 0.95))
truncate_weights <- function(weights, probs = c(0.01, 0.99)) {
  if (!is.numeric(weights)) stop("weights must be numeric.", call. = FALSE)
  if (length(probs) != 2L || any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("probs must be two increasing probabilities between 0 and 1.", call. = FALSE)
  }
  caps <- stats::quantile(weights, probs = probs, na.rm = TRUE, names = FALSE)
  truncated <- pmin(pmax(weights, caps[1]), caps[2])
  list(
    original = weights,
    weights = truncated,
    lower_cap = caps[1],
    upper_cap = caps[2],
    probs = probs,
    n_changed = sum(truncated != weights, na.rm = TRUE)
  )
}

#' Fit adjusted survival with adjustedCurves
#'
#' Optional adjustedCurves backend complementing fit_survival.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @param conf_int Whether to compute or display confidence intervals.
#' @param method Method passed to the weighting or adjusted-survival backend; see Details.
#' @param bootstrap Whether to bootstrap adjusted survival estimates.
#' @param n_boot Number of bootstrap replicates.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return An `ekb_iptw_survival` list with the backend fit and analysis metadata.
#' @details Uses adjustedCurves::adjustedsurv with default method iptw_km and supplied numeric weights. Additional arguments go to that function. No weights are estimated here.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' \dontrun{
#' w <- manual_ate_weights(d, "treatment", "age")$weights
#' fit_iptw_survival(d, "time", "event", "treatment", w, conf_int = FALSE)
#' }
fit_iptw_survival <- function(data,
                              time,
                              event,
                              treatment,
                              weights,
                              conf_int = TRUE,
                              method = "iptw_km",
                              bootstrap = FALSE,
                              n_boot = 500,
                              ...) {
  .ekb_require("adjustedCurves")
  .ekb_assert_columns(data, c(time, event, treatment))
  dat <- as.data.frame(data)
  dat[[event]] <- .ekb_binary_event(dat[[event]], event)
  dat[[treatment]] <- .ekb_binary_factor(dat[[treatment]], name = treatment)
  w <- if (inherits(weights, "ekb_weighting")) weights$weights else .ekb_resolve_weights(dat, weights)

  w <- .ekb_resolve_weights(dat, w)
  .ekb_validate_time(dat, time)
  keep <- stats::complete.cases(dat[, c(time, event, treatment), drop = FALSE]) & !is.na(w)
  dat <- droplevels(dat[keep, , drop = FALSE])
  w <- w[keep]
  dat[[treatment]] <- .ekb_binary_factor(dat[[treatment]], name = treatment)

  fit <- adjustedCurves::adjustedsurv(
    data = dat,
    variable = treatment,
    ev_time = time,
    event = event,
    method = method,
    treatment_model = w,
    conf_int = conf_int,
    bootstrap = bootstrap,
    n_boot = n_boot,
    ...
  )
  structure(
    list(fit = fit, data = dat, time = time, event = event, treatment = treatment, weights = w, method = method),
    class = "ekb_iptw_survival"
  )
}

#' Plot an adjustedCurves survival fit
#'
#' Plots an ekb_iptw_survival object through the adjustedCurves plot method.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param conf_int Whether to compute or display confidence intervals.
#' @param risk_table Whether to display a number-at-risk table.
#' @param risk_table_use_weights Whether the adjustedCurves risk table uses weights.
#' @param risk_table_stratify Whether to stratify the risk table by treatment.
#' @param xlab Horizontal axis label.
#' @param ylab Vertical axis label.
#' @param legend.title Optional legend title.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return A ggplot or composed plot returned by the backend.
#' @details Additional arguments go to the adjustedCurves plot method. Confidence bands require pammtools; risk tables require cowplot.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' \dontrun{
#' w <- manual_ate_weights(d, "treatment", "age")$weights
#' a <- fit_iptw_survival(d, "time", "event", "treatment", w, conf_int = FALSE)
#' plot_iptw_survival(a, conf_int = FALSE, risk_table = FALSE)
#' }
plot_iptw_survival <- function(object,
                               conf_int = TRUE,
                               risk_table = TRUE,
                               risk_table_use_weights = FALSE,
                               risk_table_stratify = TRUE,
                               xlab = "Time (months)",
                               ylab = "Adjusted survival probability",
                               legend.title = NULL,
                               ...) {
  if (!inherits(object, "ekb_iptw_survival")) stop("object must be created by fit_iptw_survival().", call. = FALSE)
  if (isTRUE(conf_int)) .ekb_require("pammtools")
  if (isTRUE(risk_table)) .ekb_require("cowplot")
  plot(
    object$fit,
    conf_int = conf_int,
    risk_table = risk_table,
    risk_table_use_weights = risk_table_use_weights,
    risk_table_stratify = risk_table_stratify,
    linetype = TRUE,
    color = TRUE,
    xlab = xlab,
    ylab = ylab,
    legend.title = legend.title %||% object$treatment,
    gg_theme = theme_ekbmed(),
    ...
  )
}
