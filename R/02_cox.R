# Cox proportional hazards regression

#' Fit a Cox proportional-hazards model
#'
#' Fits a prespecified outcome model with optional reusable weights.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param covariates Character vector of prespecified outcome-model covariates. Include treatment explicitly for treatment-effect models.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @param robust NULL delegates variance selection to [survival::coxph()]; non-integer weights normally trigger robust SE. TRUE explicitly requests sandwich variance, including with integer weights; FALSE requests model-based variance.
#' @param start Optional start-time column for counting-process intervals; every start must precede stop.
#' @param ties Tie-handling method passed to [survival::coxph()].
#' @param model Whether to retain the fitted model frame.
#' @param x Whether to retain the model matrix.
#' @param y Whether to retain the survival response in the fitted Cox object.
#' @return An `ekb_cox_fit` list with fit, data, covariates, weights and weighting metadata.
#' @details Factor coding and reference levels are retained from the input. No covariate selection is performed. Missing outcome-model values follow the current coxph/session na.action policy. Weights are resolved once and are aligned by row position, not by patient ID. A weighted model with additional outcome covariates is an additionally adjusted conditional model; it is not automatically doubly robust.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit <- fit_cox(d, "time", "event", c("treatment", "age"))
#' tidy_cox(fit)
fit_cox <- function(data,
                    time,
                    event,
                    covariates,
                    weights = NULL,
                    robust = NULL,
                    start = NULL,
                    ties = "efron",
                    model = TRUE,
                    x = TRUE,
                    y = TRUE) {
  .ekb_require("survival")
  .ekb_assert_columns(data, c(time, event, start, covariates))
  if (!length(covariates)) stop("At least one covariate is required.", call. = FALSE)

  if (inherits(weights, "ekb_weighting")) {
    message(
      "Using propensity-score weights from fit_weighting() ",
      "(estimand = ", weights$estimand, ")."
    )
  }

  dat <- as.data.frame(data)
  dat[[event]] <- .ekb_binary_event(dat[[event]], event)
  .ekb_validate_time(dat, time, start)
  w <- .ekb_resolve_weights(dat, weights)
  formula <- .ekb_surv_formula(time, event, covariates, start = start)

  args <- list(
    formula = formula,
    data = dat,
    ties = ties,
    model = model,
    x = x,
    y = y
  )
  if (!is.null(w)) args$weights <- w
  if (!is.null(robust)) args$robust <- robust

  fit <- do.call(survival::coxph, args)
  structure(
    list(
      fit = fit,
      data = dat,
      time = time,
      event = event,
      start = start,
      covariates = covariates,
      weights = w,
      weighting = if (inherits(weights, "ekb_weighting")) weights else NULL,
      weighted = !is.null(w),
      robust = robust,
      call = match.call()
    ),
    class = "ekb_cox_fit"
  )
}

#' Extract Cox effects and standard errors
#'
#' Accepts an ekb_cox_fit or survival::coxph object; uses robust SE when available.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param conf.level Confidence level between zero and one.
#' @return A tibble with terms, log hazard ratios, SE, HR, confidence limits, p-values, n and events.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit <- fit_cox(d, "time", "event", c("treatment", "age"))
#' tidy_cox(fit)
tidy_cox <- function(object, conf.level = 0.95) {
  .ekb_require("broom")
  fit <- if (inherits(object, "ekb_cox_fit")) object$fit else object
  if (!inherits(fit, "coxph")) stop("object must be an ekb_cox_fit or survival::coxph object.", call. = FALSE)

  z <- broom::tidy(fit, conf.int = TRUE, conf.level = conf.level, exponentiate = FALSE)
  if ("robust.se" %in% names(z)) z$std.error <- z$robust.se
  z$log_hr <- z$estimate
  z$hr <- exp(z$estimate)
  z$hr_conf.low <- exp(z$conf.low)
  z$hr_conf.high <- exp(z$conf.high)
  z$n <- fit$n
  z$events <- fit$nevent
  z[, c("term", "log_hr", "std.error", "statistic", "p.value", "hr", "hr_conf.low", "hr_conf.high", "n", "events")]
}

#' Extract the global Cox Wald test
#'
#' Accepts an ekb_cox_fit or survival::coxph object.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @return A data frame containing statistic, df and p.value.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit <- fit_cox(d, "time", "event", c("treatment", "age"))
#' cox_global_wald(fit)
cox_global_wald <- function(object) {
  fit <- if (inherits(object, "ekb_cox_fit")) object$fit else object
  if (!inherits(fit, "coxph")) stop("object must be an ekb_cox_fit or survival::coxph object.", call. = FALSE)
  w <- summary(fit)$waldtest
  data.frame(
    statistic = unname(w[1]),
    df = unname(w[2]),
    p.value = unname(w[3]),
    stringsAsFactors = FALSE
  )
}

#' Check the proportional-hazards assumption
#'
#' Applies survival::cox.zph to an ekb_cox_fit or survival::coxph object.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param transform Time transformation for [survival::cox.zph()].
#' @param global Whether to include the global proportional-hazards test.
#' @return A cox.zph object.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit <- fit_cox(d, "time", "event", c("treatment", "age"))
#' check_cox_ph(fit)
check_cox_ph <- function(object, transform = "km", global = TRUE) {
  .ekb_require("survival")
  fit <- if (inherits(object, "ekb_cox_fit")) object$fit else object
  if (!inherits(fit, "coxph")) stop("object must be an ekb_cox_fit or survival::coxph object.", call. = FALSE)
  survival::cox.zph(fit, transform = transform, global = global)
}

#' Format a Cox coefficient table
#'
#' Adds publication-oriented formatted values to tidy_cox output.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param digits Number of decimal places for estimates.
#' @param p_digits Number of decimal places for p-values.
#' @return A coefficient tibble with additional HR (95% CI) and p columns.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit <- fit_cox(d, "time", "event", c("treatment", "age"))
#' format_cox(fit)
format_cox <- function(object, digits = 2, p_digits = 3) {
  z <- tidy_cox(object)
  z$`HR (95% CI)` <- sprintf(
    paste0("%.", digits, "f (%.", digits, "f\u2013%.", digits, "f)"),
    z$hr, z$hr_conf.low, z$hr_conf.high
  )
  z$p <- .ekb_format_p(z$p.value, digits = p_digits)
  z
}
