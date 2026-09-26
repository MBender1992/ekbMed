# Publication-ready clinical tables
# gtsummary backend with sparse, journal-style flextable rendering.

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.ekb_table_label <- function(x, labels = NULL) {
  x <- as.character(x)
  if (is.null(labels)) return(x)
  if (is.null(names(labels))) {
    stop("labels must be a named character vector.", call. = FALSE)
  }
  ifelse(x %in% names(labels), unname(labels[x]), x)
}

.ekb_table_level_label <- function(variable, level, level_labels = NULL) {
  level <- as.character(level)
  if (is.null(level_labels)) return(level)

  mapping <- level_labels[[variable]]
  if (is.null(mapping)) return(level)
  if (is.null(names(mapping))) {
    stop(
      sprintf("level_labels[['%s']] must be a named character vector.", variable),
      call. = FALSE
    )
  }

  if (level %in% names(mapping)) {
    return(as.character(mapping[[level]]))
  }
  level
}

.ekb_gtsummary_labels <- function(labels = NULL, variables = NULL) {
  if (is.null(labels)) return(NULL)
  if (is.null(names(labels))) {
    stop("labels must be a named character vector.", call. = FALSE)
  }

  if (is.null(variables)) {
    variables <- names(labels)
  } else {
    variables <- intersect(variables, names(labels))
  }

  lapply(
    variables,
    function(v) rlang::new_formula(rlang::sym(v), as.character(labels[[v]]))
  )
}

.ekb_table_pvalue <- function(x, digits = 3) {
  x <- as.numeric(x)
  out <- rep("", length(x))
  ok <- !is.na(x)

  threshold <- 10^(-digits)
  out[ok & x < threshold] <- paste0("<", formatC(threshold, format = "f", digits = digits))
  out[ok & x >= threshold] <- formatC(x[ok & x >= threshold], format = "f", digits = digits)

  # Journal-style removal of unnecessary trailing zeroes.
  idx <- ok & x >= threshold
  out[idx] <- sub("0+$", "", out[idx])
  out[idx] <- sub("\\.$", "", out[idx])
  out
}

.ekb_table_ratio <- function(x, digits = 2) {
  x <- as.numeric(x)
  out <- rep("", length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

.ekb_table_smd <- function(x, digits = 2) {
  x <- abs(as.numeric(x))
  out <- rep("", length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

.ekb_table_group_levels <- function(x) {
  if (is.factor(x)) {
    lev <- levels(droplevels(x))
  } else {
    lev <- unique(as.character(stats::na.omit(x)))
  }
  lev
}

.ekb_table_group_headers <- function(tbl,
                                     data,
                                     by,
                                     by_labels = NULL,
                                     weights = NULL,
                                     weighted_digits = 0) {
  lev <- .ekb_table_group_levels(data[[by]])

  if (is.null(weights)) {
    counts <- vapply(
      lev,
      function(x) sum(!is.na(data[[by]]) & as.character(data[[by]]) == x),
      numeric(1)
    )
    n_label <- "N"
  } else {
    w <- if (inherits(weights, "ekb_weighting")) {
      weights$weights
    } else {
      .ekb_resolve_weights(data, weights)
    }
    w <- .ekb_resolve_weights(data, w)

    counts <- vapply(
      lev,
      function(x) {
        idx <- !is.na(data[[by]]) & as.character(data[[by]]) == x
        sum(w[idx], na.rm = TRUE)
      },
      numeric(1)
    )
    n_label <- "Weighted N"
  }

  display <- vapply(
    lev,
    function(x) .ekb_table_level_label(by, x, by_labels),
    character(1)
  )

  formatted_counts <- if (is.null(weights)) {
    format(round(counts), big.mark = ",", scientific = FALSE, trim = TRUE)
  } else {
    formatC(
      counts,
      format = "f",
      digits = weighted_digits,
      big.mark = ","
    )
  }

  headers <- setNames(
    as.list(
      sprintf(
        "**%s**  \n%s = %s",
        display,
        n_label,
        formatted_counts
      )
    ),
    paste0("stat_", seq_along(lev))
  )

  rlang::exec(gtsummary::modify_header, tbl, !!!headers)
}

.ekb_table_total_header <- function(data,
                                    weights = NULL,
                                    weighted_digits = 0,
                                    label = "Total") {
  if (is.null(weights)) {
    n <- nrow(data)
    n_text <- format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
    return(paste0("**", label, "**  \nN = ", n_text))
  }

  w <- if (inherits(weights, "ekb_weighting")) {
    weights$weights
  } else {
    .ekb_resolve_weights(data, weights)
  }
  w <- .ekb_resolve_weights(data, w)

  n_text <- formatC(
    sum(w, na.rm = TRUE),
    format = "f",
    digits = weighted_digits,
    big.mark = ","
  )

  paste0("**", label, "**  \nWeighted N = ", n_text)
}

.ekb_table_caption <- function(tbl, caption = NULL) {
  if (is.null(caption) || !nzchar(caption)) return(tbl)
  attr(tbl, "ekb_caption") <- as.character(caption)
  tbl
}


# -------------------------------------------------------------------------
# Baseline / descriptive tables
# -------------------------------------------------------------------------

#' Create a baseline characteristics table
#'
#' Creates a gtsummary baseline table with optional weighting, Total and balance diagnostics.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param by Name of the table grouping column.
#' @param include Optional character vector of variables to include.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @param labels Optional named character vector mapping variable names to display labels.
#' @param by_labels Optional named list keyed by the grouping variable, containing a named vector mapping its levels to display labels.
#' @param continuous Continuous summaries: `"median_iqr"` or `"mean_sd"`.
#' @param missing Missing-category display policy; see the function default and Details.
#' @param missing_text Label for missing or unknown values.
#' @param add_p Whether to add comparison p-values.
#' @param add_smd Whether to display absolute standardized mean differences instead of p-values.
#' @param overall Whether to append a Total column.
#' @param p_digits Number of decimal places for p-values.
#' @param smd_digits Number of decimal places for standardized mean differences.
#' @param caption Optional plain-text table caption.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return A gtsummary object with ekb_table_type and ekb_weighted attributes.
#' @details Default unweighted tables use p-values; default weighted tables use absolute SMDs, requiring smd. Both cannot be requested simultaneously. Weighted N is the sum of weights (pseudo-population), not unique patients. Missing weights are rejected. Additional arguments go to gtsummary::tbl_summary or tbl_svysummary. Missing-category options follow gtsummary.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("gtsummary", quietly = TRUE)) {
#'   tbl_baseline(d, "treatment", include = c("age", "sex"), add_p = FALSE)
#' }
tbl_baseline <- function(data,
                         by,
                         include = NULL,
                         weights = NULL,
                         labels = NULL,
                         by_labels = NULL,
                         continuous = c("median_iqr", "mean_sd"),
                         missing = "ifany",
                         missing_text = "Missing/unknown",
                         add_p = is.null(weights),
                         add_smd = !is.null(weights),
                         overall = TRUE,
                         p_digits = 3,
                         smd_digits = 2,
                         caption = NULL,
                         ...) {
  .ekb_require("gtsummary")
  .ekb_require("rlang")
  .ekb_require("dplyr")

  continuous <- match.arg(continuous)
  dat <- as.data.frame(data)
  .ekb_assert_columns(dat, c(by, include))
  # Keep backend group order aligned with the explicit header order.
  dat[[by]] <- factor(dat[[by]], levels = .ekb_table_group_levels(dat[[by]]))

  if (isTRUE(add_p) && isTRUE(add_smd)) {
    stop(
      "Choose either add_p = TRUE or add_smd = TRUE for a publication table, not both.",
      call. = FALSE
    )
  }

  by_sym <- rlang::sym(by)
  include_sel <- if (is.null(include)) {
    rlang::expr(dplyr::all_of(!!setdiff(names(dat), by)))
  } else {
    rlang::expr(dplyr::all_of(!!include))
  }

  variables <- include %||% setdiff(names(dat), by)
  label_spec <- .ekb_gtsummary_labels(labels, variables)

  continuous_stat <- switch(
    continuous,
    median_iqr = "{median} ({p25}, {p75})",
    mean_sd = "{mean} ({sd})"
  )

  statistic <- list(
    gtsummary::all_continuous() ~ continuous_stat,
    gtsummary::all_categorical() ~ "{n} ({p}%)"
  )

  if (is.null(weights)) {
    tbl <- gtsummary::tbl_summary(
      data = dat,
      by = !!by_sym,
      include = !!include_sel,
      label = label_spec,
      statistic = statistic,
      missing = missing,
      missing_text = missing_text,
      ...
    )
  } else {
    .ekb_require("survey")

    w <- if (inherits(weights, "ekb_weighting")) {
      weights$weights
    } else {
      .ekb_resolve_weights(dat, weights)
    }

    w <- .ekb_resolve_weights(dat, w)
    if (anyNA(w)) {
      stop("Table weights must be complete.", call. = FALSE)
    }

    dat$.ekb_weight <- w
    des <- survey::svydesign(
      ids = ~1,
      weights = ~.ekb_weight,
      data = dat
    )

    tbl <- gtsummary::tbl_svysummary(
      data = des,
      by = !!by_sym,
      include = !!include_sel,
      label = label_spec,
      statistic = statistic,
      missing = missing,
      missing_text = missing_text,
      ...
    )
  }

  if (!is.null(weights)) {
    if (inherits(weights, "ekb_weighting")) {
      message(
        "IPTW-weighted baseline table enabled (estimand = ",
        weights$estimand,
        "). Group and Total headers report weighted N."
      )
    } else {
      message(
        "Weighted baseline table enabled. Group and Total headers report weighted N."
      )
    }
  }

  tbl <- .ekb_table_group_headers(
    tbl,
    dat,
    by,
    by_labels,
    weights = weights
  )

  if (isTRUE(overall)) {
    tbl <- gtsummary::add_overall(
      tbl,
      last = TRUE,
      col_label = .ekb_table_total_header(
        dat,
        weights = weights,
        label = "Total"
      )
    )
  }

  if (isTRUE(add_p)) {
    tbl <- gtsummary::add_p(
      tbl,
      pvalue_fun = function(x) .ekb_table_pvalue(x, digits = p_digits)
    )
    tbl <- gtsummary::modify_header(tbl, p.value = "**P value**")
  }

  if (isTRUE(add_smd)) {
    .ekb_require("smd")

    tbl <- gtsummary::add_difference(
      tbl,
      test = list(
        gtsummary::all_continuous() ~ "smd",
        gtsummary::all_categorical() ~ "smd"
      ),
      estimate_fun = list(
        gtsummary::all_continuous() ~ function(x) .ekb_table_smd(x, smd_digits),
        gtsummary::all_categorical() ~ function(x) .ekb_table_smd(x, smd_digits)
      )
    )

    tbl <- gtsummary::modify_column_hide(
      tbl,
      columns = dplyr::any_of(c("conf.low", "conf.high", "p.value"))
    )
    tbl <- gtsummary::modify_header(tbl, estimate = "**|SMD|**")
  }

  tbl <- gtsummary::modify_header(tbl, label = "**Characteristic**")
  tbl <- gtsummary::bold_labels(tbl)

  if (!is.null(weights)) {
    tbl <- gtsummary::modify_source_note(
      tbl,
      "Weighted N is the sum of analysis weights and represents the IPTW pseudo-population, not the number of unique patients."
    )
  }

  tbl <- .ekb_table_caption(tbl, caption)

  attr(tbl, "ekb_table_type") <- "baseline"
  attr(tbl, "ekb_weighted") <- !is.null(weights)
  tbl
}


# -------------------------------------------------------------------------
# Cox regression tables
# Supports ordinary Cox and multiply-imputed Cox objects.
# -------------------------------------------------------------------------

#' Create a Cox regression table
#'
#' Accepts ekb_cox_fit, ekb_mi_cox or survival::coxph objects.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param labels Optional named character vector mapping variable names to display labels.
#' @param include Optional character vector of variables to include.
#' @param conf.level Confidence level between zero and one.
#' @param hr_digits Number of decimal places for hazard ratios and confidence limits.
#' @param p_digits Number of decimal places for p-values.
#' @param caption Optional plain-text table caption.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return A gtsummary regression table with weighting and imputation attributes.
#' @details Additional arguments go to gtsummary::tbl_regression. MI objects are converted to mice::mira from their stored fitted models without repeating imputation. gtsummary may require broom.helpers for regression formatting.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit <- fit_cox(d, "time", "event", c("treatment", "age"))
#' if (requireNamespace("gtsummary", quietly = TRUE) &&
#'     requireNamespace("broom.helpers", quietly = TRUE)) tbl_cox(fit)
tbl_cox <- function(object,
                    labels = NULL,
                    include = NULL,
                    conf.level = 0.95,
                    hr_digits = 2,
                    p_digits = 3,
                    caption = NULL,
                    ...) {
  .ekb_require("gtsummary")
  .ekb_require("rlang")

  if (inherits(object, "ekb_mi_cox")) {
    .ekb_require("mice")
    model <- mice::as.mira(object$fits)
    variables <- object$covariates
  } else if (inherits(object, "ekb_cox_fit")) {
    model <- object$fit
    variables <- object$covariates
  } else if (inherits(object, "coxph")) {
    model <- object
    variables <- all.vars(stats::delete.response(stats::terms(object)))
  } else {
    stop(
      "object must be an ekb_cox_fit, ekb_mi_cox, or survival::coxph object.",
      call. = FALSE
    )
  }

  label_spec <- .ekb_gtsummary_labels(labels, variables)

  args <- list(
    x = model,
    exponentiate = TRUE,
    conf.level = conf.level,
    label = label_spec,
    estimate_fun = function(x) .ekb_table_ratio(x, digits = hr_digits),
    pvalue_fun = function(x) .ekb_table_pvalue(x, digits = p_digits)
  )

  if (!is.null(include)) {
    args$include <- include
  }

  args <- c(args, list(...))
  tbl <- do.call(gtsummary::tbl_regression, args)

  tbl <- gtsummary::modify_header(
    tbl,
    label = "**Characteristic**",
    estimate = "**HR**",
    conf.low = paste0("**", round(100 * conf.level), "% CI**"),
    p.value = "**P value**"
  )

  tbl <- gtsummary::bold_labels(tbl)
  tbl <- .ekb_table_caption(tbl, caption)

  attr(tbl, "ekb_table_type") <- "cox"
  attr(tbl, "ekb_weighted") <- isTRUE(object$weighted)
  attr(tbl, "ekb_imputed") <- inherits(object, "ekb_mi_cox")
  tbl
}


# -------------------------------------------------------------------------
# Subgroup treatment-effect tables
# Supports ordinary and multiply-imputed subgroup objects.
# -------------------------------------------------------------------------

#' Create a subgroup treatment-effect table
#'
#' Accepts ekb_subgroup_cox or ekb_subgroup_mi_cox objects.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param labels Optional named character vector mapping variable names to display labels.
#' @param level_labels Optional named list of named character vectors mapping variable levels to display labels.
#' @param show_n Whether to show observation counts.
#' @param hr_digits Number of decimal places for hazard ratios and confidence limits.
#' @param p_digits Number of decimal places for p-values.
#' @param caption Optional plain-text table caption.
#' @return A gtsummary object with within-level treatment effects and an optional Overall row.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' sg <- subgroup_cox(d, "time", "event", "treatment", "sex")
#' if (requireNamespace("gtsummary", quietly = TRUE)) tbl_subgroup(sg)
tbl_subgroup <- function(object,
                         labels = NULL,
                         level_labels = NULL,
                         show_n = TRUE,
                         hr_digits = 2,
                         p_digits = 3,
                         caption = NULL) {
  .ekb_require("gtsummary")
  .ekb_require("dplyr")

  if (!inherits(object, c("ekb_subgroup_cox", "ekb_subgroup_mi_cox"))) {
    stop(
      "object must be an ekb_subgroup_cox or ekb_subgroup_mi_cox object.",
      call. = FALSE
    )
  }

  d <- as.data.frame(object$results)
  .ekb_assert_columns(
    d,
    c("subgroup", "level", "n", "hr", "conf.low", "conf.high", "p.value")
  )

  # Keep Overall first, followed by the original subgroup order.
  subgroup_order <- object$subgroups %||% unique(d$subgroup[d$subgroup != "Overall"])
  order_levels <- unique(c("Overall", subgroup_order, d$subgroup))
  d$.ekb_order <- match(d$subgroup, order_levels)
  d$.ekb_row <- seq_len(nrow(d))
  d <- d[order(d$.ekb_order, d$.ekb_row), , drop = FALSE]

  subgroup_display <- vapply(
    d$subgroup,
    function(x) {
      if (identical(x, "Overall")) return("Overall")
      .ekb_table_label(x, labels)
    },
    character(1)
  )

  first_in_block <- !duplicated(d$subgroup)
  subgroup_display[!first_in_block] <- ""

  category_display <- mapply(
    function(sg, lev) {
      if (identical(sg, "Overall")) return("")
      .ekb_table_level_label(sg, lev, level_labels)
    },
    d$subgroup,
    d$level,
    USE.NAMES = FALSE
  )

  estimable <- is.finite(d$hr) & is.finite(d$conf.low) & is.finite(d$conf.high)
  hr_ci <- rep("Not estimable", nrow(d))
  hr_ci[estimable] <- sprintf(
    paste0("%.", hr_digits, "f [%.", hr_digits, "f, %.", hr_digits, "f]"),
    d$hr[estimable],
    d$conf.low[estimable],
    d$conf.high[estimable]
  )

  body <- data.frame(
    subgroup = subgroup_display,
    category = category_display,
    n = as.integer(d$n),
    hr_ci = hr_ci,
    p_value = .ekb_table_pvalue(d$p.value, digits = p_digits),
    is_heading = first_in_block,
    stringsAsFactors = FALSE
  )

  if (!isTRUE(show_n)) {
    body$n <- NULL
  }

  tbl <- gtsummary::as_gtsummary(body)

  tbl <- gtsummary::modify_header(
    tbl,
    subgroup = "**Subgroup**",
    category = "**Category**",
    hr_ci = "**HR (95% CI)**",
    p_value = "**P value**"
  )

  if (isTRUE(show_n)) {
    tbl <- gtsummary::modify_header(tbl, n = "**N**")
    tbl <- gtsummary::modify_column_alignment(tbl, columns = n, align = "center")
  }

  tbl <- gtsummary::modify_column_hide(tbl, columns = is_heading)
  tbl <- gtsummary::modify_column_alignment(
    tbl,
    columns = c(hr_ci, p_value),
    align = "right"
  )

  tbl <- gtsummary::modify_bold(
    tbl,
    columns = subgroup,
    rows = is_heading
  )

  tbl <- .ekb_table_caption(tbl, caption)

  attr(tbl, "ekb_table_type") <- "subgroup"
  attr(tbl, "ekb_weighted") <- isTRUE(object$weighted)
  attr(tbl, "ekb_imputed") <- inherits(object, "ekb_subgroup_mi_cox")
  tbl
}


# -------------------------------------------------------------------------
# Survival summary table
# -------------------------------------------------------------------------

#' Format survival summary data
#'
#' Formats the output of summarize_survival as a rectangular table.
#'
#' @param summary_object The list returned by [summarize_survival()].
#' @param digits Number of decimal places for estimates.
#' @param percent_digits Number of decimal places for percentages.
#' @return A tibble; use flextable::flextable for Word output. This function does not return gtsummary.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' tbl_survival(summarize_survival(fit_survival(d, "time", "event")))
tbl_survival <- function(summary_object,
                         digits = 1,
                         percent_digits = 1) {
  .ekb_require("dplyr")
  .ekb_require("tidyr")

  if (!is.list(summary_object) || is.null(summary_object$median)) {
    stop("summary_object must be created by summarize_survival().", call. = FALSE)
  }

  med <- as.data.frame(summary_object$median)
  med$median_ci <- ifelse(
    is.na(med$median),
    "Not reached",
    sprintf(
      paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
      med$median,
      med$conf.low,
      med$conf.high
    )
  )

  out <- med[, c("group", "n", "events", "median_ci"), drop = FALSE]
  names(out)[names(out) == "median_ci"] <- sprintf(
    "Median (%g%% CI)",
    100 * summary_object$conf.level
  )

  lm <- as.data.frame(summary_object$landmarks)
  if (nrow(lm)) {
    lm$value <- sprintf(
      paste0(
        "%.", percent_digits, "f%% [%.", percent_digits,
        "f, %.", percent_digits, "f%%]"
      ),
      100 * lm$estimate,
      100 * lm$conf.low,
      100 * lm$conf.high
    )
    lm$label <- paste0(lm$time, "-month")
    wide <- tidyr::pivot_wider(
      lm[, c("group", "label", "value")],
      names_from = "label",
      values_from = "value"
    )
    out <- dplyr::left_join(out, wide, by = "group")
  }

  tibble::as_tibble(out)
}



# -------------------------------------------------------------------------
# Clinical outcome table
# Combines best response, dichotomous response rates, and optional
# Kaplan-Meier medians / landmark estimates in one publication-style table.
# -------------------------------------------------------------------------

.ekb_outcome_categorical_p <- function(x, group) {
  keep <- !is.na(x) & !is.na(group)
  x <- droplevels(factor(x[keep]))
  group <- droplevels(factor(group[keep]))

  if (nlevels(x) < 2L || nlevels(group) < 2L) {
    return(NA_real_)
  }

  tab <- table(x, group)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))

  if (any(chi$expected < 5)) {
    return(tryCatch(
      stats::fisher.test(tab)$p.value,
      error = function(e) chi$p.value
    ))
  }

  chi$p.value
}

.ekb_outcome_positive <- function(x, positive_level = "Yes") {
  if (is.logical(x)) return(x)
  if (is.numeric(x) && identical(positive_level, "Yes")) return(x == 1)
  as.character(x) == as.character(positive_level)
}

.ekb_outcome_format_n_pct <- function(n, denom, digits = 1) {
  if (is.na(n) || is.na(denom) || denom <= 0) return("")
  pct <- 100 * n / denom
  paste0(
    format(n, big.mark = ",", scientific = FALSE, trim = TRUE),
    " (",
    formatC(pct, format = "f", digits = digits),
    "%)"
  )
}

.ekb_outcome_format_survival <- function(est, low, high, digits = 1, percent = FALSE) {
  if (is.na(est)) return("Not reached")

  if (isTRUE(percent)) {
    return(sprintf(
      paste0("%.", digits, "f%% [%.", digits, "f, %.", digits, "f]"),
      100 * est, 100 * low, 100 * high
    ))
  }

  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    est, low, high
  )
}

.ekb_outcome_survival_spec <- function(x, name = NULL) {
  if (is.character(x)) {
    if (is.null(names(x)) || !all(c("time", "event") %in% names(x))) {
      stop(
        "Each survival specification must contain named 'time' and 'event' entries.",
        call. = FALSE
      )
    }
    return(list(
      time = unname(x[["time"]]),
      event = unname(x[["event"]]),
      label = name %||% unname(x[["time"]])
    ))
  }

  if (is.list(x)) {
    if (is.null(x$time) || is.null(x$event)) {
      stop(
        "Each survival specification must contain 'time' and 'event'.",
        call. = FALSE
      )
    }
    return(list(
      time = x$time,
      event = x$event,
      label = x$label %||% name %||% x$time
    ))
  }

  stop("Invalid survival specification.", call. = FALSE)
}

#' Create a response and survival outcome table
#'
#' Combines observed response descriptions with optional weighted or unweighted survival summaries.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param by Name of the table grouping column.
#' @param best_response Optional name of the best-overall-response column. Factor levels determine order.
#' @param response_rates Optional named character vector of binary response-column names, for example `c(ORR = "orr", DCR = "dcr")`. Create these columns explicitly upstream.
#' @param survival Optional named list of endpoints, each containing `time` and `event`, for example `list(OS = c(time = "months", event = "death"))`.
#' @param landmarks Numeric follow-up times in the same unit as the time column. NULL or `numeric()` requests no landmarks.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @param labels Optional named character vector mapping variable names to display labels.
#' @param by_labels Optional named list keyed by the grouping variable, containing a named vector mapping its levels to display labels.
#' @param total Whether to include the Total column.
#' @param positive_level Positive category for binary response columns. Logical and numeric columns are also supported.
#' @param add_p Whether to add comparison p-values.
#' @param percent_digits Number of decimal places for percentages.
#' @param survival_digits Number of decimal places for median survival.
#' @param p_digits Number of decimal places for p-values.
#' @param conf.level Confidence level between zero and one.
#' @param survival_heading Optional survival-section heading; defaults reflect weighting.
#' @param caption Optional plain-text table caption.
#' @return A gtsummary object with ekb_weighted_survival and ekb_survival_weighting attributes.
#' @details BOR and binary response columns remain observed unweighted descriptions even when survival is weighted. Binary response columns must be prepared upstream; they are not inferred from best_response. Response denominators exclude NA; nonmissing coded unknown values remain included unless converted explicitly. Group and Total headers show observed N. Survival estimates and survival log-rank tests use the supplied weights. Landmarks are optional and are labelled in months. The weighted pooled Total survival curve is descriptive, not a treatment contrast.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("gtsummary", quietly = TRUE)) {
#'   tbl_outcomes(d, "treatment", survival = list(OS = c(time = "time", event = "event")),
#'     landmarks = c(12, 24), add_p = FALSE)
#' }
tbl_outcomes <- function(data,
                         by,
                         best_response = NULL,
                         response_rates = NULL,
                         survival = NULL,
                         landmarks = NULL,
                         weights = NULL,
                         labels = NULL,
                         by_labels = NULL,
                         total = TRUE,
                         positive_level = "Yes",
                         add_p = TRUE,
                         percent_digits = 1,
                         survival_digits = 1,
                         p_digits = 3,
                         conf.level = 0.95,
                         survival_heading = NULL,
                         caption = NULL) {
  .ekb_require("gtsummary")

  dat <- as.data.frame(data)

  survival_requested <- !is.null(survival) && length(survival) > 0L
  survival_weighted <- survival_requested && !is.null(weights)
  survival_weights <- NULL
  weighting_label <- NULL

  if (survival_weighted) {
    if (inherits(weights, "ekb_weighting")) {
      survival_weights <- weights$weights
      weighting_label <- paste0("IPTW-adjusted (", weights$estimand, ")")
      message(
        "IPTW-adjusted survival estimates enabled (estimand = ",
        weights$estimand,
        "). Categorical response outcomes remain unweighted descriptive summaries."
      )
    } else {
      survival_weights <- .ekb_resolve_weights(dat, weights)
      weighting_label <- "weighted"
      message(
        "Weighted survival estimates enabled. ",
        "Categorical response outcomes remain unweighted descriptive summaries."
      )
    }

    survival_weights <- .ekb_resolve_weights(dat, survival_weights)
  }

  if (is.null(survival_heading)) {
    survival_heading <- if (survival_weighted) {
      if (inherits(weights, "ekb_weighting")) {
        "Survival (IPTW-adjusted)"
      } else {
        "Survival (weighted)"
      }
    } else {
      "Survival (unadjusted)"
    }
  }

  response_vars <- unname(response_rates %||% character())
  survival_vars <- character()

  if (!is.null(survival)) {
    if (!is.list(survival)) {
      stop("survival must be a named list of time/event specifications.", call. = FALSE)
    }
    if (is.null(names(survival)) || any(!nzchar(names(survival)))) {
      stop("survival must be a named list.", call. = FALSE)
    }

    survival_specs <- lapply(
      seq_along(survival),
      function(i) .ekb_outcome_survival_spec(survival[[i]], names(survival)[i])
    )
    survival_vars <- unique(unlist(lapply(
      survival_specs,
      function(x) c(x$time, x$event)
    )))
  } else {
    survival_specs <- list()
  }

  .ekb_assert_columns(
    dat,
    c(by, best_response, response_vars, survival_vars)
  )

  group_levels <- .ekb_table_group_levels(dat[[by]])
  if (length(group_levels) < 2L) {
    stop("by must contain at least two observed groups.", call. = FALSE)
  }

  rows <- list()
  row_index <- 1L

  add_row <- function(variable,
                      row_type,
                      label,
                      values = rep("", length(group_levels)),
                      total_value = "",
                      p.value = "",
                      is_heading = FALSE) {
    z <- data.frame(
      variable = variable,
      row_type = row_type,
      label = label,
      p.value = p.value,
      is_heading = is_heading,
      stringsAsFactors = FALSE
    )
    for (j in seq_along(group_levels)) {
      z[[paste0("stat_", j)]] <- values[j]
    }
    if (isTRUE(total)) z$stat_total <- total_value
    z
  }

  # Best overall response -------------------------------------------------
  if (!is.null(best_response)) {
    x <- dat[[best_response]]
    response_label <- .ekb_table_label(best_response, labels)

    p <- if (isTRUE(add_p)) {
      .ekb_outcome_categorical_p(x, dat[[by]])
    } else {
      NA_real_
    }

    rows[[row_index]] <- add_row(
      variable = best_response,
      row_type = "label",
      label = response_label,
      p.value = if (isTRUE(add_p)) .ekb_table_pvalue(p, p_digits) else "",
      is_heading = TRUE
    )
    row_index <- row_index + 1L

    levels_response <- if (is.factor(x)) {
      levels(droplevels(x))
    } else {
      unique(as.character(stats::na.omit(x)))
    }

    for (lev in levels_response) {
      vals <- vapply(
        group_levels,
        function(g) {
          idx <- !is.na(dat[[by]]) &
            as.character(dat[[by]]) == g &
            !is.na(x)
          denom <- sum(idx)
          n <- sum(idx & as.character(x) == lev)
          .ekb_outcome_format_n_pct(n, denom, percent_digits)
        },
        character(1)
      )

      idx_total <- !is.na(dat[[by]]) & !is.na(x)
      total_value <- .ekb_outcome_format_n_pct(
        sum(idx_total & as.character(x) == lev),
        sum(idx_total),
        percent_digits
      )

      rows[[row_index]] <- add_row(
        variable = best_response,
        row_type = "level",
        label = lev,
        values = vals,
        total_value = total_value
      )
      row_index <- row_index + 1L
    }
  }

  # ORR / DCR / other dichotomous response rates ------------------------
  if (length(response_vars)) {
    rate_names <- names(response_rates)
    if (is.null(rate_names)) rate_names <- rep("", length(response_vars))

    for (i in seq_along(response_vars)) {
      v <- response_vars[i]
      x <- dat[[v]]
      positive <- .ekb_outcome_positive(x, positive_level)

      display_label <- if (nzchar(rate_names[i])) {
        rate_names[i]
      } else {
        .ekb_table_label(v, labels)
      }

      vals <- vapply(
        group_levels,
        function(g) {
          idx <- !is.na(dat[[by]]) &
            as.character(dat[[by]]) == g &
            !is.na(x)
          denom <- sum(idx)
          n <- sum(idx & positive, na.rm = TRUE)
          .ekb_outcome_format_n_pct(n, denom, percent_digits)
        },
        character(1)
      )

      p <- if (isTRUE(add_p)) {
        .ekb_outcome_categorical_p(x, dat[[by]])
      } else {
        NA_real_
      }

      idx_total <- !is.na(dat[[by]]) & !is.na(x)
      total_value <- .ekb_outcome_format_n_pct(
        sum(idx_total & positive, na.rm = TRUE),
        sum(idx_total),
        percent_digits
      )

      rows[[row_index]] <- add_row(
        variable = v,
        row_type = "label",
        label = display_label,
        values = vals,
        total_value = total_value,
        p.value = if (isTRUE(add_p)) .ekb_table_pvalue(p, p_digits) else "",
        is_heading = TRUE
      )
      row_index <- row_index + 1L
    }
  }

  # Survival --------------------------------------------------------------
  if (length(survival_specs)) {
    if (!is.null(survival_heading) && nzchar(survival_heading)) {
      rows[[row_index]] <- add_row(
        variable = ".survival_heading",
        row_type = "label",
        label = survival_heading,
        is_heading = TRUE
      )
      row_index <- row_index + 1L
    }

    landmarks_use <- if (is.null(landmarks)) numeric() else sort(unique(as.numeric(landmarks)))

    for (spec in survival_specs) {
      fit <- fit_survival(
        dat,
        time = spec$time,
        event = spec$event,
        group = by,
        weights = survival_weights,
        conf.level = conf.level
      )

      sm <- summarize_survival(
        fit,
        landmarks = landmarks_use,
        conf.level = conf.level,
        extend = FALSE
      )

      p <- if (isTRUE(add_p)) {
        survival_logrank(
          dat,
          time = spec$time,
          event = spec$event,
          group = by,
          weights = survival_weights
        )$p.value
      } else {
        NA_real_
      }

      med <- as.data.frame(sm$median)
      med_vals <- vapply(
        group_levels,
        function(g) {
          z <- med[as.character(med$group) == g, , drop = FALSE]
          if (!nrow(z)) return("")
          .ekb_outcome_format_survival(
            z$median[1],
            z$conf.low[1],
            z$conf.high[1],
            digits = survival_digits,
            percent = FALSE
          )
        },
        character(1)
      )

      total_med_value <- ""
      total_lm <- data.frame()
      if (isTRUE(total)) {
        fit_total <- fit_survival(
          dat,
          time = spec$time,
          event = spec$event,
          group = NULL,
          weights = survival_weights,
          conf.level = conf.level
        )
        sm_total <- summarize_survival(
          fit_total,
          landmarks = landmarks_use,
          conf.level = conf.level,
          extend = FALSE
        )
        med_total <- as.data.frame(sm_total$median)
        if (nrow(med_total)) {
          total_med_value <- .ekb_outcome_format_survival(
            med_total$median[1],
            med_total$conf.low[1],
            med_total$conf.high[1],
            digits = survival_digits,
            percent = FALSE
          )
        }
        total_lm <- as.data.frame(sm_total$landmarks)
      }

      rows[[row_index]] <- add_row(
        variable = spec$time,
        row_type = "level",
        label = paste0("Median ", spec$label, ", months"),
        values = med_vals,
        total_value = total_med_value,
        p.value = if (isTRUE(add_p)) .ekb_table_pvalue(p, p_digits) else ""
      )
      row_index <- row_index + 1L

      if (length(landmarks_use)) {
        lm <- as.data.frame(sm$landmarks)

        for (tt in landmarks_use) {
          lm_vals <- vapply(
            group_levels,
            function(g) {
              z <- lm[
                as.character(lm$group) == g &
                  abs(as.numeric(lm$time) - tt) < 1e-8,
                ,
                drop = FALSE
              ]
              if (!nrow(z)) return("")
              .ekb_outcome_format_survival(
                z$estimate[1],
                z$conf.low[1],
                z$conf.high[1],
                digits = survival_digits,
                percent = TRUE
              )
            },
            character(1)
          )

          total_value <- ""
          if (isTRUE(total) && nrow(total_lm)) {
            z_total <- total_lm[
              abs(as.numeric(total_lm$time) - tt) < 1e-8,
              ,
              drop = FALSE
            ]
            if (nrow(z_total)) {
              total_value <- .ekb_outcome_format_survival(
                z_total$estimate[1],
                z_total$conf.low[1],
                z_total$conf.high[1],
                digits = survival_digits,
                percent = TRUE
              )
            }
          }

          rows[[row_index]] <- add_row(
            variable = spec$time,
            row_type = "landmark",
            label = paste0(format(tt, trim = TRUE, scientific = FALSE), "-month ", spec$label),
            values = lm_vals,
            total_value = total_value
          )
          row_index <- row_index + 1L
        }
      }
    }
  }

  if (!length(rows)) {
    stop(
      "Specify at least one of best_response, response_rates, or survival.",
      call. = FALSE
    )
  }

  body <- do.call(rbind, rows)

  # Put displayed columns into journal-style order.
  stat_cols <- paste0("stat_", seq_along(group_levels))
  display_stat_cols <- c(
    stat_cols,
    if (isTRUE(total)) "stat_total" else character()
  )
  body <- body[, c(
    "variable", "row_type", "label",
    display_stat_cols,
    "p.value", "is_heading"
  ), drop = FALSE]

  tbl <- gtsummary::as_gtsummary(body)

  tbl <- gtsummary::modify_column_hide(
    tbl,
    columns = dplyr::any_of(c("variable", "row_type", "is_heading"))
  )

  tbl <- gtsummary::modify_header(
    tbl,
    label = "**Characteristic**",
    p.value = if (isTRUE(add_p)) "**P value**" else ""
  )

  tbl <- .ekb_table_group_headers(tbl, dat, by, by_labels)

  if (isTRUE(total)) {
    tbl <- gtsummary::modify_header(
      tbl,
      stat_total = .ekb_table_total_header(dat, label = "Total")
    )
  }

  tbl <- gtsummary::modify_column_alignment(
    tbl,
    columns = dplyr::all_of(display_stat_cols),
    align = "center"
  )

  if (isTRUE(add_p)) {
    tbl <- gtsummary::modify_column_alignment(
      tbl,
      columns = p.value,
      align = "right"
    )
  } else {
    tbl <- gtsummary::modify_column_hide(tbl, columns = p.value)
  }

  tbl <- gtsummary::modify_bold(
    tbl,
    columns = label,
    rows = is_heading
  )

  survival_note <- if (survival_weighted) {
    if (inherits(weights, "ekb_weighting")) {
      paste0(
        "Categorical outcomes are shown as observed n (%). Survival results are IPTW-adjusted Kaplan-Meier estimates with 95% CIs (",
        weights$estimand,
        " estimand)."
      )
    } else {
      "Categorical outcomes are shown as observed n (%). Survival results are weighted Kaplan-Meier estimates with 95% CIs."
    }
  } else {
    "Categorical outcomes are shown as n (%). Survival results are unweighted Kaplan-Meier estimates with 95% CIs."
  }

  survival_note <- sub("95%", paste0(100 * conf.level, "%"), survival_note, fixed = TRUE)
  tbl <- gtsummary::modify_source_note(tbl, survival_note)

  if (isTRUE(add_p)) {
    survival_test_note <- if (survival_weighted) {
      if (inherits(weights, "ekb_weighting")) {
        "the IPTW-adjusted log-rank test"
      } else {
        "the weighted log-rank test"
      }
    } else {
      "the log-rank test"
    }

    tbl <- gtsummary::modify_source_note(
      tbl,
      paste0(
        "P values use Pearson's chi-squared test or Fisher's exact test for categorical outcomes, as appropriate, and ",
        survival_test_note,
        " for survival outcomes."
      )
    )
  }

  tbl <- .ekb_table_caption(tbl, caption)

  attr(tbl, "ekb_table_type") <- "outcomes"
  attr(tbl, "ekb_group_levels") <- group_levels
  attr(tbl, "ekb_weighted_survival") <- survival_weighted
  attr(tbl, "ekb_survival_weighting") <- weighting_label
  tbl
}



# -------------------------------------------------------------------------
# Treatment course / discontinuation / toxicity table
# -------------------------------------------------------------------------

.ekb_clean_table_category <- function(x) {
  x_chr <- as.character(x)
  blank <- !is.na(x_chr) & trimws(x_chr) == ""
  x_chr[blank] <- NA_character_
  x_chr
}

.ekb_course_levels <- function(x,
                               variable = NULL,
                               level_order = NULL,
                               missing = c("ifany", "no"),
                               missing_text = "Missing/unknown") {
  missing <- match.arg(missing)
  x_chr <- .ekb_clean_table_category(x)

  lev <- unique(stats::na.omit(x_chr))

  explicit_order <- NULL
  if (!is.null(level_order) && !is.null(variable)) {
    explicit_order <- level_order[[variable]]
  }

  if (!is.null(explicit_order)) {
    explicit_order <- as.character(explicit_order)
    lev <- c(
      explicit_order[explicit_order %in% lev],
      setdiff(lev, explicit_order)
    )
  } else if (length(lev) && all(grepl("^[0-9]+(?:\\\\.[0-9]+)?$", lev))) {
    lev <- lev[order(as.numeric(lev))]
  } else if (is.factor(x)) {
    factor_order <- levels(droplevels(x))
    factor_order <- factor_order[factor_order %in% lev]
    lev <- c(factor_order, setdiff(lev, factor_order))
  }

  if (missing == "ifany" && anyNA(x_chr) && !missing_text %in% lev) {
    lev <- c(lev, missing_text)
  }

  lev
}

.ekb_course_value <- function(x, level, missing_text = "Missing/unknown") {
  x_chr <- .ekb_clean_table_category(x)
  if (identical(level, missing_text)) {
    is.na(x_chr)
  } else {
    !is.na(x_chr) & x_chr == as.character(level)
  }
}

#' Describe treatment discontinuation and toxicity
#'
#' Creates ordered categorical summaries, treating blank strings as missing.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param by Name of the table grouping column.
#' @param end_reason Optional end-of-treatment reason column name.
#' @param toxicity Optional toxicity yes/no column name.
#' @param toxicity_grade Optional toxicity grade column name.
#' @param toxicity_type Optional toxicity type column name.
#' @param labels Optional named character vector mapping variable names to display labels.
#' @param by_labels Optional named list keyed by the grouping variable, containing a named vector mapping its levels to display labels.
#' @param total Whether to include the Total column.
#' @param level_order Optional named list mapping column names to ordered vectors of category labels.
#' @param toxicity_positive_level Category identifying patients with toxicity.
#' @param toxicity_details Use `"toxicity_positive"` to describe grades/types among patients with toxicity, or `"all"` to use the full cohort.
#' @param missing Missing-category display policy; see the function default and Details.
#' @param missing_text Label for missing or unknown values.
#' @param add_p Whether to add comparison p-values.
#' @param percent_digits Number of decimal places for percentages.
#' @param p_digits Number of decimal places for p-values.
#' @param caption Optional plain-text table caption.
#' @return A gtsummary object with treatment_course metadata and an optional Total column.
#' @details Blank strings are treated as missing. Category order follows level_order, then factor levels or natural numeric ordering. Toxicity grades/types default to the documented-toxicity subgroup. Percentages use the corresponding group size, including missing-category rows; missing = "no" hides the row without changing that denominator.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' d$reason <- rep(c("Completed", "Toxicity"), 40)
#' if (requireNamespace("gtsummary", quietly = TRUE)) {
#'   tbl_treatment_course(d, "treatment", end_reason = "reason", add_p = FALSE)
#' }
tbl_treatment_course <- function(data,
                                 by,
                                 end_reason = NULL,
                                 toxicity = NULL,
                                 toxicity_grade = NULL,
                                 toxicity_type = NULL,
                                 labels = NULL,
                                 by_labels = NULL,
                                 total = TRUE,
                                 level_order = NULL,
                                 toxicity_positive_level = "Yes",
                                 toxicity_details = c("toxicity_positive", "all"),
                                 missing = c("ifany", "no"),
                                 missing_text = "Missing/unknown",
                                 add_p = TRUE,
                                 percent_digits = 1,
                                 p_digits = 3,
                                 caption = "Treatment discontinuation and toxicity") {
  .ekb_require("gtsummary")

  toxicity_details <- match.arg(toxicity_details)
  missing <- match.arg(missing)

  dat <- as.data.frame(data)
  .ekb_assert_columns(
    dat,
    c(by, end_reason, toxicity, toxicity_grade, toxicity_type)
  )

  group_levels <- .ekb_table_group_levels(dat[[by]])
  if (length(group_levels) < 2L) {
    stop("by must contain at least two observed groups.", call. = FALSE)
  }

  rows <- list()
  row_index <- 1L

  add_row <- function(variable,
                      row_type,
                      label,
                      values = rep("", length(group_levels)),
                      total_value = "",
                      p.value = "",
                      is_heading = FALSE) {
    z <- data.frame(
      variable = variable,
      row_type = row_type,
      label = label,
      p.value = p.value,
      is_heading = is_heading,
      stringsAsFactors = FALSE
    )

    for (j in seq_along(group_levels)) {
      z[[paste0("stat_", j)]] <- values[j]
    }
    if (isTRUE(total)) z$stat_total <- total_value
    z
  }

  add_categorical_block <- function(variable,
                                    display_label,
                                    analysis_data = dat,
                                    row_type = "level") {
    x_raw <- analysis_data[[variable]]
    x <- .ekb_clean_table_category(x_raw)
    group <- analysis_data[[by]]

    p <- if (isTRUE(add_p)) {
      .ekb_outcome_categorical_p(x, group)
    } else {
      NA_real_
    }

    rows[[row_index]] <<- add_row(
      variable = variable,
      row_type = "label",
      label = display_label,
      p.value = if (isTRUE(add_p)) .ekb_table_pvalue(p, p_digits) else "",
      is_heading = TRUE
    )
    row_index <<- row_index + 1L

    levs <- .ekb_course_levels(
      x_raw,
      variable = variable,
      level_order = level_order,
      missing = missing,
      missing_text = missing_text
    )

    for (lev in levs) {
      vals <- vapply(
        group_levels,
        function(g) {
          idx_group <- !is.na(group) & as.character(group) == g
          denom <- sum(idx_group)

          if (!denom) return("")

          hit <- .ekb_course_value(
            x,
            lev,
            missing_text = missing_text
          )

          n <- sum(idx_group & hit, na.rm = TRUE)
          .ekb_outcome_format_n_pct(n, denom, percent_digits)
        },
        character(1)
      )

      idx_total <- !is.na(group)
      hit_total <- .ekb_course_value(
        x,
        lev,
        missing_text = missing_text
      )
      total_value <- .ekb_outcome_format_n_pct(
        sum(idx_total & hit_total, na.rm = TRUE),
        sum(idx_total),
        percent_digits
      )

      rows[[row_index]] <<- add_row(
        variable = variable,
        row_type = row_type,
        label = as.character(lev),
        values = vals,
        total_value = total_value
      )
      row_index <<- row_index + 1L
    }
  }

  # End-of-treatment reasons ------------------------------------------------
  if (!is.null(end_reason)) {
    add_categorical_block(
      end_reason,
      .ekb_table_label(end_reason, labels)
    )
  }

  # Any treatment-related toxicity -----------------------------------------
  if (!is.null(toxicity)) {
    add_categorical_block(
      toxicity,
      .ekb_table_label(toxicity, labels)
    )
  }

  # Toxicity details are usually most interpretable among patients with
  # documented toxicity. The full-cohort option is available if desired.
  detail_data <- dat
  detail_suffix <- ""

  if (
    toxicity_details == "toxicity_positive" &&
      !is.null(toxicity)
  ) {
    keep <- !is.na(dat[[toxicity]]) &
      as.character(dat[[toxicity]]) == as.character(toxicity_positive_level)

    detail_data <- dat[keep, , drop = FALSE]
    detail_suffix <- " (among patients with toxicity)"
  }

  if (!is.null(toxicity_grade)) {
    add_categorical_block(
      toxicity_grade,
      paste0(.ekb_table_label(toxicity_grade, labels), detail_suffix),
      analysis_data = detail_data,
      row_type = "detail"
    )
  }

  if (!is.null(toxicity_type)) {
    add_categorical_block(
      toxicity_type,
      paste0(.ekb_table_label(toxicity_type, labels), detail_suffix),
      analysis_data = detail_data,
      row_type = "detail"
    )
  }

  if (!length(rows)) {
    stop(
      "Specify at least one of end_reason, toxicity, toxicity_grade, or toxicity_type.",
      call. = FALSE
    )
  }

  body <- do.call(rbind, rows)
  stat_cols <- paste0("stat_", seq_along(group_levels))
  display_stat_cols <- c(
    stat_cols,
    if (isTRUE(total)) "stat_total" else character()
  )

  body <- body[, c(
    "variable", "row_type", "label",
    display_stat_cols, "p.value", "is_heading"
  ), drop = FALSE]

  tbl <- gtsummary::as_gtsummary(body)

  tbl <- gtsummary::modify_column_hide(
    tbl,
    columns = dplyr::any_of(c("variable", "row_type", "is_heading"))
  )

  tbl <- gtsummary::modify_header(
    tbl,
    label = "**Characteristic**",
    p.value = if (isTRUE(add_p)) "**P value**" else ""
  )

  tbl <- .ekb_table_group_headers(tbl, dat, by, by_labels)

  if (isTRUE(total)) {
    tbl <- gtsummary::modify_header(
      tbl,
      stat_total = .ekb_table_total_header(dat, label = "Total")
    )
  }

  tbl <- gtsummary::modify_column_alignment(
    tbl,
    columns = dplyr::all_of(display_stat_cols),
    align = "center"
  )

  if (isTRUE(add_p)) {
    tbl <- gtsummary::modify_column_alignment(
      tbl,
      columns = p.value,
      align = "right"
    )
  } else {
    tbl <- gtsummary::modify_column_hide(tbl, columns = p.value)
  }

  tbl <- gtsummary::modify_bold(
    tbl,
    columns = label,
    rows = is_heading
  )

  tbl <- gtsummary::modify_source_note(
    tbl,
    "Values are shown as n (%). Toxicity grade/type distributions are restricted to patients with documented toxicity unless toxicity_details = 'all'."
  )

  if (isTRUE(add_p)) {
    tbl <- gtsummary::modify_source_note(
      tbl,
      "P values use Pearson's chi-squared test or Fisher's exact test, as appropriate."
    )
  }

  tbl <- .ekb_table_caption(tbl, caption)

  attr(tbl, "ekb_table_type") <- "treatment_course"
  attr(tbl, "ekb_group_levels") <- group_levels
  tbl
}

# -------------------------------------------------------------------------
# Journal-style flextable renderer
# -------------------------------------------------------------------------

#' Render a clinical table for Word
#'
#' Converts a gtsummary object into an editable black-and-white journal-style flextable.
#'
#' @param tbl A gtsummary object.
#' @param font Font family for table text.
#' @param textsize Body text size in points.
#' @param headsize Header text size in points.
#' @param footsize Footnote text size in points.
#' @param autofit Whether to fit cell widths and heights to content.
#' @return An editable flextable object.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("gtsummary", quietly = TRUE) &&
#'     requireNamespace("flextable", quietly = TRUE) &&
#'     requireNamespace("officer", quietly = TRUE)) {
#'   as_ekb_flextable(tbl_baseline(d, "treatment", "age", add_p = FALSE))
#' }
as_ekb_flextable <- function(tbl,
                             font = "Arial",
                             textsize = 9,
                             headsize = 9,
                             footsize = 8,
                             autofit = TRUE) {
  .ekb_require("gtsummary")
  .ekb_require("flextable")
  .ekb_require("officer")

  if (!inherits(tbl, "gtsummary")) {
    stop("tbl must be a gtsummary object.", call. = FALSE)
  }

  ft <- gtsummary::as_flex_table(tbl)

  # Render the caption natively in flextable. This avoids literal Markdown
  # markers (e.g. **Title**) appearing in Word/RStudio output.
  caption <- attr(tbl, "ekb_caption", exact = TRUE)
  if (!is.null(caption) && nzchar(caption)) {
    caption_par <- flextable::as_paragraph(
      flextable::as_chunk(
        caption,
        props = officer::fp_text(
          font.family = font,
          font.size = headsize,
          bold = TRUE
        )
      )
    )
    ft <- flextable::set_caption(
      ft,
      caption = caption_par
    )
  }

  # Sparse black-and-white journal style: no vertical grid, strong outer
  # rules, light rule below the header, compact spacing.
  border_outer <- officer::fp_border(color = "black", width = 1.0)
  border_header <- officer::fp_border(color = "black", width = 0.6)

  ft <- flextable::border_remove(ft)
  ft <- flextable::font(ft, fontname = font, part = "all")
  ft <- flextable::fontsize(ft, size = textsize, part = "body")
  ft <- flextable::fontsize(ft, size = headsize, part = "header")
  ft <- flextable::fontsize(ft, size = footsize, part = "footer")
  ft <- flextable::bold(ft, bold = TRUE, part = "header")

  ft <- flextable::padding(
    ft,
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 2,
    padding.right = 2,
    part = "all"
  )

  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::align(ft, align = "center", part = "body")
  ft <- flextable::align(ft, j = 1, align = "left", part = "header")
  ft <- flextable::align(ft, j = 1, align = "left", part = "body")
  ft <- flextable::align(ft, align = "left", part = "footer")

  # Hierarchical indentation for publication-style outcome/subgroup tables.
  # Level rows are indented once; survival landmarks are indented one step
  # further so they visually sit underneath the corresponding median row.
  table_type <- attr(tbl, "ekb_table_type", exact = TRUE)
  if (isTRUE(table_type %in% c("outcomes", "subgroup", "treatment_course")) && "row_type" %in% names(tbl$table_body)) {
    level_rows <- which(tbl$table_body$row_type == "level")
    if (length(level_rows)) {
      ft <- flextable::padding(
        ft,
        i = level_rows,
        j = 1,
        padding.left = 12,
        part = "body"
      )
    }

    landmark_rows <- which(tbl$table_body$row_type == "landmark")
    if (length(landmark_rows)) {
      ft <- flextable::padding(
        ft,
        i = landmark_rows,
        j = 1,
        padding.left = 24,
        part = "body"
      )
    }

    detail_rows <- which(tbl$table_body$row_type == "detail")
    if (length(detail_rows)) {
      ft <- flextable::padding(
        ft,
        i = detail_rows,
        j = 1,
        padding.left = 12,
        part = "body"
      )
    }
  }

  ft <- flextable::hline_top(ft, border = border_outer, part = "header")
  ft <- flextable::hline_bottom(ft, border = border_header, part = "header")
  ft <- flextable::hline_bottom(ft, border = border_outer, part = "body")

  ft <- flextable::valign(ft, valign = "center", part = "all")
  ft <- flextable::set_table_properties(ft, layout = "autofit", width = 1)

  if (isTRUE(autofit)) {
    ft <- flextable::autofit(ft)
  }

  ft
}
