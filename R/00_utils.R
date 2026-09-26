# Internal helpers for ekbMed scripts

`%||%` <- function(x, y) if (is.null(x)) y else x

.ekb_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is required for this function. Please install it first.", package),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.ekb_assert_columns <- function(data, columns) {
  columns <- unique(stats::na.omit(columns))
  missing_cols <- setdiff(columns, names(data))
  if (length(missing_cols)) {
    stop(
      sprintf("Missing column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.ekb_bt <- function(x) {
  paste0("`", gsub("`", "\\`", x, fixed = TRUE), "`")
}

.ekb_surv_formula <- function(time, event, covariates = NULL, start = NULL) {
  rhs <- if (length(covariates)) paste(.ekb_bt(covariates), collapse = " + ") else "1"
  response <- if (is.null(start)) {
    sprintf("survival::Surv(%s, %s)", .ekb_bt(time), .ekb_bt(event))
  } else {
    sprintf("survival::Surv(%s, %s, %s)", .ekb_bt(start), .ekb_bt(time), .ekb_bt(event))
  }
  stats::as.formula(paste(response, "~", rhs), env = parent.frame())
}

.ekb_model_formula <- function(response, covariates) {
  rhs <- if (length(covariates)) paste(.ekb_bt(covariates), collapse = " + ") else "1"
  stats::as.formula(paste(.ekb_bt(response), "~", rhs), env = parent.frame())
}

.ekb_resolve_weights <- function(data, weights = NULL) {
  if (is.null(weights)) return(NULL)

  if (inherits(weights, "ekb_weighting")) {
    out <- weights$weights
  } else if (is.character(weights) && length(weights) == 1L) {
    .ekb_assert_columns(data, weights)
    out <- data[[weights]]
  } else {
    out <- weights
  }

  if (!is.numeric(out) || length(out) != nrow(data)) {
    stop(
      "weights must be an ekb_weighting object, a numeric vector of nrow(data), or a single column name.",
      call. = FALSE
    )
  }
  if (any(!is.finite(out) & !is.na(out)) || any(out <= 0, na.rm = TRUE)) {
    stop("All non-missing weights must be finite and > 0.", call. = FALSE)
  }
  out
}

.ekb_binary_event <- function(x, name = "event") {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(as.integer(x))
  vals <- sort(unique(stats::na.omit(x)))
  if (!all(vals %in% c(0, 1))) {
    stop(sprintf("%s must be coded 0/1 (or logical).", name), call. = FALSE)
  }
  as.integer(x)
}

.ekb_binary_factor <- function(x, reference = NULL, name = "treatment") {
  f <- droplevels(factor(x))
  if (nlevels(f) != 2L) {
    stop(sprintf("%s must have exactly two observed levels.", name), call. = FALSE)
  }
  if (!is.null(reference)) {
    if (!reference %in% levels(f)) {
      stop(sprintf("Reference level '%s' not found in %s.", reference, name), call. = FALSE)
    }
    f <- stats::relevel(f, ref = reference)
  }
  f
}

.ekb_drop_constant_covariates <- function(data, covariates) {
  covariates[vapply(covariates, function(v) {
    length(unique(stats::na.omit(data[[v]]))) >= 2L
  }, logical(1))]
}

.ekb_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

.ekb_extract_term <- function(tbl, variable) {
  if (!"term" %in% names(tbl)) stop("Result table has no 'term' column.", call. = FALSE)
  idx <- startsWith(tbl$term, variable) | startsWith(tbl$term, .ekb_bt(variable))
  if (sum(idx) != 1L) {
    stop(
      sprintf("Could not uniquely identify the coefficient for '%s'. Matching terms: %s",
              variable, paste(tbl$term[idx], collapse = ", ")),
      call. = FALSE
    )
  }
  tbl[idx, , drop = FALSE]
}

.ekb_format_p <- function(p, digits = 3, threshold = 0.001) {
  out <- rep(NA_character_, length(p))
  ok <- !is.na(p)
  out[ok & p < threshold] <- paste0("<", format(threshold, scientific = FALSE))
  out[ok & p >= threshold] <- formatC(p[ok & p >= threshold], format = "f", digits = digits)
  out
}

.ekb_weighted_ess <- function(weights) {
  if (!length(weights) || all(is.na(weights))) return(NA_real_)
  w <- weights[is.finite(weights) & !is.na(weights)]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

.ekb_as_tibble <- function(x) {
  .ekb_require("tibble")
  tibble::as_tibble(x)
}

.ekb_validate_time <- function(data, time, start = NULL) {
  for (v in c(time, start)) {
    x <- data[[v]]
    if (!is.numeric(x) || any(!is.finite(x) & !is.na(x)) || any(x < 0, na.rm = TRUE))
      stop(sprintf("%s must contain finite, non-negative numeric times or NA.", v), call. = FALSE)
  }
  if (!is.null(start) && any(data[[start]] >= data[[time]], na.rm = TRUE))
    stop("Counting-process intervals must have start < stop.", call. = FALSE)
}

.ekb_complete_ps <- function(data, treatment, covariates) {
  cols <- unique(c(treatment, covariates))
  if (any(!stats::complete.cases(data[, cols, drop = FALSE])))
    stop("Missing treatment/PS covariates: explicitly impute or select complete cases before weighting.", call. = FALSE)
  for (v in cols) {
    if (is.numeric(data[[v]]) && any(!is.finite(data[[v]])))
      stop("Non-finite treatment/PS covariate values are not allowed.", call. = FALSE)
  }
}

.ekb_pool_summary <- function(pooled, conf.level = 0.95) {
  # tidy.mipo offers stable CI names across mice releases.
  .ekb_require("broom")
  broom::tidy(pooled, conf.int = TRUE, conf.level = conf.level, exponentiate = TRUE)
}

.ekb_treatment_term <- function(fit, variable) {
  j <- match(.ekb_bt(variable), attr(stats::terms(fit), "term.labels"))
  if (is.na(j)) j <- match(variable, attr(stats::terms(fit), "term.labels"))
  mm <- stats::model.matrix(fit)
  cols <- colnames(mm)[attr(mm, "assign") == j]
  if (length(cols) != 1L) stop("Treatment must correspond to exactly one model coefficient.", call. = FALSE)
  cols
}
