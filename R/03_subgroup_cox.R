# Subgroup Cox analyses (legacy emR meta-analysis logic, renamed more accurately)

.ekb_subgroup_single <- function(data,
                                 time,
                                 event,
                                 treatment,
                                 subgroup,
                                 level,
                                 covariates = NULL,
                                 weights = NULL,
                                 treatment_reference = NULL,
                                 robust = NULL,
                                 conf.level = 0.95) {
  dat <- as.data.frame(data)
  idx <- !is.na(dat[[subgroup]]) & as.character(dat[[subgroup]]) == as.character(level)
  dat <- dat[idx, , drop = FALSE]
  if (!is.null(weights) && !is.character(weights)) weights <- weights[idx]

  required <- unique(c(time, event, treatment, setdiff(covariates, subgroup)))
  keep <- stats::complete.cases(dat[, required, drop = FALSE])
  if (!is.null(weights)) {
    ww <- .ekb_resolve_weights(dat, weights)
    keep <- keep & !is.na(ww)
    weights <- ww[keep]
  }
  dat <- droplevels(dat[keep, , drop = FALSE])
  if (!nrow(dat) || sum(dat[[event]] == 1, na.rm = TRUE) == 0L) {
    return(data.frame(subgroup = subgroup, level = as.character(level), n = nrow(dat),
      events = 0L, term = NA_character_, hr = NA_real_, conf.low = NA_real_,
      conf.high = NA_real_, p.value = NA_real_, note = "No analyzable observations/events."))
  }
  dat[[treatment]] <- tryCatch(
    .ekb_binary_factor(dat[[treatment]], reference = treatment_reference, name = treatment),
    error = function(e) NULL
  )
  if (is.null(dat[[treatment]])) {
    return(data.frame(
      subgroup = subgroup, level = as.character(level), n = nrow(dat), events = sum(dat[[event]] == 1, na.rm = TRUE),
      term = NA_character_, hr = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_,
      note = "Treatment has fewer than two observed levels in this subgroup.", stringsAsFactors = FALSE
    ))
  }

  covs <- setdiff(covariates %||% character(), c(subgroup, treatment))
  covs <- .ekb_drop_constant_covariates(dat, covs)
  fit <- fit_cox(
    dat,
    time = time,
    event = event,
    covariates = c(treatment, covs),
    weights = weights,
    robust = robust
  )
  tt <- tidy_cox(fit, conf.level = conf.level)
  tr <- tt[tt$term == .ekb_treatment_term(fit$fit, treatment), , drop = FALSE]
  data.frame(
    subgroup = subgroup,
    level = as.character(level),
    n = fit$fit$n,
    events = fit$fit$nevent,
    term = tr$term,
    hr = tr$hr,
    conf.low = tr$hr_conf.low,
    conf.high = tr$hr_conf.high,
    p.value = tr$p.value,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Estimate treatment effects within subgroups
#'
#' Fits the binary treatment contrast separately within each observed subgroup level.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param subgroups Character vector of subgroup-variable names.
#' @param covariates Character vector of prespecified outcome-model covariates. Include treatment explicitly for treatment-effect models.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @param treatment_reference Optional reference level. Otherwise the first existing factor level (or default factor ordering) is used.
#' @param robust NULL delegates variance selection to [survival::coxph()]; non-integer weights normally trigger robust SE. TRUE explicitly requests sandwich variance, including with integer weights; FALSE requests model-based variance.
#' @param include_overall Whether to include a full-cohort treatment-effect row.
#' @param interaction Whether to request optional unweighted interaction likelihood-ratio tests. FALSE by default; weighted interaction LRTs are not returned.
#' @param conf.level Confidence level between zero and one.
#' @return An `ekb_subgroup_cox` list with results, optional interaction_tests and treatment metadata.
#' @details Adjustment excludes the grouping variable and locally constant covariates. Global input weights are subset without refitting. P-values test the treatment effect within each level, not interaction. Sparse levels without observations/events or both treatments return missing effects with a note. Additional treatment adjustment and reference coding follow fit_cox.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' sg <- subgroup_cox(d, "time", "event", "treatment", "sex")
#' sg$results
subgroup_cox <- function(data,
                         time,
                         event,
                         treatment,
                         subgroups,
                         covariates = NULL,
                         weights = NULL,
                         treatment_reference = NULL,
                         robust = NULL,
                         include_overall = TRUE,
                         interaction = FALSE,
                         conf.level = 0.95) {
  .ekb_assert_columns(data, c(time, event, treatment, subgroups, covariates))

  if (inherits(weights, "ekb_weighting")) {
    message(
      "Using propensity-score weights from fit_weighting() for subgroup Cox analysis ",
      "(estimand = ", weights$estimand, ")."
    )
  }

  dat <- as.data.frame(data)
  dat[[event]] <- .ekb_binary_event(dat[[event]], event)
  dat[[treatment]] <- .ekb_binary_factor(dat[[treatment]], reference = treatment_reference, name = treatment)
  w <- .ekb_resolve_weights(dat, weights)

  out <- list()
  k <- 1L
  for (sg in subgroups) {
    f <- droplevels(factor(dat[[sg]]))
    if (nlevels(f) < 2L) {
      warning(sprintf("Skipping subgroup '%s': fewer than two observed levels.", sg), call. = FALSE)
      next
    }
    for (lev in levels(f)) {
      out[[k]] <- .ekb_subgroup_single(
        data = dat,
        time = time,
        event = event,
        treatment = treatment,
        subgroup = sg,
        level = lev,
        covariates = covariates,
        weights = w,
        treatment_reference = treatment_reference,
        robust = robust,
        conf.level = conf.level
      )
      k <- k + 1L
    }
  }
  res <- if (length(out)) do.call(rbind, out) else data.frame()

  overall <- NULL
  if (isTRUE(include_overall)) {
    covs <- setdiff(covariates %||% character(), treatment)
    covs <- .ekb_drop_constant_covariates(dat, covs)
    fit <- fit_cox(dat, time, event, c(treatment, covs), weights = w, robust = robust)
    tt <- tidy_cox(fit, conf.level = conf.level)
    tr <- tt[tt$term == .ekb_treatment_term(fit$fit, treatment), , drop = FALSE]
    overall <- data.frame(
      subgroup = "Overall",
      level = "Overall",
      n = fit$fit$n,
      events = fit$fit$nevent,
      term = tr$term,
      hr = tr$hr,
      conf.low = tr$hr_conf.low,
      conf.high = tr$hr_conf.high,
      p.value = tr$p.value,
      note = NA_character_,
      stringsAsFactors = FALSE
    )
    res <- rbind(res, overall)
  }

  interaction_tests <- NULL
  if (isTRUE(interaction)) {
    if (!is.null(w)) {
      warning("Interaction LRTs with IPTW are not returned by default because the likelihood-ratio test does not use the robust sandwich variance. Subgroup HR/p-values are still returned.", call. = FALSE)
    } else {
      valid_sg <- subgroups[vapply(subgroups, function(sg) length(unique(stats::na.omit(dat[[sg]]))) >= 2L, logical(1))]
      its <- lapply(valid_sg, function(sg) {
        covs <- setdiff(covariates %||% character(), c(treatment, sg))
        main_terms <- c(treatment, sg, covs)
        int_term <- paste0(.ekb_bt(treatment), "*", .ekb_bt(sg))
        response <- sprintf("survival::Surv(%s, %s)", .ekb_bt(time), .ekb_bt(event))
        f0 <- stats::as.formula(paste(response, "~", paste(.ekb_bt(main_terms), collapse = " + ")))
        f1 <- stats::as.formula(paste(response, "~", paste(c(int_term, .ekb_bt(covs)), collapse = " + ")))
        m0 <- survival::coxph(f0, data = dat)
        m1 <- survival::coxph(f1, data = dat)
        a <- stats::anova(m0, m1, test = "LRT")
        data.frame(subgroup = sg, interaction_p = a[["P(>|Chi|)"]][2], stringsAsFactors = FALSE)
      })
      interaction_tests <- do.call(rbind, its)
    }
  }

  structure(
    list(
      results = .ekb_as_tibble(res),
      interaction_tests = if (is.null(interaction_tests)) NULL else .ekb_as_tibble(interaction_tests),
      treatment = treatment,
      treatment_reference = levels(dat[[treatment]])[1],
      treatment_comparator = levels(dat[[treatment]])[2],
      weighted = !is.null(w),
      weighting = if (inherits(weights, "ekb_weighting")) weights else NULL,
      call = match.call()
    ),
    class = "ekb_subgroup_cox"
  )
}
