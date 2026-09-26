# Multiple-imputation subgroup Cox analyses
# Historical filename retained for source-order compatibility.

.ekb_subgroup_mi_diag_row <- function(diagnostic,
                                      imputation,
                                      subgroup,
                                      level,
                                      scope) {
  data.frame(
    imputation = imputation,
    subgroup = subgroup,
    level = as.character(level),
    scope = scope,
    max_abs_smd = diagnostic$max_abs_smd,
    min_ess = min(diagnostic$ess$ess, na.rm = TRUE),
    max_weight = max(diagnostic$weight_quantiles$weight, na.rm = TRUE),
    positivity_flag = diagnostic$positivity_flag,
    stringsAsFactors = FALSE
  )
}

#' Pool subgroup treatment effects across imputations
#'
#' Reuses an upstream imputation with optional global or within-subgroup propensity-score weighting.
#'
#' @param imputation An existing `ekb_imputation` or `mice::mids` object. It is reused and never re-imputed downstream.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param treatment Name of the binary treatment column. The second factor level is compared with the first.
#' @param subgroups Character vector of subgroup-variable names.
#' @param covariates Character vector of prespecified outcome-model covariates. Include treatment explicitly for treatment-effect models.
#' @param ps_covariates Character vector of baseline propensity-score predictors. In MI functions, NULL disables weighting; supplying predictors enables within-imputation weighting.
#' @param estimand Target estimand passed to [WeightIt::weightit()], default `"ATE"`.
#' @param weight_method Propensity-score estimation method passed to WeightIt.
#' @param stabilize Whether to request stabilized weights from WeightIt.
#' @param weight_scope `"global"` estimates weights on the full cohort within each imputation, then subsets them; `"within_subgroup"` re-estimates weights within each subgroup level and imputation. The overall row always uses full-cohort weights.
#' @param treatment_reference Optional reference level. Otherwise the first existing factor level (or default factor ordering) is used.
#' @param robust NULL delegates variance selection to [survival::coxph()]; non-integer weights normally trigger robust SE. TRUE explicitly requests sandwich variance, including with integer weights; FALSE requests model-based variance.
#' @param include_overall Whether to include a full-cohort treatment-effect row.
#' @param conf.level Confidence level between zero and one.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return An `ekb_subgroup_mi_cox` list with results, imputation, weighting and diagnostic_summary. Counts include n_mean and events_mean across imputations.
#' @details Each level is analysed within every completed imputation and treatment effects are pooled with mice::pool. Global weights are estimated once per imputation; within_subgroup removes the grouping variable and constant PS predictors before fitting within each level. At least one varying PS predictor must remain. Empty levels or loss of either treatment level raise an error. No interaction tests are produced. Additional arguments are passed to WeightIt::weightit. Returned N is the rounded mean across imputations; n_mean retains the unrounded mean. Robust pooling follows fit_mi_cox.
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
#'   subgroup_mi_cox(mi, "time", "event", "treatment", "sex")
#' }
subgroup_mi_cox <- function(imputation,
                            time,
                            event,
                            treatment,
                            subgroups,
                            covariates = NULL,
                            ps_covariates = NULL,
                            estimand = "ATE",
                            weight_method = "glm",
                            stabilize = FALSE,
                            weight_scope = c("global", "within_subgroup"),
                            treatment_reference = NULL,
                            robust = NULL,
                            include_overall = TRUE,
                            conf.level = 0.95,
                            ...) {
  .ekb_require("mice")
  mi <- .ekb_get_mids(imputation)
  weight_scope <- match.arg(weight_scope)

  d1 <- mice::complete(mi, 1)
  .ekb_assert_columns(d1, c(time, event, treatment, subgroups, covariates, ps_covariates))
  weighted <- length(ps_covariates %||% character()) > 0L

  if (weighted) {
    message(
      "IPTW enabled: propensity scores and weights will be estimated ",
      "within each imputed dataset using weight_scope = '",
      weight_scope,
      "' (estimand = ", estimand, ")."
    )
  }

  completed <- lapply(seq_len(mi$m), function(i) {
    d <- mice::complete(mi, i)
    d[[treatment]] <- .ekb_binary_factor(d[[treatment]], reference = treatment_reference, name = treatment)
    d
  })

  treatment_reference_final <- levels(completed[[1]][[treatment]])[1]
  treatment_comparator_final <- levels(completed[[1]][[treatment]])[2]

  # One full-cohort PS model per imputation. These weights are used for
  # weight_scope = "global" and for the overall row in both weighting scopes.
  global_weighting <- NULL
  global_diagnostics <- NULL
  diagnostic_rows <- list()
  d_index <- 1L

  if (weighted) {
    global_weighting <- lapply(seq_len(mi$m), function(i) {
      d <- completed[[i]]
      psc <- .ekb_drop_constant_covariates(d, ps_covariates)
      suppressMessages(
        fit_weighting(
          d,
          treatment = treatment,
          ps_covariates = psc,
          method = weight_method,
          estimand = estimand,
          stabilize = stabilize,
          treatment_reference = treatment_reference_final,
          ...
        )
      )
    })
    global_diagnostics <- lapply(global_weighting, diagnose_weighting)
    for (i in seq_len(mi$m)) {
      diagnostic_rows[[d_index]] <- .ekb_subgroup_mi_diag_row(
        global_diagnostics[[i]], i, "Overall", "Overall", "global"
      )
      d_index <- d_index + 1L
    }
  }

  subgroup_levels <- function(variable) {
    x <- completed[[1]][[variable]]
    if (is.factor(x)) return(levels(x))
    unique(unlist(lapply(completed, function(d) as.character(stats::na.omit(d[[variable]])))))
  }

  results <- list()
  r_index <- 1L

  for (sg in subgroups) {
    for (lev in subgroup_levels(sg)) {
      fits <- vector("list", mi$m)
      ns <- events <- numeric(mi$m)

      for (i in seq_len(mi$m)) {
        d <- completed[[i]]
        idx <- !is.na(d[[sg]]) & as.character(d[[sg]]) == as.character(lev)
        dsub <- d[idx, , drop = FALSE]

        if (!nrow(dsub)) {
          stop(sprintf("Subgroup '%s', level '%s' is empty in imputation %s.", sg, lev, i), call. = FALSE)
        }
        dsub[[treatment]] <- .ekb_binary_factor(
          dsub[[treatment]], reference = treatment_reference_final, name = treatment
        )

        covs <- setdiff(covariates %||% character(), c(sg, treatment))
        covs <- .ekb_drop_constant_covariates(dsub, covs)
        w <- NULL

        if (weighted) {
          if (weight_scope == "global") {
            w <- global_weighting[[i]]$weights[idx]
          } else {
            psc <- setdiff(ps_covariates, sg)
            psc <- .ekb_drop_constant_covariates(dsub, psc)
            wfit <- suppressMessages(
              fit_weighting(
                dsub,
                treatment = treatment,
                ps_covariates = psc,
                method = weight_method,
                estimand = estimand,
                stabilize = stabilize,
                treatment_reference = treatment_reference_final,
                ...
              )
            )
            w <- wfit$weights
            diagnostic_rows[[d_index]] <- .ekb_subgroup_mi_diag_row(
              diagnose_weighting(wfit), i, sg, lev, "within_subgroup"
            )
            d_index <- d_index + 1L
          }
        }

        fit <- fit_cox(
          dsub,
          time = time,
          event = event,
          covariates = c(treatment, covs),
          weights = w,
          robust = robust
        )
        fits[[i]] <- fit$fit
        ns[i] <- fit$fit$n
        events[i] <- fit$fit$nevent
      }

      pooled <- mice::pool(mice::as.mira(fits))
      sm <- .ekb_pool_summary(pooled, conf.level)
      tr <- sm[sm$term == .ekb_treatment_term(fits[[1]], treatment), , drop = FALSE]
      if (nrow(tr) != 1L) {
        stop(sprintf("Could not identify one treatment coefficient for '%s' = '%s'.", sg, lev), call. = FALSE)
      }

      results[[r_index]] <- data.frame(
        subgroup = sg,
        level = as.character(lev),
        n = as.integer(round(mean(ns))),
        events = as.integer(round(mean(events))),
        n_mean = mean(ns),
        events_mean = mean(events),
        term = tr$term,
        hr = tr$estimate,
        conf.low = tr$conf.low,
        conf.high = tr$conf.high,
        p.value = tr$p.value,
        stringsAsFactors = FALSE
      )
      r_index <- r_index + 1L
    }
  }

  if (isTRUE(include_overall)) {
    fits <- vector("list", mi$m)
    ns <- events <- numeric(mi$m)

    for (i in seq_len(mi$m)) {
      d <- completed[[i]]
      covs <- setdiff(covariates %||% character(), treatment)
      covs <- .ekb_drop_constant_covariates(d, covs)
      w <- if (weighted) global_weighting[[i]]$weights else NULL

      fit <- fit_cox(
        d,
        time = time,
        event = event,
        covariates = c(treatment, covs),
        weights = w,
        robust = robust
      )
      fits[[i]] <- fit$fit
      ns[i] <- fit$fit$n
      events[i] <- fit$fit$nevent
    }

    pooled <- mice::pool(mice::as.mira(fits))
    sm <- .ekb_pool_summary(pooled, conf.level)
    tr <- sm[sm$term == .ekb_treatment_term(fits[[1]], treatment), , drop = FALSE]

    results[[r_index]] <- data.frame(
      subgroup = "Overall",
      level = "Overall",
      n = as.integer(round(mean(ns))),
      events = as.integer(round(mean(events))),
      n_mean = mean(ns),
      events_mean = mean(events),
      term = tr$term,
      hr = tr$estimate,
      conf.low = tr$conf.low,
      conf.high = tr$conf.high,
      p.value = tr$p.value,
      stringsAsFactors = FALSE
    )
  }

  diagnostic_summary <- if (length(diagnostic_rows)) {
    .ekb_as_tibble(do.call(rbind, diagnostic_rows))
  } else {
    NULL
  }

  structure(
    list(
      results = .ekb_as_tibble(do.call(rbind, results)),
      interaction_tests = NULL,
      imputation = imputation,
      weighted = weighted,
      weighting = global_weighting,
      diagnostics = global_diagnostics,
      diagnostic_summary = diagnostic_summary,
      weight_scope = if (weighted) weight_scope else NA_character_,
      estimand = if (weighted) estimand else NA_character_,
      weight_method = if (weighted) weight_method else NA_character_,
      stabilize = if (weighted) stabilize else NA,
      treatment = treatment,
      treatment_reference = treatment_reference_final,
      treatment_comparator = treatment_comparator_final,
      subgroups = subgroups,
      covariates = covariates,
      ps_covariates = ps_covariates,
      call = match.call()
    ),
    class = "ekb_subgroup_mi_cox"
  )
}
