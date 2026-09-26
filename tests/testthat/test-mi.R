test_that("imputation is explicit deterministic and reusable", {
  skip_if_not_installed("mice")
  mi <- make_mi()
  other <- make_mi()
  expect_s3_class(mi, "ekb_imputation")
  expect_equal(mi$mids$imp, other$mids$imp)
  expect_false(anyNA(mice::complete(mi$mids, 1)$age))
  expect_true(".ekb_nelson_aalen" %in% names(mi$imputation_data))
  expect_equal(mi$predictorMatrix["age", "event"], 1)
  a <- fit_mi_cox(mi, "time", "event", c("treatment", "age"), treatment = "treatment")
  expect_identical(a$imputation, mi)
  expect_false(a$weighted)
  expect_null(a$weighting)
  direct <- lapply(seq_len(mi$m), function(i) {
    d <- mice::complete(mi$mids, i)
    survival::coxph(survival::Surv(time, event) ~ treatment + age, data = d)
  })
  b <- mice::pool(mice::as.mira(direct))
  expect_equal(a$pooled$pooled$estimate, b$pooled$estimate)
  expect_equal(a$pooled$pooled$t, b$pooled$t)
})

test_that("MI refuses implicit imputation of outcomes and predictors", {
  skip_if_not_installed("mice")
  d <- simulate_clinical()[, c("time", "event", "treatment", "age", "sex")]
  expect_error(impute_clinical(d, "event", time = "time", event = "event"), "should not")
  d$sex[1] <- NA
  expect_error(impute_clinical(d, "age", predictor_vars = "sex"), "Missing predictors")
  expect_error(fit_mi_cox(d, "time", "event", "age"), "imputation must")
})

test_that("MI IPTW estimates weights per imputation and pools robust variance", {
  skip_if_not_installed("mice")
  need_weighting()
  mi <- make_mi()
  a <- capture_messages(fit_mi_cox(mi, "time", "event", c("treatment", "age"),
    treatment = "treatment", ps_covariates = c("age", "sex"), robust = TRUE))
  fit <- a$value
  expect_length(a$messages, 1)
  expect_match(a$messages, "separately within each")
  expect_identical(fit$imputation, mi)
  expect_length(fit$weighting, mi$m)
  expect_equal(nrow(fit$diagnostic_summary), mi$m)
  for (i in seq_len(mi$m)) {
    d <- mice::complete(mi$mids, i)
    w <- WeightIt::weightit(treatment ~ age + sex, data = d, method = "glm", estimand = "ATE", stabilize = FALSE)$weights
    expect_equal(fit$weighting[[i]]$weights, as.numeric(w))
    direct <- survival::coxph(survival::Surv(time, event) ~ treatment + age,
      data = d, weights = w, robust = TRUE)
    expect_equal(coef(fit$fits[[i]]), coef(direct))
    expect_equal(vcov(fit$fits[[i]]), vcov(direct))
  }
  robust_vars <- vapply(fit$fits, function(f) unname(diag(vcov(f))[1]), numeric(1))
  estimates <- vapply(fit$fits, function(f) unname(coef(f)[1]), numeric(1))
  pooled <- fit$pooled$pooled[fit$pooled$pooled$term == "treatmentB", ]
  expect_equal(pooled$ubar, mean(robust_vars))
  expect_equal(pooled$t, mean(robust_vars) + (1 + 1 / mi$m) * var(estimates))
})

test_that("MI subgroup scopes match direct pooled models", {
  skip_if_not_installed("mice")
  need_weighting()
  mi <- make_mi(impute_sex = TRUE)
  for (scope in c("unweighted", "global", "within_subgroup")) {
    weighted <- scope != "unweighted"
    cap <- capture_messages(subgroup_mi_cox(mi, "time", "event", "treatment", "sex",
      covariates = "age", ps_covariates = if (weighted) c("age", "sex") else NULL,
      weight_scope = if (weighted) scope else "global", robust = TRUE))
    a <- cap$value
    expect_identical(a$imputation, mi)
    expect_null(a$interaction_tests)
    expect_length(cap$messages, as.integer(weighted))
    if (weighted) expect_match(cap$messages, scope, fixed = TRUE)
    for (lev in c("F", "M")) {
      fits <- lapply(seq_len(mi$m), function(i) {
        d <- mice::complete(mi$mids, i)
        idx <- d$sex == lev
        w <- rep(1, sum(idx))
        if (scope == "global") {
          w <- WeightIt::weightit(treatment ~ age + sex, data = d, method = "glm", estimand = "ATE")$weights[idx]
        }
        d <- droplevels(d[idx, ])
        if (scope == "within_subgroup") {
          w <- WeightIt::weightit(treatment ~ age, data = d, method = "glm", estimand = "ATE")$weights
        }
        survival::coxph(survival::Surv(time, event) ~ treatment + age, data = d, weights = w, robust = TRUE)
      })
      direct <- broom::tidy(mice::pool(mice::as.mira(fits)), conf.int = TRUE, exponentiate = TRUE)
      z <- direct[direct$term == "treatmentB", ]
      row <- a$results[a$results$level == lev, ]
      expect_equal(row$hr, z$estimate)
      expect_equal(row$p.value, z$p.value)
      expect_equal(row$n_mean, mean(vapply(fits, function(f) f$n, numeric(1))))
    }
    expect_equal(a$results$n[a$results$level == "Overall"], nrow(mi$original_data))
  }
})
