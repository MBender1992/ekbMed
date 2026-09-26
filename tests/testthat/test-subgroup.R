test_that("every subgroup matches direct treatment effects", {
  d <- simulate_clinical()
  for (weighted in c(FALSE, TRUE)) {
    w <- if (weighted) seq(0.8, 2, length.out = nrow(d)) else NULL
    a <- subgroup_cox(d, "time", "event", "treatment", "sex",
      covariates = c("age", "constant", "sex"), weights = w, robust = TRUE)
    expect_s3_class(a, "ekb_subgroup_cox")
    expect_null(a$interaction_tests)
    expect_equal(nrow(a$results), 3)
    for (lev in levels(d$sex)) {
      idx <- d$sex == lev
      ds <- d[idx, ]
      ww <- if (weighted) w[idx] else rep(1, nrow(ds))
      b <- survival::coxph(survival::Surv(time, event) ~ treatment + age,
        data = ds, weights = ww, robust = TRUE)
      z <- a$results[a$results$level == lev, ]
      expect_equal(z$hr, unname(exp(coef(b)[1])))
      expect_equal(z$p.value, unname(summary(b)$coefficients[1, "Pr(>|z|)"]))
      expect_equal(z$n, b$n)
    }
    expect_equal(a$results$n[a$results$level == "Overall"], nrow(d))
  }
})

test_that("subgroup weight object handling is quiet internally", {
  need_weighting()
  d <- simulate_clinical()
  w <- fit_weighting(d, "treatment", "age")
  a <- capture_messages(subgroup_cox(d, "time", "event", "treatment", "sex", weights = w))
  expect_length(a$messages, 1)
  b <- subgroup_cox(d, "time", "event", "treatment", "sex", weights = w$weights)
  expect_equal(a$value$results, b$results)
})

test_that("reference changes invert subgroup treatment effects", {
  d <- simulate_clinical()
  a <- subgroup_cox(d, "time", "event", "treatment", "sex", treatment_reference = "A")
  b <- subgroup_cox(d, "time", "event", "treatment", "sex", treatment_reference = "B")
  expect_equal(a$results$hr, 1 / b$results$hr)
})

test_that("sparse subgroups return a diagnostic note", {
  d <- simulate_clinical()
  d$sg <- ifelse(d$treatment == "A", "A_only", "B_only")
  a <- subgroup_cox(d, "time", "event", "treatment", "sg", include_overall = FALSE)
  expect_true(all(is.na(a$results$hr)))
  expect_true(all(!is.na(a$results$note)))
})
