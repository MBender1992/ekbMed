test_that("legacy wrappers delegate to current survival and Cox APIs", {
  d <- simulate_clinical()
  direct <- fit_cox(d, "time", "event", c("treatment", "age"))
  legacy <- cox_output(d, "time", "event", c("treatment", "age"), output = "fit")
  expect_equal(coef(legacy), coef(direct$fit))
  expect_warning(cox_output(d, "time", "event", "age", modeltype = "backward"), "full model")
  a <- survival_time(d, "time", "event", "treatment", times = 12)
  b <- summarize_survival(fit_survival(d, "time", "event", "treatment"), 12)$landmarks
  expect_equal(a$estimate, b$estimate)
  expect_equal(add_median_survival(d, "time", "event", "treatment", statistics = FALSE)$median,
    summarize_survival(fit_survival(d, "time", "event", "treatment"), numeric())$median$median)
  expect_equal(calc_survival("2020-01-01", "2020-03-01"),
    calc_survival_interval("2020-01-01", "2020-03-01"))
  sg <- coxph_meta_analysis(d, "time", "event", "age", "sex", "treatment")
  expect_equal(sg, subgroup_cox(d, "time", "event", "treatment", "sex", covariates = "age")$results)
})

test_that("legacy ATE and MI wrappers reuse public implementations", {
  need_weighting()
  skip_if_not_installed("mice")
  d <- simulate_clinical()[, c("time", "event", "treatment", "age", "sex")]
  expect_equal(ate_weights(d, c("treatment", "age"), "treatment"),
    fit_weighting(d, "treatment", "age")$weights)
  d$age[seq(5, nrow(d), 11)] <- NA_real_
  expect_warning(a <- suppressMessages(fit_mi_iptw_cox(d, "time", "event", "treatment",
    "age", c("age", "sex"), outcome_covariates = "age", m = 2, maxit = 2, seed = 71)),
    "retained for compatibility")
  expect_s3_class(a, "ekb_mi_cox")
  b <- suppressMessages(fit_mi_cox(a$imputation, "time", "event", c("treatment", "age"),
    treatment = "treatment", ps_covariates = c("age", "sex")))
  expect_equal(a$pooled_summary, b$pooled_summary)
})

test_that("optional adjustedCurves pathway fits and plots", {
  skip_if_not_installed("adjustedCurves")
  d <- simulate_clinical()
  w <- manual_ate_weights(d, "treatment", "age")$weights
  a <- fit_iptw_survival(d, "time", "event", "treatment", w, conf_int = FALSE)
  expect_s3_class(a, "ekb_iptw_survival")
  expect_equal(a$weights, w)
  expect_s3_class(plot_iptw_survival(a, conf_int = FALSE, risk_table = FALSE), "ggplot")
})
