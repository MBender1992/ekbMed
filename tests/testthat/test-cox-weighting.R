test_that("unweighted and weighted Cox match direct fits", {
  d <- simulate_clinical()
  a <- fit_cox(d, "time", "event", c("treatment", "age"))
  b <- survival::coxph(survival::Surv(time, event) ~ treatment + age, data = d)
  expect_equal(coef(a$fit), coef(b))
  expect_equal(vcov(a$fit), vcov(b))
  w <- seq(0.7, 2, length.out = nrow(d))
  aw <- fit_cox(d, "time", "event", c("treatment", "age"), weights = w, robust = TRUE)
  bw <- survival::coxph(survival::Surv(time, event) ~ treatment + age, data = d, weights = w, robust = TRUE)
  expect_equal(coef(aw$fit), coef(bw))
  expect_equal(vcov(aw$fit), vcov(bw))
  expect_equal(tidy_cox(aw)$std.error, unname(sqrt(diag(vcov(bw)))))
  expect_equal(tidy_cox(aw)$hr, unname(exp(coef(bw))))
})

test_that("weight objects produce the same Cox fit and one message", {
  need_weighting()
  d <- simulate_clinical()
  w <- fit_weighting(d, "treatment", c("age", "sex"))
  out <- capture_messages(fit_cox(d, "time", "event", "treatment", weights = w))
  expect_length(out$messages, 1)
  expect_match(out$messages, "estimand = ATE", fixed = TRUE)
  b <- fit_cox(d, "time", "event", "treatment", weights = w$weights)
  expect_equal(coef(out$value$fit), coef(b$fit))
  expect_identical(out$value$weighting, w)
})

test_that("WeightIt ATE agrees with direct and manual implementations", {
  need_weighting()
  d <- simulate_clinical()
  a <- fit_weighting(d, "treatment", c("age", "sex"), treatment_reference = "B")
  d$treatment <- relevel(d$treatment, "B")
  b <- WeightIt::weightit(treatment ~ age + sex, data = d, method = "glm", estimand = "ATE", stabilize = FALSE)
  manual <- manual_ate_weights(d, "treatment", c("age", "sex"), "B")
  expect_s3_class(a, "ekb_weighting")
  expect_identical(levels(a$data$treatment), c("B", "A"))
  expect_equal(a$weights, as.numeric(b$weights), tolerance = 1e-7)
  expect_equal(a$weights, unname(manual$weights), tolerance = 1e-6)
  expect_equal(a$propensity_score, as.numeric(b$ps), tolerance = 1e-7)
  expect_true(all(a$weights > 0))
})

test_that("weighting diagnostics expose overlap balance and ESS", {
  need_weighting()
  d <- simulate_clinical()
  w <- fit_weighting(d, "treatment", c("age", "sex"))
  a <- diagnose_weighting(w)
  expect_true(all(c("balance", "ess", "weight_quantiles", "positivity_flag", "max_abs_smd") %in% names(a)))
  expect_equal(a$ess$ess[1], sum(w$weights[d$treatment == "A"])^2 / sum(w$weights[d$treatment == "A"]^2))
  w$propensity_score[1] <- 0.001
  expect_true(diagnose_weighting(w)$positivity_flag)
})

test_that("weighting rejects missing predictors and invalid treatment", {
  need_weighting()
  d <- simulate_clinical()
  d$age[1] <- NA
  expect_error(fit_weighting(d, "treatment", "age"), "Missing treatment/PS")
  d$age[1] <- 60
  d$treatment <- "A"
  expect_error(fit_weighting(d, "treatment", "age"), "exactly two")
  expect_error(fit_weighting(d, "treatment", character()), "at least one")
})

test_that("truncation is explicit and matches quantile capping", {
  w <- c(1, 2, 3, 4, 20)
  a <- truncate_weights(w, c(0.1, 0.9))
  cap <- quantile(w, c(0.1, 0.9), names = FALSE)
  expect_equal(a$weights, pmin(pmax(w, cap[1]), cap[2]))
  expect_identical(a$original, w)
  expect_error(truncate_weights(w, c(0.9, 0.1)), "increasing")
})

test_that("PH and Wald diagnostics return standard objects", {
  d <- simulate_clinical()
  a <- fit_cox(d, "time", "event", c("treatment", "age"))
  expect_s3_class(check_cox_ph(a), "cox.zph")
  expect_equal(cox_global_wald(a)$p.value, unname(summary(a$fit)$waldtest[3]))
  expect_true("HR (95% CI)" %in% names(format_cox(a)))
})
