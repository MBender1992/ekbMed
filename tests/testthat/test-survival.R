test_that("unweighted KM and log-rank match survival", {
  d <- simulate_clinical()
  a <- fit_survival(d, "time", "event", "treatment")
  b <- survival::survfit(survival::Surv(time, event) ~ treatment, data = d, conf.type = "log-log")
  expect_s3_class(a, "ekb_survival_fit")
  for (nm in c("time", "surv", "lower", "upper", "n.risk")) expect_equal(a$fit[[nm]], b[[nm]])
  lr <- survival_logrank(d, "time", "event", "treatment")
  direct <- survival::survdiff(survival::Surv(time, event) ~ treatment, data = d)
  expect_equal(lr$statistic, unname(direct$chisq))
  expect_equal(lr$p.value, pchisq(direct$chisq, 1, lower.tail = FALSE))
})

test_that("weighted KM retains raw counts and matches survfit", {
  d <- simulate_clinical()
  w <- seq(0.8, 2.2, length.out = nrow(d))
  a <- fit_survival(d, "time", "event", "treatment", weights = w)
  b <- survival::survfit(survival::Surv(time, event) ~ treatment, data = d, weights = w, conf.type = "log-log")
  expect_equal(a$fit$surv, b$surv)
  expect_equal(a$fit$std.err, b$std.err)
  sm <- summarize_survival(a, c(12, 24))
  expect_equal(sum(sm$median$n), nrow(d))
  expect_equal(sum(sm$median$events), sum(d$event))
  expect_true(sm$weighted)
  direct <- summary(b, times = c(12, 24))
  expect_equal(sm$landmarks$estimate, direct$surv)
  expect_equal(sm$median$median, unname(summary(b)$table[, "median"]))
})

test_that("weight objects can be reused for survival", {
  need_weighting()
  d <- simulate_clinical()
  w <- fit_weighting(d, "treatment", "age")
  expect_equal(fit_survival(d, "time", "event", weights = w)$fit$surv,
    fit_survival(d, "time", "event", weights = w$weights)$fit$surv)
})

test_that("weighted log-rank matches RISCA", {
  skip_if_not_installed("RISCA")
  d <- simulate_clinical()
  w <- seq(1, 2, length.out = nrow(d))
  a <- survival_logrank(d, "time", "event", "treatment", weights = w)
  b <- RISCA::ipw.log.rank(times = d$time, failures = d$event,
    variable = as.integer(d$treatment == "B"), weights = w)
  expect_equal(a$p.value, unname(b$p.value))
  expect_equal(a$statistic, unname(b$statistic))
})

test_that("landmarks and missingness are handled consistently", {
  d <- simulate_clinical()
  d$time[1] <- NA
  a <- fit_survival(d, "time", "event", "treatment")
  sm <- summarize_survival(a, landmarks = c(12, 1e6))
  expect_equal(sum(sm$median$n), nrow(d) - 1)
  expect_false(any(sm$landmarks$time == 1e6))
  expect_error(summarize_survival(a, -1), "non-negative")
  expect_equal(nrow(summarize_survival(a, numeric())$landmarks), 0)
  sm90 <- summarize_survival(a, 12, conf.level = 0.9)
  direct <- fit_survival(d, "time", "event", "treatment", conf.level = 0.9)
  expect_equal(sm90$landmarks$conf.low, summary(direct$fit, times = 12)$lower)
})
