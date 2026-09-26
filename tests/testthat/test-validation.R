test_that("missing columns and invalid events fail early", {
  d <- simulate_clinical()
  expect_error(fit_cox(d, "absent", "event", "age"), "Missing column")
  d$event[1] <- 2
  expect_error(fit_survival(d, "time", "event"), "coded 0/1")
  expect_error(fit_cox(d, "time", "event", "age"), "coded 0/1")
})

test_that("weight shape and values are validated", {
  d <- simulate_clinical()
  for (w in list(1:2, rep(-1, nrow(d)), rep(0, nrow(d)), rep(Inf, nrow(d)), "absent")) {
    expect_error(fit_cox(d, "time", "event", "age", weights = w))
  }
  d$w <- rep(1.2, nrow(d))
  a <- fit_cox(d, "time", "event", "age", weights = "w")
  b <- fit_cox(d, "time", "event", "age", weights = d$w)
  expect_equal(coef(a$fit), coef(b$fit))
})

test_that("binary treatment and follow-up are validated", {
  d <- simulate_clinical()
  d$treatment <- rep(c("A", "B", "C"), length.out = nrow(d))
  expect_error(subgroup_cox(d, "time", "event", "treatment", "sex"), "exactly two")
  d$time[1] <- -1
  expect_error(fit_survival(d, "time", "event"), "non-negative")
  d$time[1] <- 1
  d$start <- d$time
  expect_error(fit_cox(d, "time", "event", "age", start = "start"), "start < stop")
})

test_that("logical and factor events retain 0/1 coding", {
  expect_identical(ekbMed:::.ekb_binary_event(c(TRUE, FALSE, NA)), c(1L, 0L, NA_integer_))
  expect_identical(ekbMed:::.ekb_binary_event(factor(c("1", "0"))), c(1L, 0L))
})

test_that("ESS matches the direct formula", {
  w <- c(1, 2, 3, 4)
  expect_equal(ekbMed:::.ekb_weighted_ess(w), sum(w)^2 / sum(w^2))
  expect_true(is.na(ekbMed:::.ekb_weighted_ess(numeric())))
})
