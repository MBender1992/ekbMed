test_that("ORR DCR and unknown denominators preserve conventions", {
  x <- c("CR", "PR", "SD", "PD", NA, "Unknown")
  expect_equal(calc_ORR(x)$estimate, 0.5)
  expect_equal(calc_DCR(x)$estimate, 0.75)
  expect_equal(calc_ORR(x)$n, 4)
  expect_equal(calc_ORR(x)$conf.low, unname(binom.test(2, 4)$conf.int[1]))
  expect_true(is.na(calc_ORR(x, na.rm = FALSE)$estimate))
  expect_equal(calc_ORR(c(NA, "Unknown"))$n, 0)
})

test_that("response conversions retain the strict mixed-response cutoff", {
  a <- convert_response(c("NC", "NED", "MR", "MR", "MR", "NB"),
    trtdur = c(NA, NA, 90, 91, NA, NA))
  expect_identical(a, c("SD", "CR", "PD", "SD", NA_character_, NA_character_))
  expect_error(convert_response("MR"), "trtdur")
})

test_that("date conversion and inclusive intervals are explicit", {
  a <- convert_date(c("2020", "2020-05", "21.05.2020", "", "bad"), return_details = TRUE)
  expect_equal(a$date[1:3], as.Date(c("2020-06-15", "2020-05-15", "2020-05-21")))
  expect_identical(a$imputed[1:3], c(TRUE, TRUE, FALSE))
  expect_true(a$invalid[5])
  expect_error(convert_date("bad", invalid = "error"), "Could not parse")
  expect_equal(calc_survival_interval("2020-01-01", "2020-01-02", unit = "days"), 2)
  expect_equal(calc_survival_interval("2020-01-01", "2020-01-02", unit = "days", inclusive = FALSE), 1)
  expect_error(calc_survival_interval("2020-01-02", "2020-01-01"), "before")
})

test_that("clinical QC does not mutate input and detects temporal errors", {
  d <- data.frame(id = c(1, 1, 3), age = c(50, -2, NA),
    start = as.Date(c("2020-01-01", "2020-02-01", "2020-01-01")),
    end = as.Date(c("2020-01-02", "2020-01-01", NA)))
  before <- d
  q <- qc_clinical_data(d, id = "id", date_vars = c("start", "end"),
    nonnegative_vars = "age", required_complete_vars = "age",
    temporal_rules = data.frame(check = "start_before_end", earlier = "start", later = "end"), verbose = FALSE)
  expect_identical(d, before)
  expect_s3_class(q, "ekb_clinical_qc")
  expect_gt(nrow(q$issues), 0)
  expect_true("start_before_end" %in% q$issues$check)
  expect_output(print(q), "QC|quality|Clinical|clinical")
})
