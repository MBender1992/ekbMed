test_that("baseline tables provide Total and weighted headers", {
  skip_if_not_installed("gtsummary")
  skip_if_not_installed("survey")
  need_weighting()
  d <- simulate_clinical()
  a <- tbl_baseline(d, "treatment", c("age", "sex"), add_p = FALSE)
  expect_s3_class(a, "gtsummary")
  expect_true("stat_0" %in% names(a$table_body))
  w <- fit_weighting(d, "treatment", "age")
  cap <- capture_messages(tbl_baseline(d, "treatment", c("age", "sex"), weights = w, add_smd = FALSE))
  b <- cap$value
  expect_length(cap$messages, 1)
  headers <- b$table_styling$header$label
  expect_true(any(grepl("Weighted N", headers)))
  expect_true(any(grepl("Total", headers)))
  expect_true(attr(b, "ekb_weighted"))
  expect_error(tbl_baseline(d, "treatment", "age", add_p = TRUE, add_smd = TRUE), "either")
})

test_that("default weighted baseline supports SMD", {
  skip_if_not_installed("gtsummary")
  skip_if_not_installed("survey")
  skip_if_not_installed("smd")
  need_weighting()
  d <- simulate_clinical()
  w <- fit_weighting(d, "treatment", "age")
  a <- suppressMessages(tbl_baseline(d, "treatment", c("age", "sex"), weights = w))
  expect_true("estimate" %in% names(a$table_body))
  expect_true(any(is.finite(a$table_body$estimate)))
})

test_that("outcome weighting affects survival only and emits one message", {
  skip_if_not_installed("gtsummary")
  need_weighting()
  d <- simulate_clinical()
  w <- fit_weighting(d, "treatment", "age")
  args <- list(data = d, by = "treatment", best_response = "response",
    response_rates = c(ORR = "orr", DCR = "dcr"),
    survival = list(OS = c(time = "time", event = "event")), landmarks = c(12, 24), add_p = FALSE)
  a <- do.call(tbl_outcomes, args)
  cap <- capture_messages(do.call(tbl_outcomes, c(args, list(weights = w))))
  b <- cap$value
  expect_length(cap$messages, 1)
  expect_true(attr(b, "ekb_weighted_survival"))
  expect_match(attr(b, "ekb_survival_weighting"), "IPTW")
  expect_true("stat_total" %in% names(b$table_body))
  response_rows <- a$table_body$variable %in% c("response", "orr", "dcr")
  expect_equal(a$table_body[response_rows, ], b$table_body[response_rows, ])
  expect_equal(sum(b$table_body$row_type == "landmark"), 2)
  direct <- summarize_survival(fit_survival(d, "time", "event", "treatment", weights = w$weights), c(12, 24))
  z <- direct$median[direct$median$group == "A", ]
  expected <- ekbMed:::.ekb_outcome_format_survival(z$median, z$conf.low, z$conf.high)
  expect_equal(b$table_body$stat_1[b$table_body$variable == "time" & b$table_body$row_type == "level"], expected)
})

test_that("toxicity tables order levels and count blanks as missing", {
  skip_if_not_installed("gtsummary")
  d <- simulate_clinical()
  a <- tbl_treatment_course(d, "treatment", end_reason = "reason", toxicity = "toxicity",
    toxicity_grade = "grade", toxicity_details = "all", add_p = FALSE,
    level_order = list(grade = c("3", "2", "1")))
  expect_s3_class(a, "gtsummary")
  expect_true("stat_total" %in% names(a$table_body))
  labels <- a$table_body$label[a$table_body$variable == "grade" & a$table_body$row_type != "label"]
  expect_equal(labels[1:3], c("3", "2", "1"))
  expect_true("Missing/unknown" %in% labels)
  expect_false("" %in% labels)
  expect_match(a$table_body$stat_total[a$table_body$variable == "reason" & a$table_body$label == "Missing/unknown"], "60")
})

test_that("all table adapters return expected classes", {
  skip_if_not_installed("gtsummary")
  skip_if_not_installed("flextable")
  skip_if_not_installed("officer")
  skip_if_not_installed("broom.helpers")
  d <- simulate_clinical()
  cox <- fit_cox(d, "time", "event", c("treatment", "age"))
  sg <- subgroup_cox(d, "time", "event", "treatment", "sex")
  for (tbl in list(tbl_cox(cox), tbl_subgroup(sg), tbl_baseline(d, "treatment", "age", add_p = FALSE))) {
    expect_s3_class(tbl, "gtsummary")
    expect_s3_class(as_ekb_flextable(tbl), "flextable")
  }
  expect_s3_class(tbl_survival(summarize_survival(fit_survival(d, "time", "event"))), "tbl_df")
})

test_that("Cox and subgroup forest plots build without pixel snapshots", {
  skip_if_not_installed("patchwork")
  d <- simulate_clinical()
  cox <- fit_cox(d, "time", "event", c("treatment", "age"))
  sg <- subgroup_cox(d, "time", "event", "treatment", "sex")
  a <- plot_cox_forest(cox)
  b <- plot_subgroup_forest(sg)
  expect_s3_class(a, "ggplot")
  expect_s3_class(b, "ggplot")
  expect_s3_class(ggplot2::ggplotGrob(a), "gtable")
  expect_s3_class(ggplot2::ggplotGrob(b), "gtable")
})

test_that("survival plots and risk panels build", {
  skip_if_not_installed("ggsurvfit")
  skip_if_not_installed("patchwork")
  d <- simulate_clinical()
  a <- fit_survival(d, "time", "event", "treatment")
  expect_s3_class(plot_survival(a, risk.table = FALSE), "ggplot")
  p <- plot_survival(a, risk.table = TRUE)
  expect_s3_class(p, "patchwork")
  expect_s3_class(patchwork::patchworkGrob(p), "gtable")
})

test_that("MI forest and regression table reuse stored fits", {
  skip_if_not_installed("mice")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("gtsummary")
  skip_if_not_installed("broom.helpers")
  mi <- make_mi()
  a <- fit_mi_cox(mi, "time", "event", c("treatment", "age"))
  expect_s3_class(plot_cox_forest(a), "ggplot")
  expect_s3_class(tbl_cox(a), "gtsummary")
})

test_that("implicit weighted baseline selection excludes internal weights", {
  skip_if_not_installed("gtsummary")
  skip_if_not_installed("survey")
  d <- simulate_clinical()[, c("treatment", "age", "sex")]
  a <- suppressMessages(tbl_baseline(d, "treatment", weights = rep(1.5, nrow(d)), add_smd = FALSE))
  expect_false(".ekb_weight" %in% a$table_body$variable)
})

test_that("plain gtsummary tables can be converted without ekb metadata", {
  skip_if_not_installed("gtsummary")
  skip_if_not_installed("flextable")
  skip_if_not_installed("officer")
  a <- gtsummary::tbl_summary(data.frame(age = c(40, 50, 60)))
  expect_s3_class(as_ekb_flextable(a), "flextable")
})

test_that("character group headers match the displayed statistics", {
  skip_if_not_installed("gtsummary")
  d <- data.frame(group = c("B", "B", "A", "A"), age = c(80, 80, 20, 20))
  a <- tbl_baseline(d, "group", "age", add_p = FALSE, type = list(age ~ "continuous"))
  header <- a$table_styling$header
  first <- header$label[header$column == "stat_1"]
  body <- a$table_body$stat_1[a$table_body$variable == "age" & a$table_body$row_type == "label"]
  if (grepl("B", first, fixed = TRUE)) expect_match(body, "80") else expect_match(body, "20")
})

test_that("outcome footnotes respect the requested confidence level", {
  skip_if_not_installed("gtsummary")
  d <- simulate_clinical()
  a <- tbl_outcomes(d, "treatment", survival = list(OS = c(time = "time", event = "event")),
    conf.level = 0.9, add_p = FALSE)
  notes <- unlist(a$table_styling$source_note, use.names = FALSE)
  expect_true(any(grepl("90%", notes, fixed = TRUE)))
  expect_false(any(grepl("95%", notes, fixed = TRUE)))
})
