simulate_clinical <- function(n = 240, seed = 410) {
  set.seed(seed)
  age <- rnorm(n, 60, 9)
  sex <- factor(rep(c("F", "M"), length.out = n))
  treatment <- factor(ifelse(runif(n) < plogis(0.025 * (age - 60)), "B", "A"), levels = c("A", "B"))
  failure <- rexp(n, 0.035 * exp(0.3 * (treatment == "B") + 0.02 * (age - 60)))
  censor <- runif(n, 18, 72)
  response <- factor(sample(c("CR", "PR", "SD", "PD"), n, TRUE), levels = c("CR", "PR", "SD", "PD"))
  data.frame(time = pmin(failure, censor), event = as.integer(failure <= censor),
    treatment = treatment, age = age, sex = sex, constant = 1,
    response = response, orr = response %in% c("CR", "PR"),
    dcr = response %in% c("CR", "PR", "SD"),
    reason = rep(c("Completed", "Toxicity", "Progression", ""), length.out = n),
    toxicity = rep(c("Yes", "No", "Yes", "No"), length.out = n),
    grade = rep(c("3", "", "1", "2"), length.out = n))
}

make_mi <- function(impute_sex = FALSE) {
  d <- simulate_clinical()
  d$age[seq(5, nrow(d), 11)] <- NA_real_
  vars <- "age"
  if (impute_sex) {
    d$sex[seq(8, nrow(d), 17)] <- NA
    vars <- c(vars, "sex")
  }
  impute_clinical(d[, c("time", "event", "treatment", "age", "sex")],
    impute_vars = vars, predictor_vars = "sex", treatment = "treatment",
    time = "time", event = "event", m = 2, maxit = 2, seed = 71)
}

capture_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, message = function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(value = value, messages = messages)
}

need_weighting <- function() {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("cobalt")
}
