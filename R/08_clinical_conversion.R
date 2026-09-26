# SAP-driven clinical conversion helpers

#' Convert complete and partial clinical dates
#'
#' Applies explicit partial-date rules while retaining optional audit details.
#'
#' @param x Input vector or object.
#' @param year_only_month Month assigned to year-only dates under the explicit conversion rule.
#' @param month_only_day Day assigned to incomplete dates under the explicit conversion rule.
#' @param return_details Whether to return a tibble containing values and applied conversion rules.
#' @param invalid Return NA or raise an error for unparseable dates.
#' @return A Date vector, or a tibble containing original, date, imputed, imputed_component, rule and invalid.
#' @details Year-only dates use the specified month and day; year-month dates use the specified day. Parsing retains the legacy dmy-first, ymd-second priority. Confirm these explicit deterministic rules against the study analysis plan; this is not multiple imputation.
#' @export
#' @examples
#' convert_date(c("2020", "2020-05", "21.05.2020"), return_details = TRUE)
convert_date <- function(x,
                         year_only_month = 6,
                         month_only_day = 15,
                         return_details = FALSE,
                         invalid = c("NA", "error")) {
  .ekb_require("lubridate")
  .ekb_require("tibble")
  invalid <- match.arg(invalid)
  x_chr <- as.character(x)
  work <- trimws(x_chr)
  work[work == ""] <- NA_character_

  component <- rep(NA_character_, length(work))
  rule <- rep(NA_character_, length(work))

  year_only <- !is.na(work) & grepl("^\\d{4}$", work)
  if (any(year_only)) {
    work[year_only] <- sprintf("%s-%02d-%02d", work[year_only], as.integer(year_only_month), as.integer(month_only_day))
    component[year_only] <- "month+day"
    rule[year_only] <- sprintf("year only -> month %02d and day %02d", as.integer(year_only_month), as.integer(month_only_day))
  }

  month_only <- !is.na(work) & grepl("^\\d{4}-\\d{2}$", work)
  if (any(month_only)) {
    work[month_only] <- sprintf("%s-%02d", work[month_only], as.integer(month_only_day))
    component[month_only] <- "day"
    rule[month_only] <- sprintf("year-month -> day %02d", as.integer(month_only_day))
  }

  # Preserve legacy emR parsing priority: dmy first, then ymd for entries dmy cannot parse.
  dmy <- suppressWarnings(lubridate::dmy(work, quiet = TRUE))
  ymd <- suppressWarnings(lubridate::ymd(work, quiet = TRUE))
  out <- dmy
  out[is.na(out)] <- ymd[is.na(out)]

  invalid_idx <- !is.na(work) & is.na(out)
  if (invalid == "error" && any(invalid_idx)) {
    stop(sprintf("Could not parse date value(s): %s", paste(unique(x_chr[invalid_idx]), collapse = ", ")), call. = FALSE)
  }

  if (!isTRUE(return_details)) return(as.Date(out))
  tibble::tibble(
    original = x_chr,
    date = as.Date(out),
    imputed = !is.na(component),
    imputed_component = component,
    rule = rule,
    invalid = invalid_idx
  )
}

#' Convert clinical response categories
#'
#' Applies explicit response mappings and treatment-duration rules for mixed responses.
#'
#' @param x Input vector or object.
#' @param trtdur Numeric treatment duration in days, aligned with the response vector.
#' @param na_values Response values explicitly mapped to NA.
#' @param mapping Named character vector of source-to-target response mappings.
#' @param mixed_response Category identifying mixed responses.
#' @param mixed_cutoff_days Treatment-duration cutoff for the mixed-response conversion.
#' @param mixed_if_long Response assigned when treatment duration exceeds the cutoff.
#' @param mixed_if_short Response assigned when treatment duration is at or below the cutoff.
#' @param return_details Whether to return a tibble containing values and applied conversion rules.
#' @return A character vector, or a tibble with original, converted and rule.
#' @details Default mappings retain the validated legacy conventions NC to SD and NED to CR. Mixed response becomes SD only when duration is strictly above 90 days, otherwise PD; missing duration becomes NA. These are configurable study conventions, not universal response criteria.
#' @export
#' @examples
#' convert_response(c("NC", "NED", "MR"), trtdur = c(NA, NA, 100))
convert_response <- function(x,
                             trtdur = NULL,
                             na_values = c("", "NB"),
                             mapping = c(NC = "SD", NED = "CR"),
                             mixed_response = "MR",
                             mixed_cutoff_days = 90,
                             mixed_if_long = "SD",
                             mixed_if_short = "PD",
                             return_details = FALSE) {
  .ekb_require("tibble")
  original <- as.character(x)
  out <- original
  rule <- rep(NA_character_, length(out))

  na_idx <- !is.na(out) & out %in% na_values
  out[na_idx] <- NA_character_
  rule[na_idx] <- "mapped to NA by SAP rule"

  if (length(mapping)) {
    for (src in names(mapping)) {
      idx <- !is.na(out) & out == src
      out[idx] <- unname(mapping[[src]])
      rule[idx] <- sprintf("%s -> %s", src, unname(mapping[[src]]))
    }
  }

  mr_idx <- !is.na(out) & out == mixed_response
  if (any(mr_idx)) {
    if (is.null(trtdur)) stop("trtdur is required when mixed-response values are present.", call. = FALSE)
    if (length(trtdur) != length(out)) stop("trtdur must have the same length as x.", call. = FALSE)
    long <- mr_idx & !is.na(trtdur) & trtdur > mixed_cutoff_days
    short <- mr_idx & !is.na(trtdur) & trtdur <= mixed_cutoff_days
    unknown <- mr_idx & is.na(trtdur)
    out[long] <- mixed_if_long
    out[short] <- mixed_if_short
    out[unknown] <- NA_character_
    rule[long] <- sprintf("%s and treatment duration > %s days -> %s", mixed_response, mixed_cutoff_days, mixed_if_long)
    rule[short] <- sprintf("%s and treatment duration <= %s days -> %s", mixed_response, mixed_cutoff_days, mixed_if_short)
    rule[unknown] <- sprintf("%s with missing treatment duration -> NA", mixed_response)
  }

  if (!isTRUE(return_details)) return(out)
  tibble::tibble(original = original, converted = out, rule = rule)
}

#' Calculate a date-based follow-up interval
#'
#' Calculates duration with an explicit inclusive-day convention.
#'
#' @param start_date Start date or vector coercible to Date.
#' @param event_date End date or vector coercible to Date.
#' @param inclusive Whether to include the starting day by adding one day.
#' @param unit Duration unit passed to [lubridate::time_length()].
#' @param digits Number of decimal places for estimates.
#' @return A numeric duration vector rounded to digits.
#' @export
#' @examples
#' calc_survival_interval("2020-01-01", "2020-02-01", unit = "days")
calc_survival_interval <- function(start_date,
                                   event_date,
                                   inclusive = TRUE,
                                   unit = "months",
                                   digits = 2) {
  .ekb_require("lubridate")
  start_date <- as.Date(start_date)
  event_date <- as.Date(event_date)
  if (length(start_date) != length(event_date) && length(start_date) != 1L && length(event_date) != 1L)
    stop("Date vectors must have equal lengths or one must be scalar.", call. = FALSE)
  raw_days <- as.numeric(event_date - start_date)
  if (any(raw_days < 0, na.rm = TRUE)) stop("At least one event_date is before start_date.", call. = FALSE)
  days <- raw_days + if (isTRUE(inclusive)) 1 else 0
  duration <- lubridate::duration(days, units = "days")
  round(lubridate::time_length(duration, unit = unit), digits)
}

# Backwards-compatible name. The +1 day behavior is now explicit via inclusive = TRUE.
#' Calculate follow-up using legacy argument names
#'
#' Compatibility wrapper for calc_survival_interval using months and two decimals.
#'
#' @param startDate Legacy spelling of start_date.
#' @param eventDate Legacy spelling of event_date.
#' @param inclusive Whether to include the starting day by adding one day.
#' @return A numeric vector of follow-up durations in months.
#' @family compatibility wrappers
#' @export
#' @examples
#' calc_survival("2020-01-01", "2020-02-01")
calc_survival <- function(startDate, eventDate, inclusive = TRUE) {
  calc_survival_interval(startDate, eventDate, inclusive = inclusive, unit = "months", digits = 2)
}

#' Calculate ORR or DCR with exact intervals
#'
#' Uses evaluable responses and an exact binomial confidence interval.
#'
#' @param response Character or factor vector of response categories.
#' @param type Response rate to calculate: ORR or DCR.
#' @param evaluable Response categories included in the denominator; defaults to CR, PR, SD and PD.
#' @param conf.level Confidence level between zero and one.
#' @param na.rm Whether to discard NA responses. FALSE with any NA returns an NA summary.
#' @return A data frame with type, n, events, estimate, conf.low and conf.high.
#' @details ORR counts CR and PR; DCR counts CR, PR and SD. The default denominator includes only CR/PR/SD/PD, excluding missing and other unknown codes. Zero evaluable responses return NA estimates. Confidence intervals use stats::binom.test.
#' @export
#' @examples
#' calc_response_rate(c("CR", "PR", "SD", "PD", NA, "Unknown"))
calc_response_rate <- function(response,
                               type = c("ORR", "DCR"),
                               evaluable = c("CR", "PR", "SD", "PD"),
                               conf.level = 0.95,
                               na.rm = TRUE) {
  type <- match.arg(type)
  x <- as.character(response)
  if (!isTRUE(na.rm) && anyNA(x))
    return(data.frame(type = type, n = NA_integer_, events = NA_integer_, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_))
  if (isTRUE(na.rm)) x <- x[!is.na(x)]
  x <- x[x %in% evaluable]
  success <- if (type == "ORR") x %in% c("CR", "PR") else x %in% c("CR", "PR", "SD")
  n <- length(success)
  events <- sum(success)
  if (!n) return(data.frame(type = type, n = 0L, events = 0L, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_))
  ci <- stats::binom.test(events, n, conf.level = conf.level)$conf.int
  data.frame(type = type, n = n, events = events, estimate = events / n, conf.low = ci[1], conf.high = ci[2])
}

#' Calculate objective response rate
#'
#' Thin convenience wrapper for calc_response_rate with type ORR.
#'
#' @param response Character or factor vector of response categories.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return A one-row data frame of the response estimate and exact confidence interval.
#' @details Additional arguments are passed to calc_response_rate. CR and PR count as successes; the denominator defaults to CR/PR/SD/PD.
#' @export
#' @examples
#' calc_ORR(c("CR", "PR", "SD", "PD"))
calc_ORR <- function(response, ...) calc_response_rate(response, type = "ORR", ...)
#' Calculate disease control rate
#'
#' Thin convenience wrapper for calc_response_rate with type DCR.
#'
#' @param response Character or factor vector of response categories.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return A one-row data frame of the response estimate and exact confidence interval.
#' @details Additional arguments are passed to calc_response_rate. CR, PR and SD count as successes; the denominator defaults to CR/PR/SD/PD.
#' @export
#' @examples
#' calc_DCR(c("CR", "PR", "SD", "PD"))
calc_DCR <- function(response, ...) calc_response_rate(response, type = "DCR", ...)
