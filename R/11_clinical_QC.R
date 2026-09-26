# Clinical data quality control for ekbMed
#
# Generic, read-only sanity checks for clinical datasets.
# Project-specific rules should be defined outside this file and supplied
# to qc_clinical_data(). The function never modifies the input data.

.ekb_qc_native_missing <- function(x) {
  out <- is.na(x)

  if (is.character(x) || is.factor(x)) {
    x_chr <- trimws(as.character(x))
    out <- out | x_chr == ""
  }

  out
}


.ekb_qc_coded_missing <- function(x, missing_codes) {
  if (!length(missing_codes)) {
    return(rep(FALSE, length(x)))
  }

  x_chr <- trimws(tolower(as.character(x)))
  native <- .ekb_qc_native_missing(x)

  !native & x_chr %in% tolower(missing_codes)
}


.ekb_qc_parse_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }

  suppressWarnings(
    convert_date(
      x,
      invalid = "NA"
    )
  )
}


.ekb_qc_validate_severity <- function(x) {
  x <- toupper(as.character(x))
  allowed <- c("ERROR", "WARNING", "INFO")

  bad <- !x %in% allowed

  if (any(bad)) {
    stop(
      sprintf(
        "Unknown QC severity level(s): %s. Allowed values are ERROR, WARNING, INFO.",
        paste(unique(x[bad]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  x
}


.ekb_qc_count_patients <- function(x) {
  if (!length(x)) {
    return(NA_integer_)
  }

  x <- x[!is.na(x) & x != ""]

  if (!length(x)) {
    return(NA_integer_)
  }

  length(unique(x))
}


.ekb_qc_issue_table <- function() {
  data.frame(
    row = integer(),
    id = character(),
    severity = character(),
    category = character(),
    check = character(),
    variable = character(),
    value = character(),
    related_variable = character(),
    related_value = character(),
    difference_days = integer(),
    details = character(),
    stringsAsFactors = FALSE
  )
}


.ekb_qc_bind_issues <- function(...) {
  x <- list(...)
  x <- x[vapply(x, nrow, integer(1)) > 0L]

  if (!length(x)) {
    return(.ekb_qc_issue_table())
  }

  do.call(rbind, x)
}


#' Inspect clinical data quality without modification
#'
#' Reports missingness, duplicates, date-order violations and other prespecified checks.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param id Optional patient identifier column name.
#' @param date_vars Date columns to inspect for parsing and completeness.
#' @param temporal_rules Data frame with earlier and later column names and optional rule metadata; see Details.
#' @param nonnegative_vars Columns that must not contain negative numeric values.
#' @param required_complete_vars Columns that must not contain missing values.
#' @param missing_codes Character values treated as coded missingness for QC.
#' @param missing_threshold Missingness proportion above which a variable is flagged.
#' @param check_duplicate_ids Whether to flag duplicate identifiers.
#' @param check_duplicate_rows Whether to flag duplicate rows.
#' @param verbose Whether to print a concise QC summary.
#' @return An ekb_clinical_qc list with overview, variable_summary, issues and QC summaries.
#' @details temporal_rules requires columns check, earlier and later. Optional severity is ERROR, WARNING or INFO; optional allow_equal defaults to TRUE. Dates must be parseable. Checks report issues and never modify, impute, drop or repair patient data.
#' @export
#' @examples
#' qc_clinical_data(data.frame(id = c(1, 1, 2), age = c(50, -1, NA)),
#'   id = "id", nonnegative_vars = "age", verbose = FALSE)
qc_clinical_data <- function(
    data,
    id = NULL,
    date_vars = NULL,
    temporal_rules = NULL,
    nonnegative_vars = NULL,
    required_complete_vars = NULL,
    missing_codes = c(
      "Missing",
      "Unknown",
      "Missing/Unknown",
      "Uk",
      "NA",
      "N/A"
    ),
    missing_threshold = 0.20,
    check_duplicate_ids = TRUE,
    check_duplicate_rows = TRUE,
    verbose = TRUE
) {

  dat <- as.data.frame(data)

  if (!is.numeric(missing_threshold) ||
      length(missing_threshold) != 1L ||
      is.na(missing_threshold) ||
      missing_threshold < 0 ||
      missing_threshold > 1) {
    stop("missing_threshold must be a single number between 0 and 1.", call. = FALSE)
  }

  if (!is.null(id)) {
    .ekb_assert_columns(dat, id)
  }

  requested <- unique(c(
    date_vars,
    nonnegative_vars,
    required_complete_vars,
    if (!is.null(temporal_rules)) {
      c(temporal_rules$earlier, temporal_rules$later)
    }
  ))

  requested <- requested[!is.na(requested) & requested != ""]
  .ekb_assert_columns(dat, requested)


  # ---------------------------------------------------------------------------
  # Dataset overview
  # ---------------------------------------------------------------------------

  overview <- data.frame(
    n_rows = nrow(dat),
    n_variables = ncol(dat),
    n_unique_ids = if (is.null(id)) {
      NA_integer_
    } else {
      length(unique(as.character(dat[[id]])[!.ekb_qc_native_missing(dat[[id]])]))
    },
    stringsAsFactors = FALSE
  )


  # ---------------------------------------------------------------------------
  # Variable-level missingness
  # ---------------------------------------------------------------------------

  variable_summary <- do.call(
    rbind,
    lapply(names(dat), function(v) {

      x <- dat[[v]]

      native_missing <- .ekb_qc_native_missing(x)
      coded_missing <- .ekb_qc_coded_missing(x, missing_codes)
      combined_missing <- native_missing | coded_missing

      data.frame(
        variable = v,
        class = paste(class(x), collapse = "/"),
        n_observed = sum(!combined_missing),
        n_unique_observed = length(unique(x[!combined_missing])),
        n_missing = sum(native_missing),
        pct_missing = mean(native_missing),
        n_coded_missing = sum(coded_missing),
        pct_coded_missing = mean(coded_missing),
        n_missing_combined = sum(combined_missing),
        pct_missing_combined = mean(combined_missing),
        stringsAsFactors = FALSE
      )
    })
  )

  high_missingness <- variable_summary[
    variable_summary$pct_missing_combined >= missing_threshold,
    ,
    drop = FALSE
  ]

  high_missingness <- high_missingness[
    order(-high_missingness$pct_missing_combined),
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------------------------
  # ID checks
  # ---------------------------------------------------------------------------

  missing_ids <- data.frame(
    row = integer(),
    id = character(),
    stringsAsFactors = FALSE
  )

  duplicate_ids <- data.frame(
    row = integer(),
    id = character(),
    stringsAsFactors = FALSE
  )

  if (!is.null(id)) {
    id_missing <- .ekb_qc_native_missing(dat[[id]])
    id_chr <- as.character(dat[[id]])

    if (any(id_missing)) {
      missing_ids <- data.frame(
        row = which(id_missing),
        id = id_chr[id_missing],
        stringsAsFactors = FALSE
      )
    }

    if (isTRUE(check_duplicate_ids)) {
      duplicate_idx <- !id_missing & (
        duplicated(id_chr) |
          duplicated(id_chr, fromLast = TRUE)
      )

      if (any(duplicate_idx)) {
        duplicate_ids <- data.frame(
          row = which(duplicate_idx),
          id = id_chr[duplicate_idx],
          stringsAsFactors = FALSE
        )

        duplicate_ids <- duplicate_ids[
          order(duplicate_ids$id, duplicate_ids$row),
          ,
          drop = FALSE
        ]
      }
    }
  }


  # ---------------------------------------------------------------------------
  # Completely duplicated rows
  # ---------------------------------------------------------------------------

  duplicate_rows <- data.frame(
    row = integer(),
    id = character(),
    stringsAsFactors = FALSE
  )

  if (isTRUE(check_duplicate_rows)) {
    duplicate_idx <- duplicated(dat) |
      duplicated(dat, fromLast = TRUE)

    if (any(duplicate_idx)) {
      duplicate_rows <- data.frame(
        row = which(duplicate_idx),
        id = if (is.null(id)) {
          rep(NA_character_, sum(duplicate_idx))
        } else {
          as.character(dat[[id]][duplicate_idx])
        },
        stringsAsFactors = FALSE
      )
    }
  }


  # ---------------------------------------------------------------------------
  # Required complete variables
  # ---------------------------------------------------------------------------

  required_missing <- data.frame(
    row = integer(),
    id = character(),
    variable = character(),
    value = character(),
    stringsAsFactors = FALSE
  )

  if (length(required_complete_vars)) {
    out <- lapply(required_complete_vars, function(v) {

      x <- dat[[v]]
      missing <- .ekb_qc_native_missing(x) |
        .ekb_qc_coded_missing(x, missing_codes)

      if (!any(missing)) {
        return(NULL)
      }

      data.frame(
        row = which(missing),
        id = if (is.null(id)) {
          rep(NA_character_, sum(missing))
        } else {
          as.character(dat[[id]][missing])
        },
        variable = v,
        value = as.character(x[missing]),
        stringsAsFactors = FALSE
      )
    })

    out <- Filter(Negate(is.null), out)

    if (length(out)) {
      required_missing <- do.call(rbind, out)
    }
  }


  # ---------------------------------------------------------------------------
  # Date parsing
  # ---------------------------------------------------------------------------

  invalid_dates <- data.frame(
    row = integer(),
    id = character(),
    variable = character(),
    value = character(),
    stringsAsFactors = FALSE
  )

  if (length(date_vars)) {
    out <- lapply(date_vars, function(v) {

      raw <- dat[[v]]
      parsed <- .ekb_qc_parse_date(raw)

      missing <- .ekb_qc_native_missing(raw) |
        .ekb_qc_coded_missing(raw, missing_codes)

      invalid <- !missing & is.na(parsed)

      if (!any(invalid)) {
        return(NULL)
      }

      data.frame(
        row = which(invalid),
        id = if (is.null(id)) {
          rep(NA_character_, sum(invalid))
        } else {
          as.character(dat[[id]][invalid])
        },
        variable = v,
        value = as.character(raw[invalid]),
        stringsAsFactors = FALSE
      )
    })

    out <- Filter(Negate(is.null), out)

    if (length(out)) {
      invalid_dates <- do.call(rbind, out)
    }
  }


  # ---------------------------------------------------------------------------
  # Temporal plausibility
  #
  # Each rule tests whether:
  #   earlier <= later
  #
  # Required columns:
  #   check, earlier, later
  #
  # Optional columns:
  #   severity    default = "ERROR"
  #   allow_equal default = TRUE
  # ---------------------------------------------------------------------------

  temporal_issues <- data.frame(
    row = integer(),
    id = character(),
    check = character(),
    severity = character(),
    earlier_variable = character(),
    earlier_date = character(),
    later_variable = character(),
    later_date = character(),
    difference_days = integer(),
    stringsAsFactors = FALSE
  )

  if (!is.null(temporal_rules) && nrow(temporal_rules)) {

    temporal_rules <- as.data.frame(temporal_rules)

    required_rule_cols <- c(
      "check",
      "earlier",
      "later"
    )

    if (!all(required_rule_cols %in% names(temporal_rules))) {
      stop(
        "temporal_rules requires columns: check, earlier, later.",
        call. = FALSE
      )
    }

    if (!"severity" %in% names(temporal_rules)) {
      temporal_rules$severity <- "ERROR"
    }

    if (!"allow_equal" %in% names(temporal_rules)) {
      temporal_rules$allow_equal <- TRUE
    }

    temporal_rules$severity <- .ekb_qc_validate_severity(
      temporal_rules$severity
    )

    if (anyNA(temporal_rules$allow_equal)) {
      stop("temporal_rules$allow_equal cannot contain NA.", call. = FALSE)
    }

    out <- lapply(seq_len(nrow(temporal_rules)), function(i) {

      earlier_var <- temporal_rules$earlier[i]
      later_var <- temporal_rules$later[i]

      earlier_date <- .ekb_qc_parse_date(dat[[earlier_var]])
      later_date <- .ekb_qc_parse_date(dat[[later_var]])

      complete <- !is.na(earlier_date) & !is.na(later_date)

      if (isTRUE(temporal_rules$allow_equal[i])) {
        problem <- complete & later_date < earlier_date
      } else {
        problem <- complete & later_date <= earlier_date
      }

      if (!any(problem)) {
        return(NULL)
      }

      data.frame(
        row = which(problem),
        id = if (is.null(id)) {
          rep(NA_character_, sum(problem))
        } else {
          as.character(dat[[id]][problem])
        },
        check = temporal_rules$check[i],
        severity = temporal_rules$severity[i],
        earlier_variable = earlier_var,
        earlier_date = as.character(earlier_date[problem]),
        later_variable = later_var,
        later_date = as.character(later_date[problem]),
        difference_days = as.integer(
          later_date[problem] - earlier_date[problem]
        ),
        stringsAsFactors = FALSE
      )
    })

    out <- Filter(Negate(is.null), out)

    if (length(out)) {
      temporal_issues <- do.call(rbind, out)
    }
  }


  # ---------------------------------------------------------------------------
  # Numeric plausibility
  #
  # Variables supplied in nonnegative_vars are expected to be numeric and >= 0.
  # ---------------------------------------------------------------------------

  invalid_numeric_values <- data.frame(
    row = integer(),
    id = character(),
    variable = character(),
    value = character(),
    stringsAsFactors = FALSE
  )

  negative_values <- data.frame(
    row = integer(),
    id = character(),
    variable = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )

  if (length(nonnegative_vars)) {

    invalid_numeric_out <- list()
    negative_out <- list()

    for (v in nonnegative_vars) {

      raw <- dat[[v]]

      missing <- .ekb_qc_native_missing(raw) |
        .ekb_qc_coded_missing(raw, missing_codes)

      if (is.numeric(raw)) {
        value <- raw
        invalid_numeric <- rep(FALSE, length(raw))
      } else {
        value <- suppressWarnings(as.numeric(as.character(raw)))
        invalid_numeric <- !missing & is.na(value)
      }

      if (any(invalid_numeric)) {
        invalid_numeric_out[[v]] <- data.frame(
          row = which(invalid_numeric),
          id = if (is.null(id)) {
            rep(NA_character_, sum(invalid_numeric))
          } else {
            as.character(dat[[id]][invalid_numeric])
          },
          variable = v,
          value = as.character(raw[invalid_numeric]),
          stringsAsFactors = FALSE
        )
      }

      negative <- !is.na(value) & value < 0

      if (any(negative)) {
        negative_out[[v]] <- data.frame(
          row = which(negative),
          id = if (is.null(id)) {
            rep(NA_character_, sum(negative))
          } else {
            as.character(dat[[id]][negative])
          },
          variable = v,
          value = value[negative],
          stringsAsFactors = FALSE
        )
      }
    }

    invalid_numeric_out <- Filter(Negate(is.null), invalid_numeric_out)
    negative_out <- Filter(Negate(is.null), negative_out)

    if (length(invalid_numeric_out)) {
      invalid_numeric_values <- do.call(rbind, invalid_numeric_out)
    }

    if (length(negative_out)) {
      negative_values <- do.call(rbind, negative_out)
    }
  }


  # ---------------------------------------------------------------------------
  # Unified issue table
  # ---------------------------------------------------------------------------

  issue_missing_ids <- .ekb_qc_issue_table()
  if (nrow(missing_ids)) {
    issue_missing_ids <- data.frame(
      row = missing_ids$row,
      id = missing_ids$id,
      severity = "ERROR",
      category = "identifier",
      check = "Missing patient ID",
      variable = id %||% NA_character_,
      value = NA_character_,
      related_variable = NA_character_,
      related_value = NA_character_,
      difference_days = NA_integer_,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issue_duplicate_ids <- .ekb_qc_issue_table()
  if (nrow(duplicate_ids)) {
    issue_duplicate_ids <- data.frame(
      row = duplicate_ids$row,
      id = duplicate_ids$id,
      severity = "ERROR",
      category = "identifier",
      check = "Duplicated patient ID",
      variable = id %||% NA_character_,
      value = duplicate_ids$id,
      related_variable = NA_character_,
      related_value = NA_character_,
      difference_days = NA_integer_,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issue_duplicate_rows <- .ekb_qc_issue_table()
  if (nrow(duplicate_rows)) {
    issue_duplicate_rows <- data.frame(
      row = duplicate_rows$row,
      id = duplicate_rows$id,
      severity = "WARNING",
      category = "duplicate",
      check = "Completely duplicated row",
      variable = NA_character_,
      value = NA_character_,
      related_variable = NA_character_,
      related_value = NA_character_,
      difference_days = NA_integer_,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issue_required_missing <- .ekb_qc_issue_table()
  if (nrow(required_missing)) {
    issue_required_missing <- data.frame(
      row = required_missing$row,
      id = required_missing$id,
      severity = "ERROR",
      category = "missingness",
      check = paste0("Missing required variable: ", required_missing$variable),
      variable = required_missing$variable,
      value = required_missing$value,
      related_variable = NA_character_,
      related_value = NA_character_,
      difference_days = NA_integer_,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issue_invalid_dates <- .ekb_qc_issue_table()
  if (nrow(invalid_dates)) {
    issue_invalid_dates <- data.frame(
      row = invalid_dates$row,
      id = invalid_dates$id,
      severity = "ERROR",
      category = "date",
      check = paste0("Invalid date value: ", invalid_dates$variable),
      variable = invalid_dates$variable,
      value = invalid_dates$value,
      related_variable = NA_character_,
      related_value = NA_character_,
      difference_days = NA_integer_,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issue_temporal <- .ekb_qc_issue_table()
  if (nrow(temporal_issues)) {
    issue_temporal <- data.frame(
      row = temporal_issues$row,
      id = temporal_issues$id,
      severity = temporal_issues$severity,
      category = "temporal",
      check = temporal_issues$check,
      variable = temporal_issues$earlier_variable,
      value = temporal_issues$earlier_date,
      related_variable = temporal_issues$later_variable,
      related_value = temporal_issues$later_date,
      difference_days = temporal_issues$difference_days,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issue_invalid_numeric <- .ekb_qc_issue_table()
  if (nrow(invalid_numeric_values)) {
    issue_invalid_numeric <- data.frame(
      row = invalid_numeric_values$row,
      id = invalid_numeric_values$id,
      severity = "ERROR",
      category = "numeric",
      check = paste0(
        "Non-numeric value in numeric QC variable: ",
        invalid_numeric_values$variable
      ),
      variable = invalid_numeric_values$variable,
      value = invalid_numeric_values$value,
      related_variable = NA_character_,
      related_value = NA_character_,
      difference_days = NA_integer_,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issue_negative <- .ekb_qc_issue_table()
  if (nrow(negative_values)) {
    issue_negative <- data.frame(
      row = negative_values$row,
      id = negative_values$id,
      severity = "ERROR",
      category = "numeric",
      check = "Negative value in non-negative variable",
      variable = negative_values$variable,
      value = as.character(negative_values$value),
      related_variable = NA_character_,
      related_value = NA_character_,
      difference_days = NA_integer_,
      details = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  issues <- .ekb_qc_bind_issues(
    issue_missing_ids,
    issue_duplicate_ids,
    issue_duplicate_rows,
    issue_required_missing,
    issue_invalid_dates,
    issue_temporal,
    issue_invalid_numeric,
    issue_negative
  )


  # ---------------------------------------------------------------------------
  # Issue summary
  # ---------------------------------------------------------------------------

  summary <- data.frame(
    severity = character(),
    check = character(),
    n_issues = integer(),
    n_patients = integer(),
    stringsAsFactors = FALSE
  )

  if (nrow(issues)) {
    keys <- unique(issues[, c("severity", "check"), drop = FALSE])

    summary <- do.call(
      rbind,
      lapply(seq_len(nrow(keys)), function(i) {

        idx <- issues$severity == keys$severity[i] &
          issues$check == keys$check[i]

        data.frame(
          severity = keys$severity[i],
          check = keys$check[i],
          n_issues = sum(idx),
          n_patients = if (is.null(id)) {
            NA_integer_
          } else {
            .ekb_qc_count_patients(issues$id[idx])
          },
          stringsAsFactors = FALSE
        )
      })
    )
  }

  if (nrow(high_missingness)) {
    summary <- rbind(
      summary,
      data.frame(
        severity = "INFO",
        check = sprintf(
          "Variables with >= %.0f%% missing/coded-missing values",
          100 * missing_threshold
        ),
        n_issues = nrow(high_missingness),
        n_patients = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  }

  if (nrow(summary)) {
    severity_order <- c("ERROR", "WARNING", "INFO")

    summary$.severity_order <- match(
      summary$severity,
      severity_order
    )

    summary <- summary[
      order(
        summary$.severity_order,
        -summary$n_issues,
        summary$check
      ),
      ,
      drop = FALSE
    ]

    summary$.severity_order <- NULL
    rownames(summary) <- NULL
  }


  # ---------------------------------------------------------------------------
  # Patient-level summary
  # ---------------------------------------------------------------------------

  patient_summary <- data.frame(
    id = character(),
    n_issues = integer(),
    n_errors = integer(),
    n_warnings = integer(),
    n_info = integer(),
    checks = character(),
    stringsAsFactors = FALSE
  )

  if (!is.null(id) && nrow(issues)) {

    patient_issues <- issues[
      !is.na(issues$id) & issues$id != "",
      ,
      drop = FALSE
    ]

    if (nrow(patient_issues)) {
      ids <- unique(patient_issues$id)

      patient_summary <- do.call(
        rbind,
        lapply(ids, function(z) {

          d <- patient_issues[
            patient_issues$id == z,
            ,
            drop = FALSE
          ]

          data.frame(
            id = z,
            n_issues = nrow(d),
            n_errors = sum(d$severity == "ERROR"),
            n_warnings = sum(d$severity == "WARNING"),
            n_info = sum(d$severity == "INFO"),
            checks = paste(unique(d$check), collapse = "; "),
            stringsAsFactors = FALSE
          )
        })
      )

      patient_summary <- patient_summary[
        order(
          -patient_summary$n_errors,
          -patient_summary$n_warnings,
          -patient_summary$n_issues,
          patient_summary$id
        ),
        ,
        drop = FALSE
      ]

      rownames(patient_summary) <- NULL
    }
  }


  # ---------------------------------------------------------------------------
  # Return object
  # ---------------------------------------------------------------------------

  result <- structure(
    list(
      overview = .ekb_as_tibble(overview),
      summary = .ekb_as_tibble(summary),
      issues = .ekb_as_tibble(issues),
      patient_summary = .ekb_as_tibble(patient_summary),
      variable_summary = .ekb_as_tibble(variable_summary),
      high_missingness = .ekb_as_tibble(high_missingness),
      missing_ids = .ekb_as_tibble(missing_ids),
      duplicate_ids = .ekb_as_tibble(duplicate_ids),
      duplicate_rows = .ekb_as_tibble(duplicate_rows),
      required_missing = .ekb_as_tibble(required_missing),
      invalid_dates = .ekb_as_tibble(invalid_dates),
      temporal_issues = .ekb_as_tibble(temporal_issues),
      invalid_numeric_values = .ekb_as_tibble(invalid_numeric_values),
      negative_values = .ekb_as_tibble(negative_values),
      settings = list(
        id = id,
        date_vars = date_vars,
        temporal_rules = temporal_rules,
        nonnegative_vars = nonnegative_vars,
        required_complete_vars = required_complete_vars,
        missing_codes = missing_codes,
        missing_threshold = missing_threshold,
        check_duplicate_ids = check_duplicate_ids,
        check_duplicate_rows = check_duplicate_rows
      ),
      call = match.call()
    ),
    class = "ekb_clinical_qc"
  )

  if (isTRUE(verbose)) {
    print(result)
  }

  invisible(result)
}


#' Print a clinical quality-control summary
#'
#' Displays the summary of an ekb_clinical_qc object.
#'
#' @param x An ekb_clinical_qc object.
#' @param ... Additional arguments passed to the underlying function; see Details.
#' @return The input object, invisibly.
#' @details Additional arguments are reserved for compatibility with print.
#' @export
#' @examples
#' q <- qc_clinical_data(data.frame(id = 1:3), id = "id", verbose = FALSE)
#' print(q)
print.ekb_clinical_qc <- function(x, ...) {

  cat("\nClinical data QC\n")
  cat("================\n")

  cat(
    "Rows:      ",
    format(x$overview$n_rows, big.mark = ","),
    "\n",
    sep = ""
  )

  cat(
    "Variables: ",
    format(x$overview$n_variables, big.mark = ","),
    "\n",
    sep = ""
  )

  if (!is.na(x$overview$n_unique_ids)) {
    cat(
      "Unique IDs:",
      format(x$overview$n_unique_ids, big.mark = ","),
      "\n"
    )
  }

  cat("\n")

  if (!nrow(x$summary)) {
    cat("No predefined QC issues detected.\n")
    return(invisible(x))
  }

  print(
    as.data.frame(x$summary),
    row.names = FALSE
  )

  invisible(x)
}
